"""
This script is designed to help with the processing of bgen files for WGS data.
"""
import json
from pathlib import Path
from typing import Union

import pandas as pd


def run_splitter(coordinate_path: Union[pd.DataFrame, Path], gene_dict: Union[dict, Path], chunk_size: int = 30):
    """
    Run the splitter

    :param coordinate_path: file with the original coordinates
    :param gene_dict: dictionary object with gene positions
    :param chunk_size: size chunk for chunking
    """

    # read in the original file
    original_file = pd.read_csv(coordinate_path, sep='\t')
    # do the splitting
    modified_file = split_coordinates_file(original_file, gene_dict, chunk_size)
    # output the coordinate file with chunking
    modified_file.to_csv(coordinate_path, sep='\t')


def split_coordinates_file(coordinates_file: pd.DataFrame, gene_dict: Union[dict, Path],
                           chunk_size: int = 30) -> pd.DataFrame:
    """
    Splits the coordinates file into chunks based on the specified chunk size and gene dictionary.

    :param coordinates_file: DataFrame containing the coordinates data
    :param gene_dict: Dictionary containing gene information with genomic locations
    :param chunk_size: Size of each chunk in base pairs. Default is 30 (interpreted as 30Mb)
    :return: DataFrame with the coordinates split into chunks
    """

    # read in the gene dictionary
    if gene_dict == Path:
        with open(gene_dict, 'r') as f:
            gene_dict = json.load(f)

    # first, if we are using a 30Mb chunk then we need to convert it to Mb
    if chunk_size == 30:
        print("Creating 30Mb chunks for bgens")
        # convert chunk size to mega base
        chunk_size = chunk_size * 1000000
    # if we are not (for testing purposes, I imagine) then use that instead
    else:
        print(f"Using a custom chunk size of {chunk_size}")

    # order by position to make them more readable
    coordinates_file = sort_coordinates_by_position(coordinates_file).reset_index(drop=True)

    # re-arrange our gene dictionary for easy processing
    gene_records = []
    for gene, details in gene_dict.items():
        loc = details['genomic_location']
        gene_records.append({
            'gene': gene,
            'chrom': f"chr{loc['chromosome']}",
            'start': loc['start'],
            'end': loc['end']
        })
    gene_df = pd.DataFrame(gene_records)
    gene_df = sort_coordinates_by_position(gene_df).reset_index(drop=True)

    # start the chunking process
    chunks = []

    # we want to do this separately for each chromosome
    for chrom in sorted(coordinates_file['chrom'].unique(), key=lambda x: int(x.replace('chr', ''))):
        # get the chromosome data that we are working on
        chrom_df = coordinates_file[coordinates_file['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)
        # also subset the gene data to get the working chromosome
        chrom_genes = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

        # if both are empty then something is wrong
        if chrom_df.empty or chrom_genes.empty:
            # throw an error
            raise ValueError("Chromosome data or gene data is empty for chromosome: " + chrom)

        # work out the very start of the chromosome
        current_start = chrom_df['start'].min()
        # this will be our first chunk
        chunk_number = 1
        # also work out the end of the chromosome
        chrom_max_end = chrom_df['end'].max()

        # while loop ensures that we stay within the chromosome boundary
        while current_start < chrom_max_end:
            # add our chunk to the current position
            position_to_check = current_start + chunk_size

            # check if within gene
            gene_context = find_gene_context(chrom, position_to_check, gene_df)

            # if it's not within a gene, then we can cut the chunk here
            if gene_context['status'] != 'within':
                # fallback: use chunk_size increment ensuring it's smaller than the
                # end of the chromosome
                current_end = min(current_start + chunk_size, chrom_max_end)
            else:
                # if we are within a gene, then let's get that gene
                gene_name = gene_context['genes'][0]
                # then let's find a space downstream between this gene and the
                # next gene - directly in-between the two genes
                downstream = find_position_after_within_gene(chrom, gene_name, gene_df)
                # if there is no such gene (shouldn't happen until we run out of genes)
                # then use the chromosome end position as the end
                if downstream is None:
                    current_end = chrom_max_end
                # otherwise, let's use this new in-between position as the end of our chunk
                else:
                    current_end = min(downstream['between_position'], chrom_max_end)

            # Find all rows in coordinates_file that overlap with this chunk
            overlapping_rows = chrom_df[(chrom_df['end'] >= current_start) & (chrom_df['start'] <= current_end)]

            # For each overlapping row, create a new entry with this chunk label
            for index, overlap_row in overlapping_rows.iterrows():
                chunks.append({
                    'chrom': f"{chrom}_chunk{chunk_number}",
                    'start': overlap_row['start'],
                    'end': overlap_row['end'],
                    'chunk_start': current_start,
                    'chunk_end': current_end,
                    'vcf_prefix': overlap_row['vcf_prefix'],
                    'output_bcf': overlap_row['output_bcf'],
                    'output_bcf_idx': overlap_row['output_bcf_idx'],
                    'output_vep': overlap_row['output_vep'],
                    'output_vep_idx': overlap_row['output_vep_idx']
                })

            # move to next chunk
            current_start = current_end + 1
            chunk_number += 1

    # make a new dataframe with the new chunks
    chunk_df = pd.DataFrame(chunks)

    return chunk_df


def find_chunk_for_position(chrom_df: pd.DataFrame, position: int) -> pd.Series:
    """
    Finds the chunk for a given position within our per-chromosome DataFrame. This is just a way
    to locate a chunk in a given dataframe that pinpoints the position that we are working with.
    In other words: 1 < our position (3) > 5 = you are in the chunk that's 1-5

    :param chrom_df: DataFrame containing chromosome data with 'start' and 'end' columns.
    :param position: The position to find the chunk for.
    :return: A Series representing the chunk that contains the given position, or None if no chunk is found.
    """
    # find the chunk (i.e. row) that matches the position that we are working with
    match = chrom_df[(chrom_df['start'] <= position) & (chrom_df['end'] >= position)]
    if not match.empty:
        return match.iloc[0]
    else:
        # I believe we should always have a matching chunk (as there is already something in place for when
        # we reach the limit of a chromosome), so this should throw an error
        raise ValueError("No chunk found for the given position: " + str(position))


def find_position_after_within_gene(chrom: str, gene_name: str, gene_df: pd.DataFrame) -> dict:
    """
    If our position falls within a gene, this is a way to find a position between this gene and the next gene that doesn't
    overlap, and take the half-way point as our new break-point.
    In other words, we want to create a break in a "junk" region between two non-overlapping genes and this should
    help us do so.

    :param chrom: Chromosome identifier (e.g., 'chr1')
    :param gene_name: Name of the gene to find the position after (i.e. downstream of gene_name)
    :param gene_df: DataFrame containing gene information with 'chrom', 'gene', 'start', and 'end' columns
    :return: A dictionary with the position between the current gene and the next non-overlapping gene, or None if no such position is found
    """
    # Get all genes on this chromosome sorted by the start position (this data comes from our gene dictionary)
    genes_on_chrom = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

    # Find the index of the current gene that we are working with (i.e. our position is within a gene, what's the
    # index of that gene?
    gene_idx = genes_on_chrom[genes_on_chrom['gene'] == gene_name].index
    if gene_idx.empty:
        # this should never happen so add error catching here
        raise ValueError("Our position did not match a gene... something weird going on with: " + str(gene_name))

    # find the end of this gene
    gene_idx = gene_idx[0]
    current_end = genes_on_chrom.loc[gene_idx, 'end']

    # Look ahead to find the first downstream gene that doesn't overlap
    for i in range(gene_idx + 1, len(genes_on_chrom)):
        next_start = genes_on_chrom.loc[i, 'start']
        next_gene = genes_on_chrom.loc[i, 'gene']

        # once we find a downstream gene, we want to record the new position that is in-between these
        # two genes and also maybe the names of our genes...
        if next_start > current_end:
            midpoint = (current_end + next_start) // 2
            return {
                'between_position': midpoint,
                'within_gene': genes_on_chrom.loc[gene_idx, 'gene'],
                'next_gene': next_gene
            }

        # Update current_end in case of overlapping gene
        current_end = max(current_end, genes_on_chrom.loc[i, 'end'])

    # If no non-overlapping downstream gene found
    raise ValueError("Couldn't find any downstream genes... something odd going on with: " + str(
        genes_on_chrom.loc[gene_idx, 'gene']))


def find_gene_context(chrom: str, position: int, gene_df: pd.DataFrame) -> dict:
    """
    Finds the gene context for a given position within a chromosome. In other words, we have a given position (which
    might be a possible break point), but we want to know whether this position falls within a gene, or whether it
    falls between two genes.

    :param chrom: Chromosome identifier (e.g., 'chr1')
    :param position: The position to find the gene context for.
    :param gene_df: DataFrame containing gene information with 'chrom', 'gene', 'start', and 'end' columns.
    :return: A dictionary with the gene context, indicating whether the position is within a gene or between genes, and the closest upstream and downstream genes if applicable.
    """

    # Filter to just the chromosome of interest
    genes_on_chrom = gene_df[gene_df['chrom'] == chrom]

    # Check if position falls within any gene
    within = genes_on_chrom[(genes_on_chrom['start'] <= position) & (genes_on_chrom['end'] >= position)]
    if not within.empty:
        return {'status': 'within', 'genes': within['gene'].tolist()}

    # If not within, find closest upstream and downstream genes
    # let's get all the upstream genes
    upstream = genes_on_chrom[genes_on_chrom['end'] < position]
    # let's get all the downstream genes
    downstream = genes_on_chrom[genes_on_chrom['start'] > position]
    # then let's find the closest upstream gene
    closest_up = upstream.iloc[upstream['end'].sub(position).abs().argmin()] if not upstream.empty else None
    # and then let's find the closest downstream gene
    closest_down = downstream.iloc[downstream['start'].sub(position).abs().argmin()] if not downstream.empty else None

    # if we are already between two genes then return this as our status, but also let us know what are the
    # two genes that we are sandwiched between
    return {
        'status': 'between',
        'closest_upstream': closest_up['gene'] if closest_up is not None else None,
        'closest_downstream': closest_down['gene'] if closest_down is not None else None
    }


def sort_coordinates_by_position(df: pd.DataFrame) -> pd.DataFrame:
    """
    Simple function to sort a dataframe based on chromosome number and start position
    """

    df['chrom_num'] = df['chrom'].replace('chr', '')
    return df.sort_values(by=['chrom_num', 'start']).drop(columns='chrom_num')
