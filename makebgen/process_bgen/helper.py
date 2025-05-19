"""
This script is designed to help with the processing of bgen files for WGS data.
"""
import argparse
import json
from pathlib import Path
from typing import Union

import pandas as pd
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process bgen files for WGS data.")
    parser.add_argument("--coordinate_path", type=str, required=True, help="Path to the coordinate file.")
    parser.add_argument("--gene_dict", type=str, required=True, help="Path to the gene dictionary file.")
    parser.add_argument("--chunk_size", type=int, default=3, required=False, help="Chunk size in Mb (default: 3Mb).")
    parser.add_argument("--output_path", type=str, default="chunked_files", required=False, help="Path to the output directory (default: 'chunked_files').")
    return parser.parse_args()


def run_splitter(coordinate_path: [Path], gene_dict: Union[dict, Path, str], chunk_size: int,
                 output_path: Path) -> None:
    """
    Run the splitter.

    The aim of this function is to take the coordinate file and split it into chunks. This is because the WGS bgen
    files are very large, so we can only merge chunks of certain sizes.

    The default chunk size is 3Mb, but this can be changed to a custom size.

    The way that chunking works is we take the start position of a file, apply the 3Mb chunk size, and then make sure that
    the chunk does not overlap with any genes. If it does, we adjust the chunk to avoid the gene, but setting the chunk
    end between two genes. We also ensure that  the chunk end does not overlap with the next file.
    If it does, we adjust the chunk to end at the end position of the file.

    :param coordinate_path: file with the original coordinates
    :param gene_dict: dictionary object with gene positions
    :param chunk_size: size chunk for chunking
    :param output_path: path to the output file
    """

    # read in the original file
    coordinate_path = InputFileHandler(coordinate_path).get_file_handle()
    original_file = pd.read_csv(coordinate_path, sep='\t')
    # read in the gene dictionary
    gene_dict = InputFileHandler(gene_dict).get_file_handle()
    with open(gene_dict, 'r') as f:
        gene_dict = json.load(f)
    # do the splitting
    modified_file = split_coordinates_file(original_file, gene_dict, chunk_size)
    # output the chunks to separate files into a directory
    # create a new directory for the chunked files
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    output_chunked_file(chunked_file=modified_file, outout_dir=output_path)
    # print the output path
    print(f"Chunked file saved to: {output_path}")

def output_chunked_file(chunked_file: pd.DataFrame, outout_dir: Path) -> None:
    """
    Outputs the chunked file to separate CSV files, one for each unique chromosome.

    :param chunked_file: DataFrame containing the chunked coordinates.
    :param outout_dir: Directory to save the chunked files.
    :return: None
    """
    # select unique values in the 'chrom' column
    unique_chroms = chunked_file['chrom'].unique()
    # create a new DataFrame to store the chunked coordinates for each unique_chrom
    for chrom in unique_chroms:
        chrom_df = chunked_file[chunked_file['chrom'] == chrom].reset_index(drop=True)
        # name the new output file
        output_file = f"{chrom}.txt"
        # write the new DataFrame to a file with the header
        chrom_df.to_csv(outout_dir / output_file, sep='\t', index=False, header=True)


def split_coordinates_file(coordinates_file: pd.DataFrame, gene_dict: Union[dict, Path],
                           chunk_size: int = 3) -> pd.DataFrame:
    """
    Splits the coordinates file into chunks based on the specified chunk size and gene dictionary.

    :param coordinates_file: DataFrame containing the coordinates data
    :param gene_dict: Dictionary containing gene information with genomic locations
    :param chunk_size: Size of each chunk in base pairs. Default is 3 (interpreted as 3Mb)
    :return: DataFrame with the coordinates split into chunks
    """

    if chunk_size == 3:
        print(f"Creating {chunk_size}Mb chunks for bgens")
        chunk_size = chunk_size * 1000000
    else:
        print(f"Using a custom chunk size of {chunk_size}")

    coordinates_file = sort_coordinates_by_position(coordinates_file).reset_index(drop=True)

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

    chunks = []
    log_entries = []

    for chrom in sorted(coordinates_file['chrom'].unique(), key=lambda x: int(x.replace('chr', ''))):
        chrom_df = coordinates_file[coordinates_file['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)
        chrom_genes = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

        if chrom_df.empty or chrom_genes.empty:
            raise ValueError("Chromosome data or gene data is empty for chromosome: " + chrom)

        current_start = chrom_df['start'].min()
        chrom_max_end = chrom_df['end'].max()
        chunk_number = 1
        sub_chunk_index = 1
        splitting_logical_chunk = False
        original_chunk_start = current_start

        while current_start <= chrom_max_end:
            log_entry = {
                'chrom': chrom,
                'chunk_number': chunk_number,
                'start': current_start,
                'proposed_end_in_gene': None,
                'adjusted_to': 'not_applicable',
                'file_boundary_adjustment': 'not_applicable',
                'final_end_in_gene': None,
                'final_adjusted_to': 'not_applicable',
                'chunk_end': None,
                'file_limit_adjustment': 'not_applicable'
            }

            proposed_end = original_chunk_start + chunk_size
            gene_context = find_gene_context(chrom, proposed_end, gene_df)
            if gene_context['status'] == 'within':
                gene_name = gene_context['genes'][0]
                log_entry['proposed_end_in_gene'] = gene_name
                downstream = find_position_after_within_gene(chrom, gene_name, gene_df)
                if downstream is None:
                    proposed_end = chrom_max_end
                    log_entry['adjusted_to'] = 'end_of_chromosome'
                else:
                    proposed_end = downstream['between_position']
                    log_entry['adjusted_to'] = f"between {gene_name} and {downstream['next_gene']}"
            else:
                log_entry['proposed_end_in_gene'] = 'not_within_gene'

            overlapping_at_boundary = chrom_df[
                (chrom_df['start'] <= proposed_end) & (chrom_df['end'] > proposed_end)
            ]
            if not overlapping_at_boundary.empty:
                current_end = overlapping_at_boundary['end'].max()
                log_entry['file_boundary_adjustment'] = f"moved to {current_end} to avoid file split"
            else:
                current_end = proposed_end

            gene_context = find_gene_context(chrom, current_end, gene_df)
            if gene_context['status'] == 'within':
                gene_name = gene_context['genes'][0]
                log_entry['final_end_in_gene'] = gene_name
                downstream = find_position_after_within_gene(chrom, gene_name, gene_df)
                if downstream is None:
                    current_end = chrom_max_end
                    log_entry['final_adjusted_to'] = 'end_of_chromosome'
                else:
                    current_end = downstream['between_position']
                    log_entry['final_adjusted_to'] = f"between {gene_name} and {downstream['next_gene']}"
            else:
                log_entry['final_end_in_gene'] = 'not_within_gene'

            overlapping_rows = chrom_df[(chrom_df['end'] >= current_start) & (chrom_df['start'] <= current_end)]

            needs_split = False
            if len(overlapping_rows) > 800:
                needs_split = True
                splitting_logical_chunk = True
                limited_rows = overlapping_rows.iloc[:800]
                current_end = limited_rows['end'].max()
                log_entry['file_limit_adjustment'] = f"limited to 800 files, adjusted to {current_end}"

                gene_context = find_gene_context(chrom, current_end, gene_df)
                if gene_context['status'] == 'within':
                    gene_name = gene_context['genes'][0]
                    downstream = find_position_after_within_gene(chrom, gene_name, gene_df)
                    log_entry['final_end_in_gene'] = gene_name
                    if downstream is None:
                        current_end = chrom_max_end
                        log_entry['final_adjusted_to'] = 'end_of_chromosome (after file limit)'
                    else:
                        current_end = downstream['between_position']
                        log_entry['final_adjusted_to'] = f"between {gene_name} and {downstream['next_gene']} (after file limit)"
                else:
                    log_entry['final_end_in_gene'] = 'not_within_gene'
                    log_entry['final_adjusted_to'] = 'not_applicable'

                overlapping_rows = chrom_df[(chrom_df['end'] >= current_start) & (chrom_df['start'] <= current_end)]
            else:
                log_entry['file_limit_adjustment'] = 'not_applicable'

            if splitting_logical_chunk:
                chunk_label = f"{chrom}_chunk{chunk_number}_v{sub_chunk_index}"
            else:
                chunk_label = f"{chrom}_chunk{chunk_number}"

            for index, overlap_row in overlapping_rows.iterrows():
                chunks.append({
                    'chrom': chunk_label,
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

            log_entry['chunk_end'] = current_end
            log_entries.append(log_entry)

            next_files = chrom_df[chrom_df['start'] > current_end]
            if needs_split and not next_files.empty:
                current_start = next_files['start'].min()
                sub_chunk_index += 1
            else:
                if not next_files.empty:
                    current_start = next_files['start'].min()
                else:
                    break
                chunk_number += 1
                sub_chunk_index = 1
                splitting_logical_chunk = False
                original_chunk_start = current_start

    chunk_df = pd.DataFrame(chunks)
    log_df = pd.DataFrame(log_entries)
    log_file_path = Path("chunking_log_" + chrom + ".txt")
    log_df.to_csv(log_file_path, sep='\t', index=False)

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


if __name__ == "__main__":
    args = parse_arguments()
    run_splitter(
        coordinate_path=Path(args.coordinate_path),
        gene_dict=Path(args.gene_dict),
        chunk_size=args.chunk_size,
        output_path=Path(args.output_path)
    )
