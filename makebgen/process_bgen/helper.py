"""
This script is designed to help with the processing of bgen files for WGS data.
"""

import pandas as pd


def split_coordinates_file(coordinates_file: pd.DataFrame, gene_dict: dict, chunk_size: int = 30) -> pd.DataFrame:
    """
    Splits the coordinates file into chunks based on the specified chunk size and gene dictionary.

    :param coordinates_file: DataFrame containing the coordinates data.
    :param gene_dict: Dictionary containing gene information with genomic locations.
    :param chunk_size: Size of each chunk in base pairs. Default is 30 (interpreted as 30Mb).
    :return: DataFrame with the coordinates split into chunks.
    """

    # first, if we are using a 30Mb chunk then we need to convert it to Mb
    if chunk_size == 30:
        print("Creating 30Mb chunks for bgens")
        # convert chunk size to megabase
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

    # make a new dataframe
    chunk_df = pd.DataFrame(chunks)

    return chunk_df


def find_chunk_for_position(chrom_df, position):
    match = chrom_df[(chrom_df['start'] <= position) & (chrom_df['end'] >= position)]
    if not match.empty:
        return match.iloc[0]
    else:
        return None


def find_position_after_within_gene(chrom, gene_name, gene_df):
    # Get all genes on this chromosome sorted by start
    genes_on_chrom = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

    # Find the index of the current gene
    gene_idx = genes_on_chrom[genes_on_chrom['gene'] == gene_name].index
    if gene_idx.empty:
        return None  # Gene not found

    gene_idx = gene_idx[0]
    current_end = genes_on_chrom.loc[gene_idx, 'end']

    # Look ahead to find the first downstream gene that doesn't overlap
    for i in range(gene_idx + 1, len(genes_on_chrom)):
        next_start = genes_on_chrom.loc[i, 'start']
        next_gene = genes_on_chrom.loc[i, 'gene']

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
    return None



def find_gene_context(chrom, position, gene_df):
    # Filter to just the chromosome of interest
    genes_on_chrom = gene_df[gene_df['chrom'] == chrom]

    # Check if position falls within any gene
    within = genes_on_chrom[(genes_on_chrom['start'] <= position) & (genes_on_chrom['end'] >= position)]
    if not within.empty:
        return {'status': 'within', 'genes': within['gene'].tolist()}

    # If not within, find closest upstream and downstream genes
    upstream = genes_on_chrom[genes_on_chrom['end'] < position]
    downstream = genes_on_chrom[genes_on_chrom['start'] > position]

    closest_up = upstream.iloc[upstream['end'].sub(position).abs().argmin()] if not upstream.empty else None
    closest_down = downstream.iloc[downstream['start'].sub(position).abs().argmin()] if not downstream.empty else None

    return {
        'status': 'between',
        'closest_upstream': closest_up['gene'] if closest_up is not None else None,
        'closest_downstream': closest_down['gene'] if closest_down is not None else None
    }


def sort_coordinates_by_position(df):
    df['chrom_num'] = df['chrom'].replace('chr', '')
    return df.sort_values(by=['chrom_num', 'start']).drop(columns='chrom_num')
