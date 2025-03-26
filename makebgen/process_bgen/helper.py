"""
This script is designed to help with the processing of bgen files for WGS data.
"""

import pandas as pd


def split_coordinates_file(coordinates_file, gene_dict, chunk_size: int = 30):
    if chunk_size == 30:
        print("Creating 30Mb chunks for bgens")
        # convert chunk size to megabase
        chunk_size = chunk_size * 1000000
    else:
        print(f"Using a custom chunk size of {chunk_size}")

    # order by position
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

    for chrom in sorted(coordinates_file['chrom'].unique(), key=lambda x: int(x.replace('chr', ''))):
        chrom_df = coordinates_file[coordinates_file['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)
        chrom_genes = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

        if chrom_df.empty or chrom_genes.empty:
            # throw an error
            raise ValueError("Chromosome data or gene data is empty for chromosome: " + chrom)

        current_start = chrom_df['start'].min()
        chunk_number = 1
        chrom_max_end = chrom_df['end'].max()

        while current_start < chrom_max_end:
            position_to_check = current_start + chunk_size

            # check if within gene
            gene_context = find_gene_context(chrom, position_to_check, gene_df)

            if gene_context['status'] != 'within':
                # fallback: use chunk_size increment
                current_end = min(current_start + chunk_size, chrom_max_end)
            else:
                # ind downstream_pos
                gene_name = gene_context['genes'][0]
                downstream = find_position_after_within_gene(chrom, gene_name, gene_df)
                if downstream is None:
                    current_end = chrom_max_end
                else:
                    current_end = min(downstream['between_position'], chrom_max_end)

            # Find all rows in coordinates_file that overlap with this chunk
            overlapping_rows = chrom_df[(chrom_df['end'] >= current_start) & (chrom_df['start'] <= current_end)]

            # For each overlapping row, create a new entry with this chunk label
            for index, overlap_row in overlapping_rows.iterrows():
                chunks.append({
                    'chrom': f"{chrom}_chunk{chunk_number}",
                    'start': current_start,
                    'end': current_end,
                    'vcf_prefix': overlap_row['vcf_prefix'],
                    'output_bcf': overlap_row['output_bcf'],
                    'output_bcf_idx': overlap_row['output_bcf_idx'],
                    'output_vep': overlap_row['output_vep'],
                    'output_vep_idx': overlap_row['output_vep_idx']
                })

            # move to next chunk
            current_start = current_end + 1
            chunk_number += 1

    chunk_df = pd.DataFrame(chunks)

    chunk_df = chunk_df.drop_duplicates()

    return chunk_df


def find_chunk_for_position(chrom_df, position):
    match = chrom_df[(chrom_df['start'] <= position) & (chrom_df['end'] >= position)]
    if not match.empty:
        return match.iloc[0]
    else:
        return None


def find_position_after_within_gene(chrom, gene_name, gene_df):
    # Filter for this chromosome and sort by start
    genes_on_chrom = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

    # Find the index of the gene we're within
    gene_idx = genes_on_chrom[genes_on_chrom['gene'] == gene_name].index
    if gene_idx.empty:
        return None  # Gene not found

    gene_idx = gene_idx[0]

    # Make sure there's a next gene after this one
    if gene_idx + 1 >= len(genes_on_chrom):
        return None  # No next gene

    current_gene_end = genes_on_chrom.loc[gene_idx, 'end']
    next_gene_start = genes_on_chrom.loc[gene_idx + 1, 'start']

    # Calculate midpoint
    midpoint = (current_gene_end + next_gene_start) // 2

    return {
        'between_position': midpoint,
        'within_gene': genes_on_chrom.loc[gene_idx, 'gene'],
        'next_gene': genes_on_chrom.loc[gene_idx + 1, 'gene']
    }


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
