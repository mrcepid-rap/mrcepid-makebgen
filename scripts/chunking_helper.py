"""
This script is designed to help with the processing of bgen files for WGS data.
It splits large genomic coordinate files into chunks, while ensuring that chunks do not break across genes or file boundaries,
and that each chunk contains no more than a specified number of files (default: 750).
"""

import argparse
import json
from pathlib import Path
from typing import Union, Tuple

import pandas as pd
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler


def parse_arguments():
    """Parse command-line arguments for input file paths, chunk size, and output directory."""
    parser = argparse.ArgumentParser(description="Process bgen files for WGS data.")
    parser.add_argument("--coordinate_path", type=str, required=True, help="Path to the coordinate file.")
    parser.add_argument("--gene_dict", type=str, required=True, help="Path to the gene dictionary file.")
    parser.add_argument("--chunk_size", type=int, default=3, required=False, help="Chunk size in Mb (default: 3Mb).")
    parser.add_argument("--output_path", type=str, default="chunked_files", required=False,
                        help="Path to the output directory (default: 'chunked_files').")
    return parser.parse_args()


def output_chunked_file(chunked_file: pd.DataFrame, output_dir: Path) -> None:
    """
    Write each chromosome's chunked coordinates to a separate file in the output directory.

    :param chunked_file: DataFrame containing chunked coordinates
    :param output_dir: Directory where the chunked files will be saved.
    :return: None
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for chrom in chunked_file['chrom'].unique():
        chrom_df = chunked_file[chunked_file['chrom'] == chrom].reset_index(drop=True)
        output_file = f"{chrom}.txt"
        chrom_df.to_csv(output_dir / output_file, sep='\t', index=False, header=True)


def split_coordinates_file(coordinates_file: pd.DataFrame, gene_dict: Union[dict, Path],
                           chunk_size: int = 3) -> pd.DataFrame:
    """
    Orchestrates the splitting of coordinates into chunks based on chunk size and gene boundaries.
    Returns a DataFrame of all chunked records.
    """
    print(f"Creating {chunk_size}Mb chunks for bgens")
    chunk_size = chunk_size * 1_000_000
    gene_df = parse_gene_dict_to_df(gene_dict)

    chunks, logs = process_chromosome_chunks(coordinates_file, gene_df, chunk_size)

    chunk_df = pd.DataFrame(chunks)
    chunk_df = sort_coordinates_by_position(chunk_df).reset_index(drop=True)
    log_df = pd.DataFrame(logs)
    log_file_path = Path(f"chunking_log_{chunks[-1]['chrom'].split('_')[0]}.txt")
    log_df.to_csv(log_file_path, sep='\t', index=False)

    return chunk_df


def process_chromosome_chunks(coordinates_file: pd.DataFrame, gene_df: pd.DataFrame, chunk_size: int) -> Tuple[list, list]:
    """
    Iterates over each chromosome and splits its coordinates into chunks.

    :param coordinates_file: DataFrame containing genomic coordinates
    :param gene_df: DataFrame containing gene information.
    :param chunk_size: Size of each chunk in base pairs.
    :return: A tuple containing a list of dictionaries representing chunks and a list of logs for each chunk.
    """
    chunks = []
    logs = []

    for chrom in coordinates_file['chrom'].unique():
        chrom_df = coordinates_file[coordinates_file['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)
        chrom_genes = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)

        if chrom_df.empty or chrom_genes.empty:
            raise ValueError(f"Chromosome data or gene data is empty for {chrom}")

        chrom_chunks, chrom_logs = generate_chunks_for_chromosome(chrom_df, chrom_genes, chrom, chunk_size)
        chunks.extend(chrom_chunks)
        logs.extend(chrom_logs)

    return chunks, logs


def generate_chunks_for_chromosome(chrom_df: pd.DataFrame, chrom_genes: pd.DataFrame, chrom: str, chunk_size: int) -> Tuple[list, list]:
    """
    Splits a single chromosome's coordinate data into valid chunks.

    :param chrom_df: DataFrame containing chromosome coordinates.
    :param chrom_genes: DataFrame containing genes for the current chromosome.
    :param chrom: Chromosome identifier (e.g., 'chr1').
    :param chunk_size: Size of each chunk in base pairs.
    :return: A list of dictionaries representing chunks and a list of logs for each chunk.
    """
    chunks = []
    logs = []

    current_start = chrom_df['start'].min()
    chrom_max_end = chrom_df['end'].max()
    chunk_number = 1
    sub_chunk_index = 1

    while current_start <= chrom_max_end:
        original_chunk_start = current_start
        chunk_info, overlapping_rows = adjust_chunk_end(
            chrom_df, chrom_genes, chrom, current_start, original_chunk_start + chunk_size, chunk_number
        )

        # Create label: v1 if we hit file limit, otherwise just chunk number
        chunk_label = (
            f"{chrom}_chunk{chunk_number}_v{sub_chunk_index}"
            if chunk_info['file_limit_adjustment'] != 'not_applicable'
            else f"{chrom}_chunk{chunk_number}"
        )

        for row in overlapping_rows.itertuples(index=False):
            chunks.append({
                'chrom': chunk_label,
                'start': row.start,
                'end': row.end,
                'chunk_start': chunk_info['chunk_start'],
                'chunk_end': chunk_info['chunk_end'],
                'vcf_prefix': row.vcf_prefix,
                'output_bcf': row.output_bcf,
                'output_bcf_idx': row.output_bcf_idx,
                'output_vep': row.output_vep,
                'output_vep_idx': row.output_vep_idx
            })

        logs.append(chunk_info)

        # Prepare for next chunk
        next_files = chrom_df[chrom_df['start'] > chunk_info['chunk_end']]
        if chunk_info['file_limit_adjustment'] != 'not_applicable' and not next_files.empty:
            current_start = next_files['start'].min()
            sub_chunk_index += 1
        else:
            if not next_files.empty:
                current_start = next_files['start'].min()
            else:
                break
            chunk_number += 1
            sub_chunk_index = 1

    return chunks, logs


def adjust_chunk_end(chrom_df: pd.DataFrame, chrom_genes: pd.DataFrame, chrom: str, current_start: int,
                     proposed_end: int, chunk_number: int) -> Tuple[dict, pd.DataFrame]:
    """
    Adjust the chunk end to avoid falling within a gene and avoid splitting files.
    """
    log_entry = {
        'chrom': chrom,
        'chunk_number': f"{chrom}_chunk{chunk_number}",
        'chunk_start': current_start,
        'file_boundary_adjustment': 'not_applicable',
        'final_adjusted_to': 'not_applicable',
        'final_end_in_gene': None,
        'chunk_end': None,
        'file_limit_adjustment': 'not_applicable'
    }

    max_end = chrom_df['end'].max()
    current_end = min(proposed_end, max_end)

    # Step 1: Make sure we don't end within a gene
    while True:
        gene_context = find_gene_context(current_end, chrom_genes)
        if gene_context['status'] == 'within':
            gene_name = gene_context['genes'][0]
            downstream = find_position_after_within_gene(chrom, gene_name, chrom_genes)
            log_entry['final_end_in_gene'] = gene_name

            if downstream is None:
                current_end = max_end
                break
            else:
                current_end = downstream['between_position']
        else:
            break

    # Step 2: Adjust to file boundary if needed
    overlapping_at_boundary = chrom_df[
        (chrom_df['start'] <= current_end) & (chrom_df['end'] > current_end)
    ]
    if not overlapping_at_boundary.empty:
        current_end = overlapping_at_boundary['end'].max()
        log_entry['file_boundary_adjustment'] = f"moved to {current_end} to avoid file split"

        # Re-check: if this new boundary falls within a gene, restart step 1
        gene_context = find_gene_context(current_end, chrom_genes)
        if gene_context['status'] == 'within':
            gene_name = gene_context['genes'][0]
            downstream = find_position_after_within_gene(chrom, gene_name, chrom_genes)
            log_entry['final_end_in_gene'] = gene_name
            if downstream is None:
                current_end = max_end
            else:
                current_end = downstream['between_position']

    # Step 3: Final context note
    final_context = find_gene_context(current_end, chrom_genes)
    if final_context['status'] == 'between':
        up = final_context.get('closest_upstream')
        down = final_context.get('closest_downstream')
        if up and down:
            log_entry['final_adjusted_to'] = f"between {up} and {down}"
        elif up:
            log_entry['final_adjusted_to'] = f"after {up}"
        elif down:
            log_entry['final_adjusted_to'] = f"before {down}"
        else:
            log_entry['final_adjusted_to'] = 'no nearby genes found'

    # Step 4: Extract overlapping rows and enforce file count
    overlapping_rows = chrom_df[
        (chrom_df['end'] >= current_start) & (chrom_df['start'] <= current_end)
    ]
    if len(overlapping_rows) > 750:
        limited_rows = overlapping_rows.iloc[:750]
        current_end = limited_rows['end'].max()
        overlapping_rows = limited_rows
        log_entry['file_limit_adjustment'] = f"limited to 750 files, adjusted to {current_end}"

    log_entry['chunk_end'] = current_end
    return log_entry, overlapping_rows



def parse_gene_dict_to_df(gene_dict: dict) -> pd.DataFrame:
    """
    Converts a gene dictionary into a sorted DataFrame with columns ['gene', 'chrom', 'start', 'end'].

    :param gene_dict: Dictionary with gene names as keys and genomic locations as values.
    :return: DataFrame with sorted gene coordinates.
    """
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
    return sort_coordinates_by_position(gene_df).reset_index(drop=True)


def find_position_after_within_gene(chrom: str, gene_name: str, gene_df: pd.DataFrame) -> Union[dict, None]:
    """
    Finds a midpoint between the current gene and the next non-overlapping gene on the same chromosome.
    Returns None if no downstream gene exists.

    :param chrom: Chromosome identifier (e.g., 'chr1')
    :param gene_name: Name of the current gene
    :param gene_df: DataFrame containing gene information with columns ['gene', 'chrom', 'start', 'end']
    :return: Dictionary with 'between_position', 'within_gene', and 'next_gene' keys, or None if no next gene.
    """
    genes_on_chrom = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)
    gene_idx = genes_on_chrom[genes_on_chrom['gene'] == gene_name].index
    if gene_idx.empty:
        raise ValueError("Gene not found: " + str(gene_name))

    gene_idx = gene_idx[0]
    current_end = genes_on_chrom.loc[gene_idx, 'end']

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
        current_end = max(current_end, genes_on_chrom.loc[i, 'end'])

    return None


def find_gene_context(position: int, genes_on_chrom: pd.DataFrame) -> dict:
    """
    Determines whether the given position falls within a gene or between genes.
    Returns a dictionary indicating status and gene context.

    :param position: Position to check
    :param genes_on_chrom: DataFrame of genes on the chromosome
    """
    within = genes_on_chrom[(genes_on_chrom['start'] <= position) & (genes_on_chrom['end'] >= position)]
    if not within.empty:
        return {'status': 'within', 'genes': within['gene'].tolist()}

    upstream = genes_on_chrom[genes_on_chrom['end'] < position]
    downstream = genes_on_chrom[genes_on_chrom['start'] > position]

    closest_up = upstream.iloc[-1] if not upstream.empty else None
    closest_down = downstream.iloc[0] if not downstream.empty else None

    return {
        'status': 'between',
        'closest_upstream': closest_up['gene'] if closest_up is not None else None,
        'closest_downstream': closest_down['gene'] if closest_down is not None else None
    }


def sort_coordinates_by_position(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sorts a coordinate dataframe by chromosome number (after stripping 'chr') and start position.

    :param df: coordinate dataframe
    :return: sorted dataframe with an additional 'chrom_num' column for sorting
    """
    df['chrom_num'] = df['chrom'].str.replace('chr', '')
    return df.sort_values(by=['chrom_num', 'start']).drop(columns='chrom_num')


if __name__ == "__main__":
    # Parse CLI args
    args = parse_arguments()

    # Load coordinate file
    coordinate_path = InputFileHandler(Path(args.coordinate_path)).get_file_handle()
    original_file = pd.read_csv(coordinate_path, sep='\t')
    original_file = sort_coordinates_by_position(original_file).reset_index(drop=True)

    # Load gene dictionary
    gene_dict_path = InputFileHandler(Path(args.gene_dict)).get_file_handle()
    with open(gene_dict_path, 'r') as f:
        gene_dict = json.load(f)

    # Create output directory
    output_path = Path(args.output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    # Perform chunking
    chunked_df = split_coordinates_file(original_file, gene_dict, args.chunk_size)

    # Save output
    output_chunked_file(chunked_file=chunked_df, output_dir=output_path)
    print(f"Chunked file saved to: {output_path}")
