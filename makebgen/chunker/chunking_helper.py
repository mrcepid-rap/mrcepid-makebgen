"""
This is a helper script to chunk the BGEN files into chunks of 3Mb or less than 750 files, so that they can
be concatenated into a single BGEN file. Note that a key part of this script is that the ends of each chunk
cannot be inside a gene, instead we must find a chunk end that is safe, i.e. not overlapping with any gene.
"""

import json
import dxpy
import pandas as pd
from pathlib import Path
from intervaltree import IntervalTree
from typing import Union, Tuple, List, Dict, Any

from general_utilities.import_utils.file_handlers.dnanexus_utilities import generate_linked_dx_file, LOGGER


def parse_gene_dict(gene_dict_path: Union[str, Path]) -> pd.DataFrame:
    """
    Parse a gene dictionary JSON file into a DataFrame.

    :param gene_dict_path: Path to the gene dictionary JSON file.
    :return: DataFrame containing gene information.
    """
    # load the gene dictionary from the JSON file
    with open(gene_dict_path, 'r') as genes:
        raw = json.load(genes)

    records = []

    # Extract relevant information from the raw gene dictionary
    for gene, data in raw.items():
        loc = data['genomic_location']
        records.append({
            'gene': gene,
            'chrom': f"chr{loc['chromosome']}",
            'start': loc['start'],
            'end': loc['end']
        })

    # Create a DataFrame from the records
    return pd.DataFrame(records)


def build_interval_tree(df: pd.DataFrame, chrom: str, start_col: str, end_col: str,
                        label_col: str = None) -> IntervalTree:
    """
    Builds an IntervalTree from a DataFrame for a specific chromosome.

    :param df: DataFrame containing genomic intervals.
    :param chrom: Chromosome to filter the DataFrame on.
    :param start_col: Column name for the start positions.
    :param end_col: Column name for the end positions.
    :param label_col: Optional column name for labels to associate with intervals.
    :return: IntervalTree containing intervals for the specified chromosome.
    """

    # create an IntervalTree for the specified chromosome
    tree = IntervalTree()

    # For the specified chromosome, iterate through each row of the DataFrame
    # and insert an interval into the interval tree using the start and end positions.
    # If a label column is specified, associate that label with the interval.
    for index, row in df[df['chrom'] == chrom].iterrows():
        label = row[label_col] if label_col else None
        tree[row[start_col]:row[end_col]] = label
    return tree


def is_position_within_gene(tree: IntervalTree, pos: int) -> Union[str, None]:
    """
    Checks if a given position is within any gene interval in the IntervalTree.

    :param tree: IntervalTree containing gene intervals.
    :param pos: Position to check.
    :return: The label of the first gene that contains the position, or None if no gene contains it.
    """
    hits = list(tree[pos])
    return hits[0].data if hits else None


def get_safe_chunk_ends(chrom_df, gene_intervals):
    """
    Filters out positions in chrom_df['end'] that overlap with any interval in gene_intervals.

    :param chrom_df: DataFrame containing chromosome data with 'end' positions.
    :param gene_intervals: IntervalTree containing gene intervals.
    :return: List of 'end' positions that do not overlap with any gene intervals.
    """
    safe_ends = []
    for pos in chrom_df['end']:
        # Check if the position overlaps with any interval in the tree
        if not gene_intervals.overlaps(pos):
            safe_ends.append(pos)
    return safe_ends


def find_next_safe_end(safe_chunk_ends: List[int], proposed_end: int) -> int:
    """
    Finds the next safe end position that is less than the proposed end.
    To be used when the proposed end is within a gene.

    :param safe_chunk_ends: List of safe chunk end positions.
    :param proposed_end: The proposed end position to check against.
    :return: The maximum safe end position that is less than the proposed end.
    """
    next_safe_ends = [end for end in safe_chunk_ends if end < proposed_end]
    if not next_safe_ends:
        raise RuntimeError(f"No safe chunk end found after {proposed_end}")
    return max(next_safe_ends)


def chunk_chromosome(chrom_df: pd.DataFrame, gene_df: pd.DataFrame, chrom: str, chunk_size_bp: int) -> Tuple[
    pd.DataFrame, List[Dict]]:
    """
    Splits a chromosome DataFrame into chunks of specified size, ensuring that chunks do not overlap with genes.

    :param chrom_df: DataFrame containing chromosome data with 'start' and 'end' positions.
    :param gene_df: DataFrame containing gene information.
    :param chrom: Chromosome identifier to filter the DataFrame.
    :param chunk_size_bp: Size of each chunk in base pairs.
    :return: Tuple containing a DataFrame of chunks and a log entry for each chunk.
    """

    # build an interval tree for the genes on this chromosome
    gene_tree = build_interval_tree(gene_df, chrom, 'start', 'end', 'gene')
    # sort the chromosome DataFrame by 'start' position
    chrom_df = chrom_df.sort_values(by='start').reset_index(drop=True)

    chunks = []
    log_entries = []

    # set the initial start position and maximum end position for the chromosome
    current_start = chrom_df['start'].min()
    chrom_max = chrom_df['end'].max()
    # get the safe chunk ends that do not overlap with genes
    safe_chunk_ends = get_safe_chunk_ends(chrom_df, gene_tree)
    chunk_number = 0

    while current_start <= chrom_max:
        # Calculate the ideal end position for the current chunk
        ideal_end = current_start + chunk_size_bp
        # note how many files are remaining in the DataFrame
        files_remaining = chrom_df[chrom_df['start'] >= current_start]
        # if there are no files remaining, break the loop
        if files_remaining.empty:
            break

        # Filter to candidate ends <= ideal_end and safe
        safe_within_limit = [end for end in safe_chunk_ends if current_start < end <= ideal_end]
        # if no safe end is found within the limit, raise an error
        if not safe_within_limit:
            raise RuntimeError(f"No safe chunk end found near {ideal_end} on {chrom}")

        # Set the proposed end to the maximum of the safe ends within the limit
        proposed_end = max(safe_within_limit)

        # add a look-ahead check to see if we are near the end of the chromosome
        next_df = chrom_df[chrom_df['start'] > proposed_end]
        if not next_df.empty and len(next_df) < 300:
            # absorb the remaining files into this chunk
            proposed_end = chrom_df['end'].max()

        # see if the proposed end is within a gene (it absolutely should not be)
        within_gene = is_position_within_gene(gene_tree, proposed_end)
        # add a check to ensure the proposed end is not within a gene (none of them should be)
        if within_gene:
            # this happens for chr6, which ends within a non-coding gene
            LOGGER.warning(
                f"File limit reached but proposed_end {proposed_end} is within a gene, something went wrong. "
                f"If you are on chromosome 6, this is likely due to the chromosome ending within a non-coding gene. "
                f"Also worth checking - did you set the chunk size too small? "
                f"Please look at the chunking logs and the concatenated output to ensure everything is correct.")
        # calculate the chunk rows that fall within the proposed end
        chunk_rows = files_remaining[files_remaining['start'] <= proposed_end]

        # Log info
        log_entry = {
            'chrom': chrom,
            'chunk_number': f"{chrom}_chunk{chunk_number}",
            'chunk_start': current_start,
            'proposed_end': ideal_end,
            'adjusted_for_file': proposed_end,
            'adjusted_for_gene': 'safe',
            'final_chunk_end': proposed_end,
            'not_within_gene': within_gene is None
        }
        # Add gene context (for logging purposes)
        nearby_genes = gene_df[gene_df['chrom'] == chrom].sort_values(by='start').reset_index(drop=True)
        upstream_gene = nearby_genes[nearby_genes['end'] < proposed_end]
        downstream_gene = nearby_genes[nearby_genes['start'] > proposed_end]
        gene_before = upstream_gene.iloc[-1]['gene'] if not upstream_gene.empty else 'None'
        gene_after = downstream_gene.iloc[0]['gene'] if not downstream_gene.empty else 'None'
        log_entry['between_genes'] = f"{gene_before} and {gene_after}"

        # Create the chunk entry
        for _, row in chunk_rows.iterrows():
            chunks.append({
                'chrom': log_entry['chunk_number'],
                'start': row['start'],
                'end': row['end'],
                'chunk_start': current_start,
                'chunk_end': proposed_end,
                'vcf_prefix': row['vcf_prefix'],
                'output_bcf': row['output_bcf'],
                'output_bcf_idx': row['output_bcf_idx'],
                'output_vep': row['output_vep'],
                'output_vep_idx': row['output_vep_idx']
            })

        # Add the log entry to the list
        log_entries.append(log_entry)
        # Update the current start position to the next start after the proposed end
        next_df = chrom_df[chrom_df['start'] > proposed_end]
        # if there are no more rows, break the loop
        if next_df.empty:
            break
        # ensure the next start is greater than the proposed end
        current_start = next_df['start'].min()
        # increment the chunk number
        chunk_number += 1

    # convert the chunks list to a DataFrame
    chunk_df = pd.DataFrame(chunks)

    # return the chunk DataFrame and the log entries
    return chunk_df, log_entries


def split_coordinates_file(coordinates_file: pd.DataFrame, gene_df: pd.DataFrame, chunk_size: int):
    """
    Splits the coordinates file into chunks based on the specified chunk size,
    ensuring that chunks do not overlap with genes.

    :param coordinates_file: Dataframe with the genetic file coordinates.
    :param gene_df: Dataframe with gene information.
    :param chunk_size: Size of each chunk in megabases (default: 3Mb).
    :return: Tuple containing a DataFrame of chunks and a list of log entries.
    """
    # Convert chunk size from megabases to base pairs
    chunk_size_bp = int(chunk_size * 1_000_000)

    # Initialize lists to hold all chunks and logs
    all_chunks = []
    all_logs = []

    # Iterate through each unique chromosome in the coordinates file
    for chrom in coordinates_file['chrom'].unique():
        # Filter the coordinates file for the current chromosome
        chrom_df = coordinates_file[coordinates_file['chrom'] == chrom]
        # run the chunking function for the current chromosome
        chunk_df, log_entries = chunk_chromosome(chrom_df, gene_df, chrom, chunk_size_bp)
        # Add the chunk DataFrame and log entries to the output
        all_chunks.append(chunk_df)
        all_logs.extend(log_entries)

    # Concatenate all chunks into a single DataFrame and reset the index
    # return the logs
    return pd.concat(all_chunks).reset_index(drop=True), all_logs


def get_chunk_number(path: Path) -> int:
    """
    Extracts the chunk number from the filename.
    :param path: Path object representing the file.
    :return: int: The chunk number extracted from the filename.
    """
    return int(path.name.split('chunk')[-1])


def chunking_helper(gene_dict: Path, coordinate_path: Path, chunk_size: int) -> Tuple[List[Any], List[Dict]]:
    """
    Main function to create BGEN chunk coordinates

    :param gene_dict: Path to the JSON file containing gene dictionary.
    :param coordinate_path: Path to the input coordinate file.
    :param chunk_size: Size of each chunk in megabases.
    :return: List of file paths for the chunked output files.
    """

    # Read the coordinate file and prase gene dictionary
    coord_df = pd.read_csv(coordinate_path, sep='\t')
    gene_df = parse_gene_dict(gene_dict)

    # remove duplicates (if they exist)
    coord_df = coord_df.drop_duplicates(subset=["chrom", "start", "end", "vcf_prefix"], keep="first")

    # chrom
    unique_chrom = coord_df["chrom"].unique()[0]
    output_path = Path(f"chunked_output_{unique_chrom}")
    output_path.mkdir(parents=True, exist_ok=True)

    # Split the coordinates file into chunks
    chunked_df, log_entries = split_coordinates_file(coord_df, gene_df, chunk_size)

    # Save each chunk to a separate file
    for chunk_label in chunked_df['chrom'].unique():
        chunk_df = chunked_df[chunked_df['chrom'] == chunk_label]
        chunk_df.to_csv(output_path / f"{chunk_label}.txt", sep='\t', index=False)

    # Save the log entries to a text file
    log_entires_path = f"chunking_log_{unique_chrom}.txt"
    pd.DataFrame(log_entries).to_csv(log_entires_path, sep='\t', index=False)
    # also save all chunks combined into a single file
    chunks_combined_path = f"all_chunks_combined_{unique_chrom}.txt"
    chunked_df.to_csv(chunks_combined_path, sep='\t', index=False)

    output_files = sorted(Path(output_path).iterdir(), key=lambda p: int(p.name.split('chunk')[-1].replace('.txt', '')))

    # # we should upload the helper files for reference & safe-keeping
    log_files = [dxpy.dxlink(generate_linked_dx_file(file=log_entires_path, delete_on_upload=False)),
                 dxpy.dxlink(
                     generate_linked_dx_file(file=chunks_combined_path, delete_on_upload=False))]

    return output_files, log_files


if __name__ == "__main__":

    output_files, log_files = chunking_helper(
        gene_dict=Path("/Users/alish.palmos/PycharmProjects/mrcepid-makebgen/test/test_data/final_dict.json"),
        coordinate_path=Path("/Users/alish.palmos/PycharmProjects/mrcepid-makebgen/test/test_data/chr6_concatenated_output.txt"),
        chunk_size=9
    )
