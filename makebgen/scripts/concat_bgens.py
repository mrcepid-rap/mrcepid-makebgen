import gzip
import re
from pathlib import Path
from typing import Iterator, List, Dict, Union, Iterable, Any

from general_utilities.association_resources import bgzip_and_tabix
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import CommandExecutor, build_default_command_executor

CMD_EXEC = build_default_command_executor()


def process_chunk_group(batch_index: int, chunk: List[Dict[str, Union[str, Path, dict]]],
                        cmd_exec: CommandExecutor = CMD_EXEC) -> Dict[str, Path]:
    """
    Process a group of chunks:
    - Download files
    - Concatenate BGEN, SAMPLE, VEP
    - Index BGEN and VEP
    - Clean up originals

    :param batch_index: Index used in naming output files.
    :param chunk: List of dictionaries, each containing keys: 'bgen', 'sample', 'vep', and 'chrom'.
    :param cmd_exec: CommandExecutor instance to run shell commands. Defaults to CMD_EXEC.
    :return: Dictionary of output file paths with keys: 'bgen', 'bgen_index', 'sample', 'vep', 'vep_index'.
    """

    files = {
        'bgen': [],
        'sample': [],
        'vep': []
    }

    # download files
    for i, row in enumerate(chunk):
        bgen = InputFileHandler(row['bgen']).get_file_handle()
        sample = InputFileHandler(row['sample']).get_file_handle()
        vep = InputFileHandler(row['vep']).get_file_handle()

        files['bgen'].append(bgen)
        files['sample'].append(sample)
        files['vep'].append(vep)

    # chrom label is the same for all rows in the chunk
    # (assuming they are all from the same chromosome)
    bgen_files = [row['bgen'] for row in chunk]
    chrom_label = check_all_same_chrom(bgen_files)

    # create output files
    out_bgen = Path(f"{chrom_label}_mergedchunk_{batch_index}.bgen")
    out_sample = Path(f"{chrom_label}_mergedchunk_{batch_index}.sample")
    out_vep = Path(f"{chrom_label}_mergedchunk_{batch_index}.vep.tsv")

    # Assume you're inside the correct working directory in Docker
    input_bgens = ' '.join(f'-g /test/{bgen.name}' for bgen in files['bgen'])
    cmd = f'cat-bgen {input_bgens} -og /test/{out_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # index the bgen
    cmd = f'bgenix -index -g /test/{out_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # ensure the sample files are identical
    with files['sample'][0].open("r") as in_f:
        with out_sample.open("w") as out_f:
            for line in in_f:
                out_f.write(line)

    # concatenate and sort VEP
    with out_vep.open("w") as out_f:
        for i, vep in enumerate(files['vep']):
            with gzip.open(vep, "rt") as in_f:
                for j, line in enumerate(in_f):
                    if i == 0 and j == 0:
                        out_f.write(line)  # write header once
                    elif j > 0:
                        out_f.write(line)

    # sort VEP
    final_vep, final_vep_idx = bgzip_and_tabix(out_vep, comment_char='C', end_row=2)

    # clean up original files
    for f in files['bgen'] + files['sample'] + files['vep']:
        if f.exists():
            f.unlink()
    out_vep.unlink()

    return {
        'bgen': out_bgen,
        'bgen_index': Path(f"{out_bgen}.bgi"),
        'sample': out_sample,
        'vep': final_vep,
        'vep_index': final_vep_idx
    }


def chunk_dict_reader(reader: Union[Iterable[Dict[str, Any]], Dict[str, Any]], chunk_size: int) -> Iterator[
    List[Dict[str, Any]]]:
    """
    Yield batches of rows from an iterable of dictionaries or a single dictionary in chunks of `chunk_size`.

    :param reader: An iterable of dicts (e.g. csv.DictReader, list of dicts), or a plain dictionary
    :param chunk_size: Number of rows per chunk
    """
    if isinstance(reader, dict):
        # Convert dictionary to list of {"key": ..., "value": ...} format
        reader = [{"key": k, "value": v} for k, v in reader.items()]

    chunk = []
    for row in reader:
        chunk.append(row)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


def extract_chrom_from_filename(filename: str) -> str:
    """
    Extract the chromosome from a filename formatted as "<chromosome>_<rest_of_filename>".

    :param filename: Filename from which to extract the chromosome.
    :return: The extracted chromosome part of the filename.
    """
    return Path(filename).stem.split('_')[0]


def check_all_same_chrom(file_list: List[str]) -> str:
    """
    Check if all files in the list belong to the same chromosome based on their filenames.
    :param file_list: List of filenames to check.
    :return: The chromosome name if all files are from the same chromosome.
    """
    chroms = {extract_chrom_from_filename(f) for f in file_list}
    if len(chroms) != 1:
        raise ValueError(f"Inconsistent chromosomes in chunk: {chroms}")
    return chroms.pop()

