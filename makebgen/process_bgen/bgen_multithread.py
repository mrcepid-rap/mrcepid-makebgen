import csv
from pathlib import Path
from typing import List, Dict, Any

import dxpy
from general_utilities.association_resources import check_gzipped
from general_utilities.import_utils.file_handlers.dnanexus_utilities import generate_linked_dx_file, \
    download_dxfile_by_name
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.subjob_utility import SubjobUtility, prep_current_image
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

from makebgen.process_bgen.process_bgen import make_final_bgen, make_bgen_from_vcf

LOGGER = MRCLogger().get_logger()


# loaded_module = check_subjob_decorator()
# if loaded_module:

def split_batch_files(batch: Path, max_rows: int) -> List[Path]:
    """
    Splits a TSV file into smaller files with a maximum number of rows.
    :param batch: Path to the TSV file.
    :param max_rows: Maximum number of rows per split file.
    :return: List of file paths for the split files.
    """
    split_files = []

    # Open the batch file and read its contents using DictReader
    with open(batch, 'r') as infile:
        reader = csv.DictReader(infile, delimiter='\t')

        # Extract the header (need to have one to make sure the output is valid)
        fieldnames = reader.fieldnames
        if fieldnames is None:
            raise ValueError(f"No header found in file: {batch}")

        # Initialize variables to keep track of rows and file index
        rows = []
        count = 0
        file_index = 0

        # Iterate through each row in the file
        for row in reader:
            # Append the current row
            rows.append(row)
            count += 1

            # If max_rows is reached, write to a new split file
            if count == max_rows:
                out_path = batch.parent / f"{batch.stem}_split_{file_index}.tsv"
                with open(out_path, 'w') as outfile:
                    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
                    # write the header
                    writer.writeheader()
                    # write the rows
                    writer.writerows(rows)
                # Append the output path to the list of split files
                split_files.append(out_path)

                # Reset the rows and count for the next split file
                rows = []
                count = 0
                file_index += 1

        # Write any remaining rows to a final split file
        if rows:
            out_path = batch.parent / f"{batch.stem}_split_{file_index}.tsv"
            with open(out_path, 'w', newline='') as outfile:
                writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
                writer.writeheader()
                writer.writerows(rows)
            split_files.append(out_path)

    # Return the list of split file paths
    return split_files


def process_batches(chunked_files: List, make_bcf: bool, output_prefix: str) -> Dict[int, List[Dict[str, Any]]]:
    """
    A function to process a batch of chunked files, converting BCF files to BGEN format and merging them.
    :param chunked_files: A list of chunked files to process.
    :param make_bcf: Should a concatenated BCF be made in addition to the BGEN?
    :param output_prefix: The prefix for the output files.
    :return: A list of dictionaries containing the output files for the batch.
    """
    # Log the start of processing for this batch

    subjob_launcher = SubjobUtility(log_update_time=600, incrementor=5, download_on_complete=False)

    for batch_index, batch_file in enumerate(chunked_files):

        LOGGER.info(f"Processing batch {batch_index} with {sum(1 for _ in open(batch_file)) - 1} chunked files...")

        # Split the batch file into smaller chunks if it exceeds the maximum number of rows
        # The maximum number of rows per chunk has to be 750 so that we can merge 1 BGEN
        split_batch = split_batch_files(batch_file, max_rows=750)

        # Launch subjobs for each chunk file in the batch
        for i, chunk_file in enumerate(split_batch):
            chunk_file_link = generate_linked_dx_file(chunk_file)

            subjob_launcher.launch_job(
                function=process_single_chunk,
                inputs={
                    'chunk_file': {"$dnanexus_link": chunk_file_link.get_id()},
                    'chunk_index': i,
                    'batch_index': batch_index,
                    'make_bcf': make_bcf,
                    'output_prefix': output_prefix
                },
                outputs=['bgen', 'bgen_index', 'sample', 'vep', 'vep_index', 'vcfprefix', 'start', 'batch_index'],
                instance_type='mem2_ssd1_v2_x16',
                name=f'makebgen_batch{batch_index}_chunk{i}',
            )

    subjob_launcher.submit_queue()

    bgen_prefixes: Dict[int, List[Dict[str, Any]]] = {}

    for subjob_output in subjob_launcher:
        if subjob_output['bgen'] is not None:
            # Each output is a dictionary with keys 'bgen', 'index', 'sample', 'vep', 'vep_idx'
            # Collect each output dictionary into a list
            if subjob_output['batch_index'] not in bgen_prefixes:
                bgen_prefixes[subjob_output['batch_index']] = []

            bgen_prefixes[subjob_output['batch_index']].append(subjob_output)

    return bgen_prefixes

def process_subjob_outputs(bgen_chunks: Dict[int, List[Dict[str, Any]]], make_bcf: bool,
                           output_prefix: str) -> dict[str, list[Any]]:

    output = {
        'bgen': [],
        'index': [],
        'sample': [],
        'vep': [],
        'vep_idx': []
    }

    for batch_index, batch_files in bgen_chunks.items():
        LOGGER.info(f"Processing batch {batch_index} with {len(batch_files)} chunk files...")

        file_handles = []
        for file_info in batch_files:
            for prefix in ['bgen', 'bgen_index', 'sample', 'vep', 'vep_index']:
                file_handles.append(download_dxfile_by_name(file_info[prefix]))

        merged = make_final_bgen(bgen_prefixes=batch_files, output_prefix=f"{output_prefix}_{batch_index}",
                                 make_bcf=make_bcf)

        # Set output as a list of dxlinks to the final files
        output['bgen'].append(merged['bgen']['file'])
        output['index'].append(merged['bgen']['index'])
        output['sample'].append(merged['bgen']['sample'])
        output['vep'].append(merged['vep']['file'])
        output['vep_idx'].append(merged['vep']['index'])

        # Now that the final files are gone, delete all chunk-level temporary files
        # Specify the suffix to delete
        # Delete all chunk-level temporary files except the final output files
        for file in file_handles:
            if file is not None:
                try:
                    file.unlink()
                    print(f"Deleted: {file}")
                except Exception as e:
                    print(f"Failed to delete {file}: {e}")

    return output


@dxpy.entry_point('process_single_chunk')
def process_single_chunk(chunk_file: dict, chunk_index: int,
                         batch_index: int, make_bcf: bool, output_prefix: str) -> dict:
    """
    Convert all VCFs in a chunk file to BGENs, then merge into one BGEN for that chunk.

    :param chunk_file: The file containing the chunk of coordinates to process.
    :param chunk_index: The index of the chunk being processed.
    :param batch_index: The index of the batch this chunk belongs to.
    :param make_bcf: Should a concatenated BCF be made in addition to the BGEN?
    :param output_prefix: The prefix for the output files.
    :return: A dictionary containing the output files for the chunk.
    """
    LOGGER.info(f"Starting processing of chunk {chunk_index} in batch {batch_index}")

    cmd_executor = prep_current_image([chunk_file])
    chunk_file = InputFileHandler(chunk_file, download_now=True).get_file_handle()

    with check_gzipped(chunk_file) as coord_file:
        coord_reader = list(csv.DictReader(coord_file, delimiter="\t"))

        thread_utility = ThreadUtility(
            incrementor=100,
            thread_factor=2,
            error_message='bcf to bgen thread failed'
        )

        previous_vep_id = None
        bgen_inputs = []

        for row in coord_reader:
            thread_utility.launch_job(
                make_bgen_from_vcf,
                vcf_id=row['output_bcf'],
                vep_id=row['output_vep'],
                previous_vep_id=previous_vep_id,
                start=row['start'],
                make_bcf=make_bcf,
                cmd_exec=cmd_executor
            )
            previous_vep_id = row['output_vep']

        for result in thread_utility:
            bgen_inputs.append(result)

    output_prefix = f"{output_prefix}_batch{batch_index}_chunk{chunk_index}"
    final_files = make_final_bgen(bgen_inputs, output_prefix, make_bcf)

    output = {
        'bgen': generate_linked_dx_file(final_files['bgen']['file']),
        'bgen_index': generate_linked_dx_file(final_files['bgen']['index']),
        'sample': generate_linked_dx_file(final_files['bgen']['sample']),
        'vep': generate_linked_dx_file(final_files['vep']['file']),
        'vep_index': generate_linked_dx_file(final_files['vep']['index']),
        'batch_index': batch_index
    }

    if final_files['bcf']['file'] is not None:
        output['bcf'] = dxpy.dxlink(generate_linked_dx_file(final_files['bcf']['file']))
        output['bcf_idx'] = dxpy.dxlink(generate_linked_dx_file(final_files['bcf']['index']))

    output['vcfprefix'] = output_prefix
    output['start'] = row['start']

    LOGGER.info(f"Finished chunk {chunk_index} in batch {batch_index}")

    return output
