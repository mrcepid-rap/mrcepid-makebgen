import csv
from pathlib import Path
from typing import List, Dict, Any

import dxpy
from general_utilities.association_resources import check_gzipped
from general_utilities.import_utils.file_handlers.dnanexus_utilities import generate_linked_dx_file
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.subjob_utility import SubjobUtility, prep_current_image
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger

from makebgen.process_bgen.process_bgen import make_final_bgen, make_bgen_from_vcf

LOGGER = MRCLogger().get_logger()


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
    A function that handles processing of batches of chunked files, by converting BCF files to BGEN format and
    merging them. This method uses a SobjobUtility to launch subjobs and then collect them on completion.

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
                           output_prefix: str) -> Dict[str, List[Any]]:

    """Collect batches of merged files from :func:`process_batches` and merge them into final BGEN files. This is
    necessary as the `cat-bgen` command can only handle ≤750 files at a time.

    Note that this method is not parallelized. This is by design, as the process of merging BGEN files is storage
    intensive, with a typical batch requiring around 1Tb of disk space to complete when a 9Mb ideal_chunk_size is
    used. This means that parallelizing this step would require a lot of disk space, which is not available on the
    DNANexus platform without high cost.

    :param bgen_chunks: A dictionary where keys are batch indices and values are lists of dictionaries containing file
        information for all files in a batch.
    :param make_bcf: Should a concatenated BCF be made in addition to the BGEN?
    :param output_prefix: The prefix for the output files. The final BGEN files will be named like
        :param: output_prefix_batch_index.bgen.
    :return: A dictionary containing lists of dxlinks to the final BGEN, index, sample, VEP, and VEP index files for each batch.
    """

    output = {
        'bgen': [],
        'index': [],
        'sample': [],
        'vep': [],
        'vep_idx': []
    }

    subjob_launcher = SubjobUtility(log_update_time=600, incrementor=5, download_on_complete=False)

    for batch_index, batch_files in bgen_chunks.items():
        LOGGER.info(f"Processing batch {batch_index} with {len(batch_files)} chunk files...")

        prefixes = ['bgen', 'bgen_index', 'sample', 'vep', 'vep_index', 'vcfprefix', 'start']
        inputs = {prefix: [] for prefix in prefixes}
        inputs['batch_index'] = batch_index
        inputs['make_bcf'] = make_bcf
        inputs['output_prefix'] = f"{output_prefix}_{batch_index}"

        for file_info in batch_files:
            for prefix in prefixes:
                inputs[prefix].append(file_info[prefix])

        subjob_launcher.launch_job(
            function=process_final_chunk,
            inputs=inputs,
            outputs=['bgen', 'bgen_index', 'sample', 'vep', 'vep_index'],
            instance_type='mem2_ssd1_v2_x32',
            name=f'makebgen_batch{batch_index}',
        )

    subjob_launcher.submit_queue()

    for subjob_output in subjob_launcher:
        # Set output as a list of dxlinks to the final files

        output['bgen'].append(subjob_output['bgen'])
        output['index'].append(subjob_output['bgen_index'])
        output['sample'].append(subjob_output['sample'])
        output['vep'].append(subjob_output['vep'])
        output['vep_idx'].append(subjob_output['vep_index'])

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
        coord_reader = csv.DictReader(coord_file, delimiter="\t")

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
                start=int(row['start']),
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

@dxpy.entry_point('process_final_chunk')
def process_final_chunk(bgen: List[dict], bgen_index: List[dict], sample: List[dict], vep: List[dict], vep_index: List[dict],
                        vcfprefix: List[str], start: List[int],
                        batch_index: int, make_bcf: bool, output_prefix: str) -> dict:
    """Convert all VCFs in a chunk file to BGENs, then merge into one BGEN for that batch.

    This function is very similar to process_single_chunk, but merges all the BGEN files for a batch into a single BGEN file.

    :param bgen: A list containing file pointers for all BGEN files in the batch.
    :param bgen_index: A list containing file pointers for all BGEN index files in the batch.
    :param sample: A list containing file pointers for all sample files in the batch.
    :param vep: A list containing file pointers for all VEP files in the batch.
    :param vep_index: A list containing file pointers for all VEP index files in the batch.
    :param vcfprefix: A list of VCF prefixes for the batch, used to ensure sorting when merging.
    :param start: A list of start coordinates for each VCF file in the batch, used to ensure sorting when merging.
    :param batch_index: The index of the batch.
    :param make_bcf: Should a concatenated BCF be made in addition to the BGEN?
    :param output_prefix: The prefix for the output files.
    :return: A dictionary containing the output files for the batch.
    """
    LOGGER.info(f"Starting processing of batch {batch_index}.")

    batch_files = []
    for i in range(len(bgen)):
        batch_files.append({
            'bgen': InputFileHandler(bgen[i], download_now=True),
            'bgen_index': InputFileHandler(bgen_index[i], download_now=True),
            'sample': InputFileHandler(sample[i], download_now=True),
            'vep': InputFileHandler(vep[i], download_now=True),
            'vep_index': InputFileHandler(vep_index[i], download_now=True),
            'vcfprefix': vcfprefix[i],
            'start': int(start[i])  # We do this to be 200% sure...
        })

    merged = make_final_bgen(bgen_prefixes=batch_files, output_prefix=f"{output_prefix}_{batch_index}",
                             make_bcf=make_bcf)

    output = {
        'bgen': generate_linked_dx_file(merged['bgen']['file']),
        'bgen_index': generate_linked_dx_file(merged['bgen']['index']),
        'sample': generate_linked_dx_file(merged['bgen']['sample']),
        'vep': generate_linked_dx_file(merged['vep']['file']),
        'vep_index': generate_linked_dx_file(merged['vep']['index']),
    }

    return output