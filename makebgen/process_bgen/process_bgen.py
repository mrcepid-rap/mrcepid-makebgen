import gzip
import os
from pathlib import Path
from typing import Dict, List, Any

import dxpy
from general_utilities.association_resources import replace_multi_suffix, bgzip_and_tabix
from general_utilities.import_utils.file_handlers.dnanexus_utilities import generate_linked_dx_file
from general_utilities.import_utils.file_handlers.input_file_handler import FileType
from general_utilities.import_utils.import_lib import InputFileHandler
from general_utilities.job_management.command_executor import CommandExecutor, build_default_command_executor
from general_utilities.job_management.subjob_utility import SubjobUtility, prep_current_image

from makebgen.deduplication.deduplication import deduplicate_variants

# Small note to explain logic – When CMD EXEC is passed to functions, it is to enable functional testing rather than
# end-to-end testing outside the DNANexus environment. During 'normal' execution the global CMD_EXEC is used.
CMD_EXEC = build_default_command_executor()


def make_bgen_from_vcf(vcf_id: Path, vep_id: str, previous_vep_id: str, start: int, make_bcf: bool,
                       cmd_exec: CommandExecutor = CMD_EXEC) -> Dict[str, int]:
    """Downloads the BCF/VEP for a single chunk and processes it.

    Returns a dict of the chunk name and start coordinate. This is so the name can be provided for merging, sorted by
    start coordinate so sorting doesn't have to happen twice.

    :param vcf_id: A string dxid in the form of file-12345... pointing to the VCF to be processed
    :param vep_id: A string dxid in the form of file-12345... pointing to the VEP annotations for the VCF
    :param previous_vep_id: A string dxid in the form of file-12345... pointing to the VEP annotations for the PREVIOUS VCF
    :param start: Start coordinate for this chunk
    :param make_bcf: Is a bcf being made for this chromosome?
    :param cmd_exec: A command executor object to run commands on the docker instance. Default is the global CMD_EXEC.
    :return: A dictionary with key of processed prefix and value of the start coordinate for that bgen
    """

    ### NOTE:
    # try and convert this function to work with pysam & bgen modules if possible

    # download the file
    vcf_path = InputFileHandler(vcf_id).get_file_handle()

    # Set names and DXPY files for bcf/vep file
    # Get a prefix name for all files, the 1st element of the suffixes is ALWAYS the chunk number.
    # Extract this part
    suffix_from_name = "." + vcf_path.name.split(".")[1]  # e.g., '.chunk1'
    # Apply the function
    vcf_prefix = replace_multi_suffix(vcf_path, suffix_from_name).name

    # Download and remove duplicate sites (in both the VEP and BCF) due to erroneous multi-allelic processing by UKBB
    deduplicate_variants(vep_id=vep_id, previous_vep_id=previous_vep_id, vcf_prefix=vcf_prefix, vcf_path=vcf_path)

    # And convert processed bcf into bgenv1.2
    cmd = f'plink2 --threads 2 --memory 10000 ' \
          f'--bcf /test/{vcf_path.name} ' \
          f'--export bgen-1.2 \'bits=\'8 \'ref-first\' ' \
          f'--vcf-half-call r ' \
          f'--out /test/{vcf_prefix} ' \
          f'--output-chr chrMT ' \
          f'--new-id-max-allele-len 1500'
    cmd_exec.run_cmd_on_docker(cmd)

    # Delete the original .bcf from the instance to save space (if we aren't making a bcf later)
    # Note, for now let's only do this if we are working on a VM (i.e. not locally)
    if InputFileHandler(vcf_id).get_file_type() != FileType.LOCAL_PATH:

        if not make_bcf:
            vcf_path.unlink()

    # Return both the processed prefix and the start coordinate to enable easy merging/sorting
    return {'vcfprefix': vcf_prefix,
            'start': start}


def correct_sample_file(template_sample: Path, output_prefix: str) -> Path:
    """Correct the incorrect sample from plink2

    Plink to cannot write a correct sample file, so we need to correct it here.

    plink2 outputs a file like the following:

        ID_1 ID_2 missing sex
        0 0 0 D
        0 1000001 0 NA
        0 1000002 0.0008 NA
        0 1000003 0.0008 NA

    We need it to look like this:

        ID_1 ID_2 missing sex
        0 0 0 D
        1000001 1000001 0 NA
        1000002 1000002 0 NA
        1000003 1000003 0 NA

    The primary change is the addition of the sample ID to ID_1. We also overwrite the missingness column to 0 since
    the information is not calculated correctly by plink2 for all bgen files being merged.

    :param template_sample: The sample file to correct
    :param output_prefix: The prefix for the output file
    :return: The Path object pointing to the corrected sample file
    """

    # Create a sample file for the final bgen
    # plink writes an incorrectly formatted sample file, so we need to create a new one
    final_sample = Path(f'{output_prefix}.sample')
    with final_sample.open('w') as final_writer, \
            template_sample.open('r') as template_reader:

        for line in template_reader:
            split_sample = line.rstrip().split(' ')
            if split_sample[0] == 'ID_1':
                final_writer.write('ID_1 ID_2 missing sex\n')
            elif split_sample[1] == '0':
                final_writer.write('0 0 0 D\n')
            else:
                final_writer.write(f'{split_sample[1]} {split_sample[1]} 0 NA\n')

    return final_sample


def make_final_bgen(bgen_prefixes: List[Dict[str, Any]], output_prefix: str, make_bcf: bool,
                    cmd_exec: CommandExecutor = CMD_EXEC) -> Dict[str, Dict[str, Path]]:
    """Concatenate the final per-chromosome bgen file while processing the .vep.gz annotations into a single .tsv file.

    VEP annotation concatenation ASSUMES that file prefixes are sorted by coordinate. This is a safe assumption as the
    previous steps in the pipeline should ensure coordinate sorted input.

    The output format of this method is a named dictionary wuth three keys, 'bgen', 'vep', and 'bcf' each containing a
    dictionary with keys:
        * 'file' - points to the final bgen, vep, and (if requested) bcf files respectively.
        * 'index' - points to the index of the final bgen, vep, and (if requested) bcf files respectively.

    :param bgen_prefixes: a dictionary of form {vcfprefix: file_info} Where 'file_info' is a dictionary containing AT
        LEAST a key for 'start' indicating the start position of this particular file.
    :param output_prefix: The prefix for the final bgen file. Will be named <output_prefix>.bgen.
    :param make_bcf: Should a bcf be made in addition to bgen? This can lead to VERY long runtimes if the number of
        sites is large.
    :param cmd_exec: A command executor object to run commands on the docker instance. Default is the global CMD_EXEC.
    :return: A named dict containing the final bgen, vep, (and if requested) bcf file paths with associated indices.
    """

    # Set output names here:
    final_bgen = Path(f'{output_prefix}.bgen')
    final_bgen_idx = Path(f'{output_prefix}.bgen.bgi')

    # sort the bgen files by the start position to make sure that they are truly in the correct order, but store
    # the original prefix for the final bgen file - note the test data should be sorted in the opposite order
    # (i.e. ['test_input2', 'test_input1'])
    sorted_bgen_prefixes = [bgens['vcfprefix'] for bgens in sorted(bgen_prefixes, key=lambda x: x['start'])]

    # Create a correct sample file:
    final_sample = correct_sample_file(Path(f'{sorted_bgen_prefixes[0]}.sample'), output_prefix)

    # Create a command line for concatenating bgen files and execute
    print('STARTING BGEN CONCATENATION')
    cmd = ' '.join([f'-g /test/{file}.bgen' for file in sorted_bgen_prefixes])
    cmd = f'cat-bgen {cmd} -og /test/{final_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # And index the bgen for random query later
    print('STARTING BGEN INDEXING')
    cmd = f'bgenix -index -g /test/{final_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # Collect & concat VEP annotations at this time
    print('STARTING VEP CONCATENATION')
    concat_vep = Path(f'{output_prefix}.vep.tsv')
    with concat_vep.open('w') as vep_writer:
        for file_n, bgen_prefix in enumerate(sorted_bgen_prefixes):

            vep_reader = None
            current_vep = Path(f'{bgen_prefix}.vep.tsv')
            if current_vep.exists():
                vep_reader = current_vep.open("r")
            else:
                current_vep = Path(f'{bgen_prefix}.vep.tsv.gz')
                if current_vep.exists():
                    vep_reader = gzip.open(current_vep, "rt")

            with vep_reader:
                for line_n, line in enumerate(vep_reader):
                    if file_n == 0 and line_n == 0:  # Only write header of first file
                        vep_writer.write(line)
                    elif line_n != 0:
                        vep_writer.write(line)

        current_vep.unlink()

    print('STARTING BZIP AND TABIX INDEXING OF VEP ANNOTATIONS')
    # bgzip and tabix index the resulting annotations
    final_vep, final_vep_idx = bgzip_and_tabix(concat_vep, comment_char='C', end_row=2)

    # If make_bcf == True, then we actually do the bcf concatenation
    if make_bcf:
        final_bcf = Path(f'{output_prefix}.bcf')
        final_bcf_idx = Path(f'{output_prefix}.bcf.csi')
        bcf_cmd = ' '.join([f'-g /test/{file}.bcf' for file in sorted_bgen_prefixes])
        bcf_cmd = f'bcftools concat --threads {os.cpu_count()} -a -D -Ob -o /test/{output_prefix}.bcf {bcf_cmd}'
        cmd_exec.run_cmd_on_docker(bcf_cmd)

        # And index:
        bcf_idx = f'bcftools index /test/{output_prefix}.bcf'
        cmd_exec.run_cmd_on_docker(bcf_idx)

    else:
        final_bcf = None
        final_bcf_idx = None

    return {'bgen': {'file': final_bgen, 'index': final_bgen_idx, 'sample': final_sample},
            'vep': {'file': final_vep, 'index': final_vep_idx},
            'bcf': {'file': final_bcf, 'index': final_bcf_idx}}


def make_final_bgen_chunk(bgen_prefixes: List[Dict[str, Any]], output_prefix: str, make_bcf: bool,
                          cmd_exec: CommandExecutor = CMD_EXEC) -> Dict[str, Dict[str, Path]]:
    """Concatenate the final per-chromosome bgen file while processing the .vep.gz annotations into a single .tsv file.

    VEP annotation concatenation ASSUMES that file prefixes are sorted by coordinate. This is a safe assumption as the
    previous steps in the pipeline should ensure coordinate sorted input.

    The output format of this method is a named dictionary wuth three keys, 'bgen', 'vep', and 'bcf' each containing a
    dictionary with keys:
        * 'file' - points to the final bgen, vep, and (if requested) bcf files respectively.
        * 'index' - points to the index of the final bgen, vep, and (if requested) bcf files respectively.

    :param bgen_prefixes: a dictionary of form {vcfprefix: file_info} Where 'file_info' is a dictionary containing AT
        LEAST a key for 'start' indicating the start position of this particular file.
    :param output_prefix: The prefix for the final bgen file. Will be named <output_prefix>.bgen.
    :param make_bcf: Should a bcf be made in addition to bgen? This can lead to VERY long runtimes if the number of
        sites is large.
    :param cmd_exec: A command executor object to run commands on the docker instance. Default is the global CMD_EXEC.
    :return: A named dict containing the final bgen, vep, (and if requested) bcf file paths with associated indices.
    """

    # Set output names here:
    final_bgen = Path(f'{output_prefix}.bgen')
    final_bgen_idx = Path(f'{output_prefix}.bgen.bgi')

    # sort the bgen files by the start position to make sure that they are truly in the correct order, but store
    # the original prefix for the final bgen file - note the test data should be sorted in the opposite order
    # (i.e. ['test_input2', 'test_input1'])
    sorted_bgen_prefixes = [bgens['vcfprefix'] for bgens in sorted(bgen_prefixes, key=lambda x: x['start'])]

    # Create a correct sample file:
    final_sample = correct_sample_file(Path(f'{sorted_bgen_prefixes[0]}.sample'), output_prefix)

    # Create a command line for concatenating bgen files and execute
    print('STARTING BGEN CONCATENATION')
    cmd = ' '.join([f'-g /test/{file}.bgen' for file in sorted_bgen_prefixes])
    cmd = f'cat-bgen {cmd} -og /test/{final_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # Start BGEN indexing
    print('STARTING BGEN INDEXING')
    output = {}
    subjob_launcher = SubjobUtility(log_update_time=600, incrementor=5, download_on_complete=False)
    final_bgen_link = generate_linked_dx_file(final_bgen)
    print("About to start a subjob to index the final BGEN file...")
    subjob_launcher.launch_job(
        function=index_bgen_file,
        inputs={
            'bgen_file': {"$dnanexus_link": final_bgen_link.get_id()},
        },
        outputs=['bgen_index'],
        instance_type='mem3_ssd3_x12',
        name=f'indexing_{final_bgen}',
    )
    subjob_launcher.submit_queue()
    for subjob_output in subjob_launcher:
        output['bgen_index'] = subjob_output['bgen_index']
    # then let's get the actual file
    final_bgen_idx = InputFileHandler(output['bgen_index']).get_file_handle()

    # Collect & concat VEP annotations at this time
    print('STARTING VEP CONCATENATION')
    concat_vep = Path(f'{output_prefix}.vep.tsv')
    with concat_vep.open('w') as vep_writer:
        for file_n, bgen_prefix in enumerate(sorted_bgen_prefixes):

            vep_reader = None
            current_vep = Path(f'{bgen_prefix}.vep.tsv')
            if current_vep.exists():
                vep_reader = current_vep.open("r")
            else:
                current_vep = Path(f'{bgen_prefix}.vep.tsv.gz')
                if current_vep.exists():
                    vep_reader = gzip.open(current_vep, "rt")

            with vep_reader:
                for line_n, line in enumerate(vep_reader):
                    if file_n == 0 and line_n == 0:  # Only write header of first file
                        vep_writer.write(line)
                    elif line_n != 0:
                        vep_writer.write(line)

        current_vep.unlink()

    print('STARTING BZIP AND TABIX INDEXING OF VEP ANNOTATIONS')
    # bgzip and tabix index the resulting annotations
    final_vep, final_vep_idx = bgzip_and_tabix(concat_vep, comment_char='C', end_row=2)

    # If make_bcf == True, then we actually do the bcf concatenation
    if make_bcf:
        final_bcf = Path(f'{output_prefix}.bcf')
        final_bcf_idx = Path(f'{output_prefix}.bcf.csi')
        bcf_cmd = ' '.join([f'-g /test/{file}.bcf' for file in sorted_bgen_prefixes])
        bcf_cmd = f'bcftools concat --threads {os.cpu_count()} -a -D -Ob -o /test/{output_prefix}.bcf {bcf_cmd}'
        cmd_exec.run_cmd_on_docker(bcf_cmd)

        # And index:
        bcf_idx = f'bcftools index /test/{output_prefix}.bcf'
        cmd_exec.run_cmd_on_docker(bcf_idx)

    else:
        final_bcf = None
        final_bcf_idx = None

    return {'bgen': {'file': final_bgen, 'index': final_bgen_idx, 'sample': final_sample},
            'vep': {'file': final_vep, 'index': final_vep_idx},
            'bcf': {'file': final_bcf, 'index': final_bcf_idx}}


@dxpy.entry_point('index_bgen_file')
def index_bgen_file(bgen_file: dict) -> dict:
    """
    Index a BGEN file using bgenix.

    :param bgen_file: A dictionary with a key 'bgen_file' pointing to the BGEN file to be indexed.
    :return: A dictionary with a key 'bgen_index' pointing to the indexed BGEN file.
    """

    cmd_exec = prep_current_image([bgen_file])
    final_bgen = InputFileHandler(bgen_file).get_file_handle()
    print('Indexing BGEN file:', final_bgen)
    cmd = f'bgenix -index -g /test/{final_bgen.name}'

    cmd_exec.run_cmd_on_docker(cmd)

    output = {}
    generate_linked_dx_file(final_bgen)
    output['bgen_index'] = generate_linked_dx_file(final_bgen)
    print('Indexed BGEN file:', output['bgen_index'])

    return output
