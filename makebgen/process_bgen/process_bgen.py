import os
from pathlib import Path
from typing import Dict, Union

import dxpy
from general_utilities.association_resources import download_dxfile_by_name, replace_multi_suffix, bgzip_and_tabix
from general_utilities.import_utils.import_lib import input_filetype_parser
from general_utilities.job_management.command_executor import CommandExecutor, build_default_command_executor

from makebgen.deduplication.deduplication import deduplicate_variants

# Small note to explain logic – When CMD EXEC is passed to functions, it is to enable functional testing rather than
# end-to-end testing outside the DNANexus environment. During 'normal' execution the global CMD_EXEC is used.
CMD_EXEC = build_default_command_executor()


def make_bgen_from_vcf(vcf_id: Union[str, Path], vep_id: str, previous_vep_id: str, start: int, make_bcf: bool,
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

    if isinstance(input_filetype_parser(vcf_id), dxpy.DXFile):

        vcf_path = download_dxfile_by_name(vcf_id, print_status=False)

        # Set names and DXPY files for bcf/vep file
        # Get a prefix name for all files, the 1st element of the suffixes is ALWAYS the chunk number.
        vcf_prefix = replace_multi_suffix(vcf_path, vcf_path.suffixes[0])

        # Download and remove duplicate sites (in both the VEP and BCF) due to erroneous multi-allelic processing by UKBB
        deduplicate_variants(vep_id=vep_id, previous_vep_id=previous_vep_id, vcf_prefix=vcf_prefix, vcf_path=vcf_path)

    else:

        # if not using DNA Nexus we can assign a local path to the vcf_path variable
        vcf_path = Path(vcf_id)

        # also, we can change the filepath so that we keep the prefix only
        vcf_prefix = replace_multi_suffix(vcf_path, vcf_path.suffixes[0]).stem

    # And convert processed bcf into bgenv1.2
    cmd = f'plink2 --threads 2 --memory 10000 ' \
          f'--bcf /test/{vcf_path} ' \
          f'--export bgen-1.2 \'bits=\'8 \'ref-first\' ' \
          f'--vcf-half-call r ' \
          f'--out /test/{vcf_prefix} ' \
          f'--output-chr chrMT ' \
          f'--new-id-max-allele-len 1500'
    cmd_exec.run_cmd_on_docker(cmd)

    # Delete the original .bcf from the instance to save space (if we aren't making a bcf later)
    # Note, for now let's only do this if we are running on DNA Nexus

    if isinstance(input_filetype_parser(vcf_id), dxpy.DXFile):

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


def make_final_bgen(bgen_prefixes: dict, output_prefix: str, make_bcf: bool,
                    cmd_exec: CommandExecutor = CMD_EXEC) -> Dict[str, Dict[str, Path]]:
    """Concatenate the final per-chromosome bgen file while processing the .vep.gz annotations into a single .tsv file.

    VEP annotation concatenation ASSUMES that file prefixes are sorted by coordinate. This is a safe assumption as the
    previous steps in the pipeline should ensure coordinate sorted input.

    The output format of this method is a named dictionary wuth three keys, 'bgen', 'vep', and 'bcf' each containing a
    dictionary with keys:
        * 'file' - points to the final bgen, vep, and (if requested) bcf files respectively.
        * 'index' - points to the index of the final bgen, vep, and (if requested) bcf files respectively.

    :param bgen_prefixes: a dictionary of form {vcfprefix: start_coordinate}
    :param output_prefix: The prefix for the final bgen file. Will be named <output_prefix>.bgen.
    :param make_bcf: Should a bcf be made in addition to bgen? This can lead to VERY long runtimes if the number of
        sites is large.
    :param cmd_exec: A command executor object to run commands on the docker instance. Default is the global CMD_EXEC.
    :return: A named dict containing the final bgen, vep, (and if requested) bcf file paths with associated indices.
    """

    # Set output names here:
    final_bgen = Path(f'{output_prefix}.bgen')
    final_bgen_idx = Path(f'{output_prefix}.bgen.bgi')

    # Sort the bgen files according to coordinate
    sorted_bgen_prefixes = sorted(bgen_prefixes)

    # Create a correct sample file:
    final_sample = correct_sample_file(Path(f'{sorted_bgen_prefixes[0]}.sample'), output_prefix)

    # Create a command line for concatenating bgen files and execute
    cmd = ' '.join([f'-g /test/{file}.bgen' for file in sorted_bgen_prefixes])
    cmd = f'cat-bgen {cmd} -og /test/{final_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # And index the bgen for random query later
    cmd = f'bgenix -index -g /test/{final_bgen}'
    cmd_exec.run_cmd_on_docker(cmd)

    # Collect & concat VEP annotations at this time
    concat_vep = Path(f'{output_prefix}.vep.tsv')
    with concat_vep.open('w') as vep_writer:
        all_lines = []
        for file_n, bgen_prefix in enumerate(sorted_bgen_prefixes):
            current_vep = Path(f'{bgen_prefix}.vep.tsv')
            with current_vep.open('r') as vep_reader:
                for line_n, line in enumerate(vep_reader):
                    if file_n == 0 and line_n == 0:  # Only write header of first file
                        vep_writer.write(line)
                    elif line_n != 0:
                        all_lines.append(line)

            current_vep.unlink()

        # Sort the collected lines (excluding the header)
        sorted_lines = sorted(all_lines, key=lambda vep_line: vep_line.split('\t')[1])
        vep_writer.writelines(sorted_lines)

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
