#!/usr/bin/env python
# mrcepid-collapsevariants 0.0.1
# Generated by dx-app-wizard.
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/
import os
import csv
import dxpy
import shutil
from typing import Dict
from pathlib import Path

from makebgen.deduplication.deduplication import deduplicate_variants
from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.job_management.command_executor import build_default_command_executor, CommandExecutor
from general_utilities.mrc_logger import MRCLogger
from general_utilities.association_resources import (
    generate_linked_dx_file,
    bgzip_and_tabix,
    download_dxfile_by_name,
    replace_multi_suffix,
    check_gzipped
)

LOGGER = MRCLogger().get_logger()
# Small note to explain logic – When CMD EXEC is passed to functions, it is to enable functional testing rather than
# end-to-end testing outside the DNANexus environment. During 'normal' execution the global CMD_EXEC is used.
CMD_EXEC = build_default_command_executor()


def make_bgen_from_vcf(vcf_id: str, vep_id: str, previous_vep_id: str, start: int, make_bcf: bool,
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

    vcf_path = download_dxfile_by_name(vcf_id, print_status=False)

    # Set names and DXPY files for bcf/vep file
    # Get a prefix name for all files, the 1st element of the suffixes is ALWAYS the chunk number.
    vcf_prefix = replace_multi_suffix(vcf_path, vcf_path.suffixes[0])

    # Download and remove duplicate sites (in both the VEP and BCF) due to erroneous multi-allelic processing by UKBB
    deduplicate_variants(vep_id, previous_vep_id, vcf_prefix)

    # And convert processed bcf into bgenv1.2
    cmd = f'plink2 --threads 2 --memory 10000 ' \
          f'--bcf /test/{vcf_path} ' \
          f'--export bgen-1.2 \'bits=\'8 ' \
          f'--vcf-half-call r ' \
          f'--out /test/{vcf_prefix} ' \
          f'--new-id-max-allele-len 1500'
    cmd_exec.run_cmd_on_docker(cmd)

    # Delete the original .bcf from the instance to save space (if we aren't making a bcf later)
    if not make_bcf:
        vcf_path.unlink()

    # Return both the processed prefix and the start coordinate to enable easy merging/sorting
    return {'vcfprefix': vcf_prefix,
            'start': start}


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

    # Create a sample file for the final bgen
    final_sample = Path(f'{output_prefix}.sample')
    shutil.copy(Path(f'{sorted_bgen_prefixes[0]}.sample'),
                final_sample)

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
        for file_n, bgen_prefix in enumerate(sorted_bgen_prefixes):

            current_vep = Path(f'{bgen_prefix}.vep.tsv')
            with current_vep.open('r') as vep_reader:
                for line_n, line in enumerate(vep_reader):
                    if file_n == 0 and line_n == 0:  # Only write header of first file
                        vep_writer.write(line)
                    elif line_n != 0:
                        vep_writer.write(line)

        current_vep.unlink()

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


@dxpy.entry_point('main')
def main(output_prefix: str, coordinate_file: dict, make_bcf: bool) -> dict:
    """Main entry point into this applet. This function initiates the conversion of all bcf files for a given chromosome
    into a single .bgen file.

    Coordinate file must have the following columns:

        chrom   start   end     vcf_prefix      output_bcf      output_bcf_idx  output_vep      output_vep_idx

    :param output_prefix: Output prefix. Output file will be named <output_prefix>.bgen
    :param coordinate_file: A file containing the coordinates of all bcf files to be processed.
    :param make_bcf: Should a concatenated bcf be made in addition to the bgen?
    :return: An output dictionary following DNANexus conventions.
    """

    # Get the processed coordinate file
    total_bcf = 0
    coordinate_path = download_dxfile_by_name(coordinate_file)
    with check_gzipped(coordinate_path) as coord_file:
        coord_file_reader = csv.DictReader(coord_file, delimiter="\t")

        thread_utility = ThreadUtility(incrementor=100,
                                       thread_factor=2,
                                       error_message='A bcf to bgen thread failed')
        previous_vep_id = None
        for row in coord_file_reader:
            total_bcf+=1
            thread_utility.launch_job(make_bgen_from_vcf,
                                      vcf_id=row['output_bcf'],
                                      vep_id=row['output_vep'],
                                      previous_vep_id=previous_vep_id,
                                      start=row['start'],
                                      make_bcf=make_bcf)
            previous_vep_id = row['output_vep']

        # And gather the resulting futures which are returns of all bgens we need to concatenate:
        bgen_prefixes = {}
        for result in thread_utility:
            bgen_prefixes[result['vcfprefix']] = result['start']

    # Convert all input bcfs to bgen
    LOGGER.info(f'Converting {total_bcf} bcf(s) to single bgen')

    # Now mash all the bgen files together
    LOGGER.info(f'Merging bgen files into {output_prefix}.bgen...')
    final_files = make_final_bgen(bgen_prefixes, output_prefix, make_bcf)

    # Set output
    output = {'bgen': dxpy.dxlink(generate_linked_dx_file(final_files['bgen']['file'])),
              'index': dxpy.dxlink(generate_linked_dx_file(final_files['bgen']['index'])),
              'sample': dxpy.dxlink(generate_linked_dx_file(final_files['bgen']['sample'])),
              'vep': dxpy.dxlink(generate_linked_dx_file(final_files['vep']['file'])),
              'vep_idx': dxpy.dxlink(generate_linked_dx_file(final_files['vep']['index']))}

    if final_files['bcf']['file'] is not None:
        output['bcf'] = dxpy.dxlink(generate_linked_dx_file(final_files['bcf']['file']))
        output['bcf_idx'] = dxpy.dxlink(generate_linked_dx_file(final_files['bcf']['index']))

    return output


dxpy.run()
