#!/usr/bin/env python
# mrcepid-collapsevariants 0.0.1
# Generated by dx-app-wizard.
#
# Author: Eugene Gardner (eugene.gardner at mrc.epid.cam.ac.uk)
#
# DNAnexus Python Bindings (dxpy) documentation:
#   http://autodoc.dnanexus.com/bindings/python/current/

import re
import csv
import dxpy
import gzip
import pandas as pd
import numpy as np
from pathlib import Path

from general_utilities.job_management.thread_utility import ThreadUtility
from general_utilities.mrc_logger import MRCLogger
from general_utilities.association_resources import run_cmd, generate_linked_dx_file


LOGGER = MRCLogger().get_logger()


# Download necessary resources for this applet
def ingest_resources() -> None:

    cmd = 'docker pull egardner413/mrcepid-burdentesting:latest'
    run_cmd(cmd, is_docker=False)


# Returns a DNANexus file-pointer for the previous VEP index in a file
def get_previous_vep(vcfprefix: str) -> str:

    coord_file_reader = pd.read_csv('coordinates.files.tsv.gz', delimiter="\t")

    # Need to find the correct file
    # Doing this regex just to make sure that weird prefixes/suffices don't screw up our search space
    vep_match = re.search('(ukb\\d+_c[\dXY]{1,2}_b\\d+_v\\d+_chunk\\d+)', vcfprefix)
    search_string = vep_match.group(1)

    # Get the index of the current file
    current_index = coord_file_reader[coord_file_reader['fileprefix'] == search_string].index.to_list()[0]

    # subtract one to get previous vep index and search...
    # This is going to do something VERY dumb for the very first file in the index (chr1_b1_chunk1), which is return the
    # very last vep index (chrY) due to pandas indexing conventions. I don't think this matters as none of the variants
    # will be found in the previous index, so similar to what will be done for the first VEP on each chromosome.
    matched_vep = coord_file_reader.iloc[current_index - 1]['vep_dxpy']

    return matched_vep


# Code will *stream* a vep annotation for the given chromosome rather than download, and perform processing to
# get it into the correct format for merging.
def load_vep(vep: str) -> pd.DataFrame:

    vep_header = ['CHROM', 'POS', 'REF', 'ALT', 'ogVarID', 'FILTER', 'AF', 'F_MISSING', 'AN', 'AC', 'MANE',
                  'ENST', 'ENSG', 'BIOTYPE', 'SYMBOL', 'CSQ', 'gnomAD_AF', 'CADD', 'REVEL', 'SIFT', 'POLYPHEN',
                  'LOFTEE', 'AA', 'AApos', 'PARSED_CSQ', 'MULTI', 'INDEL', 'MINOR', 'MAJOR', 'MAF', 'MAC']

    # We stream the .vep annotation from dnanexus as it is faster for smaller files like this, and ensures that we don't
    # generate race conditions when loading the same file twice (which we do to check for duplicates...)
    current_df = pd.read_csv(gzip.open(dxpy.open_dxfile(vep, mode='rb'), 'rt'), sep="\t", header=None, names=vep_header)

    # Need to add a column with a modifed variant ID that is required during downstream processing
    current_df['newCHROM'] = current_df['CHROM'].apply(lambda x: x.strip('chr'))
    current_df['varID'] = current_df.apply(lambda row: f'{row["newCHROM"]}:{row["POS"]}:{row["REF"]}:{row["ALT"]}',
                                           axis=1)
    current_df = current_df.drop(columns={'newCHROM'})

    # Make sure varID is in the correct place in the table when writing
    vep_header.insert(4, 'varID')
    current_df = current_df[vep_header]

    # For some reason AF is not corrected in the previous step (filterbcf) when every genotype is missing
    # (i.e. F_MISSING == 1).
    current_df['AF'] = current_df['AF'].apply(lambda x: x if x != '.' else 0)
    current_df['AF'] = current_df['AF'].astype(np.float64)

    return current_df


# Helper function for process_vep() to remove duplicate variants
def deduplicate_vep(current_df: pd.DataFrame, previous_df: pd.DataFrame, vcfprefix: str) -> pd.DataFrame:

    # We have two types of duplicate variants...
    removed_variants = []

    # 1. Variants that are duplicates INTERNAL to this file
    var_counts = current_df['varID'].value_counts()
    duplicates = var_counts[var_counts == 2]

    # This iterates over all duplicates in this file...
    # [0] = the variant ID
    # [1] = the number of times the variant ID appears
    for dup in duplicates.iteritems():
        dup_index = current_df[current_df['varID'] == dup[0]].index.to_list()
        dup_rows = current_df.iloc[dup_index]

        # Then we iterate through each duplicate and select the 'best' variant based on the following hierarchy:
        # 1. PASS/FAIL
        # 2. Lower missingness
        # 3. Higher MAF (logic here is that the 'better' call is on the correct genotype and will happen more often.
        #    I'm not going to mess around with trying to merge genotypes, which is likely to be the best solution.
        keep_index = None
        current_filter = None
        current_missing = None
        current_af = None
        for row in dup_rows.iterrows():
            index = row[0]
            table = row[1]
            if keep_index is None:
                keep_index = index
                current_filter = table['FILTER']
                current_missing = table['F_MISSING']
                current_af = table['AF']
            else:
                if current_filter == 'FAIL' and table['FILTER'] == 'PASS':
                    keep_index = index
                    current_filter = table['FILTER']
                    current_missing = table['F_MISSING']
                    current_af = table['AF']
                elif current_missing > table['F_MISSING']:
                    keep_index = index
                    current_filter = table['FILTER']
                    current_missing = table['F_MISSING']
                    current_af = table['AF']
                elif current_af < table['AF']:
                    keep_index = index
                    current_filter = table['FILTER']
                    current_missing = table['F_MISSING']
                    current_af = table['AF']

        dup_index.remove(keep_index)
        removed_variants.append(current_df[current_df.index.isin(dup_index)])
        current_df = current_df[current_df.index.isin(dup_index) == False]
        current_df = current_df.reset_index(drop=True)

    # Need to create an empty pd.DataFrame so that we can return something for deduplicate_bcf() to use.
    # deduplicate_bcf() will just return immediately for an empty DataFrame, so no harm here.
    if len(removed_variants) == 0:
        removed_variants = pd.DataFrame()
    else:
        removed_variants = pd.concat(removed_variants)

    # 2. Variants that were represented in the previous file
    #    Just a note on this: Unfortunately, I don't think I can retrospectively go back and remove variants from the
    #    previous file by identical criteria. I just have to go with 'remove the second instance' for this type of
    #    duplicate.
    duplicates = current_df[current_df['varID'].isin(previous_df['varID'].to_list())]
    current_df = current_df[current_df.index.isin(duplicates.index.to_list()) == False]

    # Merge the duplicate variants from both mode 1. & 2.:
    removed_variants = pd.concat([removed_variants, duplicates])

    # And write the final vep file
    current_df.to_csv(path_or_buf=f'{vcfprefix}.vep.tsv', na_rep='NA', index=False, sep="\t")

    return removed_variants


# Helper function for make_bgen_from_vcf() to remove sites identified as duplicate in process_vep() from the bcf
def deduplicate_bcf(vcfprefix: str, removed_df: pd.DataFrame) -> None:

    # We are going to try to do this as one very large bcftools query if we find multiple variants
    # This will hopefully speed up rare instances of multiple variants needing to be removed from a given bcf
    query_builder = []

    for row in removed_df.iterrows():
        table = row[1]
        position = table['POS']
        ref = table['REF']
        alt = table['ALT']
        og_id = table['ogVarID']

        # Generate a query that matches on position, ref, alt, and og_id. I am fairly certain (I have checked multiple
        # files) that all duplicates come from different multi-allelics, so this _should_ be safe...
        query = f'(POS={position} && REF="{ref}" && ALT="{alt}" && ID="{og_id}")'
        query_builder.append(query)

    query = ' || '.join(query_builder)

    cmd = f'bcftools view --threads 2 -e \'{query}\' -Ob -o /test/{vcfprefix}.deduped.bcf /test/{vcfprefix}.bcf'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

    # And rename the deduped bcf to match the final one we want to output
    old_vcf = Path(f'{vcfprefix}.bcf')
    old_vcf.unlink()
    Path(f'{vcfprefix}.deduped.bcf').rename(old_vcf)


# Downloads the BCF/VEP for a single chunk and processes it.
# Returns a dict of the chunk name and start coordinate. This is so the name can be provided for merging, sorted by
# start coordinate so sorting doesn't have to happen twice
def make_bgen_from_vcf(vcf_id: str, vep_id: str, start: int) -> dict:

    vcf = dxpy.DXFile(vcf_id.rstrip())

    # Set names and DXPY files for bcf/vep file
    LOGGER.info(f'Processing bcf: {vcf.describe()["name"]}')
    vcfprefix = vcf.describe()['name'].rstrip('.bcf')  # Get a prefix name for all files
    dxpy.download_dxfile(vcf.get_id(), f'{vcfprefix}.bcf')

    # Download and remove duplicate sites (in both the VEP and BCF) due to erroneous multi-allelic processing by UKBB
    current_df = load_vep(vep_id)
    previous_df = load_vep(get_previous_vep(vcfprefix))
    removed_df = deduplicate_vep(current_df, previous_df, vcfprefix)
    if len(removed_df) != 0:  # Only do the bcf if we have ≥ 1 variant to exclude
        deduplicate_bcf(vcfprefix, removed_df)

    # And convert processed bcf into bgenv1.2
    cmd = f'plink2 --threads 2 --memory 10000 ' \
          f'--bcf /test/{vcfprefix}.bcf ' \
          f'--export bgen-1.2 \'bits=\'8 \'sample-v2\' ' \
          f'--vcf-half-call r ' \
          f'--out /test/{vcfprefix} ' \
          f'--set-all-var-ids @:#:\$r:\$a ' \
          f'--new-id-max-allele-len 500'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

    # Delete the original .bcf from the instance to save space
    Path(f'{vcfprefix}.bcf').unlink()

    # Return both the processed prefix and the start coordinate to enable easy merging/sorting
    return {'vcfprefix': vcfprefix,
            'start': start}


# Concatenate the final per-chromosome bgen file while processing the .vep.gz annotations
# bgen_prefixes is a dictionary of form {vcfprefix: start_coordinate}
def make_final_bgen(bgen_prefixes: dict, chromosome: str) -> None:

    vep_files = []  # A list containing pandas DataFrames that will be concatenated together

    # How this works: we just keep adding each chunk .bgen file to the end of the cat-bgen command since it doesn't
    # appear to accept a file-list as an input, rather than command-line
    cmd = "cat-bgen"
    sorted_bgen_prefixes = sorted(bgen_prefixes)

    for i in range(0, len(sorted_bgen_prefixes)):
        if i == 0:
            cp_cmd = f'cp {sorted_bgen_prefixes[i]}.sample {chromosome}.filtered.sample'
            run_cmd(cp_cmd, is_docker=False)
        cmd += f' -g /test/{sorted_bgen_prefixes[i]}.bgen'

        # Collect VEP annotations at this time
        vep_files.append(pd.read_csv(sorted_bgen_prefixes[i] + ".vep.tsv", sep="\t"))
        Path(f'{sorted_bgen_prefixes[i]}.vep.tsv').unlink()

    # Add the output filename and concatenate
    cmd += f' -og /test/ {chromosome}.filtered.bgen'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

    # And index the bgen for random query later
    cmd = f'bgenix -index -g /test/{chromosome}.filtered.bgen'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')

    # Mash the vep files together and write to .tsv format:
    vep_index = pd.concat(vep_files)
    vep_index = vep_index.sort_values('POS')  # Make sure sorted
    vep_index.to_csv(path_or_buf=f'{chromosome}.filtered.vep.tsv', na_rep='NA', index=False, sep="\t")

    # bgzip and tabix index the resulting annotations
    cmd = f'bgzip /test/{chromosome}.filtered.vep.tsv'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')
    cmd = f'tabix -c C -s 1 -b 2 -e -2 /test/{chromosome}.filtered.vep.tsv.gz'
    run_cmd(cmd, is_docker=True, docker_image='egardner413/mrcepid-burdentesting:latest')


@dxpy.entry_point('main')
def main(chromosome, coordinate_file):

    # Grab resources required for this applet
    LOGGER.info(f'Ingesting applet resources')
    ingest_resources()

    # Convert all input bcfs to bgen
    LOGGER.info(f'Converting chromosome {chromosome} bcf(s) to bgen(s)')
    # Get the processed coordinate file
    coordinate_file = dxpy.DXFile(coordinate_file)
    dxpy.download_dxfile(coordinate_file.get_id(), "coordinates.files.tsv.gz")
    coord_file_reader = csv.DictReader(gzip.open('coordinates.files.tsv.gz', mode='rt'), delimiter="\t")

    thread_utility = ThreadUtility(incrementor=10,
                                   thread_factor=2,
                                   error_message='A bcf to bgen thread failed')
    for row in coord_file_reader:
        if row['#chrom'] == chromosome:
            thread_utility.launch_job(make_bgen_from_vcf,
                                      vcf_id=row['bcf_dxpy'],
                                      vep_id=row['vep_dxpy'],
                                      start=row['start'])

    # And gather the resulting futures which are returns of all bgens we need to concatenate:
    bgen_prefixes = {}
    for result in thread_utility:
        bgen_prefixes[result['vcfprefix']] = result['start']

    # Now mash all the bgen files together
    LOGGER.info(f'Merging bgen files for chromosome {chromosome} together...')
    make_final_bgen(bgen_prefixes, chromosome)

    # Set output
    output = {"bgen": dxpy.dxlink(generate_linked_dx_file(f'{chromosome}.filtered.bgen')),
              "index": dxpy.dxlink(generate_linked_dx_file(f'{chromosome}.filtered.bgen.bgi')),
              "sample": dxpy.dxlink(generate_linked_dx_file(f'{chromosome}.filtered.sample')),
              "vep": dxpy.dxlink(generate_linked_dx_file(f'{chromosome}.filtered.vep.tsv.gz')),
              "vep_idx": dxpy.dxlink(generate_linked_dx_file(f'{chromosome}.filtered.vep.tsv.gz.tbi'))}

    return output


dxpy.run()
