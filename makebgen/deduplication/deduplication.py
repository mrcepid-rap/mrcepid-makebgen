import gzip
from pathlib import Path
from typing import Tuple, Optional, Union

import dxpy
import pandas as pd
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import build_default_command_executor, CommandExecutor
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger().get_logger()
CMD_EXEC = build_default_command_executor()


def deduplicate_variants(vep_id: Union[str, Path], previous_vep_id: Union[str, Path], vcf_prefix: Union[str, Path],
                         vcf_path: Path, cmd_exec: CommandExecutor = CMD_EXEC) -> Path:
    """Entry point into the various deduplication methods in this module. This method only handles the logic flow of
    running the other methods.

    :param vep_id: A string dxid in the form of file-12345... pointing to the VEP annotations for the VCF
    :param previous_vep_id: A string dxid in the form of file-12345... pointing to the VEP annotations for the PREVIOUS VCF
    :param vcf_prefix: A Path object pointing to the VCF prefix of the file to deduplicate
    :param vcf_path: A Path object pointing to the VCF file to deduplicate
    :param cmd_exec: A command executor object to run commands on the docker instance. Default is the global CMD_EXEC.
    :return: Path to the written VEP annotation file.
    """

    # Load relevant VEP annotations
    current_df = load_vep(vep_id)
    previous_df = load_vep(previous_vep_id)
    deduped_df, removed_df = remove_vep_duplicates(current_df, previous_df)

    if len(removed_df) != 0:  # Only do the bcf if we have ≥ 1 variant to exclude
        LOGGER.warning(f'BCF with prefix {vcf_prefix} has {len(removed_df)} duplicate variants. Removing...')
        remove_query_string = build_query_string(removed_df)
        remove_bcf_duplicates(query_string=remove_query_string, vcf_prefix=vcf_prefix, vcf_path=vcf_path, cmd_exec=cmd_exec)

    return write_vep_table(deduped_df, vcf_prefix)


def load_vep(vep_id: str) -> Optional[pd.DataFrame]:
    """Read a VEP annotation into a pandas DataFrame.

    This method is a simple wrapper for :func:`pd.read_csv` which will *stream* a vep annotation for the given
    chromosome rather than download using :func:`dxpy.open_dxfile`.

    For the first VEP annotation file for a given chromosome, this method will return an empty DataFrame as it is not
    possible for duplicates to exist.

    To be clear, the reason this method has an optional return, rather than a set column definition, is because the
    names of the columns are not guaranteed to be the same between different runs. This is due to additional annotations
    being parameterized.

    :param vep_id: The dxid of the VEP annotation file.
    tooling) or not. The default is set to True.
    :return: An optional pandas DataFrame of the VEP annotation. If vep_id is None, return None.
    """

    if vep_id is None:
        return None
    else:
        vep_id = InputFileHandler(vep_id).get_file_handle()
        # We stream the .vep annotation from dnanexus as it is faster for smaller files like this, and ensures that we don't
        # generate race conditions when loading the same file twice (which we do to check for duplicates...)
        #
        # Note to future devs – DO NOT remove gzip even though pandas can direct read gzip. It is not compatible with
        # dxpy.open_dxfile and will error out.

        current_df = pd.read_csv(gzip.open(Path(vep_id)), sep="\t", index_col=False)

        # This is legacy naming of the ID column and too difficult to refactor in the downstream pipeline
        current_df.rename(columns={'ID': 'varID'}, inplace=True)

        return current_df


# Helper function for process_vep() to remove duplicate variants
def remove_vep_duplicates(current_vep_df: pd.DataFrame, previous_vep_df: Optional[pd.DataFrame]) \
        -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Check for duplicate variants in the VEP annotation and remove them.

    See in-line comments for specific details on how duplicates are removed. In brief, we remove two types of duplicates:

    1. Those in the current VEP annotation file – these are typically InDels that are multi-allelic with a SNP but
        were stored as separate variants twice in the UKBB-provided VCF.

    2. Those in the previous VEP annotation file – the same as (1), but one multi-allelic variant was stored
        in the previous VCF and the other in the current VCF.

    This method will also write the final VEP annotation file to disk to be used by later methods in this workflow. It
    is a flat .tsv file and NOT gzipped like the original VEP annotation file.

    :param current_vep_df: A pandas DataFrame of the current VEP annotation.
    :param previous_vep_df: A pandas DataFrame of the previous VEP annotation.
    :return: A pandas DataFrame of the removed duplicate variants.
    """

    # We have two types of duplicate variants...
    removed_variants = []

    # 1. Variants that are duplicates INTERNAL to this file
    var_counts = current_vep_df['varID'].value_counts()
    duplicates = var_counts[var_counts == 2]

    # This iterates over all duplicates in this file...
    # [0] = the variant ID
    # [1] = the number of times the variant ID appears
    for dup in duplicates.items():
        dup_index = current_vep_df[current_vep_df['varID'] == dup[0]].index.to_list()
        dup_rows = current_vep_df.iloc[dup_index]

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
        removed_variants.append(current_vep_df[current_vep_df.index.isin(dup_index)])
        current_vep_df = current_vep_df[current_vep_df.index.isin(dup_index) == False]
        current_vep_df = current_vep_df.reset_index(drop=True)

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
    if previous_vep_df is not None:
        duplicates = current_vep_df[current_vep_df['varID'].isin(previous_vep_df['varID'].to_list())]
        current_vep_df = current_vep_df[current_vep_df.index.isin(duplicates.index.to_list()) == False]

    # Merge the duplicate variants from both mode 1. & 2.:
    removed_variants = pd.concat([removed_variants, duplicates])

    return current_vep_df, removed_variants


def write_vep_table(deduplicated_vep: pd.DataFrame, vcf_prefix: Path) -> Path:
    """Write a vep file to disk.

    :param deduplicated_vep: A pandas DataFrame of the deduplicated VEP annotation.
    :param vcf_prefix: A Path object pointing to the VCF prefix of the file to deduplicate.
    :return: Path to the written VEP annotation file.
    """

    # ensure that vcf_prefix is a path-like object
    vcf_prefix = Path(vcf_prefix)
    # Add the suffix to the stem
    vcf_path = vcf_prefix.parent / (vcf_prefix.name + '.vep.tsv')
    # export the data using the new filepath
    deduplicated_vep.to_csv(path_or_buf=vcf_path, na_rep='NA', index=False, sep="\t")

    # vcf_path = vcf_prefix.with_suffix('.'.join(vcf_prefix.suffixes + ['vep','tsv']))
    # deduplicated_vep.to_csv(path_or_buf=vcf_path, na_rep='NA', index=False, sep="\t")

    return vcf_path


def build_query_string(removed_df: pd.DataFrame) -> str:
    """Build a query string for BCFTools view of duplicate variants in the BCF to remove.

    This method is a loop over the dataframe to concatenate individual query strings into one long string.

    :param removed_df: A pandas DataFrame of duplicate variants to remove.
    :return: A query string compatible with bcftools view.
    """

    # We are going to try to do this as one very large bcftools query if we find multiple variants
    # This will hopefully speed up rare instances of multiple variants needing to be removed from a given bcf
    query_builder = []

    for row in removed_df.iterrows():
        table = row[1]
        position = table['POS']
        ref = table['REF']
        alt = table['ALT']
        ma_id = table['MA']

        # Generate a query that matches on position, ref, alt, and og_id. I am fairly certain (I have checked multiple
        # files) that all duplicates come from different multi-allelics (e.g., MA), so this _should_ be safe...
        query = f'(POS={position} && REF="{ref}" && ALT="{alt}" && MA="{ma_id}")'
        query_builder.append(query)

    query = ' || '.join(query_builder)

    return query


# Helper function for make_bgen_from_vcf() to remove sites identified as duplicate in process_vep() from the bcf
def remove_bcf_duplicates(query_string: str, vcf_prefix: Path, vcf_path: Path,
                          cmd_exec: CommandExecutor = CMD_EXEC) -> Path:
    """Remove duplicate variants using bcftools view.

    This method is just a wrapper for the bcftools view command to remove duplicate variants from the BCF file. Following
    deduplication, the original BCF file is deleted and the deduplicated BCF is renamed to the original BCF name.

    :param query_string: The query string to pass to bcftools view from :func:`build_query_string`.
    :param vcf_prefix: A Path object pointing to the VCF prefix of the file to deduplicate.
    :param vcf_path: A Path object pointing to the VCF file to deduplicate.
    :param cmd_exec: A command executor object to run commands on the docker instance. Default is the global CMD_EXEC.
    :return: Path to the deduplicated BCF file.
    """
    print(vcf_path)

    cmd = f'bcftools view --threads 2 -e \'{query_string}\' -Ob -o /test/{vcf_prefix}.deduped.bcf /test/{vcf_path.name}'
    cmd_exec.run_cmd_on_docker(cmd)

    # And rename the deduped bcf to match the final one we want to output
    old_vcf = Path(f'{vcf_path}')
    old_vcf.unlink()
    return Path(f'{vcf_prefix}.deduped.bcf').rename(old_vcf)
