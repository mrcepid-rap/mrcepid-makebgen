import gzip
import dxpy
import pandas as pd

from pathlib import Path
from general_utilities.job_management.command_executor import build_default_command_executor
from general_utilities.mrc_logger import MRCLogger

LOGGER = MRCLogger().get_logger()
CMD_EXEC = build_default_command_executor()


def deduplicate_variants(vep_id: str, previous_vep_id: str, vcf_prefix: Path) -> None:
    """Entry point into the various deduplication methods in this module. This method only handles the logic flow of
    running the other methods.

    :param vep_id: A string dxid in the form of file-12345... pointing to the VEP annotations for the VCF
    :param previous_vep_id: A string dxid in the form of file-12345... pointing to the VEP annotations for the PREVIOUS VCF
    :param vcf_prefix: A Path object pointing to the VCF prefix of the file to deduplicate
    """

    # Load relevant VEP annotations
    current_df = load_vep(vep_id)
    previous_df = load_vep(previous_vep_id)
    removed_df = remove_vep_duplicates(current_df, previous_df, vcf_prefix)

    if len(removed_df) != 0:  # Only do the bcf if we have ≥ 1 variant to exclude
        LOGGER.warning(f'BCF with prefix {vcf_prefix} has {len(removed_df)} duplicate variants. Removing...')
        LOGGER.warning(removed_df)
        remove_query_string = build_query_string(removed_df)
        remove_bcf_duplicates(remove_query_string, vcf_prefix)


def load_vep(vep_id: str) -> pd.DataFrame:
    """Read a VEP annotation into a pandas DataFrame.

    This method is a simple wrapper for :func:`pd.read_csv` which will *stream* a vep annotation for the given
    chromosome rather than download using :func:`dxpy.open_dxfile`.

    For the first VEP annotation file for a given chromosome, this method will return an empty DataFrame as it is not
    possible for duplicates to exist.

    :param vep_id: The dxid of the VEP annotation file.
    :return: An optional pandas DataFrame of the VEP annotation. If vep_id is None, return None.
    """
    vep_header = ['CHROM', 'POS', 'REF', 'ALT', 'ID', 'FILTER', 'AF', 'F_MISSING', 'AN', 'AC', 'MANE',
                  'ENST', 'ENSG', 'BIOTYPE', 'SYMBOL', 'CSQ', 'gnomAD_AF', 'CADD', 'REVEL', 'SIFT', 'POLYPHEN',
                  'LOFTEE', 'AA', 'AApos', 'PARSED_CSQ', 'MULTI', 'INDEL', 'MINOR', 'MAJOR', 'MAF', 'MAC']

    if vep_id is None:
        return pd.DataFrame(columns = vep_header)
    else:

        # We stream the .vep annotation from dnanexus as it is faster for smaller files like this, and ensures that we don't
        # generate race conditions when loading the same file twice (which we do to check for duplicates...)
        #
        # Note to future devs – DO NOT remove gzip eventhough pandas can direct read gzip. It is not compatible with
        # dxpy.open_dxfile and will error out.
        current_df = pd.read_csv(gzip.open(dxpy.open_dxfile(vep_id, mode='rb'), mode='rt'), sep="\t", header=None, names=vep_header)

        return current_df


# Helper function for process_vep() to remove duplicate variants
def remove_vep_duplicates(current_vep_df: pd.DataFrame, previous_vep_df: pd.DataFrame,
                            vcf_prefix: Path) -> pd.DataFrame:
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
    :param vcf_prefix: A Path object pointing to the VCF prefix of the file to deduplicate.
    :return: A pandas DataFrame of the removed duplicate variants.
    """

    # We have two types of duplicate variants...
    removed_variants = []

    # 1. Variants that are duplicates INTERNAL to this file
    var_counts = current_vep_df['ID'].value_counts()
    duplicates = var_counts[var_counts == 2]

    # This iterates over all duplicates in this file...
    # [0] = the variant ID
    # [1] = the number of times the variant ID appears
    for dup in duplicates.items():
        dup_index = current_vep_df[current_vep_df['ID'] == dup[0]].index.to_list()
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
    duplicates = current_vep_df[current_vep_df['ID'].isin(previous_vep_df['ID'].to_list())]
    current_vep_df = current_vep_df[current_vep_df.index.isin(duplicates.index.to_list()) == False]

    # Merge the duplicate variants from both mode 1. & 2.:
    removed_variants = pd.concat([removed_variants, duplicates])

    # And write the final vep file
    current_vep_df.to_csv(path_or_buf=f'{vcf_prefix}.vep.tsv', na_rep='NA', index=False, sep="\t")

    return removed_variants


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
        og_id = table['ID']

        # Generate a query that matches on position, ref, alt, and og_id. I am fairly certain (I have checked multiple
        # files) that all duplicates come from different multi-allelics, so this _should_ be safe...
        query = f'(POS={position} && REF="{ref}" && ALT="{alt}" && ID="{og_id}")'
        query_builder.append(query)

    query = ' || '.join(query_builder)

    return query


# Helper function for make_bgen_from_vcf() to remove sites identified as duplicate in process_vep() from the bcf
def remove_bcf_duplicates(query_string: str, vcf_prefix: Path) -> None:
    """Remove duplicate variants using bcftools view.

    This method is just a wrapper for the bcftools view command to remove duplicate variants from the BCF file. Following
    deduplication, the original BCF file is deleted and the deduplicated BCF is renamed to the original BCF name.

    :param query_string: The query string to pass to bcftools view from :func:`build_query_string`.
    :param vcf_prefix: A Path object pointing to the VCF prefix of the file to deduplicate.
    """

    cmd = f'bcftools view --threads 2 -e \'{query_string}\' -Ob -o /test/{vcf_prefix}.deduped.bcf /test/{vcf_prefix}.bcf'
    CMD_EXEC.run_cmd_on_docker(cmd)

    # And rename the deduped bcf to match the final one we want to output
    old_vcf = Path(f'{vcf_prefix}.bcf')
    old_vcf.unlink()
    Path(f'{vcf_prefix}.deduped.bcf').rename(old_vcf)
