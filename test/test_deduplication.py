import csv
import filecmp
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest
from general_utilities.association_resources import check_gzipped, replace_multi_suffix
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor

from makebgen.deduplication.deduplication import (deduplicate_variants, remove_vep_duplicates,
                                                  remove_bcf_duplicates, build_query_string)
from test_process_bgen import delete_test_files

test_data_dir = Path(__file__).parent / 'test_data'


@pytest.mark.parametrize(
    "coordinate_path",
    [
        Path('test_coords.txt'),
    ]
)
def test_deduplicate_variants(coordinate_path):
    """
    Test the deduplicate_variants function to ensure it correctly processes and
    deduplicates variants.

    This test performs the following steps:
    1. Verifies the existence of the coordinate file.
    2. Reads and modifies the coordinate file to update file paths for local testing.
    3. Saves the modified coordinate file.
    4. Reads the modified coordinate file and simulates the deduplication process.
    5. Compares the output files with expected results to ensure correctness.

    :param coordinate_path: Path to the coordinate file used for deduplication.
    """

    # let's deal with the coords file, which needs to
    # be re-formatted to work locally
    # make sure the path to the original coord file exists
    full_path = test_data_dir / coordinate_path
    assert full_path.exists(), f"Path does not exist: {full_path}"

    # read it in so we can change the filepaths
    df = pd.read_csv(full_path, sep='\t')
    columns_to_modify = ['output_bcf', 'output_bcf_idx',
                         'output_vep', 'output_vep_idx']
    # for each column, add the paths of our local files
    for column in columns_to_modify:
        df[column] = str('test_data/') + df[column].astype(str)
    # save this file as 'formatted_coords.txt'
    df.to_csv(test_data_dir / 'formatted_coords.txt', sep='\t')
    # make sure this file exists
    new_coords = test_data_dir / 'formatted_coords.txt'
    assert new_coords

    # now let's read the new coords file into the dict reader
    # here we are simulating the run for make_bgen_from_vcf (which is where the
    # deduplication would occur)
    with check_gzipped(new_coords) as coord_file:
        coord_file_reader = csv.DictReader(coord_file, delimiter="\t")

        previous_vep_id = None
        for row in coord_file_reader:
            # if not using DNA Nexus we can assign a local path to the vcf_path variable
            vcf_path = Path(row['output_bcf'])
            # also, we can change the filepath so that we keep the prefix only
            vcf_prefix = replace_multi_suffix(vcf_path, vcf_path.suffixes[0]).stem

            deduped = deduplicate_variants(
                vep_id=row['output_vep'],
                previous_vep_id=previous_vep_id,
                vcf_prefix=vcf_prefix,
                dna_nexus_run=False
            )

        assert filecmp.cmp('test_input1.vep.tsv',
                           test_data_dir / 'expected_output/test_input1.vep.tsv',
                           shallow=False)
        assert filecmp.cmp('test_input2.vep.tsv',
                           test_data_dir / 'expected_output/test_input2.vep.tsv',
                           shallow=False)


@pytest.mark.parametrize(
    "current_vep_df, previous_vep_df, expected_deduped_df, expected_removed_df",
    [
        (
                pd.DataFrame({
                    'varID': ['var1', 'var1', 'var2', 'var3', 'var4', 'var4'],
                    'FILTER': ['PASS', 'FAIL', 'PASS', 'PASS', 'FAIL', 'PASS'],
                    'F_MISSING': [0.1, 0.2, 0.1, 0.1, 0.3, 0.1],
                    'AF': [0.2, 0.1, 0.3, 0.4, 0.1, 0.2]
                }),
                pd.DataFrame({
                    'varID': ['var3', 'var5'],
                    'FILTER': ['PASS', 'PASS'],
                    'F_MISSING': [0.1, 0.1],
                    'AF': [0.4, 0.5]
                }),
                pd.DataFrame({
                    'varID': ['var1', 'var2', 'var4'],
                    'FILTER': ['PASS', 'PASS', 'PASS'],
                    'F_MISSING': [0.1, 0.1, 0.1],
                    'AF': [0.2, 0.3, 0.2]
                }).reset_index(drop=True),
                pd.DataFrame({
                    'varID': ['var1', 'var4', 'var3'],
                    'FILTER': ['FAIL', 'FAIL', 'PASS'],
                    'F_MISSING': [0.2, 0.3, 0.1],
                    'AF': [0.1, 0.1, 0.4]
                }).reset_index(drop=True)
        )
    ]
)
def test_remove_vep_duplicates(current_vep_df, previous_vep_df, expected_deduped_df, expected_removed_df):
    # Call the function
    deduped_df, removed_df = remove_vep_duplicates(current_vep_df, previous_vep_df)

    # Reset index for comparison
    deduped_df = deduped_df.reset_index(drop=True)
    removed_df = removed_df.reset_index(drop=True)

    # Assert the results
    pd.testing.assert_frame_equal(deduped_df, expected_deduped_df)
    pd.testing.assert_frame_equal(removed_df, expected_removed_df)


@pytest.mark.parametrize(
    "removed_df, expected_query",
    [
        (
                pd.DataFrame({
                    'POS': [12345, 67890],
                    'REF': ['A', 'G'],
                    'ALT': ['T', 'C'],
                    'MA': ['MA1', 'MA2']
                }),
                '(POS=12345 && REF="A" && ALT="T" && MA="MA1") || (POS=67890 && REF="G" && ALT="C" && MA="MA2")'
        ),
        (
                pd.DataFrame({
                    'POS': [11111],
                    'REF': ['C'],
                    'ALT': ['G'],
                    'MA': ['MA3']
                }),
                '(POS=11111 && REF="C" && ALT="G" && MA="MA3")'
        ),
        (
                pd.DataFrame({
                    'POS': [],
                    'REF': [],
                    'ALT': [],
                    'MA': []
                }),
                ''
        )
    ]
)
def test_build_query_string(removed_df, expected_query):
    """
    Test the build_query_string function to ensure it correctly constructs a query string
    (which marks duplicated variants) from the provided DataFrame.

    This test performs the following steps:
    1. Calls the build_query_string function with the provided DataFrame.
    2. Compares the generated query string with the expected query string to ensure correctness.

    :param removed_df: DataFrame containing the variants to be included in the query string.
    :param expected_query: The expected query string.
    """

    query = build_query_string(removed_df)
    assert query == expected_query


@pytest.mark.parametrize(
    "query_string, vcf_file, vcf_prefix, query_df",
    [
        (
                '(POS=100679512 && REF="A" && ALT="C")',
                Path("test_data/test_input1.vcf.filtered.bcf"),
                Path('test_input1.vcf.filtered'),
                pd.DataFrame({
                    'POS': [100679512],
                    'REF': ['A'],
                    'ALT': ['C'],
                }),
        ),
        (
                '(POS=36432507 && REF="C" && ALT="T")',
                Path("test_data/test_input2.vcf.filtered.bcf"),
                Path('test_input2.vcf.filtered'),
                pd.DataFrame({
                    'POS': [36432507],
                    'REF': ['C'],
                    'ALT': ['T'],
                }),
        ),
    ]
)
def test_remove_bcf_duplicates(query_string, vcf_file, vcf_prefix, query_df):
    """
    Test the remove_bcf_duplicates function to ensure it correctly removes duplicate variants from a BCF file.
    This test is a little convoluted, but we want to make sure that we are deleting the duplicates from out files.

    This test performs the following steps:
    1. Copies the source BCF file to the current working directory.
    2. Creates a DockerMount and CommandExecutor for running the bcftools command in a Docker container.
    3. Calls the remove_bcf_duplicates function with the provided query string, VCF prefix, and command executor.
    4. Verifies that the resulting BCF file exists.
    5. Extracts the relevant data from the resulting BCF file and checks that the specified duplicate variants are not present.

    :param query_string: The query string to pass to bcftools view for removing duplicates.
    :param vcf_file: Path to the source BCF file.
    :param vcf_prefix: Path prefix for the VCF file.
    :param query_df: DataFrame containing the expected variants to be removed.
    """

    source_file = vcf_file
    destination_dir = Path('.')
    shutil.copy(source_file, destination_dir)

    # we need to create a docker image locally
    test_mount = DockerMount(Path(os.getcwd()), Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    # Call the function
    result = remove_bcf_duplicates(query_string, vcf_prefix, cmd_exec)
    # make sure the file exists
    assert result.exists()

    # subset the data from the newly created file
    docker_command = f'docker run -v $(pwd):/test egardner413/mrcepid-burdentesting:latest bcftools view test/{result} | grep -v "^#" | cut -f1-5 > test_input.txt'
    subprocess.run(docker_command, shell=True, check=True)

    # read this dataframe in and make sure that our query string is no longer there
    test_df = pd.read_csv('test_input.txt', sep='\t', header=None)
    test_df.columns = ['CHR', 'POS', 'ID', 'REF', 'ALT']

    # Select only the relevant columns for comparison
    test_df = test_df[['POS', 'REF', 'ALT']]

    # Check if the combination is present in the DataFrame
    is_present = ((test_df['POS'] == query_df['POS'].values[0]) &
                  (test_df['REF'] == query_df['REF'].values[0]) &
                  (test_df['ALT'] == query_df['ALT'].values[0])).any()
    # Use assert to ensure the combination is not present
    assert not is_present, "The combination is present in the DataFrame."

    # delete the test files if we don't need them
    delete_test_files(test_data_dir.parent)
