"""
If you want to keep the output of these tests, change the flag KEEP_TEMP to True.
"""
import csv
import filecmp
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pysam
import pytest
from general_utilities.association_resources import check_gzipped, replace_multi_suffix
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor

from makebgen.deduplication.deduplication import (deduplicate_variants, remove_vep_duplicates,
                                                  remove_bcf_duplicates, build_query_string)

test_data_dir = Path(__file__).parent / 'test_data'

# Set this flag to True if you want to keep (copy) the temporary output files
KEEP_TEMP = False


@pytest.fixture
def temporary_path(tmp_path, monkeypatch):
    """
    Prepare a temporary working directory that contains a copy of the test_data
    directory, then change the working directory to it.

    If KEEP_TEMP is True, after the test the entire temporary directory will be copied
    to a folder 'temp_test_outputs' in the project root.
    """
    # Determine where the original test_data directory is located.
    # (Assumes it is at <project_root>/test_data)
    test_data_source = Path(__file__).parent / "test_data"

    # Create the destination folder inside the tmp_path.
    destination = tmp_path / "test_data"
    destination.parent.mkdir(parents=True, exist_ok=True)

    # Copy the entire test_data directory into the temporary directory.
    shutil.copytree(test_data_source, destination)

    # Change the current working directory to the temporary directory.
    monkeypatch.chdir(tmp_path)

    # Yield the temporary directory to the test.
    yield tmp_path

    # After the test, if KEEP_TEMP is True, copy the temporary directory to a persistent location.
    if KEEP_TEMP:
        persistent_dir = Path(__file__).parent / "temp_test_outputs" / tmp_path.name
        persistent_dir.parent.mkdir(exist_ok=True)
        shutil.copytree(tmp_path, persistent_dir, dirs_exist_ok=True)
        print(f"Temporary output files have been copied to: {persistent_dir}")


@pytest.mark.parametrize(
    "coordinate_path",
    [
        Path('test_coords.txt'),
    ]
)
def test_deduplicate_variants(temporary_path, coordinate_path):
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
    assert new_coords.exists(), f"File was not created: {new_coords}"

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

            deduplicate_variants(
                vep_id=Path(row['output_vep']),
                previous_vep_id=previous_vep_id,
                vcf_path=vcf_path,
                vcf_prefix=vcf_prefix
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
def test_remove_vep_duplicates(temporary_path, current_vep_df, previous_vep_df, expected_deduped_df,
                               expected_removed_df):
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
def test_build_query_string(temporary_path, removed_df, expected_query):
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
    "query_string, vcf_file, vcf_prefix, vcf_path, query_df, should_fail",
    [
        (
                '(POS=100679512 && REF="A" && ALT="C")',
                Path("test_data/test_input1.chunk1.vcf.gz"),
                Path('test_input1.chunk1.vcf.filtered'),
                Path('test_data/test_input1.chunk1.vcf.filtered.bcf'),
                pd.DataFrame({'POS': [100679512], 'REF': ['A'], 'ALT': ['C']}),
                False
        ),
        (
                '(POS=36432507 && REF="C" && ALT="T")',
                Path("test_data/test_input2.chunk2.vcf.gz"),
                Path('test_input2.chunk2.vcf.filtered'),
                Path('test_data/test_input2.chunk2.vcf.filtered.bcf'),
                pd.DataFrame({'POS': [36432507], 'REF': ['C'], 'ALT': ['T']}),
                False
        ),
    ]
)
def test_remove_bcf_duplicates(temporary_path, query_string, vcf_file, vcf_prefix, vcf_path, query_df, should_fail):
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
    # Add test data prefix to the input data
    vcf_path = InputFileHandler(vcf_path).get_file_handle()
    result = remove_bcf_duplicates(query_string, vcf_prefix, vcf_path, cmd_exec)
    # make sure the file exists
    assert result.exists()

    # subset the data from the newly created file
    docker_command = f'docker run -v $(pwd):/test egardner413/mrcepid-burdentesting:latest bcftools query -f "%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\n" test/{result.name} > test_input.txt'
    subprocess.run(docker_command, shell=True, check=True)

    # read this dataframe in and make sure that our query string is no longer there
    test_df = pd.read_csv('test_input.txt', sep='\t', header=None)
    test_df.columns = ['CHR', 'POS', 'ID', 'REF', 'ALT']

    # Select only the relevant columns for comparison
    test_df = test_df[['POS', 'REF', 'ALT']]

    # Check if the combination is present in the DataFrame using df.query()
    is_present = not test_df.query(
        "POS == @query_df['POS'].values[0] and REF == @query_df['REF'].values[0] and ALT == @query_df['ALT'].values[0]"
    ).empty
    # Use assert to ensure the combination is not present
    assert not is_present, "The combination is present in the DataFrame."

    # Open the VCF or BCF file
    vcf = pysam.VariantFile(vcf_file)

    # Count number of variants (records) minus the one we have removed
    num_variants = len(list(vcf)) - 1

    # Assert expected number of variants
    expected_length = len(test_df)
    assert num_variants == expected_length, f"Unexpected number of variants removed! Expected {expected_length}, got {num_variants}"
