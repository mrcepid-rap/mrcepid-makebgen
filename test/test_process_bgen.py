"""
Note: these tests are designed to be run sequentially (first test, then second, then third etc.)
This is due to the fact that the functions are run sequentially, but also so that we can look at
the output files as they get generated. To run them all, please run tests for 'Current File'.
"""

import csv
import filecmp
import glob
import gzip
import os
import shutil
from pathlib import Path

import pandas as pd
import pytest
from general_utilities.association_resources import check_gzipped
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor

from makebgen.process_bgen.process_bgen import make_bgen_from_vcf, correct_sample_file, make_final_bgen

test_data_dir = Path(__file__).parent / 'test_data'


@pytest.mark.parametrize(
    "coordinate_path",
    [
        Path('test_coords.txt'),
    ]
)
def test_make_bgen_from_vcf(coordinate_path, make_bcf=False):
    """
    Test the `make_bgen_from_vcf` function.

    This test function reads a coordinate file, reformats the file paths to work locally,
    and then processes each row to create BGEN files from VCF files. It also verifies
    that the generated BGEN files match the expected output.

    :param coordinate_path: Path to the coordinate file containing the coordinates of all BCF files to be processed.
    :param make_bcf: Boolean flag indicating whether a concatenated BCF should be made in addition to the BGEN.
    """

    total_bcf = 0

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

    # we need to create a docker image locally
    test_mount = DockerMount(Path(os.getcwd()), Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    # now let's read the new coords file into the dict reader
    with check_gzipped(new_coords) as coord_file:
        coord_file_reader = csv.DictReader(coord_file, delimiter="\t")

        previous_vep_id = None
        for row in coord_file_reader:
            total_bcf += 1
            make_bgen_from_vcf(
                vcf_id=row['output_bcf'],
                vep_id=row['output_vep'],
                previous_vep_id=previous_vep_id,
                start=row['start'],
                make_bcf=make_bcf,
                cmd_exec=cmd_exec,
                dna_nexus_run=False)
            previous_vep_id = row['output_vep']

        assert filecmp.cmp('test_input1.bgen',
                           test_data_dir / 'expected_output/test_input1.bgen',
                           shallow=False)
        assert filecmp.cmp('test_input2.bgen',
                           test_data_dir / 'expected_output/test_input2.bgen',
                           shallow=False)


@pytest.mark.parametrize(
    "sample_file",
    [
        Path('test_input1.sample'),
        Path('test_input2.sample'),
    ]
)
def test_correct_sample_file(sample_file, tmp_path):
    """
    Test the `correct_sample_file` function.

    This test function takes a sample file and a temporary path, corrects the sample file,
    and verifies that the corrected sample file exists and contains the expected content.

    :param sample_file: Path to the sample file to be corrected.
    :param tmp_path: Temporary path for storing the corrected sample file.
    """

    output_prefix = tmp_path / "test_output"
    corrected_sample = correct_sample_file(sample_file, output_prefix)
    assert corrected_sample.exists()

    with corrected_sample.open('r') as file:
        lines = file.readlines()

    assert lines[0].strip() == 'ID_1 ID_2 missing sex'
    assert lines[1].strip() == '0 0 0 D'
    assert lines[2].strip() == 'HG00096 HG00096 0 NA'
    assert lines[3].strip() == 'HG00097 HG00097 0 NA'
    assert lines[4].strip() == 'HG00099 HG00099 0 NA'


@pytest.mark.parametrize(
    "bgen_prefixes, output_prefix, make_bcf",
    [
        ({'test_input1', 'test_input2'},
         'test_bgen',
         False)
    ]
)
def test_make_final_bgen(bgen_prefixes, output_prefix, make_bcf):
    """
    Test the `make_final_bgen` function.

    This test function sets up the necessary environment, including creating a Docker image locally,
    copying required files, and then calls the `make_final_bgen` function to generate the final BGEN files.
    It verifies that the generated BGEN files match the expected output.

    :param bgen_prefixes: A set of prefixes for the BGEN files to be processed.
    :param output_prefix: The prefix for the output BGEN files.
    :param make_bcf: Boolean flag indicating whether a concatenated BCF should be made in addition to the BGEN.
    """

    # we need to create a docker image locally
    test_mount = DockerMount(Path(os.getcwd()), Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    # for this to work we need to copy some files into the current working directory

    for file in glob.glob(str(test_data_dir / "*.vep.tsv.gz")):
        new_filename = os.path.basename(file).replace('.vcf', '').replace('.gz', '')
        shutil.copy(file, new_filename + '.gz')

        with gzip.open(new_filename + '.gz', 'rb') as f_in:
            with open(new_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    final_bgens = make_final_bgen(bgen_prefixes, output_prefix, make_bcf,
                                  dna_nexus_run=False, cmd_exec=cmd_exec)

    assert filecmp.cmp(final_bgens['bgen']['file'],
                       test_data_dir / 'expected_output/test_bgen.bgen',
                       shallow=False)
    assert filecmp.cmp(final_bgens['bgen']['sample'],
                       test_data_dir / 'expected_output/test_bgen.sample',
                       shallow=False)

    delete_test_files(test_data_dir.parent)


def delete_test_files(directory):
    """
    Delete all the files after we are done testing them
    """
    # Use glob to find files starting with "test_input" or "test_bgen"
    files_to_delete = glob.glob(os.path.join(directory, 'test_input*')) + \
                      glob.glob(os.path.join(directory, 'test_bgen*'))

    # Iterate and delete each file
    for file_path in files_to_delete:
        try:
            os.remove(file_path)
        except Exception as e:
            print(f"Error deleting {file_path}: {e}")

    print("Test output files have been deleted")
