"""
Note: these tests are designed to be run all at once, as the input of a downstream test relies on the output of an
upstream test. The best way to do this is to hit "Current File" in the top right corner, if using PyCharm; or, run
pytest in the directory via CLI. If you want to keep the output of these tests, change the flag KEEP_TEMP to True.
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
from bgen import BgenReader
from general_utilities.association_resources import check_gzipped
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler
from general_utilities.job_management.command_executor import DockerMount, CommandExecutor

from makebgen.process_bgen.process_bgen import make_bgen_from_vcf, correct_sample_file, make_final_bgen

test_data_dir = Path(__file__).parent / 'test_data'


@pytest.fixture(scope="module")
def pipeline_data():
    """
    Provides a shared dictionary for storing intermediate file names
    produced by each stage of the pipeline.
    """
    return {}


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


@pytest.fixture
def formatted_coords_file(coordinate_path: Path, test_data_directory: Path = test_data_dir) -> Path:
    """
    Fixture to create a formatted coordinates file for testing.
    """
    full_path = test_data_directory / coordinate_path
    # Ensure the coordinate file exists
    assert full_path.exists(), f"Path does not exist: {full_path}"

    # Read the coordinate file into a DataFrame
    df = pd.read_csv(full_path, sep='\t')

    # Modify specific columns to include the 'test_data/' prefix
    columns_to_modify = ['output_bcf', 'output_bcf_idx', 'output_vep', 'output_vep_idx']
    df[columns_to_modify] = df[columns_to_modify].apply(lambda col: col.apply(lambda x: f'test_data/{x}'))

    # Save the modified DataFrame to a new file
    formatted_coords = test_data_dir / 'formatted_coords.txt'
    df.to_csv(formatted_coords, sep='\t')

    # Ensure the new formatted coordinates file was created
    assert formatted_coords.exists(), "Formatted coordinates file was not created"

    return formatted_coords


@pytest.mark.parametrize(
    "coordinate_path",
    [
        Path('test_coords.txt'),
    ]
)
def test_make_bgen_from_vcf(temporary_path, pipeline_data, formatted_coords_file, make_bcf=False):
    """
    Test the `make_bgen_from_vcf` function.

    This test function reads a coordinate file, reformats the file paths to work locally,
    and then processes each row to create BGEN files from VCF files. It also verifies
    that the generated BGEN files match the expected output.

    :param coordinate_path: Path to the coordinate file containing the coordinates of all BCF files to be processed.
    :param make_bcf: Boolean flag indicating whether a concatenated BCF should be made in addition to the BGEN.
    """

    total_bcf = 0

    # Use the formatted coordinates file from the fixture
    new_coords = formatted_coords_file
    input_coords = InputFileHandler(new_coords)

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
                vcf_id=Path(row['output_bcf']),
                vep_id=row['output_vep'],
                previous_vep_id=previous_vep_id,
                start=row['start'],
                make_bcf=make_bcf,
                cmd_exec=cmd_exec
            )
            previous_vep_id = row['output_vep']

        assert filecmp.cmp('test_input1.bgen',
                           test_data_dir / 'expected_output/test_input1.bgen',
                           shallow=False)
        assert filecmp.cmp('test_input2.bgen',
                           test_data_dir / 'expected_output/test_input2.bgen',
                           shallow=False)

    pipeline_data['test_input1.sample'] = temporary_path / 'test_input1.sample'
    pipeline_data['test_input2.sample'] = temporary_path / 'test_input2.sample'
    pipeline_data['test_input1.bgen'] = temporary_path / 'test_input1.bgen'
    pipeline_data['test_input2.bgen'] = temporary_path / 'test_input2.bgen'


@pytest.mark.parametrize(
    "sample_key",
    [
        'test_input1.sample',
        'test_input2.sample',
    ]
)
def test_correct_sample_file(temporary_path: Path, pipeline_data: dict, sample_key: str, tmp_path: Path):
    """
    Test the `correct_sample_file` function.

    This test retrieves a sample file path from pipeline_data (using the provided key),
    corrects the sample file using the correct_sample_file function, and verifies that the
    corrected file exists and contains the expected content.
    """
    # Retrieve the actual sample file from pipeline_data using the key.
    sample_file = pipeline_data.get(sample_key)
    assert sample_file is not None, f"No entry found in pipeline_data for key '{sample_key}'"

    print(f"Testing sample file: {sample_file}")

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
        ({ 'test_input1': 100679512,
           'test_input2': 36432507},
         'test_bgen',
         False)
    ]
)
def test_make_final_bgen(temporary_path, pipeline_data, bgen_prefixes, output_prefix, make_bcf):
    """
    Test the `make_final_bgen` function.

    This test function sets up the necessary environment, including creating a Docker image locally,
    copying required files, and then calls the `make_final_bgen` function to generate the final BGEN files.
    It verifies that the generated BGEN files match the expected output.

    :param bgen_prefixes: A set of prefixes for the BGEN files to be processed.
    :param output_prefix: The prefix for the output BGEN files.
    :param make_bcf: Boolean flag indicating whether a concatenated BCF should be made in addition to the BGEN.
    """

    # first we need to retrieve the files from our pipeline_data object
    try:
        bgen_file1 = pipeline_data['test_input1.bgen']
        bgen_file2 = pipeline_data['test_input2.bgen']
        sample_file1 = pipeline_data['test_input1.sample']
        sample_file2 = pipeline_data['test_input2.sample']
    except KeyError:
        pytest.skip("BGEN files not available in pipeline_data from previous stage")

    # let's copy these into the current directory, otherwise the function won't work
    current_dir = Path(os.getcwd())
    for file in (bgen_file1, bgen_file2, sample_file1, sample_file2):
        target = current_dir / file.name
        if not target.exists():
            shutil.copy(file, target)
            print(f"Copied {file} to {target}")

    # we need to create a docker image locally
    test_mount = DockerMount(Path(os.getcwd()), Path('/test/'))
    cmd_exec = CommandExecutor(docker_image='egardner413/mrcepid-burdentesting', docker_mounts=[test_mount])

    # for this to work we also need to copy some other test files into the current working directory
    for file in glob.glob(str(test_data_dir / "*.vep.tsv.gz")):
        new_filename = os.path.basename(file).replace('.vcf', '').replace('.gz', '')
        shutil.copy(file, new_filename + '.gz')

        with gzip.open(new_filename + '.gz', 'rb') as f_in:
            with open(new_filename, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

    # finally let's run the function
    final_bgens = make_final_bgen(bgen_prefixes, output_prefix, make_bcf,
                                  cmd_exec=cmd_exec)

    # check for number of samples
    with BgenReader(final_bgens['bgen']['file'], delay_parsing=True) as bfile:
        assert len(bfile.samples) == 3202
    # check for number of variants
    with BgenReader(final_bgens['bgen']['file'], delay_parsing=True) as bfile:
        assert len(bfile.positions()) == 1465

    # make sure the output is as expected
    assert filecmp.cmp(final_bgens['bgen']['file'],
                       test_data_dir / 'expected_output/test_bgen.bgen',
                       shallow=False)
    assert filecmp.cmp(final_bgens['bgen']['sample'],
                       test_data_dir / 'expected_output/test_bgen.sample',
                       shallow=False)
