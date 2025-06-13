import json
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from scripts.chunking_helper import split_coordinates_file

# a bit of setting up
test_data_dir = Path(__file__).parent / 'test_data'
input_file = test_data_dir / 'test_coords_v2.txt'
with open(test_data_dir / 'final_dict.json', 'r') as f:
    gene_dict_local = json.load(f)


@pytest.mark.parametrize("input_data, chunk_size, expected", [
    # Test case 1: two regions on chr7 with one fitting in chunk1, the other in chunk2
    (
            pd.DataFrame({
                'chrom': ['chr7', 'chr7'],
                'start': [36432507, 100679512],
                'end': [36442739, 100694238],
                'vcf_prefix': ['test_input2', 'test_input1'],
                'output_bcf': ['test_input2.vcf.filtered.bcf', 'test_input1.vcf.filtered.bcf'],
                'output_bcf_idx': ['test_input2.vcf.filtered.bcf.csi', 'test_input1.vcf.filtered.bcf.csi'],
                'output_vep': ['test_input2.vcf.vep.tsv.gz', 'test_input1.vcf.vep.tsv.gz'],
                'output_vep_idx': ['test_input2.vcf.vep.tsv.gz.tbi', 'test_input1.vcf.vep.tsv.gz.tbi'],
            }),
            3,
            [
                {'chrom': 'chr7_chunk1', 'start': 36432507, 'end': 36442739,
                 'chunk_start': 36432507, 'chunk_end': 39519370, 'vcf_prefix': 'test_input2'},
                {'chrom': 'chr7_chunk2', 'start': 100679512, 'end': 100694238,
                 'chunk_start': 100679512, 'chunk_end': 104057999, 'vcf_prefix': 'test_input1'},
            ]
    ),
    # Test case 2: single record spans entire chunk
    (
            pd.DataFrame({
                'chrom': ['chr1'],
                'start': [1000000],
                'end': [3000000],
                'vcf_prefix': ['sample1'],
                'output_bcf': ['sample1.bcf'],
                'output_bcf_idx': ['sample1.bcf.csi'],
                'output_vep': ['sample1.vep.tsv.gz'],
                'output_vep_idx': ['sample1.vep.tsv.gz.tbi'],
            }),
            2.5,
            [
                {'chrom': 'chr1_chunk1', 'start': 1000000, 'end': 3000000,
                 'chunk_start': 1000000, 'chunk_end': 3618261, 'vcf_prefix': 'sample1'},
            ]
    ),
    # Test case 3: record overlaps multiple chunks
    (
            pd.DataFrame({
                'chrom': ['chr2', 'chr3'],
                'start': [1000000, 7500000],
                'end': [7000000, 10000000],
                'vcf_prefix': ['sample2', 'sample3'],
                'output_bcf': ['sample2.bcf', 'sample3.bcf'],
                'output_bcf_idx': ['sample2.bcf.csi', 'sample3.bcf.csi'],
                'output_vep': ['sample2.vep.tsv.gz', 'sample3.vep.tsv.gz'],
                'output_vep_idx': ['sample2.vep.tsv.gz.tbi', 'sample3.vep.tsv.gz.tbi'],
            }),
            3,
            [
                {'chrom': 'chr2_chunk1', 'start': 1000000, 'end': 7000000,
                 'chunk_start': 1000000, 'chunk_end': 7104756, 'vcf_prefix': 'sample2'},
                {'chrom': 'chr3_chunk1', 'start': 7500000, 'end': 10000000,
                 'chunk_start': 7500000, 'chunk_end': 10733678, 'vcf_prefix': 'sample3'},
            ]
    ),
    # Additional test cases:
    (
            pd.DataFrame({
                'chrom': ['chr1', 'chr1', 'chr2'],
                'start': [1000000, 3010000, 6010000],
                'end': [3000000, 6000000, 9000000],
                'vcf_prefix': ['sample1', 'sample1', 'sample2'],
                'output_bcf': ['sample1.bcf', 'sample1.bcf', 'sample2.bcf'],
                'output_bcf_idx': ['sample1.bcf.csi', 'sample1.bcf.csi', 'sample2.bcf.csi'],
                'output_vep': ['sample1.vep.tsv.gz', 'sample1.vep.tsv.gz', 'sample2.vep.tsv.gz'],
                'output_vep_idx': ['sample1.vep.tsv.gz.tbi', 'sample1.vep.tsv.gz.tbi', 'sample2.vep.tsv.gz.tbi'],
            }),
            3,
            [
                {'chrom': 'chr1_chunk1', 'start': 1000000, 'end': 3000000,
                 'chunk_start': 1000000, 'chunk_end': 6101489, 'vcf_prefix': 'sample1'},
                {'chrom': 'chr1_chunk1', 'start': 3010000, 'end': 6000000,
                 'chunk_start': 1000000, 'chunk_end': 6101489, 'vcf_prefix': 'sample1'},
                {'chrom': 'chr2_chunk1', 'start': 6010000, 'end': 9000000,
                 'chunk_start': 6010000, 'chunk_end': 9010000, 'vcf_prefix': 'sample2'},
            ]
    ),
])
def test_split_coordinates_pytest(input_data: pd.DataFrame, chunk_size: int, expected: pd.DataFrame):
    """
    Tests for the chunking method in BGEN files, using mark parametrize as our input data
    """

    # run the function
    result = split_coordinates_file(coordinates_file=input_data, gene_dict=gene_dict_local, chunk_size=chunk_size)

    # just for display purposes, so we can see the whole df
    pd.set_option('display.max_columns', None)
    print(result)

    # subset the data to test the columns that we need
    expected_df = pd.DataFrame(expected)
    result_subset = result[['chrom', 'start', 'end', 'chunk_start', 'chunk_end', 'vcf_prefix']].reset_index(drop=True)

    # run the pytest
    pd.testing.assert_frame_equal(result_subset, expected_df)


@pytest.mark.parametrize(
    argnames=['input_data', 'gene_dict', 'chunk_size', 'output_path', 'expected_chunks'],
    argvalues=[
        (Path("test_coords_v2.txt"), Path('final_dict.json'), 3, 'chunked_files', 35),
    ]
)
def test_command_line(input_data: Path, gene_dict: Path, chunk_size: int, output_path: str,
                      expected_chunks: int) -> None:
    """
    Test the command line interface of the script.
    """
    # Construct paths dynamically
    coordinate_path = test_data_dir / input_data
    gene_dict_path = test_data_dir / gene_dict
    output_path = test_data_dir / output_path

    try:
        result = subprocess.run(
            [
                "python",
                str(Path(__file__).parent.parent / "makebgen/process_bgen/helper.py"),
                "--coordinate_path", str(coordinate_path),
                "--gene_dict", str(gene_dict_path),
                "--chunk_size", str(chunk_size),
                "--output_path", str(output_path)
            ],
            capture_output=True,
            text=True,
            check=True
        )
        print(result.stdout)
        print(result.stderr)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with return code {e.returncode}")
        print(f"Error output: {e.stderr}")

    # Check if the output directory exists
    output_dir = Path('chunked_files')
    assert output_dir.exists(), f"Output directory {output_dir} does not exist."

    # Check there are 70 files in the output directory
    output_files = list(output_dir.glob("*"))
    assert len(
        output_files) == expected_chunks, f"Expected 35 files, but found {len(output_files)} files in {output_dir}."
