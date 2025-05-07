import json
from pathlib import Path

import pandas as pd
import pytest

from makebgen.process_bgen.helper import split_coordinates_file, run_splitter

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
                 'chunk_start': 100679512, 'chunk_end': 100694238, 'vcf_prefix': 'test_input1'},
            ]
    ),
    # Test case 2: single record spans entire chunk
    (
            pd.DataFrame({
                'chrom': ['chr1'],
                'start': [100],
                'end': [300],
                'vcf_prefix': ['sample1'],
                'output_bcf': ['sample1.bcf'],
                'output_bcf_idx': ['sample1.bcf.csi'],
                'output_vep': ['sample1.vep.tsv.gz'],
                'output_vep_idx': ['sample1.vep.tsv.gz.tbi'],
            }),
            250,
            [
                {'chrom': 'chr1_chunk1', 'start': 100, 'end': 300,
                 'chunk_start': 100, 'chunk_end': 300, 'vcf_prefix': 'sample1'},
            ]
    ),
    # Test case 3: record overlaps multiple chunks
    (
            pd.DataFrame({
                'chrom': ['chr2', 'chr3'],
                'start': [100, 750],
                'end': [700, 1000],
                'vcf_prefix': ['sample2', 'sample3'],
                'output_bcf': ['sample2.bcf', 'sample3.bcf'],
                'output_bcf_idx': ['sample2.bcf.csi', 'sample3.bcf.csi'],
                'output_vep': ['sample2.vep.tsv.gz', 'sample3.vep.tsv.gz'],
                'output_vep_idx': ['sample2.vep.tsv.gz.tbi', 'sample3.vep.tsv.gz.tbi'],
            }),
            300,
            [
                {'chrom': 'chr2_chunk1', 'start': 100, 'end': 700,
                 'chunk_start': 100, 'chunk_end': 700, 'vcf_prefix': 'sample2'},
                {'chrom': 'chr3_chunk1', 'start': 750, 'end': 1000,
                 'chunk_start': 750, 'chunk_end': 1000, 'vcf_prefix': 'sample3'},
            ]
    ),
    # Additional test cases:
    (
            pd.DataFrame({
                'chrom': ['chr1', 'chr1', 'chr2'],
                'start': [100, 301, 601],
                'end': [300, 600, 900],
                'vcf_prefix': ['sample1', 'sample1', 'sample2'],
                'output_bcf': ['sample1.bcf', 'sample1.bcf', 'sample2.bcf'],
                'output_bcf_idx': ['sample1.bcf.csi', 'sample1.bcf.csi', 'sample2.bcf.csi'],
                'output_vep': ['sample1.vep.tsv.gz', 'sample1.vep.tsv.gz', 'sample2.vep.tsv.gz'],
                'output_vep_idx': ['sample1.vep.tsv.gz.tbi', 'sample1.vep.tsv.gz.tbi', 'sample2.vep.tsv.gz.tbi'],
            }),
            300,
            [
                {'chrom': 'chr1_chunk1', 'start': 100, 'end': 300,
                 'chunk_start': 100, 'chunk_end': 600, 'vcf_prefix': 'sample1'},
                {'chrom': 'chr1_chunk1', 'start': 301, 'end': 600,
                 'chunk_start': 100, 'chunk_end': 600, 'vcf_prefix': 'sample1'},
                {'chrom': 'chr2_chunk1', 'start': 601, 'end': 900,
                 'chunk_start': 601, 'chunk_end': 900, 'vcf_prefix': 'sample2'},
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
    argnames=['input_data', 'gene_dict', 'chunk_size', 'expected'],
    argvalues=[
        (Path("test_coords.txt"), 'final_dict.json', 3, 'test_coords_with_chunking.csv'),
        (Path("test_coords_v2.txt"), 'final_dict.json', 3, 'test_coords_v2_with_chunking.csv'),
        (Path("formatted_coords.txt"), 'final_dict.json', 3, 'formatted_coords_with_chunking.csv'),
        # (Path("test_coordinates.txt"), 'test/test_data/final_dict.json', 3, 'test_coordinates_with_chunking.csv'),
    ]
)
def test_run_splitter(input_data, gene_dict, chunk_size, expected):
    """
    This is a pytest using more of a real-life example
    """

    data = test_data_dir / input_data
    gene_dict = test_data_dir / gene_dict
    result = run_splitter(coordinate_path=data, gene_dict=gene_dict, chunk_size=chunk_size)

    # show all pandas columns
    pd.set_option('display.max_columns', None)

    # assert dataframes are equal
    actual = pd.read_csv(result, sep='\t')
    expected_df = pd.read_csv(test_data_dir / 'expected_output' / expected, sep="\t")

    pd.testing.assert_frame_equal(actual, expected_df)
