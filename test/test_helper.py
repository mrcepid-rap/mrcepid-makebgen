import json
from pathlib import Path

import pandas as pd
import pytest

from makebgen.process_bgen.helper import split_coordinates_file

test_data_dir = Path(__file__).parent / 'test_data'

with open(test_data_dir / 'final_dict.json', 'r') as f:
    gene_dict = json.load(f)


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
            24257005,
            [
                {'chrom': 'chr7_chunk1', 'start': 36432507, 'end': 60689512, 'vcf_prefix': 'test_input2'},
                {'chrom': 'chr7_chunk3', 'start': 84946519, 'end': 100694238, 'vcf_prefix': 'test_input1'},
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
                {'chrom': 'chr1_chunk1', 'start': 100, 'end': 300, 'vcf_prefix': 'sample1'},
            ]
    ),
    # Test case 3: record overlaps multiple chunks
    (
            pd.DataFrame({
                'chrom': ['chr2'],
                'start': [100],
                'end': [700],
                'vcf_prefix': ['sample2'],
                'output_bcf': ['sample2.bcf'],
                'output_bcf_idx': ['sample2.bcf.csi'],
                'output_vep': ['sample2.vep.tsv.gz'],
                'output_vep_idx': ['sample2.vep.tsv.gz.tbi'],
            }),
            300,
            [
                {'chrom': 'chr2_chunk1', 'start': 100, 'end': 400, 'vcf_prefix': 'sample2'},
                {'chrom': 'chr2_chunk2', 'start': 401, 'end': 700, 'vcf_prefix': 'sample2'},
            ]
    ),
    # Additional tests
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
            250,
            [
                {'chrom': 'chr1_chunk1', 'start': 100, 'end': 350, 'vcf_prefix': 'sample1'},
                {'chrom': 'chr1_chunk2', 'start': 351, 'end': 600, 'vcf_prefix': 'sample1'},
                {'chrom': 'chr2_chunk1', 'start': 601, 'end': 851, 'vcf_prefix': 'sample2'},
                {'chrom': 'chr2_chunk2', 'start': 852, 'end': 900, 'vcf_prefix': 'sample2'},
            ]
    ),
])
def test_split_coordinates_file(input_data, chunk_size, expected):
    result = split_coordinates_file(coordinates_file=input_data, gene_dict=gene_dict, chunk_size=chunk_size)

    pd.set_option('display.max_columns', None)

    print(result)

    expected_df = pd.DataFrame(expected)
    result_subset = result[['chrom', 'start', 'end', 'vcf_prefix']].reset_index(drop=True)

    pd.testing.assert_frame_equal(result_subset, expected_df)
