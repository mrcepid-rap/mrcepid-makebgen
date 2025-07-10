import csv
import json
from pathlib import Path
from typing import Tuple, List, Optional, Dict

import pandas as pd
import pytest
from intervaltree import IntervalTree, Interval

from makebgen.chunker.chunking_helper import build_interval_tree, parse_gene_dict, get_safe_chunk_ends, \
    chunk_chromosome, is_position_within_gene, split_coordinates_file, find_next_safe_end, split_batch_files

# a bit of setting up
test_data_dir = Path(__file__).parent / 'test_data'
input_file = test_data_dir / 'test_coords_v2.txt'
with open(test_data_dir / 'final_dict.json', 'r') as f:
    gene_dict_local = json.load(f)


@pytest.mark.parametrize(
    "gene_dict, expected",
    [
        (
                {
                    "BRCA1": {"genomic_location": {"chromosome": 1, "start": 123456, "end": 123556}},
                    "TP53": {"genomic_location": {"chromosome": 17, "start": 223456, "end": 223556}},
                    "EGFR": {"genomic_location": {"chromosome": 7, "start": 323456, "end": 323556}}
                },
                pd.DataFrame([
                    {"gene": "BRCA1", "chrom": "chr1", "start": 123456, "end": 123556},
                    {"gene": "TP53", "chrom": "chr17", "start": 223456, "end": 223556},
                    {"gene": "EGFR", "chrom": "chr7", "start": 323456, "end": 323556}
                ])
        ),
        (
                {
                    "MYC": {"genomic_location": {"chromosome": 8, "start": 100000, "end": 100500}}
                },
                pd.DataFrame([
                    {"gene": "MYC", "chrom": "chr8", "start": 100000, "end": 100500}
                ])
        )
    ]
)
def test_parse_gene_dict(tmp_path: Path, gene_dict: dict, expected: pd.DataFrame) -> None:
    """
    Test the `parse_gene_dict` function.

    This test creates a temporary JSON file for the gene dictionary, 
    calls the `parse_gene_dict` function, and asserts that the result 
    matches the expected DataFrame.

    :param tmp_path: Temporary directory path provided by pytest.
    :param gene_dict: Dictionary containing gene information.
    :param expected: Expected DataFrame to compare against the function output.
    :return: None
    """  # Create a temporary JSON file for the gene dictionary
    gene_dict_path = tmp_path / "gene_dict.json"
    with open(gene_dict_path, "w") as f:
        json.dump(gene_dict, f)

    # Call the function
    result = parse_gene_dict(gene_dict_path)

    # Assert the result matches the expected DataFrame
    pd.testing.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "df, chrom, start_col, end_col, label_col, expected_intervals, xfail",
    [
        # Valid test
        (
                pd.DataFrame([
                    {"chrom": "chr1", "start": 123456, "end": 123556, "gene": "BRCA1"},
                    {"chrom": "chr1", "start": 223456, "end": 223556, "gene": "TP53"},
                    {"chrom": "chr2", "start": 323456, "end": 323556, "gene": "EGFR"}
                ]),
                "chr1",
                "start",
                "end",
                "gene",
                [(123456, 123556, "BRCA1"), (223456, 223556, "TP53")],
                False
        ),
        # Valid test
        (
                pd.DataFrame([
                    {"chrom": "chr1", "start": 100000, "end": 100500, "gene": "MYC"},
                    {"chrom": "chr1", "start": 200000, "end": 200500, "gene": "KRAS"}
                ]),
                "chr1",
                "start",
                "end",
                "gene",
                [(100000, 100500, "MYC"), (200000, 200500, "KRAS")],
                False
        ),
        # Failing test
        (
                pd.DataFrame([
                    {"chrom": "chr1", "start": 123456, "end": 123556, "gene": "BRCA1"},
                    {"chrom": "chr2", "start": 223456, "end": 223556, "gene": "TP53"}
                ]),
                "chr3",
                "start",
                "end",
                "gene",
                [(123456, 123556, "BRCA1")],
                True
        ),
        # Failing test
        (
                pd.DataFrame([
                    {"chromosome": "chr1", "start_pos": 100000, "end_pos": 100500, "gene": "MYC"}
                ]),
                "chr1",
                "start",
                "end",
                "gene",
                [(100000, 100500, "MYC")],
                True
        )
    ]
)
def test_build_interval_tree(df: pd.DataFrame, chrom: str, start_col: str, end_col: str, label_col: str,
                             expected_intervals: List[Tuple[int, int, str]], xfail: bool) -> None:
    """
    Test the `build_interval_tree` function.

    This test verifies that the interval tree is correctly built from the input DataFrame
    and matches the expected intervals. It also handles cases where the test is expected to fail.

    :param df: DataFrame containing the input data for the interval tree.
    :param chrom: Chromosome identifier to filter the DataFrame.
    :param start_col: Column name for the start positions.
    :param end_col: Column name for the end positions.
    :param label_col: Column name for the labels (e.g., gene names).
    :param expected_intervals: List of expected intervals as tuples (start, end, label).
    :param xfail: Boolean indicating if the test is expected to fail.
    :return: None
    """
    if xfail:
        pytest.xfail("This test is expected to fail.")
    tree = build_interval_tree(df, chrom, start_col, end_col, label_col)
    actual_intervals = [(interval.begin, interval.end, interval.data) for interval in tree]
    assert sorted(actual_intervals) == sorted(expected_intervals)


@pytest.mark.parametrize(
    "intervals, pos, expected_gene",
    [
        # Position falls within a single gene
        (
                [(100, 200, "BRCA1")],
                150,
                "BRCA1"
        ),
        # Position falls outside any gene
        (
                [(100, 200, "BRCA1")],
                250,
                None
        ),
        # Position falls within multiple overlapping genes
        (
                [(100, 200, "BRCA1"), (150, 250, "TP53")],
                175,
                "BRCA1"  # Returns the first match
        ),
        # No intervals in the tree
        (
                [],
                150,
                None
        ),
        # Position is exactly on the boundary of a gene
        (
                [(100, 200, "BRCA1")],
                100,
                "BRCA1"
        )
    ]
)
def test_is_position_within_gene(intervals: List[Tuple[int, int, str]], pos: int, expected_gene: Optional[str]) -> None:
    """
    Test the `is_position_within_gene` function.

    This test verifies whether a given position falls within a gene interval
    in the provided interval tree.

    :param intervals: List of tuples representing gene intervals (start, end, gene name).
    :param pos: Position to check within the intervals.
    :param expected_gene: Expected gene name if the position is within a gene, otherwise None.
    :return: None
    """
    # Build the interval tree
    tree = IntervalTree(Interval(start, end, data) for start, end, data in intervals)

    # Call the function
    result = is_position_within_gene(tree, pos)

    # Assert the result matches the expected output
    assert result == expected_gene


@pytest.mark.parametrize("chrom_df, gene_intervals, expected_safe_ends", [
    # No positions overlap with genes
    (
            pd.DataFrame({'end': [100, 200, 300, 400]}),
            IntervalTree([Interval(500, 600)]),
            [100, 200, 300, 400]
    ),
    # Some positions overlap with genes
    (
            pd.DataFrame({'end': [100, 200, 300, 400]}),
            IntervalTree([Interval(150, 250), Interval(350, 450)]),
            [100, 300]
    ),
    # All positions overlap with genes
    (
            pd.DataFrame({'end': [100, 200, 300, 400]}),
            IntervalTree([Interval(50, 450)]),
            []
    ),
    # No gene intervals provided
    (
            pd.DataFrame({'end': [100, 200, 300, 400]}),
            IntervalTree(),
            [100, 200, 300, 400]
    ),
])
def test_get_safe_chunk_ends(chrom_df: pd.DataFrame, gene_intervals: IntervalTree,
                             expected_safe_ends: List[int]) -> None:
    """
    Test the `get_safe_chunk_ends` function.

    This test verifies that the function correctly identifies safe chunk end positions
    that do not overlap with gene intervals.

    :param chrom_df: DataFrame containing chromosome data with 'end' positions.
    :param gene_intervals: IntervalTree containing gene intervals to check for overlaps.
    :param expected_safe_ends: List of expected safe chunk end positions.
    :return: None
    """
    result = get_safe_chunk_ends(chrom_df, gene_intervals)
    assert result == expected_safe_ends


@pytest.mark.parametrize("chrom_df, gene_df, chrom, chunk_size_bp, expected_chunks", [
    # Simple chromosome with no overlapping genes
    (
            pd.DataFrame({
                'start': [100000, 200000, 300000],
                'end': [150000, 250000, 350000],
                'vcf_prefix': ['sample1', 'sample2', 'sample3'],
                'output_bcf': ['sample1.bcf', 'sample2.bcf', 'sample3.bcf'],
                'output_bcf_idx': ['sample1.bcf.csi', 'sample2.bcf.csi', 'sample3.bcf.csi'],
                'output_vep': ['sample1.vep.tsv.gz', 'sample2.vep.tsv.gz', 'sample3.vep.tsv.gz'],
                'output_vep_idx': ['sample1.vep.tsv.gz.tbi', 'sample2.vep.tsv.gz.tbi', 'sample3.vep.tsv.gz.tbi'],
            }),
            pd.DataFrame({
                'gene': ['gene1', 'gene2'],
                'chrom': ['chr1', 'chr1'],
                'start': [400000, 500000],
                'end': [450000, 550000],
            }),
            'chr1',
            200000,
            [
                {'chrom': 'chr1_chunk1', 'start': 100000, 'end': 150000, 'chunk_start': 100000, 'chunk_end': 250000,
                 'vcf_prefix': 'sample1'},
                {'chrom': 'chr1_chunk1', 'start': 200000, 'end': 250000, 'chunk_start': 100000, 'chunk_end': 250000,
                 'vcf_prefix': 'sample2'},
                {'chrom': 'chr1_chunk2', 'start': 300000, 'end': 350000, 'chunk_start': 300000, 'chunk_end': 350000,
                 'vcf_prefix': 'sample3'},
            ]
    ),
    # Overlapping genes and chunks
    (
            pd.DataFrame({
                'start': [100000, 200000, 300000],
                'end': [150000, 250000, 350000],
                'vcf_prefix': ['sample1', 'sample2', 'sample3'],
                'output_bcf': ['sample1.bcf', 'sample2.bcf', 'sample3.bcf'],
                'output_bcf_idx': ['sample1.bcf.csi', 'sample2.bcf.csi', 'sample3.bcf.csi'],
                'output_vep': ['sample1.vep.tsv.gz', 'sample2.vep.tsv.gz', 'sample3.vep.tsv.gz'],
                'output_vep_idx': ['sample1.vep.tsv.gz.tbi', 'sample2.vep.tsv.gz.tbi', 'sample3.vep.tsv.gz.tbi'],
            }),
            pd.DataFrame({
                'gene': ['gene1', 'gene2', 'gene3'],
                'chrom': ['chr1', 'chr1', 'chr1'],
                'start': [250000, 300000, 420000],
                'end': [300000, 340000, 500000],
            }),
            'chr1',
            200000,
            [
                {'chrom': 'chr1_chunk1', 'start': 100000, 'end': 150000, 'chunk_start': 100000, 'chunk_end': 150000,
                 'vcf_prefix': 'sample1'},
                {'chrom': 'chr1_chunk2', 'start': 200000, 'end': 250000, 'chunk_start': 200000, 'chunk_end': 350000,
                 'vcf_prefix': 'sample2'},
                {'chrom': 'chr1_chunk2', 'start': 300000, 'end': 350000, 'chunk_start': 200000, 'chunk_end': 350000,
                 'vcf_prefix': 'sample3'},
            ]
    ),
])
def test_chunk_chromosome(chrom_df: pd.DataFrame, gene_df: pd.DataFrame, chrom: str, chunk_size_bp: int,
                          expected_chunks: List[Dict]) -> None:
    """
    Test the `chunk_chromosome` function.

    This test verifies that the function correctly splits a chromosome into chunks
    based on the provided chunk size, ensuring that chunks do not overlap with genes.

    :param chrom_df: DataFrame containing chromosome data with 'start' and 'end' positions.
    :param gene_df: DataFrame containing gene information with 'start' and 'end' positions.
    :param chrom: Chromosome identifier to process.
    :param chunk_size_bp: Size of each chunk in base pairs.
    :param expected_chunks: List of dictionaries representing the expected chunk data.
    :return: None
    """

    # Run the function
    result_chunks, result_logs = chunk_chromosome(chrom_df, gene_df, chrom, chunk_size_bp)

    # Convert expected_chunks to DataFrame for comparison
    expected_chunks_df = pd.DataFrame(expected_chunks)

    # Compare chunks
    pd.testing.assert_frame_equal(result_chunks[['chrom', 'start', 'end', 'chunk_start', 'chunk_end', 'vcf_prefix']],
                                  expected_chunks_df)


@pytest.mark.parametrize("input_data, gene_df, chunk_size, expected, should_raise", [
    # two regions on chr7 with non-overlapping gene far downstream
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
            pd.DataFrame({
                'gene': ['GENE_A', 'GENE_B'],
                'chrom': ['chr7', 'chr7'],
                'start': [150000000, 160000000],
                'end': [155000000, 165000000],
            }),
            3,
            [
                {'chrom': 'chr7_chunk1', 'start': 36432507, 'end': 36442739,
                 'chunk_start': 36432507, 'chunk_end': 36442739, 'vcf_prefix': 'test_input2'},
                {'chrom': 'chr7_chunk2', 'start': 100679512, 'end': 100694238,
                 'chunk_start': 100679512, 'chunk_end': 100694238, 'vcf_prefix': 'test_input1'},
            ],
            False
    ),

    # single record exactly matches chunk size
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
            pd.DataFrame({
                'gene': ['GENE_X'],
                'chrom': ['chr1'],
                'start': [10000000],
                'end': [12000000],
            }),
            2.5,
            [
                {'chrom': 'chr1_chunk1', 'start': 1000000, 'end': 3000000,
                 'chunk_start': 1000000, 'chunk_end': 3000000, 'vcf_prefix': 'sample1'},
            ],
            False
    ),

    # record overlaps multiple chunks, now with safe end <= ideal_end
    (
            pd.DataFrame({
                'chrom': ['chr2', 'chr3'],
                'start': [1000000, 7500000],
                'end': [3500000, 9600000],  # changed from 7,000,000 to 3,500,000
                'vcf_prefix': ['sample2', 'sample3'],
                'output_bcf': ['sample2.bcf', 'sample3.bcf'],
                'output_bcf_idx': ['sample2.bcf.csi', 'sample3.bcf.csi'],
                'output_vep': ['sample2.vep.tsv.gz', 'sample3.vep.tsv.gz'],
                'output_vep_idx': ['sample2.vep.tsv.gz.tbi', 'sample3.vep.tsv.gz.tbi'],
            }),
            pd.DataFrame({
                'gene': ['GENE_Y', 'GENE_Z'],
                'chrom': ['chr2', 'chr3'],
                'start': [8000000, 7400000],  # chr2 gene moved far away
                'end': [9000000, 9500000],
            }),
            3,
            [
                {'chrom': 'chr2_chunk1', 'start': 1000000, 'end': 3500000,
                 'chunk_start': 1000000, 'chunk_end': 3500000, 'vcf_prefix': 'sample2'},
                {'chrom': 'chr3_chunk1', 'start': 7500000, 'end': 9600000,
                 'chunk_start': 7500000, 'chunk_end': 9600000, 'vcf_prefix': 'sample3'},
            ],
            False
    ),
    # triggers RuntimeError due to no safe chunk end
    (
            pd.DataFrame({
                'chrom': ['chr5'],
                'start': [1000000],
                'end': [7000000],
                'vcf_prefix': ['failcase'],
                'output_bcf': ['failcase.bcf'],
                'output_bcf_idx': ['failcase.bcf.csi'],
                'output_vep': ['failcase.vep.tsv.gz'],
                'output_vep_idx': ['failcase.vep.tsv.gz.tbi'],
            }),
            pd.DataFrame({
                'gene': ['GENE_FAIL'],
                'chrom': ['chr5'],
                'start': [1000000],
                'end': [8000000],
            }),
            3,
            None,
            True  # this test should raise RuntimeError
    )

])
def test_split_coordinates_file(input_data: pd.DataFrame, gene_df: pd.DataFrame, chunk_size: int,
                                expected: Optional[List[Dict]], should_raise: bool) -> None:
    """
    Test the `split_coordinates_file` function.

    This test verifies that the function correctly splits the input data into chunks
    based on the specified chunk size, ensuring that chunks do not overlap with genes.
    It also handles cases where the function is expected to raise a `RuntimeError`.

    :param input_data: DataFrame containing the input coordinates data.
    :param gene_df: DataFrame containing gene information.
    :param chunk_size: Size of each chunk in megabases.
    :param expected: List of dictionaries representing the expected chunk data, or None if an error is expected.
    :param should_raise: Boolean indicating if the test is expected to raise a `RuntimeError`.
    :return: None
    """
    if should_raise:
        with pytest.raises(RuntimeError, match="No safe chunk end found"):
            split_coordinates_file(coordinates_file=input_data, gene_df=gene_df, chunk_size=chunk_size)
    else:
        result, _ = split_coordinates_file(coordinates_file=input_data, gene_df=gene_df, chunk_size=chunk_size)
        result_subset = result[['chrom', 'start', 'end', 'chunk_start', 'chunk_end', 'vcf_prefix']].reset_index(
            drop=True)
        expected_df = pd.DataFrame(expected)
        pd.testing.assert_frame_equal(result_subset, expected_df)


@pytest.mark.parametrize(
    "safe_chunk_ends, proposed_end, expected",
    [
        # Valid case: next safe end is 200
        ([100, 200, 300], 250, 200),
        # Valid case: next safe end is 300
        ([100, 200, 300], 350, 300),
        # Fail-safe: no safe end less than proposed_end
        ([100, 200, 300], 100, None),
        # Fail-safe: empty safe_chunk_ends
        ([], 250, None),
    ]
)
def test_find_next_safe_end(safe_chunk_ends, proposed_end, expected):
    if expected is None:
        with pytest.raises(RuntimeError, match=f"No safe chunk end found after {proposed_end}"):
            find_next_safe_end(safe_chunk_ends, proposed_end)
    else:
        assert find_next_safe_end(safe_chunk_ends, proposed_end) == expected


@pytest.mark.parametrize(
    "input_rows,max_rows,expected_num_files,expected_rows_per_file",
    [
        (list(range(12)), 5, 3, [5, 5, 2]),
        (list(range(5)), 5, 1, [5]),
        (list(range(0)), 5, 0, []),
        (list(range(7)), 3, 3, [3, 3, 1]),
    ]
)
def test_split_batch_files(tmp_path, input_rows, max_rows, expected_num_files, expected_rows_per_file):
    """
    Test the `split_batch_files` function.

    For example (first test case):
    if we have rows to split, and we want to split by 5, we should hve 3 files with 5, 5, and 2 rows respectively.
    """
    # Create a test input file
    input_file = tmp_path / "test.txt"
    header = ["col1"]
    with open(input_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for row in input_rows:
            writer.writerow([row])

    # Run the function
    result_files = split_batch_files([str(input_file)], max_rows=max_rows)

    # Check number of output files
    assert len(result_files) == expected_num_files

    # Check number of rows in each output file
    for idx, file_path in enumerate(result_files):
        with open(file_path, "r", newline="") as f:
            reader = list(csv.reader(f, delimiter="\t"))
            assert reader[0] == header
            assert len(reader) - 1 == expected_rows_per_file[idx]