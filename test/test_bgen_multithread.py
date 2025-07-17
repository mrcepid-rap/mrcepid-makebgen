import csv

import pytest

from makebgen.process_bgen.bgen_multithread import split_batch_files


@pytest.mark.parametrize(
    "input_rows, max_rows, expected_num_files, expected_rows_per_file",
    [
        (list(range(12)), 5, 3, [5, 5, 2]),
        (list(range(5)), 5, 1, [5]),
        (list(range(0)), 5, 0, []),
        (list(range(7)), 3, 3, [3, 3, 1]),
        # Edge case: max_rows is 1, each row should be in its own file
        (list(range(4)), 1, 4, [1, 1, 1, 1]),
        # Edge case: input_rows is empty, max_rows is 1
        (list(range(0)), 1, 0, []),
        # Edge case: input_rows is less than max_rows
        (list(range(2)), 5, 1, [2]),
        # Edge case: input_rows is exactly a multiple of max_rows
        (list(range(10)), 2, 5, [2, 2, 2, 2, 2]),
        # Edge case: max_rows is greater than input_rows
        (list(range(3)), 10, 1, [3]),
    ]
)
def test_split_batch_files(tmp_path, input_rows, max_rows, expected_num_files, expected_rows_per_file):
    """
    Test the `split_batch_files` function.

    For example (first test case):
    if we have rows to split, and we want to split by 5, we should have 3 files with 5, 5, and 2 rows respectively.
    """
    # Create a test input file
    input_file = tmp_path / "test.txt"
    header = ["col1"]
    with open(input_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for row in input_rows:
            writer.writerow([row])

    # Run the function with a single Path
    result_files = split_batch_files(input_file, max_rows=max_rows)

    # Check number of output files
    assert len(result_files) == expected_num_files

    # Check number of rows in each output file
    for idx, file_path in enumerate(result_files):
        with open(file_path, "r", newline="") as f:
            reader = list(csv.reader(f, delimiter="\t"))
            assert reader[0] == header
            assert len(reader) - 1 == expected_rows_per_file[idx]
