from pathlib import Path

import pytest
from bgen import BgenReader

from scripts.concat_bgens import process_chunk_group, chunk_dict_reader

# Directory containing test data files
test_data_dir = Path(__file__).parent / 'test_data'


@pytest.mark.parametrize(
    ["input_chunks", "expected_output"],
    [
        (
                [
                    {
                        'bgen': 'test_data/expected_output/chr1_test_input1.bgen',
                        'index': 'test_data/expected_output/chr1_test_input1.bgen.bgi',
                        'sample': 'test_data/expected_output/chr1_test_input1.sample',
                        'vep': 'test_data/expected_output/chr1_test_input1.vep.tsv.gz',
                        'vep_idx': 'test_data/expected_output/chr1_test_input1.vep.tsv.gz.tbi'
                    },
                    {
                        'bgen': 'test_data/expected_output/chr1_test_input2.bgen',
                        'index': 'test_data/expected_output/chr1_test_input2.bgen.bgi',
                        'sample': 'test_data/expected_output/chr1_test_input2.sample',
                        'vep': 'test_data/expected_output/chr1_test_input2.vep.tsv.gz',
                        'vep_idx': 'test_data/expected_output/chr1_test_input2.vep.tsv.gz.tbi'
                    }
                ],
                {
                    'bgen': 'chr1_test_input1_mergedchunk_1.bgen',
                    'bgen_index': 'chr1_test_input1_mergedchunk_1.bgen.bgi',
                    'sample': 'chr1_test_input1_mergedchunk_1.sample',
                    'vep': 'chr1_test_input1_mergedchunk_1.vep.tsv.gz',
                    'vep_index': 'chr1_test_input1_mergedchunk_1.vep.tsv.gz.tbi'
                }
        )
    ]
)
def test_process_chunk_group(input_chunks, expected_output):
    # Simulate chunks of 3 (or fewer) items
    for chunk_index, chunk in enumerate(chunk_dict_reader(input_chunks, chunk_size=3), start=1):
        output = process_chunk_group(chunk=chunk, batch_index=chunk_index)

        # Assert all expected output files exist
        for key, expected_file in expected_output.items():
            actual_file = Path(output[key])
            assert actual_file.exists(), f"{key} file not found: {actual_file}"

        # Count total variants in original files
        raw_total_variants = 0
        for row in chunk:
            raw_bgen_path = Path(row["bgen"])
            with BgenReader(raw_bgen_path, delay_parsing=True) as raw_bgen:
                raw_total_variants += len(raw_bgen.positions())

        # Count variants in merged BGEN
        with BgenReader(Path(output['bgen']), delay_parsing=True) as merged_bgen:
            merged_positions = merged_bgen.positions()
            assert len(merged_positions) == raw_total_variants, (
                f"Merged BGEN variant count ({len(merged_positions)}) "
                f"does not match sum of inputs ({raw_total_variants})"
            )

        # Cleanup created files
        for file in output.values():
            path = Path(file)
            if path.exists():
                path.unlink()
