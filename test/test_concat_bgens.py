import csv
from pathlib import Path
from typing import Iterator, List, Dict

import pytest
from general_utilities.association_resources import check_gzipped
from general_utilities.import_utils.file_handlers.input_file_handler import InputFileHandler

from scripts.concat_bgens import process_chunk_group, chunk_dict_reader

# Directory containing test data files
test_data_dir = Path(__file__).parent / 'test_data'


@pytest.mark.parametrize(
    ["input_coords"],
    [
        ("concat_coords.txt",)
    ]
)
def test_process_chunk_group(input_coords):
    data = test_data_dir / input_coords
    coords = InputFileHandler(data).get_file_handle()

    with check_gzipped(coords) as coord_file:
        coord_file_reader = csv.DictReader(coord_file, delimiter="\t")
        for chunk_index, chunk in enumerate(chunk_dict_reader(coord_file_reader, chunk_size=3), start=1):
            # process the chunk group
            output = process_chunk_group(chunk=chunk, batch_index=chunk_index)

            for file in output.values():
                file = Path(file)
                if file.exists():
                    assert file.is_file(), f"Path {file} is not a file"
                    assert file.stat().st_size > 0, f"File {file} is empty"
                    file.unlink()
                else:
                    print(f"Warning: File {file} does not exist and will not be deleted.")

