import hashlib
import gzip

import pandas as pd

from BEstimate.datafiles import DataFiles


def test_DataFiles__get_home_sapiens_interfaces_path__exists() -> None:
    # Given
    hash_md5_compressed = hashlib.md5()
    hash_md5_decompressed = hashlib.md5()
    expected_compressed_md5_sum = DataFiles.MD5_H_SAPIENS_INTERFACES_GZIP
    expected_decompressed_md5_sum = DataFiles.MD5_H_SAPIENS_INTERFACES_TXT

    # When
    gzip_file = DataFiles.get_home_sapiens_interfaces_path()

    # Then
    assert gzip_file.exists()
    assert gzip_file.is_file()
    hash_md5_compressed.update(gzip_file.read_bytes())
    assert hash_md5_compressed.hexdigest() == expected_compressed_md5_sum

    # Finally - round-robin - decompress to ensure it is a valid gzip file
    with gzip.open(gzip_file, mode="rb") as f:
        file_as_bytes = f.read()
        hash_md5_decompressed.update(file_as_bytes)
    assert hash_md5_decompressed.hexdigest() == expected_decompressed_md5_sum


def test_DataFiles__get_home_sapiens_interfaces_as_dataframe():
    # Given
    expected_columns = DataFiles.H_SAPIENS_INTERFACES_COLUMNS
    expected_row_count = 121575
    expected_shape = (expected_row_count, len(expected_columns))

    # When
    df = DataFiles.get_home_sapiens_interfaces_as_dataframe()

    # Then
    assert not df.empty
    assert list(df.columns) == list(expected_columns)
    assert df.shape == expected_shape
