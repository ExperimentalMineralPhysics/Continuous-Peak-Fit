import numpy as np
import pytest
from cpf.util.io import json_numpy_serializer


@pytest.mark.parametrize(
    "test_params",
    # NumPy object | Expected value | Expected type
    (
        # Individual objects
        (np.bool(True), True, bool),
        (np.bool(False), False, bool),
        (np.bool(-2), True, bool),
        (np.bool(-1), True, bool),
        (np.bool(0), False, bool),
        (np.bool(1), True, bool),
        (np.bool(2), True, bool),
        (np.float16(12), 12.0, float),
        (np.float32(12), 12.0, float),
        (np.float64(12), 12.0, float),
        (np.int8(-1), -1, int),
        (np.int8(0), 0, int),
        (np.int8(1), 1, int),
        (np.int16(-1), -1, int),
        (np.int16(0), 0, int),
        (np.int16(1), 1, int),
        (np.int32(-1), -1, int),
        (np.int32(0), 0, int),
        (np.int32(1), 1, int),
        (np.int64(-1), -1, int),
        (np.int64(0), 0, int),
        (np.int64(1), 1, int),  #
        (np.uint8(0), 0, int),
        (np.uint8(1), 1, int),
        (np.uint16(0), 0, int),
        (np.uint16(1), 1, int),
        (np.uint32(0), 0, int),
        (np.uint32(1), 1, int),
        (np.uint64(0), 0, int),
        (np.uint64(1), 1, int),
        (np.str_(123), "123", str),
        (np.timedelta64(123), 123, int),
        (np.void(123), bytes(123), bytes),
        # Arrays
        (
            np.array([list(range(-2, 3))] * 5, dtype=np.int8),
            [list(range(-2, 3))] * 5,
            list,
        ),
        (
            np.array([list(range(-2, 3))] * 5, dtype=np.int16),
            [list(range(-2, 3))] * 5,
            list,
        ),
        (
            np.array([list(range(-2, 3))] * 5, dtype=np.int32),
            [list(range(-2, 3))] * 5,
            list,
        ),
        (
            np.array([list(range(-2, 3))] * 5, dtype=np.int64),
            [list(range(-2, 3))] * 5,
            list,
        ),
        (
            np.array([list(range(5))] * 5, dtype=np.uint8),
            [list(range(5))] * 5,
            list,
        ),
        (
            np.array([list(range(5))] * 5, dtype=np.uint16),
            [list(range(5))] * 5,
            list,
        ),
        (
            np.array([list(range(5))] * 5, dtype=np.uint32),
            [list(range(5))] * 5,
            list,
        ),
        (
            np.array([list(range(5))] * 5, dtype=np.uint64),
            [list(range(5))] * 5,
            list,
        ),
        (
            np.array([[i for i in range(5)]] * 5, dtype=np.float16),
            [[float(i) for i in range(5)]] * 5,
            list,
        ),
        (
            np.array([[i for i in range(5)]] * 5, dtype=np.float32),
            [[float(i) for i in range(5)]] * 5,
            list,
        ),
        (
            np.array([[i for i in range(5)]] * 5, dtype=np.float64),
            [[float(i) for i in range(5)]] * 5,
            list,
        ),
    ),
)
def test_json_numpy_serializer(
    test_params,
):
    # Unpack test params
    input, expected_value, expected_type = test_params
    output = json_numpy_serializer(input)
    assert output == expected_value
    assert isinstance(output, expected_type)


def test_file_list():
    pass


def test_image_list():
    pass


def test_get_file_keys():
    pass


def test_any_terms_null():
    pass


def test_replace_null_terms():
    pass


def test_any_errors_huge():
    pass


def test_peak_string():
    pass


def test_peak_hkl():
    pass


def test_peak_phase():
    pass


def test_title_file_names():
    pass


def test_make_outfile_name():
    pass


def test_lmfit_fix_int_data_type():
    pass


def test_number_to_string():
    pass


def test_licit_filename():
    pass


def test_figure_suptitle_space():
    pass
