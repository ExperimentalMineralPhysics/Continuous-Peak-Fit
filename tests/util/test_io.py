from typing import Any

import numpy as np
import pandas as pd
import pytest
from cpf.util.io import has_value, numpy_to_json


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
def test_numpy_to_json(
    test_params,
):
    # Unpack test params
    input, expected_value, expected_type = test_params
    output = numpy_to_json(input)
    assert output == expected_value
    assert isinstance(output, expected_type)


def test_file_list():
    pass


def test_image_list():
    pass


def test_get_file_keys():
    pass


@pytest.mark.parametrize(
    "test_params",
    (  # Object | Value | Expected result
        # =============================================================================
        # True cases (value is present)
        # =============================================================================
        # Flat dict
        (
            {
                key: value
                for key, value in (
                    (0, None),
                    (1, 1),
                    (2, 2),
                )
            },
            None,
            True,
        ),
        # Nested dict
        (
            {
                key: {
                    key: value
                    for key, value in (
                        (0, None),
                        (1, 1),
                        (2, 2),
                    )
                }
                for key in (0, 1, 2)
            },
            None,
            True,
        ),
        # Very nested dict
        (
            {0: {0: {0: {0: None}}}},
            None,
            True,
        ),
        # Flat list
        (
            [i for i in range(5)] + [None],
            None,
            True,
        ),
        # Nested list
        (
            [([i for i in range(5)] + [None]) * 5],
            None,
            True,
        ),
        # Very nested list
        (
            [
                0,
                [
                    0,
                    [
                        0,
                        [
                            0,
                            [
                                0,
                                None,
                            ],
                        ],
                    ],
                ],
            ],
            None,
            True,
        ),
        # Pandas dataframe (None)
        (
            pd.DataFrame(
                {
                    "alpha": [1.2, np.nan, 3.4, 4.5, None],
                    "beta": [5.1, 6.2, None, np.nan, 9.5],
                    "gamma": [np.nan, 2.3, 3.3, None, 5.5],
                    "delta": [7.7, 8.8, np.nan, 1.1, None],
                    "epsilon": [None, 0.5, 1.5, np.nan, 4.4],
                }
            ),
            None,
            True,
        ),
        # Pandas dataframe (Nan)
        (
            pd.DataFrame(
                {
                    "alpha": [1.2, np.nan, 3.4, 4.5, None],
                    "beta": [5.1, 6.2, None, np.nan, 9.5],
                    "gamma": [np.nan, 2.3, 3.3, None, 5.5],
                    "delta": [7.7, 8.8, np.nan, 1.1, None],
                    "epsilon": [None, 0.5, 1.5, np.nan, 4.4],
                }
            ),
            np.nan,
            True,
        ),
        # =============================================================================
        # False cases (value is not present)
        # =============================================================================
        # Flat dict
        (
            {
                key: value
                for key, value in (
                    (0, 0),
                    (1, 1),
                    (2, 2),
                )
            },
            None,
            False,
        ),
        # Nested dict
        (
            {
                key: {
                    key: value
                    for key, value in (
                        (0, 0),
                        (1, 1),
                        (2, 2),
                    )
                }
                for key in (0, 1, 2)
            },
            None,
            False,
        ),
        # Very nested dict
        (
            {0: {0: {0: {0: 0}}}},
            None,
            False,
        ),
        # Flat list
        (
            [i for i in range(5)],
            None,
            False,
        ),
        # Nested list
        (
            [[i for i in range(5)] * 5],
            None,
            False,
        ),
        # Very nested list
        (
            [
                0,
                [
                    0,
                    [
                        0,
                        [
                            0,
                            [
                                0,
                            ],
                        ],
                    ],
                ],
            ],
            None,
            False,
        ),
        # Pandas dataframe (None)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)],
                    "beta": [i for i in range(5)],
                    "gamma": [i for i in range(5)],
                    "delta": [i for i in range(5)],
                    "epsilon": [i for i in range(5)],
                }
            ),
            None,
            False,
        ),
        # Pandas dataframe (Nan)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)],
                    "beta": [i for i in range(5)],
                    "gamma": [i for i in range(5)],
                    "delta": [i for i in range(5)],
                    "epsilon": [i for i in range(5)],
                }
            ),
            np.nan,
            False,
        ),
    ),
)
def test_has_value(
    test_params: tuple[Any, Any, bool],
):
    # Unpack test params
    obj, val, result = test_params
    # Check that the result is as expected
    assert has_value(obj, val) == result


def test_replace_value():
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
