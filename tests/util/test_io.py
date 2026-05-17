from typing import Any, TypeVar

import numpy as np
import pandas as pd
import pytest
from cpf.util.io import (
    has_value,
    numpy_to_json,
    peak_hkl,
    peak_phase,
    peak_string,
    replace_value,
)

T = TypeVar("T", dict, list, pd.DataFrame)


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


def test_get_file_indices():
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


@pytest.mark.parametrize(
    "test_params",
    (  # Initial object | Old value | New value | Expected output
        # Simple dict
        (
            {
                0: None,
                1: 1,
                2: None,
                3: 3,
            },
            None,
            0,
            {
                0: 0,
                1: 1,
                2: 0,
                3: 3,
            },
        ),
        # Nested dict
        (
            {
                0: {
                    0: None,
                    1: 1,
                    2: None,
                    3: 3,
                },
                1: {
                    0: None,
                    1: 1,
                    2: None,
                    3: 3,
                },
            },
            None,
            0,
            {
                0: {
                    0: 0,
                    1: 1,
                    2: 0,
                    3: 3,
                },
                1: {
                    0: 0,
                    1: 1,
                    2: 0,
                    3: 3,
                },
            },
        ),
        # Simple list
        (
            [None, 1, 2, None, 4],
            None,
            0,
            [0, 1, 2, 0, 4],
        ),
        # Nested list
        (
            [[None, 1, 2, None, 4] * 5],
            None,
            0,
            [[0, 1, 2, 0, 4] * 5],
        ),
        # Pandas DataFrame (None)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [None],
                    "beta": [i for i in range(5)] + [None],
                }
            ),
            None,
            0,
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [0],
                    "beta": [i for i in range(5)] + [0],
                }
            ),
        ),
        # Pandas DataFrame (NaN)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [np.nan],
                    "beta": [i for i in range(5)] + [np.nan],
                }
            ),
            np.nan,
            0,
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [0],
                    "beta": [i for i in range(5)] + [0],
                }
            ),
        ),
        # Pandas DataFrame (check that NaN and None are interchangeable)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [None],
                    "beta": [i for i in range(5)] + [None],
                }
            ),
            np.nan,
            0,
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [0],
                    "beta": [i for i in range(5)] + [0],
                }
            ),
        ),
        # Pandas DataFrame (check that NaN and None are interchangeable)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [np.nan],
                    "beta": [i for i in range(5)] + [np.nan],
                }
            ),
            None,
            np.nan,
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [None],
                    "beta": [i for i in range(5)] + [None],
                }
            ),
        ),
        # Pandas DataFrame (check that replacing with None works)
        (
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [np.nan],
                    "beta": [i for i in range(5)] + [np.nan],
                }
            ),
            None,
            None,
            pd.DataFrame(
                {
                    "alpha": [i for i in range(5)] + [None],
                    "beta": [i for i in range(5)] + [None],
                }
            ),
        ),
    ),
)
def test_replace_value(test_params: tuple[T, Any, Any, T]):
    # Unpack test params
    obj, old, new, out = test_params
    if isinstance(obj, pd.DataFrame) and isinstance(out, pd.DataFrame):
        assert replace_value(obj, old, new).astype(float).equals(out.astype(float))
    else:
        assert replace_value(obj, old, new) == out


def test_has_huge_errors():
    pass


@pytest.mark.parametrize(
    "test_params",
    (  # Fit order | Use file name? | Peak | Expected output
        (
            {
                "peak": [
                    {
                        "phase": "Other",
                        "hkl": "000",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": 110,
                    },
                ],
            },
            False,
            "all",
            "Other (000) & Fe-BCC (110)",
        ),
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                        "hkl": "110",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "220",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "330",
                    },
                ]
            },
            True,
            "all",
            "Fe-BCC-110_Fe-BCC-220_Fe-BCC-330",
        ),
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                        "hkl": "110",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "220",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "330",
                    },
                ]
            },
            False,
            0,
            "Fe-BCC (110)",
        ),
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                        "hkl": "110",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "220",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "330",
                    },
                ]
            },
            True,
            [0, 1],
            "Fe-BCC-110_Fe-BCC-220",
        ),
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                        "hkl": "110",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "220",
                    },
                    {
                        "phase": "Fe-BCC",
                        "hkl": "330",
                    },
                ]
            },
            False,
            "1,2,",
            "Fe-BCC (220) & Fe-BCC (330)",
        ),
        # If 'phase' key is absent
        (
            {
                "peak": [
                    {
                        "hkl": "110",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "330",
                    },
                ]
            },
            False,
            "1,2",
            "Peak (220) & Peak (330)",
        ),
        # If 'hkl' key is absent
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                    },
                    {
                        "phase": "Fe-BCC",
                    },
                    {
                        "phase": "Fe-BCC",
                    },
                ]
            },
            False,
            "1,2",
            "Fe-BCC (2) & Fe-BCC (3)",
        ),
    ),
)
def test_peak_string(
    test_params: tuple[dict[str, Any], bool, str | int | list[int], str],
):
    # Unpack the test params
    fit_order, as_filename, peak, output = test_params
    assert peak_string(fit_order, peak=peak, fname=as_filename) == output


@pytest.mark.parametrize(
    "test_params",
    (  # Fit order | Peak | As a string? | Expected output
        (
            {
                "peak": [
                    {
                        "hkl": "110",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "330",
                    },
                ]
            },
            "all",
            True,
            ["110", "220", "330"],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": "110",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "330",
                    },
                ]
            },
            "all",
            False,
            [[1, 1, 0], [2, 2, 0], [3, 3, 0]],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": "-1-10",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "-3-30",
                    },
                ]
            },
            "all",
            True,
            ["-1-10", "220", "-3-30"],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": [1, 1, 0],
                    },
                    {
                        "hkl": [2, 2, 0],
                    },
                    {
                        "hkl": [3, 3, 0],
                    },
                ]
            },
            "all",
            False,
            [[1, 1, 0], [2, 2, 0], [3, 3, 0]],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": [-1, -1, 0],
                    },
                    {
                        "hkl": [2, 2, 0],
                    },
                    {
                        "hkl": [-3, -3, 0],
                    },
                ]
            },
            "all",
            True,
            ["-1-10", "220", "-3-30"],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": "-1-10",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "-3-30",
                    },
                ]
            },
            "all",
            False,
            [[-1, -1, 0], [2, 2, 0], [-3, -3, 0]],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": "110",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "330",
                    },
                ]
            },
            "0,1",
            True,
            ["110", "220"],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": "110",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "330",
                    },
                ]
            },
            [1, 2],
            False,
            [[2, 2, 0], [3, 3, 0]],
        ),
        (
            {
                "peak": [
                    {
                        "hkl": "-1-10",
                    },
                    {
                        "hkl": "220",
                    },
                    {
                        "hkl": "-3-30",
                    },
                ]
            },
            0,
            False,
            [[-1, -1, 0]],
        ),
        # If no 'hkl' key is provided
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                    },
                    {
                        "phase": "Fe-BCC",
                    },
                    {
                        "phase": "Fe-BCC",
                    },
                ]
            },
            [1, 2],
            False,
            [[0, 0, 0], [0, 0, 0]],
        ),
        (
            {
                "peak": [
                    {
                        "phase": "Fe-BCC",
                    },
                    {
                        "phase": "Fe-BCC",
                    },
                    {
                        "phase": "Fe-BCC",
                    },
                ]
            },
            [1, 2],
            True,
            ["000", "000"],
        ),
    ),
)
def test_peak_hkl(
    test_params: tuple[
        dict[str, Any], str | int | list[int], bool, list[str | list[int]]
    ],
):
    # Unpack test params
    orders, peak, as_string, output = test_params
    assert output == peak_hkl(orders, peak=peak, string=as_string)


@pytest.mark.parametrize(
    "test_params",
    (  # Orders | Peak to analyse | Expected output
        (
            {
                "peak": [
                    {"phase": "Fe-BCC"},
                    {"phase": "CaO-FCC"},
                    {"phase": "Cr-BCC"},
                ],
            },
            "all",
            ["Fe-BCC", "CaO-FCC", "Cr-BCC"],
        ),
        (
            {
                "peak": [
                    {"phase": "Fe-BCC"},
                    {"phase": "CaO-FCC"},
                    {"phase": "Cr-BCC"},
                ],
            },
            0,
            ["Fe-BCC"],
        ),
        (
            {
                "peak": [
                    {"phase": "Fe-BCC"},
                    {"phase": "CaO-FCC"},
                    {"phase": "Cr-BCC"},
                ],
            },
            [0, 1],
            ["Fe-BCC", "CaO-FCC"],
        ),
        (
            {
                "peak": [
                    {"phase": "Fe-BCC"},
                    {"phase": "CaO-FCC"},
                    {"phase": "Cr-BCC"},
                ],
            },
            "0,1,2",
            ["Fe-BCC", "CaO-FCC", "Cr-BCC"],
        ),
        (
            {
                "peak": [
                    {"hkl": "110"},
                    {"hkl": "220"},
                    {"hkl": "330"},
                ],
            },
            "all",
            ["Unknown", "Unknown", "Unknown"],
        ),
    ),
)
def test_peak_phase(
    test_params: tuple[dict[str, Any], str | int | list[int], list[str]],
):
    # Unpack test params
    orders, peak, output = test_params
    assert peak_phase(orders, peak) == output


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
