from typing import Callable

import numpy as np
import pytest
from cpf.output_formatters.crystallographic_operations import (
    plane_indices_3_to_4,
    plane_indices_4_to_3,
    vector_indices_3_to_4,
    vector_indices_4_to_3,
)
from numpy.testing import assert_equal

plane_indices_test_matrix = (
    # Miller indices | Miller-Bravais indices | Input Type
    ((0, 0, 1), (0, 0, 0, 1), tuple),
    ((0, 1, 0), (0, 1, -1, 0), np.ndarray),
    ((0, 1, 1), (0, 1, -1, 1), list),
    ((1, 0, 0), (1, 0, -1, 0), tuple),
    ((1, 0, 1), (1, 0, -1, 1), np.ndarray),
    ((1, 1, 0), (1, 1, -2, 0), list),
    ((1, 1, 1), (1, 1, -2, 1), tuple),
    ((0, 0, 2), (0, 0, 0, 2), np.ndarray),
    ((0, 1, 2), (0, 1, -1, 2), list),
    ((0, 2, 0), (0, 2, -2, 0), tuple),
    ((0, 2, 1), (0, 2, -2, 1), np.ndarray),
    ((0, 2, 2), (0, 2, -2, 2), list),
    ((1, 0, 2), (1, 0, -1, 2), tuple),
    ((1, 1, 2), (1, 1, -2, 2), np.ndarray),
    ((1, 2, 0), (1, 2, -3, 0), list),
    ((1, 2, 1), (1, 2, -3, 1), tuple),
    ((1, 2, 2), (1, 2, -3, 2), np.ndarray),
    ((2, 0, 0), (2, 0, -2, 0), list),
    ((2, 0, 1), (2, 0, -2, 1), tuple),
    ((2, 0, 2), (2, 0, -2, 2), np.ndarray),
    ((2, 1, 0), (2, 1, -3, 0), list),
    ((2, 1, 1), (2, 1, -3, 1), tuple),
    ((2, 1, 2), (2, 1, -3, 2), np.ndarray),
    ((2, 2, 0), (2, 2, -4, 0), list),
    ((2, 2, 1), (2, 2, -4, 1), tuple),
    ((2, 2, 2), (2, 2, -4, 2), np.ndarray),
    # Placeholder for copying for if we want to add more
    # ((), (), list),
    # ((), (), tuple),
    # ((), (), np.ndarray),
)


@pytest.mark.parametrize("test_params", plane_indices_test_matrix)
def test_plane_indices_3_to_4(
    test_params: tuple[tuple[int, ...], tuple[int, ...], type],
):
    # Unpack test parameters
    m_indices, mb_indices, input_type = test_params

    # Convert input indices (M indices) to correct type
    if input_type == np.ndarray:
        input_value = np.array(m_indices)
    else:
        input_value = input_type(m_indices)

    # Check that object conversion has worked
    assert isinstance(input_value, input_type)

    # Run the function
    assert_equal(
        plane_indices_3_to_4(input_value),
        np.array(mb_indices),
    )


@pytest.mark.parametrize("test_params", plane_indices_test_matrix)
def test_plane_indices_4_to_3(
    test_params: tuple[tuple[int, ...], tuple[int, ...], type],
):
    # Unpack test parameters
    m_indices, mb_indices, input_type = test_params

    # Convert input indices (M indices) to correct type
    if input_type == np.ndarray:
        input_value = np.array(mb_indices)
    else:
        input_value = input_type(mb_indices)

    # Check that object conversion has worked
    assert isinstance(input_value, input_type)

    # Run the function
    assert_equal(
        plane_indices_4_to_3(input_value),
        np.array(m_indices),
    )


vector_indices_test_matrix = (
    # Miller indices | Miller-Bravais indices | Input Type
    ((0, 0, 1), (0, 0, 0, 1), list),
    ((0, 0, 2), (0, 0, 0, 1), tuple),
    ((0, 0, 3), (0, 0, 0, 1), np.ndarray),
    ((0, 0, 4), (0, 0, 0, 1), list),
    ((0, 1, 0), (-1, 2, -1, 0), tuple),
    ((0, 1, 1), (-1, 2, -1, 3), np.ndarray),
    ((0, 1, 2), (-1, 2, -1, 6), list),
    ((0, 1, 3), (-1, 2, -1, 9), tuple),
    ((0, 1, 4), (-1, 2, -1, 12), np.ndarray),
    ((0, 2, 0), (-1, 2, -1, 0), list),
    ((0, 2, 1), (-2, 4, -2, 3), tuple),
    ((0, 2, 2), (-1, 2, -1, 3), np.ndarray),
    ((0, 2, 3), (-2, 4, -2, 9), list),
    ((0, 2, 4), (-1, 2, -1, 6), tuple),
    ((0, 3, 0), (-1, 2, -1, 0), np.ndarray),
    ((0, 3, 1), (-1, 2, -1, 1), list),
    ((0, 3, 2), (-1, 2, -1, 2), tuple),
    ((0, 3, 3), (-1, 2, -1, 3), np.ndarray),
    ((0, 3, 4), (-1, 2, -1, 4), list),
    ((0, 4, 0), (-1, 2, -1, 0), tuple),
    ((0, 4, 1), (-4, 8, -4, 3), np.ndarray),
    ((0, 4, 2), (-2, 4, -2, 3), list),
    ((0, 4, 3), (-4, 8, -4, 9), tuple),
    ((0, 4, 4), (-1, 2, -1, 3), np.ndarray),
    ((1, 0, 0), (2, -1, -1, 0), list),
    ((1, 0, 1), (2, -1, -1, 3), tuple),
    ((1, 0, 2), (2, -1, -1, 6), np.ndarray),
    ((1, 0, 3), (2, -1, -1, 9), list),
    ((1, 0, 4), (2, -1, -1, 12), tuple),
    ((1, 1, 0), (1, 1, -2, 0), np.ndarray),
    ((1, 1, 1), (1, 1, -2, 3), list),
    ((1, 1, 2), (1, 1, -2, 6), tuple),
    ((1, 1, 3), (1, 1, -2, 9), np.ndarray),
    ((1, 1, 4), (1, 1, -2, 12), list),
    ((1, 2, 0), (0, 1, -1, 0), tuple),
    ((1, 2, 1), (0, 1, -1, 1), np.ndarray),
    ((1, 2, 2), (0, 1, -1, 2), list),
    ((1, 2, 3), (0, 1, -1, 3), tuple),
    ((1, 2, 4), (0, 1, -1, 4), np.ndarray),
    ((1, 3, 0), (-1, 5, -4, 0), list),
    ((1, 3, 1), (-1, 5, -4, 3), tuple),
    ((1, 3, 2), (-1, 5, -4, 6), np.ndarray),
    ((1, 3, 3), (-1, 5, -4, 9), list),
    ((1, 3, 4), (-1, 5, -4, 12), tuple),
    ((1, 4, 0), (-2, 7, -5, 0), np.ndarray),
    ((1, 4, 1), (-2, 7, -5, 3), list),
    ((1, 4, 2), (-2, 7, -5, 6), tuple),
    ((1, 4, 3), (-2, 7, -5, 9), np.ndarray),
    ((1, 4, 4), (-2, 7, -5, 12), list),
    ((2, 0, 0), (2, -1, -1, 0), tuple),
    ((2, 0, 1), (4, -2, -2, 3), np.ndarray),
    ((2, 0, 2), (2, -1, -1, 3), list),
    ((2, 0, 3), (4, -2, -2, 9), tuple),
    ((2, 0, 4), (2, -1, -1, 6), np.ndarray),
    ((2, 1, 0), (1, 0, -1, 0), list),
    ((2, 1, 1), (1, 0, -1, 1), tuple),
    ((2, 1, 2), (1, 0, -1, 2), np.ndarray),
    ((2, 1, 3), (1, 0, -1, 3), list),
    ((2, 1, 4), (1, 0, -1, 4), tuple),
    ((2, 2, 0), (1, 1, -2, 0), np.ndarray),
    ((2, 2, 1), (2, 2, -4, 3), list),
    ((2, 2, 2), (1, 1, -2, 3), tuple),
    ((2, 2, 3), (2, 2, -4, 9), np.ndarray),
    ((2, 2, 4), (1, 1, -2, 6), list),
    ((2, 3, 0), (1, 4, -5, 0), tuple),
    ((2, 3, 1), (1, 4, -5, 3), np.ndarray),
    ((2, 3, 2), (1, 4, -5, 6), list),
    ((2, 3, 3), (1, 4, -5, 9), tuple),
    ((2, 3, 4), (1, 4, -5, 12), np.ndarray),
    ((2, 4, 0), (0, 1, -1, 0), list),
    ((2, 4, 1), (0, 2, -2, 1), tuple),
    ((2, 4, 2), (0, 1, -1, 1), np.ndarray),
    ((2, 4, 3), (0, 2, -2, 3), list),
    ((2, 4, 4), (0, 1, -1, 2), tuple),
    ((3, 0, 0), (2, -1, -1, 0), np.ndarray),
    ((3, 0, 1), (2, -1, -1, 1), list),
    ((3, 0, 2), (2, -1, -1, 2), tuple),
    ((3, 0, 3), (2, -1, -1, 3), np.ndarray),
    ((3, 0, 4), (2, -1, -1, 4), list),
    ((3, 1, 0), (5, -1, -4, 0), tuple),
    ((3, 1, 1), (5, -1, -4, 3), np.ndarray),
    ((3, 1, 2), (5, -1, -4, 6), list),
    ((3, 1, 3), (5, -1, -4, 9), tuple),
    ((3, 1, 4), (5, -1, -4, 12), np.ndarray),
    ((3, 2, 0), (4, 1, -5, 0), list),
    ((3, 2, 1), (4, 1, -5, 3), tuple),
    ((3, 2, 2), (4, 1, -5, 6), np.ndarray),
    ((3, 2, 3), (4, 1, -5, 9), list),
    ((3, 2, 4), (4, 1, -5, 12), tuple),
    ((3, 3, 0), (1, 1, -2, 0), np.ndarray),
    ((3, 3, 1), (1, 1, -2, 1), list),
    ((3, 3, 2), (1, 1, -2, 2), tuple),
    ((3, 3, 3), (1, 1, -2, 3), np.ndarray),
    ((3, 3, 4), (1, 1, -2, 4), list),
    ((3, 4, 0), (2, 5, -7, 0), tuple),
    ((3, 4, 1), (2, 5, -7, 3), np.ndarray),
    ((3, 4, 2), (2, 5, -7, 6), list),
    ((3, 4, 3), (2, 5, -7, 9), tuple),
    ((3, 4, 4), (2, 5, -7, 12), np.ndarray),
    ((4, 0, 0), (2, -1, -1, 0), list),
    ((4, 0, 1), (8, -4, -4, 3), tuple),
    ((4, 0, 2), (4, -2, -2, 3), np.ndarray),
    ((4, 0, 3), (8, -4, -4, 9), list),
    ((4, 0, 4), (2, -1, -1, 3), tuple),
    ((4, 1, 0), (7, -2, -5, 0), np.ndarray),
    ((4, 1, 1), (7, -2, -5, 3), list),
    ((4, 1, 2), (7, -2, -5, 6), tuple),
    ((4, 1, 3), (7, -2, -5, 9), np.ndarray),
    ((4, 1, 4), (7, -2, -5, 12), list),
    ((4, 2, 0), (1, 0, -1, 0), tuple),
    ((4, 2, 1), (2, 0, -2, 1), np.ndarray),
    ((4, 2, 2), (1, 0, -1, 1), list),
    ((4, 2, 3), (2, 0, -2, 3), tuple),
    ((4, 2, 4), (1, 0, -1, 2), np.ndarray),
    ((4, 3, 0), (5, 2, -7, 0), list),
    ((4, 3, 1), (5, 2, -7, 3), tuple),
    ((4, 3, 2), (5, 2, -7, 6), np.ndarray),
    ((4, 3, 3), (5, 2, -7, 9), list),
    ((4, 3, 4), (5, 2, -7, 12), tuple),
    ((4, 4, 0), (1, 1, -2, 0), np.ndarray),
    ((4, 4, 1), (4, 4, -8, 3), list),
    ((4, 4, 2), (2, 2, -4, 3), tuple),
    ((4, 4, 3), (4, 4, -8, 9), np.ndarray),
    ((4, 4, 4), (1, 1, -2, 3), list),
    # Placeholder for copying for if we want to add more
    # ((), (), list),
    # ((), (), tuple),
    # ((), (), np.ndarray),
)


@pytest.mark.parametrize("test_params", vector_indices_test_matrix)
def test_vector_indices_3_to_4(
    test_params: tuple[tuple[int, ...], tuple[int, ...], type],
):
    # Unpack test parameters
    m_indices, mb_indices, input_type = test_params

    # Convert input indices (M indices) to correct type
    if input_type == np.ndarray:
        input_value = np.array(m_indices)
    else:
        input_value = input_type(m_indices)

    # Check that object conversion has worked
    assert isinstance(input_value, input_type)

    # Run the function
    assert_equal(
        vector_indices_3_to_4(input_value),
        np.array(mb_indices) / np.gcd.reduce(mb_indices),
    )


@pytest.mark.parametrize("test_params", vector_indices_test_matrix)
def test_vector_indices_4_to_3(
    test_params: tuple[tuple[int, ...], tuple[int, ...], type],
):
    # Unpack test parameters
    m_indices, mb_indices, input_type = test_params

    # Convert input indices (MB indices) to correct type
    if input_type == np.ndarray:
        input_value = np.array(mb_indices)
    else:
        input_value = input_type(mb_indices)

    # Check that object conversion has worked
    assert isinstance(input_value, input_type)

    # Run the function
    assert_equal(
        vector_indices_4_to_3(input_value),
        (np.array(m_indices) / np.gcd.reduce(m_indices)),
    )
