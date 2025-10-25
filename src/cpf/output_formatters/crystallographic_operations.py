#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__all__ = ["plane_indices_4_to_3", "vector_indices_3_to_4"]

import numpy as np
from numpy.typing import ArrayLike


def plane_indices_3_to_4(m_indices: ArrayLike) -> np.ndarray:
    """
    Converts set of Miller plane indices (hkl) into Miller-Bravais indices (hkil).

    Parameters
    ----------
    m_indices : ArrayLike
        (...,3) array of Miller indices.

    Raises
    ------
    ValueError
        If there are not 3 indices.

    Returns
    -------
    m_indices : np.array
        (...,4) array of Miller-Bravais indices.

    """

    # check inputs are the right size
    m_indices = np.atleast_1d(m_indices)
    if m_indices.shape[-1] != 3:
        raise ValueError("This function takes 3 value Miller indices as inputs.")

    # do conversion
    mb_indices = np.zeros(m_indices.shape[:-1] + (4,), dtype=int)
    mb_indices[..., 0] = m_indices[..., 0]
    mb_indices[..., 1] = m_indices[..., 1]
    mb_indices[..., 2] = -(m_indices[..., 0] + m_indices[..., 1])
    mb_indices[..., 3] = m_indices[..., 2]

    return mb_indices.squeeze()


def plane_indices_4_to_3(mb_indices: ArrayLike) -> np.ndarray:
    """
    Converts set of Miller-Bravais plane indices (hkil) into Miller (hkl) indices

    Parameters
    ----------
    mb_indices : ArrayLike
        (...,4) array of Miller-Bravais indices.

    Raises
    ------
    ValueError
        If there are not 4 indices.
        If the sum of the first three indices is not 0.

    Returns
    -------
    m_indices : np.array
        (...,3) array of Miller indices.

    """

    # check inputs are the right size
    mb_indices = np.atleast_1d(mb_indices)
    if mb_indices.shape[-1] != 4:
        raise ValueError(
            "This function takes 4 value Miller-Bravais indices as inputs."
        )
    if not np.allclose(mb_indices[..., :3].sum(axis=-1), 0.0):
        raise ValueError("Invalid inputs: the first three values must sum to zero")

    # do conversion
    m_indices = np.zeros(mb_indices.shape[:-1] + (3,), dtype=int)
    m_indices[..., 0] = mb_indices[..., 0]
    m_indices[..., 1] = mb_indices[..., 1]
    m_indices[..., 2] = mb_indices[..., 3]

    return m_indices.squeeze()


def vector_indices_3_to_4(m_indices: ArrayLike) -> np.ndarray:
    """
    Converts Miller vector indices (uvw) into Miller-Bravais indices (UVTW)

    Parameters
    ----------
    m_indices : ArrayLike
        (...,3) array of Miller vector indices.

    Raises
    ------
    ValueError
        If there are not 3 indices.

    Returns
    -------
    mb_indices : np.array
        (...,4) array of Miller-Bravais indices.
    """
    # check inputs are the right size
    m_indices = np.atleast_1d(m_indices)
    if m_indices.shape[-1] != 3:
        raise ValueError("This function takes 3 value Miller indices as inputs.")

    # Convert from uvw to UVTW
    mb_indices = np.zeros(m_indices.shape[:-1] + (4,), dtype=int)
    mb_indices[..., 0] = 2 * m_indices[..., 0] - m_indices[..., 1]
    mb_indices[..., 1] = 2 * m_indices[..., 1] - m_indices[..., 0]
    mb_indices[..., 2] = -(mb_indices[..., 0] + mb_indices[..., 1])
    mb_indices[..., 3] = m_indices[..., 2] * 3

    # Divide by highest common factor
    mb_indices = mb_indices / np.gcd.reduce(mb_indices)

    return mb_indices.squeeze()


def vector_indices_4_to_3(mb_indices: ArrayLike):
    """
    Converts Miller-Bravais vector indices (UVTW) into Miller indices (uvw)

    Parameters
    ----------
    mb_indices : ArrayLike
        (...,4) array of Miller-Bravais indices.

    Raises
    ------
    ValueError
        If there are not 4 indices.
        If the sum of the first three indices is not 0.

    Returns
    -------
    m_indices : np.array
        (...,3) array of Miller indices.
    """

    # Check inputs are the right size and are valid Miller-Bravais indices
    mb_indices = np.atleast_1d(mb_indices)
    if mb_indices.shape[-1] != 4:
        raise ValueError(
            "This function takes 4 value Miller-Bravais indices as inputs."
        )
    if not np.allclose(mb_indices[..., :3].sum(axis=-1), 0.0):
        raise ValueError(
            "This is not a valid Miller-Bravais index. "
            "The sum of the first three indices must be equal to 0"
        )

    # Convert from UVTW to uvw
    m_indices = np.zeros(mb_indices.shape[:-1] + (3,), dtype=int)
    m_indices[0] = 2 * mb_indices[..., 0] + mb_indices[..., 1]
    m_indices[1] = 2 * mb_indices[..., 1] + mb_indices[..., 0]
    m_indices[2] = mb_indices[..., 3]
    m_indices = m_indices / np.gcd.reduce(m_indices)

    return m_indices.squeeze()
