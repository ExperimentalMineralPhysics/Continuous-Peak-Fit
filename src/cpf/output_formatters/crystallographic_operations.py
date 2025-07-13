#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__all__ = ["indices4to3", "indices3to4"]

import numpy as np
from numpy.typing import ArrayLike


def indices4to3(mb_indices: ArrayLike) -> np.ndarray:
    """
    Converts set of Miller-Bravais (hkil) indices into Miller (hkl) indices

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
    m_indices[..., 0] = int(mb_indices[..., 0] - mb_indices[..., 2])
    m_indices[..., 1] = int(mb_indices[..., 1] - mb_indices[..., 2])
    m_indices[..., 2] = int(mb_indices[..., 3])

    return m_indices.squeeze()


def indices3to4(m_indices: ArrayLike) -> np.ndarray:
    """
    Cnoverts Miller (hkl) indices into Miller-Bravais (hkil) indices

    Parameters
    ----------
    m_indices : ArrayLike
        (...,3) array of Miller indices.

    Raises
    ------
    ValueError
        If there are not 4 indices.
        If the sum of the first three indices is not 0.

    Returns
    -------
    mb_indices : np.array
        (...,4) array of Miller-Bravais indices.

    """
    # check inputs are the right size
    m_indices = np.atleast_1d(m_indices)
    if m_indices.shape[-1] != 3:
        raise ValueError("This function takes 3 value Miller indices as inputs.")

    # do conversion
    mb_indices = np.zeros(m_indices.shape[:-1] + (4,), dtype=int)
    mb_indices[..., 0] = 1 / 3 * (2 * m_indices[..., 0] - m_indices[..., 1])
    mb_indices[..., 1] = 1 / 3 * (2 * m_indices[..., 1] - m_indices[..., 0])
    mb_indices[..., 2] = -(mb_indices[..., 0] + mb_indices[..., 1])
    mb_indices[..., 3] = m_indices[..., 2]

    return mb_indices.squeeze()
