#!/usr/bin/env python3
# -*- coding: utf-8 -*-



__all__ = ["indicies4to3", "indicies3to4"]

import numpy as np



def indicies4to3(MBindicies: np.typing.ArrayLike) -> np.ndarray:
    """
    Cnoverts set of Miller-Bravais (hkil) indicies into Miller (hkl) indicies

    Parameters
    ----------
    Mindicies : np.typing.ArrayLike
        (...,4) array of Miller indicies.

    Raises
    ------
    ValueError
        If there are not 4 indicies.
        If the sum of the first three indicies is not 0.

    Returns
    -------
    MBindicies : np.array
        (...,3) array of Miller-Bravais indicies.

    """
    """
    
    Returns
    -------
    None.

    """

    #check inputs are the right size
    MBindicies = np.atleast_1d(MBindicies)
    if MBindicies.shape[-1] != 4:
        raise ValueError("This function takes 4 value Miller-Bravais indicies as inputs.")
    if not np.allclose(MBindicies[...,:3].sum(axis=-1), 0.0):
        raise ValueError('Invalid inputs: the first three values must sum to zero')

    #do conversion
    Mindicies = np.zeros(MBindicies.shape[:-1]+ (3,), dtype=int)
    Mindicies[..., 0] = int(MBindicies[..., 0] - MBindicies[..., 2])
    Mindicies[..., 1] = int(MBindicies[..., 1] - MBindicies[..., 2])
    Mindicies[..., 2] = int(MBindicies[..., 3])
    
    return Mindicies.squeeze()

    
    
def indicies3to4(Mindicies: np.typing.ArrayLike) -> np.ndarray:
    """
    Cnoverts Miller (hkl) indicies into Miller-Bravais (hkil) indicies 

    Parameters
    ----------
    Mindicies : np.typing.ArrayLike
        (...,3) array of Miller indicies.

    Raises
    ------
    ValueError
        If there are not 4 indicies.
        If the sum of the first three indicies is not 0.

    Returns
    -------
    MBindicies : np.array
        (...,4) array of Miller-Bravais indicies.

    """
    #check inputs are the right size
    Mindicies = np.atleast_1d(Mindicies)
    if Mindicies.shape[-1] != 3:
        raise ValueError("This function takes 3 value Miller indicies as inputs.")

    #do conversion
    MBindicies = np.zeros(Mindicies.shape[:-1]+ (4,), dtype=int)
    MBindicies[..., 0] = 1/3 * (2* Mindicies[..., 0] - Mindicies[..., 1])
    MBindicies[..., 1] = 1/3 * (2* Mindicies[..., 1] - Mindicies[..., 0])
    MBindicies[..., 2] = - (MBindicies[..., 0] + MBindicies[..., 1])
    MBindicies[..., 3] = Mindicies[..., 2]
    
    return MBindicies.squeeze()