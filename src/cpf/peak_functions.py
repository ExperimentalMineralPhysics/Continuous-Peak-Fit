#!/usr/bin/env python

"""
Functions describing the peak shape and properties.

Ideally if a new peak shape is needed it would be
possible to define it here and it would propagate through the code.
But this is not implemented and might not be possible...
"""

__all__ = [
    "peak_components",
    "expand_component_string",
    "gaussian_peak",
    "lorentzian_peak",
    "pseudo_voigt_peak",
]


__doc__ = "Functions for Pseudo-Voigt peak shape. "


import numpy as np

from cpf.util.logging import get_logger

logger = get_logger("cpf.peak_functions")



def peak_components(full=False, include_profile=True):
    """
    Lists the parameters needed for the peak shape function.
    If 'full' is True returns all the parameters needed for the fit (including background and symmetry)
    Otherwise returns just the peak parameters.
    

    Parameters
    ----------
    full : bool, optional
        Return just the parameters for the peak shape (if False) or with additional 
        parameters (if True). The default is False.
    include_profile : bool, optional
        Include the peak shape 'profile' in the output. The default is True.

    Returns
    -------
    comp_list : list
        List of peak component 1 letter abbreviations.
    comp_names : list
        List of peak component full names.

    """
    comp_list = ["h", "d", "w"]
    comp_names = ["height", "d-space", "width"]
    if include_profile:
        comp_list.extend("p")
        comp_names.extend(["profile"])
    if full:
        comp_list.append("bg")
        comp_names.append("background")
        comp_list.append("s")
        comp_names.append("symmetry")

    return comp_list, comp_names


def expand_component_string(comp):
    """
    Exapand he compnent name from its short representation to its full name.
    e.g. "d" --> "d-spance"
    
    Conversion is performed using lookup of outputs from peak_functions.peak_components()

    Parameters
    ----------
    comp : str
        Single letter peak shape component.

    Raises
    ------
    ValueError
        Unrecognised peak property type.

    Returns
    -------
    out : str
        Full name for peak shape component in input.

    """
    comp_list, comp_names = peak_components(full=True)
    if comp in comp_names:
        out = comp
    elif comp in comp_list:
        out = comp_names[comp_list.index(comp)]
    else:
        raise ValueError("Unrecognised peak property type")
    return out


# Gaussian shape
def gaussian_peak(two_theta, two_theta_0, w_all, h_all):
    """
    Calculates intensiities for Gaussian peak shape

    Parameters
    ----------
    two_theta : np.array
        Array of all two theat values to calculate intensity for.
    two_theta_0 : np.array
        Centroid of the peak at each two theta.
    w_all : np.array
        Width of the peak at each two theta value.
    h_all : np.array
        Height of the peak at each two theta value.

    Returns
    -------
    gauss_peak : np.array
        Intensity of the Gaussian peak at each two theta value.

    """
    w_all = w_all / np.sqrt(np.log(4))
    gauss_peak = h_all * np.exp((-((two_theta - two_theta_0) ** 2)) / (2 * w_all**2))
    return gauss_peak


def lorentzian_peak(two_theta, two_theta_0, w_all, h_all):
    """
    Calculates intensiities for Lorentz peak shape

    Parameters
    ----------
    two_theta : np.array
        Array of all two theat values to calculate intensity for.
    two_theta_0 : np.array
        Centroid of the peak at each two theta.
    w_all : np.array
        Width of the peak at each two theta value.
    h_all : np.array
        Height of the peak at each two theta value.

    Returns
    -------
    lorentz_peak : np.array
        Intensity of the Lorentzian peak at each two theta value.

    """
    lorentz_peak = h_all * w_all**2 / ((two_theta - two_theta_0) ** 2 + w_all**2)
    return lorentz_peak


def pseudo_voigt_peak(two_theta, two_theta_0, w_all, h_all, l_g_ratio):
    """
    Calculates intensiities for Lorentz peak shape

    Parameters
    ----------
    two_theta : np.array
        Array of all two theat values to calculate intensity for.
    two_theta_0 : np.array
        Centroid of the peak at each two theta.
    w_all : np.array
        Width of the peak at each two theta value.
    h_all : np.array
        Height of the peak at each two theta value.
    l_g_ratio : np.array
        Proportions of Gauss and Lorentz peak at each two theta value.
        if l_g_ratio = 1 returns Gaussian peak 
        if l_g_ratio = 0 returns Lorentzian peak 
        
    Returns
    -------
    PesudoVoigt_peak : np.array
        Intensity of the Pesudo-Voigt peak at each two theta value.

    """
    PesudoVoigt_peak = l_g_ratio * gaussian_peak(two_theta, two_theta_0, w_all, h_all) + (
        1 - l_g_ratio
    ) * lorentzian_peak(two_theta, two_theta_0, w_all, h_all)
    return PesudoVoigt_peak
