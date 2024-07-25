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

import numpy as np
# from cpf.XRD_FitPattern import logger
from cpf.logger_functions import logger


def peak_components(full=False, include_profile=True):
    """
    Lists the parameters needed for each peak.
    If full is True returns all the parameters needed for the fit (including background and symetry)
    Otherwise returns just the peak parameters

    :param full:
    :return comp_list:
    :return comp_names:
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
    Exapand he compnent name from its short representation to its long name.
    e.g. "d" --> "d-spance"

    :param component name:
    :return: expanded component name
    """
    comp_list, comp_names = peak_components(full=True)
    if comp in comp_names:
        out = comp
    elif comp in comp_list:
        out = comp_names[comp_list.index(comp)]
    else:
        raise ValueError("Unrecognised peak property type")

    # if comp == "d" or comp == "d-space":
    #     out = "d-space"
    # elif comp == "h" or comp == "height":
    #     out = "height"
    # elif comp == "w" or comp == "width":
    #     out = "width"
    # elif comp == "p" or comp == "profile":
    #     out = "profile"
    # elif comp == "bg" or "f" or comp == "background":
    #     out = "background"
    # else:
    #     raise ValueError("Unrecognised peak property type")
    return out


# Gaussian shape
def gaussian_peak(two_theta, two_theta_0, w_all, h_all):
    """
    Gaussian shape
    :param two_theta:
    :param two_theta_0:
    :param w_all:
    :param h_all:
    :return:
    """
    w_all = w_all / np.sqrt(np.log(4))
    gauss_peak = h_all * np.exp((-((two_theta - two_theta_0) ** 2)) / (2 * w_all**2))
    return gauss_peak


def lorentzian_peak(two_theta, two_theta_0, w_all, h_all):
    """
    Lorentz shape
    :param two_theta:
    :param two_theta_0:
    :param w_all:
    :param h_all:
    :return:
    """
    l_peak = h_all * w_all**2 / ((two_theta - two_theta_0) ** 2 + w_all**2)
    return l_peak


def pseudo_voigt_peak(two_theta, two_theta_0, w_all, h_all, l_g_ratio):
    """
    Pseudo-Voigt
    :param two_theta:
    :param two_theta_0:
    :param w_all:
    :param h_all:
    :param l_g_ratio:
    :return:
    """
    p_v_peak = l_g_ratio * gaussian_peak(two_theta, two_theta_0, w_all, h_all) + (
        1 - l_g_ratio
    ) * lorentzian_peak(two_theta, two_theta_0, w_all, h_all)
    return p_v_peak
