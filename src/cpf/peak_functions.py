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
import uncertainties.unumpy as unp

from cpf.util.logging import get_logger

logger = get_logger("cpf.peak_functions")



def peak_components(full=False, include_profile=True, include_combined=False):
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
    if include_combined:
        comp_list.append("a")
        comp_names.append("area")

    return comp_list, comp_names


def expand_component_string(comp):
    """
    Exapand he compnent name from its short representation to its full name.
    e.g. "d" --> "d-space"
    
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
    comp_list, comp_names = peak_components(full=True,include_combined=True)
    if comp in comp_names:
        out = comp
    elif comp in comp_list:
        out = comp_names[comp_list.index(comp)]
    else:
        raise ValueError("Unrecognised peak property type")
    return out


def compress_component_string(comp):
    """
    Compress the compnent name from its long name to its short representation.
    e.g.  "d-space" --> "d"
    
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
    comp_list, comp_names = peak_components(full=True,include_combined=True)
    if comp in comp_list:
        out = comp
    elif comp in comp_names:
        out = comp_list[comp_names.index(comp)]
    else:
        raise ValueError("Unrecognised peak property type")
    return out



# Gaussian shape
def gaussian_peak(two_theta, two_theta_0, w_all, h_all):
    """
    Calculates intensiities (height) at each position for Gaussian peak shape

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
    Calculates intensity (height) at each position for Lorentz peak shape

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
    Calculates intensities (heights) of Pseudo-Voigt peak shape

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
        height of the Pesudo-Voigt peak at each two theta value.

    """
    PesudoVoigt_peak = l_g_ratio * gaussian_peak(two_theta, two_theta_0, w_all, h_all) + (
        1 - l_g_ratio
    ) * lorentzian_peak(two_theta, two_theta_0, w_all, h_all)
    return PesudoVoigt_peak



def area(w_all, h_all, l_g_ratio):
    """
    Calculates area intensiities for Pseudo-Voigt peak shape

    Parameters
    ----------
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
        Area of the Pesudo-Voigt peak at each two theta value.

    """
    # Gauss sum
    sumG = h_all * unp.sqrt(np.pi*(2 * (w_all/unp.sqrt(unp.log(4)))**2))
    
    # lotentz sum
    sumL = h_all * np.pi /unp.sqrt(1/w_all**2)
    # lorentz sum for h_all=1 and w_all = i converges on pi at infinity.
    # we just assume this here as it is simplest but should perhaps have a cut off     
    
    return l_g_ratio * sumG + (1 - l_g_ratio) * sumL
    

def plot_sum():
    """
    Plot shape functions for the peak profiles. 
    
    Calculates some sums to show self.area is correct.
    
    Returns
    -------
    None.

    """
    # make graph showing the peak functions
    x = np.linspace(-5,5,500)
    
    h=1
    w=1
    
    y_g = pseudo_voigt_peak(x,0,w,h,1)
    y_l = pseudo_voigt_peak(x,0,w,h,0)
    
    import matplotlib.pyplot as plt
    plt.plot(x,y_g,'-r', label = "Gaussian (profile=1)")
    plt.plot(x, y_l,'-b', label = "Lorentz (profile=0)")
    plt.legend()
    plt.title(f"cpf peak shapes: height={h}, width={w}$")
    
        
    # proof of sum functions.
    h = 3
    w = .234
    
    a = np.linspace(-10000,10000,500000)
    y_l = pseudo_voigt_peak(a,0,w,h,0)
    y_g = pseudo_voigt_peak(a,0,w,h,1)
    y_mix = pseudo_voigt_peak(a,0,w,h,0.5)
    
    L_sum = np.sum(y_l)*(a[1]-a[0])
    G_sum = np.sum(y_g)*(a[1]-a[0])
    mix_sum = np.sum(y_mix)*(a[1]-a[0])
    
    G_sum2 = area(w, h, 1)
    L_sum2 = area(w, h, 0)
    mix_sum2 = area(w, h, 0.5)
    
    print("Gaussian")
    print(f"Numerical area: {G_sum}")
    print(f"arithmetic area: {G_sum2}")
    print(f"difference = {G_sum2-G_sum}")
    
    print("Lorentz")
    print(f"Numerical area: {L_sum}")
    print(f"arithmetic area: {L_sum2}")
    print(f"difference = {L_sum2-L_sum}")
    
    
    print("mixed")
    print(f"Numerical area: {mix_sum}")
    print(f"arithmetic area: {mix_sum2}")
    print(f"difference = {mix_sum2-mix_sum}")
    
    
    
    
    
    
    
    
    
    
    
    