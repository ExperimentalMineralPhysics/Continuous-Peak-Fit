#!/usr/bin/env python3

"""
Functions for cpf describing the fourier or spline series used as parameter sets for the
peak shape and background properties.
"""

__all__ = [
    "coefficient_types",
    "coefficient_type_as_number",
    "coefficient_type_as_string",
    "get_series_type",
    "get_params_type",
    "get_number_coeff",
    "get_order_from_coeff",
    "get_order_from_params",
    "get_series_mean",
    "coefficient_expand",
    "spline_expand",
    "fourier_expand",
    "background_expansion",
]

import numpy as np
import numpy.ma as ma
from scipy.interpolate import CubicSpline, make_interp_spline

import cpf.peak_functions as pf
import cpf.series_constraints as sc
from cpf.util.logging import get_logger

logger = get_logger("cpf.series_functions")



def coefficient_types(full=False):
    """
    Defines the coefficient types possible and the numbers associated with them.

    The series types and their associated numbers are: 
         0: 'fourier'
         1: 'linear', 'spline_linear'
         2: 'quadratic', 'spline_quadratic'
         3: 'cubic', 'spline_cubic', 'spline', 'spline-cubic'
         4: 'linear_open', 'spline_linear_open'
         5: 'quadratic_open', 'spline_quadratic_open'
         6: 'cubic_open', 'spline_cubic_open', 'spline_open',
         7: 'independent'

    Parameters
    ----------
    full : bool, optional
        Return dirctionary with just series numbers (False) or with all ancillary 
        information for fitting (True). The default is False.

    Returns
    -------
    coeff_types : dict
        dictionary of key, value pairs.

    """
    coeff_types = {}
    # add fourier series type
    coeff_types = {"fourier": {"num": 0,  "expansion_function": "fourier_expand" },}   
    # add spline, linear, periodic series type
    coeff_types |= {"linear": {
            "num": 1,
            "expansion_function": "spline_expand",
            "boundary_conditions": "periodic",
            "spline_type": "linear",
            },}
    coeff_types |= {"spline_linear": coeff_types["linear"],}
    # add spline, quadratic, periodic series type
    coeff_types |= {"quadratic": {
            "num": 2,
            "expansion_function": "spline_expand",
            "boundary_conditions": "periodic",
            "spline_type": "quadratic",
            },}
    coeff_types |= {"spline_quadratic": coeff_types["quadratic"],}
    # add spline, cubic, periodic series type    
    coeff_types |= {"cubic": {
            "num": 3,
            "expansion_function": "spline_expand",
            "boundary_conditions": "periodic",
            "spline_type": "quadratic",
            },}
    coeff_types |= {"spline_cubic": coeff_types["cubic"],}
    coeff_types |= {"spline": coeff_types["cubic"],}
    coeff_types |= {"spline-cubic": coeff_types["cubic"],}
    
    # add spline, linear, open series type    
    coeff_types |= {"linear_open": {
            "num": 4,
            "expansion_function": "spline_expand",
            "boundary_conditions": "natural",
            "spline_type": "linear",
            },}
    coeff_types |= {"spline_linear_open": coeff_types["linear_open"],}
    # add spline, quardartic, open series type    
    coeff_types |= {"quadratic_open": {
            "num": 5,
            "expansion_function": "spline_expand",
            "boundary_conditions": "natural",
            "spline_type": "quadratic",
            },}
    coeff_types |= {"spline_quadratic_open": coeff_types["quadratic_open"],}
    # add spline, cubic, open series type    
    coeff_types |= {"cubic_open": {
            "num": 6,
            "expansion_function": "spline_expand",
            "boundary_conditions": "natural",
            "spline_type": "quadratic",
            },}
    coeff_types |= {"spline_cubic_open": coeff_types["cubic_open"],}
    coeff_types |= {"spline_open": coeff_types["cubic_open"],}    
    
    # add independent coefficients. 
    coeff_types |= {"independent": {
            "num": 7,
            "expansion_function": "spline_expand",
            "boundary_conditions": "natural",
            "spline_type": "independent",
            },}
        
    if full!=True:
        for i in range(len(coeff_types)):
            ky = list(coeff_types.keys())[i]
            coeff_types[ky] = coeff_types[ky]["num"]
    
    return coeff_types


def coefficient_type_as_number(series_type, return_error=1):
    """
    Returns serise type as a number from given string. 
    The key,value pairs are defined in series_functions.coefficient_types()

    Parameters
    ----------
    series_type : str
        series type name string.
    return_error : bool, optional
        return an error string, rather than raising an error. 

    Raises
    ------
    ValueError
        Unrecognised coefficient series type.

    Returns
    -------
    out : str or int
        Either: integer index for series type or,
            a string saying that that input string is not recognised.
    """
    types = coefficient_types()
    if series_type in types.keys():
        out = types[series_type]
    elif series_type in types.values():
        #if input is a recognised number then return a number. 
        out = series_type
    else:
        error_str = ("Unrecognised string index for series type. The valid options are "
            "defined in cpf.series_functions.coefficient_types()"
            )
        if return_error == 0:
            raise ValueError(error_str)
        else:
            out = error_str
    return out
    

def coefficient_type_as_string(series_type):
    """
    Returns series type as a string from given number. 
    The key,value pairs are defined in series_functions.coefficient_types()

    Parameters
    ----------
    series_type : float, int, 
        series type name string.

    Raises
    ------
    ValueError
        Unrecognised coefficient series type.

    Returns
    -------
    out : str
        Either: string index for series type or,
            a string saying that that input string is not recognised.
    """
    types = coefficient_types()
    if series_type in list(types.keys()):
        #if input is a recognised strong then return a strong.
        out = series_type
    elif series_type in list(types.values()):
        out = list(types.keys())[list(types.values()).index(series_type)]
    else:
        raise ValueError(
            "Unrecognised number index for series type. The valid options are "
            "defined in cpf.series_functions.coefficient_types()"
        )
    return out


def get_series_type(param, param_str, comp=None):
    """
    Get series type from parameter dictionary.

    Parameters
    ----------
    param : dict
        dictionary with multiple coefficients per component.
    param_str : str
        base string to select parameters.
    comp : str, optional
        component to add to base string to select parameters. The default is None.

    Returns
    -------
    out : list
        nested list of [parameters, errors] in alphanumerical order.

    """
    if comp is not None:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    new_str = new_str + "_tp"
    if isinstance(param, dict) and new_str in param:
        out = param[new_str]
    elif isinstance(param, dict) and new_str not in param:
        out = coefficient_types()["fourier"]
    else:
        out = None
    return out


def get_params_type(orders, comp, peak=0):
    """
    Get series type from input orders.
    
    Parameters
    ----------
    orders : dict
        order dictionary from cpf,Settings.settings or input.
    comp : str
        component string to select parameter. The default is None.
    peak : int, optional
        Choose which peak in orders to return series for. The default is 0.

    Returns
    -------
    series_type : string
        Series type for coefficients.

    """
    comp_str = pf.expand_component_string(comp) + "_type"
    if comp_str in orders:
        series_type = orders[comp_str]
    elif len(orders["peak"]) > peak and comp_str in orders["peak"][peak]:
        series_type = orders["peak"][peak][comp_str]
    else:
        series_type = "fourier"

    return series_type


def get_number_coeff(orders, comp, peak=0, azimuths=None):
    """
    Returns the expected number of coefficients from order value
    
    For a fourier series the number of coefficients is 2n+1.
    The same convention is adopted for spline series. 
    For 'independent' series, the number of unique azimuths is needed to determine 
    the number of coefficients. 
    
    Parameters
    ----------
    orders : dict
        order dictionary from cpf,Settings.settings or input.
    comp : str
        component string to select parameter. The default is None.
    peak : int, optional
        Choose which peak in orders to return series for. The default is 0.
    azimuths : np.array, optional
        array of unique azimuths. The default is None.

    Raises
    ------
    ValueError
        "Cannot define number of independent values without a number of coefficients."

    Returns
    -------
    n_param : int
        number of parameters in the series..

    """
    parm_str = get_params_type(orders, comp, peak)
    parm_num = coefficient_type_as_number(parm_str)

    # if parm_num == 5:  # independent
    if parm_num == coefficient_types(full=True)["independent"]["num"]:  # independent
        if azimuths is None:
            raise ValueError(
                "Cannot define number of independent values without a number of coefficients."
            )
        else:
            n_param = len(np.unique(azimuths[~ma.array(azimuths).mask]))

    elif comp == "bg" or comp == "background" or comp == "f":
        n_param = sc.BiggestValue(orders["background"][peak]) * 2 + 1

    else:  # everything else.
        n_param = (
            sc.BiggestValue(orders["peak"][peak][pf.expand_component_string(comp)]) * 2
            + 1
        )
    return n_param


def get_order_from_coeff(n_coeff, parm_num=0, azimuths=None):
    """
    Returns the order of a series from the number of coefficients it contains
    
    For a fourier series the order is (n-1)/2.
    The same convention is adopted for spline series. 
    For 'independent' series, the number of unique azimuths is needed to determine 
    the number of coefficients. 
    
    Parameters
    ----------
    n_coeff : num
        Number of cofficents in the series.
    parm_num : int or str, optional
        Label for series type - either a string or a numeric value, as defined in 
        series_functions.coefficient_types(). The default is 0.
    azimuths : np.array, optional
        array of unique azimuths. The default is None.

    Raises
    ------
    ValueError
        "Cannot define number of independent values without a number of coefficients."

    Returns
    -------
    order : int
        Order of the series.

    """
    parm_num = coefficient_type_as_number(parm_num)
    if parm_num == coefficient_types(full=True)["independent"]["num"]:  # independent
        if azimuths is None:
            raise ValueError(
                "Cannot define number of independent values without a number of coefficients."
            )
        else:
            order = azimuths.shape[0]
    else:  # everything else.
        order = (n_coeff - 1) / 2
    return order


def get_order_from_params(params, comp=None, peak=0):
    """
    Calculate the order of the series from a list of coefficients

    The order of the series is: (len(coefficients)-1) / 2
    
    Parameters
    ----------
    params : list
        List of series coefficients.
    comp : str
        component string to select parameter. The default is None.
    peak : int, optional
        Choose which peak in orders to return series for. The default is 0.

    Raises
    ------
    ValueError
        Params are not in a recongised type.

    Returns
    -------
    order : int
        Order fo the series.
        
    """
    # Given list of Fourier coefficients return order (n)
    if isinstance(params, (list,)):
        l = len(params)
    elif isinstance(params, (float,)):
        l = np.size(params)
    elif isinstance(params, (int,)):
        l = 1
    elif isinstance(params, (dict,)):
        parm_str = get_params_type(params, comp, peak)
        if comp == "bg" or comp == "background" or comp == "f":
            l = np.max(params["background"][peak])
        else:  # everything else.
            l = np.max(params["peak"][peak][pf.expand_component_string(comp)])
    else:
        logger.debug(" ".join(map(str, [("params", params)])))
        logger.debug(" ".join(map(str, [("type(params)", type(params))])))
        err_str = "Parameter list is not list, float or a dictionary."
        logger.critical(" ".join(map(str, [(err_str)])))
        raise ValueError(err_str)

    #convert from lenth to an order.
    order = get_order_from_coeff(l)

    return order


def get_series_mean(param, param_str, comp=None):
    """
    Calcualte the mean of the parameter series from parameter dictionary.

    Parameters
    ----------
    param : dict
        dict with multiple coefficients per component.
    param_str : str
        base string to select parameters.
    comp : str, optional
        component to add to base string to select parameters. The default is None.

    Returns
    -------
    mean : float
        weighted mean of the series.

    """
    # FIX ME: here we need to be able to discard the outliers.
    # We should use the medaian and the mean deviation from the median...
    if get_series_type(param, param_str, comp=comp) == coefficient_types()["fourier"]:
        # if it is a Fourier series just get the first value.
        mean = param[param_str + "_" + comp + "0"].value
    else:
        # get a mean of all the coefficients
        mean_tmp = []
        done = 0
        n = 0
        while done == 0:
            try:
                mean_tmp.append(param[param_str + "_" + comp + str(n)].value)
                n = n + 1
            except:
                # now we have run out of coefficients. So get the mean and then leave the loop.
                mean = np.mean(mean_tmp)
                done = 1

    return mean


def coefficient_expand(
    azimuth,
    param=None,
    coeff_type="fourier",
    comp_str=None,
    start_end=[0, 360],
    **params,
):
    """
    Calcuate the value of a series at each azimuth. 

    Parameters
    ----------
    azimuth : np.array
        array of unique azimuths.
    param : list, optional
        Series coefficients. The default is None.
    coeff_type : str or int, optional
        The type of series to be expanded. The default is "fourier".
    comp_str : str, optional
        String denoting which component of the peak is to be calculated. The default is None.
    start_end : list, optional
        Limits for the expnasion. The default is [0, 360].
    **params : dict
        parameter dictionary.

    Raises
    ------
    ValueError
        Unrecognised series type.

    Returns
    -------
    out : np.array
        Coefficient value at each azimuth.

    """
    series_name = coefficient_type_as_string(coeff_type)
    all_series = coefficient_types(full=True)
    
    # FIXME: this could be changed so that all_series[series_name]["expansion_function"]
    # is used with getattr -- allowing easier future expansion of the series types. 
    if all_series[series_name]["expansion_function"] == "fourier_expand":
        out = fourier_expand(azimuth, inp_param=param, comp_str=comp_str, **params)
        
    elif all_series[series_name]["expansion_function"] == "spline_expand":
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            bc_type=all_series[series_name]["boundary_conditions"],
            kind=all_series[series_name]["spline_type"],
            **params,
        )
    
    else:
        raise ValueError(
            "Unrecognised number index for series type. The valid options are "
            "defined in cpf.series_functions.coefficient_types()"
        )
        
    return out


def spline_expand(
    azimuth,
    inp_param=None,
    comp_str=None,
    start_end=[0, 360],
    bc_type="periodic",
    kind="cubic",
    **params,
):
    """
    Calculate series value at each azimuth for given spline coefficients

    Parameters
    ----------
    azimuths : np.array
        array of unique azimuths.
    inp_param : float, list, optional
        list of values at spline tie points or number of spline tie points. The default is None.
    comp_str : str, optional
        str to determine which coefficients to use from params. The default is None.
    start_end : list, optional
        Minimum and maximum azimuth. The default is [0, 360].
    bc_type : str, optional
        Bounding conditions type: Options are "indepeddnt", "periodic", and "natural". 
        The default is "periodic".
    kind : str, optional
        Type of spline or independent series. The default is "cubic".
    **params : dict
        lmfit dict of coefficients as parameters.

    Returns
    -------
    np.array
        series value for each azimuth.

    """
    if inp_param is not None:
        if not isinstance(inp_param, np.float64):
            inp_param = np.array(inp_param)
            if len(inp_param.shape) > 1:
                if np.any(np.array(inp_param.shape) > 1):
                    inp_param = np.squeeze(inp_param)
                elif np.all(np.array(inp_param.shape) == 1):
                    inp_param = np.squeeze(inp_param)
                    inp_param = np.array([inp_param], float)
    else:
        # create relevant list of parameters from dict
        str_keys = [
            key for key, val in params.items() if comp_str in key and "tp" not in key
        ]
        inp_param = []
        for j in range(len(str_keys)):
            inp_param.append(params[comp_str + str(j)])
    
    if kind == "independent":
        points = np.unique(azimuth)
    elif bc_type == "periodic":
        points = np.linspace(start_end[0], start_end[1], np.size(inp_param) + 1)
        inp_param = np.append(inp_param, inp_param[0])
    elif isinstance(bc_type, (list, tuple, np.ndarray)):
        points = bc_type
    else:
        points = np.linspace(start_end[0], start_end[1], np.size(inp_param))

    if kind == "cubic":
        k = 3
    elif kind == "quadratic":
        k = 2
    elif kind == "linear" or kind == "independent":
        k = 1
    else:
        raise ValueError("Unknown spline type.")

    fout = np.ones(azimuth.shape)

    if (
        azimuth.size == 1
    ):  # this line is required to catch error when out is single number.
        try:
            fout = inp_param[0]
        except IndexError:
            fout = inp_param
    else:
        fout[:] = inp_param[0]
    # essentially d_0, h_0 or w_0
    if not isinstance(inp_param, np.float64) and np.size(inp_param) > 1:
        if k == 3:
            spl = CubicSpline(
                points, inp_param, bc_type=bc_type, extrapolate="periodic"
            )
        else:
            spl = make_interp_spline(points, inp_param, k=k)

        fout = spl(azimuth)

    return np.squeeze(fout)


def fourier_expand(
    azimuth, inp_param=None, comp_str=None, start_end=[0, 360], **params
):
    """
    Calculate series value at each azimuth for given fourier coefficients

    Parameters
    ----------
    azimuth : np.array
        array of unique azimuths.
    inp_param : float, optional
        list of Fourier coefficients. The default is None.
    comp_str : str, optional
        str to determine which coefficients to use from params. The default is None.
    start_end : list, optional
        Minimum and maximum azimuth. The default is [0, 360].
    **params : dict
        lmfit dict of coefficients as parameters.

    Returns
    -------
    np.array
        series value for each azimuth.

    """
    if inp_param is not None:
        # FIX ME: Need to check the fourier is a licit length
        if not isinstance(inp_param, np.float64):
            inp_param = np.array(inp_param)
            if len(inp_param.shape) > 1:
                if np.any(np.array(inp_param.shape) > 1):
                    inp_param = np.squeeze(inp_param)
                elif np.all(np.array(inp_param.shape) == 1):
                    inp_param = np.squeeze(inp_param)
                    inp_param = np.array([inp_param], float)
    else:
        # create relevant list of parameters from dict
        str_keys = [
            key for key, val in params.items() if comp_str in key and "tp" not in key
        ]
        inp_param = []
        for j in range(len(str_keys)):
            inp_param.append(params[comp_str + str(j)])
    fout = np.ones(azimuth.shape)
    # this line is required to catch error when out is single number.
    if azimuth.size == 1:
        try:
            fout = inp_param[0]
        except IndexError:
            fout = inp_param
    else:
        fout[:] = inp_param[0]
    # essentially d_0, h_0 or w_0

    azm_tmp = np.deg2rad(
        (azimuth - start_end[0]) / (start_end[-1] - start_end[0]) * 360
    )
    if not isinstance(inp_param, np.float64) and np.size(inp_param) > 1:
        for i in range(1, int((len(inp_param) - 1) / 2) + 1):
            # len(param)-1 should never be odd because of initial a_0 parameter
            # try:
            # try/except is a ctch for expanding a fourier series that has failed and has nones as coefficient values.
            # azm_tmp stretches the azimuths between the max and min of start_end. In XRD cases this should have no effect
            # but is added for consistency with the spline functions
            fout = (
                fout
                + inp_param[(2 * i) - 1] * np.sin((azm_tmp) * i)
                + inp_param[2 * i] * np.cos((azm_tmp) * i)
            )
    return np.squeeze(fout)


def background_expansion(azimuth_two_theta, orders, params):
    """
    Calculate background value at each azimuth / two theta pair for given series 
    coefficients.

    Parameters
    ----------
    azimuth_two_theta : list of np.array
        list of arrays for azimuth and two theta.
    orders : dict
        order dictionary from cpf,Settings.settings or input.
    params : list
        list of input parameters.

    Returns
    -------
    bg_all : np.array
        Background intensity at each azimuth / two theta pair.

    """
    azimuth, two_theta = azimuth_two_theta
    two_theta_prime = two_theta - orders["range"][0]

    # recreate background array to pass to coefficient_expand
    background_keys = [
        key for key, val in params.items() if "bg_c" in key and "tp" not in key
    ]
    n_term_fourier = []
    i = 0
    backg = []
    backg_tp = []
    while sum(n_term_fourier) < len(background_keys):
        f = sum("bg_c" + str(i) in L for L in background_keys)
        n_term_fourier.append(f)
        fourier_background = []
        for k in range(f):
            fourier_background.append(params["bg_c" + str(i) + "_f" + str(k)])
            if "bg_c" + str(i) + "_f_tp" in params:
                b_tp = params["bg_c" + str(i) + "_f_tp"]
            else:
                b_tp = coefficient_types()["fourier"]

        backg.append(np.array(fourier_background))
        backg_tp.append(b_tp)
        i = i + 1

    bg_all = np.zeros(azimuth.shape)
    for i in range(len(backg)):
        out = coefficient_expand(azimuth, backg[i], backg_tp[i])
        bg_all = bg_all + (out * (two_theta_prime ** float(i)))
    return bg_all
