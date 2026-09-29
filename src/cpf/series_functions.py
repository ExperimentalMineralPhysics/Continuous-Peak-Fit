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

import re
import numpy as np
import numpy.ma as ma
from scipy.interpolate import make_interp_spline, PPoly

import uncertainties.unumpy as unp
from uncertainties import UFloat, ufloat
 
import matplotlib.pyplot as plt
from lmfit import Parameters
    
import cpf.peak_functions as pf
import cpf.series_constraints as sc
from cpf.lmfit_model import coefficient_fit, initiate_params, gather_param_errs_to_list
from cpf.util.io import replace_value
from cpf.util.logging import get_logger

logger = get_logger("cpf.series_functions")



# TO DO: 
# 1. rename coefficient_expand as series_expand.
# 2. reorder the coefficient types so that the higher orderones have a greater numerical value.
#    then I can use the highest value as a check when combining the series. 


def coefficient_types(full=False):
    """
    Defines the coefficient types possible and the numbers associated with them.

    The series types and their associated numbers are:
         0: 'fourier'
         1: 'spline_linear', 'linear'
         2: 'spline_quadratic', 'quadratic'
         3: 'spline_cubic', 'cubic', 'spline', 'spline-cubic'
         4: 'spline_linear_open', 'linear_open'
         5: 'spline_quadratic_open', 'quadratic_open'
         6: 'spline_cubic_open', 'cubic_open', 'spline_open'
         7: 'independent'
    When cycled through the code the first value in each of these options is the one that
    will be added to the json files. So 'linear' will be written as 'spline_linear'

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
    coeff_types = {
        "fourier": {"num": 0, "expansion_function": "fourier_expand"},
    }
    # add spline, linear, periodic series type
    coeff_types |= {
        "spline_linear": {
            "num": 1,
            "expansion_function": "spline_expand",
            "boundary_conditions": "periodic",
            "spline_type": "linear",
        },
    }
    coeff_types |= {
        "linear": coeff_types["spline_linear"],
    }
    # add spline, quadratic, periodic series type
    coeff_types |= {
        "spline_quadratic": {
            "num": 2,
            "expansion_function": "spline_expand",
            "boundary_conditions": "periodic",
            "spline_type": "quadratic",
        },
    }
    coeff_types |= {
        "quadratic": coeff_types["spline_quadratic"],
    }
    # add spline, cubic, periodic series type
    coeff_types |= {
        "spline_cubic": {
            "num": 3,
            "expansion_function": "spline_expand",
            "boundary_conditions": "periodic",
            "spline_type": "cubic",
        },
    }
    coeff_types |= {
        "cubic": coeff_types["spline_cubic"],
    }
    coeff_types |= {
        "spline": coeff_types["spline_cubic"],
    }
    coeff_types |= {
        "spline-cubic": coeff_types["spline_cubic"],
    }

    # add spline, linear, open series type
    coeff_types |= {
        "spline_linear_open": {
            "num": 4,
            "expansion_function": "spline_expand",
            "boundary_conditions": "not-a-knot",
            "spline_type": "linear",
        },
    }
    coeff_types |= {
        "linear_open": coeff_types["spline_linear_open"],
    }
    # add spline, quardartic, open series type
    coeff_types |= {
        "spline_quadratic_open": {
            "num": 5,
            "expansion_function": "spline_expand",
            "boundary_conditions": "not-a-knot",
            "spline_type": "quadratic",
        },
    }
    coeff_types |= {
        "quadratic_open": coeff_types["spline_quadratic_open"],
    }
    # add spline, cubic, open series type
    coeff_types |= {
        "spline_cubic_open": {
            "num": 6,
            "expansion_function": "spline_expand",
            "boundary_conditions": "not-a-knot",
            "spline_type": "cubic",
        },
    }
    coeff_types |= {
        "cubic_open": coeff_types["spline_cubic_open"],
    }
    coeff_types |= {
        "spline_open": coeff_types["spline_cubic_open"],
    }

    # add independent coefficients.
    coeff_types |= {
        "independent": {
            "num": 7,
            "expansion_function": "spline_expand",
            "boundary_conditions": "natural",
            "spline_type": "independent",
        },
    }

    if full != True:
        for i in range(len(coeff_types)):
            ky = list(coeff_types.keys())[i]
            coeff_types[ky] = coeff_types[ky]["num"]

    return coeff_types


def coefficient_type_as_number(series_type, return_error=1):
    """
    Returns series type as a number from given string.
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
        # if input is a recognised number then return a number.
        out = series_type
    else:
        error_str = (
            "Unrecognised string index for series type. The valid options are "
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
        # if input is a recognised strong then return a strong.
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
            n_param = len(ma.unique(azimuths).compressed())

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
    
        if np.int64(order) == order:
            order = np.int64(order)
        else:
            raise ValueError(
                "Cannot have a non-integer order value. The number of coefficients provided must be incorrect"
            )
        
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

    # convert from lenth to an order.
    order = get_order_from_coeff(l)

    return order


def get_series_mean(param, param_str, comp=None):
    """
    Calcualte the mean of the parameter series from lmfit parameter dictionary or
    settings coefficient dictionary.

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
    if param_str in param:
        # then this is not a lmfit parameter dictionary but a coefficient dictionary.
        if param_str + "_type" in param:  #
            series_type = param[param_str + "_type"]
        else:
            series_type = "fourier"
        if series_type == "fourier":
            mean = param[param_str][0]
        else:
            # get a mean of all the coefficients
            mean = np.nanmean(param[param_str])

    else:
        if (
            get_series_type(param, param_str, comp=comp)
            == coefficient_types()["fourier"]
        ):
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
                    done = 1
            # now we have run out of coefficients. So get the mean and then leave the loop.
            mean = np.mean(mean_tmp)

    return mean


def coefficient_expand(
    azimuth,
    param=None,
    coeff_type="fourier",
    comp_str=None,
    start_end=[0, 360],
    no_negatives = False,
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
        out = fourier_expand(azimuth, 
                             inp_param=param, 
                             comp_str=comp_str, 
                             no_negatives=no_negatives,
                             **params)

    elif all_series[series_name]["expansion_function"] == "spline_expand":
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            no_negatives=no_negatives,
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


def construct_spline(
        coeffs,
        azimuth = None,
        start_end=[0, 360],
        bc_type="periodic",
        kind=None,
        **params):
    """
    Constructs a scipy B-spline from the paratmers given. 
    
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
    spl : scipy B-spline
        spline defined by the inputs.
    
    """
    
    if kind == "independent":
        if ma.isMaskedArray(azimuth):
            points = ma.unique(azimuth).compressed()
        else:
            points = np.unique(azimuth)
    elif bc_type == "periodic":
        points = np.linspace(start_end[0], start_end[1], np.size(coeffs) + 1)
        coeffs = np.append(coeffs, coeffs[0])
    elif isinstance(bc_type, (list, tuple, np.ndarray)):
        points = bc_type
    else:
        points = np.linspace(start_end[0], start_end[1], np.size(coeffs))
        
    if kind == "cubic":
        k = 3
    elif kind == "quadratic":
        k = 2
    elif kind == "linear" or kind == "independent":
        k = 1
        bc_type = None  # catch error feeding into make_interp_spline
    else:
        raise ValueError("Unknown spline type.")   

    if k >= len(points):
        # catch if the spline is underconstrained
        k = len(points) - 1
        if k < 0:
            k=0
            
    spl = make_interp_spline(points, unp.nominal_values(coeffs), k=k, bc_type=bc_type )   
    
    return spl, points


def spline_expand(
    azimuth,
    inp_param=None,
    comp_str=None,
    start_end=[0, 360],
    bc_type="periodic",
    kind=None,
    no_negatives=False,
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

    if (
        azimuth.size == 1
    ):  # this line is required to catch error when out is single number.
        try:
            fout = inp_param[0]
        except IndexError:
            fout = inp_param
    else:
        # fout[:] = inp_param[0]
        fout = np.ones(azimuth.shape) * inp_param[0]
    # essentially d_0, h_0 or w_0
    if not isinstance(inp_param, np.float64) and np.size(inp_param) > 1:
        spl, _ = construct_spline(
                coeffs=inp_param,
                azimuth = azimuth,
                start_end=start_end,
                bc_type=bc_type,
                kind=kind,
                **params)
        fout = spl(azimuth)

    if no_negatives and np.any(fout<0):
        fout[fout<0] = np.finfo(fout.dtype).eps
        
    if isinstance(inp_param[0], UFloat):
        # then the input is an array of values with errors. 
        # these errors will be greater than the formal errors on any fit.
        # Used when calculating combined series.
        if kind == "independent":
            inp = unp.std_devs(inp_param)
        else:
            # run to end-1 because have to cut value added by spline_expand.
            inp = unp.std_devs(inp_param)#[:-1]
        if no_negatives:
            # should prevent negative errors of itself
            kind = "linear"
        errs = spline_expand(
            azimuth,
            inp_param=inp,
            comp_str=None,
            start_end=start_end,
            bc_type=bc_type,
            kind=kind,
            no_negatives=True,
            **params,
        )
        # prevent negative errors
        if np.any(errs<0):
            errs[errs<0]=np.min(np.array([fout[errs<0], np.abs(errs[errs<0])]), axis=0)
        
        fout = unp.uarray(unp.nominal_values(fout), errs)
        
    return np.squeeze(fout)


def fourier_expand(
    azimuth, inp_param=None, comp_str=None, start_end=[0, 360], no_negatives=True, **params
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
            key for key, val in params.items() if comp_str in key and "tp" not in key and "err" not in key
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
        # fout[:] = inp_param[0]
        fout = np.ones(azimuth.shape) * inp_param[0]
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
    if no_negatives and np.any(fout<0):
        fout[fout<0] = np.finfo(fout.dtype).eps
    return np.squeeze(fout)


def background_expansion(azimuth_two_theta, orders, params,
                         no_negatives=False):
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
        out = coefficient_expand(azimuth, backg[i], backg_tp[i], 
                                 no_negatives = no_negatives)
        bg_all = bg_all + (out * (two_theta_prime ** float(i)))
    return bg_all

    
def combine_series(
        param_dict,
        azimuth = None,
        start_end=[0, 360],
        **kwargs
    ):
    """
    Numerically combine series from peak fits and return values for each Azimuth. 
    
    The default combined series is peak_functions.area(). But alternatives can be set 
    via kwargs if desired.

    Parameters
    ----------
    param_dict : dict
        dictionary of fitted peak parameters, with peak series values, errors and types in
    azimuth : list, np.array, optional
        azimuths to calculate the series at. The default is None.
    start_end : list, numpy.ndarray, optional
        start and end azimuths to combine series over. The default is [0, 360].

    Keyword Parameters
    ------------------
    combined_series_name : string
        Name of the combined series, and the attribute name for the function.
        Default is "area"
    combination_function : method
        Method of peak_functions that performs the combination/integration of the data.
        Default is "getattr(pf, combined_series_name)"    
    num_azimuths : int
        Number of azimuths to calculate peak at if azimuth is not set. 
        Default is 360

    Returns
    -------
    combined : numpy.ndarray
        Value of combined series at each given Azimuth.
    """
    
    # get kwargs that are needed
    combined_series_name = kwargs.get("combined_series_name", "area")
    combination_function = kwargs.get("combination_function", getattr(pf, combined_series_name))
    num_azimuths = kwargs.get("num_azimuths", 360)

    if azimuth is None:
        azimuth = np.linspace(start_end[0], start_end[1], num_azimuths)
        
    # expand series around the azimuth values
    h = unp.uarray(replace_value(param_dict["height"]), replace_value(param_dict["height_err"]))
    height = coefficient_expand(azimuth, 
                              param=h, 
                              coeff_type=param_dict["height_type"],
                              comp_str="height",
                              start_end=start_end)
    w = unp.uarray(replace_value(param_dict["width"]), replace_value(param_dict["width_err"]))
    width = coefficient_expand(azimuth, 
                              param=w, 
                              coeff_type=param_dict["width_type"],
                              comp_str="width",
                              start_end=start_end)
    p = unp.uarray(replace_value(param_dict["profile"]), replace_value(param_dict["profile_err"]))
    profile = coefficient_expand(azimuth, 
                              param=p,
                              coeff_type=param_dict["profile_type"],
                              comp_str="profile",
                              start_end=start_end)

    # get combined values.    
    combined = combination_function(width, height, profile)
    
    if 0:
        import matplotlib.pyplot as plt
        plt.plot(azimuth, combined, '.')
        
    return combined
    

def get_combined_series(
        param_dict,
        azimuth = None,
        start_end=[0, 360],
        **kwargs,
    ):
    """
    Combines series into new series. The default is for area of the peak, but custom integrations 
    can be specififed in kwargs.
    
    The integreated series is a spline type series and has an order that is at least equal to the 
    sum of the originating series' orders.
    
    Thie combining of the series is done numerically by calling peak_functions.integrated (or defined function)
    at the spline knot points of the combined series. This works in testing for both the values and their errors.
    
    An old algorthim that fits the new series to an expanded set of values is still in the code but disabled by default.
    If needed it can be acces by adding "fit_series_new" = True to the kwargs.

    Parameters
    ----------
    param_dict : dict
        dictionary of fitted peak parameters, with peak series values, errors and types in
    azimuth : list, np.array, optional
        azimuths to calculate the series at. The default is None.
    start_end : list, numpy.ndarray, optional
        start and end azimuths to combine series over. The default is [0, 360].

    Keyword Parameters
    ------------------
    combined_series_name : string
        Name of the combined series, and the attribute name for the function.
        Default is "area"
    combination_function : method
        Method of peak_functions that performs the combination/integration of the data.
        Default is "getattr(pf, combined_series_name)"    
    order_scale : int
        Value to multiply the order of the combined series by. Testing shows that 
        the errors sometimes need a factor of two to be able to replicate the shape of the errors 
        proucted by combining the series.
        Default is 2. 
    fit_series_new" : bool
        Switch to access an old algorthim. Default is False

    Returns
    -------
    combined_series : dict
        dictionary of combined series parameters, with peak series values, errors and type in.

    """
    
    # get kwargs that are needed
    combined_series_name = kwargs.get("combined_series_name", "area")
    order_scale = kwargs.get("order_scale", 2)
        
    # type of new series.
    # use maximum of series types as numbers
    i_types = [coefficient_type_as_number(param_dict["width_type"]),
                    coefficient_type_as_number(param_dict["height_type"]),
                    coefficient_type_as_number(param_dict["profile_type"])
                    ]
    i_type = np.max(i_types)
    if i_type == 0:
        i_type = 3
        logger.moreinfo("Change combined series type from 'fourier' to 'cubic spline' because fouriers do not propagate series errors correctly.")
    i_bc = coefficient_types(full=True)[coefficient_type_as_string(i_type)]['boundary_conditions']
    
    # get order and knot points of new series
    if i_type == coefficient_types()["independent"]:
        # then independent and new series needs same order as old because order 
        # is the number of independent values.
        # In this case none of the values can be larger than number of values so use max
        i_order = np.max([len(param_dict["width"]),
                        len(param_dict["height"]),
                        len(param_dict["profile"])
                        ])
        
        # make sure only have unique azimuths
        azimuth = np.unique(azimuth)
    elif (any(x == coefficient_types()["spline_linear"] for x in i_types) or
        any(x == coefficient_types()["linear_open"] for x in i_types)):
        # if any of the series are line splines then we need to make sure the knots in the combined series 
        #valign with those in the base series. --> use lowest common multiple as test of this.
        i_order = get_order_from_coeff(
                        np.lcm.reduce([(2*len(param_dict["width"])+1 if (coefficient_type_as_number(param_dict["width_type"])==0 ) else len(param_dict["width"])),
                                      (2*len(param_dict["height"])+1 if (coefficient_type_as_number(param_dict["height_type"])==0 ) else len(param_dict["height"])),
                                      (2*len(param_dict["profile"])+1 if (coefficient_type_as_number(param_dict["profile_type"])==0 ) else len(param_dict["profile"]))
                                      ])
                    )  #  * order_scale
    else:
        # for Fourier series, order of combined series is sum of orders -- 
        # i.e. sin(x) * sin(x) has order sin^2(x).
        # but this only applies if order is greater than 0.

        # because series type is forced to be spline, we have to double any series that is a fourier.
        i_order = (get_order_from_coeff(len(param_dict["width"]), coefficient_type_as_number(param_dict["width_type"])) 
                           * (2 if (coefficient_type_as_number(param_dict["width_type"])==0 and len(param_dict["width"]) > 1) else 1) + 
                   get_order_from_coeff(len(param_dict["height"]), coefficient_type_as_number(param_dict["height_type"]))
                           * (2 if (coefficient_type_as_number(param_dict["height_type"])==0 and len(param_dict["height"]) > 1) else 1) + 
                   get_order_from_coeff(len(param_dict["profile"]), coefficient_type_as_number(param_dict["profile_type"]))
                           * (2 if (coefficient_type_as_number(param_dict["profile_type"])==0 and len(param_dict["profile"]) > 1) else 1)
                        ) * order_scale
    
    
    # FIX ME : propagate symmetry properly!!!
    if "symmetry" in param_dict:
        symmetry = param_dict["symmetry"]
    else:
        symmetry = 1
    symmetry = 1
    # force symmetry to be 1 because all we care about here is getting the new series to align with the old ones. The 
    # axzimuths that this happens at are mutable. 
        

    if kwargs.get("fit_series_new", False):
        # fit combined series with new series; use series fitting code used in 
        # data fitting. 
                
        num_azimuths = kwargs.get("num_azimuths", 360)
        if i_order*2+1 > num_azimuths:
            num_azimuths = int(np.ceil( (i_order*2+1) /num_azimuths))*num_azimuths
        azimuth = np.linspace(start_end[0], start_end[1], num_azimuths)
        
        
        combined = combine_series(param_dict, azimuth=azimuth, start_end=start_end, **kwargs)
    
        master_params = Parameters()  
        param_str = "peak_0"
        component = pf.compress_component_string(combined_series_name)
        if i_type == coefficient_types()["independent"]:
            # independent so have values from combined series. 
            # parse into fit-like dictionary.
            combined_series = {}
            combined_series[combined_series_name] = list(unp.nominal_values(combined))
            combined_series[combined_series_name+"_err"] = list(unp.std_devs(combined))
            combined_series[combined_series_name+"_type"] = coefficient_type_as_string(i_type)
        elif np.all(unp.nominal_values(combined)==0):
            # all the values are zeros
            combined_series = {}
            combined_series[combined_series_name] = [0] * get_number_coeff({"peak": [{"area": i_order}]},"area")
            combined_series[combined_series_name+"_err"] = [0] * get_number_coeff({"peak": [{"area": i_order}]},"area")
            combined_series[combined_series_name+"_type"] = coefficient_type_as_string(i_type)
        
        else:
            master_params = initiate_params(
                master_params,
                param_str,
                component,
                coeff_type=i_type,
                trig_orders=i_order,
                limits=None,
                value=None,
                types=True,
            )
            fout = coefficient_fit(
                azimuth=azimuth,
                ydata=unp.nominal_values(combined),
                inp_param=master_params,
                param_str=param_str + "_" + component,
                symmetry=symmetry,
                errs=unp.std_devs(combined),
                fit_method="leastsq",
                start_end = start_end
                )
            
            if logger.is_below_level(level="DEBUG"):
                fout.plot(show_init=True)
                fout.params.pretty_print()
            """
            get errors on progagated series.
            --------------------------------
            Fitting the combined values with a series reproduces the centroid values. 
            But the errors are too small so fit the errors with a series to get the expected coefficient errors.
            """
            master_params = Parameters()  
            master_params = initiate_params(
                master_params,
                param_str,
                component,
                coeff_type=i_type,
                trig_orders=i_order,
                limits=[0, np.max(unp.std_devs(combined))],
                value=None,
                types=True,
            )
            ferrs_out = coefficient_fit(
                azimuth=azimuth,
                ydata=unp.std_devs(combined),
                inp_param=master_params,
                param_str=param_str + "_" + component,
                symmetry=symmetry,
                errs=None,#unp.std_devs(combined)*0 + 1E-6, # np.array(data_val_errors),
                fit_method="leastsq",
                start_end = start_end
            )
            combined_series = {}
            combined_series[combined_series_name] = gather_param_errs_to_list(
                                                    fout.params, "peak_0", comp=component
                                                )[0]
            combined_series[combined_series_name+"_err"] = gather_param_errs_to_list(
                                                    ferrs_out.params, "peak_0", comp=component
                                                )[0]
            combined_series[combined_series_name+"_type"] = coefficient_type_as_string(i_type)
            if logger.is_below_level(level="DEBUG"):
                """
                test the new series can reproduce the errors in the original data.
                """
                
                if i_bc == "periodic":
                    azimuth2 = np.linspace(start_end[0], start_end[1], len(combined_series[combined_series_name])+1)
                    azimuth2 = azimuth2[:-1]
                else:
                    azimuth2 = np.linspace(start_end[0], start_end[1], len(combined_series[combined_series_name]))
                
                i = unp.uarray(combined_series[combined_series_name], combined_series[combined_series_name+"_err"])
                reconstructed_series = coefficient_expand(azimuth*symmetry, 
                                          param=i,
                                          coeff_type=i_type,
                                          comp_str=component,
                                          start_end=start_end)
                plt.figure()
                plt.plot(azimuth, unp.nominal_values(combined), '.',azimuth, unp.nominal_values(reconstructed_series), '-')
                plt.plot(azimuth2, unp.nominal_values(i),'o')
                plt.title(f"series: {combined_series_name}")
                
                plt.figure()
                plt.plot(azimuth, unp.std_devs(combined), '.',azimuth, unp.std_devs(reconstructed_series), '-')
                plt.plot(azimuth2, unp.std_devs(i),'o')
                plt.title(f"errors in {combined_series_name}")
                # stop
    
    else:
        
        if azimuth is None:
            # we need to calculate the positions of the azimuths given the order. 
            #the azimuth is none logic is only so that we can calculate the combined series at any and all azimuths if desired.
            if i_bc == "periodic":
                azimuth = np.linspace(start_end[0], start_end[1], (2*i_order+1) + 1)
                azimuth = azimuth[:-1]
            else:
                azimuth = np.linspace(start_end[0], start_end[1], (2*i_order+1))
                
        combined = combine_series(param_dict, azimuth=azimuth, start_end=start_end, **kwargs)
        # get combined series from the fits above; make fit-like dictionary for it.
        combined_series = {}
        combined_series[combined_series_name] = unp.nominal_values(combined)
        combined_series[combined_series_name+"_err"] = unp.std_devs(combined)
        combined_series[combined_series_name+"_type"] = coefficient_type_as_string(i_type)
    
        if logger.is_below_level(level="DEBUG") and kwargs.get("prevent_interative_looping", False) is False:
            """
            test the new series can reproduce the values and errors in the original data.
            """
            # calculate all the combined series at all azimuths
            num_azimuths = kwargs.get("num_azimuths", 720)
            kwargs["prevent_interative_looping"] = True
            azims_all = np.linspace(start_end[0], start_end[1], num_azimuths)
            combined_all_azimuths = get_combined_series(param_dict,
                    azimuth = azims_all,
                    start_end=start_end,
                    **kwargs)
            # put values and errors together so can get errors from series. 
            i = unp.uarray(combined_series[combined_series_name], combined_series[combined_series_name+"_err"])
            reconstructed_series = coefficient_expand(azims_all, 
                                      param=i,
                                      coeff_type=i_type,
                                      comp_str=combined_series_name,
                                      start_end=start_end)
            
            plt.figure()
            plt.plot(azims_all, combined_all_azimuths[combined_series_name], 'r.', label="combined and expanded values")
            plt.plot(azimuth, unp.nominal_values(combined), 'ob', label="Series tie points")
            plt.plot( azims_all, unp.nominal_values(reconstructed_series), '-', label="reconstructed and expanded series")
            plt.legend()
            plt.title(f"series: {combined_series_name}; order={i_order}; i_type={coefficient_type_as_string(i_type)}")
            
            plt.figure()
            plt.plot(azims_all, combined_all_azimuths[combined_series_name+"_err"], 'r.', label="combined and expanded errors")
            plt.plot(azimuth, unp.std_devs(combined), 'ob', label="Series tie points errors")
            plt.plot(azims_all, unp.std_devs(reconstructed_series), '-', label="reconstructed and expanded errors")
            plt.legend()
            plt.title(f"errors in {combined_series_name}; order={i_order}; i_type={coefficient_type_as_string(i_type)}")
            # stop
    

    if 0 and logger.is_below_level(level="DEBUG"):
        """
        This loop does everything again but records the time taken.
        Only to be used for testing purposes. 
        """
        import time
        s1 = time.time()
        kwargs["prevent_interative_looping"] = True
        combined_series = get_combined_series(param_dict,
                            azimuth = None,
                            start_end=[0, 360],
                            **kwargs)
        if logger.is_below_level(level="DEBUG"):
            e1 = time.time()
            s2 = time.time()
            kwargs.update({"fit_series_new": True})
            combined_series_old = get_combined_series(param_dict,
                                azimuth = None,
                                start_end=[0, 360],
                                **kwargs)
            e2 = time.time()
            print(f"new serires method: {e1-s1}s")
            print(f"old serires method: {e2-s2}s")
            print(f"speed up: {((e2-s2)-(e1-s1))/(e2-s2)*100}%")
            print(f"times faster: {(e2-s2)/(e1-s1)}")
    
    return combined_series



def series_properties(
    coefficients,
    correlation_coeffs=None,
    subpattern=0,
    peak=0,
    param = "height",
    azm_spacing = 0.01, #merge these parameters into 1
    start_end = [0, 360], #merge these parameters into 1
    **kwargs,
):
    """
    Calcualte series properties from the coefficients.

    Parameters
    ----------
    coefficients : dict
        Coeffecient dictionary used in cpf.
    correlation_coeffs : dict, optional
        correlation coefficient dictionary as created by lmfit. The default is None.
    subpattern : list, int, optional
        Which subpattern in the coefficients to calulcate parameters for. The default is 0.
    peak : int, optional
        Which peak in the subpattern to calulcate parameters for. The default is 0.
    param: str, optional
        Peak profile parameter to calculate properties for
    azm_spacing : float or list or np.array, optional
        either:
            Precision to calulate the properties for (if required).
        or:
            list of azimuths to calculate the properties at
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    differential_coefficients : dict
        Dictionary of the calculated properties and their errors. The propertues calculated from the
        series coefficients are:
            "series_mean" -- mean d-spacing of diffraction ring, assuming 2d or 3d 'SampleGeometry'
            "series_max" -- maximum values of series
            "orientation max" -- maximum values of series
            "series_min" -- minimum values of series
            "orientation min" -- minimum values of series

    """

    # validate the inputs.
    if isinstance(coefficients, dict):
        coefficients = [coefficients]

    if not isinstance(coefficients, list):
        raise ValueError("The coefficients need to be a list of dictionaries.")

    # catch 'null' terms in fits
    coefficients = replace_value(coefficients, old=None, new=np.nan)

    if "symmetry" in coefficients[subpattern]["peak"][peak]:
        sym = coefficients[subpattern]["peak"][peak]["symmetry"]
    else:
        sym = 1
        
    # calulate maximum and minimum and their positions.
    # import time
    # start_diff = time.time()
    if coefficients[subpattern]["peak"][peak][param + "_type"] == "fourier":
        # get 
        properties_tmp = fourier_properties(unp.uarray(coefficients[subpattern]["peak"][peak][param], coefficients[subpattern]["peak"][peak][param+"_err"]),
                                        symmetry=sym, 
                                        start_end=start_end
                                        )
    else:
        # stop
        # calulate maximum and minimum and their positions.
        if len(ma.unique(azm_spacing).compressed()) != 1:
            # then need to use uniquie azimuths that were fed in
            orientations = ma.unique(azm_spacing).compressed()
        else:
            n = (start_end[0]- start_end[0]) / azm_spacing + 1
            orientations = np.linspace(start_end[0], start_end[1], int(n))
        # differential spline and get roots. 
        properties_tmp = spline_properties(unp.uarray(coefficients[subpattern]["peak"][peak][param], coefficients[subpattern]["peak"][peak][param+"_err"]),
                                       azimuth = orientations,
                                       symmetry=sym, 
                                       start_end=start_end,
                                       spline_type = coefficients[subpattern]["peak"][peak][param + "_type"],
                                       # bc_type="periodic", #*************need settings
                                       # kind=None,  #******************need settings
                                       )
      
    properties = {}
    properties["series mean"] = unp.nominal_values(properties_tmp["series mean"])
    properties["series mean err"] = unp.std_devs(properties_tmp["series mean"])
    properties["series max"] = unp.nominal_values(properties_tmp["series max"])
    properties["series max err"] = unp.std_devs(properties_tmp["series max"])
    properties["series min"] = unp.nominal_values(properties_tmp["series min"])
    properties["series min err"] = unp.std_devs(properties_tmp["series min"])
    properties["series max orientation"] = unp.nominal_values(properties_tmp["series max orientation"])
    properties["series max orientation err"] = unp.std_devs(properties_tmp["series max orientation"])
    properties["series min orientation"] = unp.nominal_values(properties_tmp["series min orientation"])
    properties["series min orientation err"] = unp.std_devs(properties_tmp["series min orientation"])
    if param is not None:
        entries = list(properties)
        for i in range(len(entries)):
            properties[re.sub("series", param, entries[i])] = properties.pop(entries[i])
    return properties



def fourier_properties(coeffs, symmetry=1, start_end=[0,360]):
    """
    Calculate the mean, global maximum and global minimum of a fourier series.

    A spline is constructed from the supplied tie-point coefficients using
    ``construct_spline()``. The spline derivative is then formed analytically
    and converted to a piecewise polynomial representation. Stationary points
    are obtained from the exact roots of the derivative and are combined with
    the spline endpoints (and knots for linear splines) to determine the
    global extrema.

    If the input coefficients are ``uncertainties.ufloat`` objects, the
    uncertainties are propagated through the spline basis functions to obtain
    uncertainties on the mean, maximum and minimum values. Uncertainties on
    the extremum orientations are estimated using first-order implicit
    differentiation of the stationary-point condition dS/dx = 0.

    Parameters
    ----------
    coeffs : array-like
        Spline tie-point coefficients. Elements may be floats or ufloats.
    azimuth : array-like or None, optional
        Azimuth coordinates corresponding to the tie points (only used for independent series type).
    symmetry : int, optional
        symmetry of azimuthal values
    start_end : list-like, optional
        Lower and upper bounds of the spline domain.
    spline_type : str, optional
        name of spline type 

    Returns
    -------
    dict
        Dictionary containing mean, max, min and orientations.
    """
    """
    edited after copilot output
    """    
    
    def fourier_prime_nominal(azm, coeffs):
        """
        Nominal derivative f'(azm).
        Used for checking roots.
        """
        N = get_order_from_coeff(len(coeffs))
        coeffs2 = coeffs*0
        for n in range(1, N + 1):
            # reverse order of the values in the series so that sin/cos are switched in the fourier expansion
            coeffs2[2*n] = coeffs[2*n - 1] * n
            coeffs2[2*n-1] = coeffs[2*n] * n * -1
        y = fourier_expand(np.rad2deg(azm), inp_param=coeffs2, comp_str=None, no_negatives=False)
        return y
    
    def fourier_second_nominal(azm, coeffs):
        """
        Nominal second derivative f''(azm).
        Used for classification and uncertainty propagation.
        """
        N = get_order_from_coeff(len(coeffs))
        coeffs2 = coeffs*0
        for n in range(1, N + 1):
            # order of sin and cos reversed twice by double differential. So the order is back in the initial order
            coeffs2[2*n - 1] = coeffs[2*n - 1] * n**2 * -1
            coeffs2[2*n] = coeffs[2*n] * n**2 * -1
        y = fourier_expand(np.rad2deg(azm), inp_param=coeffs2, comp_str=None, no_negatives=False)
        return y
    
    def derivative_polynomial_roots(coeffs, unit_tol=1e-7, residual_tol=1e-7):
        """
        Solve f'(azm)=0 algebraically by writing it as a polynomial in z = exp(i azm).
    
        Returns nominal azm roots in [0, 2*pi).
        """
        N = get_order_from_coeff(len(coeffs))
        if N == 0:
            return np.array([])
        # Polynomial P(z) = z^N f'(azm), where z = exp(i azm)
        # Powers run from 0 to 2N
        poly = np.zeros(2*N + 1, dtype=complex)
        for n in range(1, N + 1):
            a_n = unp.nominal_values(coeffs[2*n - 1])
            b_n = unp.nominal_values(coeffs[2*n])
            # f'(azm) contains:
            # n/2 * (a_n + i b_n) z^n
            # n/2 * (a_n - i b_n) z^(-n)
            c_pos = 0.5 * n * (a_n + 1j*b_n)
            c_neg = 0.5 * n * (a_n - 1j*b_n)
            power_pos = N + n
            power_neg = N - n
            # np.roots wants descending powers
            poly[2*N - power_pos] += c_pos
            poly[2*N - power_neg] += c_neg
        poly = np.trim_zeros(poly, trim="f")
        if len(poly) <= 1:
            return np.array([])
        z_roots = np.roots(poly)
        # Keep only roots close to the unit circle
        z_roots = z_roots[np.abs(np.abs(z_roots) - 1.0) < unit_tol]
        azm = np.mod(np.angle(z_roots), 2*np.pi)
        # Check derivative residual
        azm = np.array([x for x in azm if abs(fourier_prime_nominal(x, coeffs)) < residual_tol])
        # Remove duplicates
        azm = np.unique(np.round(azm, 12))
        return azm
    
    
    def uncertain_root_azm(azm0, coeffs):
        """
        First-order uncertainty propagation for the stationary point location.
    
        If g(azm, coeffs) = f'(azm) = 0, then
            dazm/dc = - (dg/dc) / (dg/dazm)
        where dg/dazm = f''(azm).
        """
        N = get_order_from_coeff(len(coeffs))
        f2 = fourier_second_nominal(azm0, coeffs)
        if abs(f2) < 1e-14:
            # Degenerate stationary point, uncertainty in azm is ill-defined
            return unp.ufloat(azm0, np.nan)
        azm_u = azm0
        for n in range(1, N + 1):
            a_n = coeffs[2*n - 1]
            b_n = coeffs[2*n]
            a0 = unp.nominal_values(a_n)
            b0 = unp.nominal_values(b_n)
            # g = f' = sum n a_n cos(n azm) - n b_n sin(n azm)
            dg_da = n * np.cos(n * azm0)
            dg_db = -n * np.sin(n * azm0)
            dazm_da = -dg_da / f2
            dazm_db = -dg_db / f2
            azm_u += dazm_da * (a_n - a0)
            azm_u += dazm_db * (b_n - b0)
        return azm_u
    
    mean = coeffs[0]
    azm_roots = derivative_polynomial_roots(coeffs)
    if len(azm_roots) == 0:
        return {
            "series mean": mean,
            "series max": None,
            "series max orientation": np.nan,
            "series min": None,
            "series min orientation": np.nan,
        }
    stationary = []
    for azm0 in azm_roots:
        azm_u = uncertain_root_azm(azm0, coeffs) # in radians
        # val_u = fourier_value(azm_u, coeffs)
        
        val_u = fourier_expand(np.rad2deg(unp.nominal_values(azm_u)), inp_param=coeffs, comp_str=None, no_negatives=False)
        # stop
        f2 = fourier_second_nominal(azm0, coeffs)
        if f2 < 0:
            kind = "local maximum"
        elif f2 > 0:
            kind = "local minimum"
        else:
            kind = "flat / degenerate"
        stationary.append({
            "azm": unp.degrees(azm_u),
            "value": val_u,
            "second_derivative": f2,
            "kind": kind
        })
    # Global max/min by nominal value
    values_nom = np.array([unp.nominal_values(p["value"]) for p in stationary])
    i_max = np.argmax(values_nom)
    i_min = np.argmin(values_nom)
    result = {
        "series mean": mean,
        "series max": stationary[i_max]["value"],
        "series max orientation": stationary[i_max]["azm"]/symmetry,
        "series min": stationary[i_min]["value"],
        "series min orientation": stationary[i_min]["azm"]/symmetry,
    }
    if 0:    
        print("Mean:")
        print(result["series mean"])
        print("\nMaximum:")
        print("azm =", result["series max orientation"])
        print("f   =", result["series max"])
        print("\nMinimum:")
        print("azm =", result["series min orientation"])
        print("f   =", result["series min"])
        if 0:
            print("\nAll stationary points:")
            for p in result["stationary_points"]:
                print(p["kind"])
                print("  azm =", p["azm"])
                print("  f   =", p["value"])
    return result
    



def spline_properties(coeffs,
                      azimuth=None,
                      symmetry=1,
                      start_end=[0, 360],
                      spline_type='spline_cubic'):
    """
    Calculate the mean, global maximum and global minimum of a spline series.

    A spline is constructed from the supplied tie-point coefficients using
    ``construct_spline()``. The spline derivative is then formed analytically
    and converted to a piecewise polynomial representation. Stationary points
    are obtained from the exact roots of the derivative and are combined with
    the spline endpoints (and knots for linear splines) to determine the
    global extrema.

    If the input coefficients are ``uncertainties.ufloat`` objects, the
    uncertainties are propagated through the spline basis functions to obtain
    uncertainties on the mean, maximum and minimum values. Uncertainties on
    the extremum orientations are estimated using first-order implicit
    differentiation of the stationary-point condition dS/dx = 0.

    Parameters
    ----------
    coeffs : array-like
        Spline tie-point coefficients. Elements may be floats or ufloats.
    azimuth : array-like or None, optional
        Azimuth coordinates corresponding to the tie points (only used for independent series type).
    symmetry : int, optional
        symmetry of azimuthal values
    start_end : list-like, optional
        Lower and upper bounds of the spline domain.
    spline_type : str, optional
        name of spline type 

    Returns
    -------
    dict
        Dictionary containing mean, max, min and orientations.
    """
    """
    edited after copilot output
    """

    # get spline type
    series_name = coefficient_type_as_string(spline_type)
    all_series = coefficient_types(full=True)

    coeffs = np.asarray(coeffs, dtype=object)

    # Nominal spline
    coeffs_nom = np.array([unp.nominal_values(c) for c in coeffs],dtype=float)
    spl, tie_points = construct_spline(
        coeffs=coeffs_nom,
        azimuth=azimuth,
        start_end=start_end,
        bc_type=all_series[series_name]["boundary_conditions"],
        kind=all_series[series_name]["spline_type"],
    )
    xmin = spl.t[spl.k]
    xmax = spl.t[-spl.k - 1]

    # Candidate extrema locations
    if spl.k == 1:
        # Piecewise-linear spline:
        # extrema occur at knots, not derivative roots
        x_candidates = np.unique(spl.t[spl.k:-spl.k])
    else:
        try:
            dspl = spl.derivative()
            pp = PPoly.from_spline(dspl)
            roots = pp.roots(extrapolate=False)
            roots = roots[ (roots >= xmin) & (roots <= xmax) ]
        except Exception:
            roots = np.array([])
        x_candidates = np.concatenate(([xmin], roots, [xmax]))
        x_candidates = np.unique(np.round(x_candidates, 12))

    # Nominal values
    y_nom = spl(x_candidates)
    imax = np.argmax(y_nom)
    imin = np.argmin(y_nom)
    xmax_nom = x_candidates[imax]
    xmin_nom = x_candidates[imin]

    # Build basis splines
    n = len(coeffs)
    eye = np.eye(n)
    basis = []
    for i in range(n):
        b, _ = construct_spline(
            coeffs=eye[i],
            azimuth=azimuth,
            start_end=start_end,
            bc_type=all_series[series_name]["boundary_conditions"],
            kind=all_series[series_name]["spline_type"],
        )
        basis.append(b)

    # Evaluate spline with uncertainty propagation
    def eval_ufloat(x):
        out = 0
        for c, b in zip(coeffs, basis):
            out += c * b(x)
        return out

    # Exact mean
    L = xmax - xmin
    mean = 0
    for c, b in zip(coeffs, basis):
        mean += c * ( b.integrate(xmin, xmax) / L )
    # Max/min values
    series_max = eval_ufloat(xmax_nom)
    series_min = eval_ufloat(xmin_nom)
    # Extremum orientation uncertainty
    def orientation_error(x0):
        # Linear splines have discontinuous derivatives
        if spl.k == 1:
            return ufloat(x0, 0.0)
        try:
            f2 = spl.derivative(2)(x0)
            if abs(f2) < 1e-12:
                return ufloat(x0, np.nan)
            xu = ufloat(x0, 0.0)
            for c, b in zip(coeffs, basis):
                dg_dc = b.derivative()(x0)
                xu += -(dg_dc / f2) * (c - unp.nominal_values(c))
            return xu
        except Exception:
            return ufloat(x0, np.nan)
    max_orientation = orientation_error(xmax_nom)
    min_orientation = orientation_error(xmin_nom)
    return {
        "series mean": mean,
        "series max": series_max,
        "series max orientation": max_orientation,
        "series min": series_min,
        "series min orientation": min_orientation,
    }