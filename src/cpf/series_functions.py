#!/usr/bin/env python3

"""
Functions describing the fourier or spline series used as parameter sets for the
peak shape and background properties.
"""

__all__ = [
    "coefficient_types",
    "coefficient_type_as_number",
    "coefficient_type_as_string",
    "get_series_type",
    "params_get_type",
    "get_number_coeff",
    "get_order_from_coeff",
    "get_order_from_params",  # not needed?
    "fourier_order",  # same as get order from params?
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
from cpf.logger_functions import logger


def coefficient_types():
    """
    Defines the coefficient types possible and the numbers associated with them.
    :return: dictionary of key, value pairs.
    """
    coeff_types = {
        "fourier": 0,
        "spline_linear": 1,
        "linear": 1,
        "spline_quadratic": 2,
        "quadratic": 2,
        "spline_cubic_closed": 3,
        "spline_cubic": 3,
        "spline-cubic": 3,
        "cubic": 3,
        "spline": 3,
        "spline_cubic_open": 4,
        "independent": 5,
    }
    return coeff_types


def coefficient_type_as_number(coeff_type, return_error=1):
    """

    :param coeff_type:
    :return: numerical index for series type
    """
    # FIX ME: DMF can refactor function as dictionary key,value pairs and call in print
    if coeff_type == "fourier" or coeff_type == 0:
        out = 0
    elif coeff_type == "spline_linear" or coeff_type == "linear" or coeff_type == 1:
        out = 1
    elif (
        coeff_type == "spline_quadratic" or coeff_type == "quadratic" or coeff_type == 2
    ):
        out = 2
    elif (
        coeff_type == "spline_cubic_closed"
        or coeff_type == "spline_cubic"
        or coeff_type == "spline-cubic"
        or coeff_type == "cubic"
        or coeff_type == "spline"
        or coeff_type == 3
    ):
        out = 3
    elif coeff_type == "spline_cubic_open" or coeff_type == 4:
        out = 4
    elif coeff_type == "independent" or coeff_type == 5:
        out = 5
    else:
        error_str = "Unrecognised coefficient series type, the valid options are fourier, etc..."
        if return_error == 0:
            raise ValueError(error_str)
            # FIX ME: write out all the licit options in the error message.
        else:
            out = error_str
    return out


def coefficient_type_as_string(series_type):
    """

    :param series_type: numeric series type.
    :return: strong index for series type
    """
    # FIX ME: DMF can refactor as dictionary
    if series_type == 0:  # or series_type=='fourier':
        out = "fourier"
    elif (
        series_type == 1
    ):  # or series_type == 'spline_linear' or series_type == 'linear':
        out = "spline_linear"
    elif (
        series_type == 2
    ):  # or series_type == 'spline_quadratic' or series_type == 'quadratic':
        out = "spline_quadratic"
    elif (
        series_type == 3
    ):  # or series_type == 'spline_cubic_closed' or series_type == 'spline_cubic'
        # or series_type == 'spline-cubic' or series_type == 'cubic'
        # or series_type == 'spline':
        out = "spline_cubic"
    elif (
        series_type == 4
    ):  # series_type == 'spline_cubic_open' or series_type == 'spline-cubic-open'
        # or series_type == 'cubic-open' or series_type == 'spline-open':
        out = "spline_cubic_open"
    elif series_type == 5:  # or series_type == 'independent':
        out = "independent"
    else:
        raise ValueError(
            "Unrecognised coefficient series type, the valid options are "
            "fourier"
            ", etc..."
        )
        # FIX ME: write out all the licit options in the error message.
    return out


def get_series_type(param, param_str, comp=None):
    """
    Make a nested list of parameters and errors from a lmfit Parameters class
    :param param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: nested list of [parameters, errors] in alphanumerical order
    """
    if comp is not None:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    new_str = new_str + "_tp"
    if isinstance(param, dict) and new_str in param:
        out = param[new_str]
    elif isinstance(param, dict) and new_str not in param:
        out = 0
    else:
        out = None
    return out


def get_series_mean(param, param_str, comp=None):
    """
    Calcualte the mean of the parameter series a lmfit Parameters class
    :param param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: weighted mean of parameters
    """

    # if comp is not None:
    #     new_str = param_str + "_" + comp
    # else:
    #     new_str = param_str

    # FIX ME: here we need to be able to discard the outliers.
    # We should use the medaian and the mean deviation from the median...

    if get_series_type(param, param_str, comp=comp) == 0:
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


def params_get_type(orders, comp, peak=0):
    """
    Parameters
    ----------
    orders : TYPE
        DESCRIPTION.
    comp : TYPE
        DESCRIPTION.
    peak : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    """
    comp_str = pf.expand_component_string(comp) + "_type"
    if comp_str in orders:
        coeff_type = orders[comp_str]
    elif len(orders["peak"]) > peak and comp_str in orders["peak"][peak]:
        coeff_type = orders["peak"][peak][comp_str]
    else:
        coeff_type = "fourier"

    return coeff_type


def get_number_coeff(orders, comp, peak=0, azimuths=None):
    """
    Gets the expected number of coefficients from orders.
    If independent this needs the azimuths as well.
    :param orders:
    :param comp:
    :param peak:
    :param azimuths:
    :return:
    """

    parm_str = params_get_type(orders, comp, peak)
    parm_num = coefficient_type_as_number(parm_str)

    if parm_num == 5:  # independent
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
    if parm_num == 5:  # independent
        if azimuths is None:
            raise ValueError(
                "Cannot define number of independent values without a number of coefficients."
            )
        else:
            order = azimuths.shape[0]
    else:  # everything else.
        order = (n_coeff - 1) / 2
    return order


# Not needed?
def get_order_from_params(params, comp=None, peak=0):
    # Given list of Fourier coefficients return order (n)
    if isinstance(params, (list,)):
        l = len(params)
    elif isinstance(params, (float,)):
        l = np.size(params)
    elif isinstance(params, (dict,)):
        parm_str = params_get_type(params, comp, peak)
        if comp == "bg" or comp == "background" or comp == "f":
            l = np.max(params["background"][peak])
        else:  # everything else.
            l = np.max(params["peak"][peak][pf.expand_component_string(comp)])
    elif isinstance(params, (int,)):
        l = 1
    else:
        logger.debug(" ".join(map(str, [("params", params)])))
        logger.debug(" ".join(map(str, [("type(params)", type(params))])))
        err_str = "Parameter list is not list, float or a dictionary."
        logger.critical(" ".join(map(str, [(err_str)])))
        raise ValueError(err_str)

    order = get_order_from_coeff(l)

    return order


def fourier_order(params):
    """

    :param params:
    :return:
    """
    # Given list of Fourier coefficients return order (n) of the Fourier series.

    if isinstance(params, (list,)):
        order = int((len(params) - 1) / 2)
    elif isinstance(params, (float,)):
        order = int((np.size(params) - 1) / 2)
    elif isinstance(params, (int,)):
        order = 0
    else:
        logger.debug(" ".join(map(str, [("params", params)])))
        logger.debug(" ".join(map(str, [("type(params)", type(params))])))
        err_str = "Parameter list is not list, float or an integer. I do not know what to do with it."
        logger.critical(" ".join(map(str, [(err_str)])))
        raise ValueError(err_str)
    return order


def coefficient_expand(
    azimuth,
    param=None,
    coeff_type="fourier",
    comp_str=None,
    start_end=[0, 360],
    **params,
):
    """
    Expand the coefficients to give one value for each azimuth.
    :param azimuth:
    :param param:
    :param coeff_type:
    :param comp_str:
    :param params:
    :return: coefficient value at each azimuth
    """

    coeff_type = coefficient_type_as_number(coeff_type)

    if coeff_type == 0:
        out = fourier_expand(azimuth, inp_param=param, comp_str=comp_str, **params)
    elif coeff_type == 1:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            bc_type="natural",
            kind="linear",
            **params,
        )
    elif coeff_type == 2:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            bc_type="natural",
            kind="quadratic",
            **params,
        )
    elif coeff_type == 3:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            bc_type="periodic",
            kind="cubic",
            **params,
        )
    elif coeff_type == 4:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            bc_type="natural",
            kind="cubic",
            **params,
        )
    elif coeff_type == 5:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=start_end,
            bc_type="natural",
            kind="independent",
            **params,
        )
    else:
        raise ValueError(
            "Unrecognised coefficient series type, the valid options are "
            "fourier"
            ", etc..."
        )
        # FIX ME: write out all the licit options in the error message.
    return out


# spline expansion function
def spline_expand(
    azimuth,
    inp_param=None,
    comp_str=None,
    start_end=[0, 360],
    bc_type="periodic",
    kind="cubic",
    **params,
):
    """Calculate Spline interpolation given input coefficients.
    :param kind:
    :param start_end:
    :param bc_type:
    :param azimuth: arr data array float
    :param inp_param: list of values at spline tie points
    :param comp_str: str to determine which coefficients to use from params
    :param params: lmfit dict of coefficients as parameters
    :return:
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

    if kind == "cubic":  # and bc_type=='periodic':
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


# fourier expansion function
def fourier_expand(
    azimuth, inp_param=None, comp_str=None, start_end=[0, 360], **params
):
    """Calculate Fourier expansion given input coefficients
    :param azimuth: arr data array float
    :param inp_param: list of Fourier coefficients float
    :param comp_str: str to determine which coefficients to use from params
    :param params: lmfit dict of coefficients as parameters
    :return: coefficient value for each azimuth
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
            )  # single col array
            # except:
            #    pass
    return np.squeeze(fout)


# fourier expansion function
def background_expansion(azimuth_two_theta, orders, params):
    """
    Calculate the Fourier expansion of the background terms
    :param coeff_type:
    :param azimuth_two_theta: list of arrays for azimuth and theta float
    :param inp_param: list of input parameters
    :return: Fourier expansion result float arr
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
                b_tp = 0

        backg.append(np.array(fourier_background))
        backg_tp.append(b_tp)
        i = i + 1

    # Not sure if this is correct, thought that if orders then array would be eg [5,0] and would want a fourier
    # expansion with 5 parms for offset and then no fourier expansion for slope i.e. 5 parms then 0.
    # But below would interpret this as effectively [0,0] instead.
    # ANSWER: I am not sure this is true any more.

    bg_all = np.zeros(azimuth.shape)
    for i in range(len(backg)):
        out = coefficient_expand(azimuth, backg[i], backg_tp[i])
        bg_all = bg_all + (out * (two_theta_prime ** float(i)))
    return bg_all
