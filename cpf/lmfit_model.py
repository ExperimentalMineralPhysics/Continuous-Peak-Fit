#!/usr/bin/env python

"""
Functions using lmfit to make a semi-bounded fitting model. 
"""

__all__ = [
    "parse_bounds",
    # "create_new_params",
    "params_to_new_params",
    # "param_fit_array",
    "gather_param_errs_to_list",
    "gather_params_to_list",
    "gather_params_from_dict",
    "initiate_all_params_for_fit",
    "initiate_params",
    "vary_params",
    "un_vary_params",
    "un_vary_part_params",
    "un_vary_single_param",
    "peaks_model",
    "fit_model",
    "coefficient_fit",
]

import numpy as np
import sys
import warnings
from lmfit import Parameters, Model
import cpf.peak_functions as pf
import cpf.series_functions as sf

np.set_printoptions(threshold=sys.maxsize)


def parse_bounds(bounds, data_as_class, ndat=None, n_peaks=1, param=None):
    """
    Turn text limits in bounds into numbers and return bounds in dictionary format.
    :param bounds:
    :param dspace:
    :param intens:
    :param two_theta:
    :param ndat:
    :param n_peaks:
    :param param:
    :return:
    """
    limits = {}
    if bounds is None:
        bounds = {
            "background": ["min", "max"],
            "d-space": ["min", "max"],
            "height": [0, "max-min"],
            "profile": [0, 1],
            "width": ["range/(ndata)", "range/2/npeaks"],
        }
    if ndat is None:
        ndat = np.size(data_as_class.dspace)

    choice_list = ["d-space", "height", "width", "profile", "background"]
    if param is not None:
        choice_list = param

    for par in choice_list:
        if par == "height" or par == "background":
            vals = data_as_class.intensity
        elif par == "width":
            vals = data_as_class.tth
        else:
            vals = data_as_class.dspace

        b = bounds[par]
        b = [str(w).replace("inf", "np.inf") for w in b]
        b = [str(w).replace("range", "(max-min)") for w in b]
        b = [w.replace("ndata", str(ndat)) for w in b]
        b = [w.replace("max", str(np.max(vals))) for w in b]
        b = [w.replace("min", str(np.min(vals))) for w in b]
        b = [w.replace("npeaks", str(n_peaks)) for w in b]
        b = [eval(w) for w in b]
        limits[par] = b

    return limits


# def create_new_params(num_peaks, dfour, hfour, wfour, pfour, bg_four, order_peak):
#     """
#     Create the new_params dictionary from results
#     :param num_peaks: number of peaks int
#     :param dfour: list of coefficients per peak float
#     :param hfour: list of coefficients per peak float
#     :param wfour: list of coefficients per peak float
#     :param pfour: list of coefficients per peak float
#     :param bg_four:
#     :param order_peak:
#     :return: new_params dictionary
#     """
#     peaks = []
#     for j in range(num_peaks):
#         if (
#                 "symmetry" in order_peak[j]
#         ):  # orders['peak'][j]:  # FIX ME: the symmetry part of this is a horrible way of doing it.
#             peaks.append(
#                 {
#                     "d-space": dfour[j][0],
#                     "height": hfour[j][0],
#                     "width": wfour[j][0],
#                     "profile": pfour[j][0],
#                     "symmetry": order_peak[j]["symmetry"],
#                 }
#             )
#             # "symmetry": orders['peak'][j]['symmetry']})
#         else:
#             peaks.append(
#                 {
#                     "d-space": dfour[j][0],
#                     "height": hfour[j][0],
#                     "width": wfour[j][0],
#                     "profile": pfour[j][0],
#                 }
#             )
#     new_params = {"background": bg_four, "peak": peaks}

#     return new_params


def params_to_new_params(params, orders=None):
    """
    Gather Parameters class to new_params dictionary format
    :param params:
    :param num_peaks:
    :param num_bg:
    :param order_peak:
    :return:
    """

    num_peaks = len(orders["peak"])
    num_bg = len(orders["background"])

    peaks = []
    for i in range(num_peaks):
        new_str = "peak_" + str(i) + "_d"
        d_space = gather_param_errs_to_list(params, new_str)
        d_space_tp = sf.get_series_type(params, new_str)
        new_str = "peak_" + str(i) + "_h"
        h_space = gather_param_errs_to_list(params, new_str)
        h_tp = sf.get_series_type(params, new_str)
        new_str = "peak_" + str(i) + "_w"
        w_space = gather_param_errs_to_list(params, new_str)
        w_tp = sf.get_series_type(params, new_str)
        new_str = "peak_" + str(i) + "_p"
        p_space = gather_param_errs_to_list(params, new_str)
        p_tp = sf.get_series_type(params, new_str)

        tmp_peaks = {}
        if "phase" in orders["peak"][i]:
            tmp_peaks["phase"] = orders["peak"][i]["phase"]
        if "hkl" in orders["peak"][i]:
            tmp_peaks["hkl"] = orders["peak"][i]["hkl"]
        tmp_peaks["d-space"] = d_space[0]
        tmp_peaks["d-space_err"] = d_space[1]
        tmp_peaks["d-space_type"] = sf.coefficient_type_as_string(d_space_tp)
        tmp_peaks["height"] = h_space[0]
        tmp_peaks["height_err"] = h_space[1]
        tmp_peaks["height_type"] = sf.coefficient_type_as_string(h_tp)
        tmp_peaks["width"] = w_space[0]
        tmp_peaks["width_err"] = w_space[1]
        tmp_peaks["width_type"] = sf.coefficient_type_as_string(w_tp)
        tmp_peaks["profile"] = p_space[0]
        tmp_peaks["profile_err"] = p_space[1]
        tmp_peaks["profile_type"] = sf.coefficient_type_as_string(p_tp)
        if "symmetry" in orders["peak"][i]:
            tmp_peaks["symmetry"] = orders["peak"][i]["symmetry"]

        peaks.append(tmp_peaks)

    # Get background parameters
    bg_space = []
    bg_space_err = []
    for b in range(num_bg):
        new_str = "bg_c" + str(b) + "_f"
        bg_spc = gather_param_errs_to_list(params, new_str)
        bg_space.append(bg_spc[0])
        bg_space_err.append(bg_spc[1])
    bg_tp = sf.coefficient_type_as_string(sf.get_series_type(params, new_str))

    new_params = {
        "background": bg_space,
        "background_err": bg_space_err,
        "background_type": bg_tp,
        "peak": peaks,
    }

    return new_params


# def flatten(li):
#     """
#     :param li:
#     :return:
#     """
#     return sum(([x] if not isinstance(x, list) else flatten(x) for x in li), [])


# def param_fit_array(num_p, len_bg, val=None, p=None, b=None, orda=None):
#     """
#     Makes array which represents all the fourier series in the fit.
#     new_arr{0} is a list of 4 element lists - one 4 element list for each peak.
#     new_arr[0] is the background representative.
#     :param num_p:
#     :param len_bg:
#     :param val:
#     :param p:
#     :param b:
#     :param orda:
#     :return:
#     """
#     if val:
#         v = val
#     else:
#         v = 0

#     if not p:
#         p = [v, v, v, v]

#     new_arr = [[], []]
#     for y in range(num_p):
#         new_arr[0].append(p[:])
#     for y in range(len_bg):
#         new_arr[1].extend([v][:])
#     return new_arr


def gather_param_errs_to_list(
    inp_param, param_str, comp=None, nterms=None
):  # fix me: Doesn't run without Comp=None.
    # Maybe it was never validated before SAH ran it.
    """
    Make a nested list of parameters and errors from a lmfit Parameters class
    :param inp_param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: nested list of [parameters, errors] in alphanumerical order
    """

    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    param_list = []
    err_list = []
    total_params = []
    str_keys = [
        key for key, val in inp_param.items() if new_str in key and "tp" not in key
    ]
    # print('yes', new_str, str_keys)
    for i in range(len(str_keys)):
        param_list.append(inp_param[new_str + str(i)].value)
        err_list.append(inp_param[new_str + str(i)].stderr)
    total_params.append(param_list)
    total_params.append(err_list)

    return total_params


def gather_params_to_list(inp_param, param_str, comp=None):
    """
    Make a list of parameters from a lmfit Parameters class
    :param inp_param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: list of parameters in alphanumerical order
    """

    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    param_list = []
    str_keys = [key for key, val in inp_param.items() if new_str in key]
    # print('yes', new_str, str_keys)
    for i in range(len(str_keys)):
        param_list.append(inp_param[new_str + str(i)].value)
    return param_list


def gather_params_from_dict(inp_param, param_str, comp):
    """
    Make a list of parameters from a dictionary
    :param inp_param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: list of parameters in alphanumerical order
    """

    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    param_list = []
    str_keys = [
        key for key, val in inp_param.items() if new_str in key and "tp" not in key
    ]
    for i in range(len(str_keys)):
        param_list.append(inp_param[new_str + str(i)])
    return param_list


def initiate_all_params_for_fit(
    settings_as_class, data_as_class, values=None, debug=False
):
    """
    Initiate all the parameters needed for the fitting with lmfit.

    :param inp_param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: list of parameters in alphanumerical order
    """

    comp_list, comp_names = pf.peak_components(full=False)

    master_params = Parameters()  # Initiate total Parameter class to add to

    lims = parse_bounds(
        settings_as_class.fit_bounds, data_as_class
    )  # d_space.flatten(), intens.flatten(), twotheta.flatten()
    # )

    # for i in range(len(comp_list)):

    #     comp_list_section = comp_list[i]
    #     comp_names_section = comp_names[i]

    #     for j in range(len(comp_list_section)):

    #         comps = comp_list_section[j]
    #         comp_name = comp_names_section[j]

    #         if comps[0] == "bg":

    # Initiate background parameters
    comps = "bg"
    coeff_type = sf.params_get_type(settings_as_class.subfit_orders, comps)
    for k in range(len(settings_as_class.subfit_orders["background"])):
        param_str = "bg_c" + str(k)
        comp = "f"
        if k == 0:
            limits = lims["background"]
        else:
            limits = None
        n_coeff = sf.get_number_coeff(
            settings_as_class.subfit_orders, comps, peak=k, azimuths=data_as_class.azm
        )
        # add values is present
        vals = None
        if values:
            vals = (values["background"][k],)
        # intiate parameters
        master_params = initiate_params(
            master_params,
            param_str,
            comp,
            coeff_type=coeff_type,
            num_coeff=n_coeff,
            value=vals,
            trig_orders=settings_as_class.subfit_orders["background"][k],
            limits=limits,
        )

    # initiate peak(s)
    peak_comp_list = comp_list
    comp_names_list = comp_names

    for j in range(len(settings_as_class.subfit_orders["peak"])):
        # defines peak string to start parameter name
        param_str = "peak_" + str(j)

        # initiate symmetry
        comp = "s"
        if "symmetry" in settings_as_class.subfit_orders["peak"][j].keys():
            symm = settings_as_class.subfit_orders["peak"][j]["symmetry"]
        else:
            symm = 1
        master_params = initiate_params(
            master_params,
            param_str,
            comp,
            0,
            limits=[10, 1],
            value=np.array(symm),
            vary=False,
        )

        comp_list = ["h", "d", "w", "p"]
        comp_names = ["height", "d-space", "width", "profile"]

        for cp in range(len(comp_list)):
            comp = comp_list[cp]
            if comp == "d":
                symmetry = 1
            else:
                symmetry = symm

            # fixed = 0
            vals = None
            if comp_names[cp] + "_fixed" in settings_as_class.subfit_orders["peak"][j]:
                # fixed = 1
                vals = settings_as_class.subfit_orders["peak"][j][
                    comp_names[cp] + "_fixed"
                ]  # values if fixed.
            elif values:
                vals = values["peak"][j][comp_names[cp]]

            coeff_type = sf.params_get_type(
                settings_as_class.subfit_orders, comp, peak=j
            )
            n_coeff = sf.get_number_coeff(
                settings_as_class.subfit_orders,
                comp,
                peak=j,
                azimuths=data_as_class.azm,
            )

            master_params = initiate_params(
                master_params,
                param_str,
                comp,
                coeff_type=coeff_type,
                num_coeff=n_coeff,
                trig_orders=settings_as_class.subfit_orders["peak"][j][comp_names[cp]],
                limits=lims[comp_names[cp]],
                value=vals,
            )

            # FIX ME. Check which params should be varying and which should not.
            # Need to incorporate vary and un-vary params as well as partial vary
    return master_params


def initiate_params(
    inp_param,
    param_str,
    comp,
    coeff_type=None,
    num_coeff=None,
    trig_orders=None,
    limits=None,
    expr=None,
    value=None,
    ind_vars=None,
    vary=True,
):
    """
    Create all required coefficients for no. terms
    :param vary:
    :param num_coeff:
    :param coeff_type:
    :param inp_param: lmfit Parameter class (can be empty)
    :param param_str: base string to create lmfit parameter
    :param comp: string to add to base string to create lmfit parameter
    :param trig_orders: number of trigonometric orders
    :param limits: list [max,min] or None
    :param expr: str or None, str expression to use or 'calc'
    :param value: int if don't know nterms, use this many orders
    :param ind_vars: array-based independent variable to be added to parameters to use in expr
    :return:
    """
    if limits:
        new_min = np.min(limits)
        new_max = np.max(limits)
        half_range = (new_max - new_min) / 2
    else:
        new_min = -np.inf
        new_max = np.inf
        half_range = np.inf
    if value is None and limits is None:
        value = 0.01
    elif value is None and limits:
        value = new_min + (new_max - new_min) * 3 / 5  # /2#*2/3
        # N.B. it seems the initiating the splines with values exactly half-way between max and min doesn't work.
        # Anywhere else in the middle of the range seems to do so.
        # FIX ME: Is this fitting-related or a code bug?
    value = np.array(value)  # force input to be array so that can be iterated over

    coeff_type = sf.coefficient_type_as_number(coeff_type)
    if (
        coeff_type == 5
    ):  # if the values at each azimuth are independent we have to have a number of coefficients given.
        if num_coeff is not None:
            num_coeff = num_coeff
        elif value is not None:
            try:
                # FIX ME: DMF is this resolved?
                # trig_orders = int((len(value) - 1)/2) #FIX ME: coefficient series length
                # -- make into separate function.
                num_coeff = len(value)
            except TypeError or AttributeError:
                # trig_orders = int((value.size - 1)/2)
                # FIX ME: I don't know why the try except is needed, but it is. it arose during the spline development
                # DMF: len is for lists, size is for arrays - need to track if this difference should exist
                # The real confusion though is len should work for arrays as well...
                num_coeff = value.size
        else:
            print(value)
            raise ValueError(
                "Cannot define independent values without a number of coefficients."
            )

    elif trig_orders is None:
        try:
            # trig_orders = int((len(value) - 1)/2) #FIX ME: coefficient series length -- make into separate function.
            num_coeff = len(value)
        except TypeError or AttributeError:
            # FIX ME: I don't know why the try except is needed, but it is. it arose during the spline development
            # DMF: len is for lists, size is for arrays - need to track if this difference should exist
            # The real confusion though is len should work for arrays as well...
            num_coeff = value.size
    else:
        num_coeff = np.max(trig_orders) * 2 + 1
        # leave the spline with the 2n+1 coefficients so that it matches the equivalent fourier series.
        # FIX ME: should we change this to make it more natural to understand?
        # DMF: how?
    # First loop to add all the parameters so can be used in expressions
    for t in range(num_coeff):
        if value.size == num_coeff and value.size > 1:
            v = value.item(t)
        elif coeff_type != 0 or t == 0:
            v = value.item(0)
        else:
            v = 0
        if t == 0 or coeff_type != 0:
            inp_param.add(
                param_str + "_" + comp + str(t),
                v,
                max=new_max,
                min=new_min,
                expr=expr,
                vary=vary,
            )
        else:
            inp_param.add(
                param_str + "_" + comp + str(t),
                v,
                max=half_range,
                min=-half_range,
                expr=expr,
                vary=vary,
            )

    if comp != "s":
        inp_param.add(
            param_str + "_" + comp + "_tp",
            sf.coefficient_type_as_number(coeff_type),
            expr=expr,
            vary=False,
        )

    return inp_param


def vary_params(inp_param, param_str, comp):
    """
    Set selected parameters to vary=True
    :param inp_param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    str_keys = [
        key for key, val in inp_param.items() if new_str in key and "tp" not in key
    ]
    for p in str_keys:
        inp_param[p].set(vary=True)
    return inp_param


def un_vary_params(inp_param, param_str, comp):
    """
    Set all but selected parameters to vary=False
    :param inp_param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    str_keys = [
        key for key, val in inp_param.items() if new_str in key and "tp" not in key
    ]
    for p in inp_param:
        if p not in str_keys:
            inp_param[p].set(vary=False)
    return inp_param


def un_vary_this_params(inp_param, param_str, comp):
    """
    Set all but selected parameters to vary=False
    :param inp_param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    str_keys = [
        key for key, val in inp_param.items() if new_str in key and "tp" not in key
    ]
    for p in str_keys:
        inp_param[p].set(vary=False)
    return inp_param


def un_vary_part_params(inp_param, param_str, comp, order=None):
    """
    Set single component of the selected parameters to vary=False
    :param inp_param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :param order: order of the coefficients. parts missing are set to vary=False
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    str_keys = [
        key for key, val in inp_param.items() if new_str in key and "tp" not in key
    ]
    new_order = int((len(str_keys) - 1) / 2)
    if isinstance(order, list):
        for i in range(new_order):
            if not np.isin(i, order):
                inp_param = un_vary_single_param(inp_param, param_str, comp, 2 * i)
                if i > 0:
                    inp_param = un_vary_single_param(
                        inp_param, param_str, comp, 2 * i - 1
                    )
    return inp_param


def un_vary_single_param(inp_param, param_str, comp, val):
    """
    Set single component of the selected parameters to vary=False
    :param inp_param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :param val: order of coefficients to set to vary=False
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + "_" + comp
    else:
        new_str = param_str
    new_str = new_str + str(val)
    str_keys = [key for key, val in inp_param.items() if new_str in key]
    for p in inp_param:
        if p in str_keys:
            inp_param[p].set(vary=False)
    return inp_param


def peaks_model(
    two_theta,
    azimuth,  # forced to exist as independent values by lmfit
    data_class=None,  # needs to contain conversion factor
    orders=None,  # orders dictionary to get minimum position of the range.
    start_end=[0, 360],
    **params
):
    """Full model of intensities at twotheta and azi given input parameters
    :param two_theta: arr values float
    :param azimuth: arr values float
    :param num_peaks: total number of peaks int
    :param nterms_back: total number of polynomial expansion components for the background int
    :param conv: inputs for the conversion call dict
    :param PenCalc: penalise background or not int
    :param params: lmfit parameter class dict
    :return lmfit model fit result
    """
    # N.B. params now doesn't persist as a parameter class, merely a dictionary, so e.g. call key/value pairs as
    # normal not with '.value'

    # expand the background
    intensity = sf.background_expansion((azimuth, two_theta), orders, params)

    peak_keys = [
        key for key, val in params.items() if "peak" in key and "tp" not in key
    ]
    num_peaks = 0
    n_term_peaks = 0
    while n_term_peaks < len(peak_keys):
        f = sum("peak_" + str(num_peaks) in L for L in peak_keys)
        n_term_peaks = n_term_peaks + f
        num_peaks = num_peaks + 1

    # loop over the number of peaks
    intens_peak = []
    for a in range(int(num_peaks)):

        if "peak_" + str(a) + "_s0" in params:
            symm = params["peak_" + str(a) + "_s0"]
        else:
            symm = 1

        # May need to fix for multiple Fourier coefficients per component
        param_str = "peak_" + str(a)
        comp = "d"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = sf.get_series_type(params, param_str, comp)
        d_all = sf.coefficient_expand(azimuth, parms, coeff_type=coeff_type, start_end=start_end)
        comp = "h"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = sf.get_series_type(params, param_str, comp)
        h_all = sf.coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type, start_end=start_end)
        comp = "w"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = sf.get_series_type(params, param_str, comp)
        w_all = sf.coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type, start_end=start_end)
        comp = "p"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = sf.get_series_type(params, param_str, comp)
        p_all = sf.coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type, start_end=start_end)

        # conversion
        two_theta_all = data_class.conversion(
            d_all, reverse=1
        )  # centroid_conversion(conv, d_all, azimuth)
        intens_peak.append(
            pf.pseudo_voigt_peak(two_theta, two_theta_all, w_all, h_all, p_all)
        )
        intensity = intensity + intens_peak[a]

    # Elapsed time for fitting
    # t_end = time.time()
    # t_elapsed = t_end - t_start
    # print(t_elapsed)

    return intensity


def fit_model(
    data_as_class,  # needs to contain intensity, tth, azi (as chunks), conversion factor
    orders,
    params,
    start_end=[0, 360],
    fit_method="leastsq",
    weights=None,
    max_n_fev=None,
):
    """Initiate model of intensities at twotheta and azi given input parameters and fit
    :param max_n_fev:
    :param intensity_fit: intensity values to fit arr
    :param two_theta: twotheta values arr
    :param azimuth: azimuth values arr
    :param params: lmfit Parameter class dict of parameters to fit
    :param num_peaks: total number of peaks in model int
    :param nterms_back: total number of polynomial component terms in background int
    :param conv: inputs for the conversion call dict
    :param fixed: unsure
    :param method: lmfit choice of fitting method e.g. a least squares
    :param weights: errors on intensity values arr of size intensity_fit
    :return: lmfit model result
    """

    # FIX ME: DMF does the above statement need addressing?
    gmodel = Model(peaks_model, independent_vars=["two_theta", "azimuth"])

    if 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            out = gmodel.fit(
                data_as_class.intensity,  # this is what we are fitting to
                params,  # parameter class to feed in
                two_theta=data_as_class.tth,
                azimuth=data_as_class.azm,
                data_class=data_as_class,  # needs to contain tth, azi, conversion factor
                orders=orders,  # orders class to get peak lengths (if needed)
                start_end=start_end, # start and end of azimuths if needed
                nan_policy="propagate",
                max_n_fev=max_n_fev,
                xtol=1e-5,
            )
    else:
        out = gmodel.fit(
            data_as_class.intensity,  # this is what we are fitting to
            params,  # parameter class to feed in
            two_theta=data_as_class.tth,
            azimuth=data_as_class.azm,
            data_class=data_as_class,  # needs to contain tth, azi, conversion factor
            orders=orders,  # orders class to get peak lengths (if needed)
            start_end=start_end, # start and end of azimuths if needed
            nan_policy="propagate",
            max_n_fev=max_n_fev,
            xtol=1e-5,
        )
    return out


def coefficient_fit(
    ydata,
    azimuth,
    inp_param,
    terms=None,
    errs=None,
    param_str="peak_0",
    symmetry=1,
    fit_method="leastsq",
    start_end = [0,360],
):
    """Fit the Fourier expansion to number of required terms
    :param ydata: Component data array to fit float
    :param azimuth: data array float
    :param inp_param: lmfit Parameter class dict
    :param terms: number of terms to fit int
    :param errs: Component data array errors of size ydata
    :param param_str: str start of parameter name
    :param symmetry: symmetry in azimuth of the fourier to be fit
    :param fit_method: lmfit method default 'leastsq'
    :return: lmfit Model result
    """

    # get NaN values.
    idx = np.isfinite(azimuth) & np.isfinite(ydata)
    if type(terms) == list:
        terms = terms[0]
    if inp_param:
        inp_param = inp_param
    else:
        # FIX ME: DMF assuming this isn't called as it shouldn't work....
        inp_param = Parameters()  # FIX ME: where is this function?
        for t in range(
            2 * terms + 1
        ):  # FIX ME: coefficient series length -- make into separate function.
            inp_param.add(param_str + "_" + str(t), 1.0)  # need limits

    if errs is None or np.array(errs).all is None:
        errs = np.ones(ydata.shape)
    else:
        errs = np.array(errs)

    ydata = np.array(ydata)
    azimuth = np.array(azimuth)

    coeff_type = inp_param.eval(param_str + "_tp")
    f_model = Model(sf.coefficient_expand, independent_vars=["azimuth"], start_end=start_end)
    # print('parameter names: {}'.format(f_model.param_names))
    # print('independent variables: {}'.format(f_model.independent_vars))

    # Attempt to mitigate failure of fit with weights containing 'None' values
    # Replace with nan and alter dtype from object to float64
    new_errs = errs[idx]
    # FIX ME: what happens if all new_errs is None?
    # DMF: I've fudged this for now to test Example 1 from start
    try:
        new_errs[new_errs is None] = 1000 * new_errs[new_errs is not None].max()
    except TypeError:
        new_errs[:] = 0.5  # Need to talk to Simon about this!!!
    new_errs = new_errs.astype("float64")
    out = f_model.fit(
        ydata[idx],
        inp_param,
        azimuth=azimuth[idx] * symmetry,
        coeff_type=coeff_type,
        start_end=start_end,
        method=fit_method,
        sigma=new_errs,
        comp_str=param_str,
        nan_policy="propagate",
    )
    return out
