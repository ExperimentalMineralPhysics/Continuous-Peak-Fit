#!/usr/bin/env python

__all__ = [
    "create_new_params",
    "params_to_new_params",
    "param_fit_array",
    "pseudo_voigt_peak",
    "gather_param_errs_to_list",
    "gather_params_to_list",
    "gather_params_from_dict",
    "initiate_params",
    "peaks_model",
    "vary_params",
    "un_vary_params",
    "fit_model",
    "fourier_expand",
    "fourier_fit",
    "change_params",
    "Fourier_backgrnd",
]

import numpy as np
import numpy.ma as ma
import sys
import warnings
from scipy.interpolate import make_interp_spline, CubicSpline
from lmfit import Parameters, Model
np.set_printoptions(threshold=sys.maxsize)

# Fitting functions required for fitting peaks to the data.
# Called by CPF_XRD_FitSubpattern.

# Known limitations/bugs
# FIX ME: need hard limits to gauss - lorentz peak profile ratios
                

def peak_parameters():
    """
    Lists the parameters needed in each peak and fit order.
    Ideally this will be called and looped over whenever the peaks are referred to
    rather than hardcoding the properties. 
    
    But I dont know how this will works with symmetry...
    
    :return comp_list:
    :return comp_names:
    """
    
    comp_list = [["bg"], ["h", "d", "w", "p", "s"]]
    comp_names = [["background"], ["height", "d-space", "width", "profile", "symmetry"]]
    
    return comp_list, comp_names



# def parse_bounds(bounds, dspace, intens, two_theta, ndat=None, n_peaks=1, param=None):
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
        if par == 'height' or par == 'background':
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


def create_new_params(num_peaks, dfour, hfour, wfour, pfour, bg_four, order_peak):
    """
    Create the new_params dictionary from results
    :param num_peaks: number of peaks int
    :param dfour: list of coefficients per peak float
    :param hfour: list of coefficients per peak float
    :param wfour: list of coefficients per peak float
    :param pfour: list of coefficients per peak float
    :param bg_four:
    :param order_peak:
    :return: new_params dictionary
    """
    peaks = []
    for j in range(num_peaks):
        if (
                "symmetry" in order_peak[j]
        ):  # orders['peak'][j]:  # FIX ME: the symmetry part of this is a horrible way of doing it.
            peaks.append(
                {
                    "d-space": dfour[j][0],
                    "height": hfour[j][0],
                    "width": wfour[j][0],
                    "profile": pfour[j][0],
                    "symmetry": order_peak[j]["symmetry"],
                }
            )
            # "symmetry": orders['peak'][j]['symmetry']})
        else:
            peaks.append(
                {
                    "d-space": dfour[j][0],
                    "height": hfour[j][0],
                    "width": wfour[j][0],
                    "profile": pfour[j][0],
                }
            )
    new_params = {"background": bg_four, "peak": peaks}

    return new_params



def params_to_new_params(params, orders=None):
    """

    Gather Parameters class to new_params dictionary format
    :param params:
    :param num_peaks:
    :param num_bg:
    :param order_peak:
    :return:
    """
    
    num_peaks=len(orders["peak"])
    num_bg   =len(orders["background"])

    peaks = []
    for i in range(num_peaks):
        new_str = "peak_" + str(i) + "_d"
        d_space = gather_param_errs_to_list(params, new_str)
        d_space_tp = get_series_type(params, new_str)
        new_str = "peak_" + str(i) + "_h"
        h_space = gather_param_errs_to_list(params, new_str)
        h_tp = get_series_type(params, new_str)
        new_str = "peak_" + str(i) + "_w"
        w_space = gather_param_errs_to_list(params, new_str)
        w_tp = get_series_type(params, new_str)
        new_str = "peak_" + str(i) + "_p"
        p_space = gather_param_errs_to_list(params, new_str)
        p_tp = get_series_type(params, new_str)

        tmp_peaks = {}
        if "phase" in orders["peak"][i]:
            tmp_peaks["phase"] = orders["peak"][i]["phase"]
        if "hkl" in orders["peak"][i]:
            tmp_peaks["hkl"] = orders["peak"][i]["hkl"]
        tmp_peaks["d-space"] = d_space[0]
        tmp_peaks["d-space_err"] = d_space[1]
        tmp_peaks["d-space-type"] = coefficient_type_as_string(d_space_tp)
        tmp_peaks["height"] = h_space[0]
        tmp_peaks["height_err"] = h_space[1]
        tmp_peaks["height-type"] = coefficient_type_as_string(h_tp)
        tmp_peaks["width"] = w_space[0]
        tmp_peaks["width_err"] = w_space[1]
        tmp_peaks["width-type"] = coefficient_type_as_string(w_tp)
        tmp_peaks["profile"] = p_space[0]
        tmp_peaks["profile_err"] = p_space[1]
        tmp_peaks["profile-type"] = coefficient_type_as_string(p_tp)
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
    bg_tp = coefficient_type_as_string(get_series_type(params, new_str))

    new_params = {
        "background": bg_space,
        "background_err": bg_space_err,
        "background-type": bg_tp,
        "peak": peaks,
    }

    return new_params

# def params_to_new_params(params, num_peaks, num_bg, order_peak):
#     """

#     Gather Parameters class to new_params dictionary format
#     :param params:
#     :param num_peaks:
#     :param num_bg:
#     :param order_peak:
#     :return:
#     """

#     peaks = []
#     for i in range(num_peaks):
#         new_str = "peak_" + str(i) + "_d"
#         d_space = gather_param_errs_to_list(params, new_str)
#         d_space_tp = get_series_type(params, new_str)
#         new_str = "peak_" + str(i) + "_h"
#         h_space = gather_param_errs_to_list(params, new_str)
#         h_tp = get_series_type(params, new_str)
#         new_str = "peak_" + str(i) + "_w"
#         w_space = gather_param_errs_to_list(params, new_str)
#         w_tp = get_series_type(params, new_str)
#         new_str = "peak_" + str(i) + "_p"
#         p_space = gather_param_errs_to_list(params, new_str)
#         p_tp = get_series_type(params, new_str)

#         tmp_peaks = {}
#         if "phase" in order_peak[i]:
#             tmp_peaks["phase"] = order_peak[i]["phase"]
#         if "hkl" in order_peak[i]:
#             tmp_peaks["hkl"] = order_peak[i]["hkl"]
#         tmp_peaks["d-space"] = d_space[0]
#         tmp_peaks["d-space_err"] = d_space[1]
#         tmp_peaks["d-space-type"] = coefficient_type_as_string(d_space_tp)
#         tmp_peaks["height"] = h_space[0]
#         tmp_peaks["height_err"] = h_space[1]
#         tmp_peaks["height-type"] = coefficient_type_as_string(h_tp)
#         tmp_peaks["width"] = w_space[0]
#         tmp_peaks["width_err"] = w_space[1]
#         tmp_peaks["width-type"] = coefficient_type_as_string(w_tp)
#         tmp_peaks["profile"] = p_space[0]
#         tmp_peaks["profile_err"] = p_space[1]
#         tmp_peaks["profile-type"] = coefficient_type_as_string(p_tp)
#         if "symmetry" in order_peak[i]:
#             tmp_peaks["symmetry"] = order_peak[i]["symmetry"]

#         peaks.append(tmp_peaks)

#     # Get background parameters
#     bg_space = []
#     bg_space_err = []
#     for b in range(num_bg):
#         new_str = "bg_c" + str(b) + "_f"
#         bg_spc = gather_param_errs_to_list(params, new_str)
#         bg_space.append(bg_spc[0])
#         bg_space_err.append(bg_spc[1])
#     bg_tp = coefficient_type_as_string(get_series_type(params, new_str))

#     new_params = {
#         "background": bg_space,
#         "background_err": bg_space_err,
#         "background-type": bg_tp,
#         "peak": peaks,
#     }

#     return new_params


def flatten(li):
    """
    :param li:
    :return:
    """
    return sum(([x] if not isinstance(x, list) else flatten(x) for x in li), [])


def param_fit_array(num_p, len_bg, val=None, p=None, b=None, orda=None):
    """
    Makes array which represents all the fourier series in the fit.
    new_arr{0} is a list of 4 element lists - one 4 element list for each peak.
    new_arr[0] is the background representative.
    :param num_p:
    :param len_bg:
    :param val:
    :param p:
    :param b:
    :param orda:
    :return:
    """
    if val:
        v = val
    else:
        v = 0

    if not p:
        p = [v, v, v, v]

    new_arr = [[], []]
    for y in range(num_p):
        new_arr[0].append(p[:])
    for y in range(len_bg):
        new_arr[1].extend([v][:])
    return new_arr


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
    gauss_peak = h_all * np.exp((-((two_theta - two_theta_0) ** 2)) / (2 * w_all ** 2))
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
    l_peak = h_all * w_all ** 2 / ((two_theta - two_theta_0) ** 2 + w_all ** 2)
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
    # FIX ME: how should the ratio be fixed?
    # The ratio is between 1 and 0. it has to be limited.
    # But need options for fixed ratio and not fitting them.
    # How do I do this?
    #
    # I have defined a new variable called profile_fixed which fixes the value if present
    # if np.any(l_g_ratio<0) or np.any(l_g_ratio>1):
    # print 'The ratio is out of range'
    # stop

    p_v_peak = l_g_ratio * gaussian_peak(two_theta, two_theta_0, w_all, h_all) + (
            1 - l_g_ratio
    ) * lorentzian_peak(two_theta, two_theta_0, w_all, h_all)
    return p_v_peak


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
    str_keys = [key for key, val in inp_param.items() if new_str in key and "tp" not in key]
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
    str_keys = [key for key, val in inp_param.items() if new_str in key and "tp" not in key]
    for i in range(len(str_keys)):
        param_list.append(inp_param[new_str + str(i)])
    return param_list



def initiate_all_params_for_fit(settings_as_class, data_as_class, values=None, debug=False):
    """
    Initiate all the parameters needed for the fitting with lmfit.
    
    :param inp_param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: list of parameters in alphanumerical order
    """
    
    comp_list, comp_names = peak_parameters()
    
    
    master_params = Parameters()  # Initiate total Parameter class to add to
    
    lims = parse_bounds(settings_as_class.fit_bounds, data_as_class)#d_space.flatten(), intens.flatten(), twotheta.flatten()
    #)
    
    # for i in range(len(comp_list)):
        
    #     comp_list_section = comp_list[i]
    #     comp_names_section = comp_names[i]
        
    #     for j in range(len(comp_list_section)):
            
    #         comps = comp_list_section[j]
    #         comp_name = comp_names_section[j]
            
    #         if comps[0] == "bg":
                
    
    # Initiate background parameters
    comps = "bg"
    coeff_type = params_get_type(settings_as_class.subfit_orders, comps)
    for k in range(len(settings_as_class.subfit_orders["background"])):
        param_str = "bg_c" + str(k)
        comp = "f"
        if k == 0:
            limits = lims["background"]
        else:
            limits = None
        n_coeff = get_number_coeff(
            settings_as_class.subfit_orders, comps, peak=k, azimuths=data_as_class.azm
        )
        # add values is present
        vals=None
        if values:
            vals = values["background"][k],
        #intiate parameters
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
        
        # master_params.pretty_print()
        # master_params = un_vary_params(
        #     master_params, param_str, comp
        # )  # set other parameters to not vary
        # master_params = vary_params(
        #     master_params, param_str, comp
        # )  # set these parameters to vary
        # if isinstance(
        #     orders["background"][b], list
        # ):  # set part of these parameters to not vary
        #     master_params = pf.un_vary_part_params(
        #         master_params, param_str, comp, orders["background"][b]
        #     )
        # print(new_bg_all[b])
        # fout = pf.coefficient_fit(
        #     azimuth=np.array(new_azi_chunks),
        #     ydata=np.array(new_bg_all[b]),
        #     inp_param=master_params,
        #     param_str=param_str + "_" + comp,
        #     symmetry=1,
        #     errs=np.array(new_bg_all_err[b]),
        #     fit_method="leastsq",
        # )
        # master_params = fout.params
            
        
    # initiate peak(s)
    
    # comp_list, comp_names = peak_parameters()
    
    peak_comp_list = comp_list[1]
    comp_names_list = comp_names[1]
    
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
        # arr_names = ["new_h_all", "new_d0", "new_w_all", "new_p_all"]
        # arr_err_names = [
        #     "new_h_all_err",
        #     "new_d0_err",
        #     "new_w_all_err",
        #     "new_p_all_err",
        # ]
    
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
                
                
            coeff_type = params_get_type(settings_as_class.subfit_orders, comp, peak=j)
            n_coeff = get_number_coeff(
                settings_as_class.subfit_orders, comp, peak=j, azimuths=data_as_class.azm
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
            # master_params = pf.un_vary_params(
            #     master_params, param_str, comp
            # )  # set other parameters to not vary
            # if fixed == 0:
            #     master_params = pf.vary_params(
            #         master_params, param_str, comp
            #     )  # set these parameters to vary
            #     if isinstance(
            #         # FIX ME: changed k to j here as think bug but needs checking!
            #         orders["peak"][j][comp_names[cp]], list
            #     ):  # set part of these parameters to not vary
            #         # FIX ME: changed k to j here as think bug but needs checking!
            #         master_params = pf.un_vary_part_params(
            #             master_params,
            #             param_str,
            #             comp,
            #             orders["peak"][j][comp_names[cp]],
            #         )
            #     if n_coeff > len(
            #         np.unique(new_azi_chunks)
            #     ):  # catch me make sure there are not more coefficients than chunks
            #         o = int(np.floor(len(np.unique(new_azi_chunks)) / 2 - 1))
            #         un_vary = [x for x in range(0, o)]
            #         master_params = pf.un_vary_part_params(
            #             master_params, param_str, comp, un_vary
            #         )
            #     fout = pf.coefficient_fit(
            #         azimuth=np.array(new_azi_chunks),
            #         ydata=vals,
            #         inp_param=master_params,
            #         param_str=param_str + "_" + comp,
            #         symmetry=symmetry,
            #         errs=vals_err,
            #         fit_method="leastsq",
            #     )
    
            # master_params = fout.params
            
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
        new_max, new_min = limits
        half_range = (new_max - new_min) / 2
    else:
        new_min = -np.inf
        new_max = np.inf
        half_range = np.inf
    if value is None and limits is None:
        value = 0.01
    elif value is None and limits:
        value = new_min + (new_max - new_min) * 3 / 5  # /2#*2/3
        # N.B. it seems the initiating the splines with values exactly half-way between man and min doesn't work.
        # Anywhere else in the middle of the range seems to do so.
        # FIX ME: Is this fitting-related or a code bug?
    value = np.array(value)  # force input to be array so that can be iterated over

    coeff_type = coefficient_type_as_number(coeff_type)
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
            coefficient_type_as_number(coeff_type),
            expr=expr,
            vary=False,
        )

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
    str_keys = [key for key, val in inp_param.items() if new_str in key and "tp" not in key]
    for p in inp_param:
        if p not in str_keys:
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
    str_keys = [key for key, val in inp_param.items() if new_str in key and "tp" not in key]
    new_order = int((len(str_keys) - 1) / 2)
    if isinstance(order, list):
        for i in range(new_order):
            if not np.isin(i, order):
                inp_param = un_vary_single_param(inp_param, param_str, comp, 2 * i)
                if i > 0:
                    inp_param = un_vary_single_param(inp_param, param_str, comp, 2 * i - 1)
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
    str_keys = [key for key, val in inp_param.items() if new_str in key and "tp" not in key]
    for p in str_keys:
        inp_param[p].set(vary=True)
    return inp_param


# def peaks_model(two_theta, azimuth, conv=None, **params):
#     """Full model of intensities at twotheta and azi given input parameters
#     :param two_theta: arr values float
#     :param azimuth: arr values float
#     :param num_peaks: total number of peaks int
#     :param nterms_back: total number of polynomial expansion components for the background int
#     :param conv: inputs for the conversion call dict
#     :param PenCalc: penalise background or not int
#     :param params: lmfit parameter class dict
#     :return lmfit model fit result
#     """
#     # N.B. params now doesn't persist as a parameter class, merely a dictionary, so e.g. call key/value pairs as
#     # normal not with '.value'

#     # recreate background array to pass to fourier_background
#     background_keys = [
#         key for key, val in params.items() if "bg_c" in key and "tp" not in key
#     ]
#     n_term_fourier = []

#     i = 0
#     backg = []
#     backg_tp = []
#     while sum(n_term_fourier) < len(background_keys):

#         f = sum("bg_c" + str(i) in L for L in background_keys)
#         n_term_fourier.append(f)
#         fourier_background = []
#         for k in range(f):
#             fourier_background.append(params["bg_c" + str(i) + "_f" + str(k)])
#             if "bg_c" + str(i) + "_f_tp" in params:
#                 b_tp = params["bg_c" + str(i) + "_f_tp"]
#             else:
#                 b_tp = 0

#         backg.append(np.array(fourier_background))
#         backg_tp.append(b_tp)
#         i = i + 1

#     # for future logging
#     # print('recreate backg to pass to fourier_backg\n', 'background_keys is ', background_keys)
#     # print('backg is ', backg)
#     intensity = background_fourier_expansion((azimuth, two_theta), backg, coeff_type=backg_tp)

#     peak_keys = [
#         key for key, val in params.items() if "peak" in key and "tp" not in key
#     ]
#     num_peaks = 0
#     n_term_peaks = 0
#     while n_term_peaks < len(peak_keys):
#         f = sum("peak_" + str(num_peaks) in L for L in peak_keys)
#         n_term_peaks = n_term_peaks + f
#         num_peaks = num_peaks + 1

#     # loop over the number of peaks
#     intens_peak = []
#     for a in range(int(num_peaks)):

#         if "peak_" + str(a) + "_s0" in params:
#             symm = params["peak_" + str(a) + "_s0"]
#         else:
#             symm = 1

#         # May need to fix for multiple Fourier coefficients per component
#         param_str = "peak_" + str(a)
#         comp = "d"
#         parms = gather_params_from_dict(params, param_str, comp)
#         coeff_type = get_series_type(params, param_str, comp)
#         d_all = coefficient_expand(azimuth, parms, coeff_type=coeff_type)
#         comp = "h"
#         parms = gather_params_from_dict(params, param_str, comp)
#         coeff_type = get_series_type(params, param_str, comp)
#         h_all = coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type)
#         comp = "w"
#         parms = gather_params_from_dict(params, param_str, comp)
#         coeff_type = get_series_type(params, param_str, comp)
#         w_all = coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type)
#         comp = "p"
#         parms = gather_params_from_dict(params, param_str, comp)
#         coeff_type = get_series_type(params, param_str, comp)
#         p_all = coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type)

#         # conversion
#         two_theta_all = centroid_conversion(conv, d_all, azimuth)
#         intens_peak.append(pseudo_voigt_peak(two_theta, two_theta_all, w_all, h_all, p_all))
#         intensity = intensity + intens_peak[a]

#     # Elapsed time for fitting
#     # t_end = time.time()
#     # t_elapsed = t_end - t_start
#     # print(t_elapsed)

#     return intensity


def peaks_model(two_theta, azimuth, #forced to exist as independent values by lmfit
        data_class = None, # needs to contain conversion factor
        orders = None, # orders dictionary to get minimum position of the range.
        **params):
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
    intensity = background_expansion((azimuth, two_theta), orders, params)

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
        coeff_type = get_series_type(params, param_str, comp)
        d_all = coefficient_expand(azimuth, parms, coeff_type=coeff_type)
        comp = "h"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = get_series_type(params, param_str, comp)
        h_all = coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type)
        comp = "w"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = get_series_type(params, param_str, comp)
        w_all = coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type)
        comp = "p"
        parms = gather_params_from_dict(params, param_str, comp)
        coeff_type = get_series_type(params, param_str, comp)
        p_all = coefficient_expand(azimuth * symm, parms, coeff_type=coeff_type)

        # conversion
        two_theta_all = data_class.conversion(d_all, reverse=1)# centroid_conversion(conv, d_all, azimuth)
        intens_peak.append(pseudo_voigt_peak(two_theta, two_theta_all, w_all, h_all, p_all))
        intensity = intensity + intens_peak[a]

    # Elapsed time for fitting
    # t_end = time.time()
    # t_elapsed = t_end - t_start
    # print(t_elapsed)

    return intensity

# # Not needed?
# def change_params(constants, *fit_param):
#     print("In change param", len(constants))
#     two_theta = constants[0]
#     azimuth = constants[1]
#     change_array = constants[2]
#     shapes = constants[3]
#     conv = constants[4]
#     symm = constants[5]
#     fixed = constants[6]

#     print("change_params:")
#     print(two_theta, azimuth)
#     print(change_array)
#     print(shapes)
#     print(conv)
#     print(symm)
#     print(fixed)
#     print(fit_param)
#     shapes = GuessApply(change_array, shapes, fit_param)
#     # print(stop)
#     print("In changes", shapes)
#     I, pen = PeaksModel_old(two_theta, azimuth, shapes, conv, symm, PenCalc=1)

#     return I  # *pen


def fit_model(
        data_as_class, # needs to contain intensity, tth, azi (as chunks), conversion factor
        orders,
        params,
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
                data_as_class.intensity, # this is what we are fitting to
                params, # parameter class to feed in
                two_theta=data_as_class.tth,
                azimuth=data_as_class.azm,
                data_class = data_as_class, # needs to contain tth, azi, conversion factor
                orders = orders, # orders class to get peak lengths (if needed)
                nan_policy="propagate",
                max_n_fev=max_n_fev,
                xtol=1e-5,
            )
    else:
        out = gmodel.fit(
            data_as_class.intensity, # this is what we are fitting to
            params, # parameter class to feed in
            two_theta=data_as_class.tth,
            azimuth=data_as_class.azm,
            data_class = data_as_class, # needs to contain tth, azi, conversion factor
            orders = orders, # orders class to get peak lengths (if needed)
            nan_policy="propagate",
            max_n_fev=max_n_fev,
            xtol=1e-5,
        )
    return out

# def fit_model(
#         intensity_fit,
#         two_theta,
#         azimuth,
#         params,
#         num_peaks,
#         nterms_back,
#         conv=None,
#         fixed=None,
#         fit_method="leastsq",
#         weights=None,
#         max_n_fev=None,
# ):
#     """Initiate model of intensities at twotheta and azi given input parameters and fit
#     :param max_n_fev:
#     :param intensity_fit: intensity values to fit arr
#     :param two_theta: twotheta values arr
#     :param azimuth: azimuth values arr
#     :param params: lmfit Parameter class dict of parameters to fit
#     :param num_peaks: total number of peaks in model int
#     :param nterms_back: total number of polynomial component terms in background int
#     :param conv: inputs for the conversion call dict
#     :param fixed: unsure
#     :param method: lmfit choice of fitting method e.g. a least squares
#     :param weights: errors on intensity values arr of size intensity_fit
#     :return: lmfit model result
#     """

#     # bounds passed but not used
#     # FIX ME: DMF does the above statement need addressing?
#     gmodel = Model(peaks_model, independent_vars=["two_theta", "azimuth"])
#     # print('parameter names: {}'.format(gmodel.param_names))
#     # print('independent variables: {}'.format(gmodel.independent_vars))

#     out = gmodel.fit(
#         intensity_fit,
#         params,
#         two_theta=two_theta,
#         azimuth=azimuth,
#         conv=conv,
#         nan_policy="propagate",
#         max_n_fev=max_n_fev,
#         xtol=1e-5,
#     )
#     return out


# def centroid_conversion(conv, args_in, azimuth):
#     """

#     :param conv:
#     :param args_in:
#     :param azimuth:
#     :return:
#     """
#     # Conv['DispersionType'] is the conversion type
#     # Conv[...] are the values required for the conversion
#     # FIX ME: these functions should be a sub function of the detector types. but need to work out how to do it.
#     if conv == 0 or conv["DispersionType"] is None:
#         args_out = args_in

#     elif conv["DispersionType"] == "AngleDispersive":
#         # args_in are d-spacing
#         # args_out is two thetas
#         args_out = 2 * np.degrees(
#             np.arcsin(conv["conversion_constant"] / 2 / args_in)
#         )  # FIX ME: check this is correct!!!

#     elif conv["DispersionType"] == "EnergyDispersive":
#         args_out = []
#         for x in range(azimuth.size):
#             # determine which detector is being used

#             if azimuth.size == 1:
#                 a = np.array(np.where(conv["azimuths"] == azimuth)[0][0])
#                 args_out.append(
#                     12.398
#                     / (
#                             2
#                             * args_in
#                             * np.sin(
#                         np.radians(conv["calibs"].mcas[a].calibration.two_theta / 2)
#                     )
#                     )
#                 )
#             else:
#                 try:
#                     a = np.array(np.where(conv["azimuths"] == azimuth[x])[0][0])
#                     args_out.append(
#                         12.398
#                         / (
#                                 2
#                                 * args_in[x]
#                                 * np.sin(
#                             np.radians(
#                                 conv["calibs"].mcas[a].calibration.two_theta / 2
#                             )
#                         )
#                         )
#                     )
#                 # FIX ME: needs handling properly
#                 except:
#                     args_out.append(0)
#         if isinstance(azimuth, np.ma.MaskedArray):
#             args_out = ma.array(args_out, mask=azimuth.mask)
#     else:
#         raise ValueError("Unrecognised conversion type")

#     return args_out


def expand_comp_string(comp):
    if comp == "d":
        out = "d-space"
    elif comp == "h":
        out = "height"
    elif comp == "w":
        out = "width"
    elif comp == "p":
        out = "profile"
    elif comp == "bg" or "f":
        out = "background"
    else:
        raise ValueError("Unrecognised peak property type")
    return out


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
    comp_str = expand_comp_string(comp) + "-type"
    if comp_str in orders:
        coeff_type = orders[comp_str]
    elif len(orders["peak"]) > peak and comp_str in orders["peak"][peak]:
        coeff_type = orders["peak"][peak][comp_str]
    else:
        coeff_type = "fourier"

    return coeff_type


def get_number_coeff(orders, comp, peak=0, azimuths=None):
    """

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
        n_param = np.max(orders["background"][peak]) * 2 + 1

    else:  # everything else.
        n_param = np.max(orders["peak"][peak][expand_comp_string(comp)]) * 2 + 1

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
            l = np.max(params["peak"][peak][expand_comp_string(comp)])
    elif isinstance(params, (int,)):
        l = 1
    else:
        print(params)
        print(type(params))
        raise ValueError("Parameter list is not list or float.")

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
        print(params)
        print(type(params))
        raise ValueError("Parameter list is not list or float.")
    return order


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
    elif coeff_type == "spline_quadratic" or coeff_type == "quadratic" or coeff_type == 2:
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
        if return_error==0:
            raise ValueError(
                error_str
            )
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
        out = 5
    else:
        raise ValueError(
            "Unrecognised coefficient series type, the valid options are "
            "fourier"
            ", etc..."
        )
        # FIX ME: write out all the licit options in the error message.
    return out


def coefficient_expand(azimuth, param=None, coeff_type="fourier", comp_str=None, **params):
    """

    :param azimuth:
    :param param:
    :param coeff_type:
    :param comp_str:
    :param params:
    :return:
    """

    coeff_type = coefficient_type_as_number(coeff_type)
    se = [0, 360]

    if coeff_type == 0:
        out = fourier_expand(azimuth, inp_param=param, comp_str=comp_str, **params)
    elif coeff_type == 1:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=se,
            bc_type="natural",
            kind="linear",
            **params
        )
    elif coeff_type == 2:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=se,
            bc_type="natural",
            kind="quadratic",
            **params
        )
    elif coeff_type == 3:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=se,
            bc_type="periodic",
            kind="cubic",
            **params
        )
    elif coeff_type == 4:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=se,
            bc_type="natural",
            kind="cubic",
            **params
        )
    elif coeff_type == 5:
        out = spline_expand(
            azimuth,
            inp_param=param,
            comp_str=comp_str,
            start_end=se,
            bc_type="natural",
            kind="independent",
            **params
        )
    else:
        raise ValueError(
            "Unrecognised coefficient series type, the valid options are "
            "fourier"
            ", etc..."
        )
        # FIX ME: write out all the licit options in the error message.
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
            params.add(param_str + "_" + str(t), 1.0)  # need limits

    if errs is None or np.array(errs).all is None:
        errs = np.ones(ydata.shape)
    else:
        errs = np.array(errs)
        
    ydata = np.array(ydata)
    azimuth = np.array(azimuth)

    coeff_type = inp_param.eval(param_str + "_tp")
    f_model = Model(coefficient_expand, independent_vars=["azimuth"])
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
        method=fit_method,
        sigma=new_errs,
        comp_str=param_str,
        nan_policy="propagate",
    )
    return out


# spline expansion function
def spline_expand(
        azimuth,
        inp_param=None,
        comp_str=None,
        start_end=[0, 360],
        bc_type="periodic",
        kind="cubic",
        **params
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
            spl = CubicSpline(points, inp_param, bc_type=bc_type, extrapolate="periodic")
        else:
            spl = make_interp_spline(points, inp_param, k=k)

        fout = spl(azimuth)

    return fout


# fourier expansion function
def fourier_expand(azimuth, inp_param=None, comp_str=None, **params):
    """Calculate Fourier expansion given input coefficients
    :param azimuth: arr data array float
    :param inp_param: list of Fourier coefficients float
    :param comp_str: str to determine which coefficients to use from params
    :param params: lmfit dict of coefficients as parameters
    :return:
    """
    # print('param, fourier expand', param)
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
        for i in range(1, int((len(inp_param) - 1) / 2) + 1):
            # len(param)-1 should never be odd because of initial a_0 parameter
            fout = (
                    fout
                    + inp_param[(2 * i) - 1] * np.sin(np.deg2rad(azimuth) * i)
                    + inp_param[2 * i] * np.cos(np.deg2rad(azimuth) * i)
            )  # single col array
    return fout


# def fourier_fit(
#         ydata,
#         azimuth,
#         param,
#         terms=None,
#         errs=None,
#         param_str="peak_0",
#         symmetry=1,
#         fit_method="leastsq",
# ):
#     """Fit the Fourier expansion to number of required terms
#     :param ydata: Component data array to fit float
#     :param azimuth: data array float
#     :param param: lmfit Parameter class dict
#     :param terms: number of terms to fit int
#     :param errs: Component data array errors of size ydata
#     :param param_str: str start of parameter name
#     :param symmetry: symmetry in azimuth of the fourier to be fit
#     :param fit_method: lmfit method default 'leastsq'
#     :return: lmfit Model result
#     """

#     # get NaN values.
#     idx = np.isfinite(azimuth) & np.isfinite(ydata)
#     if type(terms) == list:
#         terms = terms[0]
#     if param:
#         param = param
#     else:
#         # FIX ME: DMF assuming this isn't called as it shouldn't work!
#         param = Parameters()
#         for t in range(2 * terms + 1):
#             params.add(param_str + "_" + str(t), 1.0)  # need limits

#     if errs is None or errs.all is None:
#         errs = np.ones(ydata.shape)

#     # param.pretty_print()
#     f_model = Model(fourier_expand, independent_vars=["azimuth"])
#     # print('parameter names: {}'.format(f_model.param_names))
#     # print('independent variables: {}'.format(f_model.independent_vars))

#     # Attempt to mitigate failure of fit with weights containing 'None' values
#     # Replace with nan and alter dtype from object to float64
#     new_errs = errs[idx]
#     new_errs[new_errs is None] = 1000 * new_errs[new_errs is not None].max()
#     new_errs = new_errs.astype("float64")
#     out = f_model.fit(
#         ydata[idx],
#         param,
#         azimu=azimuth[idx] * symmetry,
#         method=fit_method,
#         sigma=new_errs,
#         comp_str=param_str,
#         nan_policy="propagate",
#     )
#     return out


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

    # recreate background array to pass to fourier_background
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
