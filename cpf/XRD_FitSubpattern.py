#!/usr/bin/env python

__all__ = ["fit_sub_pattern"]

# CPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from lmfit import Parameters, Model
from lmfit.model import save_modelresult  # , load_modelresult
import json
import cpf.PeakFunctions as pf
import cpf.IO_functions as io


np.set_printoptions(threshold=sys.maxsize)


# FIX ME;
# The script should not need to know what the limits are - it should only see the data that it needs to fit.
def run_initial_fit(
    azimuth,
    y_data,
    y_data_errs,
    params,
    param_str,
    comp,
    order,
    extras,
    coeff_type=0,
    method="leastsq",
    symmetry=1,
    value=None,
):
    """
    Initiate parameter class and perform first Fourier fit
    :param coeff_type:
    :param azimuth: x-data arr float
    :param y_data: y_data arr float
    :param y_data_errs: y_data error arr float
    :param params: lmfit Parameter object
    :param param_str: base string to identify components str
    :param comp: string to add to base string to identify components str
    :param order: determines number of Fourier coefficients int
    :param method: lmfit method to use str
    :param extras: list of [[min, max, expr], ind_vars], expression for bounds on parameters and independent variable
    added to params for expression
    :param symmetry:
    :param value:
    :return: lmfit Parameter class
    """
    new_limits = extras
    new_min, new_max, expr = new_limits

    if isinstance(order, list):
        new_order = max(order)
    else:
        new_order = order

    if value is None:
        value = np.zeros(new_order * 2 + 1)
        value[0] = np.mean(y_data)

    params = pf.initiate_params(
        params,
        param_str,
        comp,
        coeff_type,
        new_order,
        limits=[new_max, new_min],
        expr=expr,
        ind_vars=azimuth,
        value=value,
    )
    # set other parameters to not vary
    params = pf.un_vary_params(
        params, param_str, comp
    )
    # set part of parameter to not vary
    params = pf.un_vary_part_params(
        params, param_str, comp, order
    )
    params.pretty_print()
    fout = pf.fourier_fit(
        azimuth=azimuth,
        ydata=np.array(y_data),
        param=params,
        errs=np.array(y_data_errs),
        param_str=param_str + "_" + comp,
        symmetry=symmetry,
        fit_method=method,
    )
    params = fout.params
    return params

# def component_fit():


def update_component_fits(
    master_params,
    intens,
    azimu,
    twotheta,
    param_str,
    comp_lists,
    conversion_factor,
    bgguess,
    guesses,
    o=None,
    pfixed=None,
):
    """
    Perform component fits to data
    :param o:
    :param master_params: lmfit Parameter class object
    :param intens: intensity data array
    :param azimu: data array
    :param twotheta: data array
    :param param_str: start string to label parameters
    :param comp_lists: list of lists for parameters
    :param conversion_factor: dict inputs for conversion
    :param bgguess: initial values or background dict
    :param guesses: initial values for peak components' dict
    :param pfixed: fix p or not
    :return: list of lists with component coefficients
    """

    dfour, hfour, wfour, pfour = comp_lists

    comp = "h"
    # set other parameters to not vary
    master_params = pf.unvary_params(master_params, param_str, comp)
    # set these parameters to vary
    master_params = pf.vary_params(master_params, param_str, comp)
    # set part of these parameters to not vary
    master_params = pf.unvary_part_params(master_params, param_str, comp, o["height"])
    #master_params.pretty_print()
    out = pf.fit_model(
        intens.flatten(),
        twotheta.flatten(),
        azimu.flatten(),
        num_peaks=len(guesses["peak"]),
        nterms_back=len(bgguess),
        conv=conversion_factor,
        fixed=pfixed,
        fit_method=None,
        weights=None,
        params=master_params,
    )
    # master_params = out.params
    hfour.append(pf.gather_param_errs_to_list(master_params, param_str, comp))
    #master_params.pretty_print()

    comp = "d"
    # set other parameters to not vary
    master_params = pf.unvary_params(master_params, param_str, comp)
    # set these parameters to vary
    master_params = pf.vary_params(master_params, param_str, comp)
    # set part of these parameters to not vary
    master_params = pf.unvary_part_params(master_params, param_str, comp, o["d-space"])
    out = pf.fit_model(
        intens.flatten(),
        twotheta.flatten(),
        azimu.flatten(),
        num_peaks=len(guesses["peak"]),
        nterms_back=len(bgguess),
        conv=conversion_factor,
        fixed=pfixed,
        fit_method=None,
        weights=None,
        params=master_params,
    )
    # master_params = out.params
    dfour.append(pf.gather_param_errs_to_list(master_params, param_str, comp))

    comp = "w"
    # set other parameters to not vary
    master_params = pf.unvary_params(master_params, param_str, comp)
    # set these parameters to vary
    master_params = pf.vary_params(master_params, param_str, comp)
    # set part of these parameters to not vary
    master_params = pf.unvary_part_params(master_params, param_str, comp, o["width"])
    out = pf.fit_model(
        intens.flatten(),
        twotheta.flatten(),
        azimu.flatten(),
        num_peaks=len(guesses["peak"]),
        nterms_back=len(bgguess),
        conv=conversion_factor,
        fixed=pfixed,
        fit_method=None,
        weights=None,
        params=master_params,
    )
    # master_params = out.params
    wfour.append(pf.gather_param_errs_to_list(master_params, param_str, comp))

    comp = "p"
    if pfixed is False:
        # if 'profile_fixed' not in order_peak[0]: # orders['peak'][0]:
        # set other parameters to not vary
        master_params = pf.unvary_params(master_params, param_str, comp)
        # set these parameters to vary
        master_params = pf.vary_params(master_params, param_str, comp)
        # set part of these parameters to not vary
        master_params = pf.unvary_part_params(master_params, param_str, comp, o["profile"])
        out = pf.fit_model(
            intens.flatten(),
            twotheta.flatten(),
            azimu.flatten(),
            num_peaks=len(guesses["peak"]),
            nterms_back=len(bgguess),
            conv=conversion_factor,
            fixed=pfixed,
            fit_method=None,
            weights=None,
            params=master_params,
        )
        # master_params = out.params
    pfour.append(pf.gather_param_errs_to_list(master_params, param_str, comp))

    return master_params, dfour, hfour, wfour, pfour


def order_set_peaks(orders, peeks, bg_order):
    """

    :param bg_order:
    :param orders:
    :param peeks:
    :return:
    """

    backgnd = []

    for j in bg_order:
        backgnd.append(np.zeros(np.max(j) * 2 + 1))
    for y in range(peeks):
        # force profile orders to match the fixed orders
        # FIX ME: should this generate an error?
        if "profile_fixed" in orders["peak"][y]:
            if (
                not pf.fourier_order(orders["peak"][y]["profile_fixed"])
                    == orders["peak"][y]["profile"]
            ):
                print(
                    "Peak "
                    + str(y)
                    + ": The order of the profile does not match that of the fixed profile. "
                    "Changing profile to match fixed profile."
                )
                orders["peak"][y]["profile"] = pf.fourier_order(
                    orders["peak"][y]["profile_fixed"]
                )
            # make sure the fixed profile is a list
            if not isinstance(orders["peak"][y]["profile_fixed"], list):
                orders["peak"][y]["profile_fixed"] = [
                    orders["peak"][y]["profile_fixed"]
                ]
    return orders, backgnd


def update_previous_params_from_orders(peeks, previous_params, orders):
    """

    :param peeks:
    :param previous_params:
    :param orders:
    :return:
    """
    choice_list = ["d-space", "height", "width", "profile"]
    # FIX ME: need to check sizes of both arrays here.
    # orders take precedence over previous_params, so we add or remove coefficients as necessary
    for y in range(peeks):
        # loop over parameters
        for param in choice_list:
            coeff_type = pf.params_get_type(previous_params, param, peak=y)
            if coeff_type != 5:  # if parameters are not independent
                if np.max(orders["peak"][y][param]) > pf.fourier_order(
                        previous_params["peak"][y][param]
                ):
                    change_by = (
                            np.max(orders["peak"][y][param]) * 2
                            + 1
                            - np.size(previous_params["peak"][y][param])
                    )
                    previous_params["peak"][y][param] = (
                            previous_params["background"][y] + [0] * change_by
                    )
                elif np.max(orders["peak"][y][param]) < pf.fourier_order(
                        previous_params["peak"][y][param]
                ):
                    # print 'smaller'
                    change_by = (np.size(previous_params["peak"][y][param]) - 
                            np.max(orders["peak"][y][param]) * 2 + 1)
                    previous_params["peak"][y][param] = previous_params["peak"][y][
                                                            param
                                                        ][1: -(change_by - 1)]
        if "profile_fixed" in orders["peak"][y]:
            previous_params["peak"][y]["profile"] = orders["peak"][y][
                "profile_fixed"
            ]
        elif "profile_fixed" in previous_params["peak"][y][param]:
            del previous_params["peak"][y]["profile_fixed"]

    # loop for background orders/size
    if (
            pf.params_get_type(previous_params, "background") != 5
    ):  # if parameters are not independent
        for y in range(
                np.max([len(orders["background"]), len(previous_params["background"])])
        ):
            # loop over the degree of the background
            if (
                    len(previous_params["background"]) - 1 >= y
                    and len(orders["background"]) - 1 >= y
            ):
                # if present in both arrays make sure it is the right size.
                if np.max(orders["background"][y]) > pf.fourier_order(
                        previous_params["background"][y]
                ):
                    # print 'longer'
                    change_by = (
                            np.max(orders["background"][y]) * 2
                            + 1
                            - np.size(previous_params["background"][y])
                    )
                    previous_params["background"][y] = (
                            previous_params["background"][y] + [0] * change_by
                    )
                elif np.max(orders["background"][y]) < pf.fourier_order(
                        previous_params["background"][y]
                ):
                    # print 'smaller'
                    change_by = np.size(previous_params["background"][y]) - (
                            np.max(orders["background"][y]) * 2 + 1
                    )
                    previous_params["background"][y] = previous_params["background"][
                                                           y
                                                       ][:-change_by]
            elif (
                    len(previous_params["background"]) - 1
                    >= y
                    > len(orders["background"]) - 1
            ):
                # if previous params has too many degrees remove the higher ones
                previous_params["background"][y] = []
            elif (
                    len(previous_params["background"]) - 1
                    < y
                    <= len(orders["background"]) - 1
            ):
                # if previous parameters does not have enough degrees add more.
                previous_params["background"].append(
                    [0] * (np.max(orders["background"][y]) * 2 + 1)
                )

    # tidy up from reductions in length of PreviousParams
    previous_params["background"] = [
        x for x in previous_params["background"] if x != []
    ]
    # FIX ME: NO backg_type set here!
    return previous_params


def check_num_azimuths(peeks, azimu, orders):
    """

    :param peeks:
    :param azimu:
    :param orders:
    :return:
    """
    n_azimuths = len(np.unique(azimu))
    max_coeff = 0
    choice_list = ["d-space", "height", "width", "profile"]
    for y in range(peeks):
        # loop over parameters
        for param in choice_list:
            coeff_type = pf.params_get_type(orders, param, peak=y)
            if coeff_type != 5:  # if parameters are not independent
                max_coeff = np.max([max_coeff, np.max(orders["peak"][y][param]) * 2 + 1])
    param = "background"
    for y in range(np.max([len(orders["background"])])):
        coeff_type = pf.params_get_type(orders, param, peak=y)
        if coeff_type != 5:  # if parameters are not independent
            max_coeff = np.max([max_coeff, np.max(orders["background"][y]) * 2 + 1])
    #print(max_coeff, len(np.unique(azimu)))
    if max_coeff > len(np.unique(azimu)):
        err_str = (
                "The maximum order, %i, needs more coefficients than the number of unique azimuths, %i. "
                "It is not possible to fit."
                % (pf.get_order_from_coeff(max_coeff), len(np.unique(azimu)))
        )
        raise ValueError(err_str)


def get_manual_guesses(peeks, orders, bounds, twotheta, debug=None):
    """

    :param debug:
    :param twotheta:
    :param bounds:
    :param peeks:
    :param orders:
    :return:
    """
    t_th_guesses = np.array(orders["PeakPositionSelection"])

    #confirm there are the same number of peaks selected as are to be fit.
    if len(t_th_guesses) != peeks:
        raise ValueError("The number of peaks and the postion selection do not match.")

    # for future use in setting limits
    dist = np.max(t_th_guesses[:, 2]) - np.min(t_th_guesses[:, 2])
    width_change = dist / (2 ** peeks)
    peak_vals = np.unique(t_th_guesses[:, 2])

    # Fit Fourier series to two-theta/d-spacing for each peak
    # FIX ME: this fitting should be done to the d-spacing not the tth/energy. Probably need to import
    # the detector functions to make this happen.
    dfour = []
    for j in range(peeks):
        peek = j + 1
        t_th_guess = t_th_guesses[t_th_guesses[:, 0] == peek, 1:]
        param_str = "peak_" + str(j)
        comp = "d"
        lims = pf.parse_bounds(
            bounds, twotheta.flatten(), 0, 0, param=["d-space"]
        )
        # Confirm that peak position selections are within the bounds.
        # If not fail as gracefully as possible.
        lims_tth = np.array(lims["d-space"])
        if (
                np.min(t_th_guess[:, 1]) < lims_tth[0]
                or np.max(t_th_guess[:, 1]) > lims_tth[1]
        ):
            raise ValueError(
                "At least one of the "
                "PeakPositionSelection"
                " values is outside of the bounds. Check values in "
                "range"
                " and "
                "PeakPositionSelection"
                "."
            )

        if "d-space-type" in orders["peak"][j]:
            coeff_type = orders["peak"][j]["d-space-type"]
        else:
            coeff_type = "fourier"

        # for guesses make sure there are not too many coefficients.
        n_coeff = pf.get_number_coeff(orders, comp)
        if n_coeff > len(t_th_guess[:, 0]):
            o = int(np.floor(len(t_th_guess[:, 0]) / 2 - 1))
            # FIX ME: this should be a function!! and should check the series type.
        else:
            o = orders["peak"][j]["d-space"]

        temp_param = Parameters()
        temp_param = pf.initiate_params(
            inp_param=temp_param,
            param_str=param_str,
            coeff_type=coeff_type,
            comp=comp,
            trig_orders=o,
            limits=lims["d-space"],
            value=t_th_guess[1, 1],
        )
        fout = pf.coefficient_fit(
            azimuth=t_th_guess[:, 0],
            ydata=t_th_guess[:, 1],
            inp_param=temp_param,
            param_str=param_str + "_" + comp,
            fit_method="leastsq",
        )
        temp_param = fout.params
        if debug:
            temp_param.pretty_print()
        dfour.append(
            pf.gather_param_errs_to_list(temp_param, param_str, comp)
        )
    return dfour


def get_chunk_background_guess(twotheta, azimu, intens, orders, chunks, background_type, count, n):
    """

    :param n:
    :param count:
    :param intens:
    :param azimu:
    :param background_type:
    :param twotheta:
    :param orders:
    :param chunks:
    :return:
    """
    j = count
    # Get indices of sorted two theta values excluding the masked values
    tth_ord = ma.argsort(twotheta.flatten()[chunks[j]].compressed())
    background_guess = [[0.0] for i in range(len(orders["background"]))]
    if background_type == "coeffs":
        # Fourier expands the coefficients.
        # FIX ME! Replace fourier_expand with coeff_expand or remove background_type == "coeffs" as an option
        for k in range(len(orders["background"])):
            background_guess[j] = pf.fourier_expand(
                np.mean(azimu.flatten()[chunks[j]]),
                inp_param=orders["peak"]["background"][k],
            )
    elif background_type == "order":
        if len(orders["background"]) == 1:
            # assume background is present either at the left or right edges of the region.
            background_guess[0][0] = np.min(
                [
                    np.mean(
                        intens.flatten()[chunks[j]].compressed()[
                            tth_ord[:n]
                        ]
                    ),
                    np.mean(
                        intens.flatten()[chunks[j]].compressed()[
                            tth_ord[-n:]
                        ]
                    ),
                ]
            )
        else:  # len(orders['background']) > 1:
            # first value (offset) is mean of left-hand values
            background_guess[0][0] = np.mean(
                intens.flatten()[chunks[j]].compressed()[
                    tth_ord[:n]
                ]
            )
            # if there are more, then calculate a gradient guess.
            background_guess[1][0] = (
                                             np.mean(
                                                 intens.flatten()[chunks[j]].compressed()[
                                                     tth_ord[-n:]
                                                 ]
                                             )
                                             - np.mean(
                                                 intens.flatten()[chunks[j]].compressed()[
                                                    tth_ord[:n]
                                                 ]
                                             )
                                     ) / (
                                             twotheta.flatten()[chunks[j]].compressed()[
                                                 tth_ord[-1]
                                             ]
                                             - twotheta.flatten()[chunks[j]].compressed()[
                                                 tth_ord[0]
                                             ]
                                     )
            # leave all higher order terms as 0 -- assume they are small.
    # FIX ME: do we need a 'flat' option?
    # elif backg_type == 'flat':
    #     backg_guess = backg[0][0]
    return background_guess


def get_chunk_guesses(peeks, orders, twotheta, intens, background_guess,
                      azimu, chunks, conversion_factor, data_as_class, count, n,
                      d_space, bounds, w_guess_fraction, dfour=None, debug=None):
    """

    :param w_guess_fraction:
    :param peeks:
    :param dfour:
    :param orders:
    :param twotheta:
    :param intens:
    :param background_guess:
    :param azimu:
    :param chunks:
    :param conversion_factor:
    :param data_as_class:
    :param count:
    :param n:
    :param d_space:
    :param bounds:
    :param debug:
    :return:
    """
    d_guess = []
    h_guess = []
    w_guess = []
    p_guess = []
    p_fixed = 0
    peaks = []
    lims = []
    j = count

    for k in range(peeks):
        # FIX ME: altered from locals call as now passed as None - check!
        if dfour: # If the positions have been pre-guessed extract values from dfour.
            if debug and k == 0:
                print("dfour in locals")
            if "d-space-type" in orders["peak"][k]:
                coeff_type = orders["peak"][k]["d-space-type"]
            else:
                coeff_type = "fourier"
            t_th_guess = pf.coefficient_expand(
                np.mean(azimu.flatten()[chunks[j]]),
                dfour[k][0],
                coeff_type=coeff_type,
            )
            d_guess = data_as_class.conversion(
                t_th_guess,
                conversion_factor,
                azm=np.mean(azimu.flatten()[chunks[j]]),
            )
            # Finds the index of the closest two-theta to the d-spacing input
            idx = (
                np.abs(twotheta.flatten()[chunks[j]] - t_th_guess)
            ).argmin()
            # FIX ME: The mean of a number of the smallest values would be more stable.

            # Height is intensity of the closest pixel in d-spacing - background at that position
            h_guess = intens.flatten()[chunks[j]][idx]
            if len(background_guess) > 1:
                h_guess = h_guess - (
                        background_guess[0][0]
                        + background_guess[1][0]
                        * (
                                twotheta.flatten()[chunks[j]][idx]
                                - twotheta.flatten()[chunks[j]].min()
                        )
                )
            elif len(background_guess) == 1:
                h_guess = h_guess - background_guess[0][0]

        else:  # 'guess' guesses from data slice - assuming a single peak
            # find the brightest pixels
            # get index of nth brightest pixel.
            idx = np.argsort(
                intens.flatten()[chunks[j]].compressed()
            )[-n:][0]

            # height guess is nth highest intensity - background guess at this position
            h_guess = (intens.flatten()[chunks[j]].compressed())[idx]
            if len(background_guess) > 1:
                h_guess = h_guess - (
                        background_guess[0][0]
                        + background_guess[1][0]
                        * (
                                twotheta.flatten()[chunks[j]].compressed()[
                                    idx
                                ]
                                - twotheta.flatten()[chunks[j]].min()
                        )
                )
            elif len(background_guess) == 1:
                h_guess = h_guess - background_guess[0][0]

            # d-spacing of the highest nth intensity pixel
            d_guess = d_space.flatten()[chunks[j]].compressed()[idx]

        # w_guess is fractional width of the data range
        w_guess = (
                (
                        np.max(twotheta.flatten()[chunks[j]])
                        - np.min(twotheta.flatten()[chunks[j]])
                )
                / w_guess_fraction
                / peeks
        )
        # FIX ME: This is a bit crude. Is there a better way to do it?

        # If profile_fixed exists then set p_guess to profile_fixed values and set p_fixed to 1
        p_guess = 0.5  # Guess half if solving for profile
        p_fixed = (
            0  # Set to 0 so solving unless profile_fixed exists.
        )
        if "profile_fixed" in orders["peak"][k]:
            # Profile fixed is the list of coefficients if the profile is fixed.
            # FIX ME: profile fixed must have the same number of coefficients as required by
            # profile.
            if "symmetry" in orders["peak"][k]:
                symm = orders["peak"][k]["symmetry"]
            else:
                symm = 1
            p_guess = pf.fourier_expand(
                np.mean(azimu.flatten()[chunks[j]]) * symm,
                inp_param=orders["peak"][k]["profile_fixed"],
            )
            p_fixed = 1

        peaks.append(
            {
                "height": [h_guess],
                "d-space": [d_guess],
                "width": [w_guess],
                "profile": [p_guess],
            }
        )

        # DMF added needs checking - required to drive fit and avoid failures
        lims.append(
            pf.parse_bounds(
                bounds,
                d_space.flatten()[chunks[j]],
                intens.flatten()[chunks[j]],
                twotheta.flatten()[chunks[j]],
                param=["height", "d-space", "width", "profile"],
            )
        )

    limits = {
        "background": pf.parse_bounds(
            bounds,
            d_space.flatten()[chunks[j]],
            intens.flatten()[chunks[j]],
            twotheta.flatten()[chunks[j]],
            param=["background"],
        )["background"],
        "peak": lims,
    }
    return peaks, limits, p_fixed




def fit_sub_pattern(
    #two_theta_and_dspacings,
    #azimu,
    #intens,
    data_as_class,
    settings_as_class,
    #orders=None,
    previous_params=None,
    #bounds=None,
    save_fit=None,
    debug=False,
    refine=True,
    iterations=2,
    #f_name=None,
    #fdir=None,
    fit_method=None,
    #num=1
):
    """
    Perform the various fitting stages to the data
    :param fit_method:
    :param fdir:
    :param data_as_class:
    :param two_theta_and_dspacings:
    :param azimu:
    :param intens:
    :param orders:
    :param previous_params:
    :param bounds:
    :param save_fit:
    :param debug:
    :param refine:
    :param iterations:
    :param f_name:
    :return:
    """






# fit_sub_pattern(
#                         # [twotheta_sub, d_spacing_sub, parms_dict],
#                         # azimuth_sub,
#                         # intens_sub,
#                         sub_data,
#                         settings_for_fit, #added
#                         #orders,
#                         params,
#                         # bounds,
#                         save_fit=save_figs,
#                         debug=debug,
#                         refine=refine,
#                         iterations=iterations,
#                         # file_number=j, #added
#                         # subpattern_number=i, #added
#                         # f_name=diff_files[j],
#                         # fdir=fit_settings.Output_directory,
#                         # num=i
#                     )


    # FIXME: these lines are here to get the data class into the file and make it run with minimum reorganisaion. They should be removed eventually when the file is reorganised.
    orders = settings_as_class.subfit_orders
    bounds = settings_as_class.fit_bounds

    f_name    = settings_as_class.subfit_filename
    fdir      = settings_as_class.output_directory
    #directory = settings_as_class.output_directory
    
    twotheta = data_as_class.tth
    d_space  = data_as_class.dspace
    intens   = data_as_class.intensity
    azimu    = data_as_class.azm
    
    # FIXME: this paramter is not needed. It is carried with the data class. 
    #but getting rid of it in the code is time consuming.
    conversion_factor = data_as_class.parameters 
    
    t_th_range = orders["range"]
    
    

    # set a limit to the maximum number of function evaluations.
    # make variable in case need more iterations for other data
    default_max_f_eval = 400

    # Measure the elapsed time during fitting.
    # To help decide what is bad fit or if over fitting the data.
    t_start = time.time()

    # Sort inputs
    # two theta and d-spacing are both 'x' dimensions. Therefore, import them as single array.
    # split them up here if needed.
    
    # if np.shape(two_theta_and_dspacings)[0] == 3:
    #     twotheta = two_theta_and_dspacings[0]
    #     d_space = two_theta_and_dspacings[1]
    #     conversion_factor = two_theta_and_dspacings[2]
    # else:
    #     twotheta = two_theta_and_dspacings
        # FIX ME: need to set a value for 'conversion_factor' in here and propagate its none-use through the code.
        # FIX ME: use d-space in guesses below - need an alternative here or change reliance later

    # if debug:  # get limits for chunk fit plots.
    #     t_th_range = [np.min(twotheta.flatten()), np.max(twotheta.flatten())]

    # set data type for the intensity data.
    # This is needed for saving the fits using save_modelresult/ load_modelresult.
    # load_modelresult fails if the data is a masked integer array.
    if isinstance(intens, int) or isinstance(intens, np.int32):
        intens = float(intens)
    elif isinstance(intens, np.int64):
        intens = np.float64(intens)
    # FIX ME: this doesn't work. the data type is got by print(intens.dtype)

    # FIX ME: need to match orders of arrays to previous numbers.
    # if no parameters
    #     get orders of arrays
    # else if no orders
    #     use parameters as given
    # else
    #     change the size of the parameter arrays to match orders

    # Define the number of peaks to be fitted
    peeks = len(orders["peak"])
    # DMF: FIX ME: orders used before checked?
    if orders:
        # If exists in input file setup peaks as appropriate for fitting
        background_type = "order"
        bg_order = orders["background"]
        orders, backgnd = order_set_peaks(orders, peeks, bg_order)

    if previous_params:
        # check if the previous fit was 'good' i.e. constrains no 'null' values.
        # N.B. null values in json file are read in as None
        clean = io.any_terms_null(previous_params, val_to_find=None)
        if clean == 0:
            # the previous fit has problems so discard it
            print(
                "Propagated fit has problems so discarding it and doing fit from scratch"
            )
            previous_params = None

    if previous_params and orders:
        # If we have both, order takes precedence so update previous_params to match
        previous_params = update_previous_params_from_orders(peeks, previous_params, orders)

    if previous_params and not orders:
        backgnd = previous_params["background"]
        background_type = "coeffs"

    # FIX ME: can we have a situation with no orders or previous_params?

    # check the number of unique azimuths is greater than the number of coefficients.
    check_num_azimuths(peeks, azimu, orders)

    # set up background arrays
    # FIX ME: is this right?
    len_bg = []
    single_len_bg = []
    for val in backgnd:
        len_bg.append(len(val))
        single_len_bg.append(1)
    len_bg = np.array(len_bg)
    single_len_bg = np.array(single_len_bg)

    if not previous_params:
        step = 0
    else:
        step = 5

    # Start fitting loops
    while (step <= 100):
        # 100 is arbitrarily large number and does not (currently) reflect the num. of actual steps/stages in the process.
        # While loops are used so that if the fit is rubbish we can go back and improve it by repeating earlier steps.
        # for chunks step <= 9 and for refine <= 19
        # we are using increments of 10 so that it is possible to utilise different states after the final fit.

        if step <= 9:
            # generate chunks and initial fits
            # or parse previous fits into correct data structure

            # Measure the time taken to do the chunks, the elapsed time during fitting.
            # To help decide which is the best peak parameters.
            chunks_start = time.time()

            # If there is not a previous fit -- Fit data in azimuthal chunks for d.
            if not previous_params:
                if "PeakPositionSelection" not in orders:
                    if peeks > 1:
                        p = "peaks"
                    else:
                        p = "peak"
                    print(
                        "\nNo previous fits or selections. Assuming "
                        + str(peeks)
                        + " "
                        + p
                        + " and group fitting in azimuth...\n"
                    )
                    dfour = None
                    # FIX ME: passed as None now - check functions as expected
                else:
                    print("\nUsing manual selections for initial guesses...\n")
                    # FIX ME: Need to check that the num. of peaks for which we have parameters is the same as the
                    # number of peaks guessed at.
                    dfour = get_manual_guesses(peeks, orders, bounds, twotheta, debug=None)

                # Get chunks according to detector type.
                chunks, azichunks = data_as_class.bins(azimu, orders)

                # Final output list of azimuths with corresponding twotheta_0,h,w

                # setup arrays
                new_d0 = [[] for _ in range(peeks)]
                new_d0_err = [[] for _ in range(peeks)]
                new_h_all = [[] for _ in range(peeks)]
                new_h_all_err = [[] for _ in range(peeks)]
                new_w_all = [[] for _ in range(peeks)]
                new_w_all_err = [[] for _ in range(peeks)]
                new_p_all = [[] for _ in range(peeks)]
                new_p_all_err = [[] for _ in range(peeks)]
                new_bg_all = [[] for _ in range(len(bg_order))]
                new_bg_all_err = [[] for _ in range(len(bg_order))]
                new_azi_chunks = []

                for j in range(len(chunks)):
                    # print('\nFitting to data chunk ' + str(j + 1) + ' of ' + str(len(chunks)) + '\n')
                    if debug:
                        print("\n")
                    print(
                        "\rFitting to data chunk %s of %s "
                        % (str(j + 1), str(len(chunks))),
                        end="\r",
                    )
                    # only compute the fit and save it if there are sufficient data
                    min_dat = 21  # minimum data in each chunk
                    n = 5  # how many data to use in constructing guesses.
                    w_guess_fraction = 10.0
                    # N.B. if min_dat/n is too large the guesses will become very strange.
                    # FIX ME: These should be input variables (for use in extremes)
                    # FIX ME: possibly do this better using scipy.signal.find or scipy.signal.find_peaks_cwt

                    # FIX ME: not meeting this condition means no chunk fitting
                    # is this correct?
                    if np.sum(~intens.flatten()[chunks[j]].mask) >= min_dat:

                        # Background estimates
                        background_guess = get_chunk_background_guess(twotheta, azimu, intens, orders,
                                                                      chunks, background_type, count=j, n=n)

                        # Organise guesses to be refined.
                        peaks, limits, p_fixed = get_chunk_guesses(peeks, orders, twotheta, intens,
                                                                   background_guess,
                                                                   azimu, chunks, conversion_factor,
                                                                   data_as_class, j, n,
                                                                   d_space, bounds, w_guess_fraction,
                                                                   dfour, debug)

                        # Define parameters to pass to fit
                        params = Parameters()
                        # if orders: bg of form [num_orders_0,num_orders_1]
                        # if coeffs: bg of form [[coeff1],[coeff2, coeff3, coeff4]]

                        # Chunks limited to no Fourier expansion so only single value per polynomial order.
                        for b in range(len(background_guess)):
                            if b == 0:
                                params.add(
                                    "bg_c" + str(b) + "_f" + str(0),
                                    background_guess[b][0],
                                    min=limits["background"][0],
                                    max=limits["background"][1],
                                )
                            else:
                                params.add(
                                    "bg_c" + str(b) + "_f" + str(0), background_guess[b][0]
                                )

                        for pk in range(len(peaks)):
                            # FIX ME: do we need locals call?
                            # FIX ME: this looks identical to the else below?
                            if "dfour" in locals():
                                params.add(
                                    "peak_" + str(pk) + "_h0",
                                    peaks[pk]["height"][0],
                                    min=limits["peak"][pk]["height"][0],
                                    max=limits["peak"][pk]["height"][1],
                                )
                                params.add(
                                    "peak_" + str(pk) + "_d0",
                                    peaks[pk]["d-space"][0],
                                    min=limits["peak"][pk]["d-space"][0],
                                    max=limits["peak"][pk]["d-space"][1],
                                )
                                params.add(
                                    "peak_" + str(pk) + "_w0",
                                    peaks[pk]["width"][0],
                                    min=limits["peak"][pk]["width"][0],
                                    max=limits["peak"][pk]["width"][1],
                                )
                                if "profile_fixed" in orders["peak"][pk]:
                                    params.add(
                                        "peak_" + str(pk) + "_p0",
                                        peaks[pk]["profile"][0],
                                        min=limits["peak"][pk]["profile"][0],
                                        max=limits["peak"][pk]["profile"][1],
                                        vary=False,
                                    )
                                else:
                                    params.add(
                                        "peak_" + str(pk) + "_p0",
                                        peaks[pk]["profile"][0],
                                        min=limits["peak"][pk]["profile"][0],
                                        max=limits["peak"][pk]["profile"][1],
                                    )

                            else:
                                # If using bounds/limits should be written as below:
                                params.add(
                                    "peak_" + str(pk) + "_h0",
                                    peaks[pk]["height"][0],
                                    min=limits["peak"][pk]["height"][0],
                                    max=limits["peak"][pk]["height"][1],
                                )
                                params.add(
                                    "peak_" + str(pk) + "_d0",
                                    peaks[pk]["d-space"][0],
                                    min=limits["peak"][pk]["d-space"][0],
                                    max=limits["peak"][pk]["d-space"][1],
                                )
                                params.add(
                                    "peak_" + str(pk) + "_w0",
                                    peaks[pk]["width"][0],
                                    min=limits["peak"][pk]["width"][0],
                                    max=limits["peak"][pk]["width"][1],
                                )
                                if "profile_fixed" in orders["peak"][pk]:
                                    params.add(
                                        "peak_" + str(pk) + "_p0",
                                        peaks[pk]["profile"][0],
                                        min=limits["peak"][pk]["profile"][0],
                                        max=limits["peak"][pk]["profile"][1],
                                        vary=False,
                                    )
                                else:
                                    params.add(
                                        "peak_" + str(pk) + "_p0",
                                        peaks[pk]["profile"][0],
                                        min=limits["peak"][pk]["profile"][0],
                                        max=limits["peak"][pk]["profile"][1],
                                    )

                        if debug:
                            params.pretty_print()
                            print("\n")
                        guess = params

                        # Run actual fit
                        out = pf.fit_model(
                            intens.flatten()[chunks[j]],
                            twotheta.flatten()[chunks[j]],
                            azichunks[j],
                            num_peaks=len(peaks),
                            nterms_back=len(background_guess),
                            conv=conversion_factor,
                            fixed=p_fixed,
                            fit_method=fit_method,
                            weights=None,
                            params=params,
                        )  # , max_n_fev=peeks*default_max_f_eval)
                        params = out.params  # update lmfit parameters
                        if debug:
                            params.pretty_print()
                        for i in range(len(background_guess)):
                            new_bg_all[i].append(params["bg_c" + str(i) + "_f0"].value)
                            new_bg_all_err[i].append(
                                params["bg_c" + str(i) + "_f0"].stderr
                            )
                        for i in range(peeks):
                            new_h_all[i].append(params["peak_" + str(i) + "_h0"].value)
                            new_h_all_err[i].append(
                                params["peak_" + str(i) + "_h0"].stderr
                            )
                            new_d0[i].append(params["peak_" + str(i) + "_d0"].value)
                            new_d0_err[i].append(params["peak_" + str(i) + "_d0"].stderr)
                            new_w_all[i].append(params["peak_" + str(i) + "_w0"].value)
                            new_w_all_err[i].append(
                                params["peak_" + str(i) + "_w0"].stderr
                            )
                            # profile should now be bound by the settings in the parameter
                            new_p_all[i].append(params["peak_" + str(i) + "_p0"].value)
                            new_p_all_err[i].append(
                                params["peak_" + str(i) + "_p0"].stderr
                            )
                            # may have to set profile error when initiate params so have appropriate error value
                            # if fixed

                        new_azi_chunks.append(azichunks[j])

                        # Temp addition: append chunk azimuths and fitted peak intensities to file
                        '''
                        azm_plot = np.tile(azimu.flatten()[chunks[j]][0], 300)
                        azm_plot = ma.array(azm_plot, mask=(~np.isfinite(azm_plot)))
                        gmodel = Model(
                                pf.peaks_model, independent_vars=["two_theta", "azimuth"])
                        tth_range = np.linspace(np.min(
                            twotheta.flatten()[chunks[j]]),
                            np.max(twotheta.flatten()[chunks[j]]),
                            azm_plot.size,
                        )
                        mod_plot = gmodel.eval(
                            params=params,
                            two_theta=tth_range,
                            azimuth=azm_plot,
                            conv=conversion_factor,
                        )
                        # Get the highest peak
                        temp_int = max(mod_plot)
                        temp_peak_ints = [temp_int]
                        temp_x = tth_range[np.where(mod_plot == temp_int)][0]
                        temp_tth = [temp_x]
                        temp_widths = []
                        temp_tth_range = tth_range
                        temp_mod_plot = mod_plot
                        # Find the widest peak width
                        for i in range(peeks):
                            temp_widths.append(params["peak_" + str(i) + "_w0"].value)
                        temp_max_width = max(temp_widths)
                        # print(temp_int, temp_x, temp_max_width)
                        # Use the max width to mask out array to search for next peak
                        for pk in range(len(peaks)-1):
                            upper = temp_x+(2*temp_max_width)
                            lower = temp_x-(2*temp_max_width)
                            # print(upper, lower)
                            # temp_int = max(mod_plot[np.where(temp_tth_range)])
                            temp_mod_plot = temp_mod_plot[np.where((temp_tth_range >= upper) | (temp_tth_range <= lower))]
                            temp_int = max(temp_mod_plot)
                            temp_peak_ints.append(temp_int)
                            # temp_int = max(mod_plot[np.where((temp_tth_range >= upper) | (temp_tth_range <= lower))])
                            temp_x = temp_tth_range[np.where(temp_mod_plot == temp_int)][0]
                            temp_tth.append(temp_x)
                            temp_tth_range = (temp_tth_range[(temp_tth_range >= upper) | (temp_tth_range <= lower)])
                            # print('peak', temp_peak_ints, temp_tth)

                        # Work out which peak is which - assume in 2theta order
                        temp_tth = sorted(temp_tth)
                        # Write to file per sub-pattern per peak
                        for pk in range(len(peaks)):
                            filename = io.make_outfile_name(
                                        f_name,
                                        directory=fdir,
                                        additional_text="ChunksHeights",
                                        orders=orders,
                                        extension=".txt",
                                        overwrite=True,
                                    )
                            # print(azichunks[j], temp_peak_ints[pk])
                            if j==0:
                                with open(filename, "w") as myfile:
                                    myfile.write("%f %f\n"  % (azichunks[j], temp_peak_ints[pk]))
                            else:
                                with open(filename, "a") as myfile:
                                    myfile.write("%f %f\n"  % (azichunks[j], temp_peak_ints[pk]))

                        # print(stop)
                        # debug = True
                        '''

                        # plot the fits.
                        if debug:
                            tth_plot = twotheta.flatten()[chunks[j]]
                            int_plot = intens.flatten()[chunks[j]]
                            azm_plot = np.tile(azimu.flatten()[chunks[j]][0], 300)
                            azm_plot = ma.array(azm_plot, mask=(~np.isfinite(azm_plot)))
                            # FIX ME: required for plotting energy dispersive data.
                            gmodel = Model(
                                pf.peaks_model, independent_vars=["two_theta", "azimuth"]
                            )
                            tth_range = np.linspace(
                                np.min(twotheta.flatten()[chunks[j]]),
                                np.max(twotheta.flatten()[chunks[j]]),
                                azm_plot.size,
                            )
                            mod_plot = gmodel.eval(
                                params=params,
                                two_theta=tth_range,
                                azimuth=azm_plot,
                                conv=conversion_factor,
                            )
                            guess_plot = gmodel.eval(
                                params=guess,
                                two_theta=tth_range,
                                azimuth=azm_plot,
                                conv=conversion_factor,
                            )
                            plt.plot(tth_plot, int_plot, ".", label="data")
                            plt.plot(
                                tth_range,
                                guess_plot,
                                marker="",
                                color="green",
                                linewidth=2,
                                linestyle="dashed",
                                label="guess",
                            )
                            plt.plot(
                                tth_range,
                                mod_plot,
                                marker="",
                                color="red",
                                linewidth=2,
                                label="fit",
                            )
                            #plt.xlim(tth_range)
                            plt.legend()
                            plt.title(io.peak_string(orders) + "; Chunk " + str(j + 1))
                            plt.show()
                
                '''
                return [0, 0] 
                # print(stop)
                '''


                # Fitting stage 2:
                # Feed each d_0,h,w into coefficient function to get fit for coeficient component
                # parameters as output.
                print("\n Performing Series fits...")

                dfour = []
                hfour = []
                wfour = []
                pfour = []
                master_params = Parameters()  # Initiate total Parameter class to add to

                lims = pf.parse_bounds(
                    bounds, d_space.flatten(), intens.flatten(), twotheta.flatten()
                )
                # replace fixed array with dynamic limits

                # Initiate background parameters and perform an initial fit
                comp = "bg"
                coeff_type = pf.params_get_type(orders, comp, peak=j)
                for b in range(len(backgnd)):
                    param_str = "bg_c" + str(b)
                    comp = "f"
                    if b == 0:
                        limits = lims["background"]
                    else:
                        limits = None
                    n_coeff = pf.get_number_coeff(
                        orders, comp, peak=b, azimuths=np.array(new_azi_chunks)
                    )
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,
                        coeff_type=coeff_type,
                        num_coeff=n_coeff,
                        trig_orders=orders["background"][b],
                        limits=limits,
                    )
                    master_params = pf.un_vary_params(
                        master_params, param_str, comp
                    )  # set other parameters to not vary
                    master_params = pf.vary_params(
                        master_params, param_str, comp
                    )  # set these parameters to vary
                    if isinstance(
                        orders["background"][b], list
                    ):  # set part of these parameters to not vary
                        master_params = pf.un_vary_part_params(
                            master_params, param_str, comp, orders["background"][b]
                        )

                    fout = pf.coefficient_fit(
                        azimuth=np.array(new_azi_chunks),
                        ydata=np.array(new_bg_all[b]),
                        inp_param=master_params,
                        param_str=param_str + "_" + comp,
                        symmetry=1,
                        errs=np.array(new_bg_all_err[b]),
                        fit_method="leastsq",
                    )
                    master_params = fout.params

                # initiate peak(s)
                for j in range(peeks):
                    # defines peak string to start parameter name
                    param_str = "peak_" + str(j)

                    # initiate symmetry
                    comp = "s"
                    if "symmetry" in orders["peak"][j].keys():
                        symm = orders["peak"][j]["symmetry"]
                    else:
                        symm = 1
                    master_params = pf.initiate_params(
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
                    arr_names = ["new_h_all", "new_d0", "new_w_all", "new_p_all"]
                    arr_err_names = [
                        "new_h_all_err",
                        "new_d0_err",
                        "new_w_all_err",
                        "new_p_all_err",
                    ]

                    for cp in range(len(comp_list)):
                        comp = comp_list[cp]
                        if comp == "d":
                            symmetry = 1
                        else:
                            symmetry = symm
                        name = comp_names[cp] + "_fixed"
                        if comp_names[cp] + "_fixed" in orders["peak"][j]:
                            fixed = 1
                            vals = orders["peak"][j][
                                comp_names[cp] + "_fixed"
                            ]  # values if fixed.
                        else:
                            fixed = 0
                            vals = np.array(vars()[arr_names[cp]])[j]
                            vals_err = np.array(vars()[arr_err_names[cp]])[j]
                            # FIX ME: this was not checked properly.the values it feeds are not necessarily correct
                            # and the fixed parameters might be fit for.
                        coeff_type = pf.params_get_type(orders, comp, peak=j)
                        n_coeff = pf.get_number_coeff(
                            orders, comp, peak=j, azimuths=np.array(new_azi_chunks)
                        )
                        master_params = pf.initiate_params(
                            master_params,
                            param_str,
                            comp,
                            coeff_type=coeff_type,
                            num_coeff=n_coeff,
                            trig_orders=orders["peak"][j][comp_names[cp]],
                            limits=lims[comp_names[cp]],
                            value=vals,
                        )
                        master_params = pf.un_vary_params(
                            master_params, param_str, comp
                        )  # set other parameters to not vary
                        if fixed == 0:
                            master_params = pf.vary_params(
                                master_params, param_str, comp
                            )  # set these parameters to vary
                            if isinstance(
                                # FIX ME: changed k to j here as think bug but needs checking!
                                orders["peak"][j][comp_names[cp]], list
                            ):  # set part of these parameters to not vary
                                # FIX ME: changed k to j here as think bug but needs checking!
                                master_params = pf.un_vary_part_params(
                                    master_params,
                                    param_str,
                                    comp,
                                    orders["peak"][j][comp_names[cp]],
                                )
                            if n_coeff > len(
                                np.unique(new_azi_chunks)
                            ):  # catch me make sure there are not more coefficients than chunks
                                o = int(np.floor(len(np.unique(new_azi_chunks)) / 2 - 1))
                                un_vary = [x for x in range(0, o)]
                                master_params = pf.un_vary_part_params(
                                    master_params, param_str, comp, un_vary
                                )
                            fout = pf.coefficient_fit(
                                azimuth=np.array(new_azi_chunks),
                                ydata=vals,
                                inp_param=master_params,
                                param_str=param_str + "_" + comp,
                                symmetry=symmetry,
                                errs=vals_err,
                                fit_method="leastsq",
                            )

                        master_params = fout.params
                        # FIX ME. Check which params should be varying and which should not.
                        # Need to incorporate vary and un-vary params as well as partial vary

                if debug:
                    print("Parameters after initial Fourier fits")
                    master_params.pretty_print()

                # plot output of fourier fits....
                if debug:
                    y_lims = np.array([np.min(new_azi_chunks), np.max(new_azi_chunks)])
                    y_lims = np.around(y_lims / 180) * 180
                    azi_plot = range(np.int(y_lims[0]), np.int(y_lims[1]), 2)
                    # azi_plot = range(0, 360, 2)
                    gmodel = Model(pf.coefficient_expand, independent_vars=["azimuth"])

                    fig = plt.figure()
                    ax1 = fig.add_subplot(5, 1, 1)
                    ax1.set_title("D-spacing")
                    for j in range(peeks):
                        param_str = "peak_" + str(j)
                        comp = "d"
                        temp = pf.gather_param_errs_to_list(
                            master_params, param_str, comp
                        )
                        temp_tp = pf.get_series_type(master_params, param_str, comp)
                        if temp_tp == 5:
                            az_plt = np.array(new_azi_chunks)
                        else:
                            az_plt = azi_plot
                        gmod_plot = gmodel.eval(
                                params=master_params, 
                                azimuth=np.array(az_plt), 
                                param=temp[0], 
                                coef_type=temp_tp
                        )
                        ax1.scatter(new_azi_chunks, new_d0[j], s=10)
                        ax1.plot(az_plt, gmod_plot, )

                    # plt.subplot(512)
                    ax2 = fig.add_subplot(5, 1, 2)
                    ax2.set_title("height")
                    for j in range(peeks):
                        if "symmetry" in orders["peak"][j].keys():
                            symm = orders["peak"][j]["symmetry"]
                        else:
                            symm = 1
                        param_str = "peak_" + str(j)
                        comp = "h"
                        temp = pf.gather_param_errs_to_list(
                            master_params, param_str, comp
                        )
                        temp_tp = pf.get_series_type(master_params, param_str, comp)
                        if temp_tp == 5:
                            az_plt = np.array(new_azi_chunks)
                        else:
                            az_plt = azi_plot
                        gmod_plot = gmodel.eval(
                                params=master_params, 
                                azimuth=np.array(az_plt) * symm, 
                                param=temp[0], 
                                coef_type=temp_tp
                        )
                        ax2.scatter(new_azi_chunks, new_h_all[j], s=10)
                        ax2.plot(az_plt, gmod_plot, )

                    ax3 = fig.add_subplot(5, 1, 3)
                    # plt.subplot(513)
                    ax3.set_title("width")
                    for j in range(peeks):
                        if "symmetry" in orders["peak"][j].keys():
                            symm = orders["peak"][j]["symmetry"]
                        else:
                            symm = 1
                        param_str = "peak_" + str(j)
                        comp = "w"
                        temp = pf.gather_param_errs_to_list(
                            master_params, param_str, comp
                        )
                        temp_tp = pf.get_series_type(master_params, param_str, comp)
                        if temp_tp == 5:
                            az_plt = np.array(new_azi_chunks)
                        else:
                            az_plt = azi_plot
                        gmod_plot = gmodel.eval(
                                params=master_params, 
                                azimuth=np.array(az_plt) * symm, 
                                param=temp[0], 
                                coef_type=temp_tp
                        )
                        ax3.scatter(new_azi_chunks, new_w_all[j], s=10)
                        ax3.plot(az_plt, gmod_plot, )

                    ax4 = fig.add_subplot(5, 1, 4)
                    # plt.subplot(514)
                    ax4.set_title("Profile")
                    for j in range(peeks):
                        if "symmetry" in orders["peak"][j].keys():
                            symm = orders["peak"][j]["symmetry"]
                        else:
                            symm = 1
                        param_str = "peak_" + str(j)
                        comp = "p"
                        temp = pf.gather_param_errs_to_list(
                            master_params, param_str, comp
                        )
                        temp_tp = pf.get_series_type(master_params, param_str, comp)
                        if temp_tp == 5:
                            az_plt = np.array(new_azi_chunks)
                        else:
                            az_plt = azi_plot
                        gmod_plot = gmodel.eval(
                                params=master_params, 
                                azimuth=np.array(az_plt) * symm, 
                                param=temp[0], 
                                coef_type=temp_tp
                        )
                        ax4.scatter(new_azi_chunks, new_p_all[j], s=10)
                        ax4.plot(az_plt, gmod_plot, )

                    for k in range(len(len_bg)):
                        x_plt = len(len_bg)
                        y_plt = len(len_bg) * 4 + k + 1
                        ax5 = fig.add_subplot(5, x_plt, y_plt)
                        ax5.set_title("Background")

                        param_str = "bg_c" + str(k)
                        comp = "f"
                        temp = pf.gather_param_errs_to_list(
                            master_params, param_str, comp
                        )
                        temp_tp = pf.get_series_type(master_params, param_str, comp)
                        if temp_tp == 5:
                            az_plt = np.array(new_azi_chunks)
                        else:
                            az_plt = azi_plot
                        gmod_plot = gmodel.eval(
                                params=master_params, 
                                azimuth=np.array(az_plt),
                                param=temp[0], 
                                coef_type=temp_tp
                        )
                        ax5.scatter(new_azi_chunks, new_bg_all[k], s=10)
                        ax5.plot(az_plt, gmod_plot, )

                    fig.suptitle(io.peak_string(orders) + "; Fits to Chunks")

                    if save_fit:
                        filename = io.make_outfile_name(
                            f_name,
                            directory=fdir,
                            additional_text="ChunksFit",
                            orders=orders,
                            extension=".png",
                            overwrite=False,
                        )
                        fig.savefig(filename)
                    plt.show()
                    plt.close()

                # put fits into new_params data structure.
                # new_params = ff.create_new_params(peeks, dfour, hfour, wfour, pfour, bgfour, orders['peak'])

            else:
                # FIX ME: This should load the saved lmfit Parameter class object.
                # Need to initiate usage (save and load).
                # For now have propagated use of new_params so re-instantiate the master_params object below,
                # which is clunky.
                print("Using previously fitted parameters and propagating fit")
                new_params = previous_params

                # FIX ME: should this be some def or sub-function?
                master_params = Parameters()  # initiate total Parameter class to add to

                lims = pf.parse_bounds(
                    bounds, d_space.flatten(), intens.flatten(), twotheta.flatten()
                )

                # Initiate background parameters
                comp = "bg"
                for b in range(len(previous_params["background"])):
                    param_str = "bg_c" + str(b)
                    comp = "f"
                    if b == 0:
                        limits = lims["background"]
                    else:
                        limits = None
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,
                        value=previous_params["background"][b],
                        coeff_type=previous_params["background-type"],
                        limits=limits,
                    )

                #Initiate peaks
                for j in range(peeks):
                    # defines peak string to start parameter name
                    param_str = "peak_" + str(j)

                    # initiate symmetry
                    comp = "s"
                    if "symmetry" in orders["peak"][j].keys():
                        symm = orders["peak"][j]["symmetry"]
                    else:
                        symm = 1
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,
                        0,
                        limits=[10, 1],
                        value=np.array(symm),
                        vary=False,
                    )

                    # initiate h
                    comp = "h"
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,
                        coeff_type=previous_params["peak"][j]["height-type"],
                        value=previous_params["peak"][j]["height"],
                        limits=lims["height"],
                    )
                    # initiate d
                    comp = "d"
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,  # orders['peak'][j]['d-space'],
                        coeff_type=previous_params["peak"][j]["d-space-type"],
                        value=previous_params["peak"][j]["d-space"],
                        limits=lims["d-space"],
                    )
                    # initiate w
                    comp = "w"
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,
                        coeff_type=previous_params["peak"][j]["width-type"],
                        value=previous_params["peak"][j]["width"],
                        limits=lims["width"],
                    )
                    # initiate p
                    comp = "p"
                    if orders and "profile_fixed" in orders["peak"][j]:
                        p_vals = orders["peak"][j]["profile_fixed"]
                        p_fixed = 1
                    else:
                        p_vals = previous_params["peak"][j]["profile"]
                        p_fixed = 0
                    master_params = pf.initiate_params(
                        master_params,
                        param_str,
                        comp,
                        coeff_type=previous_params["peak"][j]["profile-type"],
                        value=p_vals,
                        limits=lims["profile"],
                    )

                # Guess for background, if orders and therefore no guess input, choose appropriate
                # copied from above to allow the code to run. FIX ME: doesn't need all this code.
                single_background = [[] for i in range(len(backgnd))]
                for i in range(len(backgnd)):
                    single_background[i] = [backgnd[i][0]]
                background_guess = single_background
                # bgguess[0][0] = 5

                # def previous_params_to_params(PreviousParams):

                # FIX ME: need to confirm the number of parameters matches the orders of the fits.

            chunks_end = time.time()
            step = step + 10

        if step >= 10:

            # if refine or step>=10 or not PreviousParams:
            if refine or step != 10 or step != 15:
                # Iterate over each parameter series in turn.
                print(
                    "\nRe-fitting for d, h, w, +/-p, bg separately... will refine %i time(s)\n"
                    % iterations
                )
                for j in range(iterations):
                    for k in range(len(backgnd)):
                        param_str = "bg_c" + str(k)
                        comp = "f"
                        # set other parameters to not vary
                        master_params = pf.un_vary_params(master_params, param_str, comp)
                        # set these parameters to vary
                        master_params = pf.vary_params(master_params, param_str, comp)
                        # set part of these parameters to not vary
                        # master_params = ff.un_vary_part_params(master_params, param_str, comp, orders['background'][k])
                        fout = pf.fit_model(
                            intens.flatten(),
                            twotheta.flatten(),
                            azimu.flatten(),
                            num_peaks=peeks,
                            nterms_back=len(backgnd),
                            conv=conversion_factor,
                            fit_method=None,
                            weights=None,
                            params=master_params,
                            max_n_fev=default_max_f_eval,
                        )
                        #if fout.success == 1:
                        master_params = fout.params

                    for k in range(peeks):
                        param_str = "peak_" + str(k)
                        # always need to iterate over d, height and width
                        comp_list = ['h', 'd', 'w']
                        comp_names = ["height", "d-space", "width", "profile"]
                        if "profile_fixed" in orders["peak"][k]:
                            p_fixed = 1
                        else:
                            p_fixed = 0
                            # if the profile is not fixed iterate over this as well.
                            comp_list.append("p")
                        for cp in range(len(comp_list)):
                            comp = comp_list[cp]
                            # print(comp_names[cp])
                            # set other parameters to not vary
                            master_params = pf.un_vary_params(master_params, param_str, comp)
                            # set these parameters to vary
                            master_params = pf.vary_params(master_params, param_str, comp)
                            # set part of these parameters to not vary
                            if isinstance(orders["peak"][k][comp_names[cp]], list):
                                master_params = pf.un_vary_part_params(
                                    master_params,
                                    param_str,
                                    comp,
                                    orders["peak"][k][comp_names[cp]],
                                )
                            if comp == "h":
                                refine_max_f_eval = 5 * peeks * default_max_f_eval
                            else:
                                refine_max_f_eval = default_max_f_eval
                            fout = pf.fit_model(
                                intens.flatten(),
                                twotheta.flatten(),
                                azimu.flatten(),
                                num_peaks=peeks,
                                nterms_back=len(backgnd),
                                conv=conversion_factor,
                                fixed=p_fixed,
                                fit_method=None,
                                weights=None,
                                params=master_params,
                                max_n_fev=refine_max_f_eval,
                            )
                            # print(fout.success, fout.n_fev)
                            master_params = fout.params

                    if debug:
                        print(" ")
                        print(
                            "Parameters after refining series fits "
                            + str(j + 1)
                            + " time(s)"
                        )
                        master_params.pretty_print()

                step = step + 10
            else:
                step = step + 11

        if step >= 20:
            # FIX ME: Need to sort out use of lmfit parameters and new_params, below will not work currently for loaded
            # parameters as master_params not created from new_params yet. If possible use the json
            # import/export within lmfit.

            # Refit through full equation with all data for d,h,w,bg independently
            print("\nFinal fit solving for all parms...\n")

            if any(x == step for x in [20, 21, 22, 25, 26, 27]):
                # if we are on the first go round of the fitting.
                max_n_f_eval = default_max_f_eval
            elif any(x == step for x in [23, 28]):
                max_n_f_eval = 2 * default_max_f_eval
            elif any(x == step for x in [24, 29]):
                max_n_f_eval = np.inf
            else:
                raise ValueError("The value of step here is not possible. Oops.")

            # set all parameters to vary
            for k in range(len(backgnd)):
                param_str = "bg_c" + str(k)
                comp = "f"
                # set these parameters to vary
                master_params = pf.vary_params(master_params, param_str, comp)
                # master_params = ff.un_vary_part_params(master_params, param_str, comp, orders['background'][k])
                # set part of these parameters to not vary
            for k in range(peeks):
                param_str = "peak_" + str(k)
                # always need to iterate over d, height and width
                comp_list = ['h', 'd', 'w']
                comp_names = ["height", "d-space", "width", "profile"]
                if "profile_fixed" in orders["peak"][k]:
                    p_fixed = 1
                else:
                    p_fixed = 0
                    # if the profile is not fixed iterate over this as well.
                    comp_list.append("p")
                for cp in range(len(comp_list)):
                    master_params = pf.vary_params(
                        master_params, param_str, comp_list[cp]
                    )
                    # set part of these parameters to not vary
                    if isinstance(orders["peak"][k][comp_names[cp]], list):
                        master_params = pf.un_vary_part_params(
                            master_params,
                            param_str,
                            comp_list[cp],
                            orders["peak"][k][comp_names[cp]],
                        )

            # master_params.pretty_print()

            # FIX ME: need to sort parameters, bgguess doesn't exist if load previous fit data
            out = pf.fit_model(
                intens.flatten(),
                twotheta.flatten(),
                azimu.flatten(),
                num_peaks=peeks,
                nterms_back=len(background_guess),
                conv=conversion_factor,
                fixed=p_fixed,
                fit_method=None,
                weights=None,
                params=master_params,
                max_n_fev=max_n_f_eval,
            )
            master_params = out.params

            if out.success == 1:
                # it worked, carry on
                step = step + 100
                master_params = out.params
            elif step == 24 and out.success == 0:
                raise ValueError(
                    "Oh Dear. It should not be possible to get here. Something has gone very wrong with the fitting."
                )
            elif step == 29 and out.success == 0:
                step = 0  # go back to the start, discard PreviousParams and do the chunks for this data set.
            elif any(x == step for x in [20, 25]) and out.success == 0:
                step = step - 7
                iterations = np.max((iterations, 3))
            elif any(x == step for x in [21, 22, 23, 26, 27, 28]) and out.success == 0:
                step = step - 9
                if any(x == step for x in [23, 28]):
                    iterations = np.max((iterations, 3))
            else:
                raise ValueError("The value of step here is not possible. Oops.")

            # FIX ME: we could make an option for output the chunks without any Fourier/global fitting.
            # Would this be useful?

    print("\nFinal Coefficients\n")
    # print(out.fit_report(min_correl=1))
    print(out.fit_report(show_correl=False))

    # Print some stats
    print("number data", out.ndata)
    print("degrees of freedom", out.nfree)
    print("ChiSquared", out.chisqr)

    # Write master_params to new_params dict object
    new_params = pf.params_to_new_params(
        master_params,
        num_peaks=peeks,
        num_bg=len(background_guess),
        order_peak=orders["peak"],
    )

    # get all correlation coefficients
    # FIX ME: this lists all the coefficients rather than just the unique half -- ie.
    # it contains corr(a,b) and corr(b,a)
    # FIX ME: should probably be a function
    # N.B. converted to a string so that it does not take up so many times in output file.
    # Also why new lines are removed.
    correl = {}
    for key in master_params.keys():
        correl[key] = master_params[key].correl
    correl_str = json.dumps(correl)
    correl_str = correl_str.replace("\n", "")
    new_params.update({"correlation_coeffs": correl_str})

    tth_min = twotheta.min()
    tth_max = twotheta.max()
    d_min = d_space.min()
    d_max = d_space.max()
    extent = [[tth_min, tth_max], [d_min, d_max]]
    # FIX ME: This is taking the maximum and the minimum of the data not the 'range' itself.
    new_params.update({"range": extent})
    if "note" in orders:
        new_params.update({"note": orders["note"]})

    # Elapsed time for fitting
    t_end = time.time()
    t_elapsed = t_end - t_start
    chunks_time = chunks_end - chunks_start
    # get the rest of the fit stats from the lmfit output.
    fit_stats = {
        "time-elapsed": t_elapsed,
        "chunks-time": chunks_time,
        "status": step,
        "sum-residuals-squared": np.sum(out.residual ** 2),
        "function-evaluations": out.nfev,
        "n-variables": out.nvarys,
        "n-data": out.ndata,
        "degree-of-freedom": out.nfree,
        "ChiSq": out.chisqr,
        "RedChiSq": out.redchi,
        "aic": out.aic,
        "bic": out.bic,
    }
    new_params.update({"FitProperties": fit_stats})

    # add peak names to new_params
    new_params.update({"PeakLabel": io.peak_string(orders)})

    # Plot results to check
    view = 1
    if save_fit == 1 or view == 1 or debug:
        print("\nPlotting results for fit...\n")

        y_lims = np.array([np.min(azimu.flatten()), np.max(azimu.flatten())])
        y_lims = np.around(y_lims / 180) * 180

        gmodel = Model(pf.peaks_model, independent_vars=["two_theta", "azimuth"])

        full_fit_intens = gmodel.eval(
            params=master_params,
            two_theta=twotheta.flatten(),
            azimuth=azimu.flatten(),
            conv=conversion_factor,
        )

        if conversion_factor["DispersionType"] == "EnergyDispersive":
            azi_plot = np.unique(azimu.flatten())
        elif conversion_factor["DispersionType"] == "AngleDispersive":
            azi_plot = np.array(list(range(np.int(y_lims[0]), np.int(y_lims[1]), 2)))
        # centroid = []
        # for i in range(peeks):
        #     param = new_params["peak"][i]["d-space"]
        #     param_type = new_params["peak"][i]["d-space-type"]
        #     centroid.append(
        #         pf.centroid_conversion(
        #             conversion_factor,
        #             pf.coefficient_expand(azi_plot, param=param, coeff_type=param_type),
        #             azi_plot,
        #         )
        #     )


        fit_centroid = []
        for i in range(peeks):
            param = new_params["peak"][i]["d-space"]
            param_type = new_params["peak"][i]["d-space-type"]
            fit_centroid.append(
                pf.centroid_conversion(
                    conversion_factor,
                    pf.coefficient_expand(azi_plot, param=param, coeff_type=param_type),
                    azi_plot,
                )
            )
        
        #plot the modelled data
        fig = plt.figure()
        data_as_class.plot_fitted(fig_plot=fig, model=full_fit_intens, fit_centroid=[azi_plot, fit_centroid] )
        title_str = io.peak_string(orders) + "; final fit"
        if "note" in orders:
            title_str = title_str + " " + orders["note"]
        plt.suptitle(title_str)
        plt.tight_layout()

        # # full_fit_intens = inp
        # # Set the maximum to the n+1th value in the array
        # # FIX ME: the max find is not necessarily efficient. A quicker way should be found if it exists.
        # max_pos = 0

        # val_max1 = -np.sort(-intens.flatten())[max_pos]  # np.max(intens.flatten())
        # val_max2 = -np.sort(-full_fit_intens.flatten())[
        #     max_pos
        # ]  # np.max(full_fit_intens.flatten())
        # val_max3 = -np.sort(-(intens.flatten() - full_fit_intens.flatten()))[max_pos]
        # # np.max((intens.flatten()-full_fit_intens.flatten()))
        # val_max = np.max([val_max1, val_max2])
        # val_min1 = np.min(intens.flatten())
        # val_min2 = np.min(full_fit_intens.flatten())
        # val_min3 = np.min((intens.flatten() - full_fit_intens.flatten()))
        # val_min = np.min([val_min1, val_min2])

        # if conversion_factor["DispersionType"] == "EnergyDispersive":
        #     xlabel_str = "Energy (keV)"
        # elif conversion_factor["DispersionType"] == "AngleDispersive":
        #     xlabel_str = r'2$\theta$ (deg)'#"Two theta (deg)"
        # else:
        #     xlabel_str = "Dispersive units"

        # dot_size = 6#0.5

        # fig = plt.figure()
        # axf = fig.add_subplot(1, 3, 1)
        # ax_o1 = plt.subplot(131)
        # # plt.scatter(twotheta, azimuth, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
        # plt.scatter(
        #     twotheta.flatten(),
        #     azimu.flatten(),
        #     s=dot_size,
        #     c=(intens.flatten()),
        #     edgecolors="none",
        #     cmap=plt.cm.magma_r,
        #     vmin=val_min,
        #     vmax=val_max,
        # )
        # ax_o1.set_title("Data")
        # ax_o1.set_xlabel(xlabel_str)
        # ax_o1.set_ylabel("Azimuth (deg)")
        # ax_o1.set_xlim([np.min(twotheta.flatten()), np.max(twotheta.flatten())])
        # ax_o1.set_ylim(y_lims)
        # locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)
        # plt.colorbar()
        # ax_o2 = plt.subplot(132)
        # # plt.scatter(twotheta, azimuth, s=1, c=np.log(full_fit_intens), edgecolors='none', cmap=plt.cm.jet)
        # plt.scatter(
        #     twotheta.flatten(),
        #     azimu.flatten(),
        #     s=dot_size,
        #     c=(full_fit_intens.flatten()),
        #     edgecolors="none",
        #     cmap=plt.cm.magma_r,
        #     vmin=val_min,
        #     vmax=val_max,
        # )
        # for i in range(peeks):
        #     plt.plot(centroid[i], azi_plot, "k--", linewidth=0.5)
        # plt.colorbar()
        # ax_o2.set_title("Model")
        # ax_o2.set_xlabel(xlabel_str)
        # ax_o2.set_ylabel("Azimuth (deg)")
        # ax_o2.set_xlim([np.min(twotheta.flatten()), np.max(twotheta.flatten())])
        # ax_o2.set_ylim(y_lims)
        # locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)
        # ax_o3 = plt.subplot(133)
        # resid_colour = residuals_colour_scheme(val_max3, val_min3)
        # print("residuals range", val_max3, val_min3)
        # # plt.scatter(twotheta, azimuth, s=1, c=np.log(intens-full_fit_intens), edgecolors='none', cmap=plt.cm.jet)

        # plt.scatter(
        #     twotheta.flatten(),
        #     azimu.flatten(),
        #     s=dot_size,
        #     c=(intens.flatten() - full_fit_intens.flatten()),
        #     edgecolors="none",
        #     cmap=resid_colour,
        #     vmin=val_min3,
        #     vmax=val_max3,
        # )
        # for i in range(peeks):
        #     plt.plot(centroid[i], azi_plot, "k--", linewidth=0.5)
        # ax_o3.set_title("Residuals (data - model)")
        # ax_o3.set_xlabel(xlabel_str)
        # ax_o3.set_ylabel("Azimuth (deg)")
        # ax_o3.set_xlim([np.min(twotheta.flatten()), np.max(twotheta.flatten())])
        # ax_o3.set_ylim(y_lims)
        # locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)
        # plt.colorbar()
        # if "note" in orders:
        #     title_str = io.peak_string(orders) + "; final fit; " + orders["note"]
        # else:
        #     title_str = io.peak_string(orders) + "; final fit"
        # plt.suptitle(title_str)
        # plt.tight_layout()

        # #        # make a second figure with the Q-Q plot of the data and fit.
        # #        #the second subplot is a cdf of the residuals with a gaussian added to it.
        # #        fig2 = plt.figure()
        # #        axQ = plt.subplots()
        # #        plt.scatter(intens.flatten(), full_fit_intens.flatten(), s=dot_size)
        # #        #axQ.set_xlabel('Observations')
        # #        #axQ.set_ylabel('Model')
        # #        plt.grid(True)
        # #        plt.tight_layout()

        # save figures without overwriting old names
        if save_fit:
            filename = io.make_outfile_name(
                f_name, directory=fdir, orders=orders, extension=".png", overwrite=False
            )
            plt.savefig(filename)

        if view == 1 or debug:
            plt.show()

        plt.close()
        print("Done with Figure")

    # Save lmfit structure
    if save_fit:
        filename = io.make_outfile_name(
            f_name, directory=fdir, orders=orders, extension=".sav", overwrite=True
        )

        try:
            save_modelresult(out, filename)
        except BaseException:
            print("Cannot save lmfit object, trying again")
            # if os.path.isfile(filename):
            #     os.remove(filename)
            #     save_modelresult(out, filename)
            # else:
            #     print("File does not exist!")

    return [new_params, out]
