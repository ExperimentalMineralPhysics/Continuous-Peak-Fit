#!/usr/bin/env python

__all__ = ["fit_sub_pattern"]

# CPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

import json
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from lmfit import Model
from lmfit.model import save_modelresult  # , load_modelresult

import cpf.lmfit_model as lmm
import cpf.peak_functions as pf
import cpf.series_constraints as sc
import cpf.series_functions as sf
from cpf.fitsubpattern_chunks import fit_chunks, fit_series
from cpf.util.io import (
    has_huge_errors,
    has_value,
    make_outfile_name,
    numpy_to_json,
    peak_string,
)
from cpf.util.logging import get_logger

logger = get_logger("cpf.XRD_FitSubpattern")


np.set_printoptions(threshold=sys.maxsize)


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
                not sf.get_order_from_params(orders["peak"][y]["profile_fixed"])
                == orders["peak"][y]["profile"]
            ):
                logger.warning(
                    " ".join(
                        map(
                            str,
                            [
                                (
                                    "Peak %i: The order of the profile does not match that of the fixed profile. Changing profile to match fixed profile."
                                    % y
                                )
                            ],
                        )
                    )
                )
                orders["peak"][y]["profile"] = sf.get_order_from_params(
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
            coeff_type = sf.get_params_type(previous_params, param, peak=y)
            coeff_type = sf.coefficient_type_as_number(coeff_type)
            if (
                coeff_type != sf.coefficient_types()["independent"]
            ):  # if parameters are not independent
                if sc.BiggestValue(orders["peak"][y][param]) > sf.get_order_from_params(
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
                elif sc.BiggestValue(
                    orders["peak"][y][param]
                ) < sf.get_order_from_params(previous_params["peak"][y][param]):
                    change_by = (
                        np.size(previous_params["peak"][y][param])
                        - np.max(orders["peak"][y][param]) * 2
                        + 1
                    )
                    previous_params["peak"][y][param] = previous_params["peak"][y][
                        param
                    ][1 : -(change_by - 1)]
        if "profile_fixed" in orders["peak"][y]:
            previous_params["peak"][y]["profile"] = orders["peak"][y]["profile_fixed"]
        elif "profile_fixed" in previous_params["peak"][y][param]:
            del previous_params["peak"][y]["profile_fixed"]

    # loop for background orders/size
    if (
        sf.coefficient_type_as_number(sf.get_params_type(previous_params, "background"))
        != sf.coefficient_types()["independent"]
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
                if np.max(orders["background"][y]) > sf.get_order_from_params(
                    previous_params["background"][y]
                ):
                    change_by = (
                        np.max(orders["background"][y]) * 2
                        + 1
                        - np.size(previous_params["background"][y])
                    )
                    previous_params["background"][y] = (
                        previous_params["background"][y] + [0] * change_by
                    )
                elif np.max(orders["background"][y]) < sf.get_order_from_params(
                    previous_params["background"][y]
                ):
                    change_by = np.size(previous_params["background"][y]) - (
                        np.max(orders["background"][y]) * 2 + 1
                    )
                    previous_params["background"][y] = previous_params["background"][y][
                        :-change_by
                    ]
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
    max_coeff = 0
    total_coeff = 0
    choice_list = ["d-space", "height", "width", "profile"]
    for y in range(peeks):
        # loop over parameters
        for param in choice_list:
            # coeff_type = sf.coefficient_type_as_number(
            #     sf.get_params_type(orders, param, peak=y)
            # )
            total_coeff += sf.get_number_coeff(orders, param, peak=y, azimuths=azimu)
            max_coeff = np.max(
                [max_coeff, sf.get_number_coeff(orders, param, peak=y, azimuths=azimu)]
            )
            # if coeff_type != sf.coefficient_types(full=True)["independent"]["num"]:
            #     # if parameters are not independent
            #     total_coeff += sf.get_number_coeff(orders, param, peak=y)
            #     max_coeff = np.max(
            #         [max_coeff, sf.get_number_coeff(orders, param, peak=y)]
            #     )                
    param = "background"
    for y in range(np.max([len(orders["background"])])):
        # coeff_type = sf.coefficient_type_as_number(
        #     sf.get_params_type(orders, param, peak=y)
        # )
        # if parameters are not independent
        total_coeff += sf.get_number_coeff(orders, param, peak=y, azimuths=azimu)
        max_coeff = np.max([max_coeff, sf.get_number_coeff(orders, "background", azimuths=azimu)])
        # if coeff_type != sf.coefficient_types(full=True)["independent"]["num"]:
        #     # if parameters are not independent
        #     total_coeff += sf.get_number_coeff(orders, param, peak=y)
        #     max_coeff = np.max([max_coeff, sf.get_number_coeff(orders, "background")])
    safe_to_fit = True
    if max_coeff > len(np.unique(azimu)):
        err_str = (
            "The maximum order, %i, needs more coefficients than the number of unique azimuths, %i. "
            "It is not possible to fit."
            % (sf.get_order_from_coeff(max_coeff), len(np.unique(azimu)))
        )
        logger.error(err_str)
        safe_to_fit = False
        # raise ValueError(err_str)
    if total_coeff > len(np.unique(azimu)):
        err_str = (
            "The number of coefficients, %i, is greater than the number of unique azimuths, %i. "
            "It is not possible to fit."
            % (total_coeff, len(np.unique(azimu)))
        )
        logger.error(err_str)
        safe_to_fit = False
        # raise ValueError(err_str)
    
    return safe_to_fit


def fit_sub_pattern(
    data_as_class,
    settings_class,
    previous_params=None,
    mode="fit",
    save_fit = False,
    **kwargs
):
    """
    Fits the selected data to the model specified in the settings. 
    
    Parameters
    ----------
    data_as_class : cpf data_class type.
        Data class as defined by cpf.input_types.
    settings_class : Settings()
        cpf Settings() class, which contains all parameters and settings to perform the fits.
    previous_params : list, optional
        list of dictionaries, containing the fits from the previous iteration. The default is None.
    mode : str
        How to process the provided data, changes path through the method.
        Options are:
        - cascade
        - fit -- fits the data.             
        The default is "fit".
    **kwargs : key, value pairs
        key, value arguments arguments. Passed though to called methods. Not used here.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # parse values that affect the fitting from settings (if they are set). 
    if "fit_options" not in settings_class.__dict__:
        settings_class.fit_options = {}
    refine = settings_class.fit_options.get("refine", True)
    iterations = settings_class.fit_options.get("iterations", 3)
    fit_method = settings_class.fit_options.get("fit_method", None)
    histogram_type = settings_class.fit_options.get("histogram_type", None)
    histogram_bins = settings_class.fit_options.get("histogram_bins", None)
    min_data_intensity = settings_class.__dict__.get("fit_min_data_intensity", 1) # ***
    min_data_intensity = settings_class.fit_options.get("min_data_intensity", min_data_intensity)
    min_peak_intensity = settings_class.__dict__.get("fit_min_peak_intensity", "std") # ***
    min_peak_intensity = settings_class.fit_options.get("min_peak_intensity", min_peak_intensity)
    large_errors = settings_class.fit_options.get("large_errors", 300)
    # save_fit = settings_class.fit_options.get("save_fit", False)
    # FIXME: rows with *** are historical and can be removed when settings are updated/standardised. 
    
    # set a limit to the maximum number of function evaluations.
    # make variable in case need more iterations for other data
    default_max_f_eval = settings_class.fit_options.get("default_max_f_eval", 400)
    max_f_eval = settings_class.fit_options.get("max_f_eval", default_max_f_eval)

    # Measure the elapsed time during fitting.
    # To help decide what is bad fit or if over fitting the data.
    t_start = time.time()

    # set data type for the intensity data.
    # This is needed for saving the fits using save_modelresult/ load_modelresult.
    # load_modelresult fails if the data is a masked integer array.
    if isinstance(data_as_class.intensity, int) or isinstance(
        data_as_class.intensity, np.int32
    ):
        data_as_class.intensity = float(data_as_class.intensity)
    elif isinstance(data_as_class.intensity, np.int64):
        data_as_class.intensity = np.float64(data_as_class.intensity)
    # FIX ME: this doesn't seem to work. the data type is got by logger.info(" ".join(map(str, [(intens.dtype)])))

    # FIX ME: need to match orders of arrays to previous numbers.
    # if no parameters
    #     get orders of arrays
    # else if no orders
    #     use parameters as given
    # else
    #     change the size of the parameter arrays to match orders

    # Define the number of peaks to be fitted
    peeks = len(settings_class.subfit_orders["peak"])

    # DMF: FIX ME: orders used before checked?
    # if settings_class.subfit_orders:
    # If exists in input file setup peaks as appropriate for fitting
    # background_type = "order"
    # bg_order = orders["background"]
    # orders, backgnd = order_set_peaks(settings_class.subfit_orders, peeks, settings_class.subfit_orders["background"])

    if previous_params:
        # check if the previous fit was 'good' i.e. constrains no 'null' values.
        # N.B. null values in json file are read in as None
        any_bad_vals = has_value(previous_params, val=None)
        any_bad_vals = has_huge_errors(
            previous_params, min_ratio=large_errors, is_huge=any_bad_vals
        )
        if any_bad_vals == True:
            # the previous fit has problems so discard it
            logger.moreinfo(
                "Propagated fit has problems so discarding it and doing fit from scratch"
            )
            previous_params = None

    if previous_params:
        # check if the previous fit d-spacings fall within the bounds.
        # if they are not then assume the previous fit must have been poor and discard it.
        fit_centroid = []
        azims = data_as_class.test_azims()
        for i in range(peeks):
            param = previous_params["peak"][i]["d-space"]
            param_type = previous_params["peak"][i]["d-space_type"]
            fit_centroid.append(
                data_as_class.conversion(
                    sf.coefficient_expand(azims, param=param, coeff_type=param_type),
                    azims,
                    reverse=1,
                )
            )
        if (
            np.min(fit_centroid) < settings_class.subfit_orders["range"][0]
            or np.min(fit_centroid) > settings_class.subfit_orders["range"][1]
        ):
            logger.moreinfo(
                "Fitted d-spacing limits are out of bounds; discarding the fit and starting again."
            )
            previous_params = None

    # initiate counter for while loop.
    if previous_params is None or previous_params == [] or previous_params is False:
        previous_params = None  # make sure expected value for later
        step = [0]
    else:
        step = [5]

    if previous_params and settings_class.subfit_orders:
        # If we have both, order takes precedence so update previous_params to match
        previous_params = update_previous_params_from_orders(
            peeks, previous_params, settings_class.subfit_orders
        )

    #initiate parameter set for fitting
    master_params = lmm.initiate_all_params_for_fit(
        settings_class,
        data_as_class,
        values=previous_params,
        # debug=debug,
    )
    
    if data_as_class.intensity.size == 0:
        # then we have no data in the range. 
        # issue warning and skip fitting. 
        
        #initiate parameter set for fitting
        master_params = lmm.initiate_all_params_for_fit(
            settings_class,
            data_as_class,
            values=previous_params,
            # debug=debug,
        )
    
        # set step to -5001 so that it is still negative at the end
        step=[-5001]  # get to the end and void the fit
        # void so send empty parameter set to out.
        prms = master_params.valuesdict()
        for i in prms.keys():
            if "_s" in i:
                pass
            else:
                 prms[i] = 0
        master_params.set(**prms)
    
    elif np.max(data_as_class.intensity) <= min_data_intensity:
        # check if the data intensity is above threshold.
        # then there is likely no determinable peak in the data
        logger.moreinfo(
            f"Not sufficient intensity in the data to proceed with fitting (I_max < {min_data_intensity})."
        )
        # set step to -21 so that it is still negative at the end
        step.append(-21)  # get to the end and void the fit
        # void so send empty parameter set to out.
        prms = master_params.valuesdict()
        for i in prms.keys():
            if "_s" in i:
                pass
            else:
                 prms[i] = 0
        master_params.set(**prms)
                     
    elif not check_num_azimuths(peeks, data_as_class.azm, settings_class.subfit_orders):
        # check the number of unique azimuths is greater than the number of coefficients.
        # logger messages added in function -- not needed here
        # set step to -31 so that it is still negative at the end
        step.append(-31)  
        # get to the end and void the fit
        prms = master_params.valuesdict()
        for i in prms.keys():
            if "_s" in i:
                pass
            else:
                 prms[i] = 0
        master_params.set(**prms)
        
    if step[-1] < 0:
        # voided but still need some numbers for the outputs
        chunks_start = time.time()
        chunks_end = time.time()   

    # Start fitting loops
    while step[-1] >= 0 and step[-1] <= 100:
        # 100 is arbitrarily large number and does not reflect the num. of actual steps/stages in the process.
        # While loops are used so that if the fit is rubbish we can go back and improve it by repeating earlier steps.
        # for chunks step <= 9 and for refine <= 19
        # we are using increments of 10 so that it is possible to record different states or routes after the final fit.

        if step[-1] <= 9:
            # generate chunks and initial fits
            # or parse previous fits into correct data structure

            # Measure the time taken to do the chunks, the elapsed time during fitting.
            # To help decide which is the best peak parameters.
            chunks_start = time.time()

            if step[-1] >= 0 and not previous_params:
                # There is no previous fit -- Fit data in azimuthal chunks
                # using manual guesses ("PeakPositionSelection") if they exist.
                chunk_fits, chunk_positions = fit_chunks(
                    data_as_class,
                    settings_class,
                    mode=mode,
                    histogram_type=histogram_type,
                    histogram_bins=histogram_bins,
                    # debug=debug,
                    fit_method=fit_method,
                )
                if mode.lower() == "cascade":
                    # some cascade option. so exit returning values.
                    return chunk_fits, chunk_positions

                # Fitting stage 2:
                logger.moreinfo(" ".join(map(str, [("Performing Series fits...")])))
                # Feed each d_0,h,w into coefficient function to get fit for coeficient component
                # get lmfit parameters as output.

                # fit the chunk values with fourier/spline series.
                # iterate over each parameter in turn
                master_params = fit_series(
                    master_params,
                    (chunk_fits, chunk_positions),
                    settings_class,
                    start_end=[data_as_class.azm_start, data_as_class.azm_end],
                    # debug=debug,
                    save_fit=save_fit,
                )

                if save_fit or logger.is_below_level("MOREINFO"):
                    # Write master_params to new_params dict object
                    new_params = lmm.params_to_new_params(
                        master_params, orders=settings_class.subfit_orders
                    )
                    # write to a file.
                    filename = make_outfile_name(
                        settings_class.subfit_filename,
                        directory=settings_class.output_directory,
                        orders=settings_class.subfit_orders,
                        additional_text="InitialSeriesFit",
                        extension=".json",
                        overwrite=True,
                    )
                    with open(filename, "w") as TempFile:
                        # Write a JSON string into the file.
                        json.dump(
                            new_params,
                            TempFile,
                            sort_keys=True,
                            indent=2,
                            default=numpy_to_json,
                        )

                # check if peak intensity is above threshold
                ave_intensity = []
                for k in range(peeks):
                    ave_intensity.append(
                        sf.get_series_mean(master_params, "peak_" + str(k), comp="h")
                    )
                if isinstance(min_peak_intensity, str):
                    # is the threshold a multiple of the data's variance of standard deviation?
                    if min_peak_intensity.find("var") != -1:
                        var = np.var(data_as_class.intensity.flatten())
                        min_peak_intensity = min_peak_intensity.strip("var").strip("*")
                        if min_peak_intensity == "":
                            min_peak_intensity = 1
                        min_peak_intensity = float(min_peak_intensity) * var

                    elif min_peak_intensity.find("std") != -1:
                        std = np.std(data_as_class.intensity.flatten())
                        min_peak_intensity = min_peak_intensity.strip("std").strip("*")
                        if min_peak_intensity == "":
                            min_peak_intensity = 1
                        min_peak_intensity = float(min_peak_intensity) * std

                    else:
                        err_str = "The 'min_peak_intensity' has to be a number or a multiple of 'var' or 'std'"
                        logger.critical(" ".join(map(str, [(err_str)])))
                        raise ValueError(err_str)

                if np.max(ave_intensity) <= min_peak_intensity:
                    # then there is no determinable peak(s) in the data
                    logger.moreinfo(
                        f"Not sufficient intensity in the chunked peaks to proceed with fitting (h_max < {min_peak_intensity})."
                    )
                    # set step to -11 so that it is still negative at the end
                    step.append(-11)  # get to the end and void the fit
                    # void so send empty parameter set to out.
                    fout = lmm.initiate_all_params_for_fit(
                        settings_class,
                        data_as_class,
                        # debug=debug,
                    )
                    prms = fout.valuesdict()
                    for i in prms.keys():
                        if "_s" in i:
                            pass
                        else:
                             prms[i] = 0
                    fout.set(**prms)

            elif step[-1] >= 0 and previous_params:
                logger.moreinfo(
                    "Using previously fitted parameters and propagating fit"
                )
                # FIX ME: This should load the saved lmfit Parameter class object.
                # Need to initiate usage (save and load).
                # For now have propagated use of new_params so re-instantiate the master_params object below,
                # which is clunky.

                # FIX ME: need to confirm the number of parameters matches the orders of the fits.

            step.append(step[-1] + 10)
            chunks_end = time.time()

        if step[-1] >= 10:
            # if refine or step>=10 or not PreviousParams:
            if (refine or step[-1] != 10 or step[-1] != 15) and iterations >= 1:
                # Iterate over each parameter series in turn.
                logger.moreinfo(
                    " ".join(
                        map(
                            str,
                            [
                                (
                                    "Re-fitting for d, h, w, +/-p, bg separately... will refine %i time(s)"
                                    % iterations
                                )
                            ],
                        )
                    )
                )
                for j in range(iterations):
                    for k in range(len(settings_class.subfit_orders["background"])):
                        param_str = "bg_c" + str(k)
                        comp = "f"
                        # set other parameters to not vary
                        master_params = lmm.un_vary_params(
                            master_params, param_str, comp
                        )
                        # set these parameters to vary
                        master_params = lmm.vary_params(master_params, param_str, comp)
                        # set part of these parameters to not vary
                        master_params = lmm.un_vary_part_params(
                            master_params,
                            param_str,
                            comp,
                            settings_class.subfit_orders["background"][k],
                        )
                        fout = lmm.fit_model(
                            data_as_class,
                            settings_class.subfit_orders,
                            master_params,
                            start_end=[data_as_class.azm_start, data_as_class.azm_end],
                            fit_method=fit_method,
                            weights=None,
                            max_n_fev=default_max_f_eval,
                        )
                        master_params = fout.params

                    # iterate over the peak components in order of size (height)
                    # get mean height of peak from parameters
                    peak_mean = []
                    for k in range(peeks):
                        peak_mean.append(
                            sf.get_series_mean(
                                master_params, "peak_" + str(k), comp="h"
                            )
                        )
                    # get order to process peaks in -- reverse order starting with the most intense.
                    peak_order = np.array(peak_mean).argsort()[::-1]

                    comp_list, comp_names = pf.peak_components(include_profile=True)

                    # iterate over peak parameters in order and then by peak.
                    # testing on triples (NaMgF3) shows this is better than iterating
                    # over peaks as the higher level loop.
                    for cp in range(len(comp_list)):
                        comp = comp_list[cp]
                        for l in range(peeks):
                            k = peak_order[l]
                            param_str = "peak_" + str(k)

                            if not (
                                comp_names[cp] + "_fixed"
                                in settings_class.subfit_orders["peak"][k]
                            ):
                                # set other parameters to not vary
                                master_params = lmm.un_vary_params(
                                    master_params, param_str, comp
                                )
                                # set these parameters to vary
                                master_params = lmm.vary_params(
                                    master_params, param_str, comp
                                )
                                # set part of these parameters to not vary
                                if isinstance(
                                    sc.SeriesValues(
                                        settings_class.subfit_orders["peak"][k][
                                            comp_names[cp]
                                        ]
                                    ),
                                    list,
                                ):  # set part of these parameters to not vary
                                    # isinstance(
                                    #     settings_class.subfit_orders["peak"][k][
                                    #         comp_names[cp]
                                    #     ],
                                    #     list,
                                    # ):
                                    master_params = lmm.un_vary_part_params(
                                        master_params,
                                        param_str,
                                        comp,
                                        sc.SeriesValues(
                                            settings_class.subfit_orders["peak"][k][
                                                comp_names[cp]
                                            ]
                                        ),
                                    )
                                if comp == "h":
                                    refine_max_f_eval = 5 * peeks * default_max_f_eval
                                else:
                                    refine_max_f_eval = default_max_f_eval
                                fout = lmm.fit_model(
                                    data_as_class,
                                    settings_class.subfit_orders,
                                    master_params,
                                    start_end=[
                                        data_as_class.azm_start,
                                        data_as_class.azm_end,
                                    ],
                                    fit_method=fit_method,
                                    weights=None,
                                    max_n_fev=refine_max_f_eval,
                                )
                                master_params = fout.params
                                logger.log_lmfit_obj(master_params, level="EFFUSIVE")

                    logger.effusive(
                        " ".join(
                            map(
                                str,
                                [
                                    (
                                        "Parameters after refining series fits %i time(s):"
                                        % (j + 1)
                                    )
                                ],
                            )
                        )
                    )
                    # master_params.pretty_print()
                    logger.log_lmfit_obj(master_params, level="EFFUSIVE")

                step.append(step[-1] + 10)

                # get mean height of chunked peaks and check if it is greater than threshold
                ave_intensity = []
                for k in range(peeks):
                    ave_intensity.append(
                        sf.get_series_mean(master_params, "peak_" + str(k), comp="h")
                    )
                if isinstance(min_peak_intensity, str):
                    # is the threshold a multiple of the data's variance of standard deviation?
                    if min_peak_intensity.find("var") != -1:
                        var = np.var(data_as_class.intensity.flatten())
                        min_peak_intensity = min_peak_intensity.strip("var").strip("*")
                        if min_peak_intensity == "":
                            min_peak_intensity = 1
                        min_peak_intensity = float(min_peak_intensity) * var

                    elif min_peak_intensity.find("std") != -1:
                        std = np.std(data_as_class.intensity.flatten())
                        min_peak_intensity = min_peak_intensity.strip("std").strip("*")
                        if min_peak_intensity == "":
                            min_peak_intensity = 1
                        min_peak_intensity = float(min_peak_intensity) * std

                    else:
                        err_str = "The 'min_peak_intensity' has to be a number or a multiple of 'var' or 'std'"
                        logger.critical(" ".join(map(str, [(err_str)])))
                        raise ValueError(err_str)

                if np.max(ave_intensity) <= min_peak_intensity:
                    # then there is no determinable peak in the data
                    logger.moreinfo(
                        "Not sufficient intensity in the chunked peaks to proceed with fitting."
                    )
                    logger.moreinfo(
                        f"most inense chunk {np.max(ave_intensity)} <= min allowed peak intensity {min_peak_intensity}"
                    )
                    # logger.moreinfo("Not sufficient intensity in the chunked peaks to proceed with fitting.")
                    # set step to -101 so that it is still negative at the end
                    step.append(-101)  # get to the end and void the fit
                    fout = master_params

            else:
                step.append(step[-1] + 11)

        if step[-1] >= 20:
            # FIX ME: Need to sort out use of lmfit parameters and new_params, below will not work currently for loaded
            # parameters as master_params not created from new_params yet. If possible use the json
            # import/export within lmfit.

            # Refit through full equation with all data for d,h,w,bg independently
            logger.moreinfo("Final fit solving for all parms...")

            # set all parameters to vary
            for k in range(len(settings_class.subfit_orders["background"])):
                param_str = "bg_c" + str(k)
                comp = "f"
                # set these parameters to vary
                master_params = lmm.vary_params(master_params, param_str, comp)
                # set part of these parameters to not vary

            comp_list, comp_names = pf.peak_components(include_profile=True)
            for k in range(peeks):
                param_str = "peak_" + str(k)

                for cp in range(len(comp_list)):
                    comp = comp_list[cp]
                    if (
                        comp_names[cp] + "_fixed"
                        in settings_class.subfit_orders["peak"][k]
                    ):
                        # set compenent not to vary
                        master_params = lmm.un_vary_this_params(
                            master_params, param_str, comp
                        )
                    else:
                        master_params = lmm.vary_params(master_params, param_str, comp)
                        # set part of these parameters to not vary
                        if isinstance(
                            settings_class.subfit_orders["peak"][k][comp_names[cp]],
                            list,
                        ):
                            master_params = lmm.un_vary_part_params(
                                master_params,
                                param_str,
                                comp,
                                sc.SeriesValues(
                                    settings_class.subfit_orders["peak"][k][
                                        comp_names[cp]
                                    ]
                                ),
                            )

            fout = lmm.fit_model(
                data_as_class,
                settings_class.subfit_orders,
                master_params,
                start_end=[data_as_class.azm_start, data_as_class.azm_end],
                fit_method=fit_method,
                weights=None,
                max_n_fev=max_f_eval,
            )
            master_params = fout.params

            if (
                fout.success is True
                and previous_params != None
                and has_huge_errors(
                    lmm.params_to_new_params(
                        master_params, orders=settings_class.subfit_orders
                    ),
                    min_ratio=large_errors,
                )
                == True
            ):
                logger.moreinfo(
                    "The fitting worked, but propagated params could have lead to rubbish fits (huge errors). Try again."
                )
                step.append(0)
                # clear previous_params so we can't get back here
                previous_params = None
            elif (
                fout.success is True
                and previous_params != None
                and has_value(master_params, val=None) == True
            ):
                logger.moreinfo(
                    "The fitting worked, but propagated params could have lead to rubbish fits (null values). Try again."
                )
                step.append(0)
                # clear previous_params so we can't get back here
                previous_params = None
            elif fout.success:  # True
                # it worked, errors are not massive, carry on
                step.append(step[-1] + 100)
                master_params = fout.params
            elif any(x == step[-1] for x in [20, 25]) and fout.success == 0:
                # increase iterations and try again
                step.append(step[-1] - 9)
                iterations = np.max((iterations, 5))
            elif any(x == step[-1] for x in [21, 26]) and fout.success == 0:
                # increase function evaluations and try again
                step.append(step[-1] - 9)
                max_f_eval = 2 * default_max_f_eval
            elif any(x == step[-1] for x in [22, 27]) and fout.success == 0:
                # just go bananas
                step.append(step[-1] - 9)
                max_f_eval = np.inf
            elif step[-1] == 28 and fout.success == 0:
                # we tried all the way with the previous fits.
                # get rid of them and start again
                step.append(0)
                previous_params = None
            else:  # 23,24,29 and other values >=30
                err_str = "The value of step here should not be achievable. Oops. \n The data is not fitting. Discard."
                logger.critical(" ".join(map(str, [(err_str)])))
                step.append(step[-1] - 200)

        if 0:  # report_status:
            print(f"    status: step = {step}")
            print(f"        time elaspsed = {time.time() - t_start}")

        if step[-1] < 0:
            # the fit is void and we need to exit.
            # make sure all the height values are nan so that we dont propagate rubbish
            comp_list, comp_names = pf.peak_components()
            for j in range(len(comp_list)):
                for i in range(peeks):
                    done = 0
                    n = 0
                    while done == 0:
                        try:
                            master_params[
                                "peak_" + str(i) + "_" + comp_list[j] + str(n)
                            ].value = None
                            n = n + 1
                        except:
                            # now we have run out of coefficients. So get the mean and then leave the loop.
                            done = 1

    # Write master_params to new_params dict object
    new_params = lmm.params_to_new_params(
        master_params, orders=settings_class.subfit_orders
    )

    # if step < 0:
    #     new_params.update({"FitProperties": fit_stats})
    # else:
    if step[-1] > 0:
        logger.effusive(" ".join(map(str, [("Final Coefficients")])))
        logger.effusive(" ".join(map(str, [(fout.fit_report(show_correl=False))])))
        # Record some stats
        logger.effusive(" ".join(map(str, [("number data %i" % fout.ndata)])))
        logger.effusive(" ".join(map(str, [("degrees of freedom %i" % fout.nfree)])))
        logger.effusive(" ".join(map(str, [("ChiSquared %f" % fout.chisqr)])))
        # get all correlation coefficients
        # FIX ME: this lists all the coefficients rather than just the unique half -- ie.
        # it contains corr(a,b) and corr(b,a)
        # FIX ME: should probably be a function
        # N.B. converted to a string so that it does not take up so many lines in output file.
        # Also why new lines are removed.
        correl = {}
        for key in master_params.keys():
            correl[key] = master_params[key].correl
        correl_str = json.dumps(correl)
        correl_str = correl_str.replace("\n", "")
        new_params.update({"correlation_coeffs": correl_str})

    # takes the maximum and minimum values to reflect data - rather then the inputs.
    new_params.update(
        {
            "range": [
                [data_as_class.tth_start, data_as_class.tth_end],
                data_as_class.conversion(
                    [data_as_class.tth_start, data_as_class.tth_end]
                ),
            ]
        }
    )
    if "note" in settings_class.subfit_orders:
        new_params.update({"note": settings_class.subfit_orders["note"]})

    # Elapsed time for fitting
    t_end = time.time()
    t_elapsed = t_end - t_start
    chunks_time = chunks_end - chunks_start
    step_str = json.dumps(step)
    step_str = step_str.replace("\n", "")
    # get the rest of the fit stats from the lmfit output.
    if step[-1] > 0:
        fit_stats = {
            "time-elapsed": t_elapsed,
            "chunks-time": chunks_time,
            "status": step_str,
            "sum-residuals-squared": np.sum(fout.residual**2),
            "function-evaluations": fout.nfev,
            "n-variables": fout.nvarys,
            "n-data": fout.ndata,
            "degree-of-freedom": fout.nfree,
            "ChiSq": fout.chisqr,
            "RedChiSq": fout.redchi,
            "aic": fout.aic,
            "bic": fout.bic,
        }

        # Store values related to the data.
        # Will be used by WriteFitMovie to keep colour scales consistent.
        new_params.update(
            {
                "DataProperties": {
                    "max": np.max(data_as_class.intensity),
                    "min": np.min(data_as_class.intensity),
                }
            }
        )
        new_params.update(
            {
                "ModelProperties": {
                    "max": np.max(fout.best_fit),
                    "min": np.min(fout.best_fit),
                }
            }
        )
        new_params.update(
            {
                "ResidualProperties": {
                    "max": np.max(data_as_class.intensity - fout.best_fit),
                    "min": np.min(data_as_class.intensity - fout.best_fit),
                }
            }
        )
        new_params.update(
            {
                "data_ranges": {
                    "data": {
                        "max": np.max(ma.filled(data_as_class.intensity, np.nan)),
                        "min": np.min(ma.filled(data_as_class.intensity, np.nan)),
                        "pt1percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan),
                            0.1,
                            method="closest_observation",
                        ),
                        "1percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan),
                            1,
                            method="closest_observation",
                        ),
                        "99percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan),
                            99,
                            method="closest_observation",
                        ),
                        "99pt9percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan),
                            99.9,
                            method="closest_observation",
                        ),
                    },
                    "model": {
                        "max": np.max(fout.best_fit),
                        "min": np.min(fout.best_fit),
                        "pt1percentile": np.nanpercentile(
                            fout.best_fit, 0.1, method="closest_observation"
                        ),
                        "1percentile": np.nanpercentile(
                            fout.best_fit, 1, method="closest_observation"
                        ),
                        "99percentile": np.nanpercentile(
                            fout.best_fit, 99, method="closest_observation"
                        ),
                        "99pt9percentile": np.nanpercentile(
                            fout.best_fit, 99.9, method="closest_observation"
                        ),
                    },
                    "residuals": {
                        "max": np.max(
                            ma.filled(data_as_class.intensity, np.nan) - fout.best_fit
                        ),
                        "min": np.min(
                            ma.filled(data_as_class.intensity, np.nan) - fout.best_fit
                        ),
                        "pt1percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan) - fout.best_fit,
                            0.1,
                            method="closest_observation",
                        ),
                        "1percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan) - fout.best_fit,
                            1,
                            method="closest_observation",
                        ),
                        "99percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan) - fout.best_fit,
                            99,
                            method="closest_observation",
                        ),
                        "99pt9percentile": np.nanpercentile(
                            ma.filled(data_as_class.intensity, np.nan) - fout.best_fit,
                            99.9,
                            method="closest_observation",
                        ),
                    },
                }
            }
        )
    else:
        # remove background estimates from the defunct model
        for i in range(len(new_params["background"])):
            new_params["background"][i] = [None] * len(new_params["background"][i])
            new_params["background_err"][i] = [None] * len(new_params["background"][i])

        fit_stats = {
            "time-elapsed": t_elapsed,
            "chunks-time": chunks_time,
            "status": step_str,
            "sum-residuals-squared": np.nan,
            "function-evaluations": np.nan,
            "n-variables": np.nan,
            "n-data": np.nan,
            "degree-of-freedom": np.nan,
            "ChiSq": np.nan,
            "RedChiSq": np.nan,
            "aic": np.nan,
            "bic": np.nan,
        }
        # Store values related to the data.
        # Will be used by WriteFitMovie to keep colour scales consistent.
        new_params.update(
            {
                "DataProperties": {
                    "max": np.max(data_as_class.intensity) if data_as_class.intensity.size > 0 else np.nan,
                    "min": np.min(data_as_class.intensity) if data_as_class.intensity.size > 0 else np.nan,
                }
            }
        )
        new_params.update({"ModelProperties": {"max": np.nan, "min": np.nan}})
        new_params.update({"ResidualProperties": {"max": np.nan, "min": np.nan}})
        new_params.update(
            {
                "data_ranges": {
                    "data": {
                        "max": np.max(data_as_class.intensity) if data_as_class.intensity.size > 0 else np.nan,
                        "min": np.min(data_as_class.intensity) if data_as_class.intensity.size > 0 else np.nan,
                        "pt1percentile": np.nanpercentile(data_as_class.intensity, 0.1) if data_as_class.intensity.size > 0 else np.nan,
                        "1percentile": np.nanpercentile(data_as_class.intensity, 1) if data_as_class.intensity.size > 0 else np.nan,
                        "99percentile": np.nanpercentile(data_as_class.intensity, 99) if data_as_class.intensity.size > 0 else np.nan,
                        "99pt9percentile": np.nanpercentile(
                            data_as_class.intensity, 99.9
                        ),
                    },
                    "model": {
                        "max": np.nan,
                        "min": np.nan,
                        "pt1percentile": np.nan,
                        "1percentile": np.nan,
                        "99percentile": np.nan,
                        "99pt9percentile": np.nan,
                    },
                    "residuals": {
                        "max": np.nan,
                        "min": np.nan,
                        "pt1percentile": np.nan,
                        "1percentile": np.nan,
                        "99percentile": np.nan,
                        "99pt9percentile": np.nan,
                    },
                }
            }
        )

    new_params.update({"FitProperties": fit_stats})

    # add peak names to new_params
    new_params.update({"PeakLabel": peak_string(settings_class.subfit_orders)})
    # add notes if they are present.
    if "note" in settings_class.subfit_orders:
        new_params.update({"note": settings_class.subfit_orders["note"]})

    # Plot results to check
    view = 0
    if (save_fit == 1 or view == 1 or logger.is_below_level("EFFUSIVE")) and step[
        -1
    ] > 0:
        logger.effusive(" ".join(map(str, [("Plotting results for fit...")])))

        orientation = data_as_class.plot_orientation
        if orientation == "vertical":
            fig = plt.figure(figsize=(4, 6))  # default figure size is [6.4, 4.8]
        else:
            fig = plt.figure()  # default figure size is [6.4, 4.8]
        fig = plot_FitAndModel(
            settings_class,
            data_as_class,
            modelresult=fout,
            param_lmfit=master_params,
            params_dict=new_params,
            figure=fig,
            orientation=orientation,
            plot_type="scatter",
        )
        title_str = peak_string(settings_class.subfit_orders) + "\n final fit"
        if "note" in settings_class.subfit_orders:
            title_str = title_str + " " + settings_class.subfit_orders["note"]

        # io.figure_suptitle_space(fig, topmargin=2)
        plt.suptitle(title_str)
        plt.tight_layout()

        if view == 1 or logger.is_below_level("EFFUSIVE"):
            plt.show()
        # save figures without overwriting old names
        if save_fit:
            filename = make_outfile_name(
                settings_class.subfit_filename,
                directory=settings_class.output_directory,
                orders=settings_class.subfit_orders,
                extension=".png",
                overwrite=False,
            )
            fig.savefig(filename)

        plt.close(fig)
        logger.effusive(" ".join(map(str, [("Done with Figure")])))

    # Save lmfit structure
    if save_fit:
        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            orders=settings_class.subfit_orders,
            extension=".sav",
            overwrite=True,
        )

        try:
            save_modelresult(fout, filename)

        except BaseException:
            logger.debug(" ".join(map(str, [("Cannot save lmfit object")])))
            # if os.path.isfile(filename):
            #     os.remove(filename)
            #     save_modelresult(out, filename)
            # else:
            #     logger.info(" ".join(map(str, [("File does not exist!")])))

    # return [new_params, fout]
    return new_params


def plot_FitAndModel(
    settings_class,
    data_as_class,
    modelresult=None,
    param_lmfit=None,
    params_dict=None,
    figure=None,
    debug=False,
    orientation="vertical",
    plot_type="default",
    plot_ColourRange=None,
):
    """
    Generates the figure for the fitted data.
    It calls the plot_fitted function of the data_class and returns the figure.
    Normally this has three axes: data, model and residuals but see Data class for specifics.

    Parameters
    ----------
    fig : TYPE
        DESCRIPTION.
    settings_class : TYPE
        DESCRIPTION.
    data_as_class : TYPE
        DESCRIPTION.
    modelresult : TYPE, optional
        DESCRIPTION. The default is None.
    param_lmfit : TYPE, optional
        DESCRIPTION. The default is None.
    params_dict : TYPE, optional
        DESCRIPTION. The default is None.
    figure : TYPE, optional
        DESCRIPTION. The default is None.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    figure : Figure
        Figure with axes generated by DATA_CLASS.plot_fitted.

    """

    # parse the input parameters
    if figure == None:
        figure = plt.figure()
    else:
        figure.clear(True)

    if modelresult:
        full_fit_intens = modelresult.best_fit
    else:
        if param_lmfit == None and params_dict == None:
            err_str = "No data has been passed."
            raise ValueError(err_str)
        elif param_lmfit == None:
            # initiate the model parameter set
            param_lmfit = lmm.initiate_all_params_for_fit(
                settings_class,
                data_as_class,
                values=params_dict,
                debug=debug,
            )
        elif params_dict == None:
            # make params dict from lmfit object
            params_dict = lmm.params_to_new_params(
                param_lmfit, orders=settings_class.subfit_orders
            )
        else:
            pass
            # all the required data are present in the required format.
            
        if not np.any(np.isfinite(list(param_lmfit.valuesdict().values())) == False):
            # if any values are not finite. Can be the case if background has a gradient (which does not have limits) AND the fit fails.
            gmodel = Model(
                lmm.peaks_model,
                independent_vars=["two_theta", "azimuth"],
                data_class=data_as_class,
                orders=settings_class.subfit_orders,
                start_end=[data_as_class.azm_start, data_as_class.azm_end],
            )
            full_fit_intens = gmodel.eval(
                params=param_lmfit,
                two_theta=data_as_class.tth.flatten(),
                azimuth=data_as_class.azm.flatten(),
            )
        else:
            full_fit_intens = np.zeros(data_as_class.tth.shape)

    azi_plot = np.unique(data_as_class.azm.flatten())
    if data_as_class.continuous_azm:
        # if there are lots and lots of data make the list shorter
        azi_plot = np.array(
            list(
                range(
                    np.int_(data_as_class.azm_start), np.int_(data_as_class.azm_end), 2
                )
            )
        )
    peeks = len(params_dict["peak"])
    fit_centroid = []
    for i in range(peeks):
        param = params_dict["peak"][i]["d-space"]
        param_type = params_dict["peak"][i]["d-space_type"]
        if not None in param:
            fit_centroid.append(
                data_as_class.conversion(
                    sf.coefficient_expand(azi_plot, param=param, coeff_type=param_type),
                    azi_plot,
                    reverse=1,
                )
            )
        else:
            fit_centroid.append(
                np.full(azi_plot.shape, np.nan)
            )  # np.zeros(azi_plot.shape))

    # plot the data and the fit
    figure = data_as_class.plot_fitted(
        fig_plot=figure,
        model=full_fit_intens,
        fit_centroid=[azi_plot, fit_centroid],
        orientation=orientation,
        plot_type=plot_type,
        plot_ColourRange=plot_ColourRange,
    )
    return figure
