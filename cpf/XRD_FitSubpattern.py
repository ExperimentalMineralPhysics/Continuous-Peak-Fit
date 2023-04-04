#!/usr/bin/env python

__all__ = ["fit_sub_pattern"]

# CPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

import sys
import time
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Model
from lmfit.model import save_modelresult  # , load_modelresult
import json
import cpf.series_functions as sf
import cpf.peak_functions as pf
import cpf.lmfit_model as lmm
import cpf.IO_functions as io
from cpf.fitsubpattern_chunks import fit_chunks, fit_series

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
                not sf.fourier_order(orders["peak"][y]["profile_fixed"])
                == orders["peak"][y]["profile"]
            ):
                print(
                    "Peak "
                    + str(y)
                    + ": The order of the profile does not match that of the fixed profile. "
                    "Changing profile to match fixed profile."
                )
                orders["peak"][y]["profile"] = sf.fourier_order(
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
            coeff_type = sf.params_get_type(previous_params, param, peak=y)
            if coeff_type != 5:  # if parameters are not independent
                if np.max(orders["peak"][y][param]) > sf.fourier_order(
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
                elif np.max(orders["peak"][y][param]) < sf.fourier_order(
                    previous_params["peak"][y][param]
                ):
                    # print 'smaller'
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
        sf.params_get_type(previous_params, "background") != 5
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
                if np.max(orders["background"][y]) > sf.fourier_order(
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
                elif np.max(orders["background"][y]) < sf.fourier_order(
                    previous_params["background"][y]
                ):
                    # print 'smaller'
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
    n_azimuths = len(np.unique(azimu))
    max_coeff = 0
    choice_list = ["d-space", "height", "width", "profile"]
    for y in range(peeks):
        # loop over parameters
        for param in choice_list:
            coeff_type = sf.params_get_type(orders, param, peak=y)
            if coeff_type != 5:  # if parameters are not independent
                max_coeff = np.max(
                    [max_coeff, np.max(orders["peak"][y][param]) * 2 + 1]
                )
    param = "background"
    for y in range(np.max([len(orders["background"])])):
        coeff_type = sf.params_get_type(orders, param, peak=y)
        if coeff_type != 5:  # if parameters are not independent
            max_coeff = np.max([max_coeff, np.max(orders["background"][y]) * 2 + 1])
    # print(max_coeff, len(np.unique(azimu)))
    if max_coeff > len(np.unique(azimu)):
        err_str = (
            "The maximum order, %i, needs more coefficients than the number of unique azimuths, %i. "
            "It is not possible to fit."
            % (sf.get_order_from_coeff(max_coeff), len(np.unique(azimu)))
        )
        raise ValueError(err_str)


def fit_sub_pattern(
    data_as_class=None,
    settings_as_class=None,
    previous_params=None,
    save_fit=False,
    debug=False,
    refine=True,
    iterations=2,
    fit_method=None,
    mode="fit",
    cascade=False,
    intensity_threshold = 0,
    large_errors = 3
):
    """
    Perform the various fitting stages to the data
    :param fit_method:
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
    :return:
    """

    # set a limit to the maximum number of function evaluations.
    # make variable in case need more iterations for other data
    default_max_f_eval = 400

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
    # FIX ME: this doesn't seem to work. the data type is got by print(intens.dtype)

    # FIX ME: need to match orders of arrays to previous numbers.
    # if no parameters
    #     get orders of arrays
    # else if no orders
    #     use parameters as given
    # else
    #     change the size of the parameter arrays to match orders

    # Define the number of peaks to be fitted
    peeks = len(settings_as_class.subfit_orders["peak"])

    # DMF: FIX ME: orders used before checked?
    # if settings_as_class.subfit_orders:
    # If exists in input file setup peaks as appropriate for fitting
    # background_type = "order"
    # bg_order = orders["background"]
    # orders, backgnd = order_set_peaks(settings_as_class.subfit_orders, peeks, settings_as_class.subfit_orders["background"])

    if previous_params:
        # check if the previous fit was 'good' i.e. constrains no 'null' values.
        # N.B. null values in json file are read in as None
        clean = io.any_terms_null(previous_params, val_to_find=None)
        clean = io.any_errors_huge(previous_params, large_errors=large_errors, clean=clean)
        if clean == 0:
            # the previous fit has problems so discard it
            print(
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
            np.min(fit_centroid) < settings_as_class.subfit_orders["range"][0] or 
            np.min(fit_centroid) > settings_as_class.subfit_orders["range"][1]
        ):
            print("fit d-spacing limits are out of bounds; discarding the fit and starting again.")
            previous_params = None


    if previous_params:
        # initiate values for while loop.
        step = 5
    else:
        step = 0

    if previous_params and settings_as_class.subfit_orders:
        # If we have both, order takes precedence so update previous_params to match
        previous_params = update_previous_params_from_orders(
            peeks, previous_params, settings_as_class.subfit_orders
        )

    # FIX ME: can we have a situation with no orders or previous_params?

    # check the number of unique azimuths is greater than the number of coefficients.
    check_num_azimuths(peeks, data_as_class.azm, settings_as_class.subfit_orders)

    # Start fitting loops
    while step >= 0 and step <= 100:
        # 100 is arbitrarily large number and does not reflect the num. of actual steps/stages in the process.
        # While loops are used so that if the fit is rubbish we can go back and improve it by repeating earlier steps.
        # for chunks step <= 9 and for refine <= 19
        # we are using increments of 10 so that it is possible to record different states or routes after the final fit.

        if step <= 9:
            # generate chunks and initial fits
            # or parse previous fits into correct data structure

            # Measure the time taken to do the chunks, the elapsed time during fitting.
            # To help decide which is the best peak parameters.
            chunks_start = time.time()

            # If there is not a previous fit -- Fit data in azimuthal chunks for d.
            if not previous_params:

                # fit the data in chunks by azimuth.
                chunk_fits, chunk_positions = fit_chunks(
                    data_as_class,
                    settings_as_class,
                    mode=mode,
                    # cascade=cascade,
                    debug=debug,
                    fit_method=fit_method
                )
                # using manual guesses ("PeakPositionSelection") if they exist.

                if mode != "fit":  # cascade==True:
                    # some cascade option. so exit returning values.
                    return chunk_fits, chunk_positions

                # Fitting stage 2:
                print("\n Performing Series fits...")
                # Feed each d_0,h,w into coefficient function to get fit for coeficient component
                # get lmfit parameters as output.

                # initiate the model parameter set
                master_params = lmm.initiate_all_params_for_fit(
                    settings_as_class, data_as_class, debug=debug
                )

                # fit the chunk values with fourier/spline series.
                # iterate over each parameter in turn
                master_params = fit_series(
                    master_params,
                    (chunk_fits, chunk_positions),
                    settings_as_class,
                    start_end = [data_as_class.azm_start, data_as_class.azm_end],
                    debug=debug,
                    save_fit=save_fit,
                )
                
                #get mean height of chunked peaks
                ave_intensity = []
                for k in range(peeks):
                    ave_intensity.append(sf.get_series_mean(master_params, "peak_"+str(k), comp="h"))
                if np.max(ave_intensity) <= intensity_threshold:
                    #then there is no determinable peak in the data
                    #set step to -11 so that it is still negative at the end
                    step = -11 #get to the end and void the fit
                    fout=master_params

            else:
                # FIX ME: This should load the saved lmfit Parameter class object.
                # Need to initiate usage (save and load).
                # For now have propagated use of new_params so re-instantiate the master_params object below,
                # which is clunky.
                print("Using previously fitted parameters and propagating fit")

                # FIX ME: need to confirm the number of parameters matches the orders of the fits.

                # initiate the model parameter set
                master_params = lmm.initiate_all_params_for_fit(
                    settings_as_class,
                    data_as_class,
                    values=previous_params,
                    debug=debug,
                )

            chunks_end = time.time()
            step = step + 10

        if step >= 10:

            # if refine or step>=10 or not PreviousParams:
            if refine or step != 10 or step != 15 and iterations >= 1:
                # Iterate over each parameter series in turn.
                print(
                    "\nRe-fitting for d, h, w, +/-p, bg separately... will refine %i time(s)\n"
                    % iterations
                )
                for j in range(iterations):
                    for k in range(len(settings_as_class.subfit_orders["background"])):
                        param_str = "bg_c" + str(k)
                        comp = "f"
                        # set other parameters to not vary
                        master_params = lmm.un_vary_params(
                            master_params, param_str, comp
                        )
                        # set these parameters to vary
                        master_params = lmm.vary_params(master_params, param_str, comp)
                        # set part of these parameters to not vary
                        # master_params = ff.un_vary_part_params(master_params, param_str, comp,
                        # orders['background'][k])
                        fout = lmm.fit_model(
                            data_as_class,
                            settings_as_class.subfit_orders,
                            master_params,
                            start_end = [data_as_class.azm_start, data_as_class.azm_end],
                            fit_method=None,
                            weights=None,
                            max_n_fev=default_max_f_eval,
                        )
                        master_params = fout.params


                    # iterate over the peak components in order of size (height)
                    # get mean height of peak from parameters
                    peak_mean = []
                    for k in range(peeks):
                        peak_mean.append(sf.get_series_mean(master_params, "peak_"+str(k), comp="h"))
                    #get order to process peaks in -- reverse order starting with the most intense.
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
                                comp_names[cp]+"_fixed" in settings_as_class.subfit_orders["peak"][k]
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
                                    settings_as_class.subfit_orders["peak"][k][
                                        comp_names[cp]
                                    ],
                                    list,
                                ):
                                    master_params = lmm.un_vary_part_params(
                                        master_params,
                                        param_str,
                                        comp,
                                        settings_as_class.subfit_orders["peak"][k][
                                            comp_names[cp]
                                        ],
                                    )
                                if comp == "h":
                                    refine_max_f_eval = 5 * peeks * default_max_f_eval
                                else:
                                    refine_max_f_eval = default_max_f_eval
                                fout = lmm.fit_model(
                                    data_as_class,
                                    settings_as_class.subfit_orders,
                                    master_params,
                                    start_end = [data_as_class.azm_start, data_as_class.azm_end],
                                    fit_method=None,
                                    weights=None,
                                    max_n_fev=refine_max_f_eval,
                                )
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
            for k in range(len(settings_as_class.subfit_orders["background"])):
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
                        comp_names[cp]+"_fixed" in settings_as_class.subfit_orders["peak"][k]
                    ):
                        #set compenent not to vary
                        master_params = lmm.un_vary_this_params(
                            master_params, param_str, comp
                        )
                    else:
                        master_params = lmm.vary_params(
                            master_params, param_str, comp
                        )
                        # set part of these parameters to not vary
                        if isinstance(
                            settings_as_class.subfit_orders["peak"][k][comp_names[cp]], list
                        ):
                            master_params = lmm.un_vary_part_params(
                                master_params,
                                param_str,
                                comp,
                                settings_as_class.subfit_orders["peak"][k][comp_names[cp]],
                            )

            fout = lmm.fit_model(
                data_as_class,
                settings_as_class.subfit_orders,
                master_params,
                start_end = [data_as_class.azm_start, data_as_class.azm_end],
                fit_method=None,
                weights=None,
                max_n_fev=default_max_f_eval,
            )
            master_params = fout.params

            if fout.success == 1:
                # it worked, carry on
                step = step + 100
                master_params = fout.params
            elif step == 24 and fout.success == 0:
                raise ValueError(
                    "Oh Dear. It should not be possible to get here. Something has gone very wrong with the fitting."
                )
            elif step == 29 and fout.success == 0:
                step = 0  # go back to the start, discard PreviousParams and do the chunks for this data set.
            elif any(x == step for x in [20, 25]) and fout.success == 0:
                step = step - 7
                iterations = np.max((iterations, 3))
            elif any(x == step for x in [21, 22, 23, 26, 27, 28]) and fout.success == 0:
                step = step - 9
                if any(x == step for x in [23, 28]):
                    iterations = np.max((iterations, 3))
            else:
                raise ValueError("The value of step here is not possible. Oops.")

            # FIX ME: we could make an option for output the chunks without any Fourier/global fitting.
            # Would this be useful?

        if step < 0:
            #the fit is void and we need to exit.
            #make sure all the height values are nan so that we dont propagate rubbish
            comp_list, comp_names = pf.peak_components()
            for j in range(len(comp_list)):
                for i in range(peeks):
                    done = 0
                    n = 0
                    while done==0:
                        try:
                            master_params["peak_"+str(i)+"_"+comp_list[j]+str(n)].value = None
                            n = n+1
                        except:
                            # now we have run out of coefficients. So get the mean and then leave the loop.
                            done = 1
    
            
    # Write master_params to new_params dict object
    new_params = lmm.params_to_new_params(
        master_params, orders=settings_as_class.subfit_orders
    )

    # if step < 0:
    #     new_params.update({"FitProperties": fit_stats})
    # else:
    if step > 0:        
        print("\nFinal Coefficients\n")
        print(fout.fit_report(show_correl=False))
    
        # Print some stats
        print("number data", fout.ndata)
        print("degrees of freedom", fout.nfree)
        print("ChiSquared", fout.chisqr)
    

    
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
                [data_as_class.tth.min(), data_as_class.tth.max()],
                [data_as_class.dspace.min(), data_as_class.dspace.max()],
            ]
        }
    )
    if "note" in settings_as_class.subfit_orders:
        new_params.update({"note": settings_as_class.subfit_orders["note"]})

    # Elapsed time for fitting
    t_end = time.time()
    t_elapsed = t_end - t_start
    chunks_time = chunks_end - chunks_start
    # get the rest of the fit stats from the lmfit output.
    if step > 0:
        fit_stats = {
            "time-elapsed": t_elapsed,
            "chunks-time": chunks_time,
            "status": step,
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
    else:
        fit_stats = {
            "time-elapsed": t_elapsed,
            "chunks-time": chunks_time,
            "status": step,
        }
        
    new_params.update({"FitProperties": fit_stats})
    
    # add peak names to new_params
    new_params.update({"PeakLabel": io.peak_string(settings_as_class.subfit_orders)})
    
    # Plot results to check
    view = 1
    if (save_fit == 1 or view == 1 or debug) and step>0:
        print("\nPlotting results for fit...\n")
        
        fig = plt.figure()
        fig = plot_FitAndModel(settings_as_class, data_as_class, param_lmfit=master_params, params_dict=new_params, figure=fig)
        title_str = io.peak_string(settings_as_class.subfit_orders) + "; final fit"
        if "note" in settings_as_class.subfit_orders:
            title_str = title_str + " " + settings_as_class.subfit_orders["note"]
        plt.suptitle(title_str)
        #plt.tight_layout()

        if view == 1 or debug:
            plt.show()
        # save figures without overwriting old names
        if save_fit:
            filename = io.make_outfile_name(
                settings_as_class.subfit_filename,
                directory=settings_as_class.output_directory,
                orders=settings_as_class.subfit_orders,
                extension=".png",
                overwrite=False,
            )
            plt.savefig(filename)

        if view == 1 or debug:
            plt.show()

        plt.close(fig)
        print("Done with Figure")

    # Save lmfit structure
    if save_fit:
        filename = io.make_outfile_name(
            settings_as_class.subfit_filename,
            directory=settings_as_class.output_directory,
            orders=settings_as_class.subfit_orders,
            extension=".sav",
            overwrite=True,
        )

        try:
            save_modelresult(fout, filename)

        except BaseException:
            print("Cannot save lmfit object")
            # if os.path.isfile(filename):
            #     os.remove(filename)
            #     save_modelresult(out, filename)
            # else:
            #     print("File does not exist!")

    return [new_params, fout]


def plot_FitAndModel(settings_as_class, data_as_class, param_lmfit=None, params_dict=None, figure=None, debug=False):
    """
    Generates the figure for the fitted data. 
    It calls the plot_fitted function of the data_class and returns the figure. 
    Normally this has three axes: data, model and residuals but see Data class for specifics.

    Parameters
    ----------
    fig : TYPE
        DESCRIPTION.
    settings_as_class : TYPE
        DESCRIPTION.
    data_as_class : TYPE
        DESCRIPTION.
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
    

    #parse the input parameters
    if figure == None:
        figure = plt.figure()
    else:
        figure.clear(True)
        
    if param_lmfit == None and params_dict == None: 
        err_str = "No data has been passed."
        raise ValueError(err_str) 
    elif param_lmfit == None:
        # initiate the model parameter set
        param_lmfit = lmm.initiate_all_params_for_fit(
            settings_as_class,
            data_as_class,
            values=params_dict,
            debug=debug,
        )
    elif params_dict == None: 
        # make params dict from lmfit object
        params_dict = lmm.params_to_new_params(
            param_lmfit, orders=settings_as_class.subfit_orders
        )
    else:
        pass
        #all the required data are present in the required format.
        
    y_lims = [data_as_class.azm_start, data_as_class.azm_end]
    gmodel = Model(
        lmm.peaks_model,
        independent_vars=["two_theta", "azimuth"],
        data_class=data_as_class,
        orders=settings_as_class.subfit_orders,
        start_end = [data_as_class.azm_start, data_as_class.azm_end],
    )
    full_fit_intens = gmodel.eval(
        params=param_lmfit,
        two_theta=data_as_class.tth.flatten(),
        azimuth=data_as_class.azm.flatten(),
    )

    azi_plot = np.unique(data_as_class.azm.flatten())
    if data_as_class.continuous_azm:
        # if there are lots and lots of data make the list shorter
        azi_plot = np.array(list(range(np.int(y_lims[0]), np.int(y_lims[1]), 2)))

    peeks = len(params_dict["peak"])
    fit_centroid = []
    for i in range(peeks):
        param = params_dict["peak"][i]["d-space"]
        param_type = params_dict["peak"][i]["d-space_type"]
        fit_centroid.append(
            data_as_class.conversion(
                sf.coefficient_expand(azi_plot, param=param, coeff_type=param_type),
                azi_plot,
                reverse=1,
            )
        )
        
    # plot the data and the fit
    data_as_class.plot_fitted(
        fig_plot=figure, model=full_fit_intens, fit_centroid=[azi_plot, fit_centroid]
    )
    return figure