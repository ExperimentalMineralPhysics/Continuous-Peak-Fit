#!/usr/bin/env python

__all__ = ["get_manual_guesses", "get_chunk_background_guess", "fit_chunks"]

# CPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

from copy import deepcopy

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from lmfit import Model, Parameters

import cpf.histograms as hist
import cpf.lmfit_model as lmm
import cpf.series_constraints as sc
import cpf.series_functions as sf
from cpf.util.io import make_outfile_name, peak_string
from cpf.util.logging import get_logger

logger = get_logger("cpf.fitsubpattern_chunks")


def get_manual_guesses(
    settings_as_class, data_as_class, return_converted=False, debug=False
):
    """
    Calculate series for set of 'manual guesses' i.e. pre-defined peak centers.

    Parameters
    ----------
    settings_as_class : cpf settings class
        DESCRIPTION.
    data_as_class : cpf data class
        DESCRIPTION.
    return_converted : bool
        Return series in d-spacing [if = True] or
        collection dimension (e.g. two theta, energy) [if = False]
        The default is False
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    dfour : list
        List of series coefficients for each peak.
    """

    # get guesses and limits in the correct parameter
    peak_pos_guesses = np.array(
        settings_as_class.subfit_orders["PeakPositionSelection"]
    )
    lims = lmm.parse_bounds(
        settings_as_class.fit_bounds, data_as_class, 0, 0, param=["d-space"]
    )  # limits alwasys in converted unit.
    limits_range = lims["d-space"]
    if return_converted == True:
        # limits alwasys in converted unit no need to chamge.
        if np.max(peak_pos_guesses[:, 2]) > np.max(limits_range) or np.min(
            peak_pos_guesses[:, 2]
        ) < np.min(limits_range):
            # then the peak_pos_guesses are not in d-spacing
            # therefore change
            peak_pos_guesses[:, 2] = data_as_class.conversion(
                peak_pos_guesses[:, 2], azm=None, reverse=False
            )

    else:  # return_converted == False:
        # limits neeed to be converted
        limits_range = data_as_class.conversion(limits_range, azm=None, reverse=True)
        # then check guesses
        if np.max(peak_pos_guesses[:, 2]) > np.max(limits_range) or np.min(
            peak_pos_guesses[:, 2]
        ) < np.min(limits_range):
            # then the peak_pos_guesses are not in collected units
            # therefore change
            peak_pos_guesses[:, 2] = data_as_class.conversion(
                peak_pos_guesses[:, 2], azm=None, reverse=True
            )

    settings_as_class.validate_position_selection(
        peak_set=settings_as_class.subfit_order_position, report=False
    )

    # Fit Fourier series to two-theta/d-spacing for each peak
    dfour = []
    for j in range(len(settings_as_class.subfit_orders["peak"])):
        peek = j + 1
        peak_pos_guess = peak_pos_guesses[peak_pos_guesses[:, 0] == peek, 1:]
        param_str = "peak_" + str(j)
        comp = "d"
        # get coefficient type
        coeff_type = sf.get_params_type(settings_as_class.subfit_orders, comp, peak=j)

        # for guesses make sure there are not too many coefficients.
        n_coeff = sf.get_number_coeff(settings_as_class.subfit_orders, comp, peak=j)

        if n_coeff > len(peak_pos_guess[:, 0]):
            o = int(np.floor(len(peak_pos_guess[:, 0]) / 2 - 1))
            # FIX ME: this should be a function!! and should check the series type.
        else:
            o = settings_as_class.subfit_orders["peak"][j]["d-space"]

        temp_param = Parameters()
        temp_param = lmm.initiate_params(
            inp_param=temp_param,
            param_str=param_str,
            coeff_type=coeff_type,
            comp=comp,
            trig_orders=sc.SeriesValues(o),
            limits=limits_range,
            value=peak_pos_guess[0, 1],
        )
        fout = lmm.coefficient_fit(
            azimuth=peak_pos_guess[:, 0],
            ydata=peak_pos_guess[:, 1],
            inp_param=temp_param,
            param_str=param_str + "_" + comp,
            fit_method="leastsq",
        )
        temp_param = fout.params

        logger.log_lmfit_obj(temp_param, level="DEBUG", space=True)
        dfour.append(lmm.gather_param_errs_to_list(temp_param, param_str, comp))
    return dfour


def get_chunk_background_guess(settings_as_class, data_chunk, n=1, debug=False):
    """
    :param data_chunk_class:
    :param settings_as_class:
    :param chunks:
    :param count:
    :param n:
    :param debug:
    :return background_guess:

    FIXME: background_type has been removed/depreciated. It is no longer needed. It should be replaced by background_fixed
    """

    cnk_i = data_chunk[0]
    cnk_tth = data_chunk[1]

    # Get indices of sorted two theta values excluding the masked values
    # tth_ord = ma.argsort(data_chunk_class.tth.compressed())
    tth_ord = np.argsort(cnk_tth)

    background_guess = [
        [0.0] for i in range(len(settings_as_class.subfit_orders["background"]))
    ]

    if "background_fixed" in settings_as_class.fit_orders:
        # this is not yet implemented but the option exists here where it is needed first.
        raise NotImplementedError

    else:
        if len(settings_as_class.subfit_orders["background"]) == 1:
            # assume background is present either at the left or right edges of the region.
            # FIXME: this assumes a positive peak. It may not work for absorption peaks. Not tested though.
            background_guess[0][0] = np.min(
                [
                    np.mean(cnk_i[tth_ord[:n]]),
                    np.mean(cnk_i[tth_ord[-n:]]),
                    # np.mean(data_chunk_class.intensity.compressed()[tth_ord[:n]]),
                    # np.mean(data_chunk_class.intensity.compressed()[tth_ord[-n:]]),
                ]
            )
        else:  # len(orders['background']) > 1:
            # first value (offset) is mean of left-hand values
            background_guess[0][0] = np.mean(cnk_i[tth_ord[:n]])
            # data_chunk_class.intensity.compressed()[tth_ord[:n]]
            # )
            # if there are more, then calculate a gradient guess.
            background_guess[1][0] = (
                np.mean(cnk_i[tth_ord[-n:]]) - np.mean(cnk_i[tth_ord[:n]])
                # np.mean(data_chunk_class.intensity.compressed()[tth_ord[-n:]])
                # - np.mean(data_chunk_class.intensity.compressed()[tth_ord[:n]])
            ) / (
                cnk_tth[tth_ord[-1]] - cnk_tth[tth_ord[0]]
                # data_chunk_class.tth.compressed()[tth_ord[-1]]
                # - data_chunk_class.tth.compressed()[tth_ord[0]]
            )
            # leave all higher order terms as 0 -- assume they are small.

    return background_guess


def get_chunk_peak_guesses_new(
    settings_as_class,
    data_chunk,
    background_guess,
    n,
    w_guess_fraction,
    dfour,
    debug=False,
):
    """


    Parameters
    ----------
    settings_as_class : TYPE
        DESCRIPTION.
    data_chunk : TYPE
        DESCRIPTION.
    background_guess : TYPE
        DESCRIPTION.
    n : TYPE
        DESCRIPTION.
    w_guess_fraction : TYPE
        DESCRIPTION.
    dfour : TYPE
        DESCRIPTION.
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    c_guess : TYPE
        DESCRIPTION.
    h_guess : TYPE
        DESCRIPTION.
    w_guess : TYPE
        DESCRIPTION.
    p_guess : TYPE
        DESCRIPTION.

    """
    """


    :param data_as_class:
    :param settings_as_class:
    :param background_guess:
    :param chunks:
    :param count:
    :param n:
    :param w_guess_fraction:
    :param dfour:
    :param debug:
    :return peaks, limits, p_fixed:
    """

    cnk_i = data_chunk[0]
    cnk_tth = data_chunk[1]
    cnk_azm = data_chunk[2]

    guess_cent = []
    guess_h = []
    guess_w = []
    guess_p = []

    for k in range(len(settings_as_class.subfit_orders["peak"])):
        # FIX ME: altered from locals call as now passed as None - check!
        if dfour:  # If the positions have been pre-guessed extract values from dfour.
            # if debug and k == 0:
            if k == 0:
                logger.debug(" ".join(map(str, [("dfour in locals")])))

            coeff_type = sf.get_params_type(
                settings_as_class.subfit_orders, "d", peak=k
            )
            pos_guess = sf.coefficient_expand(
                np.mean(cnk_azm),
                dfour[k][0],
                coeff_type=coeff_type,
            )
            # t_th_guess = data_chunk_class.conversion(
            #     d_guess, azm=np.mean(cnk_azm), reverse=True
            # )
            # Finds the index of the closest two-theta to the d-spacing input
            idx = (np.abs(cnk_tth - pos_guess)).argmin()
            # FIX ME: The mean of a number of the smallest values would be more stable.

            # Height is intensity of the closest pixel in d-spacing - background at that position
            h_guess = cnk_i[idx]
            if len(background_guess) > 1:
                h_guess = h_guess - (
                    background_guess[0][0]
                    + background_guess[1][0] * (cnk_tth[idx] - cnk_tth.min())
                )
            elif len(background_guess) == 1:
                h_guess = h_guess - background_guess[0][0]

        else:  # 'guess' guesses from data slice - assuming a single peak
            # find the brightest pixels
            # get index of nth brightest pixel.
            # idx = np.argsort(data_chunk_class.intensity.compressed())[-n:][0]
            idx = np.argsort(cnk_i)[-n:][0]

            # height guess is nth highest intensity - background guess at this position
            # h_guess = (data_chunk_class.intensity.compressed())[idx]
            h_guess = cnk_i[idx]
            if len(background_guess) > 1:
                h_guess = h_guess - (
                    background_guess[0][0]
                    + background_guess[1][0]
                    * (
                        cnk_tth[idx]
                        # data_chunk_class.tth.compressed()[idx]
                        - cnk_tth.min()
                    )
                )
            elif len(background_guess) == 1:
                h_guess = h_guess - background_guess[0][0]

            # d-spacing of the highest nth intensity pixel
            # d_guess = data_chunk_class.dspace.compressed()[idx]
            pos_guess = cnk_tth[idx]

        # w_guess is fractional width of the data range
        w_guess = (
            (np.max(cnk_tth) - np.min(cnk_tth))
            / w_guess_fraction
            / len(settings_as_class.subfit_orders["peak"])
        )
        # FIX ME: This is a bit crude. Is there a better way to do it?

        if "profile_fixed" in settings_as_class.subfit_orders["peak"][k]:
            # calculate profile at this orientation
            coeff_type = sf.get_params_type(
                settings_as_class.subfit_orders, "p", peak=k
            )
            if "symmetry" in settings_as_class.subfit_orders["peak"][k]:
                symm = settings_as_class.subfit_orders["peak"][k]["symmetry"]
            else:
                symm = 1
            p_guess = sf.coefficient_expand(
                np.mean(cnk_azm) * symm,
                param=settings_as_class.subfit_orders["peak"][k]["profile_fixed"],
                coeff_type=coeff_type,
            )
        else:
            p_guess = np.mean(settings_as_class.fit_bounds["profile"])

        guess_cent.append(pos_guess)
        guess_h.append(h_guess)
        guess_w.append(w_guess)
        guess_p.append(p_guess)

    return guess_cent, guess_h, guess_w, guess_p


def fit_chunks(
    data_as_class,
    settings_as_class,
    fit_method="least_squares",
    mode="fit",
    histogram_type=None,  # "width",
    histogram_bins=None,
    max_n_f_eval=400,
    save_fit=False,
    debug=False,
):
    """
    Take the raw data, fit the chunks and return the chunk fits
    :param fit_method:
    :param data_as_class:
    :param settings_as_class:
    :param save_fit:
    :param debug:
    :param fit_method:
    :return chunk_params:
    """

    peeks = len(settings_as_class.subfit_orders["peak"])
    bg_length = len(settings_as_class.subfit_orders["background"])

    if "PeakPositionSelection" not in settings_as_class.subfit_orders:
        if peeks > 1:
            p = "peaks"
        else:
            p = "peak"
        logger.debug(
            " ".join(
                map(
                    str,
                    [
                        (
                            "No previous fits or selections. Assuming %i %s and group fitting in azimuth..."
                            % (peeks, p)
                        )
                    ],
                )
            )
        )
        dfour = None
        # FIXME: passed as None now - check functions as expected
    else:
        logger.debug(
            " ".join(map(str, [("Using manual selections for initial guesses...")]))
        )
        # FIXME: Need to check that the num. of peaks for which we have parameters is the same as the
        # number of peaks guessed at.
        # dfour = get_manual_guesses(peeks, orders, bounds, twotheta, debug=None)
        dfour = get_manual_guesses(
            settings_as_class, data_as_class, return_converted=False, debug=debug
        )

    # Get chunks according to detector type.
    if mode != "fit":  # cascade
        chunks, chunk_bounds, azichunks = data_as_class.bins(
            settings_as_class, cascade=True
        )
    else:
        chunks, chunk_bounds, azichunks = data_as_class.bins(settings_as_class)
    # Final output list of azimuths with corresponding twotheta_0,h,w

    # setup arrays
    new_azi_chunks = []
    if mode == "fit" or mode == "cascade":
        # bunch of settings
        # only compute the fit and save it if there are sufficient data
        min_dat = 21  # minimum data in each chunk
        n = 5  # how many data to use in constructing guesses.
        w_guess_fraction = 10.0
        # N.B. if min_dat/n is too large the guesses will become very strange.
        # FIXME: These should be input variables (for use in extremes)
        # FIXME: possibly do this better using scipy.signal.find or scipy.signal.find_peaks_cwt

    else:
        pass

    out_vals = {}
    out_vals["chunks"] = []
    out_vals["d"] = [[] for _ in range(peeks)]
    out_vals["d_err"] = [[] for _ in range(peeks)]
    out_vals["h"] = [[] for _ in range(peeks)]
    out_vals["h_err"] = [[] for _ in range(peeks)]
    out_vals["w"] = [[] for _ in range(peeks)]
    out_vals["w_err"] = [[] for _ in range(peeks)]
    out_vals["p"] = [[] for _ in range(peeks)]
    out_vals["p_err"] = [[] for _ in range(peeks)]
    out_vals["bg"] = [[] for _ in range(bg_length)]
    out_vals["bg_err"] = [[] for _ in range(bg_length)]

    out_vals["peak"] = peak_string(settings_as_class.subfit_orders)

    for j in range(len(chunks)):
        # logger.info(" ".join(map(str, [('\nFitting to data chunk ' + str(j + 1) + ' of ' + str(len(chunks)) + '\n')])))
        logger.debug(f"Fitting to data chunk {j + 1} of {len(chunks)}")

        # make data class for chunks.
        # reduce data to a subset using azi_bounds
        chunk_data = data_as_class.duplicate(
            azi_bounds=chunk_bounds[j], as_masked=False
        )
        chunk_intensity = chunk_data.intensity
        chunk_tth = chunk_data.tth
        chunk_azm = chunk_data.azm

        # find other output from intensities
        if mode == "maxima":
            # get maximum from each chunk
            out_vals["h"][0].append(np.max(chunk_intensity))
            out_vals["chunks"].append(azichunks[j])
            new_azi_chunks.append(azichunks[j])

        elif mode == "range":
            # get maximum from each chunk
            out_vals["h"][0].append(np.max(chunk_intensity) - np.min(chunk_intensity))
            out_vals["chunks"].append(azichunks[j])
            new_azi_chunks.append(azichunks[j])

        elif mode == "98thpercentile":
            # FIXME: this is crude and could be done better -- i.e. with a switch for the percentile
            # get 98th percentils from each chunk
            raise NotImplementedError

        elif mode == "fit" or mode == "cascade":
            # Define parameters to pass to fit
            params = Parameters()

            # if ma.MaskedArray.count(chunk_data.intensity) >= min_dat:
            if len(chunk_intensity) >= min_dat:
                # integrate (smooth) the chunks
                if histogram_type != None:
                    chunk_tth, chunk_intensity, chunk_azm = hist.histogram1d(
                        chunk_tth,
                        chunk_intensity,
                        azi=chunk_azm,
                        histogram_type=histogram_type,
                        bin_n=histogram_bins,
                        debug=debug,
                    )

                # append azimuth to output
                new_azi_chunks.append(azichunks[j])

                # Background estimates
                background_guess = get_chunk_background_guess(
                    settings_as_class, (chunk_intensity, chunk_tth), n, debug=debug
                )

                cent_guess, h_guess, w_guess, p_guess = get_chunk_peak_guesses_new(
                    settings_as_class,
                    (chunk_intensity, chunk_tth, chunk_azm),
                    background_guess,
                    n,
                    w_guess_fraction,
                    dfour,
                    debug,
                )
                # cent_guess is in collected units. needs converting.
                lims = lmm.parse_bounds(
                    settings_as_class.fit_bounds,
                    data_as_class,
                    0,
                    len(cent_guess),
                    param=["d-space"],
                )  # limits alwasys in converted unit.
                limits_range = lims["d-space"]

                if np.max(cent_guess) > np.max(limits_range) or np.min(
                    cent_guess
                ) < np.min(limits_range):
                    # then the peak_pos_guesses are not in d-spacing
                    # therefore change
                    cent_guess = list(
                        np.atleast_1d(
                            data_as_class.conversion(
                                cent_guess, azm=None, reverse=False
                            )
                        )
                    )
                # convert to dictionary
                peaks = []
                for a, b, c, d in zip(cent_guess, h_guess, w_guess, p_guess):
                    peaks.append(
                        {
                            "d-space": [a],
                            "height": [b],
                            "width": [c],
                            "profile": [d],
                        }
                    )

                comp_list = ["h", "d", "w", "p"]
                comp_names = ["height", "d-space", "width", "profile"]

                local_orders = deepcopy(settings_as_class.subfit_orders)
                for b in range(len(local_orders["background"])):
                    local_orders["background"][b] = 0
                for pk in range(len(peaks)):
                    for i in range(len(comp_names)):
                        val = local_orders["peak"][pk][comp_names[i]]
                        # print(val)
                        if isinstance(val, int):
                            # has to be single value with no constraints.
                            if val > 0:
                                local_orders["peak"][pk][comp_names[i]] = 0

                        elif (
                            len(
                                sc.SeriesConstraints(
                                    local_orders["peak"][pk][comp_names[i]]
                                )
                            )
                            > 0
                        ):
                            # if length >0 there has to be constrains on the peak parameters

                            local_orders["peak"][pk][comp_names[i]] = [
                                0
                            ] + sc.SeriesConstraints(
                                local_orders["peak"][pk][comp_names[i]]
                            )

                        else:
                            # the orders a list of numbers - in which case stil
                            # have to be reduced to 0
                            local_orders["peak"][pk][comp_names[i]] = 0

                        if comp_names[i] + "_fixed" in local_orders["peak"][pk]:
                            # fixed = 1
                            local_orders["peak"][pk][comp_names[i] + "_fixed"] = 1

                params = lmm.initiate_all_params_for_fit(
                    settings_as_class,
                    chunk_data,
                    peak_orders=local_orders,
                    values={"background": background_guess, "peak": peaks},
                    types=False,
                )

                # unvary fixed parameters
                for pk in range(len(peaks)):
                    for i in range(len(comp_names)):
                        if comp_names[i] + "_fixed" in local_orders["peak"][pk]:
                            params = lmm.un_vary_single_param(
                                params,
                                f"peak_{pk}",
                                comp_list[i],
                                0,
                            )

                logger.debug(
                    " ".join(
                        map(str, [("Initiallised chunk; %i/%i" % (j, len(chunks)))])
                    )
                )
                logger.log_lmfit_obj(params, level="DEBUG", space=True)
                # if debug:
                if logger.is_below_level(level="DEBUG"):
                    # keep the original params for plotting afterwards
                    guess = params

                # Run actual fit
                fit = lmm.fit_model(
                    chunk_data,  # needs to contain intensity, tth, azi (as chunks), conversion factor
                    local_orders,
                    params,
                    fit_method=fit_method,
                    max_n_fev=max_n_f_eval,
                    weights=None,
                )
                params = fit.params  # update lmfit parameters

                # if debug:
                logger.debug(
                    " ".join(map(str, [("Fitted chunk; %i/%i" % (j, len(chunks)))]))
                )
                logger.log_lmfit_obj(params, level="DEBUG", space=True)

                # get values from fit and append to arrays for output
                for i in range(len(background_guess)):
                    out_vals["bg"][i].append(params["bg_c" + str(i) + "_f0"].value)
                    # catch None is errors
                    if params["bg_c" + str(i) + "_f0"].stderr != None:
                        out_vals["bg_err"][i].append(
                            params["bg_c" + str(i) + "_f0"].stderr
                        )
                    else:
                        out_vals["bg_err"][i].append(np.nan)
                comp_list = ["h", "d", "w", "p"]
                comp_names = ["height", "d-space", "width", "profile"]
                for pk in range(len(peaks)):
                    for cp in range(len(comp_list)):
                        out_vals[comp_list[cp]][pk].append(
                            params["peak_" + str(pk) + "_" + comp_list[cp] + "0"].value
                        )
                        # catch None is errors
                        if (
                            params["peak_" + str(pk) + "_" + comp_list[cp] + "0"].stderr
                            != None
                        ):
                            out_vals[comp_list[cp] + "_err"][pk].append(
                                params[
                                    "peak_" + str(pk) + "_" + comp_list[cp] + "0"
                                ].stderr
                            )
                        else:
                            out_vals[comp_list[cp] + "_err"][pk].append(np.nan)
                out_vals["chunks"].append(azichunks[j])

                # plot the fits.
                if logger.is_below_level(level="DEBUG"):
                    tth_plot = chunk_data.tth
                    int_plot = chunk_data.intensity
                    azm_plot = chunk_data.azm
                    # required for plotting energy dispersive data. - but I am not sure why
                    azm_plot = ma.array(azm_plot, mask=(~np.isfinite(azm_plot)))
                    gmodel = Model(
                        lmm.peaks_model,
                        independent_vars=["two_theta", "azimuth"],
                        data_class=chunk_data,
                        orders=settings_as_class.subfit_orders,
                    )
                    tth_model = np.linspace(
                        np.min(chunk_data.tth),
                        np.max(chunk_data.tth),
                        100,
                    )
                    azm_model = tth_model * 0 + np.mean(azm_plot)

                    mod_plot = gmodel.eval(
                        params=params,
                        two_theta=tth_model,
                        azimuth=azm_model,
                    )
                    guess_plot = gmodel.eval(
                        params=guess,
                        two_theta=tth_model,
                        azimuth=azm_model,
                    )
                    plt.plot(tth_plot, int_plot, ".", label="data")
                    plt.plot(
                        tth_model,
                        np.array(guess_plot).flatten(),
                        marker="",
                        color="green",
                        linewidth=2,
                        linestyle="dashed",
                        label="guess",
                    )
                    plt.plot(
                        tth_model,
                        np.array(mod_plot).flatten(),
                        marker="",
                        color="red",
                        linewidth=2,
                        label="fit",
                    )
                    plt.legend()
                    plt.title(
                        (
                            (
                                peak_string(settings_as_class.subfit_orders)
                                + "; azimuth = %.1f"
                            )
                            % azichunks[j]
                        )
                    )
                    plt.show()
            else:
                # there are not enough data to fit for the chunks
                pass

        else:
            raise ValueError("The mode for processing the chunks is not recognised.")

    return out_vals, new_azi_chunks


def fit_series(
    master_params,
    data,
    settings_as_class,
    start_end=[0, 360],
    debug=False,
    save_fit=False,
    max_n_f_eval=400,
):
    orders = settings_as_class.subfit_orders

    # Initiate background parameters and perform an initial fit
    comp = "bg"
    for b in range(len(orders["background"])):
        param_str = "bg_c" + str(b)
        comp = "f"

        logger.log_lmfit_obj(master_params, level="DEBUG", space=True)
        master_params = lmm.un_vary_params(
            master_params, param_str, comp
        )  # set other parameters to not vary
        master_params = lmm.vary_params(
            master_params, param_str, comp
        )  # set these parameters to vary
        if isinstance(
            orders["background"][b], list
        ):  # set part of these parameters to not vary
            master_params = lmm.un_vary_part_params(
                master_params, param_str, comp, orders["background"][b]
            )

        azimuth = data[1]
        data_vals = data[0]["bg"][b]
        data_val_errors = data[0]["bg_err"][b]
        data_val_errors = clean_errs(data_val_errors)

        fout = lmm.coefficient_fit(
            azimuth=azimuth,
            ydata=data_vals,
            inp_param=master_params,
            param_str=param_str + "_" + comp,
            symmetry=1,
            errs=np.array(data_val_errors),
            fit_method="leastsq",
        )
        # fout.plot(show_init=True)

        master_params = fout.params

    # initiate peak(s)
    for j in range(len(orders["peak"])):
        # defines peak string to start parameter name
        param_str = "peak_" + str(j)

        comp_list = ["h", "d", "w", "p"]
        comp_names = ["height", "d-space", "width", "profile"]

        for cp in range(len(comp_list)):
            comp = comp_list[cp]
            if comp == "d":
                symmetry = 1
            if "symmetry" in orders["peak"][j]:
                symmetry = orders["peak"][j]["symmetry"]
            else:
                symmetry = 1

            if comp_names[cp] + "_fixed" in orders["peak"][j]:
                fixed = 1
                # data_vals are not needed
                # data_val_errors are not needed.
            else:
                fixed = 0
                data_vals = data[0][comp][j]
                data_val_errors = data[0][comp + "_err"][j]
                data_val_errors = clean_errs(data_val_errors, outliers=2)

            #     # FIX ME: this was not checked properly.the values it feeds are not necessarily correct
            #     # and the fixed parameters might be fit for.
            # coeff_type = pf.get_params_type(orders, comp, peak=j)
            n_coeff = sf.get_number_coeff(orders, comp, peak=j, azimuths=data[1])
            master_params = lmm.un_vary_params(
                master_params, param_str, comp
            )  # set other parameters to not vary
            if fixed == 0:
                master_params = lmm.vary_params(
                    master_params, param_str, comp
                )  # set these parameters to vary
                if isinstance(
                    sc.SeriesValues(orders["peak"][j][comp_names[cp]]), list
                ):  # set part of these parameters to not vary
                    master_params = lmm.un_vary_part_params(
                        master_params,
                        param_str,
                        comp,
                        orders["peak"][j][comp_names[cp]],
                    )
                # catch to make sure there are not more coefficients than chunks
                if n_coeff > len(np.unique(data[1])):
                    o = int(np.floor(len(np.unique(data[1])) / 2 - 1))
                    un_vary = [x for x in range(0, o)]
                    master_params = lmm.un_vary_part_params(
                        master_params, param_str, comp, un_vary
                    )
                fout = lmm.coefficient_fit(
                    azimuth=data[1],
                    ydata=data_vals,
                    inp_param=master_params,
                    param_str=param_str + "_" + comp,
                    symmetry=symmetry,
                    errs=np.array(data_val_errors),
                    fit_method="leastsq",
                    max_nfev=max_n_f_eval,
                )
                # fout.plot(show_init=True)

            master_params = fout.params

    logger.debug(" ".join(map(str, [("Parameters after initial Fourier fits")])))
    logger.log_lmfit_obj(master_params, level="DEBUG", space=True)

    # plot output of series fits....
    if logger.is_below_level("DEBUG") is True:
        x_lims = start_end
        azi_plot = range(np.int_(x_lims[0]), np.int_(x_lims[1]), 2)
        gmodel = Model(sf.coefficient_expand, independent_vars=["azimuth"])

        comp_list = ["h", "d", "w", "p"]
        comp_names = ["height", "d-space", "width", "profile"]

        # loop over peak parameters and plot.
        fig = plt.figure()
        ax = []
        for i in range(len(comp_list)):
            ax.append(fig.add_subplot(5, 1, i + 1))
            ax[i].set_title(comp_names[i])
            ax[i].set_xlim(x_lims)
            for j in range(len(orders["peak"])):
                param_str = "peak_" + str(j)
                comp = comp_list[i]
                if "symmetry" in orders["peak"][j].keys() and comp != "d":
                    symm = orders["peak"][j]["symmetry"]
                else:
                    symm = 1
                temp = lmm.gather_param_errs_to_list(master_params, param_str, comp)
                temp_tp = sf.get_series_type(master_params, param_str, comp)
                if temp_tp == sf.coefficient_types()["independent"]:
                    az_plt = np.array(data[1])
                else:
                    az_plt = azi_plot
                gmod_plot = gmodel.eval(
                    params=master_params,
                    azimuth=np.array(az_plt) * symm,
                    param=temp[0],
                    coeff_type=temp_tp,
                )
                ax[i].plot(az_plt, gmod_plot)
                ax[i].scatter(data[1], data[0][comp][j], s=10)

            # set y limits ignoring the error bars
            y_lms = ax[i].get_ylim()
            ax[i].set_ylim(y_lms)

            for j in range(len(orders["peak"])):
                # add error bars to data points, over the top of the points
                ax[i].errorbar(
                    data[1],
                    data[0][comp][j],
                    yerr=data[0][comp + "_err"][j],
                    fmt="none",
                    elinewidth=0.5,
                )

            # set x-ticks by data type.
            # if notthing in the data class then continue
            try:
                label_x = settings_as_class.data_class.dispersion_ticks(
                    disp_ticks=ax[i + 1].get_xticks
                )
                ax[i + 1].set_xticks(label_x)
            except:
                pass

        # plot background
        for k in range(len(orders["background"])):
            x_plt = len(orders["background"])
            y_plt = len(orders["background"]) * 4 + k + 1
            ax.append(fig.add_subplot(5, x_plt, y_plt))
            ax[i + k + 1].set_title("Background " + str(k))

            param_str = "bg_c" + str(k)
            comp = "f"
            temp = lmm.gather_param_errs_to_list(master_params, param_str, comp)
            temp_tp = sf.get_series_type(master_params, param_str, comp)
            if temp_tp == sf.coefficient_types()["independent"]:
                az_plt = np.array(data[1])
            else:
                az_plt = azi_plot
            gmod_plot = gmodel.eval(
                params=master_params,
                azimuth=np.array(az_plt),
                param=temp[0],
                coeff_type=temp_tp,
            )
            ax[i + k + 1].scatter(data[1], data[0]["bg"][k], s=10)
            y_lms = ax[i + k + 1].get_ylim()
            # add error bars to data points, over the top of the points
            ax[i + k + 1].errorbar(
                data[1],
                data[0]["bg"][k],
                yerr=data[0]["bg_err"][k],
                fmt="none",
                elinewidth=0.5,
            )
            ax[i + k + 1].plot(
                az_plt,
                gmod_plot,
            )
            # set y limits ignoring the error bars
            ax[i + k + 1].set_ylim(y_lms)
            ax[i + k + 1].set_xlim(x_lims)
            # set x-labels by data type.
            # if notthing in the data class then continue
            try:
                label_x = settings_as_class.data_class.dispersion_ticks(
                    disp_ticks=ax[i + k + 1].get_xticks
                )
                ax[i + k + 1].set_xticks(label_x)
            except:
                pass

        fig.suptitle(peak_string(orders) + "; Fits to Chunks")

        if save_fit:
            filename = make_outfile_name(
                settings_as_class.subfit_filename,
                directory=settings_as_class.output_directory,
                additional_text="ChunksFit",
                orders=orders,
                extension=".png",
                overwrite=False,
            )
            fig.savefig(filename)
        plt.show()
        plt.close()

    return master_params


def clean_errs(error_values, outliers=5):
    """
    The chunk fits can produce errors that are either very small (<1e-10) or
    very large (>1e5). These massively skew the series fits. This functions therefore
    idientifies errors the are more than 'outliers' * stddev(error_values) from the
    median and replaces them with a value that equals median_error +/- 'outliers' * stddev(error_values)

    nan errors are accounted for.

    Parameters
    ----------
    error_values : list, np.array
        Array or list of error values.
    outliers : float, optional
        Value to filter the outliers, filters values more than this number of times standard deviaiton of values away from median.
        The default is 5.

    Returns
    -------
    error_values : np.array
        Filtered array of error values..

    """
    # N.B. (SAH) I dont know if this function is still necessary now that I have fixed the weight/errors into Model.fit.
    # [the errors needed to be square rooted before sending ... to give sesible values when calling fout.plot(show_init=True)

    error_values = np.array(error_values)

    err_stdev = np.nanmedian(error_values)
    log_err_vals = np.log10(error_values)

    dlog = np.abs(log_err_vals - np.nanmedian(log_err_vals))
    stdev_dlog = np.nanmedian(dlog)
    s = dlog / (stdev_dlog if stdev_dlog else 1.0)

    error_values[s > outliers] = np.nanmedian(error_values) + outliers * err_stdev

    # remove nan values from error_values
    error_values[np.isnan(error_values)] = (
        np.nanmedian(error_values) + outliers * err_stdev * 10
    )

    return error_values
