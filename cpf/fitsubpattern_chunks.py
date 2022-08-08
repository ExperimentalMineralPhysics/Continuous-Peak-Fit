#!/usr/bin/env python

__all__ = ["get_manual_guesses", "get_chunk_background_guess", "fit_chunks"]

# CPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from lmfit import Parameters, Model
import cpf.IO_functions as io
import cpf.series_functions as sf
import cpf.lmfit_model as lmm


def get_manual_guesses(settings_as_class, data_as_class, debug=False):
    """ 
    :param data_as_class:
    :param settings_as_class:
    :param debug:
    :return dfour:
    """
    
    peeks = len(settings_as_class.subfit_orders["peak"])
    t_th_guesses = np.array(settings_as_class.subfit_orders["PeakPositionSelection"])

    settings_as_class.validate_position_selection(peak_set=settings_as_class.subfit_order_position, report=False)

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
        lims = lmm.parse_bounds(
            settings_as_class.fit_bounds, data_as_class, 0, 0, param=["d-space"]
        )
        # get coefficient type
        coeff_type = sf.params_get_type(settings_as_class.subfit_orders, comp, peak=j)

        # for guesses make sure there are not too many coefficients.
        n_coeff = sf.get_number_coeff(settings_as_class.subfit_orders, comp)
        
        if n_coeff > len(t_th_guess[:, 0]):
            o = int(np.floor(len(t_th_guess[:, 0]) / 2 - 1))
            # FIX ME: this should be a function!! and should check the series type.
        else:
            o = settings_as_class.subfit_orders["peak"][j]["d-space"]

        temp_param = Parameters()
        temp_param = lmm.initiate_params(
            inp_param=temp_param,
            param_str=param_str,
            coeff_type=coeff_type,
            comp=comp,
            trig_orders=o,
            limits=lims["d-space"],
            value=data_as_class.conversion(t_th_guess[0, 1], azm=t_th_guess[0, 1]),
        )
        fout = lmm.coefficient_fit(
            azimuth=t_th_guess[:,0],
            ydata=data_as_class.conversion(t_th_guess[:, 1], azm=t_th_guess[:, 0]),
            inp_param=temp_param,
            param_str=param_str + "_" + comp,
            fit_method="leastsq",
        )
        temp_param = fout.params
        if debug:
            temp_param.pretty_print()
        dfour.append(
            lmm.gather_param_errs_to_list(temp_param, param_str, comp)
        )
    return dfour


def get_chunk_background_guess(settings_as_class, data_chunk_class, n, debug=False):
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
    # Get indices of sorted two theta values excluding the masked values
    tth_ord = ma.argsort(data_chunk_class.tth.compressed())
    
    background_guess = [[0.0] for i in range(len(settings_as_class.subfit_orders["background"]))]
    
    if "background_fixed" in settings_as_class.fit_orders:
        # this is not yet implemented but the option exists here where it is needed first.
        raise NotImplementedError
        
    else:
        if len(settings_as_class.subfit_orders["background"]) == 1:
            # assume background is present either at the left or right edges of the region.
            # FIXME: this assumes a positive peak. It may not work for absorption peaks. Not tested though.
            background_guess[0][0] = np.min(
                [
                    np.mean(
                        data_chunk_class.intensity.compressed()[
                            tth_ord[:n]
                        ]
                    ),
                    np.mean(
                        data_chunk_class.intensity.compressed()[
                            tth_ord[-n:]
                        ]
                    ),
                ]
            )
        else:  # len(orders['background']) > 1:
            # first value (offset) is mean of left-hand values
            background_guess[0][0] = np.mean(
                data_chunk_class.intensity.compressed()[
                    tth_ord[:n]
                ]
            )
            # if there are more, then calculate a gradient guess.
            background_guess[1][0] = (
                                              np.mean(
                                                  data_chunk_class.intensity.compressed()[
                                                      tth_ord[-n:]
                                                  ]
                                              )
                                              - np.mean(
                                                  data_chunk_class.intensity.compressed()[
                                                    tth_ord[:n]
                                                  ]
                                              )
                                      ) / (
                                              data_chunk_class.tth.compressed()[
                                                  tth_ord[-1]
                                              ]
                                              - data_chunk_class.tth.compressed()[
                                                  tth_ord[0]
                                              ]
                                      )
            # leave all higher order terms as 0 -- assume they are small.
    
    return background_guess


def get_chunk_peak_guesses(settings_as_class, data_chunk_class, 
                                                        background_guess,
                                                        n, w_guess_fraction,
                                                        dfour, debug=False):
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
    d_guess = []
    h_guess = []
    w_guess = []
    p_guess = []
    p_fixed = 0
    peaks = []
    lims = []

    for k in range(len(settings_as_class.subfit_orders["peak"])):
        # FIX ME: altered from locals call as now passed as None - check!
        if dfour: # If the positions have been pre-guessed extract values from dfour.
            if debug and k == 0:
                print("dfour in locals \n")
            
            coeff_type = sf.params_get_type(settings_as_class.subfit_orders, "d", peak=k)
            d_guess = sf.coefficient_expand(
                np.mean(data_chunk_class.azm),
                dfour[k][0],
                coeff_type=coeff_type,
            )
            t_th_guess = data_chunk_class.conversion(
                d_guess,
                azm=np.mean(data_chunk_class.azm),reverse=True
            )
            # Finds the index of the closest two-theta to the d-spacing input
            idx = (
                np.abs(data_chunk_class.tth - t_th_guess)
            ).argmin()
            # FIX ME: The mean of a number of the smallest values would be more stable.

            # Height is intensity of the closest pixel in d-spacing - background at that position
            h_guess = data_chunk_class.intensity[idx]
            if len(background_guess) > 1:
                h_guess = h_guess - (
                        background_guess[0][0]
                        + background_guess[1][0]
                        * (
                                data_chunk_class.tth[idx]
                                - data_chunk_class.tth.min()
                        )
                )
            elif len(background_guess) == 1:
                h_guess = h_guess - background_guess[0][0]

        else:  # 'guess' guesses from data slice - assuming a single peak
            # find the brightest pixels
            # get index of nth brightest pixel.
            idx = np.argsort(
                data_chunk_class.intensity.compressed()
            )[-n:][0]

            # height guess is nth highest intensity - background guess at this position
            h_guess = (data_chunk_class.intensity.compressed())[idx]
            if len(background_guess) > 1:
                h_guess = h_guess - (
                        background_guess[0][0]
                        + background_guess[1][0]
                        * (
                                data_chunk_class.tth.compressed()[
                                    idx
                                ]
                                - data_chunk_class.tth.min()
                        )
                )
            elif len(background_guess) == 1:
                h_guess = h_guess - background_guess[0][0]

            # d-spacing of the highest nth intensity pixel
            d_guess = data_chunk_class.dspace.compressed()[idx]

        # w_guess is fractional width of the data range
        w_guess = (
                (
                        np.max(data_chunk_class.tth)
                        - np.min(data_chunk_class.tth)
                )
                / w_guess_fraction
                / len(settings_as_class.subfit_orders["peak"])
        )
        # FIX ME: This is a bit crude. Is there a better way to do it?

        # If profile_fixed exists then set p_guess to profile_fixed values and set p_fixed to 1
        p_guess = 0.5  # Guess half if solving for profile
        p_fixed = (
            0  # Set to 0 so solving unless profile_fixed exists.
        )
        if "profile_fixed" in settings_as_class.subfit_orders["peak"][k]:
            # Profile fixed is the list of coefficients if the profile is fixed.
            # FIX ME: profile fixed must have the same number of coefficients as required by
            # profile.
            
            coeff_type = sf.params_get_type(settings_as_class.subfit_orders, "p", peak=k)
            if "symmetry" in settings_as_class.subfit_orders["peak"][k]:
                symm = settings_as_class.subfit_orders["peak"][k]["symmetry"]
            else:
                symm = 1
            p_guess = sf.coefficient_expand(
                np.mean(data_chunk_class.azm) * symm,
                param=settings_as_class.subfit_orders["peak"][k]["profile_fixed"],
                coeff_type=coeff_type,
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
            lmm.parse_bounds(
                settings_as_class.fit_bounds,
                data_chunk_class,
                param=["height", "d-space", "width", "profile"],
            )
        )
        data_chunk_class

    limits = {
        "background": lmm.parse_bounds(
                settings_as_class.fit_bounds,
                data_chunk_class,
            param=["background"],
        )["background"],
        "peak": lims,
    }
    return peaks, limits, p_fixed



def fit_chunks(
    data_as_class,
    settings_as_class,
    save_fit=False,
    debug=False,
    fit_method=None,
    mode="fit",
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
        print(
            "\nNo previous fits or selections. Assuming "
            + str(peeks)
            + " "
            + p
            + " and group fitting in azimuth...\n"
        )
        dfour = None
        # FIXME: passed as None now - check functions as expected
    else:
        print("\nUsing manual selections for initial guesses...\n")
        # FIXME: Need to check that the num. of peaks for which we have parameters is the same as the
        # number of peaks guessed at.
        #dfour = get_manual_guesses(peeks, orders, bounds, twotheta, debug=None)
        dfour = get_manual_guesses(settings_as_class, data_as_class, debug=debug)
                
    # Get chunks according to detector type.
    if mode != "fit": #cascade
        chunks, azichunks = data_as_class.bins(settings_as_class, cascade=True)
    else: 
        chunks, azichunks = data_as_class.bins(settings_as_class)
    # Final output list of azimuths with corresponding twotheta_0,h,w
    
    # setup arrays
    new_azi_chunks = []
    if mode == "fit" or mode=="cascade": 
        #bunch of settings
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
    
    out_vals["peak"] = io.peak_string(settings_as_class.subfit_orders)
    
    for j in range(len(chunks)):
        # print('\nFitting to data chunk ' + str(j + 1) + ' of ' + str(len(chunks)) + '\n')
        if debug:
            print("\n")
        print(
            "\rFitting to data chunk %s of %s "
            % (str(j + 1), str(len(chunks))),
            end="\r",
        )
        
        #make data class for chunks.
        chunk_data = data_as_class.duplicate()
        #reduce data to a subset
        # FIXME: this is crude but I am not convinced that it needs to be contained within the data class.
        # FEXME: maybe I need to reconstruct the data class so that dat_class.tth is a function that applies a mask when called. but this will be slower.
        chunk_data.intensity = chunk_data.intensity.flatten()[chunks[j]]
        chunk_data.tth       = chunk_data.tth.flatten()[chunks[j]]
        chunk_data.azm       = chunk_data.azm.flatten()[chunks[j]]
        chunk_data.dspace    = chunk_data.dspace.flatten()[chunks[j]]
        
        
        # find other output from intensities
        if mode == "maxima":
            # get maximum from each chunk
            out_vals["h"][0].append(np.max(chunk_data.intensity))
            out_vals["chunks"].append(azichunks[j])
            new_azi_chunks.append(azichunks[j])
            
            
        elif mode == "range":
            # get maximum from each chunk
            out_vals["h"][0].append(np.max(chunk_data.intensity)-np.min(chunk_data.intensity))
            out_vals["chunks"].append(azichunks[j])
            new_azi_chunks.append(azichunks[j])
            
        elif mode == "98thpercentile":
            # FIXME: this is crude and could be done better -- i.e. with a switch for the percentile
            # get 98th percentils from each chunk
            raise NotImplementedError
            
        elif mode == "fit" or mode=="cascade":
            
            # Define parameters to pass to fit
            params = Parameters()
            
            if ma.MaskedArray.count(chunk_data.intensity) >= min_dat:
        
                # append azimuth to output
                new_azi_chunks.append(azichunks[j])
        
                # Background estimates
                background_guess = get_chunk_background_guess(settings_as_class, chunk_data, n, debug=debug)
        
                # # Organise guesses to be refined.
                peaks, limits, p_fixed = get_chunk_peak_guesses(settings_as_class, chunk_data, 
                                                            background_guess,
                                                            n, w_guess_fraction, dfour, debug)
        
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
                comp_list = ["h", "d", "w", "p"]
                comp_names = ["height", "d-space", "width", "profile"]
                for pk in range(len(peaks)):
                    for cp in range(len(comp_list)):
                        if comp_names[cp]+"_fixed" in settings_as_class.subfit_orders["peak"][pk]:
                            vary=False
                        else:
                            vary=True
                        params.add(
                            "peak_" + str(pk) + "_"+comp_list[cp]+"0",
                            peaks[pk][comp_names[cp]][0],
                            min=limits["peak"][pk][comp_names[cp]][0],
                            max=limits["peak"][pk][comp_names[cp]][1],
                            vary=vary
                        )
    
                if debug:
                    print("Initiallised chunk fit; "+str(j)+"/"+str(len(chunks)))
                    params.pretty_print()
                    print("\n")
                    # keep the original params for plotting afterwards
                    guess = params
                    
                # Run actual fit
                fit = lmm.fit_model(
                    chunk_data, # needs to contain intensity, tth, azi (as chunks), conversion factor
                    settings_as_class.subfit_orders,
                    params,
                    fit_method=fit_method,
                    weights=None,
                )  
                params = fit.params  # update lmfit parameters
                
                if debug:
                    print("Final chunk fit; "+str(j)+"/"+str(len(chunks)))
                    params.pretty_print()
                
                # get values from fit and append to arrays for output
                for i in range(len(background_guess)):
                    out_vals["bg"][i].append(params["bg_c" + str(i) + "_f0"].value)
                    out_vals["bg_err"][i].append(
                        params["bg_c" + str(i) + "_f0"].stderr
                    )
                comp_list = ["h", "d", "w", "p"]
                comp_names = ["height", "d-space", "width", "profile"]
                for pk in range(len(peaks)):
                    for cp in range(len(comp_list)):
                        out_vals[comp_list[cp]][pk].append(params["peak_"+str(pk)+"_"+comp_list[cp]+"0"].value)
                        out_vals[comp_list[cp]+"_err"][pk].append(
                            params["peak_" + str(pk) + "_"+comp_list[cp]+"0"].stderr
                        )
                
                out_vals["chunks"].append(azichunks[j])
            
                # plot the fits.
                if debug:
                    tth_plot = chunk_data.tth
                    int_plot = chunk_data.intensity
                    azm_plot = chunk_data.azm
                    # required for plotting energy dispersive data. - but I am not sure why
                    azm_plot = ma.array(azm_plot, mask=(~np.isfinite(azm_plot)))
                    gmodel = Model(
                        lmm.peaks_model, independent_vars=["two_theta", "azimuth"], data_class = chunk_data,
                                orders = settings_as_class.subfit_orders, 
                    )
                    tth_range = np.linspace(
                        np.min(chunk_data.tth),
                        np.max(chunk_data.tth),
                        azm_plot.size,
                    )
                    mod_plot = gmodel.eval(
                        params=params,
                        two_theta=tth_range,
                        azimuth=azm_plot,
                    )
                    guess_plot = gmodel.eval(
                        params=guess,
                        two_theta=tth_range,
                        azimuth=azm_plot,
                    )
                    plt.plot(tth_plot, int_plot, ".", label="data")
                    plt.plot(
                        tth_range,
                        np.array(guess_plot).flatten(),
                        marker="",
                        color="green",
                        linewidth=2,
                        linestyle="dashed",
                        label="guess",
                    )
                    plt.plot(
                        tth_range,
                        np.array(mod_plot).flatten(),
                        marker="",
                        color="red",
                        linewidth=2,
                        label="fit",
                    )
                    plt.legend()
                    plt.title(((io.peak_string(settings_as_class.subfit_orders) + "; azimuth = %.1f" ) % azichunks[j]))
                    plt.show()
            else:
                #there are not enough data to fit for the chunks
                pass
        
        else:
            raise ValueError("The mode for processing the chunks is not recognised.")
            
    return out_vals, new_azi_chunks
    

def fit_series(master_params, data, settings_as_class, debug=False, save_fit=False):
    
    
    orders = settings_as_class.subfit_orders
    
    # Initiate background parameters and perform an initial fit
    comp = "bg"
    for b in range(len(orders["background"])):
        param_str = "bg_c" + str(b)
        comp = "f"
        
        master_params.pretty_print()
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
                            
        fout = lmm.coefficient_fit(
            azimuth=azimuth,
            ydata=data_vals,
            inp_param=master_params,
            param_str=param_str + "_" + comp,
            symmetry=1,
            errs=data_val_errors,
            fit_method="leastsq",
        )
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
            else:
                symmetry = orders["peak"][j]["symmetry"]
                
                
            if comp_names[cp] + "_fixed" in orders["peak"][j]:
                fixed = 1
                # data_vals are not needed
                # data_val_errors are not needed.
            else:
                fixed = 0
                data_vals = data[0][comp][j]
                data_val_errors = data[0][comp][j]
                
            #     # FIX ME: this was not checked properly.the values it feeds are not necessarily correct
            #     # and the fixed parameters might be fit for.
            # coeff_type = pf.params_get_type(orders, comp, peak=j)
            n_coeff = sf.get_number_coeff(
                orders, comp, peak=j, azimuths=data[1]
            )
            master_params = lmm.un_vary_params(
                master_params, param_str, comp
            )  # set other parameters to not vary
            if fixed == 0:
                master_params = lmm.vary_params(
                    master_params, param_str, comp
                )  # set these parameters to vary
                if isinstance(
                    orders["peak"][j][comp_names[cp]], list
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
                    master_params = lmm.un_vary_part_params(master_params, param_str, comp, un_vary)
                fout = lmm.coefficient_fit(
                    azimuth=data[1],
                    ydata=data_vals,
                    inp_param=master_params,
                    param_str=param_str + "_" + comp,
                    symmetry=symmetry,
                    errs=data_val_errors,
                    fit_method="leastsq",
                )
    
            master_params = fout.params
            
            # FIX ME. Check which params should be varying and which should not.
            # Need to incorporate vary and un-vary params as well as partial vary
    
    if 1:#debug:
        print("Parameters after initial Fourier fits")
        master_params.pretty_print()
        
        # plot output of fourier fits....
        x_lims = np.array([np.min(data[1]), np.max(data[1])])
        x_lims = np.around(x_lims / 90) * 90
        x_ticks = list(range(int(x_lims[0]),int(x_lims[1]+1),45))
        azi_plot = range(np.int(x_lims[0]), np.int(x_lims[1]), 2)
        gmodel = Model(sf.coefficient_expand, independent_vars=["azimuth"])
    
        comp_list = ["h", "d", "w", "p"]
        comp_names = ["height", "d-space", "width", "profile"]
        
        # loop over peak parameters and plot.
        fig = plt.figure()
        ax=[]
        for i in range(len(comp_list)):
            ax.append(fig.add_subplot(5, 1, i+1))
            ax[i].set_title(comp_names[i])
            ax[i].set_xlim(x_lims)
            for j in range(len(orders["peak"])):
            
                param_str = "peak_" + str(j)
                comp = comp_list[i]
                if "symmetry" in orders["peak"][j].keys() and comp != "d":
                    symm = orders["peak"][j]["symmetry"]
                else:
                    symm = 1
                temp = lmm.gather_param_errs_to_list(
                    master_params, param_str, comp
                )
                temp_tp = sf.get_series_type(master_params, param_str, comp)
                if temp_tp == 5:
                    az_plt = np.array(data[1])
                else:
                    az_plt = azi_plot
                gmod_plot = gmodel.eval(
                        params=master_params, 
                        azimuth=np.array(az_plt)*symm, 
                        param=temp[0], 
                        coeff_type=temp_tp
                )
                ax[i].scatter(data[1], data[0][comp][j], s=10)
                ax[i].plot(az_plt, gmod_plot, )
                # set x-labels by data type.
                label_x = settings_as_class.data_class.dispersion_ticks(disp_ticks=ax[i].get_xticks)
                ax[i].set_xticks(label_x)
    
        #plot background 
        for k in range(len(orders["background"])):
            x_plt = len(orders["background"])
            y_plt = len(orders["background"]) * 4 + k +1
            ax.append(fig.add_subplot(5, x_plt, y_plt))
            ax[i+k+1].set_title("Background "+str(k))
    
            param_str = "bg_c" + str(k)
            comp = "f"
            temp = lmm.gather_param_errs_to_list(
                master_params, param_str, comp
            )
            temp_tp = sf.get_series_type(master_params, param_str, comp)
            if temp_tp == 5:
                az_plt = np.array(data[1])
            else:
                az_plt = azi_plot
            gmod_plot = gmodel.eval(
                    params=master_params, 
                    azimuth=np.array(az_plt),
                    param=temp[0], 
                    coeff_type=temp_tp
            )
            ax[i+k+1].scatter(data[1], data[0]["bg"][k], s=10)
            ax[i+k+1].plot(az_plt, gmod_plot, )
            ax[i+k+1].set_xlim(x_lims)
            # set x-labels by data type.
            label_x = settings_as_class.data_class.dispersion_ticks(disp_ticks=ax[i+k+1].get_xticks)
            ax[i+k+1].set_xticks(label_x)
            
        fig.suptitle(io.peak_string(orders) + "; Fits to Chunks")
    
        if save_fit:
            filename = io.make_outfile_name(
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