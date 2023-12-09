#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:16:58 2021

@author: simon


This files makes the cascade plots of azimuthal intensity vis time. 
It needs to do a few things.

1. Call XRD_FitPattern to make the azimual distribution files.
    a. by fitting a peak to the data 
    b. max
    c. mean of the 90th+ percentiles.
    
    To be consistent with the philosophy of the code, the bins for the aximuths should be determined by a combination of width and the number of data contained. By keeping the number of data constant we will have more points at large two theta than at small two theta. 
    - this will need a switch in the input file somewhere. (Azibins needs another call e.g. AziNum?) 

2. Plot the outputs as a cascade plot.

3. Determine the number of peaks in each steip and plot n vs. time. 
    Peaks need to be determinable by:
        a. above or below average (e.g. He 2000)
        b. peak finding algorithm (various)
        

"""

__all__ = ["initiate", "set_range", "execute"]


import numpy as np

# import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.colors as colours
import cpf.XRD_FitPattern as XRD_FitPattern
import cpf.IO_functions as IO
import os.path
import proglog
# import copy
# import pandas as pd
import json
from cpf.settings import settings
from cpf.XRD_FitSubpattern import fit_sub_pattern
from cpf.Cosmics import image_preprocess as cosmicsimage_preprocess
from cpf.BrightSpots import SpotProcess
from cpf.IO_functions import (
    json_numpy_serializer,
    make_outfile_name,
    peak_string,
    any_terms_null,
    title_file_names,
)

# from findpeaks import findpeaks

# from scipy import ndimage as ndi
# from skimage.feature import peak_local_max
# from skimage import data, img_as_float

import re

# import plotly.graph_objects as go
# import pandas as pd
from scipy.signal import find_peaks


def initiate(*args, **kwargs):
    """
    Run checks on input files, initiate data class and check output options
    :param report:
    :param out_type:
    :param initiate_data:
    :param setting_file:
    :param inputs:
    :return fit_parameters:
    :return fit_settings:
    """

    # pass everything through to XRD_FitPattern.initiate
    settings_for_fit = XRD_FitPattern.initiate(*args, **kwargs)
    
    return settings_for_fit


def set_range(*args, **kwargs):
    """
    :param setting_file:
    :param inputs:
    :param debug:
    :param refine:
    :param save_all:
    :param propagate:
    :param iterations:
    :param track:
    :param parallel:
    :param subpattern:
    :param kwargs:
    :return:
    """

    # pass everything through to XRD_FitPattern.set_range
    XRD_FitPattern.set_range(*args, **kwargs)
    


def execute(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    save_all=False,
    parallel=True,
    subpattern="all",
    mode="cascade",
    report=False,
    **kwargs
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param setting_file:
    :param parallel:
    :param report:
    :param mode:
    :param track:
    :param propagate:
    :param save_all:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """

    if setting_class is None:
        settings_for_fit = initiate(setting_file, inputs=inputs, report=True)
    else:
        settings_for_fit = setting_class
    new_data = settings_for_fit.data_class

    # Define locally required names
    temporary_data_file = make_outfile_name(
        "PreviousChunkFit_JSON",
        directory=settings_for_fit.output_directory,
        extension=".dat",
        overwrite=True,
    )

    if settings_for_fit.calibration_data:
        data_to_fill = os.path.abspath(settings_for_fit.calibration_data)
    else:
        # data_to_fill = os.path.abspath(settings_for_fit.datafile_list[0])
        data_to_fill = settings_for_fit.image_list[0]
    new_data.fill_data(
        data_to_fill,
        settings=settings_for_fit,
        debug=debug,
    )
    
    # Get calibration parameter file
    parms_dict = new_data.calibration
    # FIXME this should be removable.

    # plot calibration file
    if debug and settings_for_fit.calibration_data is not None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        new_data.plot_collected(fig_plot=fig, axis_plot=ax)
        plt.title("Calibration data")
        plt.show()
        plt.close()

    # if parallel processing start the pool
    # if parallel is True:
    #     p = mp.Pool(processes=mp.cpu_count())
    # p = mp.Pool()

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    # Process the diffraction patterns #
    for j in range(settings_for_fit.image_number):

        print("Process ", title_file_names(image_name=settings_for_fit.image_list[j]))
        # Get diffraction pattern to process.
        new_data.import_image(settings_for_fit.image_list[j], debug=debug)

        if settings_for_fit.datafile_preprocess is not None:
            # needed because image preprocessing adds to the mask and is different for each image.
            new_data.mask_restore()
            if "cosmics" in settings_for_fit.datafile_preprocess:
                new_data = cosmicsimage_preprocess(new_data, settings_for_fit)
        else:
            # nothing is done here.
            pass

        # plot input file
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity")
            plt.title(os.path.basename(settings_for_fit.datafile_list[j]))
            plt.show()
            plt.close()

        # Get previous fit (if it exists and is required)
        if (
            os.path.isfile(temporary_data_file)
            and settings_for_fit.fit_propagate is True
        ):  # and mode == "fit":
            # Read JSON data from file
            print("Loading previous fit results from %s" % temporary_data_file)
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)

        # Switch to save the first fit in each sequence.
        if j == 0 or save_all is True:
            save_figs = 1
        else:
            save_figs = 0

        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        all_fitted_chunks = []
        all_chunk_positions = []

        for i in range(len(settings_for_fit.fit_orders)):

            # get settings for current subpattern
            settings_for_fit.set_subpattern(j, i)

            if "previous_fit" in locals():  # and mode == "fit":
                params = previous_fit
            else:
                params = []

            # Track the position of the peak centroid
            # FIXME: This is crude - the range doesn't change width. so can't account for massive change in stress.
            # But does it need to?
            tth_range = settings_for_fit.subfit_orders["range"]
            if settings_for_fit.cascade_track is True and "previous_fit" in locals():
                clean = any_terms_null(params, val_to_find=None)
                if clean == 0:
                    # the previous fit has problems so discard it
                    print(
                        "Propagated fit has problems so not sesible to track the centre of the fit."
                    )
                else:
                    print("old subpattern range", tth_range)
                    mid = []
                    for k in range(len(params["peak"])):
                        mid.append(params["peak"][k]["d-space"][0])
                        # FIXME: replace with caluculation of mean d-spacing.

                    # Replace this with call through to class structure?
                    # print(mid)
                    cent = new_data.conversion(np.mean(mid), reverse=True)
                    # print("cent", cent)
                    # print("mean tth_range", np.mean(tth_range))
                    # print((tth_range[1] + tth_range[0]) / 2)
                    move_by = cent - np.mean(tth_range)
                    print("move by", move_by)
                    tth_range = tth_range + move_by
                    # tth_range[0] = tth_range[0] - ((tth_range[1] + tth_range[0]) / 2) + cent
                    # tth_range[1] = tth_range[1] - ((tth_range[1] + tth_range[0]) / 2) + cent
                    # print("move by:", ((tth_range[1] + tth_range[0]) / 2) - cent)
                    print("new subpattern range", tth_range)

                    # copy new positions back into settings_for_fit
                    settings_for_fit.fit_orders[i]["range"] = (
                        settings_for_fit.fit_orders[i]["range"] + move_by
                    )
                    print(settings_for_fit.fit_orders[i]["range"])

                    # The PeakPositionSelections are only used if the fits are not being propagated
                    if "PeakPositionSelection" in settings_for_fit.fit_orders[i]:
                        for k in range(
                            len(settings_for_fit.fit_orders[i]["PeakPositionSelection"])
                        ):
                            settings_for_fit.fit_orders[i]["PeakPositionSelection"][k][
                                2
                            ] = (
                                settings_for_fit.fit_orders[i]["PeakPositionSelection"][
                                    k
                                ][2]
                                + move_by
                                # - ((tth_range[1] + tth_range[0]) / 2)
                                # + cent
                            )

                    # re-get settings for current subpattern
                    settings_for_fit.set_subpattern(j, i)

            sub_data = new_data.duplicate()
            sub_data.set_limits(range_bounds=tth_range)

            # Mask the subpattern by intensity if called for
            if (
                "imax" in settings_for_fit.subfit_orders
                or "imin" in settings_for_fit.subfit_orders
            ):
                sub_data = SpotProcess(sub_data, settings_for_fit)

            if mode == "set-range":

                fig_1 = plt.figure()
                sub_data.plot_masked(fig_plot=fig_1)
                plt.suptitle(peak_string(settings_for_fit.subfit_orders) + "; masking")

                filename = make_outfile_name(
                    settings_for_fit.datafile_list[j],
                    directory=settings_for_fit.output_directory,
                    additional_text="mask",
                    orders=settings_for_fit.subfit_orders,
                    extension=".png",
                    overwrite=True,
                )

                fig_1.savefig(filename)

                if debug:
                    plt.show()
                plt.close()

            else:
                tmp = fit_sub_pattern(
                    sub_data,
                    settings_for_fit,  # added
                    None,  # do not pass params they are not needed.
                    save_fit=save_figs,
                    debug=debug,
                    mode=mode,
                    # cascade=True,
                    # cascade_heights=cascade_heights,
                )
                all_fitted_chunks.append(tmp[0])
                all_chunk_positions.append(tmp[1])

        # paste the fits to an output file.
        # store the chunk fits' information as a JSON file.
        filename = make_outfile_name(
            settings_for_fit.subfit_filename,
            directory=settings_for_fit.output_directory,
            additional_text="chunks",
            # orders=settings_for_fit.subfit_orders,
            extension=".json",
            overwrite=True,
        )
        with open(filename, "w") as TempFile:
            # Write a JSON string into the file.
            json.dump(
                (all_fitted_chunks, all_chunk_positions),
                TempFile,
                sort_keys=True,
                indent=2,
                default=json_numpy_serializer,
            )

        # if propagating the fits write them to a temporary file
        if settings_for_fit.fit_propagate:
            # print json_string
            with open(temporary_data_file, "w") as TempFile:
                # Write a JSON string into the file.
                json.dump(
                    (all_fitted_chunks, all_chunk_positions),
                    TempFile,
                    sort_keys=True,
                    indent=2,
                    default=json_numpy_serializer,
                )

    # if parallel is True:
    #     p.close()

    # plot the fits
    plot_cascade_chunks(
        inputs=settings_for_fit, debug=debug, report=report, subpattern="all", **kwargs
    )

    plot_peak_count(
        inputs=settings_for_fit, debug=debug, report=report, subpattern="all", **kwargs
    )



def parallel_processing(p):
    a, kw = p
    return fit_sub_pattern(*a, **kw)


def read_saved_chunks(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    report=False,
    **kwargs
):
    """


    Parameters
    ----------
    setting_file : TYPE, optional
        DESCRIPTION. The default is None.
    setting_class : TYPE, optional
        DESCRIPTION. The default is None.
    inputs : TYPE, optional
        DESCRIPTION. The default is None.
    debug : TYPE, optional
        DESCRIPTION. The default is False.
    report : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    all_azis : TYPE
        DESCRIPTION.
    all_heights : TYPE
        DESCRIPTION.

    """

    if inputs:
        setting_class = inputs
    elif setting_class is None:
        setting_class = initiate(setting_file, report=True)
    else:
        setting_class = setting_class

    all_azis = []
    all_fits = []
    
    print("Reading chunk files")
    logger = proglog.default_bar_logger('bar')  # shorthand to generate a bar logger

    #for f in range(setting_class.image_number):
    for f in logger.iter_bar(iteration=range(setting_class.image_number)):
        setting_class.set_subpattern(f, 0)
        filename = IO.make_outfile_name(
            setting_class.subfit_filename,
            directory=setting_class.output_directory,
            additional_text="chunks",
            extension=".json",
            overwrite=True,
        )
        if os.path.isfile(filename):
            # Read JSON data from file
            with open(filename) as json_data:
                fit = json.load(json_data)
                all_azis.append(fit[1])
                all_fits.append(fit[0])

    print("Finished reading chunk files")
    return all_azis, all_fits

def get_chunks_range(chunks, series="h"):
    
    
    #make length of g
    min_all = np.inf*np.ones(len(chunks[0]))
    max_all = -np.inf*np.ones(len(chunks[0]))
    
    for e in range(len(chunks)):
        for g in range(len(chunks[e])):
            for i in range(len(chunks[e][g][series])):
                v_min = np.min(chunks[e][g][series][i])
                v_max = np.max(chunks[e][g][series][i])
                min_all[g] = np.min([min_all[g], v_min])
                max_all[g] = np.max([max_all[g], v_max])
        
    return min_all, max_all
    

def plot_cascade_chunks(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    report=False,
    plot_type = "timeseries",
    subpattern="all",
    scale="linear",
    azi_range = "all",
    vmax=np.inf,
    vmin=0,
    **kwargs
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param setting_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    if inputs:
        setting_class = inputs
    elif setting_class is None:
        setting_class = initiate(setting_file, inputs=inputs, report=True)
    else:
        setting_class = setting_class

    # get file times
    modified_time_s = []
    for i in range(setting_class.image_number):
        modified_time_s.append(os.path.getmtime(setting_class.image_list[i][0]))
    modified_time_s = np.array(modified_time_s)
    modified_time_s = modified_time_s - modified_time_s[0]
    y_label_str = r"Time (s)"
    # use file numbers if all times are the same
    if len(np.unique(modified_time_s)) == 1:
        modified_time_s = list(range(setting_class.image_number))
        y_label_str = r"Image in sequence"

    # restrict to sub-patterns listed
    setting_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_plots = len(setting_class.fit_orders)
    else:
        num_plots = len(subpattern)

    all_azis, all_data = read_saved_chunks(
        inputs=setting_class, debug=debug, report=report, subpattern=subpattern
    )
    min_all, max_all = get_chunks_range(all_data, series="h")
    print(min_all, max_all)
    
        
    for j in range(num_plots):
        # loop over the number of sets of peaks fit for (i.e. len(setting_class.fit_orders))

        setting_class.set_subpattern(0, j)

        # set plot limits
        if vmax==np.inf:
            vmax=max_all[j]
        #set colour bar ends
        if vmax < max_all[j] and vmin > min_all[j]:
            cb_extend = "both"
        elif vmax < max_all[j]:
            cb_extend = "max"
        elif vmin > min_all[j]:
            cb_extend = "min"
        else:
            cb_extend = "neither"
        #set colour scale normalisation
        norm=None
        if scale == "sqrt":
            norm = colours.PowerNorm(gamma=0.5)
        elif scale == "log":
            norm = colours.LogNorm(vmin=vmin)
        elif scale == "linear":
            norm = None
        for k in range(len(setting_class.subfit_orders["peak"])):
            
            
            if plot_type == "timeseries":
                # loop over the number of peaks in each fit_orders
                print("Making cascade plot for "+ IO.peak_string(setting_class.fit_orders[j], peak=k))
                if all_data[0][j]["h"][k]:
                    # if there is some data in the array plot it.
                    fig, ax = plt.subplots()
                    logger = proglog.default_bar_logger('bar')  # shorthand to generate a bar logger
                    for i in logger.iter_bar(iteration=range(setting_class.image_number)):
                        plt.scatter(
                            all_data[i][j]["chunks"],
                            modified_time_s[i]
                            * np.ones(np.shape(all_data[i][j]["chunks"])),
                            s=.05,
                            c=(all_data[i][j]["h"][k]),
                            #vmax= vmax,
                            #vmin= vmin,
                            cmap="YlOrBr",
                            norm=norm
                        )
                    
                    if azi_range != "all":
                        ax.set_xlim([azi_range[0],azi_range[1]])
                        
                        
                    # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                    pk = k
                    for l in range(len(setting_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = IO.peak_string(setting_class.fit_orders[j], peak=pk)
                    plt.title(ttlstr)
                    plt.xlabel(r"Azimuth (deg)")
                    plt.ylabel(y_label_str)
                    cb = plt.colorbar(extend=cb_extend)
                    cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")
                    plt.show()
                    # save the figure
                    filename = IO.make_outfile_name(
                        setting_class.datafile_basename,
                        directory=setting_class.output_directory,
                        additional_text="CascadePlot",
                        orders=setting_class.subfit_orders,
                        peak=pk,
                        extension=".png",
                        overwrite=True,
                    )
                    fig.savefig(filename, transparent=True, bbox_inches="tight")
                print("\n")

            elif plot_type == "map":
                
                leng = 500
                
                pass
            elif plot_type == "map_video":
                pass
            elif plot_type == "data_cube":
                
                #plot data cube. 
                # interactive 3D graph. (i.e. not spyder inline)
                # x, z are position in space
                # y is azimuth
                # plot points for positions with peaks greater than vmin
                # coloured by azimuth
                # size is size of symbol (and limited by vmin, vmax and normalisation)
                """
                print("Making data cube plot for "+ IO.peak_string(setting_class.fit_orders[j], peak=k))
                if all_data[0][j]["h"][k]:
                    # if there is some data in the array plot it.
                    fig, ax = plt.subplots()
                    logger = proglog.default_bar_logger('bar')  # shorthand to generate a bar logger
                    for i in logger.iter_bar(iteration=range(setting_class.image_number)):
                        points_plot = all_data[i][j]["h"][k] >= vmin
                        x = all_data[i][j]["chunks"][points_plot]
                        x_label = r"Azimuth (deg)"
                        z = modified_time_s[i] * np.ones(np.shape(all_data[i][j]["chunks"]))[points_plot]
                        z_label = y_label_str
                        y = all_data[i][j]["h"][k][points_plot]
                        y_label = 
                        plt.scatter3(
                            x, y, z,
                            s=1,
                            c=(all_data[i][j]["h"][k]),
                            vmax= vmax,
                            vmin= vmin,
                            cmap="YlOrBr",
                            norm=norm
                        )
                    # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                    pk = k
                    for l in range(len(setting_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = IO.peak_string(setting_class.fit_orders[j], peak=pk)
                    plt.title(ttlstr)
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    cb = plt.colorbar(extend=cb_extend)
                    cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")
                    plt.show()
                    # save the figure
                    filename = IO.make_outfile_name(
                        setting_class.datafile_basename,
                        directory=setting_class.output_directory,
                        additional_text="CascadePlot",
                        orders=setting_class.subfit_orders,
                        peak=pk,
                        extension=".png",
                        overwrite=True,
                    )
                    fig.savefig(filename, transparent=True, bbox_inches="tight")
                print("\n")
                """
                pass
            
            elif plot_type == "data_cube_video":
                pass
            else:
                raise ValueError(
                    "Plot type is not recognised."
                )

def peak_count(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    report=False,
    prominence=15,
    subpattern="all",
    **kwargs
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param setting_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    if inputs:
        setting_class = inputs
    elif setting_class is None:
        setting_class = initiate(setting_file, inputs=inputs, report=True)
    else:
        setting_class = setting_class

    # restrict to sub-patterns listed
    setting_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_orders = len(setting_class.fit_orders)
    else:
        num_orders = len(subpattern)

    all_azis, all_data = read_saved_chunks(
        inputs=setting_class, debug=debug, report=report, subpattern=subpattern
    )

    # get the number of peaks
    # all_peaks, all_properties, count = peak_count(setting_class=setting_class, prominence=prominence)

    all_peaks = []
    all_properties = []
    count = []
    peak_labels = []
    for j in range(num_orders):
        # loop over the number of sets of peaks fit for (i.e. len(setting_class.fit_orders))

        setting_class.set_subpattern(0, j)

        for k in range(len(setting_class.subfit_orders["peak"])):
            # loop over the number of peaks in each fit_orders

            all_peaks_tmp = []
            all_properties_tmp = []
            count_tmp = []
            if all_data[0][j]["h"][k]:
                # if there is some data in the array plot it.
                for i in range(setting_class.image_number):
                    peaks, properties = find_peaks(
                        all_data[i][j]["h"][k], prominence=prominence, width=0
                    )

                    all_peaks_tmp.append(peaks)
                    all_properties_tmp.append(properties)
                    count_tmp.append(len(peaks))

                all_peaks.append(all_peaks_tmp)
                all_properties.append(all_properties_tmp)
                count.append(count_tmp)
                # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                pk = k
                for l in range(len(setting_class.subfit_orders["peak"])):
                    if not all_data[0][j]["h"][l]:
                        pk = "all"
                # make the figure title
                ttlstr = IO.peak_string(setting_class.subfit_orders, peak=pk)
                peak_labels.append(ttlstr)

    return all_peaks, all_properties, count, peak_labels


def plot_peak_count(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    report=False,
    prominence=1,
    subpattern="all",
    rotate=False,
    **kwargs
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param setting_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    if inputs:
        setting_class = inputs
    elif setting_class is None:
        setting_class = initiate(setting_file, inputs=inputs, report=True)
    else:
        setting_class = setting_class

    # get file times
    modified_time_s = []
    for i in range(setting_class.image_number):
        modified_time_s.append(os.path.getmtime(setting_class.image_list[i][0]))
    modified_time_s = np.array(modified_time_s)
    modified_time_s = modified_time_s - modified_time_s[0]
    y_label_str = r"Time (s)"
    # use file numbers if all times are the same
    if len(np.unique(modified_time_s)) == 1:
        modified_time_s = list(range(setting_class.image_number))
        y_label_str = r"Image in sequence"

    # restrict to sub-patterns listed
    setting_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_orders = len(setting_class.fit_orders)
    else:
        num_orders = len(subpattern)

    #all_azis, all_data = read_saved_chunks(
    #    inputs=setting_class, debug=debug, report=report, subpattern=subpattern
    #)

    # get the number of peaks
    all_peaks, all_properties, count, titles = peak_count(
        setting_class=setting_class, prominence=prominence
    )

    fig, ax = plt.subplots()
    if not rotate:
        for j in range(len(all_peaks)):
            plt.plot(modified_time_s, count[j], ".-", label=titles[j])
        plt.xlabel(y_label_str)
        plt.ylabel(r"Number peaks")
    else:
        for j in range(len(all_peaks)):
            plt.plot(count[j], modified_time_s, ".-", label=titles[j])
        plt.xlabel(r"Number peaks")
        plt.ylabel(y_label_str)

    plt.legend()
    plt.show()

    filename = IO.make_outfile_name(
        "PeakCountTime",
        directory=setting_class.output_directory,
        additional_text="prominence"+str(prominence),
        extension=".png",
        overwrite=True,
    )
    fig.savefig(filename, transparent=True)


"""
# ============================
# NOTES:
# These are a set of functions that attempted to find the spots in the raw data and count them. 
# They were either very slow or didn't work.
# They are here to remind me to implement them at some time in the future because they are a 'better' way of making the measurements. 
# ============================            

def PeakCount(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):
    
    
    FitSettings, FitParameters, new_data = XRDFit.initiate(settings_file, inputs=inputs, report=True)
    
    # Get list of diffraction patterns #
    diff_files, n_diff_files = IO.FileList(FitParameters, FitSettings)
    
    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60
    
    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]
        
    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])
           
    FitSettings.fit_orders = orders_tmp
    
    
    FitSettings_files = copy.deepcopy(FitSettings.datafile_Files)
    num_peaks_all = []
    visible = []
        
    for h in range(n_diff_files):
        
        #print(FitParameters)
        #print(diff_files)
        
        #FitSettings_tmp = FitSettings
        FitSettings_tmp = FitSettings
        FitSettings_tmp.datafile_Files = [FitSettings_files[h]]
        #FitSettings_tmp.datafile_Files[0] = diff_files[h]
        
        #print(FitSettings_tmp.datafile_Files)
        #print(FitSettings_tmp)
        #print(FitParameters)
    
        outfname = IO.make_outfile_name(diff_files[h], directory=FitSettings_tmp.Output_directory, extension='csv', additional_text='peaks')
        print(outfname)
    
        if not os.path.exists(outfname):
            #This currently gets all the peaks in the image, not just those in the regios of interest
            intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data, 
                    debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                    parallel=parallel,
                    mode='ReturnImage', report=True)
            
            # initialize with default parameters. The "denoise" parameter can be of use in your case
            # import 2D example dataset
            
            plt.imshow(twotheta)
            plt.colorbar()
            
            plt.imshow(np.array(twotheta))
            plt.colorbar()    
            
            img = ma.array(intens)#, mask=False)
            plt.imshow(img)
            # #stop
        
            img.filled(0)
            img = ma.array(img, mask = None)
            plt.imshow(img)
            
            fp = findpeaks(limit=1, whitelist=['peak'], scale=False, togray=False, denoise=None)
            # #stop
            
            # make the fit
            fp.fit(img)
            # Make plot
            #fp.plot()
            #fp.plot_persistence()
            
            fp.results['persistence']
            print(fp.results['persistence'])
            print(type(fp.results['persistence']))
            
            peaks = fp.results
            #print(type(peaks))
            #print(peaks)
            #print(type(peaks))
            
            fp.results['persistence'].to_csv(outfname)
        else:
            
            peaks=pd.read_csv(outfname)
            #print(fp)
            #peaks = fp.results
            if h==0:
                #This currently gets all the peaks in the image, not just those in the regios of interest
                intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data, 
                    debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                    parallel=parallel,
                    mode='ReturnImage', report=True)
            
        # get peaks in each diffraction peak window

        num_peaks = []
        for i in range(len(FitSettings.fit_orders)):
            
            tthRange = FitSettings.fit_orders[i]['range'][0]
    
            msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])
            
            
            # find and list all peaks in unmasked area.
            #FIXME: how does this interact with the region specific masks?
            pks=[]
            for j in range(len(peaks)):
                x = peaks['x'].loc[j]
                y = peaks['y'].loc[j]
                if twotheta.mask[y][x]==False:
                    if ((twotheta[y][x] > tthRange[0]) & (twotheta[y][x] < tthRange[1]) & (twotheta.mask[y][x] == False)):
                        pks.append(peaks.loc[j])
            num_peaks.append(len(pks))
            
            #approximate fraction of ring that is masked.
            #number of pixels in two theta range
            num_pix = ((np.array(twotheta) > tthRange[0]) & (np.array(twotheta) < tthRange[1])).sum()
    
            #number of unmasked pixles in two theta range
            num_vis_pix = ((tthRange[0] < twotheta) & (twotheta < tthRange[1])).sum()
            
            
            visible.append(num_vis_pix/num_pix)
            
            # # find and list all peaks in unmasked area.
            # #FIXME: how does this interact with the region specific masks?
            # peaks = []
            # for j in range(len(fp.results['persistence'])):
            #     x = fp.results['persistence']['x'].loc[j]
            #     y = fp.results['persistence']['y'].loc[j]
            #     if twotheta.mask[y][x]==False:
            #         if ((twotheta[y][x] > tthRange[0]) & (twotheta[y][x] < tthRange[1]) & (twotheta.mask[y][x] == False)):
            #             peaks.append(fp.results['persistence'].loc[j])
            # num_peaks.append(len(peaks))
            
            # #approximate fraction of ring that is masked.
            # #number of pixels in two theta range
            # num_pix = ((np.array(twotheta) > tthRange[0]) & (np.array(twotheta) < tthRange[1])).sum()
    
            # #number of unmasked pixles in two theta range
            # num_vis_pix = ((tthRange[0] < twotheta) & (twotheta < tthRange[1])).sum()
            
            
            # visible.append(num_vis_pix/num_pix)
        num_peaks_all.append(num_peaks)
        #print(num_peaks)
        #print(visible)
     
    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60
    num_peaks_all=np.array(num_peaks_all)
    print(num_peaks_all)
    print(modified_time)
    if 1:
        fig,ax = plt.subplots()
        for i in range(len(FitSettings.fit_orders)):
            print(num_peaks_all[:,i])
            plt.plot(modified_time, num_peaks_all[:,i], '.-', label=IO.peak_string(FitSettings.fit_orders[i]))
            
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Intensity Maxima Count')
        plt.legend()
        plt.show()
            
        fig.savefig('IntensityMaximaTime2.png')  
        fig.savefig('IntensityMaximaTime2.pdf')      
        fig.savefig('IntensityMaximaTime2.eps')    
        
    #plt.scatter([*range(len(num_peaks_all))], num_peaks_all)    
    #plt.scatter
    stop
        
        # #fraction that is unmasked (as approximation for how much of ring is visible.)
        
        # tthRange = FitSettings.fit_orders[i]['range'][0]

        # msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])

        # img = ma.array(intens, mask=msk.mask)
    
        # plt.imshow(img)
        # #stop
    
        # img.filled(0)
        # img = ma.array(img, mask = None)
        # plt.imshow(img)
        
        # fp = findpeaks(limit=2, whitelist=['peak'])
        # #stop
        
        # # make the fit
        # fp.fit(img)
        # # Make plot
        
        # fp.plot()
        
        # fp.plot_persistence()
        
        # fp.results['persistence']
        
        # print(fp.results)
          
    
    
    
    
    
    
    stop
    
    
    
    for i in range(len(FitSettings.fit_orders)):
        
        
        tthRange = FitSettings.fit_orders[i]['range'][0]

        msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])

        img = ma.array(intens, mask=msk.mask)
    
        plt.imshow(img)
        #stop
    
        img.filled(0)
        img = ma.array(img, mask = None)
        plt.imshow(img)
        
        fp = findpeaks(limit=2, whitelist=['peak'])
        #stop
        
        # make the fit
        fp.fit(img)
        # Make plot
        
        fp.plot()
        
        fp.plot_persistence()
        
        fp.results['persistence']
        
        print(fp.results)
        
    stop
                

   
    
    
        
    if 1:
        fig,ax = plt.subplots()
        
        for i in range(len(FitSettings.fit_orders)):
            plt.plot(modified_time, peak_count[:][i], '.-', label=IO.peak_string(FitSettings.fit_orders[i]))
            
        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Intensity Maxima Count')
        plt.legend()
        plt.show()
            
        fig.savefig('IntensityMaximaTime.png')  
        fig.savefig('IntensityMaximaTime.pdf')      
        fig.savefig('IntensityMaximaTime.eps')    
            
        

def WatershedCount(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):
    
    
    FitSettings, FitParameters, new_data = XRDFit.initiate(settings_file, inputs=inputs, report=True)
    
    # Get list of diffraction patterns #
    diff_files, n_diff_files = IO.FileList(FitParameters, FitSettings)
    
    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60
    
    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]
        
    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])
           
    FitSettings.fit_orders = orders_tmp
    
    
    FitSettings_files = copy.deepcopy(FitSettings.datafile_Files)
    num_peaks_all = []
    visible = []
        
    for h in range(n_diff_files):
        
        #This currently gets all the peaks in the image, not just those in the regios of interest
        intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data, 
                debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                parallel=parallel,
                mode='ReturnImage', report=True)     
        
        im = img_as_float(intens)
        
        # image_max is the dilation of im with a 20*20 structuring element
        # It is used within peak_local_max function
        image_max = ndi.maximum_filter(im, size=20, mode='constant')
        
        # Comparison between image_max and im to find the coordinates of local maxima
        coordinates = peak_local_max(im, min_distance=20)
        
        # display results
        fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharex=True, sharey=True)
        ax = axes.ravel()
        ax[0].imshow(im, cmap=plt.cm.gray)
        ax[0].axis('off')
        ax[0].set_title('Original')
        
        ax[1].imshow(image_max, cmap=plt.cm.gray)
        ax[1].axis('off')
        ax[1].set_title('Maximum filter')
        
        ax[2].imshow(im, cmap=plt.cm.gray)
        ax[2].autoscale(False)
        ax[2].plot(coordinates[:, 1], coordinates[:, 0], 'r.')
        ax[2].axis('off')
        ax[2].set_title('Peak local max')
        
        fig.tight_layout()
        
        plt.show()
        




"""
