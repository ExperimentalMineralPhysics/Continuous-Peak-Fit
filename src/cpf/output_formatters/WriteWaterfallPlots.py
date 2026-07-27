__all__ = ["Requirements", "WriteOutput"]


# TODO : move to spot finding out put folder and import properly there. 


# import json
import os
# import re
# from itertools import product

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

# import cpf.peak_functions as pf
from cpf.output_formatters.fits_io import ReadSpots_to_dataframe
# from cpf.output_formatters.output_csv import make_header, write_csv
from cpf.util.io import make_outfile_name, peak_phase, peak_hkl, peak_string
from cpf.settings import get_settings
from cpf.util.logging import get_logger

import matplotlib.colors as colours
import proglog 
logger = get_logger("cpf.output_formatters.WriteSpotsPlots")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        "plot_type": "timeseries",
        "scale": "sqrt",
        "plot_orientation": "horizontal",
        "plot_against": "image position"
        
    }
    # OptionalParams = [
    #     "SampleGeometry"  # changes the strain tensor calucaltion from 2d to 3d. This determines how the cetroid and differnetial strain of the dpsaice are extracted from the fourier series.
    #     "SampleDeformation"  # changes calculation between 'compression' and 'extension'.
    #     "coefs_vals_write"  # -- pick which set of coefficients to write
    #     "ordering_of_output": None # Just leave as read -- otherwise list of dataframe headers to order by
    # ]

    return RequiredParams, OptionalParams


def WriteOutput(
    settings,
    *args,
    **kwargs,
):
    """
    Write coefficents from fits to table/csv file.


    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or
        dictionary of all the settings or
        string or path to a file with the settings in.
    :param fitStats: DESCRIPTION, defaults to True
    :type fitStats: TYPE, optional
    :param *args: DESCRIPTION
    :type *args: TYPE
    :param **kwargs: DESCRIPTION
    :type **kwargs: TYPE
    :raises ValueError: DESCRIPTION
    :return: DESCRIPTION
    :rtype: TYPE

    """

    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters

    # override with kwargs
    subpattern = kwargs.get("subpattern", "all")
    plot_type = kwargs.get("plot_type", Requirements()[1]["plot_type"])
    scale = kwargs.get("scale", Requirements()[1]["scale"])
    plot_against = kwargs.get("plot_against", Requirements()[1]["plot_against"])


    # kwargs.update(set_params)
    
    # read the data.
    spot_df = ReadSpots_to_dataframe(
        settings=settings_class,
    )
    headers = list(spot_df.columns.values)


    # plot_spot_count(settings_class, spot_df, plot_against="something")

    import cpf.spot_methods.spots_scipy_peakfind as spots_scipy_peakfind
    
    all_azis, all_data = spots_scipy_peakfind.read_intermediate_files(settings_class)

    min_all, max_all = spots_scipy_peakfind.get_chunks_range(all_data, series="h")
    
    
    if subpattern == "all":
        num_plots = len(settings_class.fit_orders)
    else:
        num_plots = len(subpattern)
    print(num_plots)
    for j in range(num_plots):
        # loop over the number of sets of peaks fit for (i.e. len(settings_class.fit_orders))
    
        settings_class.set_subpattern(0, j)
        print(j)
        # set plot limits
        if "min_all" in locals():
            vmin = min_all[j]
        if "max_all" in locals():
            vmax = max_all[j]
        # set colour bar ends
        if vmax < max_all[j] and vmin > min_all[j]:
            cb_extend = "both"
        elif vmax < max_all[j]:
            cb_extend = "max"
        elif vmin > min_all[j]:
            cb_extend = "min"
        else:
            cb_extend = "neither"
        # set colour scale normalisation
        norm = None
        if scale == "sqrt":
            norm = colours.PowerNorm(gamma=0.5,vmin=vmin, vmax=vmax)
        elif scale == "log":
            norm = colours.LogNorm(vmin=vmin, vmax=vmax)
        elif scale == "linear":
            norm = None
        for k in range(len(settings_class.subfit_orders["peak"])):
            if plot_type == "timeseries":
                # loop over the number of peaks in each fit_orders
                print(
                    "Making cascade plot for "
                    + peak_string(settings_class.fit_orders[j], peak=k)
                )
                if all_data[0][j]["h"][k]:
                    # if there is some data in the array plot it.
                    fig, ax = plt.subplots()
                    logger = proglog.default_bar_logger(
                        "bar"
                    )  # shorthand to generate a bar logger
                    for i in logger.iter_bar(
                        iteration=range(settings_class.image_number)
                    ):
                        if 0:  # len(all_data[i][j]["h"]) == 1:
                            tms = modified_time_s[j]
                            independent_label_str = "time"
                        elif 0:
                            tms = modified_time_s[i] * np.ones(
                                np.shape(all_data[i][j]["chunks"])
                            )
                            independent_label_str = "time"
                        else:
                            tms = i * np.ones(
                                np.shape(all_data[i][j]["chunks"])
                            )
                            independent_label_str = "image position in series"
    
                        plt.scatter(
                            all_data[i][j]["chunks"],
                            tms, #modified_time_s[i] * np.ones(np.shape(all_data[i][j]["chunks"])),
                            s=10,  # s=.05,
                            c=(all_data[i][j]["h"][k]),
                            # vmax=vmax,
                            # vmin=vmin,
                            # cmap="YlOrBr",
                            norm=norm,
                        )
    
                    # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                    pk: int | Literal["all"] = k
                    for l in range(len(settings_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = peak_string(settings_class.fit_orders[j], peak=pk)
                    plt.title(ttlstr)
                    plt.xlabel(r"Azimuth (deg)")
                    plt.ylabel(independent_label_str)
    
                    if 1:#azi_range == "all":
                        azi_range = [
                            np.min(all_data[i][j]["chunks"]),
                            np.max(all_data[i][j]["chunks"]),
                        ]
                    x_ticks = settings_class.data_class.dispersion_ticks(
                        disp_lims=azi_range
                    )
                    # x.set_xticks(x_ticks)
                    cb = plt.colorbar(extend=cb_extend)
                    # cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")
    
                    # Save the figure
                    filename = make_outfile_name(
                        settings_class.datafile_basename,
                        directory=settings_class.output_directory,
                        additional_text="CascadePlot",
                        orders=settings_class.subfit_orders,
                        peak=pk,
                        extension=".png",
                        overwrite=True,
                    )
                    fig.savefig(filename, transparent=True, bbox_inches="tight")
    
                    if 1:#show_plots is True:
                        plt.show()
                    else:
                        plt.close()
    
            elif plot_type == "map":
                leng = 500
    
                pass
            elif plot_type == "map_video":
                pass
            elif plot_type == "data_cube":
                # plot data cube.
                # interactive 3D graph. (i.e. not spyder inline)
                # x, z are position in space
                # y is azimuth
                # plot points for positions with peaks greater than vmin
                # coloured by azimuth
                # size is size of symbol (and limited by vmin, vmax and normalisation)
                """
                print("Making data cube plot for "+ IO.peak_string(settings_class.fit_orders[j], peak=k))
                if all_data[0][j]["h"][k]:
                    # if there is some data in the array plot it.
                    fig, ax = plt.subplots()
                    logger = proglog.default_bar_logger('bar')  # shorthand to generate a bar logger
                    for i in logger.iter_bar(iteration=range(settings_class.image_number)):
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
                    for l in range(len(settings_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = IO.peak_string(settings_class.fit_orders[j], peak=pk)
                    plt.title(ttlstr)
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    cb = plt.colorbar(extend=cb_extend)
                    cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")
                    plt.show()
                    # save the figure
                    filename = IO.make_outfile_name(
                        settings_class.datafile_basename,
                        directory=settings_class.output_directory,
                        additional_text="CascadePlot",
                        orders=settings_class.subfit_orders,
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
                raise ValueError("Plot type is not recognised.")
    
    














    # order the rows    
    
    # for prm in ParametersPlot:
        
    #     fig, ax = plt.subplots(nrows=1, ncols=1)
    
    #     plt_x = "image_position" # num when written to file
    #     plt_y = prm
    #     peaks = np.unique(df['peak'])
    #     for i in peaks:
    #         lbl = f"({i})"
    #         ax.plot(df[plt_x].loc[df['peak'] == i], df[plt_y].loc[df['peak'] == i],'.', label=lbl)
    #     ax.legend()
    #     ax.set_yscale("log", nonpositive='clip')
    #     ax.set_xlabel(plt_x)
    #     ax.set_ylabel(plt_y)
    #     ax.set_title(prm)
        

# %% 
# =============================================================================
# HERE AFTER ARE FUNCIOTNS THAT SHOULD BE MOEVED INTO OUTPUTS
# =============================================================================
        

# def plot_cascade_chunks(
#     settings_file=None,
#     settings_class=None,
#     inputs=None,
#     debug=False,
#     report: Literal[
#         "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
#     ] = "INFO",
#     plot_type="timeseries",
#     subpattern="all",
#     scale="linear",
#     azi_range="all",
#     vmax=np.inf,
#     vmin=0,
#     show_plots: bool = False,
#     **kwargs,
# ):
#     """
#     :param fit_parameters:
#     :param fit_settings:
#     :param settings_file:
#     :param parallel:
#     :param report:
#     :param mode:
#     :param save_all:
#     :param inputs:
#     :param debug:
#     :param subpattern:
#     :return:
#     """

#     if inputs:
#         settings_class = inputs
#     elif settings_class is None:
#         settings_class = initiate(settings_file, inputs=inputs, report=report)
#     else:
#         settings_class = settings_class

#     # get file times
#     modified_time_s = np.array(
#         [Path(file).stat().st_mtime for file in settings_class.image_list]
#     )
#     modified_time_s -= float(modified_time_s[0])

#     y_label_str = r"Time (s)"
#     # use file numbers if all times are the same
#     if len(np.unique(modified_time_s)) == 1:
#         modified_time_s = list(range(settings_class.image_number))
#         y_label_str = r"Image in sequence"

#     # restrict to sub-patterns listed
#     settings_class.set_subpatterns(subpatterns=subpattern)

#     if subpattern == "all":
#         num_plots = len(settings_class.fit_orders)
#     else:
#         num_plots = len(subpattern)

#     all_azis, all_data = read_saved_chunks(
#         inputs=settings_class, debug=debug, report=report, subpattern=subpattern
#     )
#     min_all, max_all = get_chunks_range(all_data, series="h")

#     for j in range(num_plots):
#         # loop over the number of sets of peaks fit for (i.e. len(settings_class.fit_orders))

#         settings_class.set_subpattern(0, j)

#         # set plot limits
#         if "min_all" in locals():
#             vmin = min_all[j]
#         if "max_all" in locals():
#             vmax = max_all[j]
#         # set colour bar ends
#         if vmax < max_all[j] and vmin > min_all[j]:
#             cb_extend = "both"
#         elif vmax < max_all[j]:
#             cb_extend = "max"
#         elif vmin > min_all[j]:
#             cb_extend = "min"
#         else:
#             cb_extend = "neither"
#         # set colour scale normalisation
#         norm = None
#         if scale == "sqrt":
#             norm = colours.PowerNorm(gamma=0.5)
#         elif scale == "log":
#             norm = colours.LogNorm(vmin=vmin)
#         elif scale == "linear":
#             norm = None
#         for k in range(len(settings_class.subfit_orders["peak"])):
#             if plot_type == "timeseries":
#                 # loop over the number of peaks in each fit_orders
#                 print(
#                     "Making cascade plot for "
#                     + peak_string(settings_class.fit_orders[j], peak=k)
#                 )
#                 if all_data[0][j]["h"][k]:
#                     # if there is some data in the array plot it.
#                     fig, ax = plt.subplots()
#                     logger = proglog.default_bar_logger(
#                         "bar"
#                     )  # shorthand to generate a bar logger
#                     for i in logger.iter_bar(
#                         iteration=range(settings_class.image_number)
#                     ):
#                         if 1:  # len(all_data[i][j]["h"]) == 1:
#                             tms = modified_time_s[j]
#                         else:
#                             tms = modified_time_s[i] * np.ones(
#                                 np.shape(all_data[i][j]["chunks"])
#                             )

#                         plt.scatter(
#                             all_data[i][j]["chunks"],
#                             modified_time_s[i]
#                             * np.ones(np.shape(all_data[i][j]["chunks"])),
#                             s=10,  # s=.05,
#                             c=(all_data[i][j]["h"][k]),
#                             vmax=vmax,
#                             vmin=vmin,
#                             # cmap="YlOrBr",
#                             norm=norm,
#                         )

#                     # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
#                     pk: int | Literal["all"] = k
#                     for l in range(len(settings_class.subfit_orders["peak"])):
#                         if not all_data[0][j]["h"][l]:
#                             pk = "all"
#                     # make the figure title
#                     ttlstr = peak_string(settings_class.fit_orders[j], peak=pk)
#                     plt.title(ttlstr)
#                     plt.xlabel(r"Azimuth (deg)")
#                     plt.ylabel(y_label_str)

#                     if azi_range == "all":
#                         azi_range = [
#                             np.min(all_data[i][j]["chunks"]),
#                             np.max(all_data[i][j]["chunks"]),
#                         ]
#                     x_ticks = settings_class.data_class.dispersion_ticks(
#                         disp_lims=azi_range
#                     )
#                     # x.set_xticks(x_ticks)
#                     cb = plt.colorbar(extend=cb_extend)
#                     # cb.set_label(r"Log$_{10}$(Intensity)")
#                     cb.set_label(r"Intensity")

#                     # Save the figure
#                     filename = make_outfile_name(
#                         settings_class.datafile_basename,
#                         directory=settings_class.output_directory,
#                         additional_text="CascadePlot",
#                         orders=settings_class.subfit_orders,
#                         peak=pk,
#                         extension=".png",
#                         overwrite=True,
#                     )
#                     fig.savefig(filename, transparent=True, bbox_inches="tight")

#                     if show_plots is True:
#                         plt.show()
#                     else:
#                         plt.close()

#             elif plot_type == "map":
#                 leng = 500

#                 pass
#             elif plot_type == "map_video":
#                 pass
#             elif plot_type == "data_cube":
#                 # plot data cube.
#                 # interactive 3D graph. (i.e. not spyder inline)
#                 # x, z are position in space
#                 # y is azimuth
#                 # plot points for positions with peaks greater than vmin
#                 # coloured by azimuth
#                 # size is size of symbol (and limited by vmin, vmax and normalisation)
#                 """
#                 print("Making data cube plot for "+ IO.peak_string(settings_class.fit_orders[j], peak=k))
#                 if all_data[0][j]["h"][k]:
#                     # if there is some data in the array plot it.
#                     fig, ax = plt.subplots()
#                     logger = proglog.default_bar_logger('bar')  # shorthand to generate a bar logger
#                     for i in logger.iter_bar(iteration=range(settings_class.image_number)):
#                         points_plot = all_data[i][j]["h"][k] >= vmin
#                         x = all_data[i][j]["chunks"][points_plot]
#                         x_label = r"Azimuth (deg)"
#                         z = modified_time_s[i] * np.ones(np.shape(all_data[i][j]["chunks"]))[points_plot]
#                         z_label = y_label_str
#                         y = all_data[i][j]["h"][k][points_plot]
#                         y_label =
#                         plt.scatter3(
#                             x, y, z,
#                             s=1,
#                             c=(all_data[i][j]["h"][k]),
#                             vmax= vmax,
#                             vmin= vmin,
#                             cmap="YlOrBr",
#                             norm=norm
#                         )
#                     # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
#                     pk = k
#                     for l in range(len(settings_class.subfit_orders["peak"])):
#                         if not all_data[0][j]["h"][l]:
#                             pk = "all"
#                     # make the figure title
#                     ttlstr = IO.peak_string(settings_class.fit_orders[j], peak=pk)
#                     plt.title(ttlstr)
#                     plt.xlabel(x_label)
#                     plt.ylabel(y_label)
#                     cb = plt.colorbar(extend=cb_extend)
#                     cb.set_label(r"Log$_{10}$(Intensity)")
#                     cb.set_label(r"Intensity")
#                     plt.show()
#                     # save the figure
#                     filename = IO.make_outfile_name(
#                         settings_class.datafile_basename,
#                         directory=settings_class.output_directory,
#                         additional_text="CascadePlot",
#                         orders=settings_class.subfit_orders,
#                         peak=pk,
#                         extension=".png",
#                         overwrite=True,
#                     )
#                     fig.savefig(filename, transparent=True, bbox_inches="tight")
#                 print("\n")
#                 """
#                 pass

#             elif plot_type == "data_cube_video":
#                 pass
#             else:
#                 raise ValueError("Plot type is not recognised.")

