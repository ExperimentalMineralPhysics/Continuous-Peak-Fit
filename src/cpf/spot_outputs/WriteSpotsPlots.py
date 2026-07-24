__all__ = ["Requirements", "WriteOutput"]


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

logger = get_logger("cpf.output_formatters.WriteSpotsPlots")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        "parameters_plot": ["d_mean", "differential", "d_max", "height mean"],
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
    fitStats=True,
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

    # kwargs.update(set_params)
    
    # read the data.
    spot_df = ReadSpots_to_dataframe(
        settings=settings_class,
    )
    headers = list(spot_df.columns.values)


    plot_spot_count(settings_class, spot_df, plot_against="something")



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
        
def spot_count(
    settings_class,
    spots_data,
    # report: Literal[
    #     "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    # ] = "INFO",
    **kwargs,
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """


    # counts = pd.DataFrame()
    
    unique_peaks = []


    combined = []
    combined_headers = []
    for fls in range(settings_class.image_number):
        
        combined_headers = ["DataFile"]
        
        combined_tmp = []
        combined_tmp.append(os.path.split(settings_class.image_list[fls])[1])
        
        if 'image_position' in spots_data:
            combined_tmp.append(spots_data["image_position"].loc[spots_data["DataFile"] == os.path.split(settings_class.image_list[fls])[1]].iloc[0])
            combined_headers += ["image_position"]
            
        #add metadata
        for md in settings_class.metadata:
            combined_headers += [md]
            combined_tmp.append(spots_data[md].loc[spots_data["DataFile"] == os.path.split(settings_class.image_list[fls])[1][0]])
        
        for i in range(len(settings_class.fit_orders)):
            for j in range(len(settings_class.fit_orders[i]["peak"])):
                phss = peak_phase(settings_class.fit_orders[i], peak = j)[0]
                pks = peak_hkl(settings_class.fit_orders[i], peak=j, as_string=True)[0]
                peak_str = peak_string(settings_class.fit_orders[i], peak=j)
                
                unique_peaks.append(peak_str)
                # f"{phss} ({pks})")
                num_spots = len(spots_data.loc[
                        ((spots_data["DataFile"] == os.path.split(settings_class.image_list[fls])[1]) & 
                         (spots_data['phases']==phss) &
                         (spots_data['peak'] == pks)
                         )]
                        )
                
                combined_headers += [peak_str]
                combined_tmp.append(num_spots)
        combined.append(combined_tmp)
        
    combined = pd.DataFrame(combined, columns=combined_headers)
    return combined, np.unique(unique_peaks)
    # stop
                
                
    # pk_count = {}
    # peaks_count = []
    # for phss in unique_phases:
    #     pk_count[phss] = pk_count.get(phss, {})
    #     for pks in unique_peaks:
    #         num_spots = []
    #         for fls in range(settings_class.image_number):
    #             print(fls, os.path.split(settings_class.image_list[fls])[1])
    #             peaks_count.append(f"{phss} {pks}")
    #             num_spots.append(
    #                 len(spots_data.loc[
    #                     ((spots_data["DataFile"] == os.path.split(settings_class.image_list[fls])[1]) & 
    #                      (spots_data['phases']==phss) &
    #                      (spots_data['peak'] == pks)
    #                      )]
    #                     )
    #                 )
            
                
    #             # print(np.array(spots_data.loc[
    #             #     ((spots_data["DataFile"] == os.path.split(settings_class.image_list[fls])[1]) & 
    #             #      (spots_data['phases']==phss) &
    #             #      (spots_data['peak'] == pks)
    #             #      )]["2theta"]))
    #         pk_count[phss][pks] = num_spots
    # stop
    # print(pk_count)
    
    # return pk_count
    # stop
    # # get the number of peaks
    # # all_peaks, all_properties, count = peak_count(settings_class=settings_class, prominence=prominence)

    # all_peaks: list[list] = []
    # all_peakAzis: list = []
    # all_properties: list[list] = []
    # count = []
    # peak_labels = []
    # for j in range(num_orders):
    #     # loop over the number of sets of peaks fit for (i.e. len(settings_class.fit_orders))

    #     settings_class.set_subpattern(0, j)

    #     for k in range(len(settings_class.subfit_orders["peak"])):
    #         # loop over the number of peaks in each fit_orders

    #         all_peaks_tmp = []
    #         all_properties_tmp = []
    #         count_tmp = []
    #         if all_data[0][j]["h"][k]:
    #             # if there is some data in the array plot it.
    #             for i in range(settings_class.image_number):
    #                 peaks, properties = find_peaks(
    #                     all_data[i][j]["h"][k], prominence=prominence, width=0
    #                 )
    #                 properties["PeakAzis"] = np.array(all_azis[i][j])[peaks]
    #                 all_peaks_tmp.append(peaks)
    #                 # all_peakAzis.append(np.array(all_azis[i][j])[peaks])
    #                 all_properties_tmp.append(properties)
    #                 count_tmp.append(len(peaks))

    #                 if 1 and i == 0:
    #                     fig = plt.figure()
    #                     plt.plot(all_data[i][j]["h"][k])
    #                     plt.plot(peaks, np.array(all_data[i][j]["h"][k])[peaks], "x")

    #                     # plt.plot(np.zeros_like(x), "--", color="gray")

    #                     if show_plots is True:
    #                         plt.show()
    #                     else:
    #                         plt.close()

    #             all_peaks.append(all_peaks_tmp)
    #             all_properties.append(all_properties_tmp)
    #             count.append(count_tmp)
    #             # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
    #             pk: int | Literal["all"] = k
    #             for l in range(len(settings_class.subfit_orders["peak"])):
    #                 if not all_data[0][j]["h"][l]:
    #                     pk = "all"
    #             # make the figure title
    #             ttlstr = peak_string(settings_class.subfit_orders, peak=pk)
    #             peak_labels.append(ttlstr)

    # return all_peaks, all_properties, count, peak_labels


def plot_spot_count(
    settings_class,
    spots_data,
    subpattern: str = "all",
    plot_against = "position",
    # report: Literal[
    #     "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    # ] = "INFO",
    **kwargs,
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    plot_orientation = settings_class.output_settings.get(
        "plot_orientation", Requirements()[1]["plot_orientation"]
    )
    plot_orientation = kwargs.get("plot_orientation", plot_orientation)

    # # get file times
    # modified_time_s = np.array(
    #     [Path(file).stat().st_mtime for file in settings_class.image_list]
    # )
    # modified_time_s -= float(modified_time_s[0])

    # y_label_str = r"Time (s)"
    # # use file numbers if all times are the same
    # if len(np.unique(modified_time_s)) == 1:
    #     modified_time_s = list(range(settings_class.image_number))
    #     y_label_str = r"Image in sequence"

    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_orders = len(settings_class.fit_orders)
    else:
        num_orders = len(subpattern)

    # get the number of peaks
    counted_peaks, peak_names = spot_count(
        settings_class, spots_data
    )

    # get axis to plot against 
    if plot_against in list(counted_peaks):
        independent_variable = counted_peaks[plot_against]
        independent_label_str = plot_against
    else:
        independent_variable = counted_peaks["image_position"]+1
        independent_label_str = "image position in series"

    fig, ax = plt.subplots()
    if plot_orientation == "vertical":
        for kys in counted_peaks.keys():
            for pks in counted_peaks[kys]:
                plt.scatter(counted_peaks[kys][pks], independent_variable, ".-", label=f"{kys} ({pks})")
        plt.xlabel(independent_label_str)
        plt.ylabel(r"Number peaks")
    else:
        for kys in peak_names:
            plt.plot(independent_variable, counted_peaks[kys], ".-", label=kys)
        plt.ylabel(r"Number spots")
        plt.xlabel(independent_label_str)
    plt.legend()

    # Save the plot
    filename = make_outfile_name(
        "PeakCountTime",
        directory=settings_class.output_directory,
        extension=".png",
        overwrite=True,
    )
    fig.savefig(filename, transparent=True)

    plt.show()
    # if show_plots is True:
    #     plt.show()
    # else:
    #     plt.close()
        
        


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

