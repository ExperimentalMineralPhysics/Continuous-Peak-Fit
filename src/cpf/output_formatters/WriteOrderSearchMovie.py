__all__ = ["Requirements", "WriteOutput"]

import glob
import json
import os
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np

# from moviepy import concatenate
from moviepy import ImageClip, VideoFileClip, concatenate_videoclips

from cpf.BrightSpots import SpotProcess
from cpf.output_formatters.fits_io import ReadFits_to_dataframe
from cpf.settings import get_settings
from cpf.util.io import (
    figure_suptitle_space,
    make_outfile_name,
    peak_string,
    title_file_names,
)
from cpf.util.logging import get_logger
from cpf.util.output_formatters import mplfig_to_npimage
from cpf.XRD_FitSubpattern import plot_FitAndModel

logger = get_logger("cpf.output_formatters.WriteOrderSearchMovie")


def Requirements():
    """
    List required and optional parameters for making outputs using
    OrderSearchMovie.

    Returns
    -------
    RequiredParams : list
        Parameters that the Output type cannot be run without, excluding the
        settings class or parameters therein.
    OptionalParams : list
        Parameters that change the output but have default settings.
    """

    # OptionalParams is for parameters that affect the contents or types of the written files.
    # Does not contain settings related to the processing.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        "fps": 10,  # frames per second
        "file_types": ["mp4"],  # movie file type
        "Irange": ["pt1percentile", "99pt9percentile"],  # range of colour scale
    }

    # OptionalParams = [ # list of parameters parsed as kwargs.
    #     "file_types"  # -- pick which type of movis to write
    #     "fps"         # -- frames per second
    # ]
    return RequiredParams, OptionalParams


def WriteOutput(
    settings,
    file_label=None,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    *args,
    **kwargs,
):
    """
    Writes a movie file of the outputs from cpf.XRD_FitPattern.order_search.

    N.B. this output requires the data files to be present.

    Parameters
    ----------
    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or
        dictionary of all the settings or
        string or path to a file with the settings in.
    file_label : str, optional
        String added to the names of the files made. The default is None.
    report : Literal[        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"    ], optional
        Sets logger level and amount of information writen to log file. The default is "INFO".
    *args : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        Raised when neither settings_class and settings_file are presuent

    Returns
    -------
    None.

    """

    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters
    fps = settings_class.output_settings.get("fps", Requirements()[1]["fps"])
    file_types = settings_class.output_settings.get(
        "file_types", Requirements()[1]["file_types"]
    )
    Irange = settings_class.output_settings.get("Irange", Requirements()[1]["Irange"])
    # override with kwargs
    fps = kwargs.get("fps", fps)
    file_types = kwargs.get("file_types", file_types)
    Irange = kwargs.get("Irange", Irange)

    if not isinstance(file_types, list):
        file_types = [file_types]
    if not (isinstance(fps, float) or isinstance(fps, int)):
        raise ValueError("The frames per second needs to be a number.")

    if file_label is not None:
        settings_class.file_label = file_label

    # make the data class.
    data_to_fill = settings_class.image_list[0]
    data_class = settings_class.data_class
    data_class.fill_data(
        data_to_fill,
        settings=settings_class,
    )

    # restrict to just the first file (as per order search)
    settings_class.set_data_files(start=0, end=1)

    # this is search data so there is a postscript in the json file label.
    # determine the label
    if "file_label" not in dir(settings_class) or settings_class.file_label is None:
        fls = glob.glob("./results/*search*.json")
        if len(fls) == 0:
            raise ValueError("There is no identified search file to plot.")
        else:
            tm = []
            for i in range(len(fls)):
                tm.append(os.path.getmtime(fls[i]))
            latest = np.argsort(tm)[-1]
            settings_class.file_label = os.path.splitext(os.path.basename(fls[latest]))[
                0
            ].split("__")[1]

    # read the data.
    df = ReadFits_to_dataframe(
        settings=settings_class,
        includeStats=True,
        includeSeriesValues=True,
        includeIntensityRanges=True,
        includePosition=True,
    )
    headers = list(df.columns.values)
    # split the notes column into columns and calculate some new values
    f = lambda x: x.split("|")[0].split("=")[1]
    df["search_peak"] = df["note"].apply(f).astype(int)
    f = lambda x: x.split("|")[1].split("=")[0]
    df["search_over"] = df["note"].apply(f)
    f = lambda x: float(x.split("|")[1].split("=")[1])
    df["search_value"] = df["note"].apply(f).astype(int)
    f = lambda x: x.split("|")[2].split("=")[1]
    df["series_type"] = df["note"].apply(f)

    peaks = df["peak"].unique()
    searches = df["series_type"].unique()
    search_value = df["search_value"].unique()

    # open data file for plotting information
    # read fit file
    json_file = make_outfile_name(
        settings_class.subfit_filename,
        directory=settings_class.output_directory,
        additional_text=settings_class.file_label,
        extension=".json",
        overwrite=True,
    )
    with open(json_file) as json_data:
        data_fit = json.load(json_data)

    # make movies
    for i in range(len(peaks)):
        # loop over the number of unique peaks

        # FIXME: here we are looping over every occurence but if there are multiple peaks in a
        # range we could be making the videos more than once. So need to check and make
        # video for only the peak that is being searched over
        # this requires filtering for a change in search value.

        # use pos_in_range and _search peak to restrict values to the peak that
        # is being searched over -- for multiple peaks in range
        df_peak = df[
            (df["peak"] == peaks[i]) & (df["pos_in_range"] == df["search_peak"])
        ]

        # make sure that we have only the peak that we are looping over.
        if len(df_peak["search_peak"].unique()) != 1:
            raise ValueError("There should be more than 1 peak in here somewhere")

        # sort the data into a sensible order
        df_peak = df_peak.sort_values(["series_type", "search_value"])

        # get the intensity ranges.
        Intensity_range = [
            np.nanmin([df_peak["data_min"], df_peak["model_min"]]),
            np.nanmax([df_peak["data_max"], df_peak["model_max"]]),
        ]
        Resid_range = [
            np.nanmin([df_peak["residual_min"]]),
            np.nanmax([df_peak["residual_max"]]),
        ]

        # setup figure
        fig = plt.figure()

        # interate over the length of the selection of peaks in df_peak
        frames = []
        for j in range(len(df_peak)):
            # get position in data frame
            position = df_peak.iloc[j]["Position in json"]
            # print(j, position)
            # position = df_peak["Position in json"].iloc[j]

            settings_class.set_subpattern(0, position)

            # print(j, position)

            sub_data = data_class.duplicate_without_detector(
                range_bounds=[
                    df_peak["range_start"].iloc[j],
                    df_peak["range_end"].iloc[j],
                ]
            )
            # sub_data = data_class.duplicate()
            # sub_data.set_limits(range_bounds=[df_peak["Range_start"].iloc[j], df_peak["Range_end"].iloc[j]])

            # Mask the subpattern by intensity if called for
            if (
                "imax" in settings_class.subfit_orders
                or "imin" in settings_class.subfit_orders
            ):
                sub_data = SpotProcess(sub_data, settings_class)

            # get position in original input file. (before order search)
            # basically set the correct subpattern
            if len(settings_class.fit_orders) == len(df):
                # the rows correspond
                settings_class.set_subpattern(0, j)
            else:
                # find the row in which the phases corrospond
                for m in range(len(settings_class.fit_orders)):
                    if peaks[i] in peak_string(settings_class.fit_orders[m]):
                        settings_class.set_subpattern(0, m)
                        # print("subpattern", m)

            # fig = plt.figure()
            # make the plot of the fits.
            fig = plot_FitAndModel(
                settings_class,
                sub_data,
                # param_lmfit=None,
                params_dict=data_fit["fits"][position],
                figure=fig,
                plot_ColourRange={
                    "max": Intensity_range[0],
                    "min": Intensity_range[-1],
                    "rmin": Resid_range[0],
                    "rmax": Resid_range[-1],
                },
            )
            title_str = (
                peak_string(settings_class.subfit_orders)
                + "; "
                + df_peak["note"].iloc[j]
            )
            # if "note" in settings_class.subfit_orders:
            #     title_str += " " + settings_class.subfit_orders["note"]
            plt.suptitle(title_str)
            figure_suptitle_space(fig, topmargin=0.4)
            # print(type(fig))
            # try:
            #     fig.show()
            # except:
            #     # print("who knows")
            #     pass
            # make the video clip
            # just addes the figure as a frame to the proto-video.
            frames.append(ImageClip(mplfig_to_npimage(fig)).with_duration(1))

        # convert to video and write
        video = concatenate_videoclips(frames, method="compose")
        for f in range(len(file_types)):
            settings_class.file_label = "search=" + df["search_over"][position]
            # make videofile name
            data_fit_tmp = data_fit["fits"][position]
            if "note" in data_fit_tmp:
                data_fit_tmp.pop("note")
            out_file = make_outfile_name(
                settings_class.subfit_filename,
                directory=settings_class.output_directory,
                orders=data_fit_tmp,
                additional_text=settings_class.file_label,
                extension=file_types[f],
                peak=df_peak["search_peak"].iloc[0],
                overwrite=True,
            )
            video.write_videofile(out_file, fps=fps)
