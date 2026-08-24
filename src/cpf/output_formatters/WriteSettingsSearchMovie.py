__all__ = ["Requirements", "WriteOutput"]

import glob
import json
import os
import pandas as pd
from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np

from copy import deepcopy

# # from moviepy import concatenate
from moviepy import ImageClip, VideoFileClip, concatenate_videoclips

from cpf.BrightSpots import SpotProcess
from cpf.output_formatters.fits_io import ReadFits_to_dataframe, ReadFits_to_list
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
        "pattern": None, # guess from files in outpuct directory
        "search_parameter": "", # guess from files in outpuct directory
        "search_over": None, # guess from files in outpuct directory
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
    pattern = settings_class.output_settings.get("pattern", Requirements()[1]["pattern"])
    search_parameter = settings_class.output_settings.get("search_parameter",Requirements()[1]["search_parameter"])
    search_over = settings_class.output_settings.get("search_over", Requirements()[1]["search_over"])
    fps = settings_class.output_settings.get("fps", Requirements()[1]["fps"])
    file_types = settings_class.output_settings.get(
        "file_types", Requirements()[1]["file_types"]
    )
    Irange = settings_class.output_settings.get("Irange", Requirements()[1]["Irange"])
    # override with kwargs
    pattern = kwargs.get("pattern", None)
    search_parameter = kwargs.get("search_parameter","")
    search_over = kwargs.get("search_over", None)
    fps = kwargs.get("fps", fps)
    file_types = kwargs.get("file_types", file_types)
    Irange = kwargs.get("Irange", Irange)
    
    # if the kwargs are not set, take a guess
    if not search_parameter:
        search_parameter = "" 
    if settings_class.image_number == 1:
        pattern = 0            
    elif pattern:
        settings_class.set_data_files(keep=pattern)
        pattern = 0
    else: # not pattern
    # if not pattern:# or not search_parameter or not search_over:
        fls = glob.glob(f"./{settings_class.output_directory}/*scan*{search_parameter}*.json")
        if len(fls) == 0:
            raise ValueError("There is no identified search file to plot.")
        tm = []
        for i in range(len(fls)):
            tm.append(os.path.getmtime(fls[i]))
        latest = np.argsort(tm)[-1]
        split = os.path.splitext(os.path.basename(fls[latest]))[0].split("__")
        # get file name
        for f_ind in range(len(settings_class.datafile_list)):
            if split[0] in str(settings_class.datafile_list[f_ind]):
                pattern = f_ind
        if settings_class.datafile_number > 1:
            # search over the first file only
            settings_class.set_data_files(keep=pattern) 
                
    if not search_parameter or search_parameter == "": #or not search_over:
        # get pattern name
        pattern_name = os.path.splitext(os.path.basename(settings_class.image_list[pattern]))[0]
        # use filenames for seaech paramter
        fls = glob.glob(f"./{settings_class.output_directory}/*{pattern_name}*scan*{search_parameter}*.json")       
        search_parameter = os.path.splitext(os.path.basename(fls[0]))[0].split("__")[1].split("=")[1]
        
    if not search_over:
        # get pattern name
        pattern_name = os.path.splitext(os.path.basename(settings_class.image_list[pattern]))[0]
        # use filenames for seaech paramter
        fls = glob.glob(f"./{settings_class.output_directory}/*{pattern_name}*scan*{search_parameter}*.json")   
        values = []
        for v in fls:
            values.append(int(os.path.splitext(os.path.basename(v))[0].split("__")[2].split("=")[1]))
        search_over = [np.min(values), np.max(values)]

    # add search values as metadata so that it can be read later. 
    if search_parameter not in settings_class.metadata:
        settings_class.metadata.append(search_parameter)
    settings_class.metadata_settings = {}

    # if len(search_over) == 2 and search_over[0] == 0:
    #     search_over[0] = 1
    search = np.arange(search_over[0],search_over[1])
    
    # get all data - read from separate json files.
    df = pd.DataFrame()
    data_fit_all = []
    for srch in search:

        settings_class_reduce = settings_class.duplicate()
        if settings_class_reduce in settings_class.__dict__:
            setattr(settings_class_reduce, search_parameter, srch) 
            
        settings_class_reduce.fit_propagate = False
        settings_class_reduce.file_label = (
            "scan="
            + search_parameter
            + "__"
            + "value="
            + str(srch)
        ) 
        
        # add search values as metadata so that it can be read later. 
        settings_class_reduce.metadata_settings[search_parameter] = srch
        
        # circulment revalidating the settings class
        settings_class_reduce._unmodified_self = settings_class_reduce._validation_copy()  

        df_tmp = ReadFits_to_dataframe(
            settings=settings_class_reduce,
            includeStats=True,
            includeSeriesValues=True,
            includeIntensityRanges=True,
            includePosition=True,
            includeIntegrated=False
        )
        df = pd.concat([df,df_tmp])    
        
        # FIXME: ReadFits_to_list burries data_fit in an unnecessary list -- hance [0] two lines down
        data_fit, metadata = ReadFits_to_list(settings_class_reduce)
        data_fit_all.append(data_fit[0])
    headers = list(df_tmp.columns.values)

    # make the data class.
    data_to_fill = settings_class.image_list[pattern]
    data_class = settings_class.data_class
    data_class.reduce_by = None
    settings_class.reduce_by = None
    data_class.fill_data(
        data_to_fill,
        settings=settings_class,
    )
       
    ranges = df["range_position"].unique()
    
    #make movies
    for i in range(len(ranges)):
        # loop over the number of ranges/subpaterns

        # use pos_in_range and _search peak to restrict values to the peak that
        # is being searched over -- for multiple peaks in range
        df_peak = df[
            (df["range_position"] == ranges[i])
        ]

        # sort the data into a sensible order
        df_peak = df_peak.sort_values([search_parameter])

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
        for j in range(len(data_fit_all)):# ddf_peak[search_parameter]:
            
            settings_class.set_subpattern(pattern, i)     

            data_copy = deepcopy(data_class)
            data_copy.reduce_by = df_peak[search_parameter].iloc[j]
           
            if data_copy.reduce_by is not None:
                data_copy.intensity = data_copy._reduce_array(data_copy.intensity)
                data_copy.tth = data_copy._reduce_array(data_copy.tth)
                data_copy.azm = data_copy._reduce_array(data_copy.azm, polar=True)

            data_copy = data_copy.duplicate_without_detector(
                range_bounds=[
                    df_peak["range_start"].iloc[0],
                    df_peak["range_end"].iloc[0],
                ]
            )

            # Mask the subpattern by intensity if called for
            if (
                "imax" in settings_class.subfit_orders
                or "imin" in settings_class.subfit_orders
            ):
                data_copy = SpotProcess(data_copy, settings_class)
       
            # make the plot of the fits.
            fig = plot_FitAndModel(
                settings_class,
                data_copy,
                params_dict = data_fit_all[j][i],
                figure=fig,
                plot_type="scatter",
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
                + search_parameter
                + " = "
                + str(df_peak[search_parameter].iloc[j])
            )
            # if "note" in settings_class.subfit_orders:
            #     title_str += " " + settings_class.subfit_orders["note"]
            plt.suptitle(title_str)
            figure_suptitle_space(fig, topmargin=0.4)
            
            # make the video clip
            # just addes the figure as a frame to the proto-video.
            frames.append(ImageClip(mplfig_to_npimage(fig)).with_duration(1))

        # convert to video and write
        video = concatenate_videoclips(frames, method="compose")
        for f in range(len(file_types)):
            settings_class.file_label = "scan=" + search_parameter
            # make videofile name
            data_fit_tmp = data_fit_all[j][i]
            if "note" in data_fit_tmp:
                data_fit_tmp.pop("note")
            out_file = make_outfile_name(
                settings_class.subfit_filename,
                directory=settings_class.output_directory,
                orders=data_fit_tmp,
                additional_text=settings_class.file_label,
                extension=file_types[f],
                peak="all",
                overwrite=True,
            )
            video.write_videofile(out_file, fps=fps)
