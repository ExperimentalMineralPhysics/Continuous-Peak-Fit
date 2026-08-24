__all__ = ["Requirements", "WriteOutput"]


import glob
import os
from typing import Literal, Optional

import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from cpf.output_formatters.fits_io import ReadFits_to_dataframe, ReadFits_to_list
from cpf.settings import get_settings
from cpf.util.io import make_outfile_name, peak_string
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteCoefficientTable")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        "pattern": None, # guess from files in outpuct directory
        "search_parameter": "", # guess from files in outpuct directory
        "search_over": None, # guess from files in outpuct directory
    }

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(
    settings,
    statistic: Literal["bic", "aic", "RedChiSq", "ChiSq"] = "bic",
    key_parameter=["d-space0", "differential"],
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    *args,
    **kwargs,
):
    """
    Plots outputs of cpf.XRD_FitPattern.order_search.
    Outputs are an indication of what is the best order to use for a fit.

    Parameters
    ----------
    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or
        dictionary of all the settings or
        string or path to a file with the settings in.
    file_label : string, optional
        Additional text in json file name added by cpf.XRD_FitPattern.order_search().
        If not present the default is to use the newest file. The default is None.
    *args : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    df : panda dataframe
        Dateframe containing all the parameters from the fits.
    """

    # make sure settings is a class
    settings_class = get_settings(settings)
    
    # Parse optional parameters
    pattern = settings_class.output_settings.get("pattern", Requirements()[1]["pattern"])
    search_parameter = settings_class.output_settings.get("search_parameter",Requirements()[1]["search_parameter"])
    search_over = settings_class.output_settings.get("search_over", Requirements()[1]["search_over"])
    # override with kwargs
    pattern = kwargs.get("pattern", None)
    search_parameter = kwargs.get("search_parameter","")
    search_over = kwargs.get("search_over", None)
    
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

    # get all data
    df = pd.DataFrame()
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
            includePosition=True,
            includeIntegrated=False
        )
        df = pd.concat([df,df_tmp])
    headers = list(df_tmp.columns.values)


    df["Fit_time_without_chunks"] = df["time-elapsed"] - df["chunks-time"]
    df["RedChiSq_per_s"] = df["RedChiSq"] / df["Fit_time_without_chunks"]
    ranges = df["range_position"].unique()
    param_plot = [
        "RedChiSq",
        "Fit_time_without_chunks",
        "RedChiSq_per_s",
        "bic",
        "aic",
        "d-space0",
        "time-elapsed",
        "chunks-time",
    ]

    param_plot = [statistic, "Fit_time_without_chunks"] + key_parameter
    cols = int(2)
    fig_scale = 1.5

    for i in range(len(ranges)):
        
        settings_class.set_subpattern(0, i)
        fig, ax = plt.subplots(
            int(np.ceil(len(param_plot)) / cols),
            cols,
            sharex=True,
            figsize=[8 * fig_scale, 6 * fig_scale],
        )
        ax = ax.flat
        title_str = (
            "Order Search - "
            + peak_string(settings_class.subfit_orders)
            + " - "
            + search_parameter
        )
        fig.suptitle(title_str)
        for h in range(len(param_plot)):
            
            df_tmp = df[
                (df["range_position"] == ranges[i])
            ]
                        
            if param_plot[h] + "_err" in headers:
                for p in df_tmp["pos_in_range"].unique():
                    ax[h].errorbar(
                        df_tmp[(df_tmp["pos_in_range"] == p)][search_parameter],
                        df_tmp[(df_tmp["pos_in_range"] == p)][param_plot[h]],
                        yerr=df_tmp[
                            param_plot[h] + "_err"
                        ][(df_tmp["pos_in_range"] == p)],
                        fmt=".-",
                        capsize=5,
                        label=f"{df_tmp['phase'].iloc[p]} ({df_tmp['peak'].iloc[p]})",
                    )
            else:
                ax[h].plot(
                    df_tmp[search_parameter],
                    df_tmp[param_plot[h]],
                    ".-",
                    label=search_parameter,
                )
                
            if (h==0  
                and np.log(np.nanmax(df_tmp[param_plot[h]])-np.nanmin(df_tmp[param_plot[h]]) ) >1
                and np.nanmin(df_tmp[param_plot[h]]) > 0
                ):
                ax[h].set_yscale("log", nonpositive='clip')
            else:
                ax[h].set_yscale("linear")

            ax[h].set_xlabel(f"{search_parameter}")
            ax[h].xaxis.set_major_locator(MaxNLocator(integer=True))
            ax[h].set_ylabel(param_plot[h])
            # ax[h].set_title(peaks[i])
            ax[h].legend()

        # write figures to file
        # position = df.iloc[i]["Position in json"]
        # data_fit_tmp = data_fit[position]
        # if "note" in data_fit_tmp:
        #     data_fit_tmp.pop("note")
        
        settings_class.file_label = (
            "scan="
            + search_parameter
        ) 
    
        lbl = settings_class.file_label
        if lbl[-3:] == "all":
            lbl = lbl[:-3] + str(i)
        settings_class.set_subpattern(pattern, i)
        lbl = peak_string(settings_class.subfit_orders, fname=True) + "__" + lbl
        out_file = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            # orders = data_fit_tmp,
            additional_text=lbl,
            extension=".png",
            peak=ranges[i],
            overwrite=True,
        )
        fig.savefig(out_file)

    return df
