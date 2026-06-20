__all__ = ["Requirements", "WriteOutput"]


import os

import glob
import numpy as np

from typing import Literal, Optional

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from  cpf.settings import get_settings
from cpf.output_formatters.fits_io import ReadFits_to_dataframe, ReadFits_to_list
from cpf.IO_functions import make_outfile_name, peak_string
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteCoefficientTable")



def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        ##"Output_directory"  # if no direcrtory is specified write to current directory.
        #"coefs_vals_write": "given"  # -- pick which set of coefficients to write
    }

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(
    settings,
    file_label = None,
    statistic: Literal["bic", "aic", "RedChiSq", "ChiSq"] = "bic",
    key_parameter = ["d-space0", "differential"],
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
    
    if file_label is not None:
        settings_class.file_label = file_label
    
    #this is search data so there is a postscript in the json file label.
    #this is search data so there is a postscript in the json file label.
    # determine the label
    if "file_label" not in dir(settings_class) or settings_class.file_label is None:
        fls = glob.glob(f"./{settings_class.output_directory}/*search*.json")
        if len(fls) == 0:
            raise ValueError("There is no identified search file to plot.")
        else:
            tm = []
            for i in range(len(fls)):
                tm.append(os.path.getmtime(fls[i]))
            latest = np.argsort(tm)[-1]
            settings_class.file_label = os.path.splitext(os.path.basename(fls[latest]))[0].split("__")[1]
            
            # make sure that we have the right file number for this set.
            possible = []
            for i, subval in enumerate(settings_class.image_list):
                if isinstance(subval, list):
                    #subval = subval[0]
                    settings_class.set_subpattern(i, 0)
                    subval = make_outfile_name(
                        settings_class.subfit_filename,
                        directory=None,
                        overwrite=True,
                    )
                if os.path.splitext(os.path.basename(fls[latest]))[0].split("__")[0] in subval:
                    possible.append(i)
            if len(possible) != 1:
                raise ValueError("There is no identified search file to plot.")
            else:
                searchdata = possible[0]

        #restrict to just the required image
        settings_class.set_data_files(keep=searchdata)
    
    # read the data.
    df = ReadFits_to_dataframe(settings=settings_class, includeStats=True, includeSeriesValues=True, includePosition=True, IncludeIntegrated=False)
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

    df["Fit_time_without_chunks"] = df["time-elapsed"]-df["chunks-time"]
    df["RedChiSq_per_s"] = df["RedChiSq"]/df["Fit_time_without_chunks"]
    
    peaks = df["peak"].unique()
    searches = df["series_type"].unique()
    param_plot = ["RedChiSq", "Fit_time_without_chunks", "RedChiSq_per_s", "bic", "aic", "d-space0", "time-elapsed", "chunks-time"]
    
    param_plot = [statistic, "Fit_time_without_chunks"] + key_parameter
    cols=int(2)
    fig_scale = 1.5

    for i in range(len(peaks)):
        fig, ax = plt.subplots(int(np.ceil(len(param_plot))/cols), cols, sharex=True, figsize=[8*fig_scale,6*fig_scale])
        ax = ax.flat
        title_str = ("Order Search - " + 
            peak_string(settings_class.subfit_orders, peak=i) +
            " - " + df["search_over"][0]
        )
        fig.suptitle(title_str)
        for h in range(len(param_plot)):
            for j in range(len(searches)):
                
                df_tmp = df[(df['peak'] == peaks[i]) & (df['series_type'] == searches[j])]
                # if there is more than 1 peak we need to filter the data to 
                # just look at the one that has been searched over.
                if len(df_tmp['search_peak'].unique()) == 1:
                    # there is only 1 peak. has to be 0.
                    k = 0
                else:
                    k= df_tmp['pos_in_range'].unique()[0]    
                if param_plot[h]+"_err" in headers:
                    ax[h].errorbar(df_tmp.loc[df_tmp['search_peak'] == k]['search_value'], 
                                   df_tmp.loc[df_tmp['search_peak'] ==k][param_plot[h]], 
                                   yerr =  df_tmp.loc[df_tmp['search_peak'] ==k][param_plot[h]+"_err"],
                                   fmt='.-', capsize=5, label=searches[j])
                else:
                    ax[h].plot(df_tmp.loc[df_tmp['search_peak'] == k]['search_value'], df_tmp.loc[df_tmp['search_peak'] ==k][param_plot[h]], '.-', label=searches[j])
        
            ax[h].set_xlabel(f"{df['search_over'][0]} order")
            ax[h].xaxis.set_major_locator(MaxNLocator(integer=True))
            ax[h].set_ylabel(param_plot[h])
            # ax[h].set_title(peaks[i])
            ax[h].legend()
        
        # write figures to file
        # position = df.iloc[i]["Position in json"]
        # data_fit_tmp = data_fit[position]
        # if "note" in data_fit_tmp:
        #     data_fit_tmp.pop("note")
        lbl = settings_class.file_label
        if lbl[-3:] == "all":
            lbl = lbl[:-3]+str(i)
        lbl = peak_string(settings_class.subfit_orders, peak=i, fname=True) + "__" + lbl
        out_file = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            # orders = data_fit_tmp,
            additional_text=lbl,
            extension='.png',
            peak=peaks[i],
            overwrite=True,
        )
        fig.savefig(out_file)
            
    return df