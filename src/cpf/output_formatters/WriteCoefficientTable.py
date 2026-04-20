__all__ = ["Requirements", "WriteOutput"]


import json
import os
from itertools import product
import re

import numpy as np
import pandas as pd

import cpf.peak_functions as pf
from  cpf.settings import get_settings
from cpf.output_formatters.fits_io import ReadFits_to_dataframe
from cpf.output_formatters.output_csv import write_csv, make_header
from cpf.IO_functions import make_outfile_name
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteCoefficientTable")



def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        ##"Output_directory"  # if no direcrtory is specified write to current directory.
        "dp": 6,  # how many decimal points to write out
        "col_width": 15,  # default column width for csv file.
        "coefs_vals_write": "all",  # -- pick which set of coefficients to write
        "ordering_of_output": None # Just leave as read -- otherwise list of dataframe headers to order by
    }

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(
    settings,
    fitStats=True,
    *args,
    **kwargs,
):
    """
    Write coefficents from fits to table/csv file. 
    
    Parameters
    ----------
    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or 
        dictionary of all the settings or 
        string or path to a file with the settings in.
    fitStats : bool, optional
        switch to include all the fit stats in the output file. The default is True.
    *args : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    """

    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters 
    dp               = settings_class.output_settings.get("dp", Requirements()[1]["dp"])
    col_width        = settings_class.output_settings.get("col_width", Requirements()[1]["col_width"])
    coefs_vals_write = settings_class.output_settings.get("coefs_vals_write", Requirements()[1]["coefs_vals_write"])
    ordering_of_output = settings_class.output_settings.get("ordering_of_output", Requirements()[1]["ordering_of_output"])
    #override with kwargs
    dp               = kwargs.get("dp", dp)
    col_width        = kwargs.get("col_width", col_width)
    coefs_vals_write = kwargs.get("coefs_vals_write", coefs_vals_write)
    ordering_of_output = kwargs.get("ordering_of_output", ordering_of_output)

    # read the data.
    df = ReadFits_to_dataframe(settings=settings_class, fitStats=fitStats)
    headers = list(df.columns.values)
    # cut data frame
    if coefs_vals_write != "all":
        df = df[coefs_vals_write]
    #order the rows
    if ordering_of_output:
        df = df.sort_values(by=ordering_of_output) 

    # make filename for output
    base = settings_class.datafile_basename
    if base is None:
        logger.info(
            " ".join(map(str, [("No base filename, using input filename instead.")]))
        )
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    out_file = make_outfile_name(
        base,
        directory=settings_class.output_directory,
        extension=".dat",
        overwrite=True,
        additional_text="all_coefficients",
    )
    
    ## outfile header
    file_header = make_header(settings_class,
                            fits="Fit coefficients", 
                            calc_options=None
                            )
    
    # write file using panda dataframe
    write_csv(out_file, df, headers, file_header, col_width=col_width, dp=dp)