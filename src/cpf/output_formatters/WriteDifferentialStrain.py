__all__ = ["Requirements", "WriteOutput"]


import json
import os
import re
from itertools import product

import numpy as np
import pandas as pd

import cpf.peak_functions as pf
from cpf.output_formatters.fits_io import ReadFits_to_dataframe
from cpf.output_formatters.output_csv import make_header, write_csv
from cpf.settings import get_settings
from cpf.util.io import make_outfile_name
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteCoefficientTable")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        ##"Output_directory"  # if no direcrtory is specified write to current directory.
        "SampleGeometry": "3d",  # -- geometry of the sample for determining the cnetres from. 2D or 3D.
        "SampleDeformation": "compression",  # changes calculation between 'compression' and 'extension'.
        "dp": 5,  # how many decimal points to write out
        "col_width": 15,  # default column width for csv file.
        "ordering_of_output": False,  # Just leave as read -- otherwise list of dataframe headers to order by
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
    SampleGeometry = settings_class.output_settings.get(
        "SampleGeometry", Requirements()[1]["SampleGeometry"]
    )
    SampleDeformation = settings_class.output_settings.get(
        "SampleDeformation", Requirements()[1]["SampleDeformation"]
    )
    dp = settings_class.output_settings.get("dp", Requirements()[1]["dp"])
    col_width = settings_class.output_settings.get(
        "col_width", Requirements()[1]["col_width"]
    )
    ordering_of_output = settings_class.output_settings.get(
        "ordering_of_output", Requirements()[1]["ordering_of_output"]
    )
    # override with kwargs
    SampleGeometry = kwargs.get("SampleGeometry", SampleGeometry)
    SampleDeformation = kwargs.get("SampleDeformation", SampleDeformation)
    dp = kwargs.get("dp", dp)
    col_width = kwargs.get("col_width", col_width)
    ordering_of_output = kwargs.get("ordering_of_output", ordering_of_output)
    # force all the kwargs that might be needed
    set_params = {
        "SampleGeometry": SampleGeometry,
        "SampleDeformation": SampleDeformation,
    }
    kwargs.update(set_params)

    # read the data.
    df = ReadFits_to_dataframe(
        settings=settings_class,
        includeSeriesValues=True,
        includeStats=fitStats,
        includeIntegrated=True,
        SampleGeometry=SampleGeometry,
        SampleDeformation=SampleDeformation,
    )
    headers = list(df.columns.values)
    # order the rows
    if ordering_of_output:
        df = df.sort_values(by=ordering_of_output)

    # limit dataframe to what we want to write.
    # order columns to be correct also
    headers_use =  ["num",
                   "DataFile",
                   'phase', 
                   'peak',   
                   ]
    headers_use += settings.metadata
    headers_use += ["d_mean",
                    "d_mean_err",  
                    'd-space4',
                    'd-space4_err',
                    'd-space3', 
                    'd-space3_err',
            # differential components
                   'differential',
                   'differential_err',
                   'orientation',  
                   'orientation_err', 
                   'd_max',     
                   'd_min',      
                   'height mean',   
                   'height mean err',
                   'width mean',   
                   'width mean err',
                   'profile mean',  
                   'profile mean err'
                   ]
    additional = list(set(pf.peak_components(include_combined=True)[1]) - set(pf.peak_components(include_combined=False)[1]))
    for add in additional:
        headers_use += [add + " mean",
                        add + " mean err",
                       ]
    headers_rename = {'d-space4':"d2cos",
                'd-space4_err':"d2cos_err" ,
                'd-space3':"d2sin",
                'd-space3_err':"d2sin_err"}
    if fitStats == True:
        extra_headers = [
                   'time-elapsed',
                   'chunks-time',
                   'sum-residuals-squared',
                   'status',
                   'function-evaluations',
                   'n-variables',
                   'n-data',
                   'degree-of-freedom',
                   'ChiSq',
                   'RedChiSq',
                   'aic',
                   'bic'
                   ]
        headers_use += extra_headers

    # check that all wanted values are present
    missing_headers = list(set(headers_use) - set(list(df)))
    # if missing headers likely means no differential strains. Add list dataframe as zeros
    for i in missing_headers:
        df[i] = 0
    df = df.loc[:, headers_use]

    # cut data frame
    df = df[headers_use]
    # rename the columns
    df.rename(columns=headers_rename, inplace=True)

    # make filename for output
    base = settings_class.datafile_basename
    if base is None:
        logger.info(
            " ".join(map(str, [("No base filename, using input filename instead.")]))
        )
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    if settings_class.file_label:
        add_text = settings_class.file_label
    else:
        add_text = ""
    add_text += "DifferentialStrains"
    out_file = make_outfile_name(
        base,
        directory=settings_class.output_directory,
        extension=".dat",
        overwrite=True,
        additional_text=add_text,
    )

    ## outfile header
    calc_options = {}
    calc_options["Sample Geometry"] = SampleGeometry
    calc_options["Sample Deformation"] = SampleDeformation
    file_header = make_header(
        settings_class,
        derived="Differnetial strains",
        calc_options=calc_options,
    )

    # write file using panda dataframe
    write_csv(out_file, df, headers, file_header, col_width=col_width, dp=dp)
