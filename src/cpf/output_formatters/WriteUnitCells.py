__all__ = ["Requirements", "WriteOutput"]


import json
import os
from itertools import product
import re

import numpy as np
import pandas as pd

import cpf.peak_functions as pf
from cpf.output_formatters.convert_fit_to_unitcell import fits_to_unitcell
from cpf.IO_functions import make_outfile_name
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteCoefficientTable")



def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        ##"Output_directory"  # if no direcrtory is specified write to current directory.
        "reflections_to_use"  # -- pick which set of reflections to use for unit cell volume
        "phase" # -- pick whick phase to fit unit cell for.
        "SampleGeometry" # -- geometry of the sample for determining the cnetres from. 2D or 3D.
        "weighted" # -- weighted fit or not. True/False
    ]

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(
    settings_class=None,
    settings_file=None,
    fitStats=True,
    *args,
    **kwargs,
):
    """
    Write unit-cell volumes derived from fitted peak centroids. Writes the values 
    to table/csv file. 
    
    :param settings_class: DESCRIPTION, defaults to None
    :type settings_class: TYPE, optional
    :param settings_file: DESCRIPTION, defaults to None
    :type settings_file: TYPE, optional
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


    # force all the kwargs that might be needed
    kwargs.pop("temperature", np.nan) # supress calculations of pressure in outputs

    # define defaults
    dp = 6  # how many decimal points to write out
    col_width = 15  # default column width for csv file.

    ## output file version
    # 1: ?
    # 2: ?
    # 3: rewritten as panda data frame - to force columns to line up.
    version = 3

    ordering_of_output = "peak"

    if settings_class is None and settings_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        from cpf.XRD_FitPattern import initiate

        settings_class = initiate(settings_file)
    
    # get the unit cells from settings.
    df = fits_to_unitcell(settings_file=settings_file,
                includeStats=False,
                includeParameters=False,
                includeSeriesValues = False,
                includeIntensityRanges = False,
                includeUnitCells = True,
                includePosition = False,
                **kwargs
                )
    
    headers = list(df.columns.values)

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
        additional_text="unit_cells",
    )

    # write file using panda dataframe    

    ## format dateframe for writing to file neatly.
    # make strings in DateFile and Peak columns all the same length
    len_datafile = np.max(df["DataFile"].str.len())
    # len_peaks = np.max(df["Peak"].str.len())
    df_tmp = df["DataFile"].str.pad(
        np.max([len_datafile, col_width]), side="left", fillchar=" "
    )
    df["DataFile"] = df_tmp
    # df_tmp = df["Peak"].str.pad(
    #     np.max([len_peaks, col_width]), side="left", fillchar=" "
    # )
    # df["Peak"] = df_tmp

    # rename the columns so that the headers are the same width as the columns
    class NewClass(object):
        pass

    columns = NewClass()
    for i in range(len(headers)):
        if headers[i] == "DataFile":
            setattr(
                columns,
                headers[i],
                headers[i].rjust(np.max([len_datafile, col_width])),
            )
        else:
            setattr(columns, headers[i], headers[i].rjust(col_width))
    columns = columns.__dict__
    df.rename(columns=columns, inplace=True)

    # make sure residual columns are saved as a single string with no line breaks.
    colms = [col for col in df.columns if 'residuals' in col]
    for i in colms:
        df[i] = df[i].apply(lambda x: np.array2string(x, separator=";", max_line_width=np.inf, formatter={"float_kind": lambda x: float_to_string_formatter(x, dp, 10) }, sign=" "))
        # df[i] = df[i].apply(lambda x: np.array2string(x, separator=";", max_line_width=np.inf, formatter={"float_kind": lambda x: f"{x:{str(col_width)}.{str(dp)}f}" if np.abs(x) > 0.1 else f"{x:{str(col_width)}.{str(dp)}e}"}))
       
    #remove hkls from data frame -- write as a header instead
    cols = [col for col in df.columns if 'hkl' in col]
    hkls = {}
    for i in cols:
        hkls[i] = df[i].iloc[0]
        df = df.drop(i, axis=1)
        
        # df[i] = df[i].apply(lambda x: np.array2string(x, separator=";", max_line_width=np.inf))

    # write data frame to csv file
    with open(out_file, "w") as f:
        logger.info(" ".join(map(str, [("Writing %s" % out_file)])))

        f.write(
            "# continuous_peak_fit : Table of unit cells from fits for input file: %s.\n"
            % settings_class.settings_file
        )
        f.write("# For more information: http://www.github.com/me/something\n")
        f.write("# File version: %i \n" % version)
        f.write("# \n")

        if len(hkls) >= 1:
            for i in list(hkls):
                f.write("# Peaks used in volume calculations: \n")
                f.write(f"#    {i.replace("hkls","").replace("hkl","").strip()} : ")
                for j in hkls[i]:
                    f.write(f"({j}) ")
                f.write("\n# \n")     

        df.to_csv(
            f,
            index=False,
            header=True,
            na_rep="".ljust(col_width),
            float_format=lambda x: f"{x:{str(col_width)}.{str(dp)}f}"
                    if np.abs(x) > 0.1
                    else f"{x:{str(col_width)}.{str(dp)}e}"
                )
        
    # rewrite the file adjusting the column widths to keep the data lined up. 
    with open(out_file, 'r') as fl: 
        in_lines = fl.readlines() 
    with open(out_file, 'w') as f:
        for line in in_lines:
            split_line = line.split(",")
            split_line_out = []
            running_length = 0
            expected_length = 0
            for i in range(len(split_line)):
                if i==0: # first column. keep narrow
                    col_here = 3
                    if "num" in split_line[i]:
                        split_line_out.append(f"{split_line[i].replace(' ',''):>{col_here}}")
                    else:
                        split_line_out.append(f"{split_line[i]:>{col_here}}")
                elif i==1: # file names
                    col_here = col_width
                    split_line_out.append(f" {split_line[i]:>{col_here}}")
                else: # data values. adjust column width to line everything up.
                    col_here = np.max([0,col_width - (running_length-expected_length)])
                    split_line_out.append(f"{split_line[i].replace(' ',''):>{col_here}}")
                    running_length += len(split_line_out[-1])
                    expected_length += col_width

            out_line = ",".join(split_line_out)
            f.write(out_line)


def float_to_string_formatter(x, dp=5, col_width=10):
    if np.abs(x) > 0.1 and np.log10(np.abs(x)) < col_width-dp-4:
        if x > 0:
            x = f" {x:{str(col_width)}.{str(dp)}f}" 
        else:
            x = f"{x:{str(col_width)}.{str(dp)}f}" 
    # elif np.log10(np.abs(x)) > col_width-dp-1:
    #     x = f"{x:{str(col_width)}.{str(dp)}e}"
    else:
        if x > 0:
            x = f" {x:{str(col_width)}.{str(dp)}e}"
        else:
            x = f"{x:{str(col_width)}.{str(dp)}e}"
    return x
    
    