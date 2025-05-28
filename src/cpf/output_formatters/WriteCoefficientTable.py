__all__ = ["Requirements", "WriteOutput"]


import json
import os
from itertools import product
import re

import numpy as np
import pandas as pd

import cpf.peak_functions as pf
from cpf.output_formatters.ReadFits import ReadFits
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
        "coefs_vals_write"  # -- pick which set of coefficients to write
    ]

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(
    settings_class=None,
    settings_file=None,
    fitStats=True,
    writefile=True,
    *args,
    **kwargs,
):
    """
    Write coefficents from fits to table/csv file. 
    
    :param settings_class: DESCRIPTION, defaults to None
    :type settings_class: TYPE, optional
    :param settings_file: DESCRIPTION, defaults to None
    :type settings_file: TYPE, optional
    :param fitStats: DESCRIPTION, defaults to True
    :type fitStats: TYPE, optional
    :param writefile: DESCRIPTION, defaults to True
    :type writefile: TYPE, optional
    :param *args: DESCRIPTION
    :type *args: TYPE
    :param **kwargs: DESCRIPTION
    :type **kwargs: TYPE
    :raises ValueError: DESCRIPTION
    :return: DESCRIPTION
    :rtype: TYPE

    """

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
    
    # read the data.
    df = ReadFits(settings_class=settings_class, fitStats=fitStats)
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
        additional_text="all_coefficients",
    )

    # write file using panda dataframe    
    if writefile == True:
        ## format dateframe for writing to file neatly.
        # make strings in DateFile and Peak columns all the same length
        len_datafile = np.max(df["DataFile"].str.len())
        len_peaks = np.max(df["Peak"].str.len())
        df_tmp = df["DataFile"].str.pad(
            np.max([len_datafile, col_width]), side="left", fillchar=" "
        )
        df["DataFile"] = df_tmp
        df_tmp = df["Peak"].str.pad(
            np.max([len_peaks, col_width]), side="left", fillchar=" "
        )
        df["Peak"] = df_tmp

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
            elif headers[i] == "Peak":
                setattr(
                    columns,
                    headers[i],
                    headers[i].rjust(np.max([len_peaks, col_width])),
                )
            else:
                setattr(columns, headers[i], headers[i].rjust(col_width))
        columns = columns.__dict__
        df.rename(columns=columns, inplace=True)

        # write data frame to csv file
        with open(out_file, "w") as f:
            logger.info(" ".join(map(str, [("Writing %s" % out_file)])))

            f.write(
                "# continuous_peak_fit : Table of all coefficients from fits for input file: %s.\n"
                % settings_class.settings_file
            )
            f.write("# For more information: http://www.github.com/me/something\n")
            f.write("# File version: %i \n" % version)
            f.write("# \n")

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
