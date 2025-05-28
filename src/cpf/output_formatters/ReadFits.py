#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ["ReadFits"]


import json
from itertools import product
import re

import numpy as np
import pandas as pd

import cpf.peak_functions as pf
from cpf.IO_functions import make_outfile_name
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.ReadFits")


def ReadFits(
    settings_class=None,
    settings_file=None,
    fitStats=True,
    *args,
    **kwargs,
):
    """
    Read coefficents from fits and return values as panda dataframe. 
    
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

    if settings_class is None and settings_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        from cpf.XRD_FitPattern import initiate

        settings_class = initiate(settings_file)

    # get what to write
    if "coefs_get" in kwargs:
        coefs_vals_get = kwargs.get("coefs_get")
    elif "coefs_vals_get" in settings_class.output_settings:
        coefs_vals_get = settings_class.output_settings["coefs_vals_get"]
    else:
        coefs_vals_get = "all"
    if isinstance(coefs_vals_get, str):
        coefs_vals_get = [coefs_vals_get]
    if coefs_vals_get == ["all"]:
        peak_properties = pf.peak_components(full=True)
        coefs_vals_get = peak_properties[1]

    # get number of coefficients.
    num_subpatterns = len(settings_class.fit_orders)
    max_coef = {}
    for w in range(len(coefs_vals_get)):
        ind = coefs_vals_get[w]
        if ind != "background":
            max_coef[coefs_vals_get[w]] = 0
        else:
            max_coef[coefs_vals_get[w]] = [0]
    for y in range(num_subpatterns):
        for x in range(len(settings_class.fit_orders[y]["peak"])):
            for w in range(len(coefs_vals_get)):
                ind = coefs_vals_get[w]
                if ind != "background" and ind != "symmetry":
                    if isinstance(settings_class.fit_orders[y]["peak"][x][ind], list):
                        coef_tmp = [
                            x
                            for x in settings_class.fit_orders[y]["peak"][x][ind]
                            if type(x) == int
                        ]
                    else:
                        coef_tmp = settings_class.fit_orders[y]["peak"][x][ind]
                    max_coef[ind] = np.max([max_coef[ind], 2 * np.max(coef_tmp) + 1])
                elif ind == "symmetry":
                    # symmetry might not be present in the fitting files.
                    try:
                        max_coef[ind] = np.max(
                            [
                                max_coef[ind],
                                np.max(2 * settings_class.fit_orders[y]["peak"][x][ind])
                                + 1,
                            ]
                        )
                    except:
                        pass
                else:
                    for v in range(len(settings_class.fit_orders[y][ind])):
                        if v > len(max_coef[ind]) - 1:
                            max_coef[ind].append(0)
                        max_coef[ind][v] = np.max(
                            [
                                max_coef[ind][v],
                                2 * settings_class.fit_orders[y][ind][v] + 1,
                            ]
                        )

    # store headers
    headers = []
    headers.append("num")
    headers.append("DataFile")
    headers.append("Peak")
    headers.append("Range_start")
    headers.append("Range_end")

    # parmeter header list
    for w in range(len(max_coef)):
        ind = coefs_vals_get[w]
        if ind != "background" and ind != "symmetry":
            for v in range(max_coef[ind]):
                headers.append(ind + str(v))
                headers.append(ind + str(v) + "err")
        elif ind == "symmetry":
            headers.append(ind)
        else:
            for u in range(len(settings_class.fit_orders[y][ind])):
                for v in range(max_coef[ind][u]):
                    headers.append(ind + str(u) + "_" + str(v))
                    headers.append(ind + str(u) + "_" + str(v) + "err")

    # include properties from the lmfit output that were passed with the fits.
    if fitStats is True:
        #read the list of parameters from the first file.
        settings_class.set_subpattern(0, 0)
        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            extension=".json",
            overwrite=True,
        )  # overwrite = True to get the file name without incrementing it.
        with open(filename) as json_data:
            tmp = json.load(json_data)
            properties = list(tmp[0]["FitProperties"].keys())
            headers += properties

    # read all the data.
    fits = []
    for z in range(settings_class.image_number):
        settings_class.set_subpattern(z, 0)

        filename = make_outfile_name(
            settings_class.subfit_filename,  # diff_files[z],
            directory=settings_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".json",
            overwrite=True,
        )  # overwrite =false to get the file name without incrlemeting it.

        # Read JSON data from file
        with open(filename) as json_data:
            fits.append(json.load(json_data))

    # make lists of the parameters to iterate over
    images = list(range(settings_class.image_number))
    subpatterns = list(range(num_subpatterns))
    max_peaks = 1
    for i in range(len(settings_class.fit_orders)):
        max_peaks = np.max([max_peaks, len(settings_class.fit_orders[i]["peak"])])
    lists = images, subpatterns, list(range(max_peaks))
    lists = np.array(list(product(*lists)))

    # sort lists by peaks
    lists = lists[np.lexsort((lists[:, 2],))]
    lists = lists[np.lexsort((lists[:, 1],))]

    RowsList = []
    # write the data to array, then append to dataframe
    for z in range(len(lists)):
        RowLst = {}
        settings_class.set_subpattern(lists[z, 0], lists[z, 1])
        data_to_write = fits[lists[z, 0]][lists[z, 1]]

        if len(data_to_write["peak"]) > lists[z, 2]:
            
            RowLst["num"] = lists[z, 0]
            
            out_name = make_outfile_name(
                settings_class.subfit_filename,
                directory="",
                extension=".json",
                overwrite=True,
            )

            out_peak = []

            # peak
            if "phase" in settings_class.fit_orders[lists[z, 1]]["peak"][lists[z, 2]]:
                out_peak = settings_class.fit_orders[lists[z, 1]]["peak"][lists[z, 2]][
                    "phase"
                ]
            else:
                out_peak = out_peak + "Peak"
            if "hkl" in settings_class.fit_orders[lists[z, 1]]["peak"][lists[z, 2]]:
                out_peak = (
                    out_peak
                    + " ("
                    + str(
                        settings_class.fit_orders[lists[z, 1]]["peak"][lists[z, 2]][
                            "hkl"
                        ]
                    )
                    + ")"
                )
            else:
                out_peak = out_peak + " " + str([lists[z, 2]])
            RowLst["DataFile"] = out_name
            RowLst["Peak"] = out_peak

            RowLst["Range_start"] = data_to_write["range"][0][0]
            RowLst["Range_end"] = data_to_write["range"][0][1]

            for w in range(len(coefs_vals_get)):
                ind = coefs_vals_get[w]
                ind_err = ind + "_err"

                if ind != "background" and ind != "symmetry":
                    for v in range(len(data_to_write["peak"][lists[z, 2]][ind])):
                        if (
                            data_to_write["peak"][lists[z, 2]][ind][v] is None
                        ):  # catch  'null' as an error
                            data_to_write["peak"][lists[z, 2]][ind][v] = np.nan
                        if (
                            data_to_write["peak"][lists[z, 2]][ind_err][v] is None
                        ):  # catch  'null' as an error
                            data_to_write["peak"][lists[z, 2]][ind_err][v] = np.nan
                        RowLst[ind + str(v)] = data_to_write["peak"][lists[z, 2]][ind][
                            v
                        ]
                        try:
                            RowLst[ind + str(v) + "err"] = data_to_write["peak"][
                                lists[z, 2]
                            ][ind_err][v]
                        except:
                            RowLst[ind + str(v) + "err"] = "None"

                elif ind == "symmetry":
                    try:
                        if (
                            data_to_write["peak"][lists[z, 2]][ind] is None
                        ):  # catch  'null' as an error
                            data_to_write["peak"][lists[z, 2]][ind] = np.nan
                        RowLst[ind] = data_to_write["peak"][lists[z, 2]][ind]
                    except:
                        pass

                else:  # background
                    for u in range(len(data_to_write[ind])):
                        for v in range(len(data_to_write[ind][u])):
                            if (
                                data_to_write[ind][u][v] is None
                            ):  # catch  'null' as an error
                                data_to_write[ind][u][v] = np.nan
                            if (
                                data_to_write[ind_err][u][v] is None
                            ):  # catch  'null' as an error
                                data_to_write[ind_err][u][v] = np.nan
                            RowLst[ind + str(u) + "_" + str(v)] = data_to_write[ind][u][
                                v
                            ]
                            try:
                                RowLst[ind + str(u) + "_" + str(v) + "err"] = (
                                    data_to_write[ind_err][u][v]
                                )
                            except:
                                RowLst[ind + str(u) + "_" + str(v) + "err"] = "None"
                                
            for w in range(len(properties)):
                ind = properties[w]

                if ind == "status":
                    # status is a list
                    try:
                        if (
                            data_to_write["FitProperties"][ind] is None
                        ):  # catch  'null' as an error
                            data_to_write["FitProperties"][ind] = np.nan
                        RowLst[ind] = re.sub(', ', ';', str(data_to_write["FitProperties"][ind]))
                    except:
                        pass

                else:
                    if (
                        data_to_write["FitProperties"][ind] is None
                    ):  # catch  'null' as an error
                        data_to_write["FitProperties"][ind] = np.nan
                    RowLst[ind] = data_to_write["FitProperties"][ind]

            RowsList.append(RowLst)

    # make data frame using headers - so columns are in good order.
    df = pd.DataFrame(RowsList, columns=headers)

    return df
