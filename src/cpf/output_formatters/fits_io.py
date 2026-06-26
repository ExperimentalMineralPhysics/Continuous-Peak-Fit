#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ["WriteFits", "ReadFits_to_list", "ReadFits_to_dataframe", "read_metadata"]


import json
import os
import re
from itertools import product

# import glob
import numpy as np
import pandas as pd

import cpf.peak_functions as pf
from cpf.output_formatters.convert_fit_to_crystallographic import (
    fourier_to_crystallographic,
)
from cpf.series_functions import series_properties, get_combined_series
from cpf.settings import get_settings
from cpf.util.io import (
    make_outfile_name,
    numpy_to_json,
    peak_hkl,
    peak_phase,
    replace_value,
)
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.fits_io")


def WriteFits(
    settings_class, fitted_param, filename_to_write=None, data_class=None, mode=None
):
    """
    Write fits and any metadata to json files.

    Writees the files as:
    {"fits"":      [list of fits for each range in settings_class.fit_orders],
     "metadata": {dict of metadata properties as named in settings class .... }
     }

    Parameters
    ----------
    settings_class : cpf settings class
        Settings class used for fitting the data.
    fitted_param : dict
        Dictionary of the fits returned by cpf. To be written to the file
    filename_to_write : string, optional
        String setting the file to be written.
        If present later optinal parameters are ignored.
    data_class : cpf data class, optional
        Data class that contains the image metadata. The default is None.
    mode : string, optional
        Switch to add string to file name (if not spedified). The default is None.

    Returns
    -------
    None.
    """
    # try and get the meta data
    # prepare output dictionary
    if filename_to_write and "PreviousFit" in filename_to_write:
        # prevent change in behaviour for now
        # FIXME: cleanup in future push
        if "fits" in fitted_param:
            out = fitted_param["fits"]
        else:
            out = fitted_param
        for i in out:
            i.pop("correlation_coeffs", None)
    elif data_class:
        metadata = data_class.get_metadata(settings_class=settings_class)
        out = {"metadata": metadata, "fits": fitted_param}
    else:
        out = {"fits": fitted_param}

    if filename_to_write is None:
        if mode == "search":
            additional_text = settings_class.file_label
        else:
            additional_text = None
        filename_to_write = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            additional_text=additional_text,
            extension=".json",
            overwrite=True,
        )

    with open(filename_to_write, "w") as TempFile:
        # Write a JSON string into the file.
        json.dump(
            out,
            TempFile,
            sort_keys=True,
            indent=2,
            default=numpy_to_json,
        )

def ReadFits_to_list(settings, replace=True, **kwargs):
    """
    Read coefficents from json fits files and return values as list.

    Parameters
    ----------
    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or
        dictionary of all the settings or
        string or path to a file with the settings in.

    Raises
    ------
    ValueError
        Raised if nether settings_class or settings_file is present.

    Returns
    -------
    fits : list
        List contiaing all the fits made when calling the settings_class/file.
    metadata : list
        List contiaing metadata for each file in settings_class/file.

    """
    
    # get kwargs that might be present 
    add_integrated = kwargs.get("add_integrated", False)
    
    if isinstance(settings, str) and "PreviousFit" in settings:
        # read previous fit
        with open(settings) as json_data:
            json_contents = json.load(json_data)
        if "fits" in json_contents:
            fits = json_contents["fits"]
        else:
            fits = json_contents
        metadata = None
    else:
        # make sure settings is a class
        settings_class = get_settings(settings)

        # read all the data.
        fits = []
        metadata = []
        for z in range(settings_class.image_number):
            settings_class.set_subpattern(z, 0)

            if settings_class.file_label:
                additional_text = settings_class.file_label
            else:
                additional_text = None
            filename = make_outfile_name(
                settings_class.subfit_filename,
                directory=settings_class.output_directory,
                extension=".json",
                additional_text=additional_text,
                overwrite=True,
            )

            # Read JSON data from file
            with open(filename) as json_data:
                json_contents = json.load(json_data)
                if isinstance(json_contents, dict):
                    # new style as dictionary with metadata
                    fits.append(json_contents["fits"])
                    metadata.append(json_contents["metadata"])
                else:
                    # old stype without metadata
                    fits.append(json_contents)
                    metadata.append([])
                if isinstance(settings_class.subfit_filename, list):
                    fnam = settings_class.subfit_filename[0]
                else:
                    fnam = settings_class.subfit_filename
                if (os.path.isfile(fnam) and
                    sorted(settings_class.metadata) != sorted(list(metadata[-1]))
                ):
                    # then we need to read the metadata from the files
                    metadata[-1] = read_metadata(settings_class)

                if add_integrated:
                    # create an integrated series and add to the fits
                    if not settings_class.data_class.continuous_azm:
                        azimuths = settings_class.data_class.azm
                        import numpy.ma as ma
                        if ma.isMaskedArray(azimuths):
                            azimuths = azimuths.compressed()
                    else:
                        azimuths = None
                    for i in range(len(fits[-1])):
                        strt_nd = fits[z][i]["range"][0]
                        settings_class.set_subpattern(z, i)
                        for j in range(len(fits[-1][i]["peak"])):
                            comb_series = get_combined_series(
                                                        fits[-1][i]["peak"][j], 
                                                        azimuth = azimuths,
                                                        start_end=strt_nd,
                                                        **kwargs)
                            fits[-1][i]["peak"][j].update(comb_series)
            
            # convert correlation coefficients into panda data frame
            for y in range(len(fits[-1])):
                if "correlation_coeffs" in fits[-1][y]:
                    try:
                        fits[-1][y]["correlation_coeffs"] = pd.DataFrame.from_dict(
                            json.loads(fits[-1][y]["correlation_coeffs"])
                        )
                    except:
                        pass

    if replace:
        # keep the null terms if we want/need.
        # used for keeting errors in the previous fits
        fits = replace_value(fits, old=None, new=0)
    return fits, metadata


def ReadFits_to_dataframe(
    settings,
    includeParameters="all",
    includeStats=False,
    includeSeriesValues = False,
    includeIntegrated=False,
    includeIntensityRanges = False,
    includeUnitCells = False,
    includePosition = False,
    *args,
    **kwargs,
):
    """
    Read coefficents from json fits files and return values as panda dataframe.

    Parameters
    ----------
    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or
        dictionary of all the settings or
        string or path to a file with the settings in.
    includeParameters : list[str], optional
        List of which peak parameters to return. The default is "all".
    includeStats : bool, optional
        Switch to include all fitting statistics in output data frame. The default is False.
    includeSeriesValues : bool or list, optional
        Switch to include values derived from the fit parameters. Either a list of parameters returned by
        cpf.output_formatters.convert_fit_to_crystallographic or a bool. The default is False.

    Raises
    ------
    ValueError
        Raised if nether settings_class or settings_file is present.

    Returns
    -------
    df : Panda data frame
        Data frame contiaing all the fits made when calling the settings_class/file.

    """
    # make sure settings is a class
    settings_class = get_settings(settings)

    # get what to write
    if includeParameters is False:
        includeParameters = []
    elif "includeParameters" in settings_class.output_settings:
        includeParameters = settings_class.output_settings["includeParameters"]
    if isinstance(includeParameters, str):
        includeParameters = [includeParameters]
    if includeParameters == ["all"]:
        peak_properties = pf.peak_components(full=True, include_combined=includeIntegrated)
        includeParameters = peak_properties[1]

    if includeSeriesValues is not False or includeUnitCells is not False:
        # Parse needed parameters
        SampleGeometry = settings_class.output_settings.get("SampleGeometry", "3d")
        SampleDeformation = settings_class.output_settings.get(
            "SampleDeformation", "compression"
        )
        # override with kwargs
        SampleGeometry = kwargs.get("SampleGeometry", SampleGeometry)
        SampleDeformation = kwargs.get("SampleDeformation", SampleDeformation)
        # force all the kwargs that might be needed
        set_params = {
            "SampleGeometry": SampleGeometry,
            "SampleDeformation": SampleDeformation,
        }
        kwargs.update(set_params)
    kwargs.update({"add_integrated": includeIntegrated})
    
    if includeIntensityRanges is not False: 

        # get the intensity maximum and minimum of the fit, model and residuals
        IntensityValues = [
            "data_max",
            "data_min",
            "model_max",
            "model_min",
            "residual_max",
            "residual_min",
        ]
    else:
        IntensityValues = []

    # read all the data.
    fits, metadata = ReadFits_to_list(settings_class, **kwargs)

    num_fits = 0
    max_peaks = 0
    for z in range(settings_class.image_number):
        settings_class.set_subpattern(z, 0)

        if includeSeriesValues is not False:
            # get converted values.
            for i in range(len(fits[z])):
                for j in range(len(fits[z][i]["peak"])):
                    crystallographic_values = fourier_to_crystallographic(
                        fits[z],
                        SampleGeometry=SampleGeometry,
                        SampleDeformation=SampleDeformation,
                        subpattern=i,
                        peak=j,
                    )
                    fits[z][i]["peak"][j]["crystallographic_values"] = (
                        crystallographic_values
                    )

                    # it not continuous_azm then need to know where the detectors are
                    if settings_class.data_class.continuous_azm == False:
                        data_class = settings_class.data_class
                        data_class.fill_data(
                            settings_class.image_list[0],
                            settings=settings_class,
                        )
                        azms = np.unique(data_class.azm)
                    else:
                        azms = 0.01  # default spacing

                    height_properties = series_properties(
                        fits[z], subpattern=i, peak=j, param="height", azm_spacing=azms
                    )
                    width_properties = series_properties(
                        fits[z], subpattern=i, peak=j, param="width", azm_spacing=azms
                    )
                    profile_properties = series_properties(
                        fits[z], subpattern=i, peak=j, param="profile", azm_spacing=azms
                    )

                    fits[z][i]["peak"][j]["crystallographic_values"] = (
                        fits[z][i]["peak"][j]["crystallographic_values"]
                        | height_properties
                    )
                    fits[z][i]["peak"][j]["crystallographic_values"] = (
                        fits[z][i]["peak"][j]["crystallographic_values"]
                        | width_properties
                    )
                    fits[z][i]["peak"][j]["crystallographic_values"] = (
                        fits[z][i]["peak"][j]["crystallographic_values"]
                        | profile_properties
                    )

                    if includeIntegrated:
                        extras = (set(
                            pf.peak_components(full=True, include_profile=True, include_combined=True)[1]) - 
                            set(pf.peak_components(full=True, include_profile=True, include_combined=False)[1])
                        )
                        for k in extras:
                            extra_properties = series_properties(
                                fits[z], subpattern = i, peak=j, param=k, azm_spacing=azms
                            )
                            fits[z][i]["peak"][j]["crystallographic_values"] = (
                                fits[z][i]["peak"][j]["crystallographic_values"] 
                                | extra_properties
                            )

        if includeSeriesValues is not False:
            # list the entries in crystallographic_values dictionary
            DerivedValues = fits[z][0]["peak"][0]["crystallographic_values"].keys()
        else:
            DerivedValues = []
        num_fits = np.max([num_fits, len(fits[z])])
        for y in range(len(fits[z])):
            max_peaks = np.max([max_peaks, len(fits[z][y]["peak"])])

    # get number of coefficients.
    # get the coefficients from the json file rather than the setting/input file
    # because the number of fits might not the same as in the settings file,
    # for example, afer running cpf.XRD_FitPatter.order_search()
    max_coef = {}
    for w in range(len(includeParameters)):
        ind = includeParameters[w]
        if ind != "background":
            max_coef[includeParameters[w]] = 0
        else:
            max_coef[includeParameters[w]] = [0]

    for w in range(len(includeParameters)):
        ind = includeParameters[w]
        for x in range(num_fits):
            # loop over num fits assumes all the json files are the same size.
            if ind == "background":
                # iterate over the length of background
                for y in range(len(fits[0][x][ind])):
                    if y >= len(max_coef[ind]):
                        max_coef[ind].append(len(fits[0][x][ind][y]))
                    else:
                        max_coef[ind][y] = np.max(
                            [max_coef[ind][y], len(fits[0][x][ind][y])]
                        )
            else:  # peak related parameter
                # iterate over the number of peaks
                for y in range(len(fits[0][x]["peak"])):
                    if ind not in fits[0][x]["peak"][y]:
                        # cannot assume the symmetry is present
                        max_coef[ind] = np.max([max_coef[ind], 0])
                    elif isinstance(fits[0][x]["peak"][y][ind], int):
                        max_coef[ind] = np.max([max_coef[ind], 1])
                    else:
                        max_coef[ind] = np.max(
                            [max_coef[ind], len(fits[0][x]["peak"][y][ind])]
                        )

    # make list of headers for panda data frame
    headers = []
    headers.append("num")
    headers.append("DataFile")
    headers.append("phase")
    headers.append("peak")
    # add metadata to list
    if settings_class.metadata:
        for i in settings_class.metadata:
            headers.append(i)
    headers.append("range_start")
    headers.append("range_end")

    # parmeter header list
    for w in range(len(max_coef)):
        ind = includeParameters[w]
        if ind == "symmetry":
            headers.append(ind)
        elif ind == "background":
            headers.append(ind + "_type")
            for u in range(len(max_coef[ind])):
                for v in range(max_coef[ind][u]):
                    headers.append(ind + str(u) + "_" + str(v))
                    headers.append(ind + str(u) + "_" + str(v) + "_err")
        else:  # ind is a peak parameter
            headers.append(ind + "_type")
            for v in range(max_coef[ind]):
                headers.append(ind + str(v))
                headers.append(ind + str(v) + "_err")
    properties = []
    if includePosition == True:
        properties.append("Position in json")
    if includeSeriesValues is not False:
        if includeSeriesValues is True:
            properties += DerivedValues
        else:
            for i in range(len(includeSeriesValues)):
                properties.append(includeSeriesValues[i])
                # properties.append(includeSeriesValues[i]+"err")
    if includeIntensityRanges is True:
        properties += IntensityValues
    if includeStats is True:
        # include properties from the lmfit output that were passed with the fits.
        # read the list of parameters from the first file.
        if "FitProperties" in fits[0][0]:
            properties += list(fits[0][0]["FitProperties"])
    if "note" in fits[0][0]:
        # additional test added by cpf.XRD_FitPattern.order_search
        properties.append("note")
        # more values neeed by cpf.Output_Formatters.WriteOrderSearchFigure.
        # which always adds notes to the json fit files.
        properties.append("pos_in_range")
    headers += properties

    # make lists of the parameters to iterate over
    # then sort them in to order by peak
    images = list(range(settings_class.image_number))
    subpatterns = list(range(num_fits))  # list(range(num_subpatterns))
    max_peaks = 1
    for i in range(len(settings_class.fit_orders)):
        max_peaks = np.max([max_peaks, len(settings_class.fit_orders[i]["peak"])])
    lists = images, subpatterns, list(range(max_peaks))
    lists = np.array(list(product(*lists)))
    # sort lists by peaks
    lists = lists[np.lexsort((lists[:, 2],))]
    lists = lists[np.lexsort((lists[:, 1],))]

    # write the data to array, then append to dataframe
    RowsList = []
    for z in range(len(lists)):
        settings_class.set_subpattern(lists[z, 0], 0)
        RowLst = {}
        data_to_write = fits[lists[z, 0]][lists[z, 1]]

        if len(data_to_write["peak"]) > lists[z, 2]:
            RowLst["num"] = lists[z, 0]
            if isinstance(settings_class.subfit_filename, list):
                # filenames have to be unique but will be a list of h5 type files
                RowLst["DataFile"] = os.path.split(settings_class.subfit_filename[0])[1]
            else:
                RowLst["DataFile"] = os.path.split(settings_class.subfit_filename)[1]

            # RowLst["Peak"] = peak_string(fits[lists[z, 0]][lists[z, 1]], peak=[lists[z, 2]], fname=False)
            RowLst["phase"] = peak_phase(
                fits[lists[z, 0]][lists[z, 1]], peak=[lists[z, 2]]
            )[0]
            RowLst["peak"] = peak_hkl(
                fits[lists[z, 0]][lists[z, 1]], peak=[lists[z, 2]]
            )[0]
            RowLst["range_start"] = data_to_write["range"][0][0]
            RowLst["range_end"] = data_to_write["range"][0][1]

            for w in settings_class.metadata:
                if "/" in w:
                    # cut to last part of h5key
                    if w[-1] == "/":
                        last = -2
                    else:
                        last = -1
                    RowLst[w] = metadata[lists[z, 0]][w.split("/")[last]]
                else:
                    RowLst[w] = metadata[lists[z, 0]][w]

            for w in range(len(includeParameters)):
                ind = includeParameters[w]
                ind_err = ind + "_err"

                if ind != "background" and ind != "symmetry":
                    RowLst[ind + "_type"] = data_to_write["peak"][lists[z, 2]][
                        ind + "_type"
                    ]
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
                            RowLst[ind + str(v) + "_err"] = data_to_write["peak"][
                                lists[z, 2]
                            ][ind_err][v]
                        except:
                            RowLst[ind + str(v) + "_err"] = "None"

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
                    RowLst[ind + "_type"] = data_to_write[ind + "_type"]
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
                                RowLst[ind + str(u) + "_" + str(v) + "_err"] = (
                                    data_to_write[ind_err][u][v]
                                )
                            except:
                                RowLst[ind + str(u) + "_" + str(v) + "_err"] = "None"

            for w in range(len(properties)):
                ind = properties[w]

                if ind == "note":
                    # note is a string
                    if "note" in fits[lists[z, 0]][lists[z, 1]]:
                        RowLst[ind] = fits[lists[z, 0]][lists[z, 1]]["note"]

                elif ind == "pos_in_range":
                    RowLst[ind] = lists[z, 2]

                elif ind in DerivedValues:
                    # in crystallographic_values dictionary
                    RowLst[ind] = data_to_write["peak"][lists[z, 2]][
                        "crystallographic_values"
                    ][ind]
                    # RowLst[ind+"err"] = data_to_write["peak"][lists[z, 2]]["crystallographic_values"][ind+"_err"]

                elif ind == "Position in json":
                    # print("here we are ")
                    # print(lists[z, 1])
                    RowLst[ind] = lists[z, 1]

                elif ind in IntensityValues:
                    if ind == "data_max":
                        RowLst[ind] = data_to_write["DataProperties"]["max"]
                    elif ind == "data_min":
                        RowLst[ind] = data_to_write["DataProperties"]["min"]
                    elif ind == "model_max":
                        RowLst[ind] = data_to_write["ModelProperties"]["max"]
                    elif ind == "model_min":
                        RowLst[ind] = data_to_write["ModelProperties"]["min"]
                    elif ind == "residual_max":
                        RowLst[ind] = data_to_write["ResidualProperties"]["max"]
                    elif ind == "residual_min":
                        RowLst[ind] = data_to_write["ResidualProperties"]["min"]
                    else:
                        raise ValueError("Unknown value to read")

                else:
                    if data_to_write["FitProperties"][ind] is None:
                        # catch 'null' as an error
                        data_to_write["FitProperties"][ind] = np.nan
                    if isinstance(data_to_write["FitProperties"][ind], str):
                        # make sure that status, or another string, does not have commass.
                        RowLst[ind] = re.sub(
                            ", ", ";", str(data_to_write["FitProperties"][ind])
                        )
                    else:
                        RowLst[ind] = data_to_write["FitProperties"][ind]

            RowsList.append(RowLst)

    # make data frame using headers - so columns are in sensible order.
    df = pd.DataFrame(RowsList, columns=headers)

    return df


def read_metadata(settings_class):
    """
    Reads the metadata for the image specified in settings_class.subfit_filename.

    Parameters
    ----------
    settings_class : cpf settings class
        Settings clsss in which settings_class.subfit_filename is set to the file
        to be read.

    Returns
    -------
    metadata : dict
        metadata of the diffraction image
    """
    if isinstance(settings_class.subfit_filename, list):
        fnam = settings_class.subfit_filename[0]
    else:
        fnam = settings_class.subfit_filename

    if settings_class.subfit_filename == None:
        raise ValueError("no file is specified")
    if not os.path.isfile(fnam):
        raise FileExistsError(
            "The file {os.path.split(fnam)[1]} does not exist on the path"
        )

    # get data from settings class
    new_data = settings_class.data_class

    new_data.fill_data(
        settings_class.subfit_filename,
        settings=settings_class,
    )

    new_data.import_image(settings=settings_class)
    metadata = new_data.get_metadata(metadata_values=settings_class.metadata)
    return metadata
