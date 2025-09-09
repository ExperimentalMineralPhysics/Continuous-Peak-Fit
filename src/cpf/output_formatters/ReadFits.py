#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ["ReadFits"]


import json
from itertools import product
import re
import glob

import numpy as np
import pandas as pd

import cpf.peak_functions as pf
import cpf.output_formatters.convert_fit_to_crystallographic as FitToCryst
from cpf.IO_functions import make_outfile_name, peak_string
from cpf.series_functions import series_properties
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.ReadFits")


def ReadFits(
    settings_class=None,
    settings_file=None,
    includeParameters = "all",
    includeStats=False,
    includeSeriesValues = False,
    includeIntensityRanges = False,
    includeUnitCells = False,
    includePosition = False,
    *args,
    **kwargs
):
    """
    Read coefficents from json fits files and return values as panda dataframe. 

    Parameters
    ----------
    settings_class : cpf.Settings.settings() Class, optional
        Class containing all the fitting parameters. The default is None.
    settings_file : *.py file, optional
        text file containing all the fitting parameters. The default is None.
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

    if settings_class is None and settings_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        from cpf.XRD_FitPattern import initiate
        settings_class = initiate(settings_file, require_datafiles=False)

    # get what to write
    if "includeParameters" in settings_class.output_settings:
        includeParameters = settings_class.output_settings["includeParameters"]
    if isinstance(includeParameters, str):
        includeParameters = [includeParameters]
    if includeParameters == ["all"]:
        peak_properties = pf.peak_components(full=True)
        includeParameters = peak_properties[1]

    if includeSeriesValues is not False:
        SampleGeometry = kwargs.get("SampleGeometry", "3d")
        SampleDeformation = kwargs.get("SampleDeformation", "compression")
        # SampleGeometry = "3d"
        # SampleDeformation = "compression"
        # if "SampleGeometry" in kwargs:
        #     SampleGeometry = kwargs["SampleGeometry"]
        # if "SampleDeformation" in kwargs:
        #     SampleDeformation = kwargs["SampleDeformation"]
    
    if includeIntensityRanges is not False: 
        # get the intensity maximum and minimum of the fit, model and residuals
        IntensityValues = ["data_max", "data_min", "model_max", "model_min", "residual_max", "residual_min"]
    else:
        IntensityValues = []
        
    # read all the data.
    fits = []
    num_fits = 0
    max_peaks = 0
    for z in range(settings_class.image_number):
        settings_class.set_subpattern(z, 0)

        if settings_class.file_label:
            additional_text = settings_class.file_label
        else:
            additional_text = None

        filename = make_outfile_name(
            settings_class.subfit_filename,  # diff_files[z],
            directory=settings_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".json",
            additional_text=additional_text,
            overwrite=True,
        )  # overwrite =false to get the file name without incrlemeting it.

        # Read JSON data from file
        with open(filename) as json_data:
            fits.append(json.load(json_data))
            if includeSeriesValues is not False:
                # get converted values.
                for i in range(len(fits[-1])):
                    for j in range(len(fits[-1][i]["peak"])):
                        crystallographic_values = FitToCryst.fourier_to_crystallographic(
                            fits[-1],
                            SampleGeometry=SampleGeometry,
                            SampleDeformation=SampleDeformation,
                            subpattern=i,
                            peak=j,
                        )
                        fits[-1][i]["peak"][j]["crystallographic_values"] = crystallographic_values
                        
                        height_properties = series_properties(fits[-1], subpattern = i, peak=j, param="height")
                        width_properties = series_properties(fits[-1], subpattern = i, peak=j, param="width")
                        profile_properties = series_properties(fits[-1], subpattern = i, peak=j, param="profile")
                        
                        fits[-1][i]["peak"][j]["crystallographic_values"] = fits[-1][i]["peak"][j]["crystallographic_values"] | height_properties
                        fits[-1][i]["peak"][j]["crystallographic_values"] = fits[-1][i]["peak"][j]["crystallographic_values"] | width_properties
                        fits[-1][i]["peak"][j]["crystallographic_values"] = fits[-1][i]["peak"][j]["crystallographic_values"] | profile_properties
                        
            if includeUnitCells is not False:
                # calculate unit cell properties and return them
                
                # get or guess phase
                if "phase" in settings_class.output_settings:
                    phase = settings_class.output_settings["phase"]
                else:
                    #list all phases in fits
                    phases = []
                    for i in range(len(fits[-1])):
                        for j in range(len(fits[-1][i]["peak"])):
                            if "phase" in fits[-1][i]["peak"][j]:
                                phases.append(fits[-1][i]["peak"][j]["phase"])
                    phase = np.unique(phases)
                    
                # get or guess jcpds file
                if "jcpds" in settings_class.output_settings:
                    jcpds = settings_class.output_settings["jcpds"]
                else:
                    jcpds = []
                    for i in range(len(phase)):
                        if glob.glob(f"*{phase[i]}*"):
                            if len(glob.glob(f"*{phase[i]}*")) != 1:
                                raise ValueError("There is more than 1 jcpds file")
                            jcpds.append(glob.glob(f"*{phase[i]}*")[0])
                    if len(phase) != len(jcpds):
                        raise ValueError("The phase and jcpds files are not certain")
                        
                for i in range (len(phase)):
                    unitcells = FitToCryst.fourier_to_unitcellvolume(
                        fits[-1],
                        SampleGeometry=SampleGeometry,
                        SampleDeformation=SampleDeformation,
                        phase = phase[i],
                        jcpds_file = jcpds[i]
                    )
                    
                    # label return with phase name and add to fits                    
                    entries = list(unitcells)
                    for j in range(len(unitcells)):
                        unitcells[re.sub(entries[j], phase[i]+" "+entries[j], entries[j])] = unitcells.pop(entries[j])
                    if not "unitcells" in fits[-1][0]:
                        fits[-1][0]["unitcells"] = unitcells
                    else:
                        fits[-1][0]["unitcells"] = fits[-1][0]["unitcells"] | unitcells
            
        if includeSeriesValues is not False:
            #list the entries in crystallographic_values dictionary
            DerivedValues = fits[-1][0]["peak"][0]["crystallographic_values"].keys()
        
        if includeUnitCells is not False:
            #list the entries in crystallographic_values dictionary
            UnitCells = fits[0][0]["unitcells"].keys()

        num_fits = np.max([num_fits, len(fits[-1])])
        for y in range(len(fits[-1])):
            max_peaks = np.max([max_peaks, len(fits[-1][y]["peak"])])
    
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
                        max_coef[ind][y] = np.max([max_coef[ind][y], len(fits[0][x][ind][y])])
            else: #peak related parameter
                # iterate over the number of peaks
                for y in range(len(fits[0][x]["peak"])):
                    if ind not in fits[0][x]["peak"][y]:
                        # cannot assume the symmetry is present
                        max_coef[ind] = np.max([max_coef[ind], 0])
                    elif isinstance(fits[0][x]["peak"][y][ind], int):
                        max_coef[ind] = np.max([max_coef[ind], 1])
                    else:
                        max_coef[ind] = np.max([max_coef[ind], len(fits[0][x]["peak"][y][ind])])
    
    # make list of headers for panda data frame
    headers = []
    headers.append("num")
    headers.append("DataFile")
    headers.append("Peak")
    headers.append("Range_start")
    headers.append("Range_end")
    
    # parmeter header list
    for w in range(len(max_coef)):
        ind = includeParameters[w]
        if ind == "symmetry":
            headers.append(ind)
        elif ind == "background":
            for u in range(len(max_coef[ind])):
                for v in range(max_coef[ind][u]):
                    headers.append(ind + str(u) + "_" + str(v))
                    headers.append(ind + str(u) + "_" + str(v) + "_err")
        else: # ind is a peak parameter
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
    if includeUnitCells is True:
        properties += UnitCells
    if includeStats is True:
        # include properties from the lmfit output that were passed with the fits.
        #read the list of parameters from the first file.
        settings_class.set_subpattern(0, 0)
        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            extension=".json",
            additional_text=settings_class.file_label,
            overwrite=True,
        )  # overwrite = True to get the file name without incrementing it.
        with open(filename) as json_data:
            tmp = json.load(json_data)
            properties += list(tmp[0]["FitProperties"].keys())
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
    subpatterns = list(range(num_fits))#list(range(num_subpatterns))
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
        settings_class.set_subpattern(lists[z,0], 0)
        RowLst = {}
        data_to_write = fits[lists[z, 0]][lists[z, 1]]
        
        if len(data_to_write["peak"]) > lists[z, 2]:
            RowLst["num"] = lists[z, 0]
            
            RowLst["DataFile"] = make_outfile_name(
                settings_class.subfit_filename,
                directory="",
                extension=".json",
                overwrite=True,
            )
            RowLst["Peak"] = peak_string(fits[lists[z, 0]][lists[z, 1]], peak=[lists[z, 2]], fname=False)
            RowLst["Range_start"] = data_to_write["range"][0][0]
            RowLst["Range_end"] = data_to_write["range"][0][1]

            for w in range(len(includeParameters)):
                ind = includeParameters[w]
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
                    RowLst[ind] = data_to_write["peak"][lists[z, 2]]["crystallographic_values"][ind]
                    # RowLst[ind+"err"] = data_to_write["peak"][lists[z, 2]]["crystallographic_values"][ind+"_err"]

                elif ind in UnitCells:
                    # in unitcells dictionary
                    if "unitcells" in data_to_write:
                        RowLst[ind] = data_to_write["unitcells"][ind]
                
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
                    if (data_to_write["FitProperties"][ind] is None):  
                        # catch 'null' as an error
                        data_to_write["FitProperties"][ind] = np.nan
                    if isinstance(data_to_write["FitProperties"][ind], str):
                        # make sure that status, or another string, does not have commass.
                        RowLst[ind] = re.sub(', ', ';', str(data_to_write["FitProperties"][ind]))
                    else:
                        RowLst[ind] = data_to_write["FitProperties"][ind]

            RowsList.append(RowLst)

    # make data frame using headers - so columns are in sensible order.
    df = pd.DataFrame(RowsList, columns=headers)

    return df
