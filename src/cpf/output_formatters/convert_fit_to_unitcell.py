__all__ = ["fits_to_unitcell"]

import glob
import json
import re

import numpy as np
import pandas as pd

from cpf.output_formatters.convert_fit_to_crystallographic import (
    fourier_to_crystallographic,
    fourier_to_unitcellvolume,
)
from cpf.output_formatters.fits_io import ReadFits_to_dataframe, ReadFits_to_list

# from uncertainties import ufloat
from cpf.output_formatters.jcpds import jcpds
from cpf.settings import get_settings
from cpf.util.io import make_outfile_name, peak_hkl, replace_null_terms
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.convert_fit_to_unitcell")


def fits_to_unitcell(settings, *args, **kwargs):
    """
    Processes all fits and returns dataframe of unit parameters calculated from
    values in fit (json) files.

    Parameters
    ----------
    settings : cpf.Settings.settings() Class, str (filename), Path
        Class containing all the fitting parameters.
    *args

    **kwargs

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

    # settings_class.output_settings["phase"]

    # Parse required parameters
    SampleGeometry = settings_class.output_settings.get("SampleGeometry", "3d")
    SampleDeformation = settings_class.output_settings.get(
        "SampleDeformation", "compression"
    )
    phase = settings_class.output_settings.get("phase", True)
    jcpds_file = settings_class.output_settings.get("jcpds", None)
    # override with kwargs
    SampleGeometry = kwargs.get("SampleGeometry", SampleGeometry)
    SampleDeformation = kwargs.get("SampleDeformation", SampleDeformation)
    phase = kwargs.get("phase", phase)
    jcpds_file = kwargs.get("jcpds", jcpds_file)

    # set the kwargs as needed.
    kwargs["includeSeriesValues"] = kwargs.get("includeSeriesValues", True)
    kwargs["includePosition"] = kwargs.get("includePosition", True)

    pressure = kwargs.get("pressure", False)
    kwargs["pressure"] = pressure

    # kwargs["phase"] = phase
    # kwargs["jcpds"] = jcpds_file

    # get the fits
    all_fits, all_metadata = ReadFits_to_list(settings=settings_class)

    # read all the data.
    all_cells = []

    for z in range(settings_class.image_number):
        settings_class.set_subpattern(z, 0)

        # if settings_class.file_label:
        #     additional_text = settings_class.file_label
        # else:
        #     additional_text = None

        # filename = make_outfile_name(
        #     settings_class.subfit_filename,  # diff_files[z],
        #     directory=settings_class.output_directory,  # directory=FitSettings.Output_directory,
        #     extension=".json",
        #     additional_text=additional_text,
        #     overwrite=True,
        # )  # overwrite =false to get the file name without incrlemeting it.

        # Read JSON data from file
        # with open(filename) as json_data:
        #     fits = json.load(json_data)

        fits = all_fits[z]
        metadata = all_metadata[z]

        # get converted values.
        for i in range(len(fits)):
            for j in range(len(fits[i]["peak"])):
                crystallographic_values = fourier_to_crystallographic(
                    fits,
                    SampleGeometry=SampleGeometry,
                    SampleDeformation=SampleDeformation,
                    subpattern=i,
                    peak=j,
                )
                fits[i]["peak"][j]["crystallographic_values"] = crystallographic_values

        # stash names for output
        cells_tmp = {}
        cells_tmp["num"] = z
        cells_tmp["DataFile"] = make_outfile_name(
            settings_class.subfit_filename,
            directory="",
            extension="",
            overwrite=True,
        )
        # add metadata
        for i in settings_class.metadata:
            cells_tmp[i] = metadata[i]

        # get or guess phase
        if not isinstance(phase, str):
            # list all phases in fits
            phases = []
            for i in range(len(fits)):
                for j in range(len(fits[i]["peak"])):
                    if "phase" in fits[i]["peak"][j]:
                        phases.append(fits[i]["peak"][j]["phase"])
            phase = np.unique(phases)

        if isinstance(phase, str):
            phase = [phase]

            # if "phase" in settings_class.output_settings:
            #     phase = settings_class.output_settings["phase"]

            # else:
            #     #list all phases in fits
            #     phases = []
            #     for i in range(len(fits)):
            #         for j in range(len(fits[i]["peak"])):
            #             if "phase" in fits[i]["peak"][j]:
            #                 phases.append(fits[i]["peak"][j]["phase"])
            #     phase = np.unique(phases)

        # get or guess jcpds file

        if not isinstance(phase, str):
            jcpds_file = []
            for i in range(len(phase)):
                if glob.glob(f"*{phase[i]}*.jcpds"):
                    if len(glob.glob(f"*{phase[i]}*.jcpds")) != 1:
                        raise ValueError("There is more than 1 jcpds file")
                    jcpds_file.append(glob.glob(f"*{phase[i]}*.jcpds")[0])
                elif glob.glob(f"*{phase[i]}*.cif"):
                    if len(glob.glob(f"*{phase[i]}*.cif")) != 1:
                        raise ValueError("There is more than 1 cif file")
                    jcpds_file.append(glob.glob(f"*{phase[i]}*.cif")[0])
            if len(jcpds_file) == 0:
                raise ValueError("There is no jcpds or cif file recognised")
            elif len(phase) != len(jcpds_file):
                raise ValueError("The phase and jcpds files do not match")

        # if "jcpds" in settings_class.output_settings:
        #     jcpds = settings_class.output_settings["jcpds"]
        # else:
        #     jcpds = []
        #     for i in range(len(phase)):
        #         if glob.glob(f"*{phase[i]}*.jcpds"):
        #             if len(glob.glob(f"*{phase[i]}*.jcpds")) != 1:
        #                 raise ValueError("There is more than 1 jcpds file")
        #             jcpds.append(glob.glob(f"*{phase[i]}*.jcpds")[0])
        #         elif glob.glob(f"*{phase[i]}*.cif"):
        #             if len(glob.glob(f"*{phase[i]}*.cif")) != 1:
        #                 raise ValueError("There is more than 1 cif file")
        #             jcpds.append(glob.glob(f"*{phase[i]}*.cif")[0])
        #     if len(jcpds) == 0:
        #         raise ValueError("There is no jcpds or cif file recognised")
        #     elif len(phase) != len(jcpds):
        #         raise ValueError("The phase and jcpds files do not match")

        # get temperature from metadata
        # get label for temperature. Should work for wild cards
        if "temperature" in settings_class.metadata_labels:
            templbl_without_wildcards = re.sub(
                r"\*", ".*", settings_class.metadata_labels["temperature"]
            )
        elif "temperature" in settings_class.data_class._default_metadata_labels:
            templbl_without_wildcards = re.sub(
                r"\*",
                ".*",
                settings_class.data_class._default_metadata_labels["temperature"],
            )
        else:
            templbl_without_wildcards = "None"
        r = re.compile(templbl_without_wildcards)
        templbl = list(filter(r.match, list(metadata)))  # Read Note below
        if len(templbl) == 1:
            templbl = templbl[0]
        elif len(templbl) > 1:
            err_str = (
                "More than one temprature has been found. Assuming the first one. "
            )
            logger.error(err_str)
            templbl = templbl[0]
        elif templbl == []:
            # empty list cause by not finding temperature
            templbl = ""
        # get temperature
        if templbl in metadata:
            temp = metadata[templbl]
        else:
            temp = 0

        # calculate unit cell properties and return them
        for i in range(len(phase)):
            kwargs_here = kwargs

            kwargs_here.pop("phase", None)
            kwargs_here.pop("jcpds", None)
            unitcells = fourier_to_unitcellvolume(
                fits,
                # SampleGeometry=SampleGeometry,
                # SampleDeformation=SampleDeformation,
                phase=phase[i],
                jcpds_file=jcpds_file[i],
                temperature=temp,
                **kwargs_here,
            )

            # label return with phase name and add to fits
            entries = list(unitcells)
            for j in range(len(unitcells)):
                unitcells[
                    re.sub(entries[j], phase[i] + " " + entries[j], entries[j])
                ] = unitcells.pop(entries[j])

            cells_tmp.update(unitcells)

        all_cells.append(cells_tmp)

    # make data frame using headers - so columns are in sensible order.
    df = pd.DataFrame(all_cells)

    return df
