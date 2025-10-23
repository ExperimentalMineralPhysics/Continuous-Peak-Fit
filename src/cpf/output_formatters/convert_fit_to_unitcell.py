__all__ = ["fits_to_unitcells"]

import json
import pandas as pd
import numpy as np
import glob
import re
# from uncertainties import ufloat

from cpf.output_formatters.jcpds import jcpds
from cpf.IO_functions import peak_hkl
from cpf.output_formatters.ReadFits import ReadFits
from cpf.IO_functions import replace_null_terms
from cpf.IO_functions import make_outfile_name
from cpf.output_formatters.convert_fit_to_crystallographic import fourier_to_crystallographic, fourier_to_unitcellvolume
from cpf.util.logging import get_logger




def fits_to_unitcell(
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
    Processes all fits and returns dataframe of unit parameters calculated from 
    values in fit (json) files. 

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
        
    # force all the kwargs that might be needed
    kwargs.pop("SampleGeometry", "3d")
    kwargs.pop("SampleDeformation", "compression")   
    
    SampleGeometry = kwargs.get("SampleGeometry", "3d")
    SampleDeformation = kwargs.get("SampleDeformation", "compression")
        
    df = ReadFits(
            settings_class=settings_class,
            settings_file=settings_file,
            includeParameters = False,
            includeStats=False,
            includeSeriesValues = True,
            includeIntensityRanges = False,
            includeUnitCells = False,
            includePosition = includePosition,
            *args,
            **kwargs)
    
    # read all the data.
    all_cells = []
    
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
            fits = json.load(json_data)
            
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
                        
            # get or guess phase
            if "phase" in settings_class.output_settings:
                phase = settings_class.output_settings["phase"]
            else:
                #list all phases in fits
                phases = []
                for i in range(len(fits)):
                    for j in range(len(fits[i]["peak"])):
                        if "phase" in fits[i]["peak"][j]:
                            phases.append(fits[i]["peak"][j]["phase"])
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
                    raise ValueError("The phase and jcpds files do not match")
                    
            # calculate unit cell properties and return them
            for i in range (len(phase)):
                unitcells = fourier_to_unitcellvolume(
                    fits,
                    # SampleGeometry=SampleGeometry,
                    # SampleDeformation=SampleDeformation,
                    phase = phase[i],
                    jcpds_file = jcpds[i],
                    Pressure=False,
                    **kwargs
                )
                
                # label return with phase name and add to fits                    
                entries = list(unitcells)
                for j in range(len(unitcells)):
                    unitcells[re.sub(entries[j], phase[i]+" "+entries[j], entries[j])] = unitcells.pop(entries[j])
                    
                cells_tmp.update(unitcells)
                    
            all_cells.append(cells_tmp)
            
    # make data frame using headers - so columns are in sensible order.
    df = pd.DataFrame(all_cells)
    
    return df
              