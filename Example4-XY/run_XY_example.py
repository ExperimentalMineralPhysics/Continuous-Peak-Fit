
from copy import deepcopy
import json
import numpy as np
import importlib.util
from pathlib import Path

import cpf
from cpf.settings import Settings

azimuth_bins = 90
tth_bins = 1500

def settings_dict_from_file(file):
    if isinstance(dioptas_settings, str):
        try:
            settingsPath = Path(dioptas_settings)
        except Exception as error:
            raise error
    if not settingsPath.suffix == ".py":
        settingsPath = settingsPath.with_suffix(".py")
    
    run_name = settingsPath.stem
    
    # store all the settings from file in a module class.
    module_name = run_name
    # store all the settings from file in a module class.
    module_name = dioptas_settings
    spec = importlib.util.spec_from_file_location(
        module_name, settingsPath
    )
    settings_as_dict = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(settings_as_dict)
    # convert to a dictionary
    return settings_as_dict.__dict__
    


# %%  make XY image files from example1

dioptas_settings = "../Example1-Fe/BCC1_Dioptas_input"
settings_as_dict = settings_dict_from_file(dioptas_settings)

settings_as_dict['run_name'] = "Dioptas-XY conversion"
settings_as_dict['datafile_directory'] = '../Example1-Fe/' + settings_as_dict['datafile_directory']
settings_as_dict['Output_directory'] = './XY_images'

settings_as_dict['Calib_data'] = '../Example1-Fe/' + settings_as_dict['Calib_data']
settings_as_dict['Calib_mask'] = '../Example1-Fe/' + settings_as_dict['Calib_mask']
settings_as_dict['Calib_param'] = '../Example1-Fe/' + settings_as_dict['Calib_param']

settings_as_dict["reduce_by"] = 1

settings = cpf.XRD_FitPattern.initiate(
    settings_as_dict)
cpf.XRD_FitPattern.write_output(
    settings, 
    out_type = "CalibratedImages", azimuth_bins=azimuth_bins, tth_bins=tth_bins)


# %% run example as XY
settings_XY = Settings()

with open(settings_as_dict['Output_directory']+'/BCC1_2GPa_10s_001_00010__calibration.json') as f:
    calib = json.load(f)

settings_as_dict_XY = {}
settings_as_dict_XY['datafile_directory'] = settings_as_dict['Output_directory']
settings_as_dict_XY['datafile_Basename'] = settings_as_dict['datafile_Basename']
settings_as_dict_XY['datafile_Ending'] = "__rebinned"+settings_as_dict['datafile_Ending']
settings_as_dict_XY['datafile_StartNum'] = settings_as_dict['datafile_StartNum']
settings_as_dict_XY['datafile_EndNum'] = settings_as_dict['datafile_EndNum']
settings_as_dict_XY['datafile_NumDigit'] = settings_as_dict['datafile_NumDigit']
settings_as_dict_XY['datafile_Step'] = settings_as_dict['datafile_Step']
settings_as_dict_XY['Calib_type'] = "XY"
settings_as_dict_XY['Calib_param'] = calib["calibration"]
settings_as_dict_XY['Calib_mask'] = {'threshold': [0,np.inf]}
settings_as_dict_XY['Output_directory'] = './results'
settings_as_dict_XY['Output_type'] = settings_as_dict['Output_type']
settings_as_dict_XY['fit_orders'] = settings_as_dict['fit_orders']
settings_as_dict_XY['run_name'] = settings_as_dict['run_name']

settings_XY = Settings()
settings_XY.populate(settings_as_dict_XY)

cpf.XRD_FitPattern.initiate(settings_XY)

cpf.XRD_FitPattern.write_output(settings_XY, out_type="CollectionMovie")

# copy settings and then set range. Do this because set_range changes settings and discards images. 
settings_XY2 = deepcopy(settings_XY)
cpf.XRD_FitPattern.set_range(settings_XY2)

cpf.XRD_FitPattern.execute(settings_XY, parallel=False)

cpf.XRD_FitPattern.write_output(settings_XY, out_type=["FitMovie"])

settings_XY.save_settings("test_save")
settings_new = cpf.XRD_FitPattern.initiate("./results/test_save")

