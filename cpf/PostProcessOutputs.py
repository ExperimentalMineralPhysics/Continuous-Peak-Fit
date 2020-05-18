#!/usr/bin/env python

__all__ = ['execute']

import os
import sys

import numpy as np
import importlib.util


# debug = True  # If True will do the plotting tests
# refine = True  # If True refine Fourier fit for d, h,w bg separately
# iterations = 1  # Number of iterations of refinement if refine=True
# SaveAll = 1 # Save all the fit images!

# FIX ME: we need a better way of loading the output options.
# import FPF_WriteMultiFit as wr
# import FPF_WritePolydefixED as ED

# FIX ME:
# need to check that dependencies are present.
# models that are required for running (and are not listed above):
# h5py   -- sudo apt-get update; sudo apt-get install python-h5py
#        -- this maybe a dependency of fabio and so installed with pyFAI but I have not checked.
# pyFAI  -- sudo apt-get install pyfai
# Fabio  -- but this appears to be installed with pyFAI.
# lmfit  -- 
# pillow -- 


# Needed for JSON to save fitted parameters. 
# copied from https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array#24375113
# on 13th November 2018
def json_numpy_serialzer(o):
    """ Serialize numpy types for json

    Parameters:
        o (object): any python object which fails to be serialized by json

    Example:

        >>> import json
        >>> a = np.array([1, 2, 3])
        >>> json.dumps(a, default=json_numpy_serializer)

    """
    numpy_types = (
        np.bool_,
        # np.bytes_, -- python `bytes` class is not json serializable     
        # np.complex64,  -- python `complex` class is not json serializable  
        # np.complex128,  -- python `complex` class is not json serializable
        # np.complex256,  -- special handling below
        # np.datetime64,  -- python `datetime.datetime` class is not json serializable
        np.float16,
        np.float32,
        np.float64,
        # np.float128,  -- special handling below
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        # np.object_  -- should already be evaluated as python native
        np.str_,
        np.timedelta64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.void,
    )

    if isinstance(o, np.ndarray):
        return o.tolist()
    elif isinstance(o, numpy_types):
        return o.item()
    elif isinstance(o, np.float128):
        return o.astype(np.float64).item()
    # elif isinstance(o, np.complex256): -- no python native for np.complex256
    #     return o.astype(np.complex128).item() -- python `complex` class is not json serializable 
    else:
        raise TypeError("{} of type {} is not JSON serializable".format(repr(o), type(o)))


def execute(settings_file=None, inputs=None, outputs=None):
    """
    :param settings_file:
    :param inputs:
    :param outputs
    :return:
    """
    # Inputs ###
    # and setup

    if settings_file:
        try:
            print('  ')
            print('The name of the settings file is: ', settings_file)

            if not str.endswith(settings_file, '.py'):
                settings_file = settings_file + '.py'
            spec = importlib.util.spec_from_file_location("settings_module", settings_file)
            FitSettings = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(FitSettings)

            # if str.endswith(settings_file, '.py'):
            #     settings_file = settings_file.strip('.py')
            # FitSettings = __import__(settings_file)
            FitParameters = dir(FitSettings)
        except:
            raise FileNotFoundError

    elif inputs:
        raise NotImplementedError

    else:
        raise ValueError('No inputs provided')

    # List all the required parameters.
    # NOTE: No doubt this will change hence denoted as a list that is worked over.
    # 1. Have option to have list of files rather than sequence. This list needs to be called diff_files.

    # Inputs:
    #   calibration type. (required)
    #   data - in form of:
    #       - file listing data files. (with or without directory)
    #       - file listing data file numbers (needs directory, file name, number digits, ending)
    #       - beginning and and numbers for files (needs directory, file name, number digits, ending)

    RequiredList = ['Calib_type',
                    # 'Calib_data',        # Removed because this is not strictly required.
                    # 'Calib_param',       # Now added to calibration function file.
                    # 'Calib_pixels',      # this is only needed for GSAS-II and should be read from image file. FIX
                    # ME: add to GSAS-II inputfile.
                    # 'Calib_mask',        # a mask file is not strictly required.
                    'datafile_directory',
                    'datafile_Basename',
                    'datafile_Ending',
                    # 'datafile_StartNum',  # now optionally replaced by datafile_Files
                    # 'datafile_EndNum',    # now optionally replaced by datafile_Files
                    'datafile_NumDigit',
                    'AziBins',
                    'fit_orders',
                    # 'Output_type',		   # should be optional
                    # 'Output_NumAziWrite',  # should be optional
                    # 'Output_directory']	   # should be optional
                    ]
    AlternativesList = [
        [['datafile_StartNum', 'datafile_EndNum'], ['datafile_Files']]
    ]

    # load the detector functions here.
    if not 'Calib_type' in FitParameters:
        sys.exit(
            'There is no ''Calib_type'' in the settings. The fitting cannot proceed until a recognised '
            'Calibration type is present.')
    else:
        # calib_mod = 'FPF_' + FitSettings.Calib_type + 'Functions'
        calib_mod = 'cpf.' + FitSettings.Calib_type + 'Functions'
    try:
        det = importlib.import_module(calib_mod)
        # det = __import__(calib_mod)
    except ImportError:
        raise ImportError(
            'The ''Calib_type'' ' + FitSettings.Calib_type + ' is not recognised; the file ' + calib_mod +
            ' does not exist. Check if the calibration type exists.')
        # FIX ME: List possible calibration types.

    # Add detector requirements to the Required List of parameters - which has been read from the Detector functions.
    RequiredList = RequiredList + det.Requirements()

    # load output function(s) and check output options are present and valid
    # FIX ME: inly works for single output should work for multiple.
    print(outputs)
    if outputs is not None:
        output_mod = 'cpf.Write' + outputs
        
    elif 'Output_type' in FitParameters:
        output_mod = 'cpf.Write' + FitSettings.Output_type
        
    else:
        print('No output stype. Raise error')
        stop
        
    try:
        wr = importlib.import_module(output_mod)
        # wr = __import__(output_mod)
    except ImportError:
        raise ImportError(
            'The ''Output_type'' ' + FitSettings.Output_type + ' is not recognised; the file ' + output_mod +
            ' does not exist.')
        # FIX ME: List possible output types.

        RequiredList = RequiredList + wr.Requirements()

    # check all the required values are present.
    all_present = 1
    # properties of the data files.
    for x in range(len(RequiredList)):
        if RequiredList[x] in FitParameters:
            print('Got: ', RequiredList[x])
        else:
            print('The settings file requires a parameter called  \'', RequiredList[x], '\'')
            all_present = 0

    possible = [[[], []] * len(AlternativesList)]
    for x in range(len(AlternativesList)):
        for y in range(2):
            for z in range(len(AlternativesList[x][y])):
                if AlternativesList[x][y][z] in FitParameters:

                    possible[x][y].append(1)
                else:
                    possible[x][y].append(0)

    # exit if all parameters are not present
    if all_present == 0:
        sys.exit(
            "The highlighted settings are missing from the input file. Fitting cannot proceed until they are all present.")



    '''
    ## Load files ##
    '''

    # calibration parameter file #
    parms_dict = det.GetCalibration(os.path.abspath(FitSettings.Calib_param))





    
    # Write the output files.
    # This is there so that the output file can be run separately from the fitting of the diffraction patterns. 
    # FIX ME: Need to add options for more than one output.
    print('\nWrite output file(s)')
    wr.WriteOutput(FitSettings, parms_dict)


