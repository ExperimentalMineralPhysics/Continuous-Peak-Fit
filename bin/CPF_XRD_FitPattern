#!/usr/bin/env python

__all__ = ['execute']

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import importlib.util

# import FPF_WriteMultiFit as wr
from cpf.XRD_FitSubpattern import FitSubpattern

np.set_printoptions(threshold=sys.maxsize)

# debug = True  # If True will do the plotting tests
# refine = True  # If True refine Fourier fit for d, h,w bg separately
# iterations = 1  # Number of iterations of refinement if refine=True


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


def execute(settings_file=None, inputs=None, debug=False, refine=True, iterations=1):
    """
    :param settings_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
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
    if 'Output_type' in FitParameters:
        # output_mod = 'FPF_Write' + FitSettings.Output_type
        output_mod = 'cpf.Write' + FitSettings.Output_type
        try:
            wr = importlib.import_module(output_mod)
            # wr = __import__(output_mod)
        except ImportError:
            raise ImportError(
                'The ''Output_type'' ' + FitSettings.Output_type + ' is not recognised; the file ' + output_mod +
                ' does not exist. Check if the calibration type exists.')
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

    # Detector check.
    # Check if the detector is required and if so if it is present and recognised.
    if 'Calib_detector' in FitParameters:
        det_name = FitSettings.Calib_detector
    else:
        det_name = None
    FitSettings.Calib_detector = det.DetectorCheck(os.path.abspath(FitSettings.Calib_data), det_name)

    # define locally required names.
    temporary_data_file = 'PreviousFit_JSON.dat'

    # backg_type = 'flat'
    # backg = [[4]]

    # Main code #

    '''
    ## Load files ##
    '''

    # calibration parameter file #
    parms_dict = det.GetCalibration(os.path.abspath(FitSettings.Calib_param))

    '''
    ### Sort background out ####
    '''
    '''
    if backg_type == 'coeffs':
        # coefficients provided, figure out height using values
        tempbackg = []
        bg_order = []
        if isinstance(backg, list):
            for val in backg:
                if not isinstance(val, list):
                    tempbackg.append([val])
                    bg_order.append(0)
                else:
                    tempbackg.append(val)
                    bg_order.append((len(val) - 1) / 2)
        else:
            tempbackg = [[backg]]
        backg = tempbackg
    
    elif backg_type == 'order':
        # orders requested, figure out height using best guess
        print('TBD')
        tempbackg = []
        bg_order = []
        if isinstance(backg, list):
            if len(backg) == 2:
                nfourier = backg[0]
                npoly = backg[1]
            else:
                nfourier = backg[0]
                npoly = 0
        else:
            nfourier = backg
            npoly = 0
        for i in range(npoly + 1):
            tempbackg_f = [1 for j in range(nfourier + 1)]
            tempbackg.append(tempbackg_f)
            bg_order.append(nfourier)
        backg = tempbackg
    
    elif backg_type == 'flat':
        backg = [[backg]]
    
    print(backg, 'backg')
    '''

    # diffraction patterns #
    if not 'datafile_Files' in FitParameters:
        n_diff_files = FitSettings.datafile_EndNum - FitSettings.datafile_StartNum + 1
        diff_files = []
        for j in range(n_diff_files):
            # make list of diffraction pattern names

            # make number of pattern
            n = str(j + FitSettings.datafile_StartNum).zfill(FitSettings.datafile_NumDigit)

            # append diffraction pattern name and directory
            diff_files.append(os.path.abspath(
                FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n +
                FitSettings.datafile_Ending))
    elif 'datafile_Files' in FitParameters:
        n_diff_files = len(FitSettings.datafile_Files)
        diff_files = []
        for j in range(n_diff_files):
            # make list of diffraction pattern names

            # make number of pattern
            n = str(FitSettings.datafile_Files[j]).zfill(FitSettings.datafile_NumDigit)

            # append diffraction pattern name and directory
            diff_files.append(os.path.abspath(
                FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n +
                FitSettings.datafile_Ending))
    else:
        n_diff_files = len(FitSettings.datafile_Files)

    # get calibration for image set.
    # get array for azimuth, two theta and d-space the same size as the image.
    azimu = det.GetAzm(os.path.abspath(FitSettings.Calib_data), parms_dict, FitSettings.Calib_pixels,
                       FitSettings.Calib_detector)
    twotheta = det.GetTth(os.path.abspath(FitSettings.Calib_data), parms_dict, FitSettings.Calib_pixels,
                          FitSettings.Calib_detector)
    dspace = det.GetDsp(os.path.abspath(FitSettings.Calib_data), parms_dict, FitSettings.Calib_pixels,
                        FitSettings.Calib_detector)
    # FIX ME: Should be able to send diffraction patter, mask and calibration parameters to a function and get back
    # three masked arrays.

    # load calibration data/image -- to get intensities for mask.
    im = det.ImportImage(os.path.abspath(FitSettings.Calib_data), debug=debug)
    intens = ma.array(im)

    # create mask from mask file if present.
    # if not make all values valid
    if 'Calib_mask' in FitParameters:
        mask_file = os.path.abspath(FitSettings.Calib_mask)
        intens = det.GetMask(mask_file, intens, twotheta, azimu, os.path.abspath(FitSettings.Calib_data),
                             FitSettings.Calib_pixels, debug=debug)
    else:
        # print np.shape(intens), type(intens)
        # print np.shape(intens.mask), type(intens)
        ImMask = np.zeros_like(intens)
        intens = ma.array(intens, mask=ImMask)
        intens = ma.masked_outside(intens, 0, np.inf)
        # print np.shape(intens), type(intens)
        # print np.shape(intens.mask), type(intens)
        # stop

    # apply mask to arrays
    azimu = ma.array(azimu, mask=intens.mask)
    twotheta = ma.array(twotheta, mask=intens.mask)
    dspace = ma.array(dspace, mask=intens.mask)

    # plot calibration file
    if debug:
        fig = plt.figure()
        if FitSettings.Calib_type == 'Med':
            # plt.scatter(twotheta, azimu, s=intens/10, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0,
            # vmax = 1000)  #Better for energy dispersive data.  #FIX ME: need variable for vmax.
            for x in range(9):
                plt.plot(twotheta[x], intens[x] + 100 * x)  # Better for energy dispersive data. #FIX ME: need variable
                # for vmax.

        else:
            # plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 1000)
            # #Better for monochromatic data.             #FIX ME: need variable for vmax.
            plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin=0)
            # Better for monochromatic data.             #FIX ME: need variable for vmax.
            plt.colorbar()
        plt.show()
        plt.close()

    # Process the diffraction patterns #

    for j in range(n_diff_files):

        print('Process ', diff_files[j])

        # Get diffraction pattern to process.
        im = det.ImportImage(diff_files[j])

        # get intensity array and mask it
        intens = ma.array(im, mask=azimu.mask)

        # plot input file
        if debug:
            fig = plt.figure()
            if FitSettings.Calib_type == 'Med':
                # plt.scatter(twotheta, azimu, s=intens/10, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0,
                # vmax = 1000)  #Better for energy dispersive data.  #FIX ME: need variable for vmax.
                for x in range(9):
                    plt.plot(twotheta[x], intens[
                        x] + 100 * x)  # Better for energy dispersive data.  #FIX ME: need variable for vmax.

            else:
                # plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0,
                # vmax = 1000)  #Better for monochromatic data.             #FIX ME: need variable for vmax.
                plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin=0,
                            vmax=4000)  # FIX ME: need variable for vmax.
                plt.colorbar()
            plt.xlabel(r'2$\theta$')
            plt.ylabel('Azimuth')
            plt.show()
            plt.close()

            # quit()

        # get previous fit (if exists)
        # load temporary fit file - if it exists.
        if os.path.isfile(temporary_data_file):
            # Read JSON data from file
            print('Loading previous fit results from %s' % temporary_data_file)
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)

        # switch to save the first fit in each sequence.
        if j == 0:
            SaveFigs = 1
        else:
            SaveFigs = 0

        # Pass each sub-pattern to FPF_Fit_Subpattern for fitting in turn.
        # get number of sub-patterns to be fitted
        n_subpats = len(FitSettings.fit_orders)
        Fitted_param = []

        for i in range(n_subpats):

            tthRange = FitSettings.fit_orders[i]['range'][0]

            # get subsections of data to pass
            subpat = np.where((twotheta >= tthRange[0]) & (twotheta <= tthRange[1]))
            twotheta_sub = twotheta[subpat]
            dspacing_sub = dspace[subpat]
            azimu_sub = azimu[subpat]
            intens_sub = intens[subpat]

            if 'previous_fit' in locals():
                params = previous_fit[i]
            else:
                params = []

            if FitSettings.fit_orders:
                orders = FitSettings.fit_orders[i]
                orders['AziBins'] = FitSettings.AziBins
            # print orders

            # print(twotheta_sub.shape, dspacing_sub.shape, parms_dict, azimu_sub.shape, intens_sub.shape, orders,
            # params, calib_mod, SaveFigs)

            # fit the subpattern
            Fitted_param.append(
                FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params,
                              DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine, iterations=iterations))

        # write output files

        # store all the fit information as a JSON file.
        # filename, file_extension = os.path.splitext(diff_files[j])
        filename = os.path.splitext(os.path.basename(diff_files[j]))[0]
        filename = filename + '.json'
        # print(filename)
        with open(filename, 'w') as TempFile:
            # Write a JSON string into the file.
            json_string = json.dump(Fitted_param, TempFile, sort_keys=True, indent=2, default=json_numpy_serialzer)

        # #Write the fits to a temporary file
        # TempFilename = open(temporary_data_file, 'w')
        # # Write a JSON string into the file.
        # json_string = json.dumps(Fitted_param, TempFilename, sort_keys=True, indent=4, default=json_numpy_serializer)
        # #json_string = json.dumps(Fitted_param, TempFilename, sort_keys=True, indent=4, separators=(',', ': '))
        # TempFilename.write(json_string)

        # print json_string
        with open(temporary_data_file, 'w') as TempFile:

            # Write a JSON string into the file.
            json_string = json.dump(Fitted_param, TempFile, sort_keys=True, indent=2, default=json_numpy_serialzer)

    # Write the output files.
    # This is there so that the output file can be run separately from the fitting of the diffraction patterns.
    # FIX ME: Need to add options for more than one output.
    print('\nWrite out file(s)')
    wr.WriteOutput(FitSettings, parms_dict)


if __name__ == '__main__':

    # Load settings fit settings file.
    sys.path.append(os.getcwd())
    settings_file = sys.argv[1]
    print(settings_file)
    execute(inputs=None, settings_file=settings_file, debug=False, refine=True, iterations=1)
