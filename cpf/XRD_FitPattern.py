#!/usr/bin/env python

__all__ = ['execute', 'writeoutput']

import json
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import multiprocessing as mp
from itertools import product
import importlib.util
from cpf.XRD_FitSubpattern import FitSubpattern, peak_string
from cpf.Cosmics import cosmicsimage

np.set_printoptions(threshold=sys.maxsize)


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


def FileList(FitParameters, FitSettings):
    # from the Settings make a list of all the data files.
    # This function is called by the output writing scripts to make sure the file names are called consistently. 
    
    #define step
    if not 'datafile_Step' in FitParameters:
        Step = 1
    else:
        Step = FitSettings.datafile_Step
    
    # diffraction patterns -- make list of files
    if not 'datafile_Files' in FitParameters:
        n_diff_files = int(np.floor((FitSettings.datafile_EndNum - FitSettings.datafile_StartNum)/Step + 1))
        diff_files = []
            
        for j in range(n_diff_files):
            # make list of diffraction pattern names

            # make number of pattern
            n = str(FitSettings.datafile_StartNum+j*Step).zfill(FitSettings.datafile_NumDigit)

            # append diffraction pattern name and directory
            diff_files.append(os.path.abspath(
                FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n +
                FitSettings.datafile_Ending))
    elif 'datafile_Files' in FitParameters:
        n_diff_files = int(len(FitSettings.datafile_Files)/Step)
        diff_files = []
        for j in range(n_diff_files):
            # make list of diffraction pattern names

            # make number of pattern
            n = str(FitSettings.datafile_Files[j*Step]).zfill(FitSettings.datafile_NumDigit)

            # append diffraction pattern name and directory
            diff_files.append(os.path.abspath(
                FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n +
                FitSettings.datafile_Ending))
    else:
        n_diff_files = int(len(FitSettings.datafile_Files)/Step + 1)
        
    return diff_files, n_diff_files




def initiate(settings_file=None, inputs=None, outtype=None, initiateData=True, **kwargs):    
    """
    :param settings_file:
    :param inputs:
    :return FitParameters:
    :return FitSettings:
    """

    if settings_file:
        try:
            print('  ')
            print('The name of the settings file is: ', settings_file)

            if not str.endswith(settings_file, '.py'):
                settings_file = settings_file + '.py'
            spec = importlib.util.spec_from_file_location("settings_module", settings_file)
            FitSettings = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(FitSettings)

            FitParameters = dir(FitSettings)

            FitSettings.inputfile = settings_file
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
                    # 'AziBins',            # required based on detector type
                    'fit_orders',
                    # 'Output_type',		   # should be optional
                    # 'Output_NumAziWrite',  # should be optional
                    # 'Output_directory']	   # should be optional
                    ]
    AlternativesList = [
        [['datafile_StartNum', 'datafile_EndNum'], ['datafile_Files']]
    ]
    OptionalList = ['datafile_Step']

    # load the detector functions here.
    if not 'Calib_type' in FitParameters:
        raise ValueError(
            'There is no ''Calib_type'' in the settings. The fitting cannot proceed until a recognised '
            'Calibration type is present.')
    else:
        calib_mod = 'cpf.' + FitSettings.Calib_type + 'Functions'
    try:
        det = importlib.import_module(calib_mod)
    except ImportError:
        raise ImportError(
            'The ''Calib_type'' ' + FitSettings.Calib_type + ' is not recognised; the file ' + calib_mod +
            ' does not exist. Check if the calibration type exists.')
        # FIX ME: List possible calibration types.

    # Add detector requirements to the Required List of parameters - which has been read from the Detector functions.
    RequiredList = RequiredList + det.Requirements()

    # load output function(s) and check output options are present and valid
    output_mod = []
    if outtype is not None:
        if isinstance(outtype, str):
            output_mod.append('cpf.Write' + outtype)
        else:
            for x in range(len(outtype)):
                output_mod.append('cpf.Write' + outtype[x])
    elif 'Output_type' in FitParameters:
        if isinstance(FitSettings.Output_type, str):
            output_mod.append('cpf.Write' + FitSettings.Output_type)
        else:
            for x in range(len(FitSettings.Output_type)):
                output_mod.append('cpf.Write' + FitSettings.Output_type[x])
    else:
        raise ValueError('No output type. Add ''Output_type'' to input file or specify ''outtype'' in command.')

    try:
        for x in range(len(output_mod)):
            wr = importlib.import_module(output_mod[x])
            RequiredList = RequiredList + wr.Requirements()
    except ImportError:
        raise ImportError(
            'The ''Output_type'' ' + FitSettings.Output_type + ' is not recognised; the file ' + output_mod +
            ' does not exist. Check if the calibration type exists.')
        # FIX ME: List possible output types.

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
            "The highlighted settings are missing from the input file. Fitting cannot proceed until they are "
            "all present.")

    # Detector check: Check if the detector is required and if so if it is present and recognised.
    if 'Calib_detector' in FitParameters:
        det_name = FitSettings.Calib_detector
    else:
        det_name = None
    if not 'Calib_data' in FitParameters:
        # this is here to catch calls for the calibration data -- which is not used as calibration data.
        # it just has to be a data file of the right size and detector. If there is no calibration data 
        # we are now just replacing it with the first data file.
        # FIX ME: the whole calib_data could be removed if wanted.
        diff_files, n_diff_files = FileList(FitParameters, FitSettings)
        FitSettings.Calib_data = diff_files[0]
    FitSettings.Calib_detector = det.DetectorCheck(os.path.abspath(FitSettings.Calib_data), det_name)
    
    if initiateData is True:
        # Data directory and file validation: check they exist.
        if os.path.exists(FitSettings.datafile_directory) is False:
                raise ImportError(
                    'The data directory ' + FitSettings.datafile_directory + ' does not exist.')
        print('The data directory exists.')
        diff_files, n_diff_files = FileList(FitParameters, FitSettings)
        for j in range(n_diff_files):
            if os.path.isfile(diff_files[j]) is False:
                raise ImportError(
                    'The data file ' + diff_files[j] + ' does not exist. Check the input file.')
        print('All the data files exist.')
    
    
    return FitSettings, FitParameters



def writeoutput(settings_file=None, FitSettings=None, parms_dict=None, outtype=None, det=None, use_bounds=False,
                differential_only=False, **kwargs):
    """
    :param settings_file:
    :param FitSettings:
    :param parms_dict:
    :param outtype:
    :param det:
    :return:
    """
    # if no params_dict then initiate. Check all the output functions are present and valid.
    if parms_dict is None:
        FitSettings, FitParameters = initiate(settings_file, outtype=outtype, initiateData=False)
            
    if det is None:
        # get parms_dict from calibration parameter file #
        # load detector here
        det = importlib.import_module('cpf.' + FitSettings.Calib_type + 'Functions')
        parms_dict = det.GetCalibration(os.path.abspath(FitSettings.Calib_param))

    output_mod = []
    if outtype is not None:
        if isinstance(outtype, str):
            output_mod.append('cpf.Write' + outtype)
        else:
            for x in range(len(outtype)):
                output_mod.append('cpf.Write' + outtype[x])
    elif 'Output_type' in FitParameters:
        if isinstance(FitSettings.Output_type, str):
            output_mod.append('cpf.Write' + FitSettings.Output_type)
        else:
            for x in range(len(FitSettings.Output_type)):
                output_mod.append('cpf.Write' + FitSettings.Output_type[x])

    for x in range(len(output_mod)):
        wr = importlib.import_module(output_mod[x])
        print('\nWrite output file(s) using', output_mod[x])
        wr.WriteOutput(FitSettings, parms_dict, differential_only=differential_only)


def execute(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, **kwargs):
    """
    :param settings_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """

    FitSettings, FitParameters = initiate(settings_file, inputs)

    # load detector
    calib_mod = 'cpf.' + FitSettings.Calib_type + 'Functions'
    det = importlib.import_module('cpf.' + FitSettings.Calib_type + 'Functions')

    # define locally required names.
    temporary_data_file = 'PreviousFit_JSON.dat'

    # calibration parameter file #
    parms_dict = det.GetCalibration(os.path.abspath(FitSettings.Calib_param))

    # get list of diffraction patterns #
    diff_files, n_diff_files = FileList(FitParameters, FitSettings)
        
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
        intens = det.GetMask(FitSettings.Calib_mask, intens, twotheta, azimu, os.path.abspath(FitSettings.Calib_data),
                             FitSettings.Calib_pixels, debug=debug)
    else:
        ImMask = np.zeros_like(intens)
        intens = ma.array(intens, mask=ImMask)
        intens = ma.masked_outside(intens, 0, np.inf)

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
            plt.xlabel('Energy (keV)')
            plt.ylabel('Intensity')

        else:
            # plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 1000)
            # #Better for monochromatic data.             #FIX ME: need variable for vmax.
            plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin=0)
            # Better for monochromatic data.             #FIX ME: need variable for vmax.
            plt.xlabel(r'2$\theta$')
            plt.ylabel('Azimuth')
            plt.colorbar()
        plt.title('Calibration data')
        plt.draw()
        #plt.show()
        plt.close()


    # if parallel processing start the pool
    if parallel is True:  
        p = mp.Pool(processes=mp.cpu_count())



    # Process the diffraction patterns #

    for j in range(n_diff_files):

        print('Process ', diff_files[j])

        # Get diffraction pattern to process.
        im = det.ImportImage(diff_files[j])

    
        # get intensity array and mask it
        intens = ma.array(im, mask=azimu.mask) 
        # FIX ME: replace mask with image_prepare mask.
        if 'Image_prepare' in FitParameters:
            if 'cosmics' in FitSettings.Image_prepare:
                print('Remove Cosmics')
                #set defaults
                gain = 2.2
                sigclip = 1.5
                objlim = 3.0
                sigfrac = 0.3
                if bool(FitSettings.Image_prepare['cosmics']) == True:
                    if ('gain' in FitSettings.Image_prepare['cosmics']):
                        gain=FitSettings.Image_prepare['cosmics']['gain']       
                    if ('sigclip' in FitSettings.Image_prepare['cosmics']):
                        sigclip=FitSettings.Image_prepare['cosmics']['sigclip']  
                    if ('objlim' in FitSettings.Image_prepare['cosmics']):
                        objlim=FitSettings.Image_prepare['cosmics']['objlim']   
                    if ('sigfrac' in FitSettings.Image_prepare['cosmics']):
                        sigfrac=FitSettings.Image_prepare['cosmics']['sigfrac'] 
                    #FIX ME: use argparse or someway of passing any aguemnt into cosmics.      
                        
                test = cosmicsimage(intens, sigclip=sigclip, objlim=objlim, gain=gain,sigfrac=sigfrac)
                num = 2
                for i in range(num):
                    test.lacosmiciteration(True)
                    #print(asdf["nnew"], asdf["niter"])
                    test.clean()
                    msk = np.logical_or(azimu.mask, np.array(test.mask, dtype='bool'))
                intens = ma.array(im, mask=msk)
                azimu = ma.array(azimu, mask=intens.mask)
                twotheta = ma.array(twotheta, mask=intens.mask)
                dspace = ma.array(dspace, mask=intens.mask)
       
        # plot input file
        if debug:
            fig = plt.figure()
            if FitSettings.Calib_type == 'Med':
                # plt.scatter(twotheta, azimu, s=intens/10, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0,
                # vmax = 1000)  #Better for energy dispersive data.  #FIX ME: need variable for vmax.
                for x in range(9):
                    plt.plot(twotheta[x], intens[
                        x] + 100 * x)  # Better for energy dispersive data.  #FIX ME: need variable for vmax.
                plt.xlabel('Energy (keV)')
                plt.ylabel('Intensity')

            else:
                # plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0,
                # vmax = 1000)  #Better for monochromatic data.             #FIX ME: need variable for vmax.
                plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin=0,
                            vmax=4000)  # FIX ME: need variable for vmax.
                plt.xlabel(r'2$\theta$')
                plt.ylabel('Azimuth')
                plt.colorbar()
            plt.title(os.path.basename(diff_files[j]))
            plt.draw()
            plt.close()

        # get previous fit (if exists and required)
        if os.path.isfile(temporary_data_file) and propagate is True:
            # Read JSON data from file
            print('Loading previous fit results from %s' % temporary_data_file)
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)

        # switch to save the first fit in each sequence.
        if j == 0 or save_all is True:
            SaveFigs = 1
        else:
            SaveFigs = 0

        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        # get number of sub-patterns to be fitted
        n_subpats = len(FitSettings.fit_orders)
        Fitted_param = []
        lmfit_models = []
        parallel_pile = []


        for i in range(n_subpats):
            tthRange = FitSettings.fit_orders[i]['range'][0]
            if hasattr(FitSettings, 'fit_orders'):  # FitSettings.fit_orders:
                orders = FitSettings.fit_orders[i]
            if hasattr(FitSettings, 'AziBins'):
                orders['AziBins'] = FitSettings.AziBins
            if 'fit_bounds' in FitParameters:
                bounds = FitSettings.fit_bounds
            else:
                bounds = None
            if 'previous_fit' in locals():
                params = previous_fit[i]                    
            else:
                params = []
                
            #Track the position of the peak centroid 
            #FIX ME; This is crude - the range doesnt change width. so cant account for massive change in stress. But does it need to?
            if track is True and 'previous_fit' in locals():
                print('tthRange', tthRange)
                mid=[]
                for k in range(len(params['peak'])):
                    mid.append(params['peak'][k]['d-space'][0])
                cent = det.Conversion(np.sum(mid)/(k+1), parms_dict, reverse=1, azm=None)
                tthRange[0] = tthRange[0] - ((tthRange[1]+tthRange[0])/2) + cent
                tthRange[1] = tthRange[1] - ((tthRange[1]+tthRange[0])/2) + cent
                print('move by:' , ((tthRange[1]+tthRange[0])/2) + cent)
                print('tthRange', tthRange)

                # The PeakPositionSelections are only used if the fits are not being propagated
                if 'PeakPositionSelection' in orders:
                    for k in range(len(orders['PeakPositionSelection'])):
                        orders['PeakPositionSelection'][k][2] = orders['PeakPositionSelection'][k][2] - ((tthRange[1]+tthRange[0])/2) + cent

            # get subsections of data to pass
            subpat = np.where((twotheta >= tthRange[0]) & (twotheta <= tthRange[1]))
            twotheta_sub = twotheta[subpat]
            dspacing_sub = dspace[subpat]
            azimu_sub = azimu[subpat]
            intens_sub = intens[subpat]

            #Mask the subpattern by intensity if called for
            if 'Imax' in orders:
                intens_sub   = ma.masked_outside(intens_sub,0,int(orders['Imax']))
                azimu_sub    = ma.array(azimu_sub, mask=intens_sub.mask)
                twotheta_sub = ma.array(twotheta_sub, mask=intens_sub.mask)
                dspacing_sub = ma.array(dspacing_sub, mask=intens_sub.mask)
                

            if 0:#debug:
                #FIX ME: this plots the input data. perhaps it should have its own switch rather than being subservient to Debug. 
                # It is not a debug it is a setup thing.
                #plot the data and the mask.
                det.plot(twotheta_sub, azimu_sub, intens_sub, dtype='mask')

            # fit the subpattern
            # tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params,
            #                     DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
            #                     iterations=iterations, fnam=diff_files[j])
            
            #tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds,
            #                    DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
            #                    iterations=iterations, fnam=diff_files[j])
            #Fitted_param.append(tmp[0])
            #lmfit_models.append(tmp[1])


            if parallel is True:#setup parallel version
                #parallel_pile.append(([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds))
            
                kwargs = {'DetFuncs': calib_mod, 'SaveFit': SaveFigs, 'debug': debug, 'refine': refine, 'iterations': iterations, 'fnam': diff_files[j]}
                args = ([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds)
                parallel_pile.append((args, kwargs))
                
                #parallel_pile.append(([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds))
                
            else: #non-parallel version.
                # fit the subpattern
                # tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params,
                #                     DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
                #                     iterations=iterations, fnam=diff_files[j])
                tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds,
                                    DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
                                    iterations=iterations, fnam=diff_files[j])
                Fitted_param.append(tmp[0])
                lmfit_models.append(tmp[1])

                

  
        if parallel is True:            
            #print(mp.cpu_count())
            #print(len(parallel_pile)) 
            
            '''
            kwargs = {'DetFuncs': calib_mod, 'SaveFit': SaveFigs, 'debug': debug, 'refine': refine, 'iterations': iterations, 'fnam': diff_files[j]}           
            result = mp.Process(target=parallel_processing, args=parallel_pile, kwargs=kwargs)
            result.start()
            result.join()
            '''
            #p = mp.Pool(processes=mp.cpu_count())#, initializer=(SubpatternParallel_init), initargs=(FitSettings, FitParameters, params, diff_files[j], calib_mod, SaveFigs, debug, refine, iterations))) 
            #result = p.map(parallel_processing1, parallel_pile) #works for just args as parallelepile
            tmp = p.map(parallel_processing, parallel_pile) 
            for i in range(n_subpats):
                Fitted_param.append(tmp[i][0])
                lmfit_models.append(tmp[i][1])

        
        # write output files

        # store the fit parameters information as a JSON file.
        filename = os.path.splitext(os.path.basename(diff_files[j]))[0]
        filename = filename + '.json'
        # print(filename)
        with open(filename, 'w') as TempFile:
            # Write a JSON string into the file.
            json_string = json.dump(Fitted_param, TempFile, sort_keys=True, indent=2, default=json_numpy_serialzer)

        # if propagating the fits write them to a temporary file
        if propagate:
            # print json_string
            with open(temporary_data_file, 'w') as TempFile:
                # Write a JSON string into the file.
                json_string = json.dump(Fitted_param, TempFile, sort_keys=True, indent=2, default=json_numpy_serialzer)

    # Write the output files.
    writeoutput(settings_file, FitSettings=FitSettings, parms_dict=parms_dict, det=det, outtype=FitSettings.Output_type)

    if parallel is True:     
        p.close()

#def parallel_processing1(p):
#    return FitSubpattern(*p)

def parallel_processing(p):
    a, kw = p
    return FitSubpattern(*a, **kw)



def setup(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, **kwargs):
    """
    :param settings_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """

    FitSettings, FitParameters = initiate(settings_file, inputs)

    # load detector
    calib_mod = 'cpf.' + FitSettings.Calib_type + 'Functions'
    det = importlib.import_module('cpf.' + FitSettings.Calib_type + 'Functions')

    # define locally required names.
    temporary_data_file = 'PreviousFit_JSON.dat'

    # calibration parameter file #
    parms_dict = det.GetCalibration(os.path.abspath(FitSettings.Calib_param))

    # get list of diffraction patterns #
    diff_files, n_diff_files = FileList(FitParameters, FitSettings)
        
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
        intens = det.GetMask(FitSettings.Calib_mask, intens, twotheta, azimu, os.path.abspath(FitSettings.Calib_data),
                             FitSettings.Calib_pixels, debug=debug)
    else:
        ImMask = np.zeros_like(intens)
        intens = ma.array(intens, mask=ImMask)
        intens = ma.masked_outside(intens, 0, np.inf)

    # apply mask to arrays
    azimu = ma.array(azimu, mask=intens.mask)
    twotheta = ma.array(twotheta, mask=intens.mask)
    dspace = ma.array(dspace, mask=intens.mask)

    # Get diffraction pattern to process.
    im = det.ImportImage(diff_files[0])

    # get intensity array and mask it
    intens = ma.array(im, mask=azimu.mask) 
    # FIX ME: replace mask with image_prepare mask.
    if 'Image_prepare' in FitParameters:
        if 'cosmics' in FitSettings.Image_prepare:
            print('Remove Cosmics')
            #set defaults
            gain = 2.2
            sigclip = 1.5
            objlim = 3.0
            sigfrac = 0.3
            if bool(FitSettings.Image_prepare['cosmics']) == True:
                if ('gain' in FitSettings.Image_prepare['cosmics']):
                    gain=FitSettings.Image_prepare['cosmics']['gain']       
                if ('sigclip' in FitSettings.Image_prepare['cosmics']):
                    sigclip=FitSettings.Image_prepare['cosmics']['sigclip']  
                if ('objlim' in FitSettings.Image_prepare['cosmics']):
                    objlim=FitSettings.Image_prepare['cosmics']['objlim']   
                if ('sigfrac' in FitSettings.Image_prepare['cosmics']):
                    sigfrac=FitSettings.Image_prepare['cosmics']['sigfrac'] 
                #FIX ME: use argparse or someway of passing any aguemnt into cosmics.      
                                
            test = cosmicsimage(intens, sigclip=sigclip, objlim=objlim, gain=gain,sigfrac=sigfrac)
            num = 2
            for i in range(num):
                test.lacosmiciteration(True)
                #print(asdf["nnew"], asdf["niter"])
                test.clean()
                msk = np.logical_or(azimu.mask, np.array(test.mask, dtype='bool'))
            intens = ma.array(im, mask=msk)
   
    # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
    # get number of sub-patterns to be fitted
    n_subpats = len(FitSettings.fit_orders)

    for i in range(n_subpats):
        tthRange = FitSettings.fit_orders[i]['range'][0]
        if hasattr(FitSettings, 'fit_orders'):  # FitSettings.fit_orders:
            orders = FitSettings.fit_orders[i]
        
        # get subsections of data to pass
        subpat = np.where((twotheta >= tthRange[0]) & (twotheta <= tthRange[1]))
        twotheta_sub = twotheta[subpat]
        dspacing_sub = dspace[subpat]
        azimu_sub = azimu[subpat]
        intens_sub = intens[subpat]
        
        #Mask the subpattern by intensity if called for
        if 'Imax' in orders:
            intens_sub   = ma.masked_outside(intens_sub,0,int(orders['Imax']))
            azimu_sub    = ma.array(azimu_sub, mask=intens_sub.mask)
            twotheta_sub = ma.array(twotheta_sub, mask=intens_sub.mask)
            dspacing_sub = ma.array(dspacing_sub, mask=intens_sub.mask)
            
            
        
        filename = os.path.splitext(os.path.basename(diff_files[0]))[0] + '_'
        filename = filename + peak_string(orders, fname=True)    
            
        #FIX ME: this plots the input data. perhaps it should have its own switch rather than being subservient to Debug. 
        # It is not a debug it is a setup thing.
        #plot the data and the mask.
        fig_1 = det.plot(twotheta_sub, azimu_sub, intens_sub, dtype='mask', name=peak_string(orders))
        #fig_1 = det.plot(twotheta_sub, azimu_sub, intens_sub, dtype='mask', name=filename)
    
    
        # save figures without overwriting old names
        i = 0
        if os.path.exists('{}.png'.format(filename)):
            i+=1
        while os.path.exists('{}_{:d}.png'.format(filename, i)):
            i += 1
        if i == 0:
            fig_1.savefig('{}_mask.png'.format(filename))
            #fig_2.savefig('{}_CDF.png'.format(filename))
        else:
            fig_1.savefig('{}_{:d}.png'.format(filename, i))
            #fig_2.savefig('{}_{:d}.png'.format(filename, i))


if __name__ == '__main__':
    # Load settings fit settings file.
    sys.path.append(os.getcwd())
    settings_file = sys.argv[1]
    # print(settings_file)
    execute(inputs=None, settings_file=settings_file, debug=False, refine=True, save_all=False, iterations=1)
