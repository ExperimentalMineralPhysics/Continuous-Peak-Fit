#!/usr/bin/env python

__all__ = ['execute', 'write_output']

import json
import os
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import multiprocessing as mp
from itertools import product
import importlib.util
import cpf.DioptasFunctions as DioptasFunctions
import cpf.GSASIIFunctions as GSASIIFunctions
import cpf.MedFunctions as MedFunctions
from cpf.XRD_FitSubpattern import FitSubpattern, peak_string
from cpf.Cosmics import cosmicsimage
from cpf.IO_functions import json_numpy_serialzer, FileList

np.set_printoptions(threshold=sys.maxsize)


# FIX ME: Need to add complexity here.
# Need option to check required functions are present if adding new output format.
# Also need to do something similar with data class
def register_default_formats():
    """
    Load all available output modules
    :return:
    """
    # FIX ME: We could add extra checks here to make sure teh required functions exist in each case.
    # FIX ME: Really these should live in their own sub-folder
    output_list = ['CoefficientTable', 'DifferentialStrain', 'MultiFit', 'Polydefix', 'PolydefixED']
    new_module = {}
    for output_module in output_list:
        module = __import__('cpf.Write'+output_module, fromlist=[None])
        new_module[output_module] = module
    return new_module


# Load potential output formats
output_methods_modules = register_default_formats()

'''
Fabio function for reference currently
def importer(module_name):
    module = __import__(module_name)
    # returns the leaf module, instead of the root module
    names = module_name.split(".")
    names.pop(0)
    for name in names:
        module = getattr(module, name)
    return module
'''


def get_output_options(output_type):
    """
    Check if input is string or list of strings
    :param output_type: string or list of strings
    :return: list of strings
    """
    output_mod_type = []
    if isinstance(output_type, str):
        output_mod_type.append(output_type)
    else:
        output_mod_type = output_type
    return output_mod_type


def detector_factory(calibration_param, calibration_type, fit_settings=None):
    """
    Factory function to provide appropriate class for data dependent on type.
    Currently supported options are Dioptas, GSASII and Med.
    :param fit_settings:
    :param calibration_param:
    :param calibration_type:
    :return:
    """
    if calibration_type == "Dioptas":
        DetectorClass = DioptasFunctions.DioptasDetector
        return DetectorClass(calibration_param, fit_settings)
    if calibration_type == "GSASII":
        DetectorClass = GSASIIFunctions.GSASIIDetector
        return DetectorClass(calibration_param, fit_settings)
    if calibration_type == "Med":
        DetectorClass = MedFunctions.MedDetector
        return DetectorClass(calibration_param, fit_settings)
    else:
        raise Exception("Unrecognized calibration type.")
    return DetectorClass




def initiate(settings_file=None, inputs=None, out_type=None, initiateData=True, **kwargs):
    """
    Run checks on input files, initiate data class and check output options
    :param settings_file:
    :param inputs:
    :return FitParameters:
    :return FitSettings:
    """
    # Check we can find and load settings file
    if settings_file:
        try:
            print('\n The name of the settings file is: ', settings_file)
            if not str.endswith(settings_file, '.py'):
                settings_file = settings_file + '.py'
            spec = importlib.util.spec_from_file_location("settings_module", settings_file)
            FitSettings = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(FitSettings)
            FitParameters = dir(FitSettings)
            FitSettings.inputfile = settings_file
        except FileNotFoundError:
            raise FileNotFoundError
    elif inputs:
        raise NotImplementedError
    else:
        raise ValueError('No inputs provided')

    # Data directory and file validation: check they exist.
    if os.path.exists(FitSettings.datafile_directory) is False:
        raise ImportError('The data directory ' + FitSettings.datafile_directory + ' does not exist.')
    else:
        print('The data directory exists.')
    diff_files, n_diff_files = FileList(FitParameters, FitSettings)
    for j in range(n_diff_files):
        if os.path.isfile(diff_files[j]) is False:
            raise ImportError('The data file ' + diff_files[j] + ' does not exist. Check the input file.')
    print('All the data files exist.')

    # Load the detector class here to access relevant functions and check required parameters are present
    if 'Calib_type' not in FitParameters:
        raise ValueError("There is no 'Calib_type' in the settings. The fitting cannot proceed until a recognised "
                         "calibration type is present.")
    else:
        data_class = detector_factory(FitSettings.Calib_param, FitSettings.Calib_type)

    # FIX ME: This doesn't seem to be used, if it should be this needs moving to class structure.
    AlternativesList = [[['datafile_StartNum', 'datafile_EndNum'], ['datafile_Files']]]
    OptionalList = ['datafile_Step']  # FIX ME: This doesn't seem to be used
    possible = [[[], []] * len(AlternativesList)]
    for x in range(len(AlternativesList)):
        for y in range(2):
            for z in range(len(AlternativesList[x][y])):
                if AlternativesList[x][y][z] in FitParameters:
                    possible[x][y].append(1)
                else:
                    possible[x][y].append(0)
    # exit if all parameters are not present

    # Load output function(s) and check output options are present and valid
    if out_type is not None:
        output_mod = get_output_options(out_type)
    elif 'Output_type' in FitParameters:
        output_mod = get_output_options(FitSettings.Output_type)
    else:
        raise ValueError("No output type. Add 'Output_type' to input file or specify 'out_type' in command.")

    # Check output format exists
    for mod in output_mod:
        if mod not in output_methods_modules.keys():
            raise ImportError("The 'Output_type' " + mod + " is not recognised; the file '" + mod + "' does not "
                                                                                                    "exist. Check if "
                                                                                                    "the calibration "
                                                                                                    "type exists.")

    return FitSettings, FitParameters, data_class



def setrange(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):
    
    
    FitSettings, FitParameters, new_data = initiate(settings_file, inputs=inputs)
    
    # search over the first file only
    # strip file list to first file
    if not 'datafile_Files' in FitParameters:
        FitSettings.datafile_EndNum = FitSettings.datafile_StartNum
    else:
        FitSettings.datafile_Files = FitSettings.datafile_Files[0]
    
    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]
        
    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])
           
    FitSettings.fit_orders = orders_tmp
    
    execute(FitSettings=FitSettings, FitParameters=FitParameters, inputs=new_data, 
            debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
            parallel=parallel,
            mode='setrange')


def ordersearch(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, search_parameter='height', search_over=[0,20], subpattern='all', search_peak=0, search_series=['fourier', 'spline'], **kwargs):
    """
    :param settings_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """
    FitSettings, FitParameters, new_data = initiate(settings_file, inputs=inputs)

    #force it to write the required output type.
    FitSettings.Output_type = 'DifferentialStrain'    

    # search over the first file only
    # strip file list to first file
    if not 'datafile_Files' in FitParameters:
        FitSettings.datafile_EndNum = FitSettings.datafile_StartNum
    else:
        FitSettings.datafile_Files = FitSettings.datafile_Files[0]
    
    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]
    # make new order search list
    if isinstance(search_over,list) and len(search_over) == 2:
        search = list(range(search_over[0], search_over[1]))
    else:
        search = [int(x) for x in str(search_over)]
    
    order_search = []
    for i in range(len(subpats)):
        
        for j in range(len(search_series)):
            tmp_order = FitSettings.fit_orders[subpats[i]]
            
            for k in range(len(search)):
                orders_s = deepcopy(tmp_order)
                if search_parameter != 'background':
                    orders_s['peak'][search_peak][search_parameter] = search[k]  
                    orders_s['peak'][search_peak][search_parameter+'-type'] = search_series[j]  
                    
                else:
                    orders_s['background'][search_peak] = search[k]
                orders_s['note'] = search_parameter+'='+str(search[k])+'; type='+search_series[j]
                order_search.append(orders_s)
    FitSettings.fit_orders = order_search
    
    execute(FitSettings=FitSettings, FitParameters=FitParameters, inputs=new_data, 
            debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
            parallel=parallel)
    
    #write a differential strain output file
    #FIX ME: using debug as a switch to get all info in file. 
    #should probably use another switch
    write_output(settings_file, FitSettings=FitSettings, FitParameters=FitParameters, debug=True, outtype='DifferentialStrain')



def write_output(settings_file=None, FitSettings=None, FitParameters=None, parms_dict=None, out_type=None,
                 det=None, use_bounds=False, differential_only=False, debug=False, **kwargs):
    """
    :param FitParameters:
    :param use_bounds:
    :param differential_only:
    :param settings_file:
    :param FitSettings:
    :param parms_dict:
    :param out_type:
    :param det:
    :return:
    """
    # Fail gracefully
    if parms_dict is None and settings_file is None and FitSettings is None and FitParameters is None:
        raise ValueError("Either the settings file or the parameter dictionary need to be specified.")
    # If no params_dict then initiate. Check all the output functions are present and valid.
    if parms_dict is None and FitSettings==None and FitParameters==None:
        FitSettings, FitParameters, new_data = initiate(settings_file, outtype=out_type)
        parms_dict = new_data.get_calibration(os.path.abspath(FitSettings.Calib_param))
    elif parms_dict is None: 
        _, _, new_data = initiate(settings_file, outtype=out_type)
        parms_dict = new_data.get_calibration(os.path.abspath(FitSettings.Calib_param))
    if out_type is not None:
        output_mod = get_output_options(out_type)
    elif 'Output_type' in FitParameters:
        output_mod = get_output_options(FitSettings.Output_type)
    else:
        raise ValueError("No output type. Add 'Output_type' to input file or specify 'out_type' in command.")

    for mod in output_mod:
        print('\nWrite output file(s) using', mod)
        wr = output_methods_modules[mod]
        wr.WriteOutput(FitSettings, parms_dict, differential_only=differential_only, debug=debug)



def execute(settings_file=None, FitSettings=None, FitParameters=None, inputs=None, debug=False, refine=True, save_all=False, 
            propagate=True, iterations=1,
            track=False, parallel=True, 
            mode='fit',
            **kwargs):
    """
    :param track:
    :param propagate:
    :param save_all:
    :param settings_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """
    if FitSettings == None:
        FitSettings, FitParameters, new_data = initiate(settings_file)
    else:
        new_data = inputs

    # Define locally required names
    temporary_data_file = 'PreviousFit_JSON.dat'

    # Get list of diffraction patterns #
    diff_files, n_diff_files = FileList(FitParameters, FitSettings)

    # Initiate data in data class with appropriate detector class
    if 'Calib_mask' in FitParameters:
        use_mask = FitSettings.Calib_mask
    else:
        use_mask = None
        
    # new_data.fill_data(os.path.abspath(FitSettings.Calib_data), FitSettings.Calib_detector, debug=debug,
    #                    calibration_mask=use_mask)
    if 'Calib_data' in FitParameters:
        new_data.fill_data(os.path.abspath(FitSettings.Calib_data), FitSettings.Calib_detector, debug=debug,
                           calibration_mask=use_mask)
    else:
        new_data.fill_data(os.path.abspath(diff_files[0]), FitSettings.Calib_detector, debug=debug,
                           calibration_mask=use_mask)
            
    # Get calibration parameter file
    parms_dict = new_data.parameters

    # Get calibration for image set.
    # Get array for azimuth, two theta and d-space the same size as the image.
    azimu = new_data.azm
    twotheta = new_data.tth
    dspace = new_data.dspace
    # print(azimu, twotheta, dspace)

    ### Replace this with call through to class structure
    # plot calibration file
    if 0:#debug and 'Calib_data' in FitParameters:  #FIX ME: SAH - I have jsut removed this so that the debug runs without the calibration file.
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
        intens = new_data.ImportImage(diff_files[j], mask=azimu.mask, debug=debug)
        ### Replace this with call through to class structure
        ## FIX ME: the mask call is needed here to pass the mask in. but should it inherit the mask from the preceeding data?

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
                intens = ma.array(intens, mask=msk)
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

        # Get previous fit (if it exists and is required)
        if os.path.isfile(temporary_data_file) and propagate is True and mode=='fit':
            # Read JSON data from file
            print('Loading previous fit results from %s' % temporary_data_file)
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)

        # Switch to save the first fit in each sequence.
        if j == 0 or save_all is True:
            SaveFigs = 1
        else:
            SaveFigs = 0

        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        # Get number of sub-patterns to be fitted
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
            if 'previous_fit' in locals() and mode=='fit':
                params = previous_fit[i]
            else:
                params = []
                
            # Track the position of the peak centroid
            # FIX ME; This is crude - the range doesnt change width. so can't account for massive change in stress.
            # But does it need to?
            if track is True and 'previous_fit' in locals():
                print('tthRange', tthRange)
                mid = []
                for k in range(len(params['peak'])):
                    mid.append(params['peak'][k]['d-space'][0])
                ### Replace this with call through to class structure?
                cent = new_data.conversion(np.sum(mid)/(k+1), parms_dict, reverse=1, azm=None)
                tthRange[0] = tthRange[0] - ((tthRange[1]+tthRange[0])/2) + cent
                tthRange[1] = tthRange[1] - ((tthRange[1]+tthRange[0])/2) + cent
                print('move by:', ((tthRange[1]+tthRange[0])/2) + cent)
                print('tthRange', tthRange)

                # The PeakPositionSelections are only used if the fits are not being propagated
                if 'PeakPositionSelection' in orders:
                    for k in range(len(orders['PeakPositionSelection'])):
                        orders['PeakPositionSelection'][k][2] = orders['PeakPositionSelection'][k][2] - \
                                                                ((tthRange[1]+tthRange[0])/2) + cent

            # DISCUSS WITH SIMON Can replace with sending new_data through to subpattern call and add functions to
            # class to find subarrays - will depend on parallelisation though....

            # Get sub-sections of data to pass
            subpat = np.where((twotheta >= tthRange[0]) & (twotheta <= tthRange[1]))
            twotheta_sub = twotheta[subpat]
            dspacing_sub = dspace[subpat]
            azimu_sub = azimu[subpat]
            intens_sub = intens[subpat]
            
            # Mask the subpattern by intensity if called for
            if 'Imax' in orders:
                intens_sub = ma.masked_outside(intens_sub, 0, int(orders['Imax']))
                azimu_sub = ma.array(azimu_sub, mask=intens_sub.mask)
                twotheta_sub = ma.array(twotheta_sub, mask=intens_sub.mask)
                dspacing_sub = ma.array(dspacing_sub, mask=intens_sub.mask)
                
            if mode=='setrange':
                #FIX ME: this plots the input data. perhaps it should have its own switch rather than being subservient to Debug. 
                # It is not a debug it is a setup thing.
                #plot the data and the mask.
                fig_1 = new_data.plot(twotheta_sub, azimu_sub, intens_sub, plottype='mask', name=peak_string(orders))
                
                # save figures without overwriting old names
                filename = os.path.splitext(os.path.basename(diff_files[j]))[0]
                filename = filename + peak_string(orders, fname=True)
                
                i = 0
                if os.path.exists('{}.png'.format(filename)):
                    i+=1
                while os.path.exists('{}_{:d}.png'.format(filename, i)):
                    i += 1
                if i == 0:
                    fig_1.savefig('{}_mask.png'.format(filename))
                else:
                    fig_1.savefig('{}_{:d}.png'.format(filename, i))
                
                if debug:
                    plt.show()
                plt.close()
                    
                step=1000

            # fit the subpattern
            # tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params,
            #                     DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
            #                     iterations=iterations, fnam=diff_files[j])
            #tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds,
            #                    DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
            #                    iterations=iterations, fnam=diff_files[j])
            #Fitted_param.append(tmp[0])
            #lmfit_models.append(tmp[1])

            else:
                if parallel is True:#setup parallel version
                    #parallel_pile.append(([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds))
                
                    kwargs = {'SaveFit': SaveFigs, 'debug': debug, 'refine': refine, 'iterations': iterations, 'fnam': diff_files[j]}
                    args = ([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, new_data, orders, params, bounds)
                    parallel_pile.append((args, kwargs))
                    
                    #parallel_pile.append(([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, bounds))
                    
                else: #non-parallel version.
                    # fit the subpattern
                    # tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params,
                    #                     DetFuncs=calib_mod, SaveFit=SaveFigs, debug=debug, refine=refine,
                    #                     iterations=iterations, fnam=diff_files[j])
                    tmp = FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, new_data, orders,
                                        params, bounds, SaveFit=SaveFigs, debug=debug, refine=refine, iterations=iterations,
                                        fnam=diff_files[j])
                    Fitted_param.append(tmp[0])
                    lmfit_models.append(tmp[1])
                
        # write output files
        if mode=='fit':
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


    if mode=='fit':
        # Write the output files.
        write_output(settings_file, FitSettings=FitSettings, FitParameters=FitParameters, parms_dict=parms_dict, outtype=FitSettings.Output_type)

    if parallel is True:     
        p.close()

#def parallel_processing1(p):
#    return FitSubpattern(*p)

def parallel_processing(p):
    a, kw = p
    return FitSubpattern(*a, **kw)



    

if __name__ == '__main__':
    # Load settings fit settings file.
    sys.path.append(os.getcwd())
    settings_file = sys.argv[1]
    # print(settings_file)
    execute(inputs=None, settings_file=settings_file, debug=False, refine=True, save_all=False, iterations=1)
