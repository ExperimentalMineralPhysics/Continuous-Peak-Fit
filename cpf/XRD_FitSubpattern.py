#!/usr/bin/env python

__all__ = ['FitSubpattern']

# FPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

import os
import sys
import time
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from lmfit import Parameters, Model
from lmfit.model import save_modelresult #, load_modelresult
import importlib

import cpf.PeakFunctions as ff

np.set_printoptions(threshold=sys.maxsize)



def AnyTermsNull(obj_to_inspect, val_to_find="Null", indexpath="", clean=None):
    ''' This function accepts a nested dictionary and list as argument
        and iterates over all values of nested dictionaries and lists.
        If any of the values are "Null" it returns 0
    '''  
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021
    if clean == None:
        clean=1
    if isinstance(obj_to_inspect, dict):
        for key, value in obj_to_inspect.items():
            clean = AnyTermsNull(value, val_to_find, indexpath + f"['{key}']", clean=clean)
    if isinstance(obj_to_inspect, list):
        for key, value in enumerate(obj_to_inspect):
            clean = AnyTermsNull(value, val_to_find, indexpath + f"[{key}]", clean=clean)
    if isinstance(obj_to_inspect, str):
        if obj_to_inspect == val_to_find:
            #print(f"Value {val_to_find} found at {indexpath}")
            clean=0
    return clean



# FIX ME;
# The script should not need to know what the limits are - it should only see the data that it needs to fit.
def peak_string(orders, fname=False):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param fname: str
    :return: string listing peak names
    """
    # FIX ME: this will be a data-type function so should probably move somewhere else in the end.
    p_str = ''
    
    for x in range(len(orders['peak'])):
        if 'phase' in orders['peak'][x]:
            p_str = p_str + orders['peak'][x]['phase']
        else:
            p_str = p_str + "Peak"
        if fname is False:
            p_str = p_str + " ("
        if 'hkl' in orders['peak'][x]:
            p_str = p_str + str(orders['peak'][x]['hkl'])
        else:
            p_str = p_str + str(x+1)
        if fname is False:
            p_str = p_str + ")"
                
        if x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
            if fname is False:
                p_str = p_str + " & "
            else:
                p_str = p_str + '_'

    return p_str


def run_initial_fit(azi, y_data, y_data_errs, params, param_str, comp, order, extras, coef_type=0, method='leastsq', symm=1,
                    value=None):
    """
    Initiate parameter class and perform first Fourier fit
    :param azi: x-data arr float
    :param y_data: y_data arr float
    :param y_data_errs: y_data error arr float
    :param params: lmfit Parameter object
    :param param_str: base string to identify components str
    :param comp: string to add to base string to identify components str
    :param order: determines number of Fourier coefficients int
    :param method: lmfit method to use str
    :param extras: list of [[min, max, expr], ind_vars], expression for bounds on parameters and independent variable
    added to params for expression
    :param symm:
    :param value:
    :return: lmfit Parameter class
    """
    new_limits = extras
    new_min, new_max, expr = new_limits
    
    if isinstance(order, list):
        orderr = max(order)
    else:
        orderr = order
    
    if value is None:
        value = np.zeros(orderr*2+1)
        value[0] = np.mean(y_data)
    # print('Initiate parms', expr)
    params = ff.initiate_params(params, param_str, comp, coef_type, orderr, limits=[new_max, new_min], expr=expr, ind_vars=azi,
                                value=value)
    params = ff.unvary_params(params, param_str, comp)  # set other parameters to not vary
    params = ff.unvary_part_params(params, param_str, comp, order)  # set part of parameter not  to vary
    params.pretty_print()
    fout = ff.Fourier_fit(azimu=azi, ydata=np.array(y_data), param=params, errs=np.array(y_data_errs),
                          param_str=param_str + '_' + comp, symm=symm, fit_method=method)
    params = fout.params

    return params

'''
def initial_component_fits(master_params, param_str, data_arrays, comp_lists, order_list, symm, limits, pfixed=None,
                           method='leastsq'):
    """
    Perform initial Fourier fits to each component data in turn
    :param master_params: lmfit Parameter class object
    :param param_str: start string to label parameters
    :param data_arrays: list of data arrays for use in component fits
    :param comp_lists: list of lists for parameters
    :param order_list: list of orders for components
    :param symm : number defining the symmetry
    :param limits: array of limits
    :param pfixed: either False or true
    :param : bounds
    :return: list of lists of each component with coefficients and errors
    """

    newAziChunks, newd0, newd0Err, newHall, newHallErr, newWall, newWallErr, newPall, newPallErr = data_arrays
    dfour, hfour, wfour, pfour = comp_lists
    d0_order, h_order, w_order, p_order = order_list
    d_lims, h_lims, w_lims, p_lims = limits
    # initiate symmetry
    comp = 's'
    master_params = ff.initiate_params(master_params, param_str, comp, 0, limits=[10, 1], value=symm)
    # Add data arrays to master_params so can use in expression
    # master_params._asteval.symtable['newAzichunks'] = newAziChunks
    # print('Master parms now:')
    # master_params.pretty_print()

    comp = 'd'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newd0),
                                    y_data_errs=np.array(newd0Err), params=master_params, param_str=param_str,
                                    comp=comp, order=d0_order, extras=d_lims, method=method,
                                    symm=1)
    dfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    
    comp = 'h'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newHall),
                                    y_data_errs=np.array(newHallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=h_order, extras=h_lims,
                                    method=method, symm=symm)
    hfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'w'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newWall),
                                    y_data_errs=np.array(newWallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=w_order, extras=w_lims,
                                    method=method, symm=symm)
    wfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'p'
    if pfixed is False:
        master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newPall),
                                        y_data_errs=np.array(newPallErr), params=master_params, param_str=param_str,
                                        comp=comp, order=p_order, extras=p_lims,
                                        method=method, symm=symm)
    else:
        master_params = ff.initiate_params(master_params, param_str=param_str, comp=comp, coef_type=0, limits=[p_lims[0],p_lims[1]], value=pfixed)
    pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()
    # print('Initial fit output parms')
    # print(dfour, hfour, wfour, pfour)
    # master_params.pretty_print()

    return master_params, dfour, hfour, wfour, pfour
'''

def update_component_fits(master_params, intens, azimu, twotheta, param_str, comp_lists, conversion_factor, bgguess,
                          Guesses, o=None, pfixed=None):
    """
    Perform component fits to data
    :param master_params: lmfit Parameter class object
    :param intens: intensity data array
    :param azimu: data array
    :param twotheta: data array
    :param param_str: start string to label parameters
    :param comp_lists: list of lists for parameters
    :param conversion_factor: dict inputs for conversion
    :param bgguess: initial values or background dict
    :param Guess: initial values for peak components dict
    :param pfixed: fix p or not
    :return: list of lists with component coefficients
    """

    dfour, hfour, wfour, pfour = comp_lists

    comp = 'd'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    master_params = ff.unvary_part_params(master_params, param_str, comp, o['d-space'])
    # set part of these parameters to not vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, fit_method=None, weights=None, params=master_params)
    # master_params = out.params
    dfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'h'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    master_params = ff.unvary_part_params(master_params, param_str, comp, o['height'])
    # set part of these parameters to not vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, fit_method=None, weights=None, params=master_params)
    #master_params = out.params
    hfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'w'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    master_params = ff.unvary_part_params(master_params, param_str, comp, o['width'])
    # set part of these parameters to not vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, fit_method=None, weights=None, params=master_params)
    # master_params = out.params
    wfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'p'
    if pfixed is False:
        # if 'profile_fixed' not in order_peak[0]: # orders['peak'][0]:
        master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
        master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
        master_params = ff.unvary_part_params(master_params, param_str, comp, o['profile'])
        # set part of these parameters to not vary
        out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                          num_peaks=len(Guesses['peak']), nterms_back=len(bgguess),
                          Conv=conversion_factor,
                          fixed=pfixed, fit_method=None, weights=None, params=master_params)
        #master_params = out.params
    pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    return master_params, dfour, hfour, wfour, pfour


def FitSubpattern(TwoThetaAndDspacings, azimu, intens, new_data, orders=None, PreviousParams=None, bounds=None,
                  SaveFit=None, debug=False, refine=True, iterations=2, fnam=None, fit_method=None):
    """
    Perform the various fitting stages to the data
    :param new_data:
    :param TwoThetaAndDspacings:
    :param azimu:
    :param intens:
    :param orders:
    :param PreviousParams:
    :param bounds:
    :param SaveFit:
    :param debug:
    :param refine:
    :param iterations:
    :param fnam:
    :return:
    """
    
    # Measure the elapsed time during fitting.
    # To help decide what is bad fit or if over fitting the data.
    t_start = time.time()

    # Sort inputs
    # two theta and d-spacing are both 'x' dimensions. therefore import them as single array.
    # split them up here if needed. 
    if np.shape(TwoThetaAndDspacings)[0] == 3:
        twotheta = TwoThetaAndDspacings[0]
        dspace = TwoThetaAndDspacings[1]
        conversion_factor = TwoThetaAndDspacings[2]
    else:
        twotheta = TwoThetaAndDspacings
        # FIX ME: need to set a value for 'conversion_factor' in here and propagate its none-use through the code.

    if debug:  # get limits for chunk fit plots.
        tthrange = [np.min(twotheta.flatten()), np.max(twotheta.flatten())]

    # set data type for the intensity data.
    # This is needed for saving the fits using save_moderesult/ load_modelresult.
    # load_modelresult fails if the data is a masked integer array.
    if isinstance(intens, int) or isinstance(intens, np.int32):
        intens = float(intens)
    elif isinstance(intens, np.int64):
        intens = np.float64(intens)
    #FIX ME: this doesnt work. the data type is got by print(intens.dtype)

    # FIX ME: need to match orders of arrays to previous numbers.
    # if no parameters
    #     get orders of arrays
    # else if no orders
    #     use parameters as given
    # else
    #     change the size of the parameter arrays to match orders

    # define the number of peaks to be fitted
    peeks = len(orders['peak'])
    # DMF: FIX ME: orders used before checked?

    # generic parameters for the limits
    # if bounds is None:
    #     bounds = {
    #           "background": ['-inf', 'inf'],
    #           "d-space":    ['-inf', 'inf'],
    #           "height":     ['-inf', 'inf'],
    #           "profile":    ['-inf', 'inf'],
    #           "width":      ['-inf', 'inf'],
    #           }
    #     '''
    #     dmin = np.min(dspace)
    #     dmax = np.max(dspace)
    #     drange = dmax - dmin
    #     imin = np.min(intens)
    #     imax = np.max(intens)
    #     ndat = np.size(dspace)
    #     '''
        
    if orders:
        # print orders
        # total = orders['AziBins']
        backg_type = 'order'
        bg_order = orders['background']
        # print(bg_order, 'order')
        #        # FIX ME: assumes a fixed number of peaks.
        #        if 'symmetry' in orders['peak'][0]:
        #            symmetry = orders['peak'][0]['symmetry']
        #        else:
        #            symmetry = 1
        # print bg_order
        backg = []
        '''
        # removed because if using orders it should be list of orders i.e. number in list = num. orders needed
        for j in range(len(bg_order)):
            backg.append(np.zeros(len(bg_order) * 2 + 1))
        '''
        for j in bg_order:
            backg.append(np.zeros(np.max(j) * 2 + 1))
        
        for y in range(peeks):
            # force profile orders to match the fixed orders
            # FIX ME: should this generate an error?
            if 'profile_fixed' in orders['peak'][y]:
                if not ff.Fourier_order(orders['peak'][y]['profile_fixed']) == orders['peak'][y]['profile']:
                    print('Peak ' + str(y) + ': The order of the profile does not match that of the fixed profile. '
                                             'Changing profile to match fixed profile.')
                    orders['peak'][y]['profile'] = ff.Fourier_order(orders['peak'][y]['profile_fixed'])
                # make sure the fixed profile is a list
                if not isinstance(orders['peak'][y]['profile_fixed'], list):
                    orders['peak'][y]['profile_fixed'] = [orders['peak'][y]['profile_fixed']]

    if PreviousParams:
        #check if the previous fit was 'good' i.e. contrains no 'null' values.
        clean = AnyTermsNull(PreviousParams, val_to_find="Null")
        if clean == 0:
            #the previous fit has problems so discard it
            print('Propagated fit has problems so discarding it and doing fit from scratch')
            PreviousParams=None
            
            
    if PreviousParams and orders:
        
        choice_list = ['d-space', 'height', 'width', 'profile']
        
        # need to check sizes of both arrays here.
        # orders takes precedence over PreviousParams so we add or remove coefficients as necessary
        for y in range(peeks):
            # loop over parameters
            for Parm in choice_list:
                '''
                for x in range(4):
                # DMF: Why not just loop for Parm in orders['peak'][y]?
                # only 4 parameters describe each peak
                if x == 0:
                    Parm = 'd-space'
                elif x == 1:
                    Parm = 'height'
                elif x == 2:
                    Parm = 'width'
                elif x == 3:
                    Parm = 'profile'
                # elif x ==4:
                #    Parm = 'profile_fixed'
                '''
                coef_type = ff.params_get_type(PreviousParams, Parm, peak=y)
                if coef_type != 5: #if parameters are not independent 
                    if np.max(orders['peak'][y][Parm]) > ff.Fourier_order(PreviousParams['peak'][y][Parm]):
                        # print 'longer'
                        change_by = np.max(orders['peak'][y][Parm]) * 2 + 1 - np.size(PreviousParams['peak'][y][Parm])
                        # PreviousParams['peak'][y][Parm] = np.pad(PreviousParams['peak'][y][Parm], (0,change_by),
                        # mode='constant', constant_values=0)
                        PreviousParams['peak'][y][Parm] = (PreviousParams['background'][y] + [0] * change_by)
                    elif np.max(orders['peak'][y][Parm]) < ff.Fourier_order(PreviousParams['peak'][y][Parm]):
                        # print 'smaller'
                        change_by = np.size(PreviousParams['peak'][y][Parm]) - (np.max(orders['peak'][y][Parm]) * 2 + 1)
                        PreviousParams['peak'][y][Parm] = PreviousParams['peak'][y][Parm][1:-(change_by - 1)]
            if 'profile_fixed' in orders['peak'][y]:
                PreviousParams['peak'][y]['profile'] = orders['peak'][y]['profile_fixed']
            elif 'profile_fixed' in PreviousParams['peak'][y][Parm]:
                del PreviousParams['peak'][y]['profile_fixed']

        # loop for background orders/size
        if ff.params_get_type(PreviousParams, 'background') != 5: #if parameters are not independent 
            for y in range(np.max([len(orders['background']), len(PreviousParams['background'])])):
                # loop over the degree of the background
                if len(PreviousParams['background']) - 1 >= y and len(orders['background']) - 1 >= y:
                    # if present in both arrays make sure it is the right size.
                    if np.max(orders['background'][y]) > ff.Fourier_order(PreviousParams['background'][y]):
                        # print 'longer'
                        change_by = np.max(orders['background'][y]) * 2 + 1 - np.size(PreviousParams['background'][y])
                        # print change_by
                        # PreviousParams['background'][y] = np.pad(PreviousParams['background'][y], (0,change_by),
                        # mode='constant', constant_values=0)
                        PreviousParams['background'][y] = (PreviousParams['background'][y] + [0] * change_by)
                    elif np.max(orders['background'][y]) < ff.Fourier_order(PreviousParams['background'][y]):
                        # print 'smaller'
                        change_by = np.size(PreviousParams['background'][y]) - (np.max(orders['background'][y]) * 2 + 1)
                        PreviousParams['background'][y] = PreviousParams['background'][y][:-change_by]
                elif len(PreviousParams['background']) - 1 >= y > len(orders['background']) - 1:
                    # if previous params has too many degrees remove the higher ones
                    PreviousParams['background'][y] = []
                elif len(PreviousParams['background']) - 1 < y <= len(orders['background']) - 1:
                    # if previous parameters does not have enough degrees add more.
                    PreviousParams['background'].append([0] * (np.max(orders['background'][y]) * 2 + 1))

        # tidy up from reductions in length of PreviousParams
        PreviousParams['background'] = [x for x in PreviousParams['background'] if x != []]
        # FIX ME: NO backg_type set here!

    if PreviousParams and not orders:
        backg = PreviousParams['background']
        backg_type = 'coeffs'

    lenbg = []
    singlelenbg = []
    for val in backg:
        lenbg.append(len(val))
        singlelenbg.append(1)
    lenbg = np.array(lenbg)
    singlelenbg = np.array(singlelenbg)


    step = 1
    while step <= 10: #10 is arbitrarily large number and does not (currently) reflect the number of actual steps/stages in nthe process
        #use while loops so that if the fit is rubbish we can go back and improve it by repeatingn earlier steps. 
        # for chunks step <= 3 and for refine <= 6
        # we are using increments of 3 so that it is possible to utilise different states after the final fit.  
        
        if step <= 3:
            # generate chunks and initial fits
            # or parse previous fits into correct data structure
        
            # If there is not previous fit -- Fit data in azimuthal chunks for d.
            if not PreviousParams:
                if 'PeakPositionSelection' not in orders:
                    # DMF changed print statement to include number peaks in inputs file
                    if peeks > 1:
                        p = 'peaks'
                    else:
                        p = 'peak'
                    print('\nNo previous fits or selections. Assuming ' + str(peeks) + ' ' + p +
                          ' and group fitting in azimuth...\n')
                else:
                    print('\nUsing manual selections for initial guesses...\n')
        
                    # FIX ME: Need to check that the number of peaks for which we have parameters is the same as the number
                    # of peaks guessed at.
        
                    tthguesses = np.array(orders['PeakPositionSelection'])
                    # for future use in setting limits
                    dist = np.max(tthguesses[:, 2]) - np.min(tthguesses[:, 2])
                    width_change = dist / (2 ** peeks)
                    peak_vals = np.unique(tthguesses[:, 2])
        
                    # Fit Fourier series to two-theta/d-spacing for each peak
                    # FIX ME: this fitting should be done to the d-spacing not the tth/energy. Probably need to import the
                    # detector functions to make this happen.
                    dfour = []
                    for j in range(peeks):
                        peek = j + 1
                        tthguess = tthguesses[tthguesses[:, 0] == peek, 1:]
                        param_str = 'peak_' + str(j)
                        comp = 'd'
                        lims = ff.parse_bounds(bounds, twotheta.flatten(), 0, 0,  param=['d-space'])
                        
                        if 'd-space-type' in orders['peak'][j]:
                            coef_type = orders['peak'][j]['d-space-type']
                        else:
                            coef_type = 'fourier'
        
                        temp_param = Parameters()
                        temp_param = ff.initiate_params(param=temp_param, param_str=param_str, coef_type=coef_type, comp=comp,
                                                        trig_orders=orders['peak'][j]['d-space'], limits = lims['d-space'], value=tthguess[1, 1])
                        fout = ff.coefficient_fit(azimu=tthguess[:, 0], ydata=tthguess[:, 1],
                                              param=temp_param, param_str=param_str + '_' + comp, fit_method='leastsq')
                        temp_param = fout.params
                        if debug:
                            temp_param.pretty_print()
                        dfour.append(ff.gather_paramerrs_to_list(temp_param, param_str, comp))
                        '''
                        dfour.append(ff.Fourier_fit_old(tthguess[:, 0], det.Conversion(tthguess, conversion_factor),
                                                        terms=orders['peak'][j]['d-space']))
                        '''
                    
                # Get chunks according to detector type.
                chunks, azichunks = new_data.bins(azimu, orders)
        
                # Final output list of azimuths with corresponding twotheta_0,h,w
        
                # wavelength=parms_dict['wavelength']
                newd0 = [[] for _ in range(peeks)]
                newd0Err = [[] for _ in range(peeks)]
                newHall = [[] for _ in range(peeks)]
                newHallErr = [[] for _ in range(peeks)]
                newWall = [[] for _ in range(peeks)]
                newWallErr = [[] for _ in range(peeks)]
                newPall = [[] for _ in range(peeks)]
                newPallErr = [[] for _ in range(peeks)]
                newBGall = [[] for _ in range(len(bg_order))]
                newBGallErr = [[] for _ in range(len(bg_order))]
                newAziChunks = []
        
                for j in range(len(chunks)):
                    # print('\nFiting to data chunk ' + str(j + 1) + ' of ' + str(len(chunks)) + '\n')
                    if debug:
                        print('\n')
                    print('\rFitting to data chunk %s of %s ' % (str(j + 1), str(len(chunks))), end='\r')
                    # print '\nFitting azimuth chunk....\n'
                    # only compute the fit and save it if there are 'enough' data
                    min_dat = 21  # minimum data in each chunk
                    n = 5  # how many data to use in constructing guesses.
                    wguess_fraction = 10.
                    # N.B. if min_dat/n is too large the guesses will become very strange.
                    # FIX ME: These should be input variables (for use in extremes)
                    # FIX ME: possibly do this better using scipy.signal.find or scipy.signal.find_peaks_cwt
        
                    if np.sum(~intens.flatten()[chunks[j]].mask) >= min_dat:
                        # Get indices of sorted two theta values excluding the masked values
                        tth_ord = ma.argsort(twotheta.flatten()[chunks[j]].compressed())
                        # Background estimates
                        backg_guess = [[0.] for i in range(len(orders['background']))]
        
                        if backg_type == 'coeffs':
                            # Fourier expand the coefficients.
                            # FIX ME!! Altered this loop from using 'j' as already outer loop, needs checking below!
                            for k in range(len(orders['background'])):
                                backg_guess[j] = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]),
                                                                   param=orders['peak']['background'][k])
                        elif backg_type == 'order':
                            # first value (offset) is mean of left hand values
                            backg_guess[0][0] = np.mean(intens.flatten()[chunks[j]].compressed()[tth_ord[:n]])
                            # if there are more then calculate a gradient guess. 
                            # leave all higher order terms as 0 -- assume they are small.
                            if len(orders['background']) > 1:
                                backg_guess[1][0] = (np.mean(intens.flatten()[chunks[j]].compressed()[tth_ord[-n:]]) - np.mean(
                                    intens.flatten()[chunks[j]].compressed()[tth_ord[:n]])) / (
                                                            twotheta.flatten()[chunks[j]].compressed()[tth_ord[-1]] -
                                                            twotheta.flatten()[chunks[j]].compressed()[tth_ord[0]])
                        # elif backg_type == 'flat':
                        #     backg_guess = backg[0][0]
                        # Organise guesses to be refined.
                        dguess = []
                        hguess = []
                        wguess = []
                        pguess = []
                        pfixed = 0
                        peaks = []
                        lims = []
        
                        for k in range(peeks):
                            if 'dfour' in locals():  # If the positions have been pre-guessed extract values from dfour.
                                if debug and k == 0:
                                    print('dfour in locals')
                                # Guess returns tthguess from dfour then converts in into d-spacing.
                                if 'd-space-type' in orders['peak'][k]:
                                    coef_type = orders['peak'][k]['d-space-type']
                                else:
                                    coef_type = 'fourier'
                                tthguess = ff.coefficient_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0], coef_type=coef_type)
                                #tthguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0])
                                dguess = new_data.conversion(tthguess, conversion_factor,
                                                             azm=np.mean(azimu.flatten()[chunks[j]]))
                                # Finds the index of the closest two-theta to the d-spacing input
                                idx = ((np.abs(twotheta.flatten()[chunks[j]] - tthguess)).argmin())
                                # FIX ME: The mean of a number of the smallest values would be more stable.
        
                                # Height is intensity of closest pixel in d-spacing - background at that position
                                hguess = intens.flatten()[chunks[j]][idx]
                                if len(backg_guess) > 1:
                                    hguess = hguess - (backg_guess[0][0] + backg_guess[1][0] * (
                                            twotheta.flatten()[chunks[j]][idx] - twotheta.flatten()[chunks[j]].min()))
                                elif len(backg_guess) == 1:
                                    hguess = hguess - backg_guess[0][0]
                                    
                                # # wguess is fractional width of the data range
                                # wguess = (det.Conversion(np.max(dspace.flatten()[chunks[j]]), conversion_factor, reverse=1,
                                # azm=np.mean(azimu.flatten()[chunks[j]])) -
                                #           det.Conversion(np.min(dspace.flatten()[chunks[j]]), conversion_factor, reverse=1,
                                #           azm=np.mean(azimu.flatten()[chunks[j]])) ) / wguess_fraction
                                # # FIX ME: This is very crude. Need a better way to do it. \
        
                                # pguess = 0.5  # guess half if solving for profile
                                # pfixed = 0  # set to 0 so solving unless profile_fixed exists.
                                # if 'profile_fixed' in orders['peak'][k]:
                                #     # FIX ME: profile fixed is the switch 1/0/true/false and profile is the profile
                                #     parameters. Make sure this is the case in the input files
                                #     pguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]),
                                #     param=orders['peak'][k]['profile_fixed'])
                                #     pfixed = 1
        
                                # peaks.append({"d-space": [dguess],
                                #               "height": [hguess],
                                #               "profile": [pguess],
                                #               "width": [wguess]
                                #               })
            
                            else: # guess guesses from data slice - assuming a single peak
                                
                                # find brightest pixels
                                
                                # get index of nth brightest pixel.
                                idx = np.argsort(intens.flatten()[chunks[j]].compressed())[-n:][0]
                                
                                # height guess is nth highest intensity - backgound guess at this position
                                hguess = (intens.flatten()[chunks[j]].compressed())[idx]                        
                                if len(backg_guess) > 1:
                                    hguess = hguess - (backg_guess[0][0] + backg_guess[1][0] *
                                                       (twotheta.flatten()[chunks[j]].compressed()[idx] -
                                                        twotheta.flatten()[chunks[j]].min() ))
                                elif len(backg_guess) == 1:
                                    hguess = hguess - backg_guess[0][0]
                                    
                                # d-spacing of highest nth intensity pixel
                                dguess = dspace.flatten()[chunks[j]].compressed()[idx]
        
                            # wguess is fractional width of the data range
                            wguess = (np.max(twotheta.flatten()[chunks[j]]) -
                                      np.min(twotheta.flatten()[chunks[j]])) / wguess_fraction / peeks
                            # FIX ME: This is a bit crude. Is there a better way to do it?
        
                            # If profile_fixed exists then set pguess to profile_fixed values and set pfixed to 1
                            pguess = 0.5  # Guess half if solving for profile
                            pfixed = 0  # Set to 0 so solving unless profile_fixed exists.
                            if 'profile_fixed' in orders['peak'][k]:
                                # Profile fixed is the list of coefficients if the profile is fixed.
                                # FIX ME: profile fixed must have the same number of coefficients as required by profile.
                                pguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]) * orders['peak'][k]['symmetry'], param=orders['peak'][k]['profile_fixed'])
                                pfixed = 1
            
                            peaks.append({"d-space": [dguess], "height": [hguess], "width": [wguess], "profile": [pguess]})
                          
                            # DMF added needs checking - required to drive fit and avoid failures
                            lims.append(ff.parse_bounds(bounds, dspace.flatten()[chunks[j]], intens.flatten()[chunks[j]],
                                                        twotheta.flatten()[chunks[j]],
                                                        param=['d-space', 'height', 'width', 'profile']))
                            # lims.append({"d-space": [np.min(dspace.flatten()[chunks[j]]), np.max(dspace.flatten()[chunks[j]])]
                            # ,
                            #              "height": [0, np.max(intens.flatten()[chunks[j]])],
                            #              "width": [(np.max(twotheta.flatten()[chunks[j]]) -
                            #              np.min(twotheta.flatten()[chunks[j]])) / 100 / peeks,
                            #                        (np.max(twotheta.flatten()[chunks[j]]) -
                            #                        np.min(twotheta.flatten()[chunks[j]])) / 2 / peeks],
                            #              "profile": [0, 1]
                            #              })
                            # lims.append({"d-space": bounds["d-space"],
                            #              "height":  bounds["height"],
                            #              "width":   bounds["width"],
                            #              "profile": bounds["profile"],
                            #              })
        
                        # # Guess for background, if orders and therefore no guess input, choose appropriate
                        # singleBackg = [[] for i in range(len(backg))]
                        # for i in range(len(backg)):
                        #     singleBackg[i] = [backg[i][0]]
                        # bgguess = singleBackg
                        # bgguess[0][0] = 5 
                        # print(bgguess)
                        
                        # FIX ME: SAH, while making refine independent of the chunk fitting:: Guesses is a superfluous variable removing it
                        #Guesses = {"background": backg_guess, #bgguess,
                        #           "peak": peaks}
                        limits = {"background": ff.parse_bounds(bounds, dspace.flatten()[chunks[j]],
                                                                intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]],
                                                                param=['background'])['background'], "peak": lims}
                        
                        # Define parameters to pass to fit
                        params = Parameters()
                        # if orders: bg of form [num_orders_0,num_orders_1]
                        # if coeffs: bg of form [[coeff1],[coeff2, coeff3, coeff4]]
        
                        # Chunks limited to no Fourier expansion so only single value per polynomial order.
                        # for b in range(len(bgguess)):
                        #     params.add('bg_c' + str(b) + '_f' + str(0), bgguess[b][0])
                        for b in range(len(backg_guess)):
                            if b == 0:
                                params.add('bg_c' + str(b) + '_f' + str(0), backg_guess[b][0], min=limits['background'][0],
                                           max=limits['background'][1])
                            else:
                                params.add('bg_c' + str(b) + '_f' + str(0), backg_guess[b][0])
        
                        for pk in range(len(peaks)):
        
                            if 'dfour' in locals():
                                params.add('peak_' + str(pk) + '_d0', peaks[pk]['d-space'][0],
                                           min=limits['peak'][pk]['d-space'][0], max=limits['peak'][pk]['d-space'][1])
                                params.add('peak_' + str(pk) + '_h0', peaks[pk]['height'][0],
                                           min=limits['peak'][pk]['height'][0], max=limits['peak'][pk]['height'][1])
                                params.add('peak_' + str(pk) + '_w0', peaks[pk]['width'][0],
                                           min=limits['peak'][pk]['width'][0], max=limits['peak'][pk]['width'][1])
                                if 'profile_fixed' in orders['peak'][pk]:
                                    params.add('peak_' + str(pk) + '_p0', peaks[pk]['profile'][0],
                                               min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1],
                                               vary=False)
                                else:
                                    params.add('peak_' + str(pk) + '_p0', peaks[pk]['profile'][0],
                                               min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1])
        
        
                            else:
                                # If using bounds/limits should be written as below:
                                params.add('peak_' + str(pk) + '_d0', peaks[pk]['d-space'][0],
                                           min=limits['peak'][pk]['d-space'][0], max=limits['peak'][pk]['d-space'][1])
                                params.add('peak_' + str(pk) + '_h0', peaks[pk]['height'][0],
                                           min=limits['peak'][pk]['height'][0], max=limits['peak'][pk]['height'][1])
                                params.add('peak_' + str(pk) + '_w0', peaks[pk]['width'][0],
                                           min=limits['peak'][pk]['width'][0], max=limits['peak'][pk]['width'][1])
                                if 'profile_fixed' in orders['peak'][pk]:
                                    params.add('peak_' + str(pk) + '_p0', peaks[pk]['profile'][0],
                                               min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1],
                                               vary=False)
                                else:
                                    params.add('peak_' + str(pk) + '_p0', peaks[pk]['profile'][0],
                                               min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1])
        
                        if debug:
                            params.pretty_print()
                            print('\n')
                            guess = params
                            
                        out = ff.FitModel(intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]], azichunks[j],
                                          num_peaks=len(peaks), nterms_back=len(backg_guess), Conv=conversion_factor,
                                          fixed=pfixed, fit_method=fit_method, weights=None, params=params)
                        params = out.params
                        if debug:
                            params.pretty_print()
        
                        for i in range(peeks):
                            newd0[i].append(params['peak_' + str(i) + '_d0'].value)
                            newd0Err[i].append(params['peak_' + str(i) + '_d0'].stderr)
                            newHall[i].append(params['peak_' + str(i) + '_h0'].value)
                            newHallErr[i].append(params['peak_' + str(i) + '_h0'].stderr)
                            newWall[i].append(params['peak_' + str(i) + '_w0'].value)
                            newWallErr[i].append(params['peak_' + str(i) + '_w0'].stderr)
                            # profile should now be bound by the settings in the parameter
                            newPall[i].append(params['peak_' + str(i) + '_p0'].value)
                            newPallErr[i].append(params['peak_' + str(i) + '_p0'].stderr)
                            # may have to set profile error when initiate params so have appropriate error value if fixed
                        for i in range(len(backg_guess)):
                            newBGall[i].append(params['bg_c' + str(i) + '_f0'].value)
                            newBGallErr[i].append(params['bg_c' + str(i) + '_f0'].stderr)
        
                        newAziChunks.append(azichunks[j])
        
                        # plot the fits.
                        if debug:
                            tth_plot = twotheta.flatten()[chunks[j]]
                            int_plot = intens.flatten()[chunks[j]]
                            azm_plot = np.tile(azimu.flatten()[chunks[j]][0], (300))
                            azm_plot = ma.array(azm_plot, mask=(~np.isfinite(azm_plot)))
                            # FIX ME: required for plotting energy dispersive data.
                            gmodel = Model(ff.PeaksModel, independent_vars=['twotheta', 'azi'])
                            tth_range = np.linspace(np.min(twotheta.flatten()[chunks[j]]),
                                                    np.max(twotheta.flatten()[chunks[j]]), azm_plot.size)
                            # mod_plot = gmodel.eval(params=params, twotheta=tth_range, azi=azm_plot,
                            #                        num_peaks=len(Guesses['peak']),
                            #                        nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed)
                            # guess_plot = gmodel.eval(params=guess, twotheta=tth_range, azi=azm_plot,
                            #                        num_peaks=len(Guesses['peak']),
                            #                        nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed)
                            mod_plot = gmodel.eval(params=params, twotheta=tth_range, azi=azm_plot, Conv=conversion_factor)
                            guess_plot = gmodel.eval(params=guess, twotheta=tth_range, azi=azm_plot, Conv=conversion_factor)
                            plt.plot(tth_plot, int_plot, '.', label="data")
                            plt.plot(tth_range, guess_plot, marker='', color='green', linewidth=2, linestyle='dashed',
                                     label='guess')
                            plt.plot(tth_range, mod_plot, marker='', color='red', linewidth=2, label='fit')
                            plt.xlim(tthrange)
                            plt.legend()
                            plt.title(peak_string(orders) + '; Chunk ' + str(j + 1))
                            plt.show()
        
                # Feed each d_0,h,w into Fourier expansion function to get fit for fourier component
                # parameters as output.
        
                print('\n Performing Series fits...')
        
                dfour = []
                hfour = []
                wfour = []
                pfour = []
                master_params = Parameters()  # Initiate total Parameter class to add to
                        
                lims = ff.parse_bounds(bounds, dspace.flatten(), intens.flatten(), twotheta.flatten())
                # replace fixed array with dynamic limits
        
                for j in range(peeks):
        
                    param_str = 'peak_' + str(j)  # defines peak string to start parameter name
        
                    # initiate symmetry
                    comp = 's'
                    if 'symmetry' in orders['peak'][j].keys():
                        symm = orders['peak'][j]['symmetry']
                    else:
                        symm = 1
                    master_params = ff.initiate_params(master_params, param_str, comp, 0, limits=[10, 1], value=np.array(symm), vary=False)

                    # initiate d
                    comp = 'd'
                    coef_type = ff.params_get_type(orders, comp, peak=j)
                    n_coef    = ff.params_get_number_coef(orders, comp, peak=j, azims=np.array(newAziChunks))
                    master_params = ff.initiate_params(master_params, param_str, comp, coef_type=coef_type, num_coef=n_coef, trig_orders=orders['peak'][j]['d-space'], limits=lims['d-space'], value=np.array(newd0[j])) 
                    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                    if isinstance(orders['peak'][k]['d-space'], list): # set part of these parameters to not vary
                        master_params = ff.unvary_part_params(master_params, param_str, comp, orders['peak'][k]['d-space'])
                    fout = ff.coefficient_fit(azimu=np.array(newAziChunks), ydata=np.array(newd0[j]), param=master_params, param_str=param_str + '_' + comp, symm=1,  errs=np.array(newd0Err[j]), fit_method='leastsq')
                    master_params = fout.params
                    # FIX ME. Check which params should be varying and which should not. Need to incorporate var y and unnvary params as well as partial vary
                    
                    
                    # initiate h
                    comp = 'h'
                    coef_type = ff.params_get_type(orders, comp, peak=j)
                    n_coef    = ff.params_get_number_coef(orders, comp, peak=j, azims=np.array(newAziChunks))
                    master_params = ff.initiate_params(master_params, param_str, comp, coef_type=coef_type, num_coef=n_coef, trig_orders=orders['peak'][j]['height'], limits=lims['height'], value=np.array(newHall[j]))
                    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                    if isinstance(orders['peak'][k]['height'], list): # set part of these parameters to not vary
                        master_params = ff.unvary_part_params(master_params, param_str, comp, orders['peak'][k]['height'])
                    fout = ff.coefficient_fit(azimu=np.array(newAziChunks), ydata=np.array(newHall[j]), param=master_params, param_str=param_str + '_' + comp, symm=symm,  errs=np.array(newHallErr[j]), fit_method='leastsq')
                    master_params = fout.params
                    
                    # initiate w
                    comp = 'w'
                    coef_type = ff.params_get_type(orders, comp, peak=j)
                    n_coef    = ff.params_get_number_coef(orders, comp, peak=j, azims=np.array(newAziChunks))
                    master_params = ff.initiate_params(master_params, param_str, comp, coef_type=coef_type, num_coef=n_coef, trig_orders=orders['peak'][j]['width'], limits=lims['width'], value=np.array(newWall[j]))
                    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                    if isinstance(orders['peak'][k]['width'], list): # set part of these parameters to not vary
                        master_params = ff.unvary_part_params(master_params, param_str, comp, orders['peak'][k]['width'])
                    fout = ff.coefficient_fit(azimu=np.array(newAziChunks), ydata=np.array(newWall[j]), param=master_params, param_str=param_str + '_' + comp, symm=symm,  errs=np.array(newWallErr[j]), fit_method='leastsq')
                    master_params = fout.params
                    
                    # initiate p
                    comp = 'p'
                    coef_type = ff.params_get_type(orders, comp, peak=j)
                    n_coef    = ff.params_get_number_coef(orders, comp, peak=j, azims=np.array(newAziChunks))
                    if 'profile_fixed' in orders['peak'][j]:
                        p_vals = orders['peak'][j]['profile_fixed']
                        pfixed = 1
                    else:
                        p_vals = None
                        pfixed = 0
                    master_params = ff.initiate_params(master_params, param_str, comp, coef_type=coef_type, num_coef=n_coef, trig_orders=orders['peak'][j]['profile'], limits=lims['profile'], value=p_vals)
                    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                    if isinstance(orders['peak'][k]['profile'], list): # set part of these parameters to not vary
                        master_params = ff.unvary_part_params(master_params, param_str, comp, orders['peak'][k]['profile'])
                    if pfixed == 0:
                        fout = ff.coefficient_fit(azimu=np.array(newAziChunks), ydata=np.array(newPall[j]), param=master_params, param_str=param_str + '_' + comp, symm=symm,  errs=np.array(newPallErr[j]), fit_method='leastsq')
                        master_params = fout.params
                    
                # Initiate background parameters and perform an initial fit
                comp = 'bg'
                coef_type = ff.params_get_type(orders, comp, peak=j)
                for b in range(len(backg)):
                    param_str= 'bg_c'+str(b)
                    comp = 'f'
                    if b==0:
                        limits = lims['background']
                    else:
                        limits = None
                    n_coef    = ff.params_get_number_coef(orders, comp, peak=b, azims=np.array(newAziChunks))
                    master_params = ff.initiate_params(master_params, param_str, comp, coef_type=coef_type, num_coef=n_coef, trig_orders=orders['background'][b], limits=limits)
                    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                    if isinstance(orders['background'][b], list): # set part of these parameters to not vary
                        master_params = ff.unvary_part_params(master_params, param_str, comp, orders['background'][b])
                    fout = ff.coefficient_fit(azimu=np.array(newAziChunks), ydata=np.array(newBGall[b]), param=master_params, param_str=param_str + '_' + comp, symm=symm,  errs=np.array(newBGallErr[b]), fit_method='leastsq')
                    master_params = fout.params
                    
                
                if debug:
                    print('Parameters after initial Fourier fits')
                    master_params.pretty_print()
                
                # plot output of fourier fits....
                if debug:
                    y_lims = np.array([np.min(newAziChunks), np.max(newAziChunks)])
                    y_lims = np.around(y_lims / 180) * 180
                    AziPlot = range(np.int(y_lims[0]), np.int(y_lims[1]), 2)
                    #AziPlot = range(0, 360, 2)
                    gmodel = Model(ff.coefficient_expand, independent_vars=['azimu'])
        
                    fig = plt.figure()
                    ax1 = fig.add_subplot(5, 1, 1)
                    ax1.set_title('D-spacing')
                    for j in range(peeks):
                        param_str = 'peak_'+str(j)
                        comp='d'
                        temp = ff.gather_paramerrs_to_list(master_params, param_str, comp)
                        temp_tp = ff.get_series_type(master_params, param_str, comp)
                        if temp_tp==5:
                            AzPlt = np.array(newAziChunks)
                        else:
                            AzPlt=AziPlot
                        gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AzPlt), param=temp[0], coef_type=temp_tp)
                        ax1.scatter(newAziChunks, newd0[j], s=10)
                        ax1.plot(AzPlt, gmod_plot)
        
                    # plt.subplot(512)
                    ax2 = fig.add_subplot(5, 1, 2)
                    ax2.set_title('height')
                    for j in range(peeks):
                        if 'symmetry' in orders['peak'][j].keys():
                            symm = orders['peak'][j]['symmetry']
                        else:
                            symm = 1
                        param_str = 'peak_'+str(j)
                        comp='h'
                        temp = ff.gather_paramerrs_to_list(master_params, param_str, comp)
                        temp_tp = ff.get_series_type(master_params, param_str, comp)
                        if temp_tp==5:
                            AzPlt = np.array(newAziChunks)
                        else:
                            AzPlt=AziPlot
                        gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AzPlt)*symm, param=temp[0], coef_type=temp_tp)
                        ax2.scatter(newAziChunks, newHall[j], s=10)
                        ax2.plot(AzPlt, gmod_plot)
        
                    ax3 = fig.add_subplot(5, 1, 3)
                    # plt.subplot(513)
                    ax3.set_title('width')
                    for j in range(peeks):
                        if 'symmetry' in orders['peak'][j].keys():
                            symm = orders['peak'][j]['symmetry']
                        else:
                            symm = 1
                        param_str = 'peak_'+str(j)
                        comp='w'
                        temp = ff.gather_paramerrs_to_list(master_params, param_str, comp)
                        temp_tp = ff.get_series_type(master_params, param_str, comp)
                        if temp_tp==5:
                            AzPlt = np.array(newAziChunks)
                        else:
                            AzPlt=AziPlot
                        gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AzPlt)*symm, param=temp[0], coef_type=temp_tp)
                        ax3.scatter(newAziChunks, newWall[j], s=10)
                        ax3.plot(AzPlt, gmod_plot)
        
                    ax4 = fig.add_subplot(5, 1, 4)
                    # plt.subplot(514)
                    ax4.set_title('Profile')
                    for j in range(peeks):
                        if 'symmetry' in orders['peak'][j].keys():
                            symm = orders['peak'][j]['symmetry']
                        else:
                            symm = 1
                        param_str = 'peak_'+str(j)
                        comp='p'
                        temp = ff.gather_paramerrs_to_list(master_params, param_str, comp)
                        temp_tp = ff.get_series_type(master_params, param_str, comp)
                        if temp_tp==5:
                            AzPlt = np.array(newAziChunks)
                        else:
                            AzPlt=AziPlot
                        gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AzPlt)*symm, param=temp[0], coef_type=temp_tp)
                        ax4.scatter(newAziChunks, newPall[j], s=10)
                        ax4.plot(AzPlt, gmod_plot)
        
                    for k in range(len(lenbg)):
                        x_plt = len(lenbg)
                        y_plt = len(lenbg) * 4 + k + 1
                        ax5 = fig.add_subplot(5, x_plt, y_plt)
                        ax5.set_title('Background')
                        
                        param_str = 'bg_c'+str(k)
                        comp='f'
                        temp = ff.gather_paramerrs_to_list(master_params, param_str, comp)
                        temp_tp = ff.get_series_type(master_params, param_str, comp)
                        if temp_tp==5:
                            AzPlt = np.array(newAziChunks)
                        else:
                            AzPlt=AziPlot
                        gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AzPlt), param=temp[0], coef_type=temp_tp)
                        ax5.scatter(newAziChunks, newBGall[k], s=10)
                        ax5.plot(AzPlt, gmod_plot)
        
                    fig.suptitle(peak_string(orders) + '; Fits to Chunks')
                    plt.show()
                    plt.close()
        
                # put fits into NewParams data structure.
                #NewParams = ff.create_newparams(peeks, dfour, hfour, wfour, pfour, bgfour, orders['peak'])
        
            else:
                # FIX ME: This should load the saved lmfit Parameter class object. Need to initiate usage (save and load).
                # For now have propagated use of NewParams so re-instantiate the master_params object below, which is clunky.
                print('Using previously fitted parameters and propagating fit')
                NewParams = PreviousParams
        
                # FIX ME: should this be some def or subfunction?
                master_params = Parameters()  # initiate total Parameter class to add to
                # FIX ME: parameters initiated with no limits currently 
        
                lims = ff.parse_bounds(bounds, dspace.flatten(), intens.flatten(), twotheta.flatten())
                # replace fixed array with dynamic limits
        
                for j in range(peeks):
        
                    param_str = 'peak_' + str(j)  # defines peak string to start parameter name
        
                    # initiate symmetry
                    comp = 's'
                    if 'symmetry' in orders['peak'][j].keys():
                        symm = orders['peak'][j]['symmetry']
                    else:
                        symm = 1
                    master_params = ff.initiate_params(master_params, param_str, comp, 0, limits=[10, 1], value=np.array(symm), vary=False)
        
                    # initiate d
                    comp = 'd'
                    master_params = ff.initiate_params(master_params, param_str, comp, #orders['peak'][j]['d-space'],
                                                       coef_type=PreviousParams['peak'][j]['d-space-type'],
                                                       value=PreviousParams['peak'][j]['d-space'], 
                                                       limits=lims['d-space'] )
        
                    # initiate h
                    comp = 'h'
                    master_params = ff.initiate_params(master_params, param_str, comp,
                                                       coef_type=PreviousParams['peak'][j]['height-type'],
                                                       value=PreviousParams['peak'][j]['height'],
                                                       limits=lims['height'] )
        
                    # initiate w
                    comp = 'w'
                    master_params = ff.initiate_params(master_params, param_str, comp, 
                                                       coef_type=PreviousParams['peak'][j]['width-type'],
                                                       value=PreviousParams['peak'][j]['width'],
                                                       limits=lims['width'] )
        
                    # initiate p
                    comp = 'p'
                    if orders and 'profile_fixed' in orders['peak'][j]:
                        pvals = orders['peak'][j]['profile_fixed']
                        pfixed = 1
                    else:
                        pvals = PreviousParams['peak'][j]['profile']
                        pfixed = 0
                    master_params = ff.initiate_params(master_params, param_str, comp,
                                                       coef_type=PreviousParams['peak'][j]['profile-type'],
                                                       value=pvals, 
                                                       limits=lims['profile'] )
        
                # Initiate background parameters
                comp = 'bg'
                for b in range(len(PreviousParams['background'])):
                    param_str= 'bg_c'+str(b)
                    comp = 'f'
                    if b==0:
                        limits = lims['background']
                    else:
                        limits = None
                    master_params = ff.initiate_params(master_params, param_str, comp,
                                                       value=PreviousParams['background'][b],
                                                       coef_type=PreviousParams['background-type'],
                                                       limits=limits)
                    
                    '''
                    
                    
                    
                    
                    
                    bg_four = PreviousParams['background'][b]
                    for f_b in range(len(bg_four)):
                        if b is 0 and f_b is 0:
                            master_params.add('bg_c' + str(b) + '_f' + str(f_b), bg_four[f_b], min=lims['background'][0],
                                              max=lims['background'][1])
                        else:
                            master_params.add('bg_c' + str(b) + '_f' + str(f_b), value=bg_four[f_b])
                    '''
                # master_params.pretty_print()
        
                # Guess for background, if orders and therefore no guess input, choose appropriate
                # copied from above to allow the code to run. FIX ME: doesn't need all this code.
                singleBackg = [[] for i in range(len(backg))]
                for i in range(len(backg)):
                    singleBackg[i] = [backg[i][0]]
                backg_guess = singleBackg
                # bgguess[0][0] = 5
        
                # def previous_params_to_params(PreviousParams):
        
                # FIX ME: need to confirm the number of parameters matches the orders of the fits.
                
            step = step + 3
            
            
        if step >= 4:
            
            if refine or step>=5 or not PreviousParams:
                # Iterate over each parameter series in turn.
                print('\nRe-fitting for d, h, w, bg separately... will refine %i time(s)\n' % iterations)
                for j in range(iterations):
                    for k in range(peeks):
                        param_str = 'peak_'+str(k)
                        comp_list = ['d', 'h', 'w'] # always need to iterate over d, height and width
                        comp_names = ['d-space','height','width','profile']
                        if 'profile_fixed' in orders['peak'][k]:
                            pfixed=1
                        else:
                            pfixed=0
                            comp_list.append('p') # if the profile is not fixed iterate over this as well.
                        for cp in range(len(comp_list)):
                            comp = comp_list[cp]
                            master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                            master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                            if isinstance(orders['peak'][k][comp_names[cp]], list): # set part of these parameters to not vary
                                master_params = ff.unvary_part_params(master_params, param_str, comp, orders['peak'][k][comp_names[cp]])
                            fout = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                                          num_peaks=peeks, nterms_back=len(backg), Conv=conversion_factor,
                                          fixed=pfixed, fit_method=None, weights=None, params=master_params)
                            master_params = fout.params
                            
                    for k in range(len(backg)):
                        param_str = 'bg_c'+str(k)
                        comp = 'f'
                        master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
                        master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                        #master_params = ff.unvary_part_params(master_params, param_str, comp, orders['background'][k])
                        # set part of these parameters to not vary
                        fout = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                                      num_peaks=peeks, nterms_back=len(backg), Conv=conversion_factor,
                                      fit_method=None, weights=None, params=master_params)
                        master_params = fout.params
                        
                    if debug:
                        print(' ')
                        print('Parameters after refining seires fits ' + str(j+1) + ' time(s)')
                        master_params.pretty_print()
                
                
            step = step + 3
        
        if step >= 7:
            # FIX ME: Need to sort out use of lmfit parameters and NewParams, below will not work currently for loaded
            # parameters as master_params not created from NewParams yet. If possible use the json import/export within lmfit.
        
            # Refit through full equation with all data for d,h,w,bg independently
            print('\nFinal fit solving for all parms...\n')
        
            if step == 7:
                max_nfev = 500
            elif step == 9:
                max_nfev = None
            else: 
                max_nfev = 1000
                
        
            # set all parameters to vary
            for k in range(peeks):
                param_str = 'peak_'+str(k)
                comp_list = ['d', 'h', 'w'] # always need to iterate over d, height and width
                comp_names = ['d-space','height','width','profile']
                if 'profile_fixed' in orders['peak'][k]:
                    pfixed=1
                else:
                    pfixed=0
                    comp_list.append('p') # if the profile is not fixed iterate over this as well.
                for cp in range(len(comp_list)):
                    master_params = ff.vary_params(master_params, param_str, comp_list[cp])
                    if isinstance(orders['peak'][k][comp_names[cp]], list): # set part of these parameters to not vary
                        master_params = ff.unvary_part_params(master_params, param_str, comp_list[cp], orders['peak'][k][comp_names[cp]])
            for k in range(len(backg)):
                param_str = 'bg_c'+str(k)
                comp = 'f'
                master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
                #master_params = ff.unvary_part_params(master_params, param_str, comp, orders['background'][k])
                # set part of these parameters to not vary
                        
            
            '''
            for p in master_params:
                master_params[p].set(vary=True)
            for x in range(peeks):
                # unvary symmetry!! This should never vary.
                new_str = 'peak_' + str(x) + '_s'
                str_keys = [key for key, val in master_params.items() if new_str in key]
                for pstr in str_keys:
                    master_params[pstr].set(vary=False)
        
                # unvary p for appropriate peaks
                if 'profile_fixed' in orders['peak'][x]:
                    new_str = 'peak_' + str(x) + '_p'
                    str_keys = [key for key, val in master_params.items() if new_str in key]
                    for pstr in str_keys:
                        master_params[pstr].set(vary=False)
                        
                # unvary specific parts of components if needed
                master_params = ff.unvary_part_params(master_params, 'peak_' + str(x), 'd', orders['peak'][x]['d-space']) 
                master_params = ff.unvary_part_params(master_params, 'peak_' + str(x), 'h', orders['peak'][x]['height']) 
                master_params = ff.unvary_part_params(master_params, 'peak_' + str(x), 'w', orders['peak'][x]['width']) 
                master_params = ff.unvary_part_params(master_params, 'peak_' + str(x), 'p', orders['peak'][x]['profile']) 
                
            for x in range(len(orders['background'])):
                master_params = ff.unvary_part_params(master_params, param_str+str(x), 'f', orders['background'][x])
            '''
            #master_params.pretty_print()
            
            # FIX ME: need to sort parameters, bgguess doesn't exist if load previous fit data
            out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), num_peaks=peeks,
                              nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed, fit_method=None,
                              weights=None, params=master_params, max_nfev=max_nfev)
            master_params = out.params
    
            if out.success == 1:
                # it worked, carry on
                step = step + 10
                master_params = out.params
            elif out.success == 0 and refine==False and PreviousParams:
                #refine the fourier series before tying again
                step = 5
                iterations = np.max((iterations, 3))
                # FIX ME: could add another option for fixing the profile if this doesn't have the desired effect
            elif out.success == 0 and refine==True and PreviousParams:
                #if we refined the previous params and it didnt work then fit the chunks again.
                step = 2
                iterations = np.max((iterations, 3))
            else:
                step = step + 1
                
            #FIX ME: we could make an option for output the chunks without any Fourier/global fitting. Would this be useful?
                
    print('\nFinal Coefficients\n')
    # print(out.fit_report(min_correl=1))
    print(out.fit_report(show_correl=False))
    
    # Print some stats
    print('number data', out.ndata)
    print('degrees of freedom', out.nfree)
    print('ChiSquared', out.chisqr)
    
    # Write master_params to NewParams dict object
    NewParams = ff.params_to_newparams(master_params, num_peaks=peeks, num_bg=len(backg_guess),
                                       order_peak=orders['peak'])

    tth_min = twotheta.min()
    tth_max = twotheta.max()
    d_min = dspace.min()
    d_max = dspace.max()
    extent = [[tth_min, tth_max], [d_min, d_max]]  
    # FIX ME: This is taking the maximum and the minimum of the data not the 'range' itself.
    NewParams.update({'range': extent})

    # Elapsed time for fitting
    t_end = time.time()
    t_elapsed = t_end - t_start
    # get the rest of the fit stats from the lmfit output.
    FitStats = {"time-elapsed": t_elapsed,
                "status": step,
                "sum-residuals-squared": np.sum(out.residual**2),
                "function-evaluations": out.nfev,
                "n-variables": out.nvarys,
                "n-data": out.ndata,
                "degree-of-freedom": out.nfree, 
                "ChiSq": out.chisqr, 
                "RedChiSq": out.redchi, 
                "aic": out.aic, 
                "bic": out.bic,
                }
    NewParams.update({'FitProperties': FitStats})


    # Plot results to check
    view = 0
    if SaveFit == 1 or view == 1 or debug:
        print('\nPlotting results for fit...\n')

        gmodel = Model(ff.PeaksModel, independent_vars=['twotheta', 'azi'])
        # fullfit_intens = gmodel.eval(params=master_params, twotheta=twotheta.flatten(), azi=azimu.flatten(),
        #                              num_peaks=num_peaks, nterms_back=len(backg_guess), Conv=conversion_factor,
        #                              fixed=pfixed)
        fullfit_intens = gmodel.eval(params=master_params, twotheta=twotheta.flatten(), azi=azimu.flatten(),
                                     Conv=conversion_factor)
        # fullfit_intens = inp
        # Set the maximum to the n+1th value in the array
        # FIX ME: the max find is not necessarily efficient. A quicker way should be found if it exists.
        max_pos = 0

        ValMax1 = -np.sort(-intens.flatten())[max_pos]  # np.max(intens.flatten())
        ValMax2 = -np.sort(-fullfit_intens.flatten())[max_pos]  # np.max(fullfit_intens.flatten())
        ValMax3 = -np.sort(-(intens.flatten() - fullfit_intens.flatten()))[max_pos]
        # np.max((intens.flatten()-fullfit_intens.flatten()))
        ValMax = np.max([ValMax1, ValMax2])
        ValMin1 = np.min(intens.flatten())
        ValMin2 = np.min(fullfit_intens.flatten())
        ValMin3 = np.min((intens.flatten() - fullfit_intens.flatten()))
        ValMin = np.min([ValMin1, ValMin2])

        y_lims = np.array([np.min(azimu.flatten()), np.max(azimu.flatten())])
        y_lims = np.around(y_lims / 180) * 180

        if conversion_factor['DispersionType'] == 'EnergyDispersive':
            xlabel_str = 'Energy (keV)'
        elif conversion_factor['DispersionType'] == 'AngleDispersive':
            xlabel_str = 'Two theta (deg)'
        else:
            xlabel_str = 'Dispersive units'

        dot_size = 16

        fig = plt.figure()
        axf = fig.add_subplot(1, 3, 1)
        axO1 = plt.subplot(131)
        # plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(intens.flatten()), edgecolors='none',
                    cmap=plt.cm.magma_r, vmin=ValMin, vmax=ValMax)
        axO1.set_title('Data')
        axO1.set_xlabel(xlabel_str)
        axO1.set_ylabel('Azimuth (deg)')
        axO1.set_xlim([np.min(twotheta.flatten()), np.max(twotheta.flatten())])
        axO1.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        plt.colorbar()
        axO2 = plt.subplot(132)
        # plt.scatter(twotheta, azimu, s=1, c=np.log(fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(fullfit_intens.flatten()), edgecolors='none',
                    cmap=plt.cm.magma_r, vmin=ValMin, vmax=ValMax)
        plt.colorbar()
        axO2.set_title('Model')
        axO2.set_xlabel(xlabel_str)
        axO2.set_ylabel('Azimuth (deg)')
        axO2.set_xlim([np.min(twotheta.flatten()), np.max(twotheta.flatten())])
        axO2.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        axO3 = plt.subplot(133)
        # plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(intens.flatten() - fullfit_intens.flatten()),
                    edgecolors='none', cmap=plt.cm.coolwarm, vmin=ValMin3, vmax=ValMax3)
        axO3.set_title('Residuals (data - model)')
        axO3.set_xlabel(xlabel_str)
        axO3.set_ylabel('Azimuth (deg)')
        axO3.set_xlim([np.min(twotheta.flatten()), np.max(twotheta.flatten())])
        axO3.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        plt.colorbar()
        plt.suptitle(peak_string(orders) + '; final fit')
        plt.tight_layout()

        #        # make a second figure with the Q-Q plot of the data and fit.
        #        #the second subplot is a cdf of the residuals with a gaussian added to it.
        #        fig2 = plt.figure()
        #        axQ = plt.subplots()
        #        plt.scatter(intens.flatten(), fullfit_intens.flatten(), s=dot_size)
        #        #axQ.set_xlabel('Observations')
        #        #axQ.set_ylabel('Model')
        #        plt.grid(True)
        #        plt.tight_layout()

        # save figures without overwriting old names
        if SaveFit:
            if not fnam == None:
                filename = os.path.splitext(os.path.basename(fnam))[0] + '_'
            else:
                filename = 'Fit2Peak_'
            filename = filename + peak_string(orders, fname=True)

            i = 0
            if os.path.exists('{}.png'.format(filename)):
                i+=1
            while os.path.exists('{}_{:d}.png'.format(filename, i)):
                i += 1
            if i == 0:
                plt.savefig('{}.png'.format(filename))
            else:
                plt.savefig('{}_{:d}.png'.format(filename, i))

        if view == 1 or debug:
            plt.show()

        plt.close()
        print('Done with Figure')

    # Save lmfit structure
    if SaveFit:
        if fnam is not None:
            filename = os.path.splitext(os.path.basename(fnam))[0] + '_'
        else:
            filename = 'Fit2Peak_'
        filename = filename + peak_string(orders, fname=True)
                
        try:
            save_modelresult(out, filename+'.sav')
        except:
            print('Cannot save lmfit object')
            if os.path.isfile(filename+'.sav'):
                os.remove(filename+'.sav')
            else:    
                print("File doesn't exists!")

    return [NewParams, out]