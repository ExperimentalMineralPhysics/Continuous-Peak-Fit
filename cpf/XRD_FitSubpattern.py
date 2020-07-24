#!/usr/bin/env python

__all__ = ['FitSubpattern']

# FPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier orders

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from lmfit import Parameters, Model
from lmfit.model import save_modelresult#, load_modelresult
import importlib

import cpf.PeakFunctions as ff

np.set_printoptions(threshold=sys.maxsize)


# FIX ME;
# The script should not need to know what the limits are - it should only see the data that it needs to fit.
def peak_string(orders, fname=False):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :return: string listing peak names
    """
    #FIX ME: this will be a data type function so should probably more somewhere else in the end. 
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


def run_initial_fit(azi, y_data, y_data_errs, params, param_str, comp, order, extras, method='leastsq', symm=1):
    """
    :param azi: x-data arr float
    :param y_data: y_data arr float
    :param y_data_errs: y_data error arr float
    :param params: lmfit Parameter object
    :param param_str: base string to identify components str
    :param comp: string to add to base string to identify components str
    :param order: determines number of Fourier coefficients int
    :param method: lmfit method to use str
    :param extras: list of min, max and expression for bounds on parameters
    :return: lmfit Parameter class
    """
    min, max, expr = extras
    params = ff.initiate_params(params, param_str, comp, order, max=max, min=min, expr=expr)
    params = ff.unvary_params(params, param_str, comp)  # set other parameters to not vary
    fout = ff.Fourier_fit(azimu=azi, ydata=np.array(y_data), param=params, errs=np.array(y_data_errs),
                          param_str=param_str + '_' + comp, symm=symm, fit_method=method)
    params = fout.params

    return params


def initial_component_fits(master_params, param_str, data_arrays, comp_lists, order_list, symm, limits, pfixed=None):
    """
    Perform initial Fourier fits to each component data in turn
    :param master_params: lmfit Parameter class object
    :param param_str: start string to label parameters
    :param data_arrays: list of data arrays for use in component fits
    :param comp_lists: list of lists for parameters
    :param order_list: list of orders for components
    :param symm : number definiing the symmetry
    :param limits: array of limits
    :param pfixed=None: either False or true
    :return: list of lists of each component with coefficients and errors
    """

    newAziChunks, newd0, newd0Err, newHall, newHallErr, newWall, newWallErr, newPall, newPallErr = data_arrays
    dfour, hfour, wfour, pfour = comp_lists
    d0_order, h_order, w_order, p_order = order_list
    d_lims, h_lims, w_lims, p_lims = limits

    # initiate symmetry
    comp = 's'
    master_params = ff.initiate_params(master_params, param_str, comp, 0, max=10, min=1, value=symm)

    # FIX ME: need to propagate method properly
    comp = 'd'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newd0),
                                    y_data_errs=np.array(newd0Err), params=master_params, param_str=param_str,
                                    comp=comp, order=d0_order, extras=d_lims, method='leastsq', symm=1)
    dfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'h'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newHall),
                                    y_data_errs=np.array(newHallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=h_order, extras=h_lims, method='leastsq', symm=symm)
    hfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'w'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newWall),
                                    y_data_errs=np.array(newWallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=w_order, extras=w_lims, method='leastsq', symm=symm)
    wfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    # FIX ME: seems to occur regardless of pfixed as checked in other places, is this correct? SAH: NO. Now Fixed.
    comp = 'p'
    if pfixed is False:
        master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newPall),
                                        y_data_errs=np.array(newPallErr), params=master_params, param_str=param_str,
                                        comp=comp, order=p_order, extras=p_lims, method='leastsq', symm=symm)
        #pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    else:
        master_params = ff.initiate_params(master_params, param_str=param_str, comp=comp, value=pfixed)    
    pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()

    #print('Initial fit output parms')
    #print(dfour, hfour, wfour, pfour)
    #master_params.pretty_print()

    return master_params, dfour, hfour, wfour, pfour


def update_component_fits(master_params, intens, azimu, twotheta, param_str, comp_lists, conversion_factor, bgguess,
                          Guesses, pfixed=None):
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
    # master_params.pretty_print()
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, method=None, weights=None, params=master_params)
    # master_params = out.params
    dfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'h'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, method=None, weights=None, params=master_params)
    master_params = out.params
    hfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'w'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, method=None, weights=None, params=master_params)
    master_params = out.params
    wfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    comp = 'p'
    if pfixed is False:
        # if 'profile_fixed' not in order_peak[0]: # orders['peak'][0]:
        master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
        master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
        out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                          num_peaks=len(Guesses['peak']), nterms_back=len(bgguess),
                          Conv=conversion_factor,
                          fixed=pfixed, method=None, weights=None, params=master_params)
        master_params = out.params
    pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))

    return master_params, dfour, hfour, wfour, pfour


def FitSubpattern(TwoThetaAndDspacings, azimu, intens, orders=None, PreviousParams=None,
                  DetFuncs='DioptasFunctions', SaveFit=None, debug=False, refine=True, iterations=1, fnam=None):
    """
    Perform the various fitting stages to the data
    :param TwoThetaAndDspacings:
    :param azimu:
    :param intens:
    :param orders:
    :param PreviousParams:
    :param DetFuncs:
    :param SaveFit:
    :param debug:
    :return:
    """
    # load detector functions to use where needed. Mostly for conversion of tth to d-spacing.
    # This is loaded from detector functions file to avoid making any assumptions about the conversion between
    # values (if indeed any is required).
    # print(DetFuncs)
    # if DetFuncs not in sys.modules:
    det = importlib.import_module(DetFuncs)
    # det = __import__(DetFuncs)

    # sort inputs

    # two theta and d-spacing are both 'x' dimensions. therefore import them as single array.
    # split them up here if needed. 
    if np.shape(TwoThetaAndDspacings)[0] == 3:
        twotheta = TwoThetaAndDspacings[0]
        dspace = TwoThetaAndDspacings[1]
        conversion_factor = TwoThetaAndDspacings[2]
    else:
        twotheta = TwoThetaAndDspacings
        # FIX ME: need to set a value for 'conversion_factor' in here and propagate its none-use through the code.

    if debug: #get limits for chunk fit plots.
        tthrange = [np.min(twotheta.flatten()), np.max(twotheta.flatten())]

    # set data type for the intensity data.
    # This is needed for saving the fits using save_modelresult/ load_modelresult. 
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

    if orders:
        # print orders
        #total = orders['AziBins']
        backg_type = 'order'
        bg_order = orders['background']
        # print(bg_order, 'order')

        #        #FIX ME: assumes a fixed number of peaks.
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
            backg.append(np.zeros(j * 2 + 1))
        # print(backg, 'backg')
        
        
        for y in range(peeks):
            
            # force profile orders to match the fixed orders
            # FIX ME: should this generate an error?
            if 'profile_fixed' in orders['peak'][y]:        
                if not ff.Fourier_order(orders['peak'][y]['profile_fixed']) == orders['peak'][y]['profile']:
                    print('Peak ' + str(y) + ': The order of the profile does not match that of the fixed profile. Changing profile to match fixed profile.')
                    orders['peak'][y]['profile'] = ff.Fourier_order(orders['peak'][y]['profile_fixed'])
        
                #make sure the fixed profile is a list
                if not isinstance(orders['peak'][y]['profile_fixed'], list):
                    orders['peak'][y]['profile_fixed'] = [orders['peak'][y]['profile_fixed']]


    if PreviousParams and orders:
        # need to check sizes of both arrays here.
        # orders takes precedence over PreviousParams so we add or remove coefficients as necessary

        for y in range(peeks):

            # loop over parameters
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

                if orders['peak'][y][Parm] > ff.Fourier_order(PreviousParams['peak'][y][Parm]):
                    # print 'longer'
                    change_by = orders['peak'][y][Parm] * 2 + 1 - np.size(PreviousParams['peak'][y][Parm])
                    # PreviousParams['peak'][y][Parm] = np.pad(PreviousParams['peak'][y][Parm], (0,change_by),
                    # mode='constant', constant_values=0)
                    PreviousParams['peak'][y][Parm] = (PreviousParams['background'][y] + [0] * change_by)

                elif orders['peak'][y][Parm] < ff.Fourier_order(PreviousParams['peak'][y][Parm]):
                    # print 'smaller'
                    change_by = np.size(PreviousParams['peak'][y][Parm]) - ((orders['peak'][y][Parm]) * 2 + 1)
                    PreviousParams['peak'][y][Parm] = PreviousParams['peak'][y][Parm][1:-(change_by - 1)]

            if 'profile_fixed' in orders['peak'][y]:
                PreviousParams['peak'][y]['profile'] = orders['peak'][y]['profile_fixed']
            elif 'profile_fixed' in PreviousParams['peak'][y][Parm]:
                del PreviousParams['peak'][y]['profile_fixed']

        # loop for background orders/size
        for y in range(np.max([len(orders['background']), len(PreviousParams['background'])])):
            # loop over the degree of the background

            if len(PreviousParams['background']) - 1 >= y and len(orders['background']) - 1 >= y:
                # if present in both arrays make sure it is the right size.

                if orders['background'][y] > ff.Fourier_order(PreviousParams['background'][y]):
                    # print 'longer'
                    change_by = orders['background'][y] * 2 + 1 - np.size(PreviousParams['background'][y])
                    # print change_by
                    # PreviousParams['background'][y] = np.pad(PreviousParams['background'][y], (0,change_by),
                    # mode='constant', constant_values=0)
                    PreviousParams['background'][y] = (PreviousParams['background'][y] + [0] * change_by)

                elif orders['background'][y] < ff.Fourier_order(PreviousParams['background'][y]):
                    # print 'smaller'
                    change_by = np.size(PreviousParams['background'][y]) - ((orders['background'][y]) * 2 + 1)
                    PreviousParams['background'][y] = PreviousParams['background'][y][:-(change_by)]

            elif len(PreviousParams['background']) - 1 >= y > len(orders['background']) - 1:
                # if previous params has too many degrees remove the higher ones
                PreviousParams['background'][y] = []

            elif len(PreviousParams['background']) - 1 < y <= len(orders['background']) - 1:
                # if previous parameters does not have enough degrees add more.
                PreviousParams['background'].append([0] * (orders['background'][y] * 2 + 1))

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

    # If there is not previous fit -- Fit data in azimuthal chunks.
    if not PreviousParams:

        if 'PeakPositionSelection' not in orders:
            # DMF changed print statement to include number peaks in inputs file
            if peeks > 1:
                p = 'peaks'
            else:
                p = 'peak'
            print('\nNo previous fits or selections. Assuming ' + str(peeks) + ' ' + p + ' and group fitting in azimuth...\n')
        else:
            print('\nUsing manual selections for initial guesses...\n')

            # FIX ME: Need to check that the number of peaks for which we have parameters is the same as the number
            # of peaks guessed at.

            tthguesses = np.array(orders['PeakPositionSelection'])
            # for future use in setting limits
            dist = np.max(tthguesses[:, 2]) - np.min(tthguesses[:, 2])
            width_change = dist / (2 ** peeks)
            peak_vals = np.unique(tthguesses[:, 2])

            # fit fourier series to two-theta/d-spacing for each peak
            # FIX ME: this fitting should be done to the d-spacing not the tth/energy. Probably need to import the
            # detector functions to make this happen.
            dfour = []
            for j in range(peeks):
                peek = j + 1
                tthguess = tthguesses[tthguesses[:, 0] == peek, 1:]
                param_str = 'peak_' + str(j)
                comp = 'd'

                temp_param = Parameters()
                temp_param = ff.initiate_params(param=temp_param, param_str=param_str, comp=comp,
                                                nterms=orders['peak'][j]['d-space'])
                # temp_param.pretty_print()
                fout = ff.Fourier_fit(azimu=tthguess[:, 0], ydata=tthguess[:, 1],
                                      param=temp_param, param_str=param_str + '_' + comp, fit_method='leastsq')
                temp_param = fout.params
                if debug:
                    temp_param.pretty_print()
                dfour.append(ff.gather_paramerrs_to_list(temp_param, param_str, comp))
                '''
                dfour.append(ff.Fourier_fit_old(tthguess[:, 0], det.Conversion(tthguess, conversion_factor),
                                                terms=orders['peak'][j]['d-space']))
                '''
        
        # get chunks according to detector type. 
        chunks, azichunks = det.bins(azimu, orders)

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

            #print('\nFiting to data chunk ' + str(j + 1) + ' of ' + str(len(chunks)) + '\n')
            print('\rFitting to data chunk %s of %s ' % (str(j + 1), str(len(chunks))), end='\r')
            if debug:
                print('\n')
            # print '\nFitting azimuth chunk....\n'
            # only compute the fit and save it if there are 'enough' data
            
            min_dat = 21  #minimum data in each chunk
            n = 5 #how many data to use in constructing guesses.
            wguess_fraction = 10.
            # N.B. if min_dat/n is too large ther guesses will become very strage. 
            #FIX ME: These should be input variables (for use in extremis)
            
            # FIX ME: possibly do this better using scipy.singal.find or scipy.signal.find_peaks_cwt
             
            if np.sum(~intens.flatten()[chunks[j]].mask) >= min_dat:
                
                # get indicies of sorted two theta values excluding the masked values
                tth_ord = ma.argsort(twotheta.flatten()[chunks[j]].compressed())
                
                                             
                # background estimates
                #make bg list
                backg_guess = [[0.] for i in range(len(orders['background']))]
                
                if backg_type == 'coeffs':
                    # Fourier expand the coefficents. 
                    for j in range(len(orders['background'])):
                        backg_guess[j] = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), param=orders['peak']['background'][j])
                    
                elif backg_type == 'order':
                    #first value (offset) is mean of left hand values
                    backg_guess[0][0] = np.mean(intens.flatten()[chunks[j]].compressed()[tth_ord[:n]])
                    # if there are more then calculate a gradient guess. 
                    # leave all higher order terms as 0 -- assume they are small.
                    if len(orders['background']) > 1:
                        backg_guess[1][0] = (np.mean(intens.flatten()[chunks[j]].compressed()[tth_ord[-n:]]) - np.mean(intens.flatten()[chunks[j]].compressed()[tth_ord[:n]])) / (twotheta.flatten()[chunks[j]].compressed()[tth_ord[-1]] - twotheta.flatten()[chunks[j]].compressed()[tth_ord[0]])
                        
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
                    if 'dfour' in locals(): # if the positions have been pre guessed extract values from dfour.
                        if debug and k == 0:
                            print('dfour in locals')
                        
                        # guess returns tthguess from dfour then converts in into d-spacing.
                        tthguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0])
                        dguess   = det.Conversion(tthguess, conversion_factor, azm=np.mean(azimu.flatten()[chunks[j]]))
                        #finds the index of the closes two theta to the d-spacing input  
                        idx = ((np.abs(twotheta.flatten()[chunks[j]] - tthguess)).argmin())
                        # FIX ME: The mean of a number of the smallest values would be more stable.

                        # height is intensity of slosest pixel in d-spacing - background at that position
                        hguess = intens.flatten()[chunks[j]][idx] 
                        if len(backg_guess) > 1:
                            hguess = hguess - (backg_guess[0][0] + backg_guess[1][0] * (twotheta.flatten()[chunks[j]][idx] - twotheta.flatten()[chunks[j]].min() ))
                        elif len(backg_guess) == 1:
                            hguess = hguess - backg_guess[0][0]
                            
                        # # wguess is fractional width of the data range
                        # wguess = (det.Conversion(np.max(dspace.flatten()[chunks[j]]), conversion_factor, reverse=1, azm=np.mean(azimu.flatten()[chunks[j]])) - 
                        #           det.Conversion(np.min(dspace.flatten()[chunks[j]]), conversion_factor, reverse=1, azm=np.mean(azimu.flatten()[chunks[j]])) ) / wguess_fraction
                        # # FIX ME: This is very crude. Need a better way to do it. \

                        # pguess = 0.5  # guess half if solving for profile
                        # pfixed = 0  # set to 0 so solving unless profile_fixed exists.
                        # if 'profile_fixed' in orders['peak'][k]:
                        #     # FIX ME: profile fixed is the switch 1/0/true/false and profile is the profile parameters. Make sure this is the case in the input files
                        #     pguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), param=orders['peak'][k]['profile_fixed']) 
                        #     pfixed = 1

                        # peaks.append({"d-space": [dguess],
                        #               "height": [hguess],
                        #               "profile": [pguess],
                        #               "width": [wguess]
                        #               })
    
                    else: # guess guesses from data slice - assuming a single peak
                        
                        # find brightest pixels
                        idx = (-intens.flatten()[chunks[j]].compressed()).argsort()[:n]
                        # height guess is nth highest intensity - backgound guess at this position
                        hguess = (intens.flatten()[chunks[j]].compressed())[idx[n-1]]
                        if len(backg_guess) > 1:
                            hguess = hguess - (backg_guess[0][0] + backg_guess[1][0] * (twotheta.flatten()[chunks[j]].compressed()[idx[n-1]] - twotheta.flatten()[chunks[j]].min() ))
                        elif len(backg_guess) == 1:
                            hguess = hguess - backg_guess[0][0]
                            
                        # d-spacing of highest nth intestity pixel
                        dguess = dspace.flatten()[chunks[j]].compressed()[idx[n - 1]]

                    # wguess is fractional width of the data range
                    wguess = (np.max(twotheta.flatten()[chunks[j]]) - np.min(twotheta.flatten()[chunks[j]])) / wguess_fraction / peeks
                    # FIX ME: This is a bit crude. Is there a better way to do it?
    
                    # if profile_fixed exists then set pguess to profile_fixed values and set pfixed to 1
                    
                    pguess = 0.5  # guess half if solving for profile
                    pfixed = 0  # set to 0 so solving unless profile_fixed exists.
                    if 'profile_fixed' in orders['peak'][k]:
                            # profile fixed is the list of coefficients if the profile is fixed. 
                            # FIX ME: profile fixed must have the same number of coefficients as required by profile.
                            pguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]) * orders['peak'][k]['symmetry'], param=orders['peak'][k]['profile_fixed']) 
                            pfixed = 1
    
                    peaks.append({"d-space": [dguess],
                                  "height": [hguess],
                                  "width": [wguess],
                                  "profile": [pguess]
                                  })
                  
                    # DMF added needs checking - required to drive fit and avoid failures
                    lims.append({"d-space": [np.min(dspace.flatten()[chunks[j]]), np.max(dspace.flatten()[chunks[j]])],
                                 "height": [0, np.max(intens.flatten()[chunks[j]])],
                                 "width": [(np.max(twotheta.flatten()[chunks[j]]) - np.min(twotheta.flatten()[chunks[j]])) / 100 / peeks,
                                           (np.max(twotheta.flatten()[chunks[j]]) - np.min(twotheta.flatten()[chunks[j]])) / 2 / peeks],
                                 "profile": [0, 1]
                                 })

                # # Guess for background, if orders and therefore no guess input, choose appropriate
                # singleBackg = [[] for i in range(len(backg))]
                # for i in range(len(backg)):
                #     singleBackg[i] = [backg[i][0]]
                # bgguess = singleBackg
                # bgguess[0][0] = 5 
                # print(bgguess)
                Guesses = {"background": backg_guess, #bgguess,
                           "peak": peaks}
                limits = {"peak": lims}
                
                 # Define parameters to pass to fit
                params = Parameters()
                # if orders bg of form [num_orders_0,num_orders_1]
                # if coeffs bg of form [[coeff1],[coeff2,coeff3,coeff4]]

                # Chunks limited to no Fourier expansion so only single value per polynomial order.
                # for b in range(len(bgguess)):
                #     params.add('bg_c' + str(b) + '_f' + str(0), bgguess[b][0])
                for b in range(len(backg_guess)):
                    params.add('bg_c' + str(b) + '_f' + str(0), backg_guess[b][0])

                for pk in range(len(Guesses['peak'])):

                    if 'dfour' in locals():
                        params.add('peak_' + str(pk) + '_d0', Guesses['peak'][pk]['d-space'][0],
                                   min=limits['peak'][pk]['d-space'][0], max=limits['peak'][pk]['d-space'][1])
                        params.add('peak_' + str(pk) + '_h0', Guesses['peak'][pk]['height'][0],
                                   min=limits['peak'][pk]['height'][0], max=limits['peak'][pk]['height'][1])
                        params.add('peak_' + str(pk) + '_w0', Guesses['peak'][pk]['width'][0],
                                   min=limits['peak'][pk]['width'][0], max=limits['peak'][pk]['width'][1])
                        if 'profile_fixed' in orders['peak'][pk]:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1],
                                       vary=False)
                        else:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1])


                    else:

                        # If using bounds/limits should be written as below:
                        params.add('peak_' + str(pk) + '_d0', Guesses['peak'][pk]['d-space'][0],
                                   min=limits['peak'][pk]['d-space'][0], max=limits['peak'][pk]['d-space'][1])
                        params.add('peak_' + str(pk) + '_h0', Guesses['peak'][pk]['height'][0],
                                   min=limits['peak'][pk]['height'][0], max=limits['peak'][pk]['height'][1])
                        params.add('peak_' + str(pk) + '_w0', Guesses['peak'][pk]['width'][0],
                                   min=limits['peak'][pk]['width'][0], max=limits['peak'][pk]['width'][1])
                        if 'profile_fixed' in orders['peak'][pk]:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1],
                                       vary=False)
                        else:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1])


                
                if debug:
                    params.pretty_print()
                    print('\n)')
                    guess = params
                out = ff.FitModel(intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]], azichunks[j],
                                  num_peaks=len(Guesses['peak']), nterms_back=len(backg_guess), Conv=conversion_factor,
                                  fixed=pfixed, method=None, weights=None, params=params)

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
                    # asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
                    # tth_plot = twotheta.flatten()[chunks[j]][asdf_ord]
                    # int_plot = intens.flatten()[chunks[j]][asdf_ord]
                    #azm_plot = azimu.flatten()[chunks[j]][asdf_ord]
                    tth_plot = twotheta.flatten()[chunks[j]]
                    int_plot = intens.flatten()[chunks[j]]
                    azm_plot = np.tile(azimu.flatten()[chunks[j]][0], (300))
                    azm_plot = ma.array(azm_plot, mask=(~np.isfinite(azm_plot))) #FIX ME: required for plotting energy dispersive data.
                    gmodel = Model(ff.PeaksModel, independent_vars=['twotheta', 'azi'])
                    tth_range = np.linspace(np.min(twotheta.flatten()[chunks[j]]), np.max(twotheta.flatten()[chunks[j]]), azm_plot.size)
                    # mod_plot = gmodel.eval(params=params, twotheta=tth_range, azi=azm_plot,
                    #                        num_peaks=len(Guesses['peak']),
                    #                        nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed)
                    # guess_plot = gmodel.eval(params=guess, twotheta=tth_range, azi=azm_plot,
                    #                        num_peaks=len(Guesses['peak']),
                    #                        nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed)
                    mod_plot = gmodel.eval(params=params, twotheta=tth_range, azi=azm_plot, Conv=conversion_factor)
                    guess_plot = gmodel.eval(params=guess, twotheta=tth_range, azi=azm_plot, Conv=conversion_factor)
                    plt.plot(tth_plot, int_plot, '.', label="data")
                    plt.plot(tth_range, guess_plot, marker='', color='green', linewidth=2, linestyle='dashed', label='guess')
                    plt.plot(tth_range, mod_plot, marker='', color='red', linewidth=2, label='fit')
                    plt.xlim(tthrange)
                    plt.legend()
                    plt.title(peak_string(orders) + '; Chunk ' + str(j+1))
                    plt.show()
                    
                    

        # Feed each d_0,h,w into fourier expansion function to get fit for fourier component
        # parameters as output.

        print('\n Doing Fourier fits...')
        # print('test0')

        dfour = []
        hfour = []
        wfour = []
        pfour = []
        master_params = Parameters()  # initiate total Parameter class to add to
        # FIX ME: parameters initiated with no limits currently 
        for j in range(peeks):

            # symmetry
            if 'symmetry' in orders['peak'][j].keys():
                symm = orders['peak'][j]['symmetry']
            else:
                symm = 1

            d0_order = orders['peak'][j]['d-space']
            h_order = orders['peak'][j]['height']
            w_order = orders['peak'][j]['width']
            p_order = orders['peak'][j]['profile']
            if 'profile_fixed' in orders['peak'][j]:
                p_fixed = orders['peak'][j]['profile_fixed']
            else:
                p_fixed = False
            # print(np.array(newAziChunks), np.array(newd0[j]), d0_order, np.array(newd0Err[j]))

            param_str = 'peak_' + str(j)  # defines peak string to start parameter name

            # Initiate parameters for each component and perform initial Fourier fit
            data_arrays = newAziChunks, newd0[j], newd0Err[j], newHall[j], newHallErr[j], newWall[j], newWallErr[j], \
                          newPall[j], newPallErr[j]
            fit_limits = [[None, None, None], [None, None, None], [None, None, None], [None, None, None]]
            # fit_limits = [[np.min(newd0),np.max(newd0),None], [0,None,None],
            #              [0,(np.max(newWall)+np.min(newWall))/2.,None], [0,1,None]] # min, max, expr in order d,
            # h, w, p

            master_params, dfour, hfour, wfour, pfour = initial_component_fits(master_params, param_str, data_arrays,
                                                                               [dfour, hfour, wfour, pfour],
                                                                               [d0_order, h_order, w_order, p_order],
                                                                               symm=symm, limits=fit_limits, pfixed=p_fixed)
            if debug:
                print('Parameters after initial Fourier fits')
                master_params.pretty_print()

        newBGall = np.array(newBGall)
        newBGallErr = np.array(newBGallErr)
        bgfour = []


        # Initiate background parameters
        comp = 'bg'
        for b in range(len(backg)):
            for f_b in range(len(backg[b])):
                master_params.add('bg_c' + str(b) + '_f' + str(f_b), backg[b][f_b])  # need limits

        # Perform initial background Fourier fit
        for i in range(len(lenbg)):
            param_str = 'bg_c' + str(i)
            master_params = ff.unvary_params(master_params, param_str, '')
            master_params = ff.vary_params(master_params, param_str, '')
            fout = ff.Fourier_fit(azimu=np.array(newAziChunks), ydata=np.array(newBGall[i]), param=master_params,
                                  errs=np.array(newBGallErr[i]), param_str=param_str, fit_method='leastsq')
            master_params = fout.params
            bgfour.append(ff.gather_paramerrs_to_list(master_params, param_str, 'f'))

        # plot output of fourier fits....
        if debug:
            y_lims = np.array([np.min(newAziChunks), np.max(newAziChunks)])
            y_lims = np.around(y_lims / 180) * 180
            AziPlot = range(np.int(y_lims[0]), np.int(y_lims[1]), 2)
            gmodel = Model(ff.Fourier_expand, independent_vars=['azimu'])

            fig = plt.figure()
            ax1 = fig.add_subplot(5, 1, 1)
            ax1.set_title('D-spacing')
            for j in range(peeks):
                ax1.scatter(newAziChunks, newd0[j], s=10)
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot), param=dfour[j][0])
                ax1.plot(AziPlot, gmod_plot)

            # plt.subplot(512)
            ax2 = fig.add_subplot(5, 1, 2)
            ax2.set_title('height')
            for j in range(peeks):
                if 'symmetry' in orders['peak'][j].keys():
                    symm = orders['peak'][j]['symmetry']
                else:
                    symm = 1
                ax2.scatter(newAziChunks, newHall[j], s=10)
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot) * symm, param=hfour[j][0])
                ax2.plot(AziPlot, gmod_plot)

            ax3 = fig.add_subplot(5, 1, 3)
            # plt.subplot(513)
            ax3.set_title('width')
            for j in range(peeks):
                if 'symmetry' in orders['peak'][j].keys():
                    symm = orders['peak'][j]['symmetry']
                else:
                    symm = 1
                ax3.scatter(newAziChunks, newWall[j], s=10)
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot) * symm, param=wfour[j][0])
                ax3.plot(AziPlot, gmod_plot)

            ax4 = fig.add_subplot(5, 1, 4)
            # plt.subplot(514)
            ax4.set_title('Profile')
            for j in range(peeks):
                if 'symmetry' in orders['peak'][j].keys():
                    symm = orders['peak'][j]['symmetry']
                else:
                    symm = 1
                ax4.scatter(newAziChunks, newPall[j], s=10)
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot) * symm, param=pfour[j][0])
                ax4.plot(AziPlot, gmod_plot)

            for k in range(len(lenbg)):
                x_plt = len(lenbg)
                y_plt = len(lenbg) * 4 + k + 1
                ax5 = fig.add_subplot(5, x_plt, y_plt)
                ax5.set_title('Background')
                ax5.scatter(newAziChunks, newBGall[k], s=10)
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot), param=bgfour[k][0])
                ax5.plot(AziPlot, gmod_plot)

            fig.suptitle(peak_string(orders) + '; Fits to Chunks')
            plt.show()
            plt.close()

        # put fits into NewParams data structure.
        NewParams = ff.create_newparams(peeks, dfour, hfour, wfour, pfour, bgfour, orders['peak'])

        # Iterate over each fourier series in turn.

        if refine:
            print('\nRe-fitting for d, h, w, bg separately...\n')
            # iterate through the variables.
            for j in range(iterations):

                dfour = []
                hfour = []
                wfour = []
                pfour = []
                # refine parameters for each peak.
                for k in range(peeks):

                    param_str = 'peak_' + str(k)
                    if 'profile_fixed' in orders['peak'][k]:
                        p_fixed = True
                    else:
                        p_fixed = False
                    if 'symmetry' in orders['peak'][k]:
                        symm = orders['peak'][k]['symmetry']
                    else:
                        symm = 1
                    master_params, dfour, hfour, wfour, pfour = update_component_fits(master_params, intens, azimu,
                                                                                      twotheta, param_str,
                                                                                      [dfour, hfour, wfour, pfour],
                                                                                      conversion_factor, backg_guess,
                                                                                      Guesses, p_fixed)

                # refine background
                param_str = 'bg_c'
                # print(param_str, lenbg, i)
                master_params = ff.unvary_params(master_params, param_str, '')
                master_params = ff.vary_params(master_params, param_str, '')
                out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), num_peaks=len(Guesses['peak']),
                                  nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed, method=None,
                                  weights=None, params=master_params)
                master_params = out.params
                # Put fits into data structure.
                NewParams = ff.create_newparams(peeks, dfour, hfour, wfour, pfour, bgfour, orders['peak'])
        
                if debug:
                    print('Parameters after refining Fourier fits ' + str(j) + 'times')
                    master_params.pretty_print()


    else:
        # FIX ME: This should load the saved lmfit Parameter class object. Need to initiate usage (save and load).
        # For now have propagated use of NewParams so re-instantiate the master_params object below, which is clunky.
        print('Using previously fitted parameters and propagating fit')
        NewParams = PreviousParams

        # FIX ME: should this be some def or subfunction?
        master_params = Parameters()  # initiate total Parameter class to add to
        # FIX ME: parameters initiated with no limits currently 

        for j in range(peeks):

            param_str = 'peak_' + str(j)  # defines peak string to start parameter name

            # initiate symmetry
            comp = 's'
            if 'symmetry' in orders['peak'][j].keys():
                symm = orders['peak'][j]['symmetry']
            else:
                symm = 1
            master_params = ff.initiate_params(master_params, param_str, comp, 0, max=10, min=1, value=np.array(symm))

            # initiate d
            comp = 'd'
            dfour = PreviousParams['peak'][j]['d-space']
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['d-space'],
                                               value=np.array(dfour))

            # initiate h
            comp = 'h'
            hfour = PreviousParams['peak'][j]['height']
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['height'], value=hfour)

            # initiate w
            comp = 'w'
            wfour = PreviousParams['peak'][j]['width']
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['width'], value=wfour)

            # initiate p
            comp = 'p'
            if orders and 'profile_fixed' in orders['peak'][j]:
                pfour = orders['peak'][j]['profile_fixed']
                pfixed = 1
            else:
                pfour = PreviousParams['peak'][j]['profile']
                pfixed = 0
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['profile'], value=pfour)

        # Initiate background parameters
        comp = 'bg'
        for b in range(len(backg)):
            # print(backg[b])
            for f_b in range(len(backg[b])):
                # print(b, f_b)
                master_params.add('bg_c' + str(b) + '_f' + str(f_b),
                                  PreviousParams['background'][b][f_b])  # need limits

        # master_params.pretty_print()

        # Guess for background, if orders and therefore no guess input, choose appropriate
        # copied from above to allow the code to run. FIX ME: doesn't need all this code.
        singleBackg = [[] for i in range(len(backg))]
        for i in range(len(backg)):
            singleBackg[i] = [backg[i][0]]
        backg_guess = singleBackg
        #bgguess[0][0] = 5

        # def previous_params_to_params(PreviousParams):

        # FIX ME: need to confirm the number of parameters matches the orders of the fits.

    # FIX ME: Need to sort out use of lmfit parameters and NewParams, below will not work currently for loaded
    # parameters as master_params not created from NewParams yet. If possible use the json import/export within lmfit.

    # Refit through full equation with all data for d,h,w,bg independently
    print('\nFinal fit solving for all parms...\n')

    num_peaks = len(orders['peak'])

    # set all parameters to vary
    for p in master_params:
        master_params[p].set(vary=True)
    for x in range(num_peaks):
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
    
    # FIX ME: need to sort parameters, bgguess doesn't exist if load previous fit data
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), num_peaks=num_peaks,
                      nterms_back=len(backg_guess), Conv=conversion_factor, fixed=pfixed, method=None,
                      weights=None, params=master_params)
    master_params = out.params
    
    print('\nFinal Coefficients\n')
    # print(out.fit_report(min_correl=1))
    print(out.fit_report())
    # print('Final Coefficients\n', NewParams)
    
    # Write master_params to NewParams dict object
    NewParams = ff.params_to_newparams(master_params, num_peaks=num_peaks, num_bg=len(backg_guess),
                                       order_peak=orders['peak'])

    # Stats now directly from lmfit info.
    # number elements
    n_points = out.ndata
    print('number data', n_points)

    # degrees of freedom (number of parameters solved for).
    deg_freedom = out.nfree
    print('degrees of freedom', deg_freedom)

    # ChiSq from fit
    ChiSq = out.chisqr
    print('ChiSquared', ChiSq)

    '''
    # Calculate SSD for fit.
    # inp, pen = ff.PeaksModel(twotheta.flatten(), azimu.flatten(), NewParams, Conv=conversion_factor)
    SSD = np.sum((intens.flatten() - inp.flatten()) ** 2)
    print('SSD', SSD)

    # https://autarkaw.org/2008/07/05/finding-the-optimum-polynomial-order-to-use-for-regression/
    # We choose the degree of polynomial for which the variance as computed by
    # Sr(m)/(n-m-1)
    # is a minimum or when there is no significant decrease in its value as the degree of polynomial is increased.
    # In the above formula, Sr(m) = sum of the square of the residuals for the mth order polynomial
    # n= number of data points
    # m=order of polynomial (so m+1 is the number of constants of the model)

    # number elements
    n_points = intens.flatten().size
    print('number data', n_points)

    # degrees of freedom (number of parameters solved for).
    deg_freedom = sum(ff.flatten(Chng))
    print('degrees of freedom', deg_freedom)

    SSD_var = SSD / (n_points + deg_freedom - 1)
    print('variance', SSD_var)

    # ChiSq is defined in the same way as GSAS. See GSAS manual page 165/166
    ChiSq = SSD / (n_points - deg_freedom)
    print('ChiSquared', ChiSq)
    '''

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
        max_pos = 5

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

    #save lmfit structure
    if SaveFit:
        if not fnam == None:
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
        
    tth_min = twotheta.min()
    tth_max = twotheta.max()
    d_min = dspace.min()
    d_max = dspace.max()
    extent = [[tth_min, tth_max], [d_min, d_max]]  # FIX ME: This is taking the maximum and the minimum of the data
    # not the 'range' itself.

    NewParams.update({'range': extent})

    FitStats = {"degree-of-freedom": deg_freedom,
                "n-points": n_points,
                "ChiSq": ChiSq
                }

    NewParams.update({'FitProperties': FitStats})

    # FIX ME: Need to put master_params into NewParams to return!

    return [ NewParams, out]
