#!/usr/bin/env python

__all__ = ['FitSubpattern']

# FPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier order

# import numpy.ma as ma
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters, Model
import importlib

import cpf.PeakFunctions as ff

np.set_printoptions(threshold=sys.maxsize)


# FIX ME;
# The script should not need to know what the limits are - it should only see the data that it needs to fit.

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
    # dfour.append(ff.gather_paramerrs_to_list(params, param_str, comp))
    # params.pretty_print()

    return params


def initial_component_fits(master_params, param_str, data_arrays, comp_lists, order_list, symm, limits):
    """
    Perform initial Fourier fits to each component data in turn
    :param master_params: lmfit Parameter class object
    :param param_str: start string to label parameters
    :param data_arrays: list of data arrays for use in component fits
    :param comp_lists: list of lists for parameters
    :param order_list: list of orders for components
    :return: list of lists of each component with coefficients and errors
    """

    newAziChunks, newd0, newd0Err, newHall, newHallErr, newWall, newWallErr, newPall, newPallErr = data_arrays
    dfour, hfour, wfour, pfour = comp_lists
    d0_order, h_order, w_order, p_order = order_list
    d_lims, h_lims, w_lims, p_lims = limits
    print(d_lims)
    print(symm)

    # initiate symmetry
    comp = 's'
    master_params = ff.initiate_params(master_params, param_str, comp, 0, max=10, min=1, value=symm)
    # master_params.pretty_print()

    # FIX ME: need to propagate method properly
    comp = 'd'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newd0),
                                    y_data_errs=np.array(newd0Err), params=master_params, param_str=param_str,
                                    comp=comp, order=d0_order, extras=d_lims, method='leastsq', symm=1)
    dfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()

    comp = 'h'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newHall),
                                    y_data_errs=np.array(newHallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=h_order, extras=h_lims, method='leastsq', symm=symm)
    hfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    master_params.pretty_print()

    comp = 'w'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newWall),
                                    y_data_errs=np.array(newWallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=w_order, extras=w_lims, method='leastsq', symm=symm)
    wfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()

    # FIX ME: seems to occur regardless of pfixed as checked in other places, is this correct?
    comp = 'p'
    master_params = run_initial_fit(azi=np.array(newAziChunks), y_data=np.array(newPall),
                                    y_data_errs=np.array(newPallErr), params=master_params, param_str=param_str,
                                    comp=comp, order=p_order, extras=p_lims, method='leastsq', symm=symm)
    pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()

    print('Initial fit output parms')
    print(dfour, hfour, wfour, pfour)
    master_params.pretty_print()

    '''
    dfour.append(ff.Fourier_fit(np.array(newAziChunks), np.array(newd0[j]), terms=d0_order,
                                errs=np.array(newd0Err[j])))
    hfour.append(ff.Fourier_fit(np.array(newAziChunks) * symm, np.array(newHall[j]), terms=h_order,
                                errs=np.array(newHallErr[j])))
    wfour.append(ff.Fourier_fit(np.array(newAziChunks) * symm, np.array(newWall[j]), terms=w_order,
                                errs=np.array(newWallErr[j]) * np.array(newWallErr[j])))
    pfour.append(ff.Fourier_fit(np.array(newAziChunks) * symm, np.array(newPall[j]), terms=p_order,
                                errs=np.array(newPallErr[j]) * np.array(newPallErr[j])))
    '''

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
    # print('d after fit')
    # master_params.pretty_print()

    comp = 'h'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, method=None, weights=None, params=master_params)
    master_params = out.params
    hfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()

    comp = 'w'
    master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
    master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                      num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                      fixed=pfixed, method=None, weights=None, params=master_params)
    master_params = out.params
    wfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    # master_params.pretty_print()

    if pfixed:
        comp = 'p'
        pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
    else:
        # if 'profile_fixed' not in order_peak[0]: # orders['peak'][0]:
        comp = 'p'
        master_params = ff.unvary_params(master_params, param_str, comp)  # set other parameters to not vary
        master_params = ff.vary_params(master_params, param_str, comp)  # set these parameters to vary
        out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(),
                          num_peaks=len(Guesses['peak']), nterms_back=len(bgguess),
                          Conv=conversion_factor,
                          fixed=pfixed, method=None, weights=None, params=master_params)
        master_params = out.params
        pfour.append(ff.gather_paramerrs_to_list(master_params, param_str, comp))
        # master_params.pretty_print()

    return master_params, dfour, hfour, wfour, pfour


def FitSubpattern(TwoThetaAndDspacings, azimu, intens, orders=None, PreviousParams=None,
                  DetFuncs='DioptasFunctions', SaveFit=None, debug=False, refine=True, iterations=1):
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
        total = orders['AziBins']
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
            print('\nNo previous fits or selections. Assuming str(peeks) ' + p + ' and group fitting in azimuth...\n')
        else:
            print('\nUsing manual selections for initial guesses...\n')

            # FIX ME: Need to check that the number of peaks for which we have parameters is the same as the number
            # of peaks guessed at.

            tthguesses = np.array(orders['PeakPositionSelection'])
            # for future use in setting limits
            print(tthguesses)
            dist = np.max(tthguesses[:, 2]) - np.min(tthguesses[:, 2])
            print(dist, peeks)
            width_change = dist / (2 ** peeks)
            peak_vals = np.unique(tthguesses[:, 2])
            print(width_change, peak_vals)

            # fit fourier series to two-theta/d-spacing for each peak
            # FIX ME: this fitting should be done to the d-spacing not the tth/energy. Probably need to import the
            # detector functions to make this happen.
            dfour = []
            for j in range(peeks):
                peek = j + 1
                tthguess = tthguesses[tthguesses[:, 0] == peek, 1:]
                print('peak', j)
                print(tthguess[:, 0])
                print(det.Conversion(tthguess, conversion_factor))
                print('terms', orders['peak'][j]['d-space'])
                param_str = 'peak_' + str(j)
                comp = 'd'

                temp_param = Parameters()
                temp_param = ff.initiate_params(param=temp_param, param_str=param_str, comp=comp,
                                                nterms=orders['peak'][j]['d-space'])
                # temp_param.pretty_print()
                fout = ff.Fourier_fit(azimu=tthguess[:, 0], ydata=det.Conversion(tthguess, conversion_factor),
                                      param=temp_param, param_str=param_str + '_' + comp, fit_method='leastsq')
                temp_param = fout.params
                temp_param.pretty_print()
                dfour.append(ff.gather_paramerrs_to_list(temp_param, param_str, comp))
                '''
                dfour.append(ff.Fourier_fit_old(tthguess[:, 0], det.Conversion(tthguess, conversion_factor),
                                                terms=orders['peak'][j]['d-space']))
                '''
                # print('tthguess:', tthguess)
                print('dfour:', dfour)

        azmax = azimu.max()
        azmin = azimu.min()
        binsize = (azmax - azmin) / total
        chunks = []
        azichunks = []

        for i in range(total):

            end = azmin + (i + 1) * binsize
            start = azmin + i * binsize

            tempazi = azimu.flatten()
            azichunk = np.where((tempazi > start) & (tempazi <= end))

            if np.unique(tempazi[azichunk]).size == 1:
                # FIXME: this is a quick way of distinguishing between angle
                # and energy dispersive data without having to know anything about hte number of detectors.... should
                # be made more robust because the azimuth is needed to convert energy to d-spacing and wont work if
                # it is wrong.
                azichunks.append(tempazi[azichunk[0][0]])
            else:
                azichunks.append(((end - start) / 2) + start)
            chunks.append(azichunk)

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

            # print '\nFitting azimuth chunk....\n'
            # only compute the fit and save it if there are 'enough' data
            if np.sum(~intens.flatten()[chunks[j]].mask) >= 20:
                # FIX ME: this switch is a crude way of finding if
                # there are less than 20 data in the bin. It should be a variable.

                # background estimates ##
                if backg_type == 'coeffs':
                    backg_base = backg[0][0]
                    # print(backg_type, backg_base, backg)
                elif backg_type == 'order':
                    backg_base = intens.flatten()[chunks[j]].argmin()
                    backg_base = intens[backg_base]
                    if not isinstance(backg_base, int):
                        backg_base = 0
                    backg[0][0] = backg_base
                    # print(backg_type, backg_base, backg)
                elif backg_type == 'flat':
                    backg_base = backg[0][0]
                    # print(backg_type, backg_base, backg)

                # Organise guesses to be refined.
                dguess = []
                hguess = []
                wguess = []
                pguess = []
                pfixed = 0
                peaks = []
                lims = []

                if 'dfour' in locals():
                    # if the positions have been pre guessed extract values from dfour.
                    print('dfour in locals')
                    for k in range(peeks):

                        # guess returns d-spacing then converted to two theta.
                        dguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0])
                        tthguess = (det.Conversion(dguess, conversion_factor, reverse=1))  # FIXME: might break here
                        # because of how changed Conversion to work with energy dispersive diffraction.

                        tthind = ((np.abs(twotheta.flatten()[chunks[j]] - tthguess)).argmin())
                        # FIX ME: this takes the closest tth value to get the intensity. The mean of a number of the
                        # smallest values would be more stable.

                        hguess = (intens.flatten()[chunks[j]][tthind] - backg_base)  # FIX ME: This is very crude.
                        # Need a better way to do it.
                        wguess = ((np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]])) / 40)
                        print(wguess)
                        # FIX ME: This is very crude. Need a better way to do it. \

                        pguess = 0.5  # guess half if solving for profile
                        pfixed = 0  # set to 0 so solving unless profile_fixed exists.
                        if 'profile_fixed' in orders['peak'][k]:
                            pguess = orders['peak'][k]['profile_fixed']  # FIX ME: this should be a Fourier expansion
                            # because it is fixed....
                            pfixed = 1

                        peaks.append({"d-space": [dguess],
                                      "height": [hguess],
                                      "profile": [pguess],
                                      "width": [wguess]
                                      })

                        # DMF added needs checking - required to drive fit and avoid failures
                        lims.append(
                            {"d-space": [np.min(dspace.flatten()[chunks[j]]), np.max(dspace.flatten()[chunks[j]])],
                             "height": [0, np.inf],
                             "profile": [0, 1],
                             "width": [0, (np.max(dspace.flatten()[chunks[j]]) +
                                           np.min(dspace.flatten()[chunks[j]])) / 2]})


                else:
                    # guess guesses from data slice - assuming a single peak
                    hguess_I = intens.flatten()[chunks[j]].argmax()
                    hguess = intens.flatten()[chunks[j]][hguess_I] * .9 - backg_base  # FIX ME: This is very crude.
                    # Need a better way to do it.

                    n = 5
                    idx = (-intens.flatten()[chunks[j]]).argsort()[:n]
                    hguess = (intens.flatten()[chunks[j]])[idx[n - 1]] - 5

                    # dguess = dspace.flatten()[chunks[j]][hguess_I]
                    dguess = dspace.flatten()[chunks[j]][idx[n - 1]]

                    wguess = 0.02  # (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/40
                    # FIX ME: This is very crude. Need a better way to do it.

                    # if profile_fixed exists then set pguess to fixed value
                    # then set pfixed, 1 if fixed else 0
                    pguess = 0.5  # guess half if solving for profile
                    pfixed = 0  # set to 0 so solving unless profile_fixed exists.

                    if 'profile_fixed' in orders['peak'][0]:
                        pguess = orders['peak'][0]['profile_fixed']
                        # FIX ME: this should be a Fourier expansion because it is fixed....
                        pfixed = 1

                    peaks.append({"d-space": [dguess],
                                  "height": [hguess],
                                  "profile": [pguess],
                                  "width": [wguess]
                                  })
                    lims.append({"d-space": [np.min(dspace.flatten()[chunks[j]]), np.max(dspace.flatten()[chunks[j]])],
                                 "height": [0, np.inf],
                                 "profile": [0, 1],
                                 "width": [0, (np.max(dspace.flatten()[chunks[j]]) +
                                               np.min(dspace.flatten()[chunks[j]])) / 2]})

                # Guess for background, if orders and therefore no guess input, choose appropriate
                singleBackg = [[] for i in range(len(backg))]
                for i in range(len(backg)):
                    singleBackg[i] = [backg[i][0]]
                bgguess = singleBackg
                bgguess[0][0] = 5

                Guesses = {"background": bgguess,
                           "peak": peaks}
                limits = {"peak": lims}

                '''
                # print('paramfitarr in', peeks, len(backg), Guesses)
                if pfixed == 1:
                    # FIX ME if True, don't fit for the peak profile?
                    Chng = ff.ParamFitArr(peeks, len(backg), 1, p=[1, 1, 1, 0])
                else:
                    Chng = ff.ParamFitArr(peeks, len(backg), 1)
                    # Fit for all - peak - dspace, height, profile, width and bg
                # print(Chng)
                '''

                # Define parameters to pass to fit
                params = Parameters()
                # if orders bg of form [num_orders_0,num_orders_1]
                # if coeffs bg of form [[coeff1],[coeff2,coeff3,coeff4]]

                # Chunks limited to no Fourier expansion so only single value per polynomial order.
                for b in range(len(bgguess)):
                    params.add('bg_c' + str(b) + '_f' + str(0), bgguess[b][0])

                for pk in range(len(Guesses['peak'])):

                    if 'dfour' in locals():
                        # SAH need to add limits to the fitter here.
                        ''' old code
                        params.add('peak_'+str(pk)+'_d0', Guesses['peak'][pk]['d-space'][0])
                        params.add('peak_'+str(pk)+'_h0', Guesses['peak'][pk]['height'][0])
                        if pfixed == 1:
                            params.add('peak_'+str(pk)+'_p0', Guesses['peak'][pk]['profile'][0],
                                       vary=False)
                        else:
                            params.add('peak_'+str(pk)+'_p0', Guesses['peak'][pk]['profile'][0])
                        params.add('peak_'+str(pk)+'_w0', Guesses['peak'][pk]['width'][0])
                        '''
                        params.add('peak_' + str(pk) + '_d0', Guesses['peak'][pk]['d-space'][0],
                                   min=limits['peak'][pk]['d-space'][0], max=limits['peak'][pk]['d-space'][1])
                        params.add('peak_' + str(pk) + '_h0', Guesses['peak'][pk]['height'][0],
                                   min=limits['peak'][pk]['height'][0], max=limits['peak'][pk]['height'][1])
                        if pfixed == 1:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1],
                                       vary=False)
                        else:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1])
                        params.add('peak_' + str(pk) + '_w0', Guesses['peak'][pk]['width'][0],
                                   min=limits['peak'][pk]['width'][0], max=limits['peak'][pk]['width'][1])


                    else:

                        # If using bounds/limits should be written as below:
                        params.add('peak_' + str(pk) + '_d0', Guesses['peak'][pk]['d-space'][0],
                                   min=limits['peak'][pk]['d-space'][0], max=limits['peak'][pk]['d-space'][1])
                        params.add('peak_' + str(pk) + '_h0', Guesses['peak'][pk]['height'][0],
                                   min=limits['peak'][pk]['height'][0], max=limits['peak'][pk]['height'][1])
                        if pfixed == 1:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1],
                                       vary=False)
                        else:
                            params.add('peak_' + str(pk) + '_p0', Guesses['peak'][pk]['profile'][0],
                                       min=limits['peak'][pk]['profile'][0], max=limits['peak'][pk]['profile'][1])
                        params.add('peak_' + str(pk) + '_w0', Guesses['peak'][pk]['width'][0],
                                   min=limits['peak'][pk]['width'][0], max=limits['peak'][pk]['width'][1])

                # print('Call fit', Chng, Guesses)
                # print('this', len(Guesses['peak']))
                # print(stop)

                params.pretty_print()
                out = ff.FitModel(intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]], azichunks[j],
                                  num_peaks=len(Guesses['peak']), nterms_back=len(bgguess), Conv=conversion_factor,
                                  fixed=pfixed, method=None, weights=None, params=params)
                # Guesses, po, pc = ff.FitModel_old(intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]],
                #                                  azichunks[j], Chng, Guesses, fixed=pfixed, Conv=conversion_factor,
                #                                  bounds=limits)
                # In above, limits passed but not used, check desired functionality

                print('\nFit to data chunk ' + str(j + 1) + ' of ' + str(len(chunks)) + '\n')
                # print(out.fit_report())
                params = out.params
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
                for i in range(len(bgguess)):
                    newBGall[i].append(params['bg_c' + str(i) + '_f0'].value)
                    newBGallErr[i].append(params['bg_c' + str(i) + '_f0'].stderr)

                    '''
                    newd0[i].append(Guesses['peak'][i]['d-space'][0])
                    newd0Err[i].append(Guesses['peak'][i]['d-space_err'][0])
                    newHall[i].append(Guesses['peak'][i]['height'][0])
                    newHallErr[i].append(Guesses['peak'][i]['height_err'][0])
                    newWall[i].append(Guesses['peak'][i]['width'][0])
                    newWallErr[i].append(Guesses['peak'][i]['width_err'][0])
                    if Guesses['peak'][i]['profile'][0] > 1:
                        # FixMe: This limits the possible widths of the peaks in the initial guess. Should be
                        # replaced by constrained fit.
                        newPall[i].append(0.9)
                    elif Guesses['peak'][i]['profile'][0] < 0:
                        newPall[i].append(0.1)
                    else:
                        newPall[i].append(Guesses['peak'][i]['profile'][0])
                    if pfixed == 1:
                        newPallErr[i].append(0.0001)
                    else:
                        newPallErr[i].append(Guesses['peak'][i]['profile_err'][0])
                    '''
                '''
                for i in range(len(Guesses['background'])):
                    newBGall[i].append(Guesses['background'][i][0])
                    newBGallErr[i].append(Guesses['background_err'][i][0])
                '''
                newAziChunks.append(azichunks[j])

                # plot the fits.
                if debug:
                    asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
                    tth_plot = twotheta.flatten()[chunks[j]][asdf_ord]
                    int_plot = intens.flatten()[chunks[j]][asdf_ord]
                    azm_plot = azimu.flatten()[chunks[j]][asdf_ord]
                    gmodel = Model(ff.PeaksModel, independent_vars=['twotheta', 'azi'])
                    mod_plot = gmodel.eval(params=params, twotheta=tth_plot, azi=azm_plot,
                                           num_peaks=len(Guesses['peak']),
                                           nterms_back=len(bgguess), Conv=conversion_factor, fixed=pfixed)
                    # mod_plot = ff.PeaksModel(tth_plot, azm_plot, Guesses, Conv=conversion_factor)

                    plt.plot(tth_plot, int_plot, '.')
                    plt.plot(tth_plot, mod_plot, 'r.-')
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
        # FIX ME: parameters initiated with no limits currently -- SAH: I think I have fixed this now. no i have not.
        # silly boy.
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
                                                                               symm=symm, limits=fit_limits)

        newBGall = np.array(newBGall)
        newBGallErr = np.array(newBGallErr)
        bgfour = []

        # print('Before initiate background params')
        # master_params.pretty_print()

        # Initiate background parameters
        comp = 'bg'
        for b in range(len(backg)):
            # print(backg[b])
            for f_b in range(len(backg[b])):
                # print(b, f_b)
                master_params.add('bg_c' + str(b) + '_f' + str(f_b), backg[b][f_b])  # need limits
        # master_params.pretty_print()

        # Perform initial background Fourier fit
        for i in range(len(lenbg)):
            param_str = 'bg_c' + str(i)
            # print(param_str, lenbg, i)
            master_params = ff.unvary_params(master_params, param_str, '')
            master_params = ff.vary_params(master_params, param_str, '')
            # master_params.pretty_print()
            fout = ff.Fourier_fit(azimu=np.array(newAziChunks), ydata=np.array(newBGall[j]), param=master_params,
                                  errs=np.array(newBGallErr[j]), param_str=param_str, fit_method='leastsq')
            # tempbg, tempbgc = ff.Fourier_fit(np.array(newAziChunks), np.array(newBGall[i]), terms=bg_order[i],
            #                                 errs=np.array(newBGallErr[i]))
            # bgfour.append(tempbg.tolist())
            master_params = fout.params
            bgfour.append(ff.gather_paramerrs_to_list(master_params, param_str, 'f'))
            # print('prettyprint after')
        master_params.pretty_print()

        # plot output of fourier fits....
        if debug:
            y_lims = np.array([np.min(newAziChunks), np.max(newAziChunks)])
            y_lims = np.around(y_lims / 180) * 180
            AziPlot = range(np.int(y_lims[0]), np.int(y_lims[1]), 2)
            gmodel = Model(ff.Fourier_expand, independent_vars=['azimu'])

            # fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(1, 5)
            fig = plt.figure()
            ax1 = fig.add_subplot(5, 1, 1)
            # plt.subplot(511)
            ax1.set_title('D-spacing')
            for j in range(peeks):
                ax1.scatter(newAziChunks, newd0[j], s=10)
                # plt.errorbar(newAziChunks,newd0[j],newd0Err[j], linestyle="None") #FIX ME: should have a different
                # colour for each peak in the subpattern
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot), param=dfour[j][0])
                # have to pass params but use of param will override
                ax1.plot(AziPlot, gmod_plot)
                # plt.plot(AziPlot, ff.Fourier_expand(np.array(AziPlot), dfour[j][0]), 'r-')

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
                ax2.plot(AziPlot, gmod_plot)  # , 'r-')
                # plt.plot(AziPlot, ff.Fourier_expand(np.array(AziPlot) * symm, hfour[j][0]), 'r-')

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
                # plt.plot(AziPlot, ff.Fourier_expand(np.array(AziPlot) * symm, np.array(wfour[j][0])), 'r-')

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
                # plt.plot(AziPlot, ff.Fourier_expand(np.array(AziPlot) * symm, np.array(pfour[j][0])), 'r-')

            for k in range(len(lenbg)):
                x_plt = len(lenbg)
                y_plt = len(lenbg) * 4 + k + 1
                ax5 = fig.add_subplot(5, x_plt, y_plt)
                ax5.set_title('Background')
                ax5.scatter(newAziChunks, newBGall[k], s=10)
                gmod_plot = gmodel.eval(params=master_params, azimu=np.array(AziPlot), param=bgfour[k][0])
                ax5.plot(AziPlot, gmod_plot)
                # plt.plot(AziPlot, ff.Fourier_expand(np.array(AziPlot), np.array(bgfour[k])), 'r-')

            plt.show()
            plt.close()

        # put fits into NewParams data structure.
        NewParams = ff.create_newparams(peeks, dfour, hfour, wfour, pfour, bgfour, orders['peak'])

        # Iterate over each fourier series in turn.

        # refine = 1  # Moved to FitPattern start
        # iterations = 1  Moved to FitPattern start
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

                    param_str = 'peak_' + str(j)
                    if 'profile_fixed' not in orders['peak'][k]:
                        p_fixed = None
                    else:
                        p_fixed = True
                    if 'symmetry' in orders['peak'][k]:
                        symm = orders['peak'][k]['symmetry']
                    else:
                        symm = 1
                    master_params, dfour, hfour, wfour, pfour = update_component_fits(master_params, intens, azimu,
                                                                                      twotheta, param_str,
                                                                                      [dfour, hfour, wfour, pfour],
                                                                                      conversion_factor, bgguess,
                                                                                      Guesses, p_fixed)

                # refine background
                param_str = 'bg_c'
                # print(param_str, lenbg, i)
                master_params = ff.unvary_params(master_params, param_str, '')
                master_params = ff.vary_params(master_params, param_str, '')
                out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), num_peaks=len(Guesses['peak']),
                                  nterms_back=len(bgguess), Conv=conversion_factor, fixed=pfixed, method=None,
                                  weights=None, params=master_params)
                master_params = out.params

                # Put fits into data structure.
                NewParams = ff.create_newparams(peeks, dfour, hfour, wfour, pfour, bgfour, orders['peak'])

                '''
                N.B. unindented one from original position
                if 'profile_fixed' in orders['peak'][0]:
                    refine_parms = 3
                else:
                    refine_parms = 4
    
                # loop over parameters
                for l in range(refine_parms):
    
                    Chng = ff.ParamFitArr(peeks, len(backg), 0)
                    Chng[0][k][l] = 1  # necessary?
                    # only 4 parameters describe each peak
                    # set parameter to change to length of array
                    # FIX ME: should this just be 1/0 as a switch or the order of the array - then we pass the
                    # order array around and match things up later?
                    if l == 0:
                        Chng[0][k][l] = orders['peak'][k]['d-space'] * 2 + 1
                        Parm = 'd-space'
                    elif l == 1:
                        Chng[0][k][l] = orders['peak'][k]['height'] * 2 + 1
                        Parm = 'height'
                    elif l == 2:
                        Chng[0][k][l] = orders['peak'][k]['width'] * 2 + 1
                        Parm = 'width'
                    elif l == 3:
                        Chng[0][k][l] = orders['peak'][k]['profile'] * 2 + 1
                        Parm = 'profile'
    
                    print(Chng)
                    NewParams, po, pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng,
                                                    NewParams, fixed=pfixed, Conv=conversion_factor)
                    # NewParams = ModParams
                    # print 'Parms out', NewParams, '\n'
                    # end of un-indent from original position   
                

                # refine background
                Chng = ff.ParamFitArr(peeks, len(backg), 0)
                for k in range(len(orders['background'])):
                    Chng[1][k] = orders['background'][k] * 2 + 1
                print(Chng)
                NewParams, po, pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng, NewParams,
                                                fixed=pfixed, Conv=conversion_factor)
                '''

    else:
        # FIX ME: This should load the saved lmfit Parameter class object. Need to initiate usage (save and load).
        # For now have propagated use of NewParams so re-instantiate the master_params object below, which is clunky.
        NewParams = PreviousParams
        pfixed = 'something'
        print('Somehow got to previousparams')
        # print(PreviousParams)
        # for ba in PreviousParams:
        #    print(PreviousParams[ba])

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
            dfour = PreviousParams['peak'][j]['d-space']
            comp = 'd'
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['d-space'],
                                               value=np.array(dfour))

            # initiate h
            hfour = PreviousParams['peak'][j]['height']
            comp = 'h'
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['height'], value=hfour)

            # initiate w
            wfour = PreviousParams['peak'][j]['width']
            comp = 'w'
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['width'], value=wfour)

            # initiate p
            pfour = PreviousParams['peak'][j]['profile']
            comp = 'p'
            master_params = ff.initiate_params(master_params, param_str, comp, orders['peak'][j]['profile'],
                                               value=pfour)

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
        bgguess = singleBackg
        bgguess[0][0] = 5

        # def previous_params_to_params(PreviousParams):

        # FIX ME: need to confirm the number of parameters matches the orders of the fits.

    # FIX ME: Need to sort out use of lmfit parameters and NewParams, below will not work currently for loaded
    # parameters as master_params not created from NewParams yet. If possible use the json import/export within lmfit.

    # Refit through full equation with all data for d,h,w,bg independently
    print('\nFinal fit solving for all parms...\n')

    '''
    Chng = [[], []]
    for x in range(len(orders['peak'])):
        if 'profile_fixed' in orders['peak'][x]:
            Chng[0].append([orders['peak'][x]['d-space'] * 2 + 1, orders['peak'][x]['height'] * 2 + 1,
                            orders['peak'][x]['width'] * 2 + 1, 0])
        else:
            Chng[0].append([orders['peak'][x]['d-space'] * 2 + 1, orders['peak'][x]['height'] * 2 + 1,
                            orders['peak'][x]['width'] * 2 + 1, orders['peak'][x]['profile'] * 2 + 1])
    for x in range(len(orders['background'])):
        Chng[1].append(orders['background'][x] * 2 + 1)
    print(Chng)
    NewParams, po, pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng, NewParams,
                                    fixed=pfixed, Conv=conversion_factor)
    '''
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
            new_str = 'peak_' + str(x) + '_' + p
            str_keys = [key for key, val in master_params.items() if new_str in key]
            for pstr in str_keys:
                master_params[pstr].set(vary=False)

    # FIX ME: need to sort parameters, bgguess doesn't exist if load previous fit data
    out = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), num_peaks=num_peaks,
                      nterms_back=len(bgguess), Conv=conversion_factor, fixed=pfixed, method=None,
                      weights=None, params=master_params)
    master_params = out.params
    print('\nFinal Coefficients\n')
    # print(out.fit_report(min_correl=1))
    print(out.fit_report())
    # print('Final Coefficients\n', NewParams)

    # Write master_params to NewParams dict object
    NewParams = ff.params_to_newparams(master_params, num_peaks=num_peaks, num_bg=len(bgguess),
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
        fullfit_intens = gmodel.eval(params=master_params, twotheta=twotheta.flatten(), azi=azimu.flatten(),
                                     num_peaks=num_peaks, nterms_back=len(bgguess), Conv=conversion_factor,
                                     fixed=pfixed)
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
        ax = fig.add_subplot(1, 3, 1)
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
        if SaveFit == 1:
            filename = ''
            for x in range(len(orders['peak'])):
                if 'phase' in orders['peak'][x]:
                    filename = filename + orders['peak'][x]['phase']
                    if 'hkl' in orders['peak'][x]:
                        filename = filename + str(orders['peak'][x]['hkl'])
                if x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
                    filename = filename + '_'
            if filename == '':
                filename = 'Fit2Peak'

            i = 1  # fix me: this doesnt make a new figure it over writes the old one.
            while os.path.exists('{}_{:d}.png'.format(filename, i)):
                i += 1
            if i == 1:
                plt.savefig('{}.png'.format(filename))
            else:
                plt.savefig('{}_{:d}.png'.format(filename, i))

        if view == 1 or debug:
            plt.show()

        plt.close()
        print('Done with Figure')

    '''
    collated = {"d-space":    fin_d,
                "height":     fin_h,
                "width":      fin_w,
                "background": backg
    }
    '''

    tth_min = twotheta.min()
    tth_max = twotheta.max()
    d_min = dspace.min()
    d_max = dspace.max()
    extent = [[tth_min, tth_max], [d_min, d_max]]  # FIX ME: This is taking the maximum and the minimum of the data
    # not the 'range' itself.

    NewParams.update({'range': extent})

    '''
    FitStats = {"degree-of-freedom": deg_freedom,
                "n-points": n_points,
                "ssd": SSD,
                "varience": SSD_var,
                "ChiSq": ChiSq
                }
    '''

    FitStats = {"degree-of-freedom": deg_freedom,
                "n-points": n_points,
                "ChiSq": ChiSq
                }

    NewParams.update({'FitProperties': FitStats})

    # FIX ME: Need to put master_params into NewParams to return!

    return NewParams
