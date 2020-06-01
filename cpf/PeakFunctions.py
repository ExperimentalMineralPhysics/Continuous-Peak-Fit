#!/usr/bin/env python

__all__ = ['create_newparams', 'params_to_newparams', 'ParamFitArr', 'PseudoVoigtPeak', 'gather_paramerrs_to_list',
           'gather_params_to_list', 'gather_params_from_dict', 'initiate_params', 'PeaksModel', 'vary_params',
           'unvary_params', 'FitModel', 'Fourier_expand', 'Fourier_fit', 'ChangeParams',
           'Fourier_backgrnd']

import numpy as np
import numpy.ma as ma
import sys
from scipy.optimize import curve_fit
from lmfit import minimize, report_fit, Model

np.set_printoptions(threshold=sys.maxsize)

# Fitting functions required for fitting peaks to the data. 
# Called by FPF_XRD_FitSubpattern.


# Known limitations/bugs
# FIX ME: need hard limits to gauss - lorentz peak profile ratios


def create_newparams(num_peaks, dfour, hfour, wfour, pfour, bgfour, order_peak):
    """
    Create the NewParams dictionary from results
    :param num_peaks: number of peaks int
    :param dfour: list of coefficients per peak float
    :param hfour: list of coefficients per peak float
    :param wfour: list of coefficients per peak float
    :param pfour: list of coefficients per peak float
    :param bgfour:
    :param order_peak:
    :return: NewParams dictionary
    """
    peaks = []
    for j in range(num_peaks):
        if 'symmetry' in order_peak[j]:  # orders['peak'][j]:  # FIX ME: the symmetry part of this is a horrible way
            # of doing it.
            peaks.append({"d-space": dfour[j][0],
                          "height": hfour[j][0],
                          "width": wfour[j][0],
                          "profile": pfour[j][0],
                          "symmetry": order_peak[j]['symmetry']})
                            # "symmetry": orders['peak'][j]['symmetry']})
        else:
            peaks.append({"d-space": dfour[j][0],
                          "height": hfour[j][0],
                          "width": wfour[j][0],
                          "profile": pfour[j][0]})
    NewParams = {"background": bgfour, "peak": peaks}

    return NewParams


def params_to_newparams(params, num_peaks, num_bg, order_peak):
    """
    Gather Parameters class to NewParams dictionary format
    :param params: lmfit Parameters class object
    :param num_peaks: number of peaks
    :param num_bg: number of background terms
    :return: NewParams dictionary object
    """

    peaks = []
    for i in range(num_peaks):
        new_str = 'peak_'+str(i)+'_d'
        dspace = gather_paramerrs_to_list(params, new_str)
        new_str = 'peak_'+str(i)+'_h'
        hspace = gather_paramerrs_to_list(params, new_str)
        new_str = 'peak_'+str(i)+'_w'
        wspace = gather_paramerrs_to_list(params, new_str)
        new_str = 'peak_'+str(i)+'_p'
        pspace = gather_paramerrs_to_list(params, new_str)

        if 'symmetry' in order_peak[i]:  # orders['peak'][j]:  # FIX ME: the symmetry part of this is a horrible
            # way of doing it.
            peaks.append({"d-space":     dspace[0],
                          "d-space_err": dspace[1],
                          "height":      hspace[0],
                          "height_err":  hspace[1],
                          "width":       wspace[0],
                          "width_err":   wspace[1],
                          "profile":     pspace[0],
                          "profile_err": pspace[1],
                          "symmetry":    order_peak[i]['symmetry']})
        else:
            peaks.append({"d-space": dspace[0],
                          "d-space_err": dspace[1],
                          "height":      hspace[0],
                          "height_err":  hspace[1],
                          "width":       wspace[0],
                          "width_err":   wspace[1],
                          "profile":     pspace[0],
                          "profile_err": pspace[1]})
    # Get background parameters
    bgspace = []
    bgspace_err = []
    for b in range(num_bg):
        new_str = 'bg_c'+str(b)+'_f'
        bgspc = gather_paramerrs_to_list(params, new_str)
        bgspace.append(bgspc[0])
        bgspace_err.append(bgspc[1])
        # print('background', ff.gather_params_to_list(params, new_str))
        # bgspace.append(ff.gather_params_to_list(params, new_str))

    # NewParams = {"background": bgspace, "peak": peaks}
    NewParams = {"background": bgspace, "background_err": bgspace_err, "peak": peaks}

    return NewParams


def flatten(li):
    return sum(([x] if not isinstance(x, list) else flatten(x)
                for x in li), [])


def ParamFitArr(numP, LenBG, val=None, p=None, b=None, orda=None):
    # Makes array which represents all the fourier series in the fit.
    # Arr{0} is a list of 4 element lists - one 4 element list for each peak.
    # Arr[0] is the background representative.

    if val:
        v = val
    else:
        v = 0

    if not p:
        p = [v, v, v, v]

    Arr = [[], []]
    for y in range(numP):
        Arr[0].append(p[:])
    for y in range(LenBG):
        Arr[1].extend([v][:])
    return Arr


# Gausian shape
def GaussianPeak(twotheta, tth0, W_all, H_all):
    W_all = W_all / np.sqrt(np.log(4))
    GaussPeak = H_all * np.exp((-(twotheta - tth0) ** 2) / (2 * W_all ** 2))
    return GaussPeak


# Lorentz shape.
def LorentzianPeak(twotheta, tth0, W_all, H_all):
    LPeak = H_all * W_all ** 2 / ((twotheta - tth0) ** 2 + W_all ** 2)
    return LPeak


# Pseudo-Voigt.
def PseudoVoigtPeak(twotheta, tth0, W_all, H_all, LGratio):
    # FIX ME: how should the ratio be fixed?
    # The ratio is between 1 and 0. it has to be limited.
    # But need options for fixed ratio and not fitting them.
    # How do I do this?
    #
    # I have defined a new variable called profile_fixed which fixes the value if present
    # if np.any(LGratio<0) or np.any(LGratio>1):
    # print 'The ratio is out of range'
    # stop

    PVpeak = LGratio * GaussianPeak(twotheta, tth0, W_all, H_all) + (1 - LGratio) * LorentzianPeak(twotheta, tth0,
                                                                                                   W_all, H_all)
    return PVpeak


def gather_paramerrs_to_list(param, param_str, comp=None, nterms=None):  # fix me: Doesn't run without Comp=None.
    # Maybe it was never validated before SAH ran it.
    """
    Make a nested list of parameters and errors from a lmfit Parameters class
    :param param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: nested list of [parameters, errors] in alphanumerical order
    """

    if comp:
        new_str = param_str + '_' + comp
    else:
        new_str = param_str
    param_list = []
    err_list = []
    total_params = []
    str_keys = [key for key, val in param.items() if new_str in key]
    # print('yes', new_str, str_keys)
    for i in range(len(str_keys)):
        param_list.append(param[new_str + str(i)].value)
        err_list.append(param[new_str + str(i)].stderr)
    total_params.append(param_list)
    total_params.append(err_list)
    return total_params


def gather_params_to_list(param, param_str, comp=None):
    """
    Make a list of parameters from a lmfit Parameters class
    :param param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: list of parameters in alphanumerical order
    """

    if comp:
        new_str = param_str + '_' + comp
    else:
        new_str = param_str
    param_list = []
    str_keys = [key for key, val in param.items() if new_str in key]
    # print('yes', new_str, str_keys)
    for i in range(len(str_keys)):
        param_list.append(param[new_str + str(i)].value)
    return param_list


def gather_params_from_dict(param, param_str, comp):
    """
    Make a list of parameters from a dictionary
    :param param: dict with multiple coefficients per component
    :param param_str: base string to select parameters
    :param comp: component to add to base string to select parameters
    :return: list of parameters in alphanumerical order
    """

    if comp:
        new_str = param_str + '_' + comp
    else:
        new_str = param_str
    param_list = []
    str_keys = [key for key, val in param.items() if new_str in key]
    for i in range(len(str_keys)):
        param_list.append(param[new_str + str(i)])
    return param_list


def initiate_params(param, param_str, comp, nterms=None, max=None, min=None, expr=None, value=1.):
    """
    Create all required coefficients for no. terms
    :param param: lmfit Parameter class (can be empty)
    :param param_str: base string to create lmfit parameter
    :param comp: string to add to base string to create lmfit parameter
    :param nterms: number of terms required
    :return:
    """
    
    value = np.array(value)  # force input to be array so that can be iterated over

    if nterms is None:
        nterms = int((len(value) - 1)/2)

    for t in range(2 * nterms + 1):
        if value.size == 1:
            ind = 0
        else:
            ind = t
        param.add(param_str + '_' + comp + str(t), value.item(ind), max=max, min=min, expr=expr)
        # param.add(param_str + '_' + comp + str(t), value, max=max, min=min, expr=expr)
    return param


def unvary_params(param, param_str, comp):
    """
    Set all but selected parameters to vary=False
    :param param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + '_' + comp
    else:
        new_str = param_str
    str_keys = [key for key, val in param.items() if new_str in key]
    # print(str_keys)
    for p in param:
        # print(p)
        if p not in str_keys:
            # print('yes')
            param[p].set(vary=False)
    return param


def vary_params(param, param_str, comp):
    """
    Set selected parameters to vary=True
    :param param: lmfit Parameter class
    :param param_str: base string to select with
    :param comp: string to add to base string to select with
    :return: updated lmfit Parameter class
    """
    if comp:
        new_str = param_str + '_' + comp
    else:
        new_str = param_str
    str_keys = [key for key, val in param.items() if new_str in key]
    # print(str_keys)
    for p in str_keys:
        # print('yes')
        param[p].set(vary=True)
    return param


#def PeaksModel(twotheta, azi, num_peaks=None, nterms_back=None, Conv=None, **params):
def PeaksModel(twotheta, azi, Conv=None, **params):
    """Full model of intensities at twotheta and azi given input parameters
    :param twotheta: arr values float
    :param azi: arr values float
    :param num_peaks: total number of peaks int
    :param nterms_back: total number of polynomial expansion components for the background int
    :param Conv: inputs for the conversion call dict
    :param PenCalc: penalise background or not int
    :param params: lmfit parameter class dict
    :return lmfit model fit result
    """

    # N.B. params now doesn't persist as a parameter class, merely a dictionary, so e.g. call key/value pairs as
    # normal not with '.value'

    # recreate backg array to pass to Fourier_backgrnd
    back_keys = [key for key, val in params.items() if 'bg_c' in key]
    nterm_fouriers = []
    
    i=0
    backg=[]
    while sum(nterm_fouriers) < len(back_keys):
        
        f = sum('bg_c' + str(i) in L for L in back_keys)
        nterm_fouriers.append(f)
        
        fbg = []
        for k in range(f):
            fbg.append(params['bg_c' + str(i) + '_f' + str(k)])
        backg.append(np.array(fbg))   
        i = i+1
    
    # for i in range(int(nterms_back)):
    #     f = sum('bg_c' + str(i) in L for L in back_keys)
    #     nterm_fouriers.append(f)
    # backg = []
    # for j in range(int(nterms_back)):
    #     fbg = []
    #     for k in range(nterm_fouriers[j]):
    #         fbg.append(params['bg_c' + str(j) + '_f' + str(k)])
    #     backg.append(np.array(fbg))
    # print(backg)
    
    # for future logging
    # print('recreate backg to pass to fourier_backg\n', 'back_keys is ', back_keys)
    # print('backg is ', backg)

    I = Fourier_backgrnd((azi, twotheta), backg)

    
    peak_keys = [key for key, val in params.items() if 'peak' in key]
    num_peaks = 0
    nterm_peaks = 0
    while nterm_peaks < len(peak_keys):
        f = sum('peak_' + str(num_peaks) in L for L in peak_keys)
        nterm_peaks = nterm_peaks + f
        num_peaks = num_peaks+1

    # loop over the number of peaks
    Ipeak = []
    for a in range(int(num_peaks)):

        if 'peak_' + str(a) + '_s0' in params:
            symm = params['peak_' + str(a) + '_s0']
        else:
            symm = 1

        # May need to fix for mutliple Fourier coefficients per component
        param_str = 'peak_' + str(a)
        comp = 'd'
        parms = gather_params_from_dict(params, param_str, comp)
        # print(parms)
        Dall = Fourier_expand(azi, param=parms)
        comp = 'h'
        parms = gather_params_from_dict(params, param_str, comp)
        Hall = Fourier_expand(azi * symm, param=parms)
        comp = 'w'
        parms = gather_params_from_dict(params, param_str, comp)
        Wall = Fourier_expand(azi * symm, param=parms)
        comp = 'p'
        parms = gather_params_from_dict(params, param_str, comp)
        Pall = Fourier_expand(azi * symm, param=parms)

        # conversion
        TTHall = CentroidConversion(Conv, Dall, azi)
        Ipeak.append(PseudoVoigtPeak(twotheta, TTHall, Wall, Hall, Pall))
        I = I + Ipeak[a]

    return I


def ChangeParams(constants, *fitparam):
    print('In changeparam', len(constants))
    twotheta = constants[0]
    azi = constants[1]
    ChangeArray = constants[2]
    Shapes = constants[3]
    Conv = constants[4]
    symm = constants[5]
    fixed = constants[6]

    print('changeparams:')
    print(twotheta, azi)
    print(ChangeArray)
    print(Shapes)
    print(Conv)
    print(symm)
    print(fixed)
    print(fitparam)
    Shapes = GuessApply(ChangeArray, Shapes, fitparam)
    # print(stop)
    print('In changes', Shapes)
    I, pen = PeaksModel_old(twotheta, azi, Shapes, Conv, symm, PenCalc=1)

    return I  # *pen


def FitModel(Intfit, twotheta, azimu, params, num_peaks, nterms_back, Conv=None, fixed=None,
             method='leastsq', weights=None):
    # ChangeArray, Shapes, Conv=None, symm=None, fixed=None, method=None,
    # weights=None, bounds=None):
    """Initiate model of intensities at twotheta and azi given input parameters and fit
    :param Intfit: intensity values to fit arr
    :param twotheta: twotheta values arr
    :param azimu: azimu values arr
    :param params: lmfit Parameter class dict of parameters to fit
    :param num_peaks: total number of peaks in model int
    :param nterms_back: total number of polynomial component terms in background int
    :param Conv: inputs for the conversion call dict
    :param fixed: unsure
    :param method: lmfit choice of fitting method e.g. least squares
    :param weights: errors on intensity values arr of size Intfit
    :return: lmfit model result
    """

    # bounds passed but not used

    # p0array = GuessGet(ChangeArray, Shapes)
    # print('here', Intfit, twotheta, azimu, ChangeArray, Shapes)
    # params.pretty_print()
    gmodel = Model(PeaksModel, independent_vars=['twotheta', 'azi'])
    # print('parameter names: {}'.format(gmodel.param_names))
    # print('independent variables: {}'.format(gmodel.independent_vars))

    # minner = Minimizer(fcn2min, params, fcn_args=(x, data))
    # result = minner.minimize()
    # out = gmodel.fit(Intfit, params, twotheta=twotheta, azi=azimu, num_peaks=num_peaks, nterms_back=nterms_back,
    #                  Conv=Conv, fixed=fixed, nan_policy='propagate')
    out = gmodel.fit(Intfit, params, twotheta=twotheta, azi=azimu, Conv=Conv, nan_policy='propagate')
    # print(out.fit_report())
    # print(out.params)
    # out.params.pretty_print()
    # result = popt.minimize()
    return out




def CentroidConversion(Conv, args_in, azi):
    # Conv['DispersionType'] is the conversion type
    # Conv[...] are the values required for the conversion
    # FIX ME: these functions should be a sub function of the detector types. but need to work out how to do it.

    if Conv['DispersionType'] is None or Conv == 0:
        args_out = args_in

    elif Conv['DispersionType'] == 'AngleDispersive':
        # args_in are d-spacing
        # args_outis two thetas
        args_out = 2 * np.degrees(
            np.arcsin(Conv['conversion_constant'] / 2 / args_in))  # FIX ME: check this is correct!!!

    elif Conv['DispersionType'] == 'EnergyDispersive':

        args_out = []
        # for x in range(len(azi)):
        for x in range(azi.size):
            # determine which detector is being used
            # a=np.asscalar(np.where(Conv['azimuths'] == azi[x])[0])

            if azi.size == 1:
                a = np.array(np.where(Conv['azimuths'] == azi)[0][0])
                # a= np.asscalar(np.where(Conv['azimuths'] == azi)[0][0])
                args_out.append(
                    12.398 / (2 * args_in * np.sin(np.radians(Conv['calibs'].mcas[a].calibration.two_theta / 2))))
            else:
                #                print type(azi)
                #                print azi.mask[x]
                if not azi.mask[x]:
                    a = (np.where(Conv['azimuths'] == azi[x])[0][0])
                    args_out.append(12.398 / (
                            2 * args_in[x] * np.sin(np.radians(Conv['calibs'].mcas[a].calibration.two_theta / 2))))
                else:
                    args_out.append(0)
        if isinstance(azi, np.ma.MaskedArray):
            args_out = ma.array(args_out, mask=azi.mask)
    else:
        raise ValueError('Unrecognised conversion type')

    return args_out


def Fourier_order(params):
    # Given list of Fourier coefficients return order (n) of the Fourier series.

    if isinstance(params, (list,)):
        order = int((len(params) - 1) / 2)
    elif isinstance(params, (float,)):
        order = int((np.size(params) - 1) / 2)
    elif isinstance(params, (int,)):
        order = 0
    else:
        print(params)
        print(type(params))
        raise ValueError('Parameter list is not list or float.')

    return order


# fourier expansion function
def Fourier_expand(azimu, param=None, comp_str=None, **params):
    """Calculate Fourier expansion given input coefficients
    :param azimu: arr data array float
    :param param: list of Fourier coefficients float
    :param comp_str: str to determine which coefficients to use from params
    :param params: lmfit dict of coefficients as parameters
    :return:
    """

    if param:
        # FIX ME: Need to check the fourier is a licit length
        if not isinstance(param, np.float64):
            param = np.array(param)
            if len(param.shape) > 1:
                if np.any(np.array(param.shape) > 1):
                    param = np.squeeze(param)
                elif np.all(np.array(param.shape) == 1):
                    param = np.squeeze(param)
                    param = np.array([param], float)
        # param = tuple(param)

    else:
        # create relevant list of parameters from dict
        str_keys = [key for key, val in params.items() if comp_str in key]
        # print(str_keys, 'str_keys')
        param = []
        for j in range(len(str_keys)):
            if str.startswith(comp_str, 'bg'):
                param.append(params[comp_str + '_f' + str(j)])
            else:
                param.append(params[comp_str + str(j)])

    fout = np.ones(azimu.shape)
    # print(azimu, fout, 'azimu', type(param))
    if azimu.size == 1:  # this line is required to catch error when out is single number.
        try:
            fout = param[0]
        except IndexError:
            fout = param
    else:
        fout[:] = param[0]
    # essentially d_0, h_0 or w_0

    if not isinstance(param, np.float64) and np.size(param) > 1:
        for i in range(1, int((len(param) - 1) / 2) + 1):
            # len(param)-1 should never be odd because of initial a_0 parameter
            fout = fout + param[(2 * i) - 1] * np.sin(np.deg2rad(azimu) * i) + param[2 * i] * np.cos(
                np.deg2rad(azimu) * i)  # single col array
    # else:
    #  azimu = np.array(azimu)
    # fout = np.ones(azimu.shape)
    # fout[:] = out
    # out = fout*param

    return fout


def Fourier_fit(ydata, azimu, param, terms=None, errs=None, param_str='peak_0', symm=1, fit_method='leastsq'):
    """Fit the Fourier expansion to number of required terms
    :param ydata: Component data array to fit float
    :param azimu: data array float
    :param param: lmfit Parameter class dict
    :param terms: number of terms to fit int
    :param errs: Component data array errors of size ydata
    :param param_str: str start of parameter name
    :param symm: symmetry in azimuth of the fourier to be fit
    :param fit_method: lmfit method default 'leastsq'
    :return: lmfit Model result
    """

    # get NaN values.
    idx = np.isfinite(azimu) & np.isfinite(ydata) 
    if type(terms) == list:
        terms = terms[0]
    if param:
        param = param
    else:
        param = Parameters()
        for t in range(2 * terms + 1):
            params.add(param_str + '_' + str(t), 1.)  # need limits

    if errs is None or errs.all is None:
        errs = np.ones(ydata.shape)

    # param.pretty_print()
    fmodel = Model(Fourier_expand, independent_vars=['azimu'])
    # print('parameter names: {}'.format(fmodel.param_names))
    # print('independent variables: {}'.format(fmodel.independent_vars))
    # Attempt to mitigate failure of fit with weights containing 'None' values
    # Replace with nan and alter dtype from object to float64
    new_errs = errs[idx]
    new_errs[new_errs == None] = 1000 * new_errs[new_errs != None].max()
    new_errs = new_errs.astype('float64')
    out = fmodel.fit(ydata[idx], param, azimu=azimu[idx]*symm, method=fit_method, sigma=new_errs, comp_str=param_str,
                     nan_policy='propagate')

    return out


# fourier expansion function
def Fourier_backgrnd(azimutheta, param):
    """
    Calculate the Fourier expansion of the background terms
    :param azimutheta: list of arrays for azimuth and theta float
    :param param: list of input parameters
    :return: Fourier expansion result float arr
    """
    azimu, twotheta = azimutheta
    twothetaprime = twotheta - twotheta.min()
    backg = param
    bg_all = np.zeros(twotheta.shape)
    # print(backg, bg_all)

    # Not sure if this is correct, thought that if orders then array would be eg [5,0] and would want a fourier
    # expansion with 5 parms for offset and then no fourier expansion for slope i.e. 5 parms then 0.
    # But below would interpret this as effectively [0,0] instead.

    for i in range(len(backg)):
        try:
            nterms = len(backg[i])
        except TypeError:
            nterms = backg[i].size
        # print((nterms - 1) / 2)
        out = backg[i][0]
        for j in range(1, int((nterms - 1) / 2) + 1):
            out = out + backg[i][(2 * j) - 1] * np.sin(np.deg2rad(azimu) * j) + backg[i][2 * j] * np.cos(
                np.deg2rad(azimu) * j)
        bg_all = bg_all + (out * (twothetaprime ** float(i)))

    return bg_all
