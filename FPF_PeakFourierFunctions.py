# /usr/bin/python

import numpy as np
import numpy.ma as ma
import copy, os, sys
from PIL import Image
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.path as mlp
from scipy.optimize import curve_fit
np.set_printoptions(threshold='nan')


# Fitting functions required for fitting peaks to the data. 
# Called by FPF_XRD_FitSubpattern.


# Known limitations/bugs
# FIX ME: need hard limits to gauss - lorentz peak profile ratios



#Gausian shape
def GaussianPeak(twotheta, tth0, W_all, H_all):
    W_all = W_all/np.sqrt(np.log(4))
    GaussPeak = H_all*np.exp((-(twotheta-tth0)**2)/(2*W_all**2))
    return GaussPeak

# Lorentz shape.
def LorentzianPeak(twotheta, tth0, W_all, H_all):
    LPeak = H_all*W_all**2/((twotheta-tth0)**2 + W_all**2)
    return LPeak


#Pseudo-Voigt.
def PseudoVoigtPeak(twotheta, tth0, W_all, H_all, LGratio):

    #FIX ME: how should the ratio be fixed?
    # The ratio is between 1 and 0. it has to be limited.
    # But need options for fixed ratio and not fitting them.
    # How do I do this?
    #
    # I have defined a new variable called profile_fixed which fixes the value if present

    if np.any(LGratio<0) or np.any(LGratio>1):
        print 'The ratio is out of range'
        #stop

    PVpeak = LGratio*GaussianPeak(twotheta, tth0, W_all, H_all) + (1-LGratio)*LorentzianPeak(twotheta, tth0, W_all, H_all)
    return PVpeak


def Fourier_order(params):
    # Given list of Fourier coefficients retrun order (n) of the Fourier series.
    
    if isinstance(params,(list,)):
        order = (len(params)-1)/2
    elif isinstance(params,(float,)):
        order = (np.size(params)-1)/2
    else:
        print 'Unknown type'
        stop

    return order


# fourier expansion function
def Fourier_expand(azimu, *param):

    param=np.array(param)
    if (len(param.shape) > 1):
        if np.any(np.array(param.shape) > 1) :
            param = np.squeeze(param)
            # print param,'1'
        elif np.all(np.array(param.shape) == 1):
            param = np.squeeze(param)
            param = np.array([param],float)
            # print param, '2'

    param=tuple(param)
    out = np.ones(azimu.shape)
    out[:] = param[0]
    # essentially d_0, h_0 or w_0

    if len(param)>1:
        for i in xrange(1, ((len(param)-1)/2)+1): # len(param)-1 should never be odd because of initial a_0 parameter
            out = out + param[(2*i)-1]*np.sin(np.deg2rad(azimu)*i) + param[2*i]*np.cos(np.deg2rad(azimu)*i) #single col array
    #else:
      #  azimu = np.array(azimu)
       # fout = np.ones(azimu.shape)
        #fout[:] = out
        #out = fout*param

    return out



def Fourier_fit(azimu,ydata,terms,param=None,errs=1):

    #print type(terms)
    if(type(terms)==list):
        terms = terms[0]
    #print terms
    #print param
    if param:
        param=param
    else:
        param = [0 for i in range((2*terms+1))]
    param=tuple(param)
    popt,pcurv = curve_fit(Fourier_expand,azimu,ydata,p0=param,sigma=errs)


    return popt,pcurv



# fourier expansion function
def Fourier_backgrnd(azimutheta, param):

    azimu,twotheta = azimutheta
    twothetaprime = twotheta-twotheta.min()
    backg=param
    bg_all = np.zeros(twotheta.shape)
    
    for i in xrange(len(backg)):
        nterms = len(backg[i])
        out = backg[i][0]
        for j in xrange(1, ((nterms-1)/2)+1): 
            out = out + backg[i][(2*j)-1]*np.sin(np.deg2rad(azimu)*j) + backg[i][2*j]*np.cos(np.deg2rad(azimu)*j)
            # print j, 'j'
        bg_all = bg_all + (out*(twothetaprime**float(i)))

    return bg_all


def FitPenaltyFunction(Azimuths, FourierValues, valid_range,Weight=1000000):

    # print Azimuths
    # print FourierValues
    # print valid_range

    # print np.size(FourierValues)
    if np.size(FourierValues) > 1:
        Vals = Fourier_expand(Azimuths, FourierValues)
    else:
        Vals = FourierValues


    Penalise = 0
    if np.min(Vals) < valid_range[0]:
        Penalise = valid_range[0] - np.min(Vals)

    if np.max(Vals) > valid_range[1]:
        Penalise = np.max(Vals) - valid_range[0]

    Penalise = Penalise**4*Weight + 1
    return Penalise


def singleInt(twothetaWL,*sparms):

    #global wavelength

    #split two theta and wavelength (if present)
    if len(twothetaWL) == 5:
        twotheta,azimu,lenbg,wavelength,profile = twothetaWL
    elif len(twothetaWL) == 4:
        twotheta,azimu,lenbg,wavelength = twothetaWL
        # twotheta = twothetaWL[0]
        # wavelength = twothetaWL[1]
    else:
        twotheta,azimu,lenbg = twothetaWL
        # twotheta = twothetaWL

    sparms=np.array(sparms)
    #print params.shape
    if (len(sparms.shape) > 1):
        if np.any(np.array(sparms.shape) > 1) :
            sparms = np.squeeze(sparms)
            #print params,'1'
        elif np.all(np.array(sparms.shape) == 1):
            sparms = np.squeeze(sparms)
            sparms = np.array([sparms],float)
            #print params, '2'
    sparms=tuple(sparms)

    ##fit d,w,h as single number with no Fourier terms
    #print 'sparams', sparms
    d_0=sparms[0]
    H_all=sparms[1]
    W_all=sparms[2]
    tmpbackg = np.array(sparms[3:3+lenbg[0]])
    backg = []
    for b in xrange(len(lenbg)):
        if b == 0:
            start = 0
            end = lenbg[b]
        else: 
            start = lenbg[b-1]
            end = lenbg[b]+start 
        backg.append(tmpbackg[start:end])
    if len(sparms) == lenbg+4:
        profile = sparms[lenbg[0]+3]
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    bg_all = Fourier_backgrnd((azimu,twotheta),backg)
    #print bg_all
    #print backg
    # bg either single num. or size of azimu array
    # print bg_all, 'bg_all'
    Int = PseudoVoigtPeak(twotheta, tth0, W_all, H_all, profile) + bg_all

    dpen  = FitPenaltyFunction(azimu, tth0, [np.min(twotheta), np.max(twotheta)])
    hpen  = FitPenaltyFunction(azimu, H_all, [0, np.max(H_all)])
    wpen  = FitPenaltyFunction(azimu, W_all, [0, (np.max(twotheta)-np.min(twotheta))])
    ppen  = FitPenaltyFunction(azimu, profile, [0, 1])
    bgpen = FitPenaltyFunction(azimu, backg, [0, np.inf])


    Int = Int * dpen * hpen * wpen * ppen * bgpen


    return Int



def singleFit(intens,twotheta,azimu,dspace,d0s,heights,widths,profile,fixed,wavelength,bg=None):


    '''
    All expansion parameters fitted
    '''
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')

    lenbg=[]
    for val in bg:
        lenbg.append(len(val))
    lenbg=np.array(lenbg)

    if fixed == 0:
        allparms = np.concatenate((d0s,heights,widths,[item for sublist in bg for item in sublist], profile), axis=None)
        #allparms.append(profile)
        #print allparms.size
        dataForFit = (twotheta, azimu,lenbg,wavelength)
        #print np.size(allparms)
        param_bounds=([-np.inf]*np.size(allparms),[np.inf]*np.size(allparms))
        param_bounds[0][-1] = 0  #force min Pseudo voigt ratio to be greater than zero
        param_bounds[1][-1] = 1  #force max Pseudo voigt ratio to be less than 1
        param_bounds[0][2] = 0   #force min height to be greater than zero
        param_bounds[0][3] = 0   #force min width to be greater than zero

    else:
        allparms = np.concatenate((d0s,heights,widths,[item for sublist in bg for item in sublist]), axis=None)
        dataForFit = (twotheta, azimu,lenbg,wavelength,profile)
        param_bounds=(-np.inf*np.ones(allparms.shape),np.inf*np.ones(allparms.shape))

    #print 'Curve_fit inputs:'
    #print allparms
    #print dataForFit
    #print param_bounds
    
    try:
        popt,pcurv = curve_fit(singleInt, dataForFit, intens, p0=allparms, bounds=param_bounds)

    except RuntimeError:
        print("Error - curve_fit failed")
        popt = np.zeros(allparms.shape)
        print popt
        print np.size(popt)
        pcurv = np.ones((popt.shape[0], popt.shape[0])) *5.

        print popt
        print pcurv

    return popt,pcurv







def ParamFit(param_change,param_fit,intens,twotheta,azimu,dspace,d0s,heights,widths,profiles,wavelength,bg=None):

    def Paramchange(fullarray,*fitparms):

        #expand input array
        twothet,azi,lenbg,tmpbackg = fullarray

        #get arrays of constants (all but background)
        D_all=d0s
        H_all=heights
        W_all=widths
        P_all=profiles

        # Organise the parameter to be fitted.
        fitparms=np.array(fitparms)
        #print params.shape
        if (len(fitparms.shape) > 1):
            if np.any(np.array(fitparms.shape) > 1) :
                fitparms = np.squeeze(fitparms)
            elif np.all(np.array(fitparms.shape) == 1):
                fitparms = np.squeeze(fitparms)
                fitparms = np.array([fitparms],float)   
        fitparms=tuple(fitparms)
        #d_0 = Fourier_expand(azi, fitparms)
        fit_all = Fourier_expand(azi, fitparms)
        # print d_0

        #Copy expanded values back into the correct place.
        #print param_change
        if param_change == 'd-space':
            D_all=np.array(fit_all)
        elif param_change == 'height':
            H_all=np.array(fit_all)
        elif param_change == 'width':
            W_all=np.array(fit_all)
        elif param_change == 'profile':
            P_all=np.array(fit_all)
        elif param_change == 'background':
            tmpbackg=np.array(fit_all)
        else:
            print 'Unknown!!!'
            stop

        #make backgorund fourier from array (either original or fitted array)
        #FIX ME: why is this not done outside this fitting routine? 
        backg = []
        for b in xrange(len(lenbg)):
            if b == 0:
                start = 0
                end = lenbg[b]
            else: 
                start = lenbg[b-1]
                end = lenbg[b]+start 
            backg.append(tmpbackg[start:end])
        
        bg_all = Fourier_backgrnd((azi,twothet),backg)

        #print wavelength
        tth_all = 2.*np.degrees(np.arcsin(wavelength/2/D_all))

        # print stop
        Int = PseudoVoigtPeak(twotheta, tth_all, W_all, H_all, P_all) + bg_all

        return Int

    '''
    Only d_0 fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
    heights=heights

    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    flatbg = np.array(np.concatenate(([item for sublist in bg for item in sublist]), axis=None))
    lenbg=[]
    for val in bg:
        lenbg.append(len(val))
    lenbg=np.array(lenbg)
    
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    
    popt,pcurv = curve_fit(Paramchange,(twotheta,azimu,lenbg,flatbg),intens,p0=param_fit)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv






def Allchange(fullarray,*allparms):

    ##fit d,w,h as single number with no Fourier terms

    #sortout allparms
    allparms=np.array(allparms)
    if (len(allparms.shape) > 1):
        if np.any(np.array(allparms.shape) > 1) :
            allparms = np.squeeze(allparms)
        elif np.all(np.array(allparms.shape) == 1):
            allparms = np.squeeze(allparms)
            allparms = np.array([allparms],float)  

    #determine if profile is beting fitted for and organise the fitted numbers
    if len(fullarray) == 6: #profile is a fitted parameter
        twothet,azi,lenbg,parmnums,wavelength,symm = fullarray
    elif len(fullarray) == 7: 
        twothet,azi,lenbg,parmnums,wavelength,symm,profiles = fullarray

    # print twothet.shape,azi.shape,parmnums.shape
    start = 0
    starth = parmnums[0]
    startw = parmnums[0:2].sum()
    if len(fullarray) == 6:
        startp = parmnums[0:3].sum()
    endp = parmnums.sum()

    d0s = allparms[0:starth]
    heights = allparms[starth:startw]
    if len(fullarray) == 6:
        widths = allparms[startw:startp]
        profiles = allparms[startp:endp]
    else:
        widths = allparms[startw:endp]


    print profiles
    

    tmpbackg = allparms[endp:]
    backg = []
    for b in xrange(len(lenbg)):
        if b == 0:
            startbg = 0
            endbg = lenbg[b]
        else: 
            startbg = lenbg[b-1]
            endbg = lenbg[b]+startbg
        backg.append(tmpbackg[startbg:endbg])

    # print d0s,heights,widths
    d_0 = Fourier_expand(azi, d0s)
    H_all = Fourier_expand(azi*symm, heights)
    W_all = Fourier_expand(azi*symm, widths)
    P_all = Fourier_expand(azi*symm, profiles)
    print np.max(d_0), np.min(d_0)
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
    bg_all = Fourier_backgrnd((azi,twothet),backg)
    # print d_0,H_all,W_all,tth0,bg

    newall =  PseudoVoigtPeak(twothet, tth0, W_all, H_all, P_all) + bg_all

    #print newall.shape

    return newall





def Allfit(intens,twotheta,azimu,d0s,heights,widths,profiles,wavelength,bg=None,symm=None, fix=None):

    print fix

    '''
    All expansion parameters fitted
    '''

    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    #allparms = np.concatenate((d0s,heights,widths,profiles,[item for sublist in bg for item in sublist]), axis=None)
    lenbg=[]
    for val in bg:
        lenbg.append(len(val))
    lenbg=np.array(lenbg)
    print 'profile', profiles
    inp = []
    if fix == None: #if profile not fixed fit for profile otherwise it is a constant.
        parmnums = np.array([len(d0s),len(heights),len(widths),len(profiles)])
        inp      = (twotheta,azimu,lenbg,parmnums,wavelength,symm)
        allparms = np.concatenate((d0s,heights,widths,profiles,[item for sublist in bg for item in sublist]), axis=None)
    else:
        parmnums = np.array([len(d0s),len(heights),len(widths)])
        inp      = (twotheta,azimu,lenbg,parmnums,wavelength,symm,profiles)
        allparms = np.concatenate((d0s,heights,widths,[item for sublist in bg for item in sublist]), axis=None)
    
    # print 'parmnums', parmnums

    # #check for bg here, won't work if try to pass to fit
    # if not bg:
    #     bg = backg
    # allparms = np.concatenate((d0s,heights,widths,profiles,[item for sublist in bg for item in sublist]), axis=None)
    # lenbg=[]
    # for val in bg:
    #     lenbg.append(len(val))
    # lenbg=np.array(lenbg)

    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')

    popt,pcurv = curve_fit(Allchange,inp,intens,p0=allparms, maxfev = 12000)
    #popt,pcurv = curve_fit(Allchange,(twotheta,azimu,lenbg,parmnums,wavelength,symm),intens,p0=allparms, maxfev = 12000)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv





def update_backgrnd(params,lenbg,lenparams):

    backg = []
    for b in xrange(len(lenbg)):
        if b == 0:
            startbg = lenparams
            endbg = startbg+lenbg[b]
        else: 
            startbg = lenparams+lenbg[b-1]
            endbg = lenbg[b]+startbg 
        backg.append(list(params[startbg:endbg]))
    return backg
