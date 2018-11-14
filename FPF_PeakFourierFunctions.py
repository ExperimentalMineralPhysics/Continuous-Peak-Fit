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


### Inputs ###
# and setup


# Fitting functions #


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
    #print param
    out = np.ones(azimu.shape)
    out[:] = param[0]
    # essentially d_0, h_0 or w_0

    if len(param)>1:
        for i in xrange(1, ((len(param)-1)/2)+1): # len(param)-1 should never be odd because of initial a_0 parameter
            out = out + param[(2*i)-1]*np.sin(np.deg2rad(azimu))**i + param[2*i]*np.cos(np.deg2rad(azimu))**i #single col array
    #else:
      #  azimu = np.array(azimu)
       # fout = np.ones(azimu.shape)
        #fout[:] = out
        #out = fout*param

    return out



def Fourier_fit(azimu,ydata,terms,param=None,errs=1):


    #print param
    if param:
        param=param
    else:
        param = [0 for i in range((2*terms+1))]
    param=tuple(param)
    popt,pcurv = curve_fit(Fourier_expand,azimu,ydata,p0=param)#,sigma=errs)

    return popt,pcurv


##warning about number of free parms e.g. if n is num. chunks/2





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
            out = out + backg[i][(2*j)-1]*np.sin(np.deg2rad(azimu))**j + backg[i][2*j]*np.cos(np.deg2rad(azimu))**j
            # print j, 'j'
        bg_all = bg_all + (out*(twothetaprime**float(i)))

    return bg_all






def singleInt(twothetaWL,*sparms):

    #global wavelength

    #split two theta and wavelength (if present)
    if len(twothetaWL) == 4:
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
    #print params, 'params'
    d_0=sparms[0]
    H_all=sparms[1]
    W_all=sparms[2]
    tmpbackg = np.array(sparms[3:])
    backg = []
    for b in xrange(len(lenbg)):
        if b == 0:
            start = 0
            end = lenbg[b]
        else: 
            start = lenbg[b-1]
            end = lenbg[b]+start 
        backg.append(tmpbackg[start:end])
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    bg_all = Fourier_backgrnd((azimu,twotheta),backg)

    # bg either single num. or size of azimu array
    # print bg_all, 'bg_all'
    Int = (H_all*np.exp((-(twotheta-tth0)**2)/(2*W_all**2))) + bg_all

    return Int



def singleFit(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):


    '''
    All expansion parameters fitted
    '''
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    allparms = np.concatenate((d0s,heights,widths,[item for sublist in bg for item in sublist]), axis=None)
    # print allparms
    lenbg=[]
    for val in bg:
        lenbg.append(len(val))
    lenbg=np.array(lenbg)
    popt,pcurv = curve_fit(singleInt, (twotheta, azimu,lenbg,wavelength), intens, p0=allparms,maxfev=8000)

    return popt,pcurv






def bgfit(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):


    def back_change(fullarray,*bgparms):

        '''
        Change background only
        '''

        twothet,azi,lenbg,d_0,H_all,W_all = fullarray
        bgparms=np.array(bgparms)
        if (len(bgparms.shape) > 1):
            if np.any(np.array(bgparms.shape) > 1) :
                bgparms = np.squeeze(bgparms)
            elif np.all(np.array(bgparms.shape) == 1):
                bgparms = np.squeeze(bgparms)
                bgparms = np.array([bgparms],float)   
        bgparms=tuple(bgparms)

        tmpbackg = bgparms
        backg = []
        for b in xrange(len(lenbg)):
            if b == 0:
                start = 0
                end = lenbg[b]
            else: 
                start = lenbg[b-1]
                end = lenbg[b]+start 
            backg.append(tmpbackg[start:end])

        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

        bg_all = Fourier_backgrnd((azi,twothet),backg)

        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg_all
        
        return Int


    heights=heights
    if not bg:
        bg = backg
    allparms = np.concatenate(([item for sublist in bg for item in sublist]), axis=None)
    lenbg=[]
    for val in bg:
        lenbg.append(len(val))
    lenbg=np.array(lenbg)
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(back_change,(twotheta,azimu,lenbg,d0s,heights,widths),intens,p0=allparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv









def dfit(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

    def dchange(fullarray,*dparms):

        '''
        Only d0s changing
        '''
        #d0s = list(d0s)
        twothet,azi,lenbg,tmpbackg = fullarray
        dparms=np.array(dparms)
        #print params.shape
        if (len(dparms.shape) > 1):
            if np.any(np.array(dparms.shape) > 1) :
                dparms = np.squeeze(dparms)
            elif np.all(np.array(dparms.shape) == 1):
                dparms = np.squeeze(dparms)
                dparms = np.array([dparms],float)   
        dparms=tuple(dparms)
        d_0 = Fourier_expand(azi, dparms)
        # print d_0
        H_all=heights
        W_all=widths
        #print wavelength
        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

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
        # print stop
        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg_all

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
    popt,pcurv = curve_fit(dchange,(twotheta,azimu,lenbg,flatbg),intens,p0=d0s)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv





def hfit(intens,twotheta,azimu,newd0s,widths,hparms,wavelength,bg=None):


    def hchange(fullarray,*hparms):

        '''
        Only h changing
        '''
        #d0s = list(d0s)
        twothet,azi,lenbg,tmpbackg = fullarray

        hparms=np.array(hparms)
        if (len(hparms.shape) > 1):
            if np.any(np.array(hparms.shape) > 1) :
                hparms = np.squeeze(hparms)
            elif np.all(np.array(hparms.shape) == 1):
                hparms = np.squeeze(hparms)
                hparms = np.array([hparms],float)   
        hparms=tuple(hparms)
        d_0 = newd0s
        H_all = Fourier_expand(azi, hparms)
        W_all = widths
        #print wavelength
        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

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
        # print stop
        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg_all

        return Int
    '''
    Only h fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
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
    popt,pcurv = curve_fit(hchange,(twotheta,azimu,lenbg,flatbg),intens,p0=hparms)

    return popt,pcurv






def wfit(intens,twotheta,azimu,newd0s,newheights,wparms,wavelength,bg=None):


    def wchange(fullarray,*wparms):

        '''
        Only w changing
        '''
        #d0s = list(d0s)
        twothet,azi,lenbg,tmpbackg = fullarray
        wparms=np.array(wparms)
        if (len(wparms.shape) > 1):
            if np.any(np.array(wparms.shape) > 1) :
                wparms = np.squeeze(wparms)
            elif np.all(np.array(wparms.shape) == 1):
                wparms = np.squeeze(wparms)
                wparms = np.array([wparms],float)   
        wparms=tuple(wparms)

        d_0 = newd0s
        H_all = newheights
        W_all = Fourier_expand(azi, wparms)
        #print wavelength
        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

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
        # print stop
        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg_all

        return Int
    '''
    Only w fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
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
    popt,pcurv = curve_fit(wchange,(twotheta,azimu,lenbg,flatbg),intens,p0=wparms)

    return popt,pcurv



def Allchange(fullarray,*allparms):

    ##fit d,w,h as single number with no Fourier terms

    twothet,azi,lenbg,parmnums,wavelength = fullarray
    # print twothet.shape,azi.shape,parmnums.shape
    allparms=np.array(allparms)
    if (len(allparms.shape) > 1):
        if np.any(np.array(allparms.shape) > 1) :
            allparms = np.squeeze(allparms)
        elif np.all(np.array(allparms.shape) == 1):
            allparms = np.squeeze(allparms)
            allparms = np.array([allparms],float)  
    start = 0
    starth = parmnums[0]
    startw = parmnums[0:2].sum()
    end = parmnums.sum()
    d0s = allparms[0:starth]
    heights = allparms[starth:startw]
    widths = allparms[startw:end]
    tmpbackg = allparms[end:]
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
    H_all = Fourier_expand(azi, heights)
    W_all = Fourier_expand(azi, widths)
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
    bg_all = Fourier_backgrnd((azi,twothet),backg)
    # print d_0,H_all,W_all,tth0,bg

    newall =  (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg_all

    #print newall.shape

    return newall


def Allfit(intens,twotheta,azimu,d0s,heights,widths,wavelength,bg=None):


    '''
    All expansion parameters fitted
    '''

    parmnums = np.array([len(d0s),len(heights),len(widths)])
    # print 'parmnums', parmnums

    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    allparms = np.concatenate((d0s,heights,widths,[item for sublist in bg for item in sublist]), axis=None)
    lenbg=[]
    for val in bg:
        lenbg.append(len(val))
    lenbg=np.array(lenbg)
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(Allchange,(twotheta,azimu,lenbg,parmnums,wavelength),intens,p0=allparms)
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
