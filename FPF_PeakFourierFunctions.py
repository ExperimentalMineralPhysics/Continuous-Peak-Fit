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

    # print type(param)
    # print param
    param=np.array(param)
    # print param.shape
    # print np.any(np.array(param.shape)>1)
    if (len(param.shape) > 1):
        if np.any(np.array(param.shape) > 1) :
            param = np.squeeze(param)
            # print param,'1'
        elif np.all(np.array(param.shape) == 1):
            param = np.squeeze(param)
            param = np.array([param],float)
            # print param, '2'
    # print len(param.shape)

    param=tuple(param)
    # print param, 'param'
    out = param[0] #make to size of azimu
    # essentially d_0, h_0 or w_0

    if len(param)>1:
        for i in xrange(1, ((len(param)-1)/2)+1): # len(param)-1 should never be odd because of initial a_0 parameter
            out = out + param[(2*i)-1]*np.sin(np.deg2rad(azimu))**i + param[2*i]*np.cos(np.deg2rad(azimu))**i #single col array
    # print out, 'out'
    else:
        fout = np.ones(azimu.shape)
        #fout[:] = out
        out = fout*param

    #print stop
    return out

    # if len(param)>1:
    #     for i in xrange(1, ((len(param)-1)/2)+1): # len(param)-1 should never be odd because of initial a_0 parameter
    #         out = out + param[(2*i)-1]*np.sin(azimu)**i + param[2*i]*np.cos(azimu)**i #single col array
    # # print out, 'out'
    # #print stop

    # #print out.shape
    # return out

def Fourier_fit(azimu,ydata,terms,param=None,errs=1):


    #print param
    if param:
        param=param
    else:
        param = [0 for i in range((2*terms+1))]
    param=tuple(param)
    # print 'fourier_fit', param
    #print errs, 'errs'

    #if errs:
    #    errs = errs
    # else:
    #    errs
    #print param.shape
    #print stop
    popt,pcurv = curve_fit(Fourier_expand,azimu,ydata,p0=param,sigma=errs)

    #print 'popt'
    #print popt
    #print 'pcurv'
    #print pcurv

    return popt,pcurv


##warning about number of free parms e.g. if n is num. chunks/2

'''
def singleInt_old(twotheta,d0s,heights,widths,backg):


    ##fit d,w,h as single number with no Fourier terms
    # print d0s
    d_0=d0s
    H_all=heights
    W_all=widths
    #print wavelength
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    # bg either single num. or size of azimu array
    #print tth0
    # print d_0,H_all,W_all,tth0,bg
    # print type(d_0),type(H_all),type(W_all),type(tth0),type(bg),type(twotheta),type(wavelength)
    #print d_0.dtype,H_all.dtype,W_all.dtype,tth0.dtype,twotheta.dtype
    Int = (H_all*np.exp((-(twotheta-tth0)**2)/(2*W_all**2))) + backg

    return Int



def singleInt_test(twotheta,*sparms):


    ##fit d,w,h as single number with no Fourier terms

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


    d_0=sparms[0]
    H_all=sparms[1]
    W_all=sparms[2]
    backg = sparms[3]
    #'
    # d_0=d0s
    # H_all=heights
    # W_all=widths
    #''
    #print wavelength
    #tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    # bg either single num. or size of azimu array
    #print tth0
    # print d_0,H_all,W_all,tth0,bg
    # print type(d_0),type(H_all),type(W_all),type(tth0),type(bg),type(twotheta),type(wavelength)
    #print d_0.dtype,H_all.dtype,W_all.dtype,tth0.dtype,twotheta.dtype
    #Int = (H_all*np.exp((-(twotheta-tth0)**2)/(2*W_all**2))) + backg

    #return Int



def singleInt_old_inV8(twotheta,*sparms):

    
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
    backg = np.array(sparms[3:])
    # backg = sparms[3]
    # print d0s
    # print d_0,H_all,W_all,backg, 'd_0,H_all,W_all,backg'
    wavelength = []
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
    #print twotheta.shape[0]
    #print type(backg)

    if isinstance(backg, (list, tuple, np.ndarray)):
        bg_all = backg[0]
        bg_out = np.ones(twotheta.shape)
        bg_out[:] = bg_all
        bg_all = bg_out
        if len(backg)>1:
            for i in xrange(1,len(backg)):
                bg_all = bg_all + (backg[i]*(twotheta**float(i)))
    else:
        bg_all = np.ones(twotheta.shape)
        bg_all[:] = backg

    # bg either single num. or size of azimu array
    # print bg_all, 'bg_all'
    Int = (H_all*np.exp((-(twotheta-tth0)**2)/(2*W_all**2))) + bg_all

    return Int
'''


def singleInt(twothetaWL,*sparms):

    #global wavelength

    #split two theta and wavelength (if present)
    if len(twothetaWL) == 2:
        twotheta = twothetaWL[0]
        wavelength = twothetaWL[1]
    else:
        twotheta = twothetaWL

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
    backg = np.array(sparms[3:])
    # backg = sparms[3]
    # print d0s
    # print d_0,H_all,W_all,backg, 'd_0,H_all,W_all,backg'
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
    #print twotheta.shape[0]
    #print type(backg)

    if isinstance(backg, (list, tuple, np.ndarray)):
        bg_all = backg[0]
        bg_out = np.ones(twotheta.shape)
        bg_out[:] = bg_all
        bg_all = bg_out
        if len(backg)>1:
            for i in xrange(1,len(backg)):
                bg_all = bg_all + (backg[i]*(twotheta**float(i)))
    else:
        bg_all = np.ones(twotheta.shape)
        bg_all[:] = backg

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
    allparms = np.concatenate((d0s,heights,widths,bg), axis=None)
    # print allparms
    popt,pcurv = curve_fit(singleInt, (twotheta, wavelength), intens, p0=allparms,maxfev=8000)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv


# def singleFit_test(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

#     '''
#     All expansion parameters fitted
#     '''
#     #check for bg here, won't work if try to pass to fit
#     if not bg:
#         bg = backg
#     # print type(d0s),d0s, intens.dtype
#     intens=np.array(intens,dtype='float64')
#     twotheta=np.array(twotheta,dtype='float64')
#     sparms = np.concatenate((d0s,heights,widths,bg), axis=None)
#     # print type(d0s),d0s, intens.dtype
#     popt,pcurv = curve_fit(singleInt, twotheta,intens,p0=sparms,maxfev=8000)
#     #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
#     return popt,pcurv



# def singleFit_old(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

#     '''
#     All expansion parameters fitted
#     '''
#     #check for bg here, won't work if try to pass to fit
#     if not bg:
#         bg = backg
#     # print type(d0s),d0s, intens.dtype
#     intens=np.array(intens,dtype='float64')
#     twotheta=np.array(twotheta,dtype='float64')
#     # print type(d0s),d0s, intens.dtype
#     popt,pcurv = curve_fit(singleInt, twotheta,intens,p0=[d0s,heights,widths, bg],maxfev=8000)
#     #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
#     return popt,pcurv





def dfit(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

    def dchange(fullarray,*dparms):

        '''
        Only d0s changing
        '''
        #d0s = list(d0s)
        twothet,azi = fullarray
        # print dparms, 'd0s 2'
        d_0 = Fourier_expand(azi, dparms)
        # print d_0
        H_all=heights
        W_all=widths
        #print wavelength
        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
        # print stop
        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg

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
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(dchange,(twotheta,azimu),intens,p0=d0s)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv





def hfit(intens,twotheta,azimu,newd0s,widths,hparms,wavelength,bg=None):


    def hchange(fullarray,*hparms):

        '''
        Only h changing
        '''
        #d0s = list(d0s)
        twothet,azi = fullarray
        # print dparms, 'd0s 2'
        d_0 = newd0s
        H_all = Fourier_expand(azi, hparms)
        W_all = widths
        #print wavelength
        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
        # print stop
        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg

        return Int
    '''
    Only h fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(hchange,(twotheta,azimu),intens,p0=hparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv






def wfit(intens,twotheta,azimu,newd0s,newheights,wparms,wavelength,bg=None):


    def wchange(fullarray,*wparms):

        '''
        Only w changing
        '''
        #d0s = list(d0s)
        twothet,azi = fullarray
        # print dparms, 'd0s 2'
        d_0 = newd0s
        H_all = newheights
        W_all = Fourier_expand(azi, wparms)
        #print wavelength
        tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))
        # print stop
        Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg

        return Int
    '''
    Only w fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(wchange,(twotheta,azimu),intens,p0=wparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv



def Allchange(fullarray,*allparms):

    ##fit d,w,h as single number with no Fourier terms

    twothet,azi,parmnums,wavelength,bg = fullarray
    # print twothet.shape,azi.shape,parmnums.shape
    allparms=np.array(allparms)
    start = 0
    starth = parmnums[0]
    startw = parmnums[0:2].sum()
    end = parmnums.sum()
    d0s = allparms[0:starth]
    heights = allparms[starth:startw]
    widths = allparms[startw:end]
    # print d0s,heights,widths
    d_0 = Fourier_expand(azi, d0s)
    H_all = Fourier_expand(azi, heights)
    W_all = Fourier_expand(azi, widths)
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    # print d_0,H_all,W_all,tth0,bg

    newall =  (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + bg

    #print newall.shape

    return newall


def Allfit(intens,twotheta,azimu,d0s,heights,widths,wavelength,bg=None):


    '''
    All expansion parameters fitted
    '''

    parmnums = np.array([len(d0s),len(heights),len(widths)])
    # print 'parmnums', parmnums
    #print d0s.shape,heights.shape,widths.shape
    allparms = np.concatenate((d0s,heights,widths), axis=None)
    #print parmnums,d0s,heights,widths
    # print allparms

    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(Allchange,(twotheta,azimu,parmnums,wavelength,bg),intens,p0=allparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv

