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
from scipy.optimize import minimize
np.set_printoptions(threshold='nan')


# Fitting functions required for fitting peaks to the data. 
# Called by FPF_XRD_FitSubpattern.


# Known limitations/bugs
# FIX ME: need hard limits to gauss - lorentz peak profile ratios


def flatten(li):
    return sum(([x] if not isinstance(x, list) else flatten(x)
        for x in li), [])



def ParamFitArr(numP, LenBG, val=None, p=None, b=None, orda=None):
    #Makes array which represents all the fourier series in the fit.
    # Arr{0} is a list of 4 element lists - one 4 element list for each peak.
    # Arr[0] is the background representative.
    
    if not val==None:
        v = val
    else:
        v = 0
        
    if p==None:
        p = [v,v,v,v]
        
        
    Arr = [[],[]]
    for y in range(numP):
        Arr[0].append(p)
    for y in range(LenBG):
        Arr[1].extend([v])
    return Arr

def ParamFitArr2(orda, *args):
    #Makes array which represents all the fourier series in the fit.
    # Arr{0} is a list of 4 element lists - one 4 element list for each peak.
    # Arr[1] is the background representative.
    # Makes all values 0 apart from those listed.
    
    #args is a list of the parts to make non-zero and to change
    #possible options are:
    # 'all', 'orders' -- pastes orders from orda
    # 'all', 1 -- makes all values 1
    # 'background'
    # 'peak', n
    # 'd-space', n
    # 'height', n
    # 'width', n
    # 'profile', n
    
    
    if isinstance(orda,dict):
        numP = len(orda['peak'])
        LenBG = len(orda['background'])
    else:
        numP = orda[0]
        LenBG = orda[1]
    
    if not val==None:
        v = val
    else:
        v = 0
        
    if p==None:
        p = [v,v,v,v]
        
        
    Arr = [[],[]]
    for y in range(numP):
        Arr[0].append(p)
    for y in range(LenBG):
        Arr[1].extend([v])
    return Arr


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
    #if np.any(LGratio<0) or np.any(LGratio>1):
        #print 'The ratio is out of range'
        #stop

    PVpeak = LGratio*GaussianPeak(twotheta, tth0, W_all, H_all) + (1-LGratio)*LorentzianPeak(twotheta, tth0, W_all, H_all)
    return PVpeak



def PeaksModel2(twotheta, azi, Shapes, Conv=None, symm=None, PenCalc=None):
    
    # Full model of intensities at twotheat and azi given parameters.    
    
    
    #make backgorund fourier from array (either original or fitted array)
    #FIX ME: why is this not done outside this fitting routine? 
#    backg = []
#    for b in xrange(len(lenbg)):
#        if b == 0:
#            start = 0
#            end = lenbg[b]
#        else: 
#            start = lenbg[b-1]
#            end = lenbg[b]+start 
#        backg.append(tmpbackg[start:end])
    
    I = Fourier_backgrnd((azi,twotheta),Shapes['background'])
    
    if not PenCalc==None:
        #Penatly Function for background
        POut = FitPenaltyFunction(azi, I, [0, np.inf])
    
    #loop over the number of peaks
    Ipeak = []
    
    #symmety 
    if symm==None:
        symm = 1
        
    for a in range(len(Shapes['peak'])):
        #print 'HERE BE DRAGONS!', Shapes,
        Dall = Fourier_expand(azi, Shapes['peak'][a]['d-space'])
        Hall = Fourier_expand(azi*symm, Shapes['peak'][a]['height'])
        Wall = Fourier_expand(azi*symm, Shapes['peak'][a]['width'])
        Pall = Fourier_expand(azi*symm, Shapes['peak'][a]['profile'])
  
        #conversion
        TTHall = CentroidConversion(Conv, Dall, azi)
        Ipeak.append(PseudoVoigtPeak(twotheta, TTHall, Wall, Hall, Pall))
        I = I + Ipeak[a]       
       
        if not PenCalc==None:
            #Penatly Function
            #apply penatly function to constrain the Fourier series
            # FIX ME: this should probably be replaced with a np.minimize function rather than the curve_fit call.
            dpen  = FitPenaltyFunction(azi, TTHall,    [np.min(twotheta), np.max(twotheta)])
            hpen  = FitPenaltyFunction(azi*symm, Hall,   [0, np.max(Hall)*1.1])
            wpen  = FitPenaltyFunction(azi*symm, Wall,   [0, (np.max(twotheta)-np.min(twotheta))])
            ppen  = FitPenaltyFunction(azi*symm, Pall, [0, 1], Weight=5)
            #print 'bgpen', POut
#            print 'dpen', dpen
#            print 'hpen', hpen
#            print 'wpen', wpen
#            print 'ppen', ppen
            POut = POut * dpen * hpen * wpen * ppen
            #print 'POut', POut
    if PenCalc==None:
        POut = Ipeak
        
    return I,POut

#def PeaksModel(twotheta, azi, d0s, Hs, Ws, Ps, bg, Conv=None, symm=None, PenCalc=None):
#    
#    # Full model of intensities at twotheat and azi given parameters.    
#    
#    
#    #make backgorund fourier from array (either original or fitted array)
#    #FIX ME: why is this not done outside this fitting routine? 
##    backg = []
##    for b in xrange(len(lenbg)):
##        if b == 0:
##            start = 0
##            end = lenbg[b]
##        else: 
##            start = lenbg[b-1]
##            end = lenbg[b]+start 
##        backg.append(tmpbackg[start:end])
#    
#    I = Fourier_backgrnd((azi,twotheta),bg)
#    
#    if not PenCalc==None:
#        #Penatly Function for background
#        POut = FitPenaltyFunction(azi, I, [0, np.inf])
#    
#    #loop over the number of peaks
#    Ipeak = []
#    
#    #symmety 
#    if not symm==None:
#        azi = azi*symm
#        
#    for a in range(len(d0s)):
#        
#        Dall = Fourier_expand(azi, d0s[a])
#        Hall = Fourier_expand(azi, Hs[a])
#        Wall = Fourier_expand(azi, Ws[a])
#        Pall = Fourier_expand(azi, Ps[a])
#  
#        #conversion
#        TTHall = CentroidConversion(Conv, Dall, azi)
#        Ipeak.append(PseudoVoigtPeak(twotheta, TTHall, Wall, Hall, Pall))
#        I = I + Ipeak[a]       
#       
#        if not PenCalc==None:
#            #Penatly Function
#            #apply penatly function to constrain the Fourier series
#            # FIX ME: this should probably be replaced with a np.minimize function rather than the curve_fit call.
#            dpen  = FitPenaltyFunction(azi, TTHall,    [np.min(twotheta), np.max(twotheta)])
#            hpen  = FitPenaltyFunction(azi, Hall,   [0, np.max(Hall)*1.1])
#            wpen  = FitPenaltyFunction(azi, Wall,   [0, (np.max(twotheta)-np.min(twotheta))])
#            ppen  = FitPenaltyFunction(azi, Pall, [0, 1], Weight=5)
#            #print 'bgpen', POut
##            print 'dpen', dpen
##            print 'hpen', hpen
##            print 'wpen', wpen
##            print 'ppen', ppen
#            POut = POut * dpen * hpen * wpen * ppen
#            #print 'POut', POut
#    if PenCalc==None:
#        POut = Ipeak
#        
#    return I,POut


def ChangeParams2(constants, *fitparam):
    
    
    twotheta    = constants[0]
    azi         = constants[1]
    ChangeArray = constants[2]
    Shapes      = constants[3]
    Conv        = constants[4]
    symm        = constants[5]
    fixed       = constants[6]
    
    Shapes = GuessApply2(ChangeArray, Shapes, fitparam)
    I,pen = PeaksModel2(twotheta, azi, Shapes, Conv, symm, PenCalc=1)
    
    return I#*pen

#def ChangeParams(constants, *fitparam):
#    
#    twotheta    = constants[0]
#    azi         = constants[1]
#    ChangeArray = constants[2]
#    dspc        = constants[3]
#    Hits        = constants[4]
#    Wets        = constants[5]
#    Profs       = constants[6]
#    bg          = constants[7]
#    Conv        = constants[8]
#    symm        = constants[9]
#    fixed       = constants[10]
#    
#    dspc, Hits, Wets, Profs, bg = GuessApply(ChangeArray, dspc, Hits, Wets, Profs, bg, fitparam)
#        
#    
#    I,pen = PeaksModel(twotheta, azi, dspc, Hits, Wets, Profs, bg, Conv, symm, PenCalc=1)
#    #print pen
#    return I#*pen


def FitModel2(Intfit, twotheta, azimu, ChangeArray, Shapes, Conv=None, symm=None, fixed=None, method=None):
    
    # FIX ME: Should the dictionary structure be passed into the curve_fit function or will be be quicker as a set of arrays?
    #print 'line 255 ChageArray', ChangeArray
    #print 'line 256, Shapes', Shapes
    p0array = GuessGet2(ChangeArray, Shapes)
    
    #print 'line 258, p0array', p0array
    #flatten bg
#    if not bgF:
#        bgF = Shapes['background']
#    flatbg = np.array(np.concatenate(([item for sublist in bgF for item in sublist]), axis=None))
#    lenbg=[]
#    for val in bgF:
#        lenbg.append(len(val))
#    lenbg=np.array(lenbg)
    
    #Intfit=np.array(Intfit,dtype='float64')
    #twotheta=np.array(twotheta,dtype='float64')
    #azimu = np.array(azimu,dtype='float64')
    
    #d0sF,HsF,WsF,PsF
    
    if not (method=='minimize'):
        
#        import pdb; pdb.set_trace()
        popt,pcurv = curve_fit(ChangeParams2,(twotheta,azimu,ChangeArray,Shapes,Conv,symm,fixed),Intfit,p0=p0array, ftol=2e-12, maxfev = 30000, method='lm')    
        
    else:
        #FIX ME: should curve fit be replaces with minimise? -- so that we can force the constraints.

        #Compose constraints array.
        cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
                {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
                {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})

        prm = minimize(MinimiseParams, p0array, args=(Intfit,twotheta,azimu,ChangeArray,d0sF[:],HsF[:],WsF[:],PsF[:],bgF[:],Conv,symm,fixed), tol=1E-10 )
        popt = prm.x
        print 'prm.x', prm.x
        pcurv=0

    #print d0sF, 'Before apply'
#    print 'Shapes \n', Shapes
    Shapes      = GuessApply2(ChangeArray, Shapes, popt.tolist(), np.sqrt(np.abs(np.diag(pcurv))).tolist())
#    Shapes_errs = Shapes.copy()
##    Shapes_errs = GuessApply2(ChangeArray, Shapes_errs, np.sqrt(np.abs(np.diag(pcurv))))
#    print 'Shapes \n', Shapes
#    
##    print 'Shapes errors\n', Shapes_errs
#    stop
#    d0sFerr, HsFerr, WsFerr, PsFerr, bgFerr = GuessApply(ChangeArray, np.zeros(np.shape(d0sF)).tolist(), np.zeros(np.shape(HsF)).tolist(), np.zeros(np.shape(WsF)).tolist(), np.zeros(np.shape(PsF)).tolist(), np.zeros(np.shape(bgF)).tolist(), np.sqrt(np.abs(np.diag(pcurv))))
#
##    print d0sF, 'After application', type(d0sF)
##    print d0sFerr, 'After application', type(d0sFerr)
#
    #print 'ChangeArr', ChangeArray
    #print 'potp', popt
    #print 'Shapes', Shapes
#    stop

    return Shapes, popt, pcurv



#def FitModel(Intfit, twotheta, azimu, ChangeArray, d0sF, HsF, WsF, PsF, bgF, Conv=None, symm=None, fixed=None, method=None):
#    
#    p0array = GuessGet(ChangeArray, d0sF, HsF, WsF, PsF, bgF)
#
#    #flatten bg
#    if not bgF:
#        bgF = backg
#    flatbg = np.array(np.concatenate(([item for sublist in bgF for item in sublist]), axis=None))
#    lenbg=[]
#    for val in bgF:
#        lenbg.append(len(val))
#    lenbg=np.array(lenbg)
#    
#    Intfit=np.array(Intfit,dtype='float64')
#    twotheta=np.array(twotheta,dtype='float64')
#    azimu = np.array(azimu,dtype='float64')
#    
#    if not (method=='minimize'):
#        
##        import pdb; pdb.set_trace()
#
#        popt,pcurv = curve_fit(ChangeParams,(twotheta,azimu,ChangeArray,d0sF,HsF,WsF,PsF,bgF,Conv,symm,fixed),Intfit,p0=p0array, ftol=2e-12, maxfev = 30000, method='lm')    
#        
#    else:
#        #FIX ME: should curve fit be replaces with minimise? -- so that we can force the constraints.
#
#        #Compose constraints array.
#        cons = ({'type': 'ineq', 'fun': lambda x:  x[0] - 2 * x[1] + 2},
#                {'type': 'ineq', 'fun': lambda x: -x[0] - 2 * x[1] + 6},
#                {'type': 'ineq', 'fun': lambda x: -x[0] + 2 * x[1] + 2})
#
#        prm = minimize(MinimiseParams, p0array, args=(Intfit,twotheta,azimu,ChangeArray,d0sF[:],HsF[:],WsF[:],PsF[:],bgF[:],Conv,symm,fixed), tol=1E-10 )
#        popt = prm.x
#        print 'prm.x', prm.x
#        pcurv=0
#
#    #print d0sF, 'Before apply'
#    
#    d0sF, HsF, WsF, PsF, bgF = GuessApply(ChangeArray, d0sF, HsF, WsF, PsF, bgF, popt)
#    d0sFerr, HsFerr, WsFerr, PsFerr, bgFerr = GuessApply(ChangeArray, np.zeros(np.shape(d0sF)).tolist(), np.zeros(np.shape(HsF)).tolist(), np.zeros(np.shape(WsF)).tolist(), np.zeros(np.shape(PsF)).tolist(), np.zeros(np.shape(bgF)).tolist(), np.sqrt(np.abs(np.diag(pcurv))))
#
##    print d0sF, 'After application', type(d0sF)
##    print d0sFerr, 'After application', type(d0sFerr)
#
#    return d0sF, HsF, WsF, PsF, bgF, d0sFerr, HsFerr, WsFerr, PsFerr, bgFerr, popt, pcurv




def MinimiseParams(p0, intens,twotheta,azimu,ChangeArray,d0s,Hs,Ws,Ps,bg,Conv,symm,fixed):
    
    d0s, Hs, Ws, Ps, bg = GuessApply(ChangeArray, d0s, Hs, Ws, Ps, bg, p0)
    
    I,pen = PeaksModel(twotheta, azimu, d0s, Hs, Ws, Ps, bg, Conv, symm)

    SSD = np.sum(((intens-I)**2).flatten())
    
    return SSD
        
    

def GuessGet2(ChangeArray, Guesses):
    
    #ran = len(Guesses['peak'])
        
    Guess = []
    
    #print 'line 380', Guesses
    for a in range(len(Guesses['peak'])):
        if not ChangeArray[0][a][0] == 0:  #d0
            e = Guesses['peak'][a]['d-space']
            Guess.extend([e])
#            if isinstance(d0s, list):
#                e = d0s[a][0:ChangeArray[0][a][0]]
#            elif ran>1:
#                e = d0s[a]#[0:ChangeArray[0][a][0]]
#                stop
#            else:
#                e = d0s
#                stop
            #Guess.extend(e[:])
        if not ChangeArray[0][a][1] == 0:  #H
            e = Guesses['peak'][a]['height']
            Guess.extend([e])
#            e = Hs[a][0:ChangeArray[0][a][1]]
            #Guess.extend(e[:])
        if not ChangeArray[0][a][2] == 0:  #W
            e = Guesses['peak'][a]['width']
            Guess.extend([e])
#            e = Ws[a][0:ChangeArray[0][a][2]]
            #Guess.extend(e[:])
        if not ChangeArray[0][a][3] == 0:  #P
            e = Guesses['peak'][a]['profile']
            Guess.extend([e])
#            e = Ps[a][0:ChangeArray[0][a][3]]
            #Guess.extend(e[:])
    for a in range(len(Guesses['background'])):
        if not ChangeArray[1][a] == 0:  #background
            e = Guesses['background'][a][0]
            Guess.extend([e])
#            e = bg[a][0:ChangeArray[1][a]]
            #Guess.extend(e[:])

    out=[]
    for i in range(len(Guess)):
        if 'list' in str(type(Guess[i])):
            for j in range(len(Guess[i])):
                out.append(Guess[i][j])
        else :
            out.append(Guess[i])
            
    return out


#def GuessGet(ChangeArray, d0s, Hs, Ws, Ps, bg):
#    
#    if isinstance(d0s,list):
#        ran = len(d0s)
#    else:
#        ran = d0s.size
#        
#    Guess = []
#    for a in range(ran):
#        if not ChangeArray[0][a][0] == 0:  #d0
#            if isinstance(d0s, list):
#                e = d0s[a][0:ChangeArray[0][a][0]]
#            elif ran>1:
#                e = d0s[a]#[0:ChangeArray[0][a][0]]
#                stop
#            else:
#                e = d0s
#                stop
#            Guess.extend(e[:])
#        if not ChangeArray[0][a][1] == 0:  #H
#            e = Hs[a][0:ChangeArray[0][a][1]]
#            Guess.extend(e[:])
#        if not ChangeArray[0][a][2] == 0:  #W
#            e = Ws[a][0:ChangeArray[0][a][2]]
#            Guess.extend(e[:])
#        if not ChangeArray[0][a][3] == 0:  #P
#            e = Ps[a][0:ChangeArray[0][a][3]]
#            Guess.extend(e[:])
#    for a in range(len(bg)):
#        if not ChangeArray[1][a] == 0:  #background
#            e = bg[a][0:ChangeArray[1][a]]
#            Guess.extend(e[:])
#                    
#    return Guess
    
def GuessApply2(ChangeArray, V, values, *args):
    
    #print type(values)
#    if isinstance(values, (tuple,)) or isinstance(values, (np.ndarray, np.generic) ):
#        values = list(values)
        #print 'changed values to tuple'
    
    UsedP = 0
    for a in range(len(V['peak'])):
        if not ChangeArray[0][a][0] == 0:  #d0
            V['peak'][a]['d-space'] = list(values[UsedP:UsedP+ChangeArray[0][a][0]])
            if args:
                V['peak'][a]['d-space_err'] = list(args[0][UsedP:UsedP+ChangeArray[0][a][0]])
            #d0s[a] = values[UsedP:UsedP+ChangeArray[0][a][0]]
            UsedP = UsedP + ChangeArray[0][a][0]
        if not ChangeArray[0][a][1] == 0:  #H
            V['peak'][a]['height'] = list(values[UsedP:UsedP+ChangeArray[0][a][1]])
            if args:
                V['peak'][a]['height_err'] = list(args[0][UsedP:UsedP+ChangeArray[0][a][1]])
            #Hs[a] = values[UsedP:UsedP+ChangeArray[0][a][1]]
            UsedP = UsedP + ChangeArray[0][a][1]
        if not ChangeArray[0][a][2] == 0:  #W
            V['peak'][a]['width'] = list(values[UsedP:UsedP+ChangeArray[0][a][2]])
            if args:
                V['peak'][a]['width_err'] = list(args[0][UsedP:UsedP+ChangeArray[0][a][2]])
            #Ws[a] = values[UsedP:UsedP+ChangeArray[0][a][2]]
            UsedP = UsedP + ChangeArray[0][a][2]
        if not ChangeArray[0][a][3] == 0:  #P
            V['peak'][a]['profile'] = list(values[UsedP:UsedP+ChangeArray[0][a][3]])
            if args:
                V['peak'][a]['profile_err'] = list(args[0][UsedP:UsedP+ChangeArray[0][a][3]])
            #Ps[a] = values[UsedP:UsedP+ChangeArray[0][a][3]]
            UsedP = UsedP + ChangeArray[0][a][3]
    for a in range(len(V['background'])):
        if not ChangeArray[1][a] == 0: #background
            V['background'][a] = list(values[UsedP:UsedP+ChangeArray[1][a]])
            if args:
                if not 'background_err' in V:
                    V['background_err'] = [[]]
                    for b in range(len(V['background'])-1):
                        V['background_err'].append([])
                V['background_err'][a] = list(args[0][UsedP:UsedP+ChangeArray[1][a]])
            #bg_[a] = values[UsedP:UsedP+ChangeArray[1][a]]
            UsedP = UsedP + ChangeArray[1][a]
    #print d0s, type(d0s), 'Guess Get'
    #return d0s, Hs, Ws, Ps, bg_
    return V


#def GuessApply(ChangeArray, d0s, Hs, Ws, Ps, bg_, values):
#    
#    
#    #print type(values)
#    if isinstance(values, (tuple,)) or isinstance(values, (np.ndarray, np.generic) ):
#        values = list(values)
#        #print 'changed values to tuple'
#        
#    #print values, type(values)
#    #print d0s, type(d0s)
#
#    UsedP = 0
#    for a in range(len(d0s)):
#        if not ChangeArray[0][a][0] == 0:  #d0
#            d0s[a] = values[UsedP:UsedP+ChangeArray[0][a][0]]
#            UsedP = UsedP + ChangeArray[0][a][0]
#        if not ChangeArray[0][a][1] == 0:  #H
#            Hs[a] = values[UsedP:UsedP+ChangeArray[0][a][1]]
#            UsedP = UsedP + ChangeArray[0][a][1]
#        if not ChangeArray[0][a][2] == 0:  #W
#            Ws[a] = values[UsedP:UsedP+ChangeArray[0][a][2]]
#            UsedP = UsedP + ChangeArray[0][a][2]
#        if not ChangeArray[0][a][3] == 0:  #P
#            Ps[a] = values[UsedP:UsedP+ChangeArray[0][a][3]]
#            UsedP = UsedP + ChangeArray[0][a][3]
#    for a in range(len(bg_)):
#        if not ChangeArray[1][a] == 0: #background
#            bg_[a] = values[UsedP:UsedP+ChangeArray[1][a]]
#            UsedP = UsedP + ChangeArray[1][a]
#    #print d0s, type(d0s), 'Guess Get'
#    return d0s, Hs, Ws, Ps, bg_


def CentroidConversion(Conv, args_in, azi):
    #Conv['DispersionType'] is the conversion type
    #Conv[...] are the vlaues required for the conversion
    #FIX ME: these functions should be a sub function of the detector types. but need to work out how to do it.
    
    # Commented out until add Energy dispersive diffraction.
    # When do this need to change wavelength parameter to 'Conversion' array with type and constant(s) in it.
    #print 'Conv', type(Conv), Conv
    if Conv['DispersionType'] == None or Conv==0:
        args_out = args_in
        
    elif Conv['DispersionType'] == 'AngleDispersive':
        #args_in are d-spacing
        #args_outis two thetas
        #wavelength = Conv[1]
        #wavelength = Conv['conversion_constant']
        args_out = 2*np.degrees(np.arcsin(Conv['conversion_constant']/2/args_in)) #FIX ME: check this is correct!!! 
        
    elif Conv['DispersionType'] == 'EnergyDispersive':
        
        args_out = []
        #for x in range(len(azi)):
        for x in range(azi.size):
            #determine which detector is being used
            #a=np.asscalar(np.where(Conv['azimuths'] == azi[x])[0])
            
            if azi.size == 1:
                a = np.array(np.where(Conv['azimuths'] == azi)[0][0])
                #a= np.asscalar(np.where(Conv['azimuths'] == azi)[0][0])
                args_out.append(12.398/(2*args_in*np.sin(np.radians(Conv['calibs'].mcas[a].calibration.two_theta/2))))
            else:
#                print type(azi)
#                print azi.mask[x]
                if not azi.mask[x]:
                    a=(np.where(Conv['azimuths'] == azi[x])[0][0])
                    args_out.append(12.398/(2*args_in[x]*np.sin(np.radians(Conv['calibs'].mcas[a].calibration.two_theta/2))))
                else:
                    args_out.append(0)
        if isinstance(azi,np.ma.MaskedArray):
            args_out = ma.array(args_out, mask=azi.mask)
    else:
        sys.exit('Unrecognised conversion type')
        
    return args_out


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

    #FIX ME: Need to check the fourier is a licit length
    param=np.array(param)
    if (len(param.shape) > 1):
        if np.any(np.array(param.shape) > 1) :
            param = np.squeeze(param)
            #print param,'1'
        elif np.all(np.array(param.shape) == 1):
            param = np.squeeze(param)
            param = np.array([param],float)
            #print param, '2'

    param=tuple(param)
    out = np.ones(azimu.shape)
    if azimu.size == 1: #this line of required to catch error when out is single number.
        out = param[0]
    else:
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



def Fourier_fit(azimu,ydata,terms,param=None,errs=None):

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
#    print ydata.shape
#    print errs.shape
#    print azimu.shape
#    print azimu
#    print errs
    popt,pcurv = curve_fit(Fourier_expand,azimu,ydata,p0=param,sigma=errs)


    return popt,pcurv



# fourier expansion function
def Fourier_backgrnd(azimutheta, param):

    azimu,twotheta = azimutheta
    twothetaprime = twotheta-twotheta.min()
    backg=param
    bg_all = np.zeros(twotheta.shape)
    
    for i in xrange(len(backg)):
        try:
            nterms = len(backg[i])
        except:
            nterms = backg[i].size
        out = backg[i][0]
        for j in xrange(1, ((nterms-1)/2)+1): 
            out = out + backg[i][(2*j)-1]*np.sin(np.deg2rad(azimu)*j) + backg[i][2*j]*np.cos(np.deg2rad(azimu)*j)
            # print j, 'j'
        bg_all = bg_all + (out*(twothetaprime**float(i)))

    return bg_all


def FitPenaltyFunction(Azimuths, FourierValues, valid_range, Weight=None, FV=None, Power=None):

    if Weight==None:
        Weight = 10#0000
    if Power==None:
        Power = 1
    
    if not np.array(Azimuths).shape==np.array(FourierValues).shape or not FV==None:
        Vals = Fourier_expand(Azimuths, FourierValues)
    else:
        Vals = FourierValues

    Penalise = 0
    if np.min(Vals) < valid_range[0]:
        Penalise = valid_range[0] - np.min(Vals)

    if np.max(Vals) > valid_range[1]:
        Penalise = np.max(Vals) - valid_range[0]

    Penalise = Penalise**Power*Weight + 1
    return Penalise







def singleInt(twothetaWL,*sparms):

    #global wavelength
    
    if len(sparms)==1:
        print 'oh dear!!!!'
        print sparms[0], type(sparms[0])
        sparms = sparms[0]
        
    #split two theta and wavelength (if present)
    if len(twothetaWL) == 5:
        twotheta,azimu,lenbg,wavelength,profile = twothetaWL
        
        Num_Peaks = np.int((len(sparms)-lenbg)/3)
        #profile = np.repeat(np.array(profile), Num_Peaks)
        
    elif len(twothetaWL) == 4:
        twotheta,azimu,lenbg,wavelength = twothetaWL
        
        Num_Peaks = np.int((len(sparms)-lenbg)/4)
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
    d_0  =np.array(sparms[0:(Num_Peaks)])
    H_all=np.array(sparms[(Num_Peaks)*1:(Num_Peaks)*2])
    W_all=np.array(sparms[(Num_Peaks)*2:(Num_Peaks)*3])
    if len(sparms) == lenbg+4*Num_Peaks:
        profile = np.array(sparms[lenbg[0]+3*Num_Peaks:])
    
    
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
    
    tth0 = CentroidConversion(wavelength, d_0)
    #tth1 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    bg_all = Fourier_backgrnd((azimu,twotheta),backg)
    #print bg_all
    #print backg
    # bg either single num. or size of azimu array
    # print bg_all, 'bg_all'
    
    #Int = PseudoVoigtPeak(twotheta, tth0, W_all, H_all, profile) + bg_all
    Int = bg_all
    for b in range(Num_Peaks):
        Int = Int + PseudoVoigtPeak(twotheta, tth0[b], W_all[b], H_all[b], profile[b])

    dpen  = FitPenaltyFunction(azimu, tth0,    [np.min(twotheta), np.max(twotheta)])
    hpen  = FitPenaltyFunction(azimu, H_all,   [0, np.max(H_all)*1.1])
    wpen  = FitPenaltyFunction(azimu, W_all,   [0, (np.max(twotheta)-np.min(twotheta))])
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
    
    if fixed == 0:  #profile is not fixed.
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

    else: # profile shape is fixed.
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
        print backg, lenbg
        bg_all = Fourier_backgrnd((azi,twothet),backg)

        print wavelength
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







    
    

