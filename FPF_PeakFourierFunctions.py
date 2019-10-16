# /usr/bin/python

import numpy as np
import numpy.ma as ma
import sys#, copy, os
#from PIL import Image
#import math
#import matplotlib.cm as cm
#import matplotlib.pyplot as plt
#import matplotlib.path as mlp
from scipy.optimize import curve_fit
from scipy.optimize import minimize
np.set_printoptions(threshold=sys.maxsize)


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
        Arr[0].append(p[:])
    for y in range(LenBG):
        Arr[1].extend([v][:])
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



def PeaksModel(twotheta, azi, Shapes, Conv=None, symm=None, PenCalc=None):
    
    # Full model of intensities at twotheat and azi given parameters.    
   
    I = Fourier_backgrnd((azi,twotheta),Shapes['background'])
    
    if not PenCalc==None:
        #Penatly Function for background
        POut = FitPenaltyFunction(azi, I, [0, np.inf])
    
    #loop over the number of peaks
    Ipeak = []
    

            
    for a in range(len(Shapes['peak'])):
        
        #symmety 
        if 'symmetry' in Shapes['peak'][a].keys():
            symm = Shapes['peak'][a]['symmetry']
        else:
            symm = 1
        
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
            POut = POut * dpen * hpen * wpen * ppen
    if PenCalc==None:
        POut = Ipeak
        
    return I,POut




def ChangeParams(constants, *fitparam):
    twotheta    = constants[0]
    azi         = constants[1]
    ChangeArray = constants[2]
    Shapes      = constants[3]
    Conv        = constants[4]
    symm        = constants[5]
    fixed       = constants[6]
    
    Shapes = GuessApply(ChangeArray, Shapes, fitparam)
    I,pen = PeaksModel(twotheta, azi, Shapes, Conv, symm, PenCalc=1)
    
    return I#*pen




def FitModel(Intfit, twotheta, azimu, ChangeArray, Shapes, Conv=None, symm=None, fixed=None, method=None, weights=None, bounds=None):
    
    p0array = GuessGet(ChangeArray, Shapes)
    
    if not (method=='minimize'):
        
        
        # FIX ME: replace with calls to LMFIT. https://lmfit.github.io/lmfit-py/model.html or equivalnet to get constrained fits.
        try:
            popt,pcurv = curve_fit(ChangeParams,(twotheta,azimu,ChangeArray,Shapes,Conv,symm,fixed),Intfit,p0=p0array, ftol=2e-12, maxfev = 12000, method='lm', sigma=weights )    
        except:
            
            popt  = [np.nan] * np.empty(len(p0array))
            pcurv = [np.nan] * np.empty((len(p0array),len(p0array),))
        
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

    Shapes      = GuessApply(ChangeArray, Shapes, popt.tolist(), np.sqrt(np.abs(np.diag(pcurv))).tolist())
    
    return Shapes, popt, pcurv






def MinimiseParams(p0, intens,twotheta,azimu,ChangeArray,d0s,Hs,Ws,Ps,bg,Conv,symm,fixed):
    
    d0s, Hs, Ws, Ps, bg = GuessApply(ChangeArray, d0s, Hs, Ws, Ps, bg, p0)
    
    I,pen = PeaksModel(twotheta, azimu, d0s, Hs, Ws, Ps, bg, Conv, symm)

    SSD = np.sum(((intens-I)**2).flatten())
    
    return SSD
        
    

def GuessGet(ChangeArray, Guesses):
    
    #ran = len(Guesses['peak'])
        
    Guess = []
    
    #print 'line 380', Guesses
    for a in range(len(Guesses['peak'])):
        if not ChangeArray[0][a][0] == 0:  #d0
            e = Guesses['peak'][a]['d-space']
            Guess.extend([e][:])
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
            Guess.extend([e][:])
#            e = Hs[a][0:ChangeArray[0][a][1]]
            #Guess.extend(e[:])
        if not ChangeArray[0][a][2] == 0:  #W
            e = Guesses['peak'][a]['width']
            Guess.extend([e][:])
#            e = Ws[a][0:ChangeArray[0][a][2]]
            #Guess.extend(e[:])
        if not ChangeArray[0][a][3] == 0:  #P
            e = Guesses['peak'][a]['profile']
            Guess.extend([e][:])
#            e = Ps[a][0:ChangeArray[0][a][3]]
            #Guess.extend(e[:])
    for a in range(len(Guesses['background'])):
        if not ChangeArray[1][a] == 0:  #background
            e = Guesses['background'][a]
            Guess.extend([e][:])
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



    
def GuessApply(ChangeArray, V, values, *args):
    
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




def CentroidConversion(Conv, args_in, azi):
    #Conv['DispersionType'] is the conversion type
    #Conv[...] are the values required for the conversion
    #FIX ME: these functions should be a sub function of the detector types. but need to work out how to do it.
    

    if Conv['DispersionType'] == None or Conv==0:
        args_out = args_in
        
    elif Conv['DispersionType'] == 'AngleDispersive':
        #args_in are d-spacing
        #args_outis two thetas
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

    #get NaN values.
    idx = np.isfinite(azimu) & np.isfinite(ydata) & (ydata > 0)
    if(type(terms)==list):
        terms = terms[0]
    if param:
        param=param
    else:
        param = [0 for i in range((2*terms+1))]
        
    if errs is None or errs.all == None:
        errs = np.ones(ydata.shape)
                    
    popt,pcurv = curve_fit(Fourier_expand,azimu[idx],ydata[idx],p0=param,sigma=errs[idx])


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
        Weight = 1000#0000
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

    Penalise = Penalise**Power*Weight #+ 1
    return Penalise



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

