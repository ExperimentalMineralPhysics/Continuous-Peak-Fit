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

calib_file = 'CeO2_Pil207_E30_2Nov2016_001.imctrl'
pix = 172 #microns ##SAH: edit 172 microns not 17.2

opt = 2

if opt==1:
    diff_file = 'CeO2_Pil207_E30_2Nov2016_001.tif'
    mask_file = 'Calibration.immask'
    tthRange = [7.6,7.8]

    total = 25

    backg = 10.

    h_order = 0  ##  no formal limit
    w_order = 5  ## no formal limit
    d0_order = 0 ##  probably limited to 2 -- generalisation might work better for radio data.

elif opt==2:
    diff_file = 'BCC1_2GPa_10s_001_00001.tif'
    mask_file = 'Diffraction.immask'
    tthRange = [16.6,17.1]

    total = 72
    not_use = 5
    backg = [4.] #have to make sure this number is the right size for the bg order.

    h_order  = 12  ##  no formal limit
    w_order  = 0   ## no formal limit
    d0_order = 2 ##  probably limited to 2 -- generalisation might work better for radio data.
    bg_order = 0 ## no limit.

# GUI needed to set inputs? ###
# SAH: A GUI would be good for the initial inputs -- most geoscientist types find it is easier to select areas with a mouse than the keyboard.
# SAH: The GUI should save an output file that is then used as the input file for batch processing lots of files.


## Functions ##

def CreateGSASIIMask(MSKfile, ImInts, ImageSize, ImTTH, ImAzi, Imy,Imx):

    #GSAS-II lisence problems later.
    #As it stands we may have to get people to copy a GSAS-II file from their repository so that it does not
    #have to be distributed by us.
    #GSAS-II lisences give them the right to incorporate wholesale this code into GSAS-II

    msks = LoadMask(MSKfile)

    #make mask array as arrayon ones and zeors for now.
    #Can be transformed into an array-mask later.

    #msk = np.ones((ImageSize[1], ImageSize[2]))
    #ImMsk = np.ones((ImageSize[1], ImageSize[2]))#, dtype=np.int32)
    ImMsk = ImInts

    #Thresholds
    # copied from pyGSAS/GSASIIimage.py Fill2ThetaAzimuthMap
    # units of mask are intensity
    IntLims = msks['Thresholds'][1]
    #print IntLims[0], IntLims[1]
    ImMsk = ma.masked_outside(ImInts,int(IntLims[0]),IntLims[1])
    #ImMsk = ma.masked_outside(ImInts.flatten(),int(IntLims[0]),IntLims[1])


    #Rings
    # copied from pyGSAS/GSASIIimage.py Fill2ThetaAzimuthMap
    # units of mask are two theta and two theta (both in degrees)
    RngLims = (msks['Rings'])
    for twoth,thickness in RngLims:
        #print max(0.01,twoth-thickness/2.), twoth+thickness/2.
        #print type(ma.masked_inside(ImTTH.flatten(),max(0.01,twoth-thickness/2.),twoth+thickness/2.))
        #print type((ImMsk))
        ImMsk.mask = ma.mask_or(ma.getmask(ImMsk),ma.getmask(ma.masked_inside(ImTTH,max(0.01,twoth-thickness/2.),twoth+thickness/2.)))


    #Arcs
    # copied from pyGSAS/GSASIIimage.py Fill2ThetaAzimuthMap
    # units of mask are two theta and azimuth (both in degrees)
    ArcLims = (msks['Arcs'])
    for twoth,azim,thickness in ArcLims:
        tamt = ma.getmask(ma.masked_inside(ImTTH,max(0.01,twoth-thickness/2.),twoth+thickness/2.))
        tama = ma.getmask(ma.masked_inside(ImAzi,azim[0],azim[1]))
        ImMsk.mask = ma.mask_or(ma.getmask(ImMsk),tamt*tama)


    #Points/Spots
    # copied from pyGSAS/GSASIIimage.py Make2ThetaAzimuthMap
    # units of mask are position (x and y) on detector (in mm)
    spots = msks['Points']
    for spX,spY,spdiam in spots:
        tamp = ma.getmask(ma.masked_less((Imx-spX)**2+(Imy-spY)**2,(spdiam/2.)**2))
        ImMsk.mask = ma.mask_or(ma.getmask(ImMsk),tamp)


    #polygon
    # GSAS-II makes a polygon mask in pyGSAS/GSASIIimage.py MakeMaskMap
    # This code though calls a Fortran script.
    # This here is therefore an equivalent code block totally in python.
    # matplotlib.path is needed here.
    # Units of mask are position (x and y) on detector (in mm)
    PolyLms = msks['Polygons']
    points = np.vstack((Imx.flatten(),Imy.flatten())).T
    #FIX ME: The points array is flattened and then in a few lines time grid is reshaped. It is possible to do this without changing the shape of the arrays?
    for polygon in PolyLms:
        path = mlp.Path(polygon)
        grid = path.contains_points(points)
        grid = np.reshape(grid, (Imx.shape[0], Imx.shape[1]))
        ImMsk.mask = ma.mask_or(ma.getmask(ImMsk),grid)


    #frames
    # GSAS-II makes a frame mask in pyGSAS/GSASIIimage.py MakeMaskMap
    # This code though calls a Fortran script.
    # This here is therefore an equivalent code block totally in python.
    # units of mask are position (x and y) on detector (in mm)

    #FIX ME: I think that a frame excludes everything outside the point list, while polgon excludes everything inside the polgon.
    # The difference in the code for the two types is a -TRUE applied to the mask.
    # This should be checked though.

    FrmLms = msks['Frames']
    print FrmLms
    #for frame in FrmLms:
    #print frame
    if FrmLms:
        path = mlp.Path(FrmLms)
        grid = path.contains_points(points)-True
        print Imx.shape[0]
        print Imx.shape[1]
        print Imx.shape[0]*Imx.shape[1]
        print grid.shape
        grid = np.reshape(grid, (Imx.shape[0], Imx.shape[1]))
        ImMsk.mask = ma.mask_or(ma.getmask(ImMsk),grid)


    if 0:
        #This is left in here for debugging.
        fig = plt.figure()
        ax = fig.add_subplot(1,2,1)
        plt.subplot(121)
        plt.scatter(Imx, Imy, s=4, c=intens, edgecolors='none', cmap=plt.cm.jet)
        ax = fig.add_subplot(1,2,2)
        plt.subplot(122)
        plt.scatter(ma.array(Imx,mask=ImMsk.mask), ma.array(Imy,mask=ImMsk.mask), s=4, c=intens, edgecolors='none', cmap=plt.cm.jet)
        plt.colorbar()
        plt.show()

        plt.close()

    return ImMsk


def pointInPolygon(pXY,xy):
    # copied from pyGSAS/GSASIIimage.py
    'Needs a doc string'
    #pXY - assumed closed 1st & last points are duplicates
    Inside = False
    N = len(pXY)
    p1x,p1y = pXY[0]
    for i in range(N+1):
        p2x,p2y = pXY[i%N]
        if (max(p1y,p2y) >= xy[1] > min(p1y,p2y)) and (xy[0] <= max(p1x,p2x)):
            if p1y != p2y:
                xinters = (xy[1]-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
            if p1x == p2x or xy[0] <= xinters:
                Inside = not Inside
        p1x,p1y = p2x,p2y
    return Inside


def Fill2ThetaAzimuthMap(masks, azim, twoth, image):
    'Needs a doc string'
    Zlim = masks['Thresholds'][1]
    rings = masks['Rings']
    arcs = masks['Arcs']
    #TA = np.dstack((ma.getdata(TA[1]),ma.getdata(TA[0]),ma.getdata(TA[2])))    #azimuth, 2-theta, dist
    #tax,tay,tad = np.dsplit(TA,3)    #azimuth, 2-theta, dist**2/d0**2
    for tth,thick in rings:
        tam = ma.mask_or(tam.flatten(),ma.getmask(ma.masked_inside(tay.flatten(),max(0.01,tth-thick/2.),tth+thick/2.)))
    for tth,azm,thick in arcs:
        tamt = ma.getmask(ma.masked_inside(tay.flatten(),max(0.01,tth-thick/2.),tth+thick/2.))
        tama = ma.getmask(ma.masked_inside(tax.flatten(),azm[0],azm[1]))
        tam = ma.mask_or(tam.flatten(),tamt*tama)
    taz = ma.masked_outside(image.flatten(),int(Zlim[0]),Zlim[1])
    tabs = np.ones_like(taz)
    tam = ma.mask_or(tam.flatten(),ma.getmask(taz))
    tax = ma.compressed(ma.array(tax.flatten(),mask=tam))   #azimuth
    tay = ma.compressed(ma.array(tay.flatten(),mask=tam))   #2-theta
    taz = ma.compressed(ma.array(taz.flatten(),mask=tam))   #intensity
    #tad = ma.compressed(ma.array(tad.flatten(),mask=tam))   #dist**2/d0**2
    tabs = ma.compressed(ma.array(tabs.flatten(),mask=tam)) #ones - later used for absorption corr.
    return tax,tay,taz,tabs
    #return tax,tay,taz,tad,tabs

def LoadMask(mskfilename):
    #copied from GSASIIImgGUI.py OnLoadMask

    File = open(mskfilename,'r')
    save = {}
    #oldThreshold = data['Thresholds'][0]
    S = File.readline()
    while S:
        if S[0] == '#':
            S = File.readline()
            continue
        [key,val] = S[:-1].split(':')
        if key in ['Points','Rings','Arcs','Polygons','Frames','Thresholds']:
            #if key == 'Thresholds':
            #    S = File.readline()
            #    continue
            save[key] = eval(val)
            #if key == 'Thresholds':
            #    save[key][0] = oldThreshold
            #    save[key][1][1] = min(oldThreshold[1],save[key][1][1])
        S = File.readline()
    File.close()

    return save


def load_inputs_file(filenam):

    '''
    Parse inputs file, create type specific inputs.
    '''

    filelines = parms_file.readlines()
    parms_dict = {}

    for item in filelines:
        newparms = item.strip('\n').split(':',1)
        parm = newparms[1]

        value = None
        # print parm
        try:
            value = int(parm)
            parms_dict[str(newparms[0])] = value
        except ValueError:
            try:
                value = float(parm)
                parms_dict[newparms[0]] = value
            except ValueError:
                if parm.startswith('['):
                    listvals = parm.strip('[').strip(']').split(',')
                    newlist = []
                    for val in listvals:
                        newValue = None
                        try:
                            newValue = int(val)
                            newlist.append(newValue)
                        except ValueError:
                            try:
                                newValue = float(val)
                                newlist.append(newValue)
                            except ValueError:
                                newlist.append(val.replace("'","").replace(" ",""))
                    parms_dict[newparms[0]] = newlist
                elif parm.startswith('{'):
                    # print parm
                    listvals = parm.strip('{').strip('}').split(',')
                    newdict = {}
                    for keyval in listvals:
                        # print keyval
                        newkey = keyval.split(':')[0].replace("'","").replace(" ","")
                        val = keyval.split(':')[1]
                        newValue = None
                        try:
                            newValue = int(val)
                            newdict[str(newkey)] = newValue
                        except ValueError:
                            try:
                                newValue = float(val)
                                newdict[str(newkey)] = newValue
                            except ValueError:
                                newdict[str(newkey)] = val.replace("'","").replace(" ","")
                    parms_dict[newparms[0]] = newdict
                elif not parm:
                    parms_dict[newparms[0]] = ''

                else:
                    parms_dict[newparms[0]] = str(parm)


    return parms_dict




## Functions below taken from GSAS-II code see https://github.com/svaksha/pyGSAS/
## Toby, B. H., & Von Dreele, R. B. (2013). "GSAS-II: the genesis of a modern open-source
## all purpose crystallography software package". Journal of Applied Crystallography,
## 46(2), 544-549. ##

#trig functions
sind = lambda x: math.sin(x*math.pi/180.)
asind = lambda x: 180.*math.asin(x)/math.pi
tand = lambda x: math.tan(x*math.pi/180.)
atand = lambda x: 180.*math.atan(x)/math.pi
atan2d = lambda y,x: 180.*math.atan2(y,x)/math.pi
cosd = lambda x: math.cos(x*math.pi/180.)
acosd = lambda x: 180.*math.acos(x)/math.pi
rdsq2d = lambda x,p: round(1.0/math.sqrt(x),p)
#numpy trig functions
npsind = lambda x: np.sin(x*np.pi/180.)
npasind = lambda x: 180.*np.arcsin(x)/np.pi
npcosd = lambda x: np.cos(x*np.pi/180.)
npacosd = lambda x: 180.*np.arccos(x)/np.pi
nptand = lambda x: np.tan(x*np.pi/180.)
npatand = lambda x: 180.*np.arctan(x)/np.pi
npatan2d = lambda y,x: 180.*np.arctan2(y,x)/np.pi

def makeMat(Angle,Axis):
    '''Make rotation matrix from Angle and Axis
    :param float Angle: in degrees
    :param int Axis: 0 for rotation about x, 1 for about y, etc.
    '''
    cs = npcosd(Angle)
    ss = npsind(Angle)
    M = np.array(([1.,0.,0.],[0.,cs,-ss],[0.,ss,cs]),dtype=np.float32)
    return np.roll(np.roll(M,Axis,axis=0),Axis,axis=1)


def peneCorr(tth,dep,tilt=0.,azm=0.):
    'Needs a doc string'
    return dep*(1.-npcosd(tth))         #best one


def GetTthAzmDsp(x,y,data): #expensive
	## SAH: lifted from: https://subversion.xray.aps.anl.gov/pyGSAS/trunk/GSASIIimage.py
    '''
    '''
    wave = data['wavelength']
    cent = data['center']
    tilt = data['tilt']
    # print tilt, cosd(tilt)
    dist = data['distance']/cosd(tilt)
    x0 = data['distance']*tand(tilt)
    phi = data['rotation']
    dep = data['DetDepth']
    LRazim = data['LRazimuth']
    azmthoff = data['azmthOff']
    dx = np.array(x-cent[0],dtype=np.float32)
    dy = np.array(y-cent[1],dtype=np.float32)
    D = ((dx-x0)**2+dy**2+data['distance']**2)      #sample to pixel distance
    X = np.array(([dx,dy,np.zeros_like(dx)]),dtype=np.float32).T
    # print np.array(([dx,dy,np.zeros_like(dx)]),dtype=np.float32).shape
    X = np.dot(X,makeMat(phi,2))
    Z = np.dot(X,makeMat(tilt,0)).T[2]
    tth = npatand(np.sqrt(dx**2+dy**2-Z**2)/(dist-Z))
    dxy = peneCorr(tth,dep,tilt,npatan2d(dy,dx))
    DX = dist-Z+dxy
    DY = np.sqrt(dx**2+dy**2-Z**2)
    tth = npatan2d(DY,DX)
    dsp = wave/(2.*npsind(tth/2.))
    azm = (npatan2d(dy,dx)+azmthoff+720.)%360.
    G = D/data['distance']**2       #for geometric correction = 1/cos(2theta)^2 if tilt=0.
    return np.array([tth,azm,G,dsp])




#related calls#
def GetTth(x,y,data):
    'Give 2-theta value for detector x,y position; calibration info in data'
    return GetTthAzmDsp(x,y,data)[0]

def GetTthAzm(x,y,data):
    'Give 2-theta, azimuth values for detector x,y position; calibration info in data'
    return GetTthAzmDsp(x,y,data)[0:2]

def GetDsp(x,y,data):
    'Give d-spacing value for detector x,y position; calibration info in data'
    return GetTthAzmDsp(x,y,data)[3]

def GetAzm(x,y,data):
    'Give azimuth value for detector x,y position; calibration info in data'
    return GetTthAzmDsp(x,y,data)[1]




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
        fout[:] = out
        out = fout

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
    print errs, 'errs'

    #if errs:
    #    errs = errs
    # else:
    #    errs
    #print 'fuck', param
    #print param.shape
    #print stop
    popt,pcurv = curve_fit(Fourier_expand,azimu,ydata,p0=param,sigma=errs)

    #print 'popt'
    #print popt
    #print 'pcurv'
    #print pcurv

    return popt,pcurv


##warning about number of free parms e.g. if n is num. chunks/2


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
    '''
    d_0=d0s
    H_all=heights
    W_all=widths
    '''
    #print wavelength
    tth0 = 2.*np.degrees(np.arcsin(wavelength/2/d_0))

    # bg either single num. or size of azimu array
    #print tth0
    # print d_0,H_all,W_all,tth0,bg
    # print type(d_0),type(H_all),type(W_all),type(tth0),type(bg),type(twotheta),type(wavelength)
    #print d_0.dtype,H_all.dtype,W_all.dtype,tth0.dtype,twotheta.dtype
    Int = (H_all*np.exp((-(twotheta-tth0)**2)/(2*W_all**2))) + backg

    return Int



def singleInt(twotheta,*sparms):


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
    popt,pcurv = curve_fit(singleInt, twotheta, intens, p0=allparms,maxfev=8000)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv


def singleFit_test(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

    '''
    All expansion parameters fitted
    '''
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    # print type(d0s),d0s, intens.dtype
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    sparms = np.concatenate((d0s,heights,widths,bg), axis=None)
    # print type(d0s),d0s, intens.dtype
    popt,pcurv = curve_fit(singleInt, twotheta,intens,p0=sparms,maxfev=8000)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv



def singleFit_old(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

    '''
    All expansion parameters fitted
    '''
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = backg
    # print type(d0s),d0s, intens.dtype
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    # print type(d0s),d0s, intens.dtype
    popt,pcurv = curve_fit(singleInt, twotheta,intens,p0=[d0s,heights,widths, bg],maxfev=8000)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv



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
    Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + backg

    return Int



def dfit(intens,twotheta,azimu,dspace,d0s,heights,widths,wavelength,bg=None):

    '''
    Only d_0 fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
    heights=heights
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = 10.
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(dchange,(twotheta,azimu),intens,p0=d0s)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv


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
    Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + backg

    return Int



def hfit(intens,twotheta,azimu,newd0s,widths,hparms,wavelength,bg=None):

    '''
    Only h fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = 10.
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(hchange,(twotheta,azimu),intens,p0=hparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv



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
    Int = (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + backg

    return Int



def wfit(intens,twotheta,azimu,newd0s,newheights,wparms,wavelength,bg=None):

    '''
    Only w fitted.
    '''
    #d0s=tuple(d0s)
    # print d0s, 'd0s 1'
    #check for bg here, won't work if try to pass to fit
    if not bg:
        bg = 10.
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(wchange,(twotheta,azimu),intens,p0=wparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv



def Allchange(fullarray,*allparms):

    ##fit d,w,h as single number with no Fourier terms

    twothet,azi,parmnums = fullarray
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

    newall =  (H_all*np.exp((-(twothet-tth0)**2)/(2*W_all**2))) + backg

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
        bg = 10.
    intens=np.array(intens,dtype='float64')
    twotheta=np.array(twotheta,dtype='float64')
    azimu = np.array(azimu,dtype='float64')
    popt,pcurv = curve_fit(Allchange,(twotheta,azimu,parmnums),intens,p0=allparms)
    #popt,pcurv = curve_fit(testInt, twotheta,intens, p0=[d0s,heights,widths])
    return popt,pcurv







#### Main code ####



## Load files ##

# calibration parameter file #
parms_file = open(calib_file,'rb')
parms_dict = load_inputs_file(parms_file)
#print parms_dict
## may be a mask file?
## SAH: yes mask files -- they will either be a list of regions (x,y) or a binary image.
## SAH: I will generate you a mask file and send it to you.



# diffraction pattern #

im = Image.open(diff_file) ##always tiff?
## SAH: No. could be propreity formats. Possible formats for GSAS-II diffracion import are: tif, ADSC, cbf, png, edf, ge, hdf5, mar, rigaku, sfrm, 'gsas-ii image file' and guess.
## SAH: GSAS-II import files are located at https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/ with the file names 'G2img_*.py'.
#im.show()
imarray = np.array(im)
#print imarray.shape,imarray.max()

## Convert data ##

#The bits of the calibration that we need to convert X,Y to 2theta, azimuth or d,azimuth  are: Wavelength (Angstroms), distance (mm), center (X,Y, mm), tilt (degrees), rotation (degrees) and DetDepth (called penetration in screen shot).
#x = np.arange(imarray.shape[1])
#y = np.arange(imarray.shape[0])[0:imarray.shape[1]]

gd = np.mgrid[0:imarray.shape[0],0:imarray.shape[1]]       ##SAH: edit
print gd.shape
y = gd[0,:,:]       ##SAH: edit
x = gd[1,:,:]       ##SAH: edit
#print y.shape, x.shape
#print x[2]
y = y * pix / 1e3   ##SAH: edit
x = x * pix / 1e3   ##SAH: edit
#print x             ##SAH: edit
#print y             ##SAH: edit
'''
## SAH: would linspace+meshgrid be a better way of making the arrays?
## DMF: Could do it this way...
xn = np.linspace(0, imarray.shape[0],num=imarray.shape[0],endpoint=False,dtype=int)#nx)
yn = np.linspace(0, imarray.shape[1],num=imarray.shape[1],endpoint=False,dtype=int)#ny)
xv, yv = np.meshgrid(xn, yn)

#print xv.shape,yv.shape,xv[1],yv[1]

'''


## Convert to twotheta and azimuth arrays ##

twotheta = GetTth(x,y,parms_dict)
azimu =  GetAzm(x,y,parms_dict)
dspace = GetDsp(x,y,parms_dict)
intens = imarray




# print twotheta.shape,azimu.shape,dspace.shape,intens.shape

#checks#
# print twotheta.shape,azimu.shape,dspace.shape,twotheta[0]
colors = ("red", "green", "blue")
#plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
#plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
#plt.colorbar()
#plt.show()




## selected limits via gui - 2theta_min, 2theta_max, n_peaks (each peak 2theta min+max)
## for each peak....


## Example introduction of mask via intensity cut. Mask created
## for intens array then applied to all arrays.
# intens = ma.asarray(intens)
# intens = ma.masked_less(intens,0)
# #print azimu[~intens.mask]
# #print intens.shape
# #newazimu = ma.asarray(azimu)
# azimu = ma.array(azimu, mask=intens.mask)
# #print azimu.shape
# twotheta = ma.array(twotheta,mask=intens.mask)
# dspace = ma.array(dspace,mask=intens.mask)



# create mask from GSAS-II mask file.
intens = CreateGSASIIMask(mask_file, intens, gd.shape, twotheta, azimu, y, x)

#apply mask to other arrays
azimu    = ma.array(azimu, mask=intens.mask)
twotheta = ma.array(twotheta,mask=intens.mask)
dspace   = ma.array(dspace,mask=intens.mask)




# plot input file


#fig = plt.figure()
#plt.scatter(twotheta, azimu, s=4, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
##plt.scatter(dspace.flatten()[tthchunk],azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
#
#plt.colorbar()
#plt.show()
#
#plt.close()


## Split data into chunks in azimuth, replace azimuths with avg.
## doesn't have to be 100

print '\nGroup fitting in azimuth...\n'


azmax = azimu.max()
azmin = azimu.min()
binsize = (azmax-azmin)/total
chunks = []
azichunks = []
for i in xrange(total):
    end = (i+1)*binsize
    start = i*binsize
    temptth = twotheta.flatten()
    tthchunk = np.where((temptth>=tthRange[0])&(temptth<=tthRange[1]))
    tempazi = azimu.flatten()
    azichunk = np.where((tempazi>start)&(tempazi<=end))
    chunkels = np.intersect1d(tthchunk,azichunk)

    azichunks.append(((end-start)/2)+start)
    chunks.append(chunkels)


##### Final output list of azimuths with corresponding twotheta_0,h,w

wavelength=parms_dict['wavelength']
newd0=[]
newd0Err=[]
newHall=[]
newHallErr=[]
newWall=[]
newWallErr=[]
newBGall=[]
newBGallErr=[]

for j in range(len(chunks)):


    # print '\nFitting azimuth chunk....\n'

    # print intens.flatten()[chunks[j]].dtype
    hguess_I=intens.flatten()[chunks[j]].argmax()
    #print hguess_I
    hguess=intens.flatten()[chunks[j]][hguess_I]-backg[0]
    #print hguess
    dguess=dspace.flatten()[chunks[j]][hguess_I]
    bgguess = backg
    wguess = 0.005
    print hguess_I,hguess,dguess,backg,'hguess'
    po,pc=singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[j]],d0s=dguess,heights=hguess,widths=wguess,wavelength=np.float(parms_dict['wavelength']),bg=bgguess)
    #print po
    AllErr = np.sqrt(np.abs(np.diag(pc)))
    #print AllErr

    #print newd0Err,'newd0Err'
    newd0.append(po[0])
    newd0Err.append(AllErr[0])
    newHall.append(po[1])
    newHallErr.append(AllErr[1])
    newWall.append(po[2])
    newWallErr.append(AllErr[2])
    newBGall.append(po[3])
    newBGallErr.append(AllErr[3])

    print po, 'p0'

    '''
    newd0.append(po[0])
    newd0Err.append(np.sqrt(pc[0,0]))
    newHall.append(po[1])
    newHallErr.append(np.sqrt(pc[1,1]))
    newWall.append(po[2])
    newWallErr.append(np.sqrt(pc[2,2]))
    newBGall.append(po[3])
    newBGallErr.append(np.sqrt(pc[3,3]))
    '''

    # #plot the fits.
    # asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
    # asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
    # asdf2 = intens.flatten()[chunks[j]][asdf_ord]
    # asdf3 = singleInt(asdf1,tuple(po))
    # #asdf3 = singleInt(asdf1,po[0],po[1],po[2],po[3])


    # plt.plot(asdf1, asdf2,'.')
    # plt.plot(asdf1, asdf3, '-')
    # plt.show()




### Feed each d_0,h,w into fourier expansion function to get fit for fourier component
### parameters as output.

# print '\n Doing d0 fourier fit...', d0_order

#print 'newd0Err', newd0Err

dfour = Fourier_fit(np.array(azichunks),np.array(newd0),terms=d0_order, errs=np.array(newd0Err))
# print '\n Doing h fourier fit...'
hfour = Fourier_fit(np.array(azichunks),np.array(newHall),terms=h_order, errs=np.array(newHallErr))
wfour = Fourier_fit(np.array(azichunks),np.array(newWall),terms=w_order, errs=np.array(newWallErr))

'''
print 'shapes'
print dfour[0].shape
print hfour[0].shape
print wfour[0].shape
#print 'hfour'
#print hfour
#print 'wfour'
#print wfour
'''

#plot output of fourier fits....

# fig = plt.figure()
# ax = fig.add_subplot(3,1,1)
# plt.subplot(311)
# plt.plot(azichunks,newd0, 'bo')
# plt.title('D-spacing')
# # print len(azichunks), 'len azichunks'
# print dfour[0]
# testout = Fourier_expand(np.array(azichunks), dfour[0])
# print testout, 'testout'
# #print stop
# plt.plot(azichunks,Fourier_expand(np.array(azichunks), dfour[0]), 'r-')
# plt.subplot(312)
# plt.plot(azichunks,newHall, 'bo')
# plt.title('height')
# plt.plot(azichunks,Fourier_expand(np.array(azichunks), hfour[0]), 'r-')
# plt.subplot(313)
# plt.plot(azichunks,newWall, 'bo')
# plt.title('width')
# plt.plot(azichunks,Fourier_expand(np.array(azichunks), np.array(wfour[0])), 'r-')
# plt.show()

# plt.close()

# print stop


### Feed each d,h,w into versions of main equation to fit only or one at a time
### including the new paramters from the Fourier fits above as initial inputs
## Slice arrays for twothetarange

print '\nRe-fitting for d, h, w separately...\n'


temptth = twotheta.flatten()
tthchunk = np.where((temptth>=tthRange[0])&(temptth<=tthRange[1]))
#print 'check point'

## Fit for d0 with existing h and w

heights = Fourier_expand(azimu.flatten()[tthchunk],hfour[0])
widths = Fourier_expand(azimu.flatten()[tthchunk],wfour[0])
#print dfour[0]
newdfour = dfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],dspace.flatten()[tthchunk],dfour[0],heights,widths,wavelength,bg=None)
print 'Old d coefficients: ', dfour[0]
print 'New d coefficients: ', newdfour[0]

## Fit for h with new d and existing w

newd0s = Fourier_expand(azimu.flatten()[tthchunk],newdfour[0])
widths = Fourier_expand(azimu.flatten()[tthchunk],wfour[0])
#print hfour[0]
newhfour = hfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],newd0s,widths,hfour[0],wavelength,bg=None)
print 'Old h coefficients: ', hfour[0]
print 'New h coefficients: ', newhfour[0]

## Fit for w with new d and h

newd0s = Fourier_expand(azimu.flatten()[tthchunk],newdfour[0])
newheights = Fourier_expand(azimu.flatten()[tthchunk],newhfour[0])
#print wfour[0]
newwfour = wfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],newd0s,newheights,wfour[0],wavelength,bg=None)
print 'Old w coefficients: ', wfour[0]
print 'New w coefficients: ', newwfour[0]



#stop


### Refit through full equation with all data for d,h,w independently

print '\nFinal fit solving for all parms...\n'

newwidths = Fourier_expand(azimu.flatten()[tthchunk],newwfour[0])

Finalparms = Allfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],newdfour[0],newhfour[0],newwfour[0],wavelength,bg=None)

starth = len(newdfour[0])
startw = len(newdfour[0])+len(newhfour[0])
end = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])

fin_d = Finalparms[0][0:starth]
fin_h = Finalparms[0][starth:startw]
fin_w = Finalparms[0][startw:end]


#print len(fin_d)

print 'Final d0 coefficients: ', Finalparms[0][0:starth]
print 'Final h coefficients: ', Finalparms[0][starth:startw]
print 'Final w coefficients: ', Finalparms[0][startw:end]





### Plot results to check

print '\nPlotting results for fit...\n'

tthchunk = np.where((twotheta.flatten()>=tthRange[0])&(twotheta.flatten()<=tthRange[1]))
# print Finalparms[0]
fit_intens = Allchange((twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],np.array([len(fin_d),len(fin_h),len(fin_w)])), *Finalparms[0].squeeze())
fullfit_intens = np.ones((intens.shape[0],intens.shape[1]),float)
oldshape = fullfit_intens.shape
fullfit_intens = fullfit_intens.flatten()
fullfit_intens[tthchunk] = fit_intens
#print fullfit_intens.max()
#print fullfit_intens.shape
#print oldshape
fullfit_intens = np.reshape(fullfit_intens,(oldshape[0],oldshape[1]))

ValMax1 = np.max(intens.flatten()[tthchunk])
ValMax2 = np.max(fullfit_intens.flatten()[tthchunk])
ValMax3 = np.max((intens.flatten()[tthchunk]-fullfit_intens.flatten()[tthchunk]))
ValMax  = 40#np.max([ValMax1, ValMax2, ValMax3])
ValMin1 = np.min(intens.flatten()[tthchunk])
ValMin2 = np.min(fullfit_intens.flatten()[tthchunk])
ValMin3 = np.min((intens.flatten()[tthchunk]-fullfit_intens.flatten()[tthchunk]))
ValMin  = -3#np.min([ValMin1, ValMin2, ValMin3])

#print fullfit_intens.shape
#print twotheta.shape
#print azimu.shape
#print intens.shape
#print tthchunk
#print twotheta[0,]
#print np.where((twotheta[0]>=tthRange[0])&(twotheta[0]<=tthRange[1]))

if 1:
    fig = plt.figure()
    ax = fig.add_subplot(1,3,1)
    plt.subplot(131)
    #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
    plt.scatter(dspace.flatten()[tthchunk],azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
    plt.subplot(132)
    #plt.scatter(twotheta, azimu, s=1, c=np.log(fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
    plt.scatter(dspace.flatten()[tthchunk], azimu.flatten()[tthchunk], s=4, c=(fullfit_intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
    plt.colorbar()
    plt.subplot(133)
    #plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
    plt.scatter(dspace.flatten()[tthchunk], azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]-fullfit_intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet)
    plt.colorbar()
    plt.show()

    plt.close()



    # print '\nFitting azimuth chunk....\n'

  #  # print intens.flatten()[chunks[j]].dtype
  #  hguess_I=intens.flatten()[chunks[j]].argmax()
  #  hguess=intens.flatten()[chunks[j]][hguess_I]
  #  dguess=dspace.flatten()[chunks[j]][hguess_I]
  #  #print hguess_I,hguess,dguess
  #  po,pc=singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[#j]],d0s=dguess,heights=hguess,widths=0.005,wavelength=np.float(parms_dict['wavelength']),bg=10.)
  #  print po,p#c
  #  newd0.append(po[0])
  #  newHall.append(po[1])
  #  newWall.append(po[2])

  #  #plot the fits.
  #  asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
  #  asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
  #  asdf2 = intens.flatten()[chunks[j]][asdf_ord]
  #  asdf3 = singleInt(asdf1,po[0],po[1],po[2])

  #  plt.plot(asdf1, asdf2,'x--')
   # plt.plot(asdf1, asdf3, '-')
  #  plt.show()












## Specify ranges to search/starting parameters ##



#Input from GUI/file?
## SAH: both. I think I would have two functions.
## SAH: A GUI setup function and a wrapper that calls this function with a parameter file or if it does not exist then calls the GUI function.

#Load plot to screen and choose regions to fit? ## SAH: yes.
#can then slice arrays and calcualte intial guesses from peaks and avgs.



## Setup/Choose function to fit ##


## make arrays. -- not line arrays n by 2. possibly better as a line array becuase then we are not carrying around an unused number.
#heights = ones(2,h_order) ## better to be intellegent guesses
#widths = ones(2,w_order)  ## better to be intellegent guesses
#d0s = ones(2,d0_order)    ## better to be intellegent guesses


## Now have starting model


## refit starting model via full intensity equation but hold other parameters fixed e.g. fit for d, not h,w,
## do for each in turn i.e. if strat with d, then h, then w


## then continue fitting with full equation, varying all parms.





'''
#start_model
Segregate into chunks in azimuth
get avg. azimuth
flatten intens and twotheta arrays to pass to fit
pass to gauss fit to get first approx of 2th_0,h,w with nterms =0 for w,h,2th for each block

#start_fourier
take array of 2th_0,h,w over azimu blocks
choose nterms for each for now
fit each with Fourier expansion of nterms

#main_fit
use all data
take fourier expansion with same nterms and input parms from start_fourier
send to allfit to allow all parms to fit

#inside main fit
redo 2th_0/d_0,h,w fit iteratively
then input to allfit to fit parms in gauss function
increase nterms?

generally think will have to pass flattened azimu,2th and intens arrays.

'''
