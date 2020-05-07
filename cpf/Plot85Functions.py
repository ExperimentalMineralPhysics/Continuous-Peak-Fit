#!/usr/bin/env python

#import sys
#sys.path.append('/Users/simon/g2conda/GSASII/')

__all__ = ['ImportImage', 'CreateGSASIIMask', 'pointInPolygon', 'Fill2ThetaAzimuthMap', 'load_inputs_file', 'LoadMask']

import numpy as np
import numpy.ma as ma
import copy, os, sys
from PIL import Image
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.path as mlp
from scipy.optimize import curve_fit
np.set_printoptions(threshold=sys.maxsize)


#import Plot85 as GsIO

# ------- above here is the intro stuff that I dont know if it is needed. 


def ImportImage(ImageName):


    #Im = GsIO.GetImageData([], ImageName)



    File = open(ImageName,'r')
    save = {}
    #oldThreshold = data['Thresholds'][0]
    S = File.readline()
    while S:
        # if the row is 10 numbers then keep otherise not data.
        if S == 'DATA:':
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




    im = Image.open(ImageName) ##always tiff?- no
    # print type(im)
    #im = im.transpose(Image.FLIP_TOP_BOTTOM)
    #im = im.transpose(Image.FLIP_LEFT_RIGHT)
    #print type(im)
    return im

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
    # print FrmLms
    #for frame in FrmLms:
    #print frame
    if FrmLms:
        path = mlp.Path(FrmLms)
        grid = path.contains_points(points)-True
        print(Imx.shape[0])
        print(Imx.shape[1])
        print(Imx.shape[0]*Imx.shape[1])
        print(grid.shape)
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

    filelines = filenam.readlines()
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


def peneCorr(tth,dep,tilt=0.,azm=0.):
    'Needs a doc string'
    return dep*(1.-npcosd(tth))         #best one


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



