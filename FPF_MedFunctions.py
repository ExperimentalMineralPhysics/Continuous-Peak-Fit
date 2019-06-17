# /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import h5py
import numpy as np
import numpy.ma as ma
from PIL import Image
from PIL.ExifTags import TAGS
import matplotlib.pyplot as plt
import Med
import Mca

# Med and related functions require the 'Numeric' and 'LinearAlgebra' packages.
# I don't know how much of this is required but for now I am just adding the packages.
# FIX ME: change to call numpy and other alternatives (if needed). 

# FIX ME: Need a dependency check.  

# For MED files from Energy dispersive diffraction patterns at X17B2, 6-BMB and others.


def ImportImage(ImageName):

    filename, file_extension = os.path.splitext(ImageName)
    
    if file_extension == '.med':
        
        dat = Med.Med()
        dat.read_file(ImageName)
        
        im_all = []
        for x in range(dat.n_detectors):
            im_all.append(dat.mcas[x].get_data())
                    
        #print im_all
        #print len(im_all)
        
    else:
        print 'Unknown file type'
        stop
    
    return im_all


    
def DetectorCheck(ImageName, detector=None):
    
    
    # FIX ME: we should try importing the image and reading the detector type from it.
    if detector == None:
        
        sys.exit('\n\nDioptas requires a detector type.\nTo list the possible detector types in the command line type:\n   import pyFAI\n   pyFAI.detectors.Detector.registry\n\n')
    # FIX ME: check the detector type is valid.

    return detector


def Conversion(tth_in, calib, reverse=0):
    #convert energy into d-spacing.
    
    print calib

    stop    
    
    #convert wavelength into angstroms.
    # this is not longer required because it is done in the GetCalibration stage. 
    wavelength = wavelength
    
    if not reverse:
        #convert tth to d-spacing
        dspc_out = wavelength/2/np.sin(np.radians(tth_in/2))
        
    else:
        #convert dspacing to tth.
        #N.B. this is the reverse finction so that labels tth and dspacing are not correct. 
        dspc_out = 2*np.degrees(np.arcsin(wavelength/2/tth_in))
        
        
    return dspc_out



def GetMask(MSKfile, ImInts, ImTTH, ImAzi, Imy, Imx):
    # Masks for 10 element detector are either for full detector or for an energy range. 
    # Currently this just works by removing  single detector.
    
    ImMsk = MSKfile

    ImMsk = ma.zeros(ImInts.shape)
    for x in range(len(MSKfile)):
        ImMsk[MSKfile[x]-1] = 1
    
    if 0:
        #Plot mask.
        #This is left in here for debugging.
        fig = plt.figure()
        ax = fig.add_subplot(1,2,1)
        plt.subplot(121)
        plt.scatter(ImTTH, ImAzi, s=4, c=ImInts, edgecolors='none', cmap=plt.cm.jet)
        ax = fig.add_subplot(1,2,2)
        plt.subplot(122)
        plt.scatter(ma.array(ImTTH,mask=ImMsk), ma.array(ImAzi,mask=ImMsk), s=4, c=ImInts, edgecolors='none', cmap=plt.cm.jet)
        plt.colorbar()
        plt.show()
    
        plt.close()
    
    
    ImInts = ma.array(ImInts, mask=ImMsk)
    
    
    return ImInts




def GetCalibration(filenam):
    'Parse inputs file, create type specific inputs.'

#    asdf = Med.Med()
#    asdf.read_file(filenam)
#    parms_dict =  asdf.get_calibration
#    print parms_dict
    parms_dict = {}


    dat = Med.Med(file=filenam)

    parms_dict['DispersionType']   = 'EnergyDispersive'
    parms_dict['calibs'] = dat
    
    if len(dat.mcas) == 10:
        
        az = np.linspace(0,180,9)
        az = np.append(az, 270)
        parms_dict['azimuths'] = az
        
    else:
        sys.exit('\n\nDetector is unknown. Edit function to add number of detectors and associated azimuths.\n\n')
  
        
    
#    filename, file_extension = os.path.splitext(filenam)
#    if file_extension == '.med':
#        t = Mca.read_ascii_file(filenam)
#        parms_dict['dets'] = t['calibration']
#        
#               
#        #define azimuths from number of detectors. This has the be assumed
#        if len(parms_dict['dets'])==10:
#            #assume 10 element detector from NSLS--X17B2 or APS--6-BMB
#        
#            l = np.shape(t['data'])
#            
#            azi = np.linspace(0,180,9)
#            azi = np.append(azi, 270)
#            #azi = np.tile(azi, [l[1],1]).T
#        
#        else:
#            print 'Unknown detector type'
#            stop
#    
#        
#        parms_dict['data_length'] = len(t['data'])
#        parms_dict['azi'] = azi
#        
#        conv = []
#        for x in range(len(parms_dict['dets'])):
#            conv.append([azi[x],parms_dict['dets'][x].two_theta]) 
#        
#        parms_dict['conversion_constant'] = conv
#        print parms_dict['conversion_constant'] 
#        stop
#    else:
#        print 'Unknown file type'
#        stop

    
#    print type(asdf)
#    print np.shape(asdf['data'])
#    print asdf['calibration'][0].two_theta
    

    
    return parms_dict






def GetTth(cali_file, parms_dict, pix=None, det=None):
    'Called get two theta but it is really getting the energy.'
    
    #dat = Med.Med()
    #dat.read_file(cali_file)
    channels = np.arange(parms_dict['calibs'].mcas[0].data.shape[0])
    
    en = []
    for x in range(parms_dict['calibs'].n_detectors):
        en.append(parms_dict['calibs'].mcas[x].channel_to_energy(channels))
        
    #print np.shape(en)
    return en


def GetAzm(cali_file, parms_dict, pix=None, det=None):
    'Give azimuth value for detector x,y position; have to assume geometry of the detectors'
    
    
    l = parms_dict['calibs'].mcas[0].data.shape
        
    azi = np.tile(parms_dict['azimuths'], [l[0],1]).T
    
#        
#    if len(poni)==10:
#        #assume 10 element detector from NSLS--X17B2 or APS--6-BMB
#    
#        l = np.shape(t['data'])
#        
#        azi = np.linspace(0,180,9)
#        azi = np.append(azi, 270)
#        azi = np.tile(azi, [l[1],1]).T
#    
#    else:
#        print 'Unknown file type'
#        stop
    
    return azi


def GetDsp(cali_file, parms_dict, pix=None, det=None):
    'Give d-spacing value for detector x,y position; calibration info in data'
        
#    dat = Med.Med()
#    dat.read_file(cali_file)
#    channels = np.arange(len(dat.data))
    channels = np.arange(parms_dict['calibs'].mcas[0].data.shape[0])
    
    dspc = []
    for x in range(parms_dict['calibs'].n_detectors):
        dspc.append(parms_dict['calibs'].mcas[x].channel_to_d(channels))
        
    #print np.shape(en)
    
    return dspc
