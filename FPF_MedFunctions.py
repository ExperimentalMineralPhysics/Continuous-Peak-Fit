# /usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import Med
import Mca

# For MED files from Energy dispersive diffraction patterns at X17B2, 6-BMB and others.

# For ease of use the *.med files are imported using Mark Rivers' Med.py and Mca.py functions.
# These were downloaded from ??? in June 2019.
# and edited to remove unnecessary dependencies e.g. 'Numeric' and 'LinearAlgebra' packages and functions that had further dependencies.
# I don't know how much of this is required but for now I am just adding the packages.

# FIX ME: Need a dependency check.  


def ImportImage(ImageName):

    filename, file_extension = os.path.splitext(ImageName)
    
    if file_extension == '.med':
        dat = Med.Med()
        dat.read_file(ImageName)
        im_all = []
        for x in range(dat.n_detectors):
            im_all.append(dat.mcas[x].get_data())
        
    else:
        sys.exit('Unknown file type. A *.med file is required')
    
    return im_all


    
def DetectorCheck(ImageName, detector=None):
    
    
    # FIX ME: we should try importing the image and reading the detector type from it.
    if detector == None:
        
        sys.exit('\n\nDioptas requires a detector type.\nTo list the possible detector types in the command line type:\n   import pyFAI\n   pyFAI.detectors.Detector.registry\n\n')
    # FIX ME: check the detector type is valid.

    return detector


def Conversion(tth_in, calib, reverse=0):
    #convert energy into d-spacing.

    
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

    return parms_dict



def GetTth(cali_file, parms_dict, pix=None, det=None):
    'Called get two theta but it is really getting the energy.'
    
    channels = np.arange(parms_dict['calibs'].mcas[0].data.shape[0])
    en = []
    for x in range(parms_dict['calibs'].n_detectors):
        en.append(parms_dict['calibs'].mcas[x].channel_to_energy(channels))
        
    return en



def GetAzm(cali_file, parms_dict, pix=None, det=None):
    'Give azimuth value for detector x,y position; have to assume geometry of the detectors'
    
    l = parms_dict['calibs'].mcas[0].data.shape
    azi = np.tile(parms_dict['azimuths'], [l[0],1]).T
    
    return azi



def GetDsp(cali_file, parms_dict, pix=None, det=None):
    'Give d-spacing value for detector x,y position; calibration info in data'
        
    channels = np.arange(parms_dict['calibs'].mcas[0].data.shape[0])
    dspc = []
    for x in range(parms_dict['calibs'].n_detectors):
        dspc.append(parms_dict['calibs'].mcas[x].channel_to_d(channels))
        
    return dspc
