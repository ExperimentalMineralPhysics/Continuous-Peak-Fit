#!/usr/bin/env python
# -*- coding: utf-8 -*-


__all__ = ['Requirements', 'ImportImage', 'DetectorCheck', 'Conversion', 'GetMask', 'GetCalibration', 'GetTth',
           'GetAzm', 'GetDsp']

import sys
import os
import h5py
import numpy as np
import numpy.ma as ma
from PIL import Image
from PIL.ExifTags import TAGS
import matplotlib.pyplot as plt

from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import fabio


# For Dioptas functions to change
#   -   Load/import data
#   -   Load mask
#   -   Load Calibration

def Requirements():
    
    RequiredParams = [
            'Calib_param',
            'Calib_detector',
            #'Calib_pixels'
            ]
    
    return RequiredParams


def ImportImage(ImageName, debug=False):
    # im = Image.open(ImageName) ##always tiff?- no

    filename, file_extension = os.path.splitext(ImageName)

    if file_extension == '.nxs':
        im_all = h5py.File(ImageName, 'r')
        # FIX ME: This assumes that there is only one image in the nxs file. If there are more then it will only read
        # 1.
        im = np.array(im_all['/entry1/instrument/detector/data'])
    else:
        im_all = fabio.open(ImageName)
        im = im_all.data

    # Dioptas flips the images to match the orientations in Plot2D. 
    # Therefore implemented here to be consistent with Dioptas.
    im = np.array(im)[::-1]
    
    if debug:
        imgplot = plt.imshow(im, vmin=0, vmax=2000)
        plt.colorbar()
        plt.show()

    #    info = i._getexif()
    #    for tag, value in info.items():
    #        decoded = TAGS.get(tag, tag)
    #        ret[decoded] = value
    #
    #    print ret
    #    stop

    # the bits of Fabio that might get most of the parameters not saved in calibration
    # print(im_all.DESCRIPTION)
    # print(im_all.factory)
    # print(im_all.header)
    # print(im_all.header_keys)

    return im


def DetectorCheck(ImageName, detector=None):
    # FIX ME: we should try importing the image and reading the detector type from it.
    if detector is None:
        sys.exit('\n\nDioptas requires a detector type.\nTo list the possible detector types in the command line '
                 'type:\n   import pyFAI\n   pyFAI.detectors.Detector.registry\n\n')
    # FIX ME: check the detector type is valid.

    return detector


def Conversion(tth_in, conv, reverse=0):
    # convert two theta values into d-spacing.

    # convert wavelength into angstroms.
    # this is not longer required because it is done in the GetCalibration stage. 
    wavelength = conv['conversion_constant']
    # print(wavelength)

    if not reverse:
        # convert tth to d_spacing
        dspc_out = wavelength / 2 / np.sin(np.radians(tth_in[:, 1] / 2))

    else:
        # convert d-spacing to tth.
        # N.B. this is the reverse function so that labels tth and d_spacing are not correct.
        # print(tth_in)
        dspc_out = 2 * np.degrees(np.arcsin(wavelength / 2 / tth_in))
        # dspc_out = 2*np.degrees(np.arcsin(wavelength/2/tth_in[:,1]))

    return dspc_out


def GetMask(MSKfile, ImInts, ImTTH, ImAzi, Imy, Imx, debug=False):
    # Dioptas mask is compressed Tiff image. 
    # Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py

    ImMsk = np.array(Image.open(MSKfile))
    ImInts = ma.array(ImInts, mask=ImMsk)

    if debug:
        # Plot mask.
        # This is left in here for debugging.
        fig = plt.figure()
        ax = fig.add_subplot(1, 2, 1)
        plt.subplot(121)
        plt.scatter(Imx, Imy, s=4, c=intens, edgecolors='none', cmap=plt.cm.jet)
        ax = fig.add_subplot(1, 2, 2)
        plt.subplot(122)
        plt.scatter(ma.array(Imx, mask=ImMsk.mask), ma.array(Imy, mask=ImMsk.mask), s=4, c=intens, edgecolors='none',
                    cmap=plt.cm.jet)
        plt.colorbar()
        plt.show()

        plt.close()

    #FIX ME: need to validate size of images vs. detector name. Otherwise the mask can be the wrong size 
    #det_size = pyFAI.detectors.ALL_DETECTORS['picam_v1'].MAX_SHAPE
    #FIX ME : this could probably all be done by using the detector class in fabio. 


    return ImInts


def GetCalibration(filenam):
    """Parse inputs file, create type specific inputs."""

    parms_file = open(filenam, 'r')

    filelines = parms_file.readlines()
    parms_dict = {}

    for item in filelines:
        newparms = item.strip('\n').split(':', 1)
        parm = newparms[1].strip()

        value = None
        # print parm
        try:
            value = int(parm)
            parms_dict[(str(newparms[0]).lower())] = value
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
                                newlist.append(val.replace("'", "").replace(" ", "").replace("\"", ""))
                    parms_dict[newparms[0]] = newlist
                elif parm.startswith('{'):
                    # print parm
                    listvals = parm.strip('{').strip('}').split(',')
                    newdict = {}
                    for keyval in listvals:
                        # print keyval
                        newkey = keyval.split(':')[0].replace("'", "").replace(" ", "").replace("\"", "")
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
                                newdict[str(newkey)] = val.replace("'", "").replace(" ", "").replace("\"", "")
                    parms_dict[newparms[0]] = newdict
                elif not parm:
                    parms_dict[newparms[0]] = ''

                else:
                    parms_dict[newparms[0]] = str(parm)

    #process 'detector_config' into pixel sizes if needed.
    if 'Detector_config' in parms_dict:
        print(parms_dict['Detector_config'])
        parms_dict['pixelsize1'] = parms_dict['Detector_config']['PixelSize1']
        parms_dict['pixelsize2'] = parms_dict['Detector_config']['PixelSize2']
                   
    # get wavelengths in Angstrom
    parms_dict['conversion_constant'] = parms_dict['Wavelength'] * 1E10

    # force all dictionary labels to be lower case -- makes s
    parms_dict = {k.lower(): v for k, v in parms_dict.items()}

    # report type
    parms_dict['DispersionType'] = 'AngleDispersive'

    # print(parms_dict)
    return parms_dict


def GetTth(cali_file, poni, pix=None, det=None):
    'Give 2-theta value for detector x,y position; calibration info in data'
    ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'], poni['rot3'],
                             pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector=det,
                             wavelength=poni['wavelength'])
    # Tth = ai.twoThetaArray()
    # print('Tth', np.min(Tth)/np.pi*180, np.max(Tth)/np.pi*180)
    return np.degrees(ai.twoThetaArray())


def GetAzm(cali_file, poni, pix=None, det=None):
    'Give azimuth value for detector x,y position; calibration info in data'
    ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'], poni['rot3'],
                             pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector=det,
                             wavelength=poni['wavelength'])
    # Chi = ai.chiArray()
    # print('Chi', np.min(Chi)/np.pi*180, np.max(Chi)/np.pi*180)
    return np.degrees(ai.chiArray())


def GetDsp(cali_file, poni, pix=None, det=None):
    'Give d-spacing value for detector x,y position; calibration info in data'
    'Get as thetheta and then convert into d-spacing'
    ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'], poni['rot3'],
                             pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector=det,
                             wavelength=poni['wavelength'])

    # dspc = poni['wavelength']*1E10/2/np.sin(ai.twoThetaArray()/2)
    # print('dspc', np.min(dspc), np.max(dspc))
    return poni['wavelength'] * 1E10 / 2 / np.sin(ai.twoThetaArray() / 2)
