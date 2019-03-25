# /usr/bin/python
# -*- coding: utf-8 -*-

import numpy as np
import numpy.ma as ma
from PIL import Image
from PIL.ExifTags import TAGS
import matplotlib.pyplot as plt

from pyFAI.azimuthalIntegrator import AzimuthalIntegrator


# For Dioptas functions to change
#   -   Load/import data
#   -   Load mask
#   -   Load Calibration



def ImportImage(ImageName):

    im = Image.open(ImageName) ##always tiff?- no
     
    # Dioptas flips the images to match the orientations in Plot2D. 
    #Therefore implemented here to be consistent with Dioptas.
    im = np.array(im)[::-1]
    
#    info = i._getexif()
#    for tag, value in info.items():
#        decoded = TAGS.get(tag, tag)
#        ret[decoded] = value
#
#    print ret
#    stop
    
    
    return im




def GetMask(MSKfile, ImInts, ImTTH, ImAzi, Imy, Imx):
    # Dioptas mask is compressed Tiff image. 
    # Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py
    
    ImMsk = np.array(Image.open(MSKfile))


    if 0:
        #Plot mask.
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
    
    
    ImInts = ma.array(ImInts, mask=ImMsk)
    
    
    return ImInts




def GetCalibration(filenam):
    'Parse inputs file, create type specific inputs.'

    parms_file = open(filenam,'rb')

    filelines = parms_file.readlines()
    parms_dict = {}

    for item in filelines:
        newparms = item.strip('\n').split(':',1)
        parm = newparms[1]

        value = None
        # print parm
        try:
            value = int(parm)
            parms_dict[lower(str(newparms[0]))] = value
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

    # force all dictionary labels to be lower case -- makes s
    parms_dict =  {k.lower(): v for k, v in parms_dict.items()}

    #print parms_dict
    return parms_dict






def GetTth(cali_file,poni, pix=None):
    'Give 2-theta value for detector x,y position; calibration info in data'   
    ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'], poni['rot3'], pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector='PILATUS1M', wavelength=poni['wavelength']) 
    #Tth = ai.twoThetaArray() 
    #print('Tth', np.min(Tth)/np.pi*180, np.max(Tth)/np.pi*180)
    return np.degrees(ai.twoThetaArray())


def GetAzm(cali_file,poni, pix=None):
    'Give azimuth value for detector x,y position; calibration info in data'
    ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'], poni['rot3'], pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector='PILATUS1M', wavelength=poni['wavelength']) 
    #Chi = ai.chiArray()
    #print('Chi', np.min(Chi)/np.pi*180, np.max(Chi)/np.pi*180)
    return np.degrees(ai.chiArray())


def GetDsp(cali_file,poni, pix=None):
    'Give d-spacing value for detector x,y position; calibration info in data'
    'Get as thetheta and then convert into d-spacing'
    ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'], poni['rot3'], pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector='PILATUS1M', wavelength=poni['wavelength']) 
    #dspc = poni['Wavelength']/2/np.sin(ai.twoThetaArray())
    #print('dspc', np.min(dspc)/np.pi*180, np.max(dspc)/np.pi*180)
    return poni['wavelength']/2/np.sin(ai.twoThetaArray())




