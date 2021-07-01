#!/usr/bin/env python
# -*- coding: utf-8 -*-


__all__ = ['DioptasDetector']
# __all__ = ['Requirements', 'ImportImage', 'DetectorCheck', 'Conversion', 'GetMask', 'GetCalibration', 'GetTth',
# 'GetAzm', 'GetDsp']

import sys
import os
import h5py
import numpy as np
import numpy.ma as ma
from PIL import Image
from PIL.ExifTags import TAGS
import matplotlib.pyplot as plt
from matplotlib import gridspec

from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import pyFAI
import fabio


class DioptasDetector:
    # For Dioptas functions to change
    #   -   Load/import data
    #   -   Load mask
    #   -   Load Calibration

    def __init__(self, calibration_parameters, fit_parameters=None):
        """
        :param calibration_parameters:
        """
        self.parameters = self.get_calibration(calibration_parameters)
        self.requirements = self.get_requirements(fit_parameters)
        self.detector = None
        self.intensity = None
        self.tth = None
        self.azm = None
        self.dspace = None

    def fill_data(self, calibration_data, detector=None, debug=None, calibration_mask=None):
        """
        :param calibration_data:
        :param detector:
        :param debug:
        :param calibration_mask:
        """
        # self.detector = self.detector_check(calibration_data, detector)
        self.get_detector(calibration_data, detector)
        self.intensity = self.get_masked_calibration(calibration_data, debug, calibration_mask)
        self.tth = self.GetTth(self.parameters, self.detector, self.intensity.mask)
        self.azm = self.GetAzm(self.parameters, self.detector, self.intensity.mask)
        self.dspace = self.GetDsp(self.parameters, self.detector, self.intensity.mask)

    def get_detector(self, calibration_data, detector=None):
        self.detector = self.detector_check(calibration_data, detector)

    def get_masked_calibration(self, calibration_data, debug, calibration_mask=None):
        """
        load calibration data/image -- to get intensities for mask
        :param calibration_data:
        :param calibration_pixels:
        :param debug:
        :param calibration_mask:
        :return:
        """
        im = self.ImportImage(calibration_data, debug=debug)
        intens = ma.array(im)
        # create mask from mask file if present. If not make all values valid
        if calibration_mask:
            intens = self.GetMask(calibration_mask, intens)
        else:
            ImMask = np.zeros_like(intens)
            intens = ma.array(intens, mask=ImMask)
            intens = ma.masked_outside(intens, 0, np.inf)
        return intens

    def get_requirements(self, parameter_settings=None):
        """
        Get the parameters required for this detector
        :return: String parameters
        """
        # List all the required parameters.
        # NOTE: No doubt this will change hence denoted as a list that is worked over.
        # 1. Have option to have list of files rather than sequence. This list needs to be called diff_files.

        # Inputs:
        #   calibration type. (required)
        #   data - in form of:
        #       - file listing data files. (with or without directory)
        #       - file listing data file numbers (needs directory, file name, number digits, ending)
        #       - beginning and and numbers for files (needs directory, file name, number digits, ending)
        RequiredList = ['Calib_type',
                        'Calib_detector'
                        # 'Calib_data',        # Removed because this is not strictly required.
                        'Calib_param',       # Now added to calibration function file.
                        # 'Calib_pixels',      # this is only needed for GSAS-II and should be read from image file. FIX
                        # ME: add to GSAS-II inputfile.
                        # 'Calib_mask',        # a mask file is not strictly required.
                        'datafile_directory',
                        'datafile_Basename',
                        'datafile_Ending',
                        # 'datafile_StartNum',  # now optionally replaced by datafile_Files
                        # 'datafile_EndNum',    # now optionally replaced by datafile_Files
                        'datafile_NumDigit',
                        'AziBins',            # required based on detector type
                        'fit_orders',
                        # 'Output_type',		   # should be optional
                        # 'Output_NumAziWrite',  # should be optional
                        # 'Output_directory']	   # should be optional
                        ]

        # Check required against inputs if given
        if parameter_settings is not None:
            # properties of the data files.
            all_present = 1
            for par in parameter_settings:
                if par in RequiredList:
                    print('Got: ', par)
                else:
                    print("The settings file requires a parameter called  \'", par, "\'")
                    all_present = 0
            if all_present == 0:
                sys.exit("The highlighted settings are missing from the input file. Fitting cannot proceed until they "
                         "are all present.")
        return RequiredList

    @staticmethod
    def ImportImage(image_name, mask=None, debug=False):
        """
        Import the data image
        :param mask:
        :param image_name: Name of file
        :param debug: extra output for debugging - plot
        :return: intensity array
        """
        # FIX ME: why is mask in here? should it ingerit from previous data if it exists?
        # im = Image.open(ImageName) ##always tiff?- no
        filename, file_extension = os.path.splitext(image_name)
        if file_extension == '.nxs':
            im_all = h5py.File(image_name, 'r')
            # FIX ME: This assumes that there is only one image in the nxs file.
            # If there are more then it will only read 1.
            im = np.array(im_all['/entry1/instrument/detector/data'])
        else:
            im_all = fabio.open(image_name)
            im = im_all.data
        # Dioptas flips the images to match the orientations in Plot2D.
        # Therefore implemented here to be consistent with Dioptas.
        im = np.array(im)[::-1]
        if debug:
            img_plot = plt.imshow(im, vmin=0, vmax=2000)
            plt.colorbar()
            plt.show()
        if mask is not None:
            return ma.array(im, mask=mask)
        else:
            return ma.array(im)

    @staticmethod
    def detector_check(calibration_data, detector=None):
        """
        Get detector information
        :param detector: Choose detector on initiation
        :param calibration_data:
        :param detector:
        :return:
        """
        # FIX ME: we should try importing the image and reading the detector type from it.
        # The image is now imported and read but this means that the detector type in settings
        # is no longer necessary and could be removed.
        if detector is None or detector == 'unknown' or detector == 'other' or detector == 'blank':
            im_all = fabio.open(calibration_data.Calib_data)
            sz = calibration_data.Calib_pixels  # Pixel_size
            if sz > 1:
                sz = sz*1e-6
            detector = pyFAI.detectors.Detector(pixel1=sz, pixel2=sz, splineFile=None, max_shape=im_all.shape)
        # FIX ME: check the detector type is valid.
        return detector

    def conversion(self, tth_in, conv, reverse=0, azm=None):
        """
        Convert two theta values into d-spacing.
        azm is needed to enable compatibility with the energy dispersive detectors
        :param tth_in:
        :param conv:
        :param reverse:
        :param azm:
        :return:
        """
        # FIX ME: Still needed?
        # Convert wavelength into angstroms.
        # this is not longer required because it is done in the GetCalibration stage.
        wavelength = conv['conversion_constant']
        # print(wavelength)
        if not reverse:
            # convert tth to d_spacing
            dspc_out = wavelength / 2 / np.sin(np.radians(tth_in / 2))
        else:
            # convert d-spacing to tth.
            # N.B. this is the reverse function so that labels tth and d_spacing are not correct.
            # print(tth_in)
            dspc_out = 2 * np.degrees(np.arcsin(wavelength / 2 / tth_in))
            # dspc_out = 2*np.degrees(np.arcsin(wavelength/2/tth_in[:,1]))
        return dspc_out

    def GetMask(self, MSKfile, ImInts):
        """
        Dioptas mask is compressed Tiff image.
        Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py
        :param MSKfile:
        :param ImInts:
        :return:
        """
        ImMsk = np.array(Image.open(MSKfile))
        ImInts = ma.array(ImInts, mask=ImMsk)

        # TALK TO SIMON ABOUT WHAT TO DO WITH THIS

        '''
        if debug:
            # N.B. The plot is a pig with even 1000x1000 pixel images and takes a long time to render.
            if ImTTH.size > 100000:
                print(' Have patience. The mask plot will appear but it can take its time to render.')
            fig_1 = plt.figure()
            ax1 = fig_1.add_subplot(1, 3, 1)
            ax1.scatter(ImTTH, ImAzi, s=1, c=(ImInts.data), edgecolors='none', cmap=plt.cm.jet, vmin=0,
                        vmax=np.percentile(ImInts.flatten(), 98))
            ax1.set_title('All data')

            ax2 = fig_1.add_subplot(1, 3, 2)
            ax2.scatter(ImTTH, ImAzi, s=1, c=ImMsk, edgecolors='none', cmap='Greys')
            ax2.set_title('Mask')
            ax2.set_xlim(ax1.get_xlim())

            ax3 = fig_1.add_subplot(1, 3, 3)
            ax3.scatter(ImTTH, ImAzi, s=1, c=(ImInts), edgecolors='none', cmap=plt.cm.jet, vmin=0,
                        vmax=np.percentile(ImInts.flatten(), 98))
            ax3.set_title('Masked data')
            ax3.set_xlim(ax1.get_xlim())
            # ax2.colorbar()
            plt.show()

            plt.close()
        '''

        # FIX ME: need to validate size of images vs. detector name. Otherwise the mask can be the wrong size
        # det_size = pyFAI.detectors.ALL_DETECTORS['picam_v1'].MAX_SHAPE
        # FIX ME : this could probably all be done by using the detector class in fabio.
        return ImInts

    def get_calibration(self, filenam):
        """
        Parse inputs file, create type specific inputs
        :param filenam:
        :return:
        """
        parms_file = open(filenam, 'r')
        filelines = parms_file.readlines()
        parms_dict = {}
        for item in filelines:
            new_parms = item.strip('\n').split(':', 1)
            parm = new_parms[1].strip()
            value = None
            # print parm
            try:
                value = int(parm)
                parms_dict[(str(new_parms[0]).lower())] = value
            except ValueError:
                try:
                    value = float(parm)
                    parms_dict[new_parms[0]] = value
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
                        parms_dict[new_parms[0]] = newlist
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
                        parms_dict[new_parms[0]] = newdict
                    elif not parm:
                        parms_dict[new_parms[0]] = ''

                    else:
                        parms_dict[new_parms[0]] = str(parm)

        # process 'detector_config' into pixel sizes if needed.
        if 'Detector_config' in parms_dict:
            # print(parms_dict['Detector_config'])
            parms_dict['pixelsize1'] = parms_dict['Detector_config']['pixel1']
            parms_dict['pixelsize2'] = parms_dict['Detector_config']['pixel2']
        # get wavelengths in Angstrom
        parms_dict['conversion_constant'] = parms_dict['Wavelength'] * 1E10
        # force all dictionary labels to be lower case -- makes s
        parms_dict = {k.lower(): v for k, v in parms_dict.items()}
        # report type
        parms_dict['DispersionType'] = 'AngleDispersive'
        return parms_dict

    def bins(self, azimu, orders):
        """
        determine bins to use in initial fitting.
        assign each data to a chunk corresponding to its azimuth value
        Returns array with indices for each bin and array of bin centroids
        :param azimu:
        :param orders:
        :return:
        """
        bins = orders['AziBins']
        azmax = azimu.max()
        azmin = azimu.min()
        binsize = (azmax - azmin) / bins
        chunks = []
        azichunks = []
        tempazi = azimu.flatten()
        for i in range(bins):
            start = azmin + i * binsize
            end = azmin + (i + 1) * binsize
            azichunk = np.where((tempazi > start) & (tempazi <= end))
            azichunks.append(((end - start) / 2) + start)
            chunks.append(azichunk)
        return chunks, azichunks

    def GetTth(self, poni, det=None, imask=None):
        """
        Give 2-theta value for detector x,y position; calibration info in data
        :param cali_file:
        :param poni:
        :param pix:
        :param det:
        :param imask:
        :return:
        """

        ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'],
                                 poni['rot3'], pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector=det,
                                 wavelength=poni['wavelength'])
        if imask is not None:
            return ma.array(np.degrees(ai.twoThetaArray()), mask=imask)
        else:
            return ma.array(np.degrees(ai.twoThetaArray()))

    def GetAzm(self, poni, det=None, imask=None):
        """
        Give azimuth value for detector x,y position; calibration info in data
        :param cali_file:
        :param poni:
        :param pix:
        :param det:
        :param imask:
        :return:
        """
        ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'],
                                 poni['rot3'], pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector=det,
                                 wavelength=poni['wavelength'])
        if imask is not None:
            return ma.array(np.degrees(ai.chiArray()), mask=imask)
        else:
            return ma.array(np.degrees(ai.chiArray()))

    def GetDsp(self, poni, det=None, imask=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        Get as thetheta and then convert into d-spacing
        :param cali_file:
        :param poni:
        :param pix:
        :param det:
        :param imask:
        :return:
        """
        ai = AzimuthalIntegrator(poni['distance'], poni['poni1'], poni['poni2'], poni['rot1'], poni['rot2'],
                                 poni['rot3'], pixel1=poni['pixelsize1'], pixel2=poni['pixelsize2'], detector=det,
                                 wavelength=poni['wavelength'])
        if imask is not None:
            return ma.array(poni['wavelength'] * 1E10 / 2 / np.sin(ai.twoThetaArray() / 2), mask=imask)
        else:
            return ma.array(poni['wavelength'] * 1E10 / 2 / np.sin(ai.twoThetaArray() / 2))




    def plot(ImDispersion, ImAzimuths, ImIntensity, dtype='data', masked=True, ImIntensity2=None, name=''):
        # Plot data.
        # possibilities:
        #   1. just plotting the data or model - with or without mask: label 'data'
        #   2. plot of data, model and differences: label 'model'
        #   3. plot of all data, mask and masked data. label 'mask'
        to_plot = []
        if dtype == 'data':
            x_plots = 1
            y_plots = 1
            to_plot.append(ImIntensity)
            plot_mask = [masked]
            plot_title = ['Data']
            plot_cmap = [plt.cm.jet]
            
            spec = gridspec.GridSpec(ncols=x_plots, nrows=y_plots,
                             width_ratios=[1], wspace=0.5,
                             hspace=0.5, height_ratios=[1])
        elif dtype == 'model':
            x_plots = 3
            y_plots = 1
            to_plot.append(ImIntensity)
            to_plot.append(ImIntensity2)
            to_plot.append(ImIntensity - ImIntensity2)
            plot_mask = [masked, masked, masked]
            plot_title = ['Data', 'Model', 'Residuals']
            plot_cmap = ['jet', 'jet', 'jet']
            
        elif dtype == 'mask':
            x_plots = 3
            y_plots = 2
            to_plot.append(ImIntensity)
            to_plot.append(np.array(ma.getmaskarray(ImIntensity), dtype='uint8')+1)
            to_plot.append(ImIntensity)
            plot_mask = [False, False, True]
            plot_title = ['All Data', 'Mask', 'Masked Data']
            plot_cmap = 'jet', 'Greys', 'jet'
            
            spec = gridspec.GridSpec(ncols=x_plots, nrows=y_plots,
                             width_ratios=[1,1,1], wspace=0.5,
                             hspace=0.5, height_ratios=[2,1])
                    
        else:
            stop    
        
        y_lims = np.array([np.min(ImAzimuths.flatten()), np.max(ImAzimuths.flatten())])
        y_lims = np.around(y_lims / 180) * 180


        # N.B. The plot is a pig with even 1000x1000 pixel images and takes a long time to render.
        if ImIntensity.size > 100000:
            print(' Have patience. The plot(s) will appear but it can take its time to render.')
            
        
        fig_1 = plt.figure()
        for i in range(x_plots):
            ax1 = fig_1.add_subplot(spec[i])
            if plot_mask[i] == True:
                ax1.scatter(ImDispersion, ImAzimuths, s=1, c=(to_plot[i]), edgecolors='none', cmap=plot_cmap[i], vmin=0) #plot all data includeing masked data.
            else:
                ax1.scatter(ImDispersion.data, ImAzimuths.data, s=1, c=(to_plot[i].data), edgecolors='none', cmap=plot_cmap[i], vmin=0)
            ax1.set_title(plot_title[i])
            ax1.set_ylim(y_lims)
            locs, labels = plt.xticks()
            plt.setp(labels, rotation=90)
            if i == 0:
                ax1.set_ylabel(r'Azimuth ($^\circ$)')
            ax1.set_xlabel(r'$2\theta$ ($^\circ$)')
            
        if y_plots > 1:
            for i in range(x_plots):
                ax1 = fig_1.add_subplot(spec[i+x_plots])
                if dtype == 'mask' and i==1:
                    #plot cdf of the intensities.
                    # sort the data in ascending order
                    x1 = np.sort(ImIntensity.data)
                    x2 = np.sort(ImIntensity)
                      
                    # get the cdf values of y
                    y1 = np.arange(np.size(x1)) / float(np.size(x1))
                    y2 = np.arange(np.size(x2)) / float(ma.count(x2))
                      
                    #ax1 = fig_1.add_subplot(1, 1, 1)   
                    ax1.plot(x1, y1, marker='.', ms=3, label='unmasked CDF')
                    ax1.plot(x2, y2, marker='.', ms=3, label='cleaned CDF')
                    #ax1.legend()
                    ax1.set_title('CDF of the intensities')
                    
                else:
                    if plot_mask[i] == True:
                        ax1.scatter(ImDispersion, to_plot[i], s=1, c=(ImAzimuths), edgecolors='none', cmap=plot_cmap[i], vmin=0) #plot all data includeing masked data.
                    else:
                        ax1.scatter(ImDispersion.data, to_plot[i].data, s=1, c=(ImAzimuths.data), edgecolors='none', cmap=plot_cmap[i], vmin=0)
                    #ax1.set_title(plot_title[i])
                    #ax1.set_ylim(y_lims)
                    locs, labels = plt.xticks()
                    plt.setp(labels, rotation=90)
                    if i == 0:
                        ax1.set_ylabel(r'Intennsity (a.u.)')
                    ax1.set_xlabel(r'$2\theta$ ($^\circ$)')
                    
        plt.suptitle(name + '; masking')
                    
        plt.show()
        plt.close()
        
        '''
        if dtype == 'mask':
        
            #plot cdf of the intensities.
            # sort the data in ascending order
            x1 = np.sort(ImIntensity)
              
            # get the cdf values of y
            y1 = np.arange(np.size(x1)) / float(np.size(x1))
              
            fig_2 = plt.figure()   
            ax1 = fig_2.add_subplot(1, 1, 1)   
            ax1.plot(x1, y1, marker='o')
            ax1.set_title('CDF of the intensities')
            
            plt.show()
        '''
        return fig_1
