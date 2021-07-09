#!/usr/bin/env python
# -*- coding: utf-8 -*-

# __all__ = ['Requirements', 'ImportImage', 'DetectorCheck', 'Conversion', 'GetMask', 'GetCalibration', 'GetTth',
#            'GetAzm', 'GetDsp']

import sys
import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from cpf import Med


# For MED files from Energy dispersive diffraction patterns at X17B2, 6-BMB and others.

# For ease of use the *.med files are imported using Mark Rivers' Med.py and Mca.py functions.
# These were downloaded from ??? in June 2019.
# and edited to remove unnecessary dependencies e.g. 'Numeric' and 'LinearAlgebra' packages and functions that had
# further dependencies.
# I don't know how much of this is required but for now I am just adding the packages.

# FIX ME: Need a dependency check.

class MedDetector:

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

    def fill_data(self, diff_file, settings=None, debug=None, calibration_mask=None):
        """
        :param calibration_data:
        :param detector:
        :param debug:
        :param calibration_mask:
        """
        # self.detector = self.detector_check(calibration_data, detector)
        self.calibration = self.get_calibration(settings.Calib_param)
        self.get_detector(diff_file, settings)
        
        self.intensity = self.get_masked_calibration(diff_file, debug, calibration_mask)
        self.tth = self.GetTth(mask=self.intensity.mask)
        self.azm = self.GetAzm(mask=self.intensity.mask)
        self.dspace = self.GetDsp(mask=self.intensity.mask)
        

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
                        # 'AziBins',            # required based on detector type
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
    def ImportImage(ImageName, mask=None, debug=False):
        """
        :param ImageName:
        :param debug:
        :return:
        """
        # FIX ME: why is mask in here? should it ingerit from previous data if it exists?
        filename, file_extension = os.path.splitext(ImageName)
        if file_extension == '.med':
            dat = Med.Med()
            dat.read_file(ImageName)
            im_all = []
            for x in range(dat.n_detectors):
                im_all.append(dat.mcas[x].get_data())
            if mask is not None:
                im_all = ma.array(im_all, mask=mask)
            else:
                im_all = ma.array(im_all)
        else:
            sys.exit('Unknown file type. A *.med file is required')
        return im_all

    @staticmethod
    def detector_check(ImageName, detector=None):
        """
        :param ImageName:
        :param detector:
        :return:
        """
        # FIX ME: we should try importing the image and reading the detector type from it.
        if detector == None:
            sys.exit('\n\nDioptas requires a detector type.\nTo list the possible detector types in the command line type:'
                     '\n   import pyFAI\n   pyFAI.detectors.Detector.registry\n\n')
        # FIX ME: check the detector type is valid.
        return detector


    def conversion(self, E_in, calib, reverse=0, azm=None):
        """
        :param E_in:
        :param calib:
        :param reverse:
        :param azm:
        :return:
        """
        # convert energy into d-spacing.
        if azm is not None:
            E_in = np.array([[azm, E_in]])
            # FIX ME: this is a horrible way to make the conversion work but it is needed for
            # processing the guesses of MED detectors.
        else:
            E_in = np.array([[0, E_in]])

        # print('E_in',E_in, 'calib,', calib)
        # print('type E_in', type(E_in))
        # print(E_in[:,1])
        # Convert Energy (keV) to wavelength (Angstroms)
        # E = hc/lamda;  h = 4.135667662(25)×10−15 eV s; c = 299792458 m/s
        wavelength = 4.135667662e-15 * (299792458 * 1E10) / (E_in[:, 1] * 1000)

        # determine which detector calibration goes with each energy.
        sort_idx = np.array(calib['azimuths']).argsort()
        de = sort_idx[np.searchsorted(calib['azimuths'], E_in[:, 0], sorter=sort_idx)]

        if not reverse:
            # convert energy to d-spacing
            dspc_out = []
            for i in range(len(de)):
                dspc_out.append(wavelength[i] / 2 / np.sin(np.radians(calib['calibs'].mcas[i].calibration.two_theta / 2)))
            dspc_out = np.squeeze(np.array(dspc_out))
        else:
            # convert dspacing to energy.
            # N.B. this is the reverse finction so that labels tth and dspacing are not correct.
            dspc_out = []
            for i in range(len(de)):
                dspc_out.append(wavelength[i] / 2 / np.sin(np.radians(calib['calibs'].mcas[i].calibration.two_theta / 2)))
                # dspc_out.append(2*np.degrees(np.arcsin(wavelength/2/calib['calibs'].mcas[i].calibration.two_theta)))
                print('dspc_out', dspc_out)
            dspc_out = np.array(dspc_out)

        return dspc_out


    def GetMask(self, MSK, ImInts=None, ImTTH=None, ImAzi=None, Imy=None, Imx=None, debug=False):
        """
        :param MSKfile:
        :param ImInts:
        :param ImTTH:
        :param ImAzi:
        :param Imy:
        :param Imx:
        :param debug:
        :return:
        """
        # Masks for 10 element detector are either for full detector or for an energy range.
        # Currently this just works by removing complete detectors.
        ImMsk = ma.zeros(ImInts.shape)
        for x in range(len(MSK)):
            ImMsk[MSK[x] - 1] = 1

        if debug:
            # Plot mask.
            # This is left in here for debugging.
            fig = plt.figure()
            ax = fig.add_subplot(1, 2, 1)
            plt.subplot(121)
            plt.scatter(ImTTH, ImAzi, s=4, c=ImInts, edgecolors='none', cmap=plt.cm.jet)
            ax.set_ylim([0, 360])
            ax = fig.add_subplot(1, 2, 2)
            plt.subplot(122)
            plt.scatter(ma.array(ImTTH, mask=ImMsk), ma.array(ImAzi, mask=ImMsk), s=4, c=ImInts, edgecolors='none',
                        cmap=plt.cm.jet)
            ax.set_ylim([0, 360])
            plt.colorbar()
            plt.show()
            plt.close()

        ImInts = ma.array(ImInts, mask=ImMsk)
        return ImInts

    def get_calibration(self, filenam):
        """
        Parse inputs file, create type specific inputs.
        :param filenam:
        :return:
        """
        parms_dict = {}
        dat = Med.Med(file=filenam)
        parms_dict['DispersionType'] = 'EnergyDispersive'
        parms_dict['calibs'] = dat
        if len(dat.mcas) == 10:
            az = np.linspace(0, 180, 9)
            az = np.append(az, 270)
            az[0] = az[0] + 2
            az[8] = az[8] - 2
            parms_dict['azimuths'] = az
        else:
            sys.exit('\n\nDetector is unknown. Edit function to add number of detectors and associated azimuths.\n\n')

        # 2 theta angle in degrees
        # parms_dict['conversion_constant'] = 6.5
        all_angle = []
        for x in range(parms_dict['calibs'].n_detectors):
            all_angle.append(parms_dict['calibs'].mcas[x].calibration.two_theta)
        parms_dict['conversion_constant'] = all_angle
        # print(parms_dict['conversion_constant'])
        return parms_dict

    def bins(self, azimu, orders):
        """
        Determine bins to use in initial fitting.
        Assign each data to a chunk corresponding to its azimuth value
        Returns array with indices for each bin and array of bin centroids
        :param azimu:
        :param orders:
        :return:
        """
        bin_vals = np.unique(azimu)
        chunks = []
        azichunks = []
        tempazi = azimu.flatten()
        for i in range(len(bin_vals)):
            azichunk = np.where((tempazi == bin_vals[i]))
            chunks.append(azichunk)
        return chunks, bin_vals

    def GetTth(self, mask=None, pix=None, det=None):
        """
        Called get two theta but it is really getting the energy.
        :param cali_file:
        :param parms_dict:
        :param pix:
        :param det:
        :return:
        """
        channels = np.arange(self.calibration['calibs'].mcas[0].data.shape[0])
        #channels = np.arange(parms_dict['calibs'].mcas[0].data.shape[0])
        en = []
        for x in range(self.calibration['calibs'].n_detectors):
            en.append(self.calibration['calibs'].mcas[x].channel_to_energy(channels))
        if mask is not None:
            en = ma.array(en, mask=mask)
        else:
            en = ma.array(en)
        return en

    def GetAzm(self, mask=None, pix=None, det=None):
        """
        Give azimuth value for detector x,y position; have to assume geometry of the detectors
        :param cali_file:
        :param parms_dict:
        :param pix:
        :param det:
        :return:
        """
        l = self.calibration['calibs'].mcas[0].data.shape
        azi = np.tile(self.calibration['azimuths'], [l[0], 1]).T
        if mask is not None:
            azi = ma.array(azi, mask=mask)
        else:
            azi = ma.array(azi, mask=None)
        return azi

    def GetDsp(self, mask=None, pix=None, det=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        :param cali_file:
        :param parms_dict:
        :param pix:
        :param det:
        :return:
        """
        channels = np.arange(self.calibration['calibs'].mcas[0].data.shape[0])
        dspc = []
        for x in range(self.calibration['calibs'].n_detectors):
            dspc.append(self.calibration['calibs'].mcas[x].channel_to_d(channels))
        if mask is not None:
            dpsc = ma.array(dspc, mask=mask)
        else:
            dpsc = ma.array(dspc)
        return dpsc
