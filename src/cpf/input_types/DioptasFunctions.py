#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["DioptasDetector"]


import re
import sys
from copy import deepcopy

import fabio
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyFAI
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from pyFAI.io import ponifile

import cpf.h5_functions as h5_functions
import cpf.logger_functions as lg
from cpf import IO_functions
from cpf.input_types._AngleDispersive_common import _AngleDispersive_common
from cpf.input_types._Masks import _masks

# from PIL import Image
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive

# from cpf.XRD_FitPattern import logger
from cpf.logger_functions import logger


class DioptasDetector:
    # For Dioptas functions to change
    #   -   Load/import data
    #   -   Load mask
    #   -   Load Calibration
    def __init__(self, settings_class=None):
        """
        :param calibration_parameters:
        """

        self.requirements = self.get_requirements()

        self.calibration = None
        self.detector = None

        self.intensity = None
        self.tth = None
        self.azm = None
        # self.dspace = None
        self.x = None
        self.y = None
        self.azm_start = -180
        self.azm_end = 180
        self.tth_start = None
        self.tth_end = None
        self.DispersionType = "AngleDispersive"
        self.continuous_azm = True

        self.azm_blocks = 45

        self.calibration = None
        self.conversion_constant = None
        self.detector = None

        if settings_class:
            self.get_calibration(settings=settings_class)
        if self.calibration:
            self.detector = self.get_detector(settings=settings_class)

    def duplicate(self):
        """
        Makes an independent copy of a Dioptas Instance

        Parameters
        ----------
        None.

        Returns
        -------
        Dioptas Instance.

        """
        new = deepcopy(self)
        return new

    def get_calibration(self, file_name=None, settings=None):
        """
        Opens the file with the calibration data in it and updates
        DioptasClass.calibration and
        DioptasClass.conversion_constant

        Either file_name or settings should be set.

        Parameters
        ----------
        file_name : string, optional
            Filename of the calibration file. In this case a *.poni file.
            The default is None.
        settings : settings class, optional
            cpf settings class containing the calibration file name.
            The default is None.

        Returns
        -------
        None.

        """
        if settings != None:
            parms_file = settings.calibration_parameters
        else:
            parms_file = file_name
        pf = ponifile.PoniFile()
        pf.read_from_file(parms_file)

        self.calibration = pf
        self.conversion_constant = pf.wavelength * 1e10  # in angstroms

    def get_detector(
        self, settings=None, calibration_file=None, diffraction_data=None, debug=False
    ):
        """
        Takes the detector information from the settings class, or the
        calibration *.json file and creates a detector-type instance which converts
        the x,y,azimith of the diffraction data pixels into two-theta vs.
        azimuth.

        Either settings or the file_name are required.

        Parameters
        ----------
        file_name : string, optional
            Filename of the calibration file. In this case a *.poni file.
            The default is None.
        settings : settings class, optional
            cpf settings class containing the calibration file name.
            The default is None.
        debug : True/False, optional
            Additional output, used for debigging.
            The default is False.

        Returns
        -------
        None.

        """
        if self.calibration == None:
            self.get_calibration(
                settings=settings, file_name=calibration_file, debug=debug
            )

        if self.calibration:
            # use the poni file/calibration to make the detector
            self.detector = AzimuthalIntegrator(
                detector=self.calibration.detector,
                dist=self.calibration.dist,
                poni1=self.calibration.poni1,
                poni2=self.calibration.poni2,
                rot1=self.calibration.rot1,
                rot2=self.calibration.rot2,
                rot3=self.calibration.rot3,
                wavelength=self.calibration.wavelength,
            )

            if (
                self.detector.detector.get_name() == "Detector"
                and "detector_config" in self.calibration.as_dict()
            ):
                config = self.calibration.as_dict()["detector_config"]

                if config["max_shape"] == None:
                    # open the file to get the shape of the data.
                    if diffraction_data is not None:
                        im_all = fabio.open(diffraction_data)
                    elif settings.calibration_data is not None:
                        im_all = fabio.open(settings.calibration_data)
                    elif settings.image_list[0] != None:
                        im_all = fabio.open(settings.image_list[0])
                    if im_all:
                        config["max_shape"] = im_all.shape

                # make sure the pixel sizes are correct
                self.detector.detector.set_config(config)

    # @staticmethod
    def import_image(
        self, image_name=None, settings=None, mask=None, dtype=None, debug=False
    ):
        """
        Import the data image into the intensity array.
        Apply new mask to the data (if given) otherwise use previous mask

        Parameters
        ----------
        image_name : string, optional
            Name of the image set to import. Either this or settings are required.
        settings : settings class, optional
            Cpf setting class that constins the image set to import.
            Either this or imagename is required.
        mask : array, optional
            Mask array to apply to data. The default is None.
        dtype : string, optional
            Data type string, to force the data type and bit depth. The default is None.
        debug : boolian, optional
            True/Flase to display debuging information. The default is False.

        Raises
        ------
        ValueError
            If there is no image_name of the settings.subpattern is not set.

        Returns
        -------
        Im : masked array
            Masked image intensity array.

        """
        # FIX ME: nxs files need to be checked. But they should be read like h5 files.

        # check inputs
        if image_name == None and settings.subpattern == None:
            raise ValueError("Settings are given but no subpattern is set.")

        if self.detector == None:
            self.get_detector(settings)

        if image_name == None:
            # load the data for the chosen subpattern.
            image_name = settings.subfit_filename

        # read image
        if isinstance(image_name, list):
            # then it is a h5 type file
            im = h5_functions.get_images(image_name)

        # elif file_extension == ".nxs":
        #     im_all = h5py.File(image_name, "r")
        #     # FIX ME: This assumes that there is only one image in the nxs file.
        #     # If there are more, then it will only read 1.
        #     im = np.array(im_all["/entry1/instrument/detector/data"])
        else:
            im_all = fabio.open(image_name)
            im = im_all.data

        # Convert the input data from integer to float because the lmfit model values
        # inherits integer properties from the data.
        #
        # Allow option for the data type to be set.
        if dtype == None:
            if self.intensity is None and np.size(self.intensity) > 2:
                # self.intensity has been set before. Inherit the dtype.
                dtype = self.intensity.dtype
            elif "int" in im[0].dtype.name:
                # the type is either int or uint - convert to float
                # using same bit precision
                precision = re.findall("\d+", im[0].dtype.name)[0]
                dtype = np.dtype("float" + precision)
        im = ma.array(im, dtype=dtype)

        # Dioptas flips the images to match the orientations in Fit2D
        # Therefore implemented here to be consistent with Dioptas.
        im = np.array(im)[::-1]

        if lg.make_logger_output(level="DEBUG"):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.imshow(im)
            plt.title(IO_functions.title_file_names(image_name=image_name))
            plt.show()
            plt.close()

        # apply mask to the intensity array
        if mask == None and ma.is_masked(self.intensity) == False:
            self.intensity = ma.array(im)
            return ma.array(im)
        elif mask is not None:
            # apply given mask
            self.intensity = ma.array(im, mask=self.fill_mask(mask, im))
            return ma.array(im, mask=mask)
        else:
            # apply mask from intensities
            self.intensity = ma.array(im, mask=self.intensity.mask)
            return ma.array(im)

    def fill_data(
        self, diff_file=None, settings=None, mask=None, make_zyx=False, debug=False
    ):
        """
        Initiates the data arrays.
        Creates the intensity, two theta and azimuth arrays from the
        Detector and the data.

        If diffraction data is provided this makes the intensity array otherwise
        the intensity array is filled with zeros.

        It must be called after setting the detector but before import_data.

        Parameters
        ----------
        diff_file : string, optional
            Filename of the diffraction data files to import.
        settings : settings class, optional
            cpf settings class containing the calibration file name.
            The default is None.
        mask : string or dirctionarry, optional
            Filename or dictionary of instructions to make data mask.
            The default is None.
        make_zyx : boolian, optional
            Switch to make x,y,z arrays for the pixels. These arrats are not
            used by default. They exist incase ever needed.
            The default is False.
        debug : True/False, optional
            Additional output, used for debugging.
            The default is False.

        Returns
        -------
        None.
        """

        # check inputs
        if diff_file != None:
            pass
        elif settings != None:  # and diff_file == None
            # load the data for the chosen subpattern.
            if settings.subfit_filename != None:
                diff_file = settings.subfit_filename
            else:
                raise ValueError("Settings are given but no subpattern is set.")
        else:
            raise ValueError("No diffraction file or settings have been given.")

        if mask == None:
            if settings.calibration_mask:  # in settings.items():
                mask = settings.calibration_mask

        if self.detector == None:
            self.get_detector(settings=settings)

        # get the intensities (without mask)
        self.intensity = self.import_image(diff_file)
        # FIXME: (June 2024) because of how self.detector is instanciated the
        # shape might not be correct (or recognised). Hence the check here and
        # inclusion of the shape in the array getting.

        if tuple(self.intensity.shape) != tuple(self.detector.detector.max_shape):
            # cast both shapes to tuples to prevent list != tuple error.
            raise ValueError(
                "The pixel size of the data and the detector are not the same"
            )

        self.tth = ma.array(
            np.rad2deg(self.detector.twoThetaArray(self.detector.detector.max_shape))
        )
        self.azm = ma.array(
            np.rad2deg(self.detector.chiArray(self.detector.detector.max_shape))
        )
        # self.dspace = self._get_d_space()
        if make_zyx:
            zyx = self.detector.calc_pos_zyx()
            self.z = zyx[0]
            self.y = zyx[1]
            self.x = zyx[2]

        # get and apply mask
        mask_array = self.get_mask(mask, self.intensity)
        self.mask_apply(mask_array, debug=debug)

        self.azm_start = (
            np.around(np.min(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )
        self.azm_end = (
            np.around(np.max(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )

    @staticmethod
    def detector_check(calibration_data, settings=None):
        """
        Get detector information
        :param settings:
        :param calibration_data:
        :return: detector:
        """
        # if "Calib_detector" in settings:
        # detector = pyFAI.detector_factory(settings.Calib_detector)
        if settings.calibration_detector is not None:
            detector = pyFAI.detector_factory(settings.calibration_detector)
        else:
            # if settings is None or detector == 'unknown' or detector == 'other' or detector == 'blank':
            im_all = fabio.open(calibration_data)
            # sz = calibration_data.Calib_pixels  # Pixel_size
            sz = calibration_data.calibration_pixel_size  # Pixel_size
            if sz > 1:
                sz = sz * 1e-6
            detector = pyFAI.detectors.Detector(
                pixel1=sz, pixel2=sz, splineFile=None, max_shape=im_all.shape
            )
        # FIX ME: check the detector type is valid.
        return detector

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
        #       - beginning and numbers for files (needs directory, file name, number digits, ending)
        required_list = [
            "Calib_type",
            "Calib_detector",
            # 'Calib_data',        # Removed because this is not strictly required.
            "Calib_param",  # Now added to calibration function file.
            # 'Calib_pixels',      # this is only needed for GSAS-II and should be read from image file. FIX
            # ME: add to GSAS-II input file.
            # 'Calib_mask',        # a mask file is not strictly required.
            "datafile_directory",
            "datafile_Basename",
            "datafile_Ending",
            # 'datafile_StartNum',  # now optionally replaced by datafile_Files
            # 'datafile_EndNum',    # now optionally replaced by datafile_Files
            "datafile_NumDigit",
            "AziBins",  # required based on detector type
            "fit_orders",
            # 'Output_type',		   # should be optional
            # 'Output_NumAziWrite',  # should be optional
            # 'Output_directory']	   # should be optional
        ]

        # Check required against inputs if given
        if parameter_settings is not None:
            # properties of the data files.
            all_present = 1
            for par in parameter_settings:
                if par in required_list:
                    logger.info(" ".join(map(str, [("Got: ", par)])))
                else:
                    logger.info(
                        " ".join(
                            map(
                                str,
                                [
                                    (
                                        "The settings file requires a parameter called  '",
                                        par,
                                        "'",
                                    )
                                ],
                            )
                        )
                    )
                    all_present = 0
            if all_present == 0:
                sys.exit(
                    "The highlighted settings are missing from the input file. Fitting cannot proceed until they "
                    "are all present."
                )
        return required_list


# add common function.
DioptasDetector._get_d_space = _AngleDispersive_common._get_d_space
DioptasDetector.conversion = _AngleDispersive_common.conversion
DioptasDetector.bins = _AngleDispersive_common.bins
DioptasDetector.set_limits = _AngleDispersive_common.set_limits
DioptasDetector.test_azims = _AngleDispersive_common.test_azims

# add masking functions to detetor class.
DioptasDetector.get_mask = _masks.get_mask
DioptasDetector.set_mask = _masks.set_mask
DioptasDetector.mask_apply = _masks.mask_apply
DioptasDetector.mask_restore = _masks.mask_restore
DioptasDetector.mask_remove = _masks.mask_remove

# these methods are all called from _Plot_AngleDispersive as they are shared with other detector types.
# Each of these methods remains here because they are called by higher-level functions:
DioptasDetector.plot_masked = _Plot_AngleDispersive.plot_masked
DioptasDetector.plot_fitted = _Plot_AngleDispersive.plot_fitted
DioptasDetector.plot_collected = _Plot_AngleDispersive.plot_collected
DioptasDetector.plot_calibrated = _Plot_AngleDispersive.plot_calibrated
# this function is added because it requires access to self:
DioptasDetector.dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks
