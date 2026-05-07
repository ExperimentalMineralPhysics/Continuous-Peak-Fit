#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["DioptasDetector"]

import sys
from copy import copy, deepcopy
from importlib.metadata import version
import os
import re
from pathlib import Path
import types
import fabio
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyFAI
from pyFAI.detectors.orientation import Orientation
from packaging.version import Version

# Logic to support multiple PyFAI versions
if Version(version("pyFAI")).major >= 2025:
    from pyFAI.integrator.azimuthal import AzimuthalIntegrator
else:
    from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from pyFAI.io import ponifile

import cpf.h5_functions as h5_functions
import cpf # need to import whole package to avoind trying to import part of incompletely iniated method (cpf.settings.issettings for _get_metadata)
from cpf import IO_functions
from cpf.input_types._AngleDispersive_common import _AngleDispersive_common
from cpf.input_types._metadata_common import _metadata_common
from cpf.input_types._Masks import _masks
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive
from cpf.util.logging import get_logger

from pyFAI.io import DefaultAiWriter

logger = get_logger("cpf.input_types.DioptasFunctions")


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

        self.plot_orientation = "vertical"

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

        self.Dispersion = "Angle"
        self.continuous_azm = True

        self.Dispersionlabel = r"2$\theta$"
        self.DispersionUnits = r"$^\circ$"
        self.Azimuthlabel = r"Azimuth"
        self.AzimuthUnits = r"$^\circ$"
        self.Observationslabel = r"Intensity"
        self.ObservationsUnits = r"counts"

        self.azm_blocks = 45

        self.reduce_by = None

        self.metadata = None
        self._default_metadata_labels_tiff  = {"time": "FILE_MODIFIED", # file creation time.
                                  }
        self._default_metadata_labels_hdf5 = {}
        # self._default_metadata_labels = {"time": "FILE_MODIFIED", # file creation time.
        #                           }
        
        self._default_h5_datakey = '/*.1/measurement/p3/'
        self._default_h5_iterate = [{"from": 0, "to": 0, "step": 1, "label":["pos"], "do":"iterate"},
                           {"do":"combine", 
                 "from": 0, 
                 "to": -1, 
                 "step": 1, 
                 "using":"position",
                 "label": ['pos'],
                 # "pos": '/*.1/measurement/azim/',
                 "dim": 0}]
        
        self.calibration = None
        self.conversion_constant = None
        self.detector = None

        if settings_class:
            if settings_class.calibration_parameters != None:
                self.get_calibration(settings=settings_class)
            if self.calibration:
                self.detector = self.get_detector(settings=settings_class)

    def duplicate(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf], with_detector=True, as_masked=True):
        """
        Makes an independent copy of a DioptasDetector Instance.

        range_bounds and azi_bounds restrict the extent of the data if needed.
        The range resturictions should be applied upon copying (if required) for memory efficieny.
        Alternatively run:
        data_class.duplicate()
        date_calss.set_limits(range_bounds=[...], azi_bounds=[...])

        Parameters
        ----------
        range_bounds : dict or array, optional
            Limits for the two theta range. The default is [-np.inf, np.inf].
        azi_bounds : dict or array, optional
            Limits for the azimuth range. The default is [-np.inf, np.inf].
        with_detector : bool, optional
            Include detector and calibrations in returned class or not. The default is True.
        as_masked : bool, optional
            Return arrays as numpy masked arrays or numpy arrays. The default is False.

        Returns
        -------
        new : DioptasDetector Instance.
            Copy of DioptasDetector with independedent data values

        """
        #FIXME: should merge the data reduction funtions here with set_limits. or call set_limits.
        #FIXME: should be able to copy the class and reduce data at the same time. rather than copying and then reducing
        # this is not memory efficient. 
        
        
        # list variables that are not just straight copied
        copy_separately = ["intensity", "tth", "azm",
                           "x", "y", "z",
                           "dspace",
                           "tth_start", "tth_end"]
        dont_copy = ["detector", "calibration"]
        
        #validate the ranges
        range_bounds, azi_bounds = self.check_bounds(range_bounds, azi_bounds)
        
        if with_detector:
            new = copy(self)
            # do not return here because likely need to cut the data down
        elif 0:
            # copy and then delete the detector and calibration, 
            # so that all other non-default values are propagated.    
            new = deepcopy(self)
            for i in dont_copy:
                setattr(new, i, None)
            # new.detector = None
            # new.calibration = None
        else:
            # make new detector instance.
            # assume the methods are not altered and 
            # add the common variables to ensure consistent behaviour
            new = DioptasDetector()
            
            # copy all the settings ignoring any methods
            for i in dir(new):
                if (type(getattr(new, i)) == types.MethodType or 
                    i[:2] == "__"):
                    # skip method copying or default object
                    continue
                elif i in dont_copy:
                    # set dont copy parameters to 0
                    setattr(new, i, None)
                elif i not in copy_separately:
                    # copy common parameters 
                    setattr(new, i, getattr(self, i))
                elif i in copy_separately:
                    # skip adding the variables in 'copy_separately'
                    # these are added below.
                    # *may* be more memory efficient than copying huge arrays and 
                    # then making them smaller.
                    continue
                else:
                    raise ValueError("Should not be possible to get here")
            
        # set new range.
        new.tth_start = range_bounds[0]
        new.tth_end = range_bounds[1]
        
        # restrict the data. 
        local_mask = np.where(
            (self.tth >= new.tth_start)
            & (self.tth <= new.tth_end)
            & (self.azm >= azi_bounds[0])
            & (self.azm <= azi_bounds[1])
        )
        new.intensity = self.intensity[local_mask]
        new.tth = self.tth[local_mask]
        new.azm = self.azm[local_mask]
        if "dspace" in dir(self):
            new.dspace = self.dspace[local_mask]
       
        if "x" in dir(new):
            if self.x is not None:
                new.x = self.x[local_mask]
        if "y" in dir(new):
            if self.y is not None:
                new.y = self.y[local_mask]
        if "z" in dir(new):
            if self.z is not None:
                new.z = self.z[local_mask]

        if as_masked == False and ma.isMaskedArray(new.intensity):
            # return flat arrays.
            new.intensity = new.intensity.compressed()
            new.tth = new.tth.compressed()
            new.azm = new.azm.compressed()
            if "dspace" in dir(new):
                new.dspace = new.dspace.compressed()
            if "x" in dir(new) and new.x is not None:
                new.x = new.x.compressed()
            if "y" in dir(new) and new.y is not None:
                new.y = new.y.compressed()
            if "z" in dir(new) and new.z is not None:
                new.z = new.z.compressed()

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
            
        if 0:
            """
            FIXME: If there is no specific detector type set in the calibration
            then is raises a warning:
              "No sensor configuration provided; using default behaviour."
            This is just unfortunate but does not cause a problem.
            the way round this is to read the poni file as a dictionary and then
            pipe it back to the ponifile object. 
            
            if not warnings are raised everytime the object is copied/duplicated
            
            The code below also raises warnings on other detectors that are specified. 
            
            Instead change how the duplication of the detector works.
            """
            import json
            
            # read a poni file into 
            pf_data = {}
            with open(parms_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    key, value = line.split(":", 1)
                    key = key.strip().lower()
                    value = value.strip()
                    # Try to parse JSON or numbers
                    try:
                        parsed_value = json.loads(value)
                    except json.JSONDecodeError:
                        parsed_value = value
                    pf_data[key] = parsed_value
            if pf_data["detector"].lower() == "detector":
                # it is a generic detector without sensor specification.
                # add this to supress warnings. 
                pf_data["detector_config"].update({'sensor': {'material': 'CdTe', 'thickness': 0.001}})
                pf_data["detector_config"].update({'sensor': pyFAI.detectors.sensors.SensorConfig.parse("CdTe, 1mm")})
            
            pf = ponifile.PoniFile()
            # pf.read_from_file(parms_file)
            pf.read_from_dict(pf_data)
        else:
            pf = ponifile.PoniFile()
            pf.read_from_file(parms_file)
            
        if pf.API_VERSION <2:
            # then there is no orientation information in the poni file
            error_str = "Support for poni v1 files has been depreciated. To proceed update your poni file to version>=2"
            logger.error(error_str)
            import sys
            sys.exit(error_str)
            
            # then there is no orientation information in the poni file
            # assume an orientation, add and update pf
            config = pf.detector.get_config()
            config["orientation"] = 2
            pf.detector.set_config(config)
            #set the orientation which means AP_VERIOSN >= 2
            pf.API_VERSION = 2.1
        if (
            not isinstance(settings.image_list[0], list) and 
            pf.as_dict().get("poni_version", 1) >= 2 and 
            "orientation" in pf.detector.get_config()
        ):
            # FIXME: this doesnt make sense to me.  
            # can't flip if using hdf5 type images (if image_list[0] is list) because 
            # the images are then upsude down relative to dioptas. 
            
            # Check orientation and patch it since pyFAI and Dioptas use different conventions:
            # - Dioptas convention: origin at the top right when looking from the sample
            # - Default pyFAI convention: origin at the bottom right when looking from the sample   
            
            """Flips the detector up-down orientation in a poni configuration dictionary. Changes the dictionary object
            in place.
            """
            """
            These functions replicate the functionality of Dioptas.
            copied from: Dioptas/dioptas/model/CalibrationModel.py
            https://github.com/Dioptas/Dioptas/blob/develop/dioptas/model/CalibrationModel.py
            """
            config = pf.detector.get_config()
            orientation = config["orientation"]
            if orientation in (Orientation.Unspecified, Orientation.BottomRight):
                config["orientation"] = Orientation.TopRight
            elif orientation == Orientation.TopRight:
                config["orientation"] = Orientation.BottomRight
            elif orientation == Orientation.BottomLeft:
                config["orientation"] = Orientation.TopLeft
            elif orientation == Orientation.TopLeft:
                config["orientation"] = Orientation.BottomLeft
            else:
                logger.error(
                    "Detector orientation is not supported: Saved .poni file is not compatible with pyFAI"
                )
            pf.detector.set_config(config)
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
                settings=settings, file_name=calibration_file,
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
                orientation=self.calibration.detector.get_config()["orientation"]
            )

            
    # @staticmethod
    def import_image(
        self, image_name=None, settings=None, mask=None, dtype=None, 
        reduce_by = None, debug=False
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
        if image_name == None and settings.subfit_filename == None:
            raise ValueError("Settings are given but no subpattern is set.")

        if self.detector == None:
            self.get_detector(settings)

        if image_name == None:
            # load the data for the chosen subpattern.
            image_name = settings.subfit_filename

        # read image
        if isinstance(image_name, list):
            # then it is a h5 type file (including *.nxs)
            im = h5_functions.get_images(image_name, settings_class=settings)
            
            #add metadata to instance.
            #done here so only need to open file once.
            self._set_metadata(None, settings=settings)

        else:
            # presume it is an image file.
            try:
                im_all = fabio.open(image_name)
                logger.moreinfo(f"This file contains {im_all.nframes} frame(s) with a combined shape of {im_all.shape}")
            except:
                err_str = "".join(
                    [
                        "The image type is not recognsed by Fabio.\n",
                        "It might be that the image type is a 'hdf5' or 'nxs' file type.\n",
                        "In which case call the images using the hdf5 functions.",
                    ]
                )
                raise ValueError(err_str)
                
            if im_all.nframes == 1:
                # squeeze to make sure 1st dimension is not 1.
                im = im_all.data.squeeze()  
            else:
                err_str = "".join(
                    [
                        "There is more than 1 image in the file.\n",
                        "Multiimage Tiffs are not currently implemented in continuous peak fit.",
                    ]
                )
                raise NotImplementedError(err_str)
                
            #add metadata to instance.
            #done here so only need to open file once.
            self._set_metadata(im_all, settings=settings)
        
        # Convert the input data from integer to float because the lmfit model values
        # inherits integer properties from the data.
        #
        # Allow option for the data type to be set.
        if dtype == None:
            if self.intensity is None and np.size(self.intensity) > 2:
                # self.intensity has been set before. Inherit the dtype.
                dtype = self.intensity.dtype
            else:
                dtype = self.GetDataType(im[0], minimumPrecision=False)
        im = ma.array(im, dtype=dtype)

        # Dioptas flips the images to match the orientations in Fit2D
        # Therefore implemented here to be consistent with Dioptas.
        # flip both the image and the calibration separately to allow subsequent 
        # loading of more images.
        im = np.array(im)[::-1]
        
        # reduce the size of the data (if called for)
        if reduce_by is not None or self.reduce_by is not None:
            if reduce_by is False:
                # used to allow the full data image to be read by data_fill as part of reading the calibrations
                pass
            elif reduce_by is not None and reduce_by != 1:
                im = self._reduce_array(im, reduce_by=reduce_by)
            elif self.reduce_by is not None and self.reduce_by != 1:
                im = self._reduce_array(im, reduce_by=self.reduce_by)
            else:
                pass

        if logger.is_below_level(level="DEBUG"):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            asdf = ax.imshow(np.log10(im+1E4))
            # print(np.log10(im))
            plt.title(IO_functions.title_file_names(image_name=image_name))
            plt.colorbar(asdf)
            plt.show()
            plt.close()

        
        # apply mask to the intensity array
        if mask == None and ma.is_masked(self.intensity) == False:
            self.intensity = ma.array(im)
            return ma.array(im)
        elif ma.is_masked(self.intensity) == True and self.intensity.mask.shape == im.shape:
            # apply mask from previous intensities and all are same size
            self.intensity = ma.array(im, mask=self.intensity.mask)
            return ma.array(im)
        else:#if mask is not None:
            # apply given mask
            self.intensity = ma.array(im, mask=self.get_mask(mask, im))
            return ma.array(im, mask=mask)
        # else:
        #     # apply new mask
        #     self.intensity = ma.array(im, mask=self.fill_mask(mask, im))
        #     return ma.array(im, mask=mask)


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

        if settings.reduce_by is not None:
            self.reduce_by = settings.reduce_by
            
        if settings.metadata_labels is not None:
            self.metadata_labels = settings.metadata_labels
        if (isinstance(diff_file, list) or 
              os.path.splitext(os.path.basename(diff_file))[1] == ".h5" or 
              os.path.splitext(os.path.basename(diff_file))[1] == ".nxs"):
            self._default_metadata_labels = self._default_metadata_labels_hdf5
            self.metadata_labels = self._default_metadata_labels_hdf5
        else:
            self._default_metadata_labels = self._default_metadata_labels_tiff
            self.metadata_labels = self._default_metadata_labels_tiff
            
        if self.detector == None:
            self.get_detector(settings=settings)

        # get the intensities (without mask) and without reduction.
        # No reduction because if reduced then the calibration is nonsence.
        self.intensity = self.import_image(diff_file, settings=settings, reduce_by=False)
        
        # FIXME: (June 2024) because of how self.detector is instanciated the
        # shape might not be correct (or recognised). Hence the check here and
        # inclusion of the shape in the array getting.

        if tuple(self.intensity.shape) != tuple(self.detector.detector.max_shape):
            # cast both shapes to tuples to prevent list != tuple error.
            raise ValueError(
                "The pixel size of the data and the detector are not the same"
            )

        self.tth = ma.array(
            self.detector.center_array(unit='2th_deg')
        )
        self.azm = ma.array(
            self.detector.center_array(unit='chi_deg')
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

        if self.reduce_by is not None:
            self.intensity = self._reduce_array(self.intensity)
            self.tth = self._reduce_array(self.tth)
            self.azm = self._reduce_array(self.azm, polar=True)
            
            if "original_mask" in dir(self):
                self.original_mask= self._reduce_array(self.original_mask)
                
            if make_zyx:
                self.z = self._reduce_array(self.z)
                self.y = self._reduce_array(self.y)
                self.x = self._reduce_array(self.x)

        self.azm_start = (
            np.floor(np.min(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )
        self.azm_end = (
            np.ceil(np.max(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )
        self.tth_start = np.min(self.tth.flatten())
        self.tth_end = np.max(self.tth.flatten())


    def _set_metadata(self, image_obj, settings=None):
        """
        Adds all possible metadata values as dictionary within the in the data_class.
        
        The values that are in settings.metadata (a list) are extracted subsequently using data_class.get_metadata 
               
        If the cpf settings class is provided and has the method 'metadata_read_func'
        then this method is used to override the internal default methods and is 
        used to create 'metadata_dictionary' which is parsed.
        In this case the settings class attribute 'metadata_labels' is still needed 
        to get required parts of the metadata. 

        For Dioptas functions the default is a fabio.open(file).header dictionary 
        
        Parameters
        ----------
        image_obj : fabio object
            image object to be parsed.
        settings : cpf settings class, optional
            If the settings class has method 'metadata_read' this overrides the 
            internal methods and is used to get the metadata. 
            The default is None.

        Returns
        -------
        metadata_dictionary
            dictionary of image metadata. 
        """
        # Defined as function to allow get_metadata to call universal image method
        metadata_dictionary = {}
        if not image_obj and not settings:
            # then nothing is provided
            # expected behaviour in some circumstances.
            self.metadata = None
            return
        elif isinstance(image_obj, list) or (not image_obj and settings) or cpf.settings.is_settings(image_obj):
             # when calling hdf5 file there is no image_obj to send (= None) and settings is
             # provided instead. 
             # 
             # (any(False if x is None else "*" in x for x in self.metadata_labels.values()) or 
             #  any(False if x is None else "/" in x for x in self.metadata_labels.values()) or 
             #  (settings and any("*" in x for x in settings.metadata)) or 
             #  (settings and any("/" in x for x in settings.metadata))
             #  ):
            # here we set pointers to the things needed when the metadata is read.
            # assuming that it is not wise (or possible) to list all the possible hdf5
            # keys which could be read as metadata. 
            
            metadata_dictionary = {}
            if settings.subfit_filename is not None:
                metadata_dictionary["image"] = settings.subfit_filename
            else:
                metadata_dictionary["image"] = settings.image_list[0][0]
            metadata_dictionary["note"] = "It is not reasonable to load all hdf5 keys into a dictionary as metadata. Instead carry file name and use keys"
            metadata_dictionary["h5_datakey"] = settings.h5_datakey       
            # add the file creation and modifications time
            # metadata_dictionary.update(self._get_file_created_modified(image_obj[0]))
        else:
            if isinstance(image_obj, str) or isinstance(image_obj, Path):
                metadata_dictionary.update(fabio.open(image_obj).header)
            else:
                metadata_dictionary.update(image_obj.header)       
            # add the file creation and modifications time
            metadata_dictionary.update(self._get_file_created_modified(image_obj))

        if settings and "metadata_read_func" in settings.__dict__:
            metadata_dictionary.update(settings.metadata_read_func(settings, image_obj=image_obj))     
        self.metadata = metadata_dictionary

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
                    logger.info(f"The settings file requires a parameter called '{par}'")
                    all_present = 0
            if all_present == 0:
                sys.exit(
                    "The highlighted settings are missing from the input file. Fitting cannot proceed until they "
                    "are all present."
                )
        return required_list
    
    def detector_description(self):
        """
        Returns a text description of the detector. 
        
        For Dioptas detectors this is called from PyFAI.  

        Returns
        -------
        description : string
            Text description of the detector and the calibration.
        """
        description = DefaultAiWriter(None, self.detector).make_headers()
        description = description.replace('\r\n', '\n')
        return description

    # add common function.
    _get_d_space = _AngleDispersive_common._get_d_space
    conversion = _AngleDispersive_common.conversion
    bins = _AngleDispersive_common.bins
    set_limits = _AngleDispersive_common.set_limits
    test_azims = _AngleDispersive_common.test_azims
    GetDataType = _AngleDispersive_common.GetDataType
    duplicate_without_detector = _AngleDispersive_common.duplicate_without_detector
    check_bounds = _AngleDispersive_common.check_bounds
    _reduce_array = _AngleDispersive_common._reduce_array
    get_metadata = _metadata_common.get_metadata
    _get_file_created_modified = _metadata_common._get_file_created_modified


    # add masking functions to detetor class.
    get_mask = _masks.get_mask
    set_mask = _masks.set_mask
    mask_apply = _masks.mask_apply
    mask_restore = _masks.mask_restore
    mask_remove = _masks.mask_remove

    # these methods are all called from _Plot_AngleDispersive as they are shared with other detector types.
    # Each of these methods remains here because they are called by higher-level functions:
    plot_masked = _Plot_AngleDispersive.plot_masked
    plot_fitted = _Plot_AngleDispersive.plot_fitted
    plot_collected = _Plot_AngleDispersive.plot_collected
    plot_calibrated = _Plot_AngleDispersive.plot_calibrated
    plot_integrated = _Plot_AngleDispersive.plot_integrated
    what_plot_type = _Plot_AngleDispersive.what_plot_type

    # this function is added because it requires access to self:
    dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks
