#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
XYFunctions.

For simple x, y data the does not require any complex conversions.
The peaks are asumed to run vertically in the image (Y direction)
If the peaks run the other way then use the YXFunctions (not yet made as a function).

The 'Calibration', such as it is, are simple finctions of x or y.

Functions to change for another imput type:
   -   get_calibration
   -   get_detector
   -   Load/import data
   -   Load mask
   -   Load Calibration



Notes (SAH, March 2023):
- The images are plotted as scatter plots. This is because I started (but didnt finish)
making the data plot as images using imshow. The switch "plot_as_image" switches between the scatter
and the imshow functions.
- Masks are either files/arrays (as in other data types) or dictionaries of polygons and thresholds, more like
the functionality of GSAS-II. But only polgons and threasholds are implemented not the other types avaliabe in
GSAS-II. An example of the new dictionary type mask is:
    Calib_mask     = {"polygon": [[(0,43),(180,43),(205,103),(127,160),(97,247),(80,393),(0,377)],
                                  [(1000,43),(820,43),(795,103),(853,160),(890,247),(910,393),(1000,377)],
                                  [(305,357),(286,375),(286,375),(305,275)],
                                 ],
                      "threshold": [0.25, 8]
                      }
- Calibrations are now also a dictionary that converts the x,y of the data set to calibrated values.
Currently only linear cnversions have been implemented. An example:
    Calib_param    = {"x_dim": 1,
                      "x": [0,1],
                      "x_label": r"pixels",
                      "y":[360, -360/1000],
                      "y_label": "Azimuth",
                      "y_start": 0,
                      "y_end": 360}
-

"""

__all__ = ["XYDetector"]


import os
import pickle
import re
import sys
from copy import copy, deepcopy

import fabio
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyFAI
from matplotlib import image

# import pyFAI
from PIL import Image

import cpf.h5_functions as h5_functions
from cpf import IO_functions
from cpf.input_types._AngleDispersive_common import _AngleDispersive_common
from cpf.input_types._Masks import _masks
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive
from cpf.util.logging import get_logger

logger = get_logger("cpf.input_types.XYFunctions")


# plot the data as an image (TRue) or a scatter plot (false).
# FIXME: the plot as image (im_show) does not work. The fitted data has to be reshaped
# into an image.
plot_as_image = True


def csv_to_image(image_name):
    # then it is a text file, assume there are less than 20 header rows.
    # assume the we do not know the delimiters.
    done = 0
    n = 0
    while done == 0:
        if os.path.splitext(image_name)[1] == ".csv":
            import csv

            im = []
            with open(image_name, "r", encoding="utf-8-sig") as f:
                reader = csv.reader(f)
                for row in reader:
                    for i in range(len(row)):
                        if row[i] == "":
                            row[i] = np.nan
                        else:
                            row[i] = float(row[i])
                        # row[row == ""] = 0
                    im.append(row)

            im = np.array(im)
            done = 1
        else:
            try:
                im = np.loadtxt(image_name, skiprows=n)
                done = 1
            except:
                done = 0
                n += 1
                if n > 20:
                    # issue an error
                    err_str = "There seems to be more than 20 header rows in the data file. Is this correct?"
                    raise ValueError(err_str)
    return im


class XYDetector:
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
        self.azm_start = np.inf
        self.azm_end = -np.inf
        self.tth_start = None
        self.tth_end = None
        
        # Image data so contonuous.
        self.Dispersion = "Angle"
        self.continuous_azm = True
        
        self.Dispersionlabel = r"Pixels"
        self.DispersionUnits = r"num"
        self.Azimuthlabel = r"Pixels"
        self.AzimuthUnits = r"num"
        self.Observationslabel = r"Intensity"
        self.ObservationsUnits = r"a.u."

        #set defualt value, do not asume is diffraction data
        self.azm_blocks = 100
            
        self.reduce_by = None
        
        self.calibration = None
        self.conversion_constant = None
        self.detector = None

        # set the start and end limits for the data
        # self.start = 0
        # self.end   = np.inf

        if settings_class:
            self.get_calibration(settings=settings_class)
        if self.calibration:
            self.detector = self.get_detector(settings=settings_class)

    def duplicate(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf], with_detector=True):
        """
        Makes an independent copy of a XYDetector Instance.

        range_bounds and azi_bounds restrict the extent of the data if needed.
        The range resturictions should be applied upon copying (if required) for memory efficieny.
        Laternatively run:
        data_class.duplicate()
        date_calss.set_limit2(range_bounds=[...], azi_bounds=[...])

        Parameters
        ----------
        range_bounds : dict or array, optional
            Limits for the two theta range. The default is [-np.inf, np.inf].
        azi_bounds : dict or array, optional
            Limits for the azimuth range. The default is [-np.inf, np.inf].

        Returns
        -------
        new : XYDetector Instance.
            Copy of XYDetector with independedent data values

        """

        if with_detector:
            new = copy(self)
        else:
            # copy and then delete the detector and calibration, 
            # so that anyother non-default values are propagated.             
            new = deepcopy(self)
            new.detector = None
            new.calibration = None
            
            # new = XYDetector()
            # # must define conversion constants here beacuce this is the only other 
            # # variable that is REQUIRED.
            # new.conversion_constant = self.conversion_constant

        local_mask = np.where(
            (self.tth >= range_bounds[0])
            & (self.tth <= range_bounds[1])
            & (self.azm >= azi_bounds[0])
            & (self.azm <= azi_bounds[1])
        )

        new.intensity = deepcopy(self.intensity[local_mask])
        new.tth = deepcopy(self.tth[local_mask])
        new.azm = deepcopy(self.azm[local_mask])
        if "dspace" in dir(self):
            new.dspace = deepcopy(self.dspace[local_mask])

        # set nee range.
        new.tth_start = range_bounds[0]
        new.tth_end = range_bounds[1]

        if "x" in dir(self):
            if self.x is not None:
                new.x = deepcopy(self.x[local_mask])
        if "y" in dir(self):
            if self.y is not None:
                new.y = deepcopy(self.y[local_mask])
        if "z" in dir(self):
            if self.z is not None:
                new.z = deepcopy(self.z[local_mask])

        return new

    def get_calibration(self, file_name=None, settings=None):
        """
        Opens the file with the calibration data in it and updates
        XYDetector.calibration and
        XYDetector.conversion_constant

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
        # make empty calibration
        self.calibration = {
            # "x_dim":   0,
            # "x":       [0,1],
            # "x_start": None, # begining of the 'tth' axis
            # "x_end":   None, # end of the 'tth' axis
            # "y":       [0,1],
            # "y_start": None, # begining of the 'azimith' axis
            # "y_end":   None, # end of the 'azimith' axis
            # "conversion_constant": 1,
            # "max_shape": [0,0]
        }

        # parms_file is file
        if settings != None:
            temp_calib = settings.calibration_parameters
        else:
            parms_file = file_name
            with open(parms_file, "wb") as f:
                temp_calib = pickle.dump("calibration", f)

        # pass temp_calib into the calibration
        for key, val in temp_calib.items():
            self.calibration[key] = val

        if "conversion_constant" in self.calibration:
            self.conversion_constant = self.calibration["conversion_constant"]
        else:
            self.conversion_constant = False
        # self.azm_start = self.calibration["y_start"]
        # self.azm_end = self.calibration["y_end"]

    def get_detector(
        self, settings=None, calibration_file=None, diffraction_data=None, debug=False
    ):
        """
        :param diff_file:
        :param settings:
        """

        if self.calibration == None:
            self.get_calibration(
                settings=settings, file_name=calibration_file, #, debug=debug
            )

        if diffraction_data is None and settings is not None:
            if settings.calibration_data is None:
                diffraction_data = settings.subfit_filename
            else:
                diffraction_data = settings.calibration_data

        self.detector = OrthogonalDetector(
            calibration=self.calibration, diffraction_data=diffraction_data, debug=debug
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
        # check inputs
        if image_name == None and settings.subpattern == None:
            raise ValueError("Settings are given but no subpattern is set.")

        if self.detector == None:
            self.get_detector(settings)

        if image_name == None:
            # load the data for the chosen subpattern.
            image_name = settings.subfit_filename

        if isinstance(image_name, list):
            # then it is a h5 type file
            im = h5_functions.get_images(image_name)
        elif (
            os.path.splitext(image_name)[1] == ".txt"
            or os.path.splitext(image_name)[1] == ".csv"
        ):
            # load csvs or txt file as an image.
            im = csv_to_image(image_name)
        else:
            im = image.imread(image_name)
            im = im.data

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

        if self.calibration["x_dim"] != 0:
            im = im.T
            
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
            ax.imshow(im)
            plt.title(IO_functions.title_file_names(image_name=image_name))
            plt.show()
            plt.close()

        """
        # FIXME: this was a smoothing filter imported to sort an LLNL data set.
        # ideally it would be set in a separate image_preprocess file along with the mask
        # and cosmic preproseccing of the images.
        # The input file would then contain a preparation function something like this:
        # data_prepare = {"smooth": {"Gaussian": 2},
        #                "mask": Calib_mask,
        #                "order": {"smooth", "mask"}}
        # end FIXME
        from skimage import filters
        # smooth_mean = ndi.correlate(bright_square, mean_kernel)
        sigma = 2
        # smooth = filters.gaussian(bright_square, sigma)
        im = filters.gaussian(im, sigma)
        """

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


    def fill_data(
        self, diff_file=None, settings=None, mask=None, make_zyx=False, debug=False
    ):
        """
        Initiates the data arrays.
        Creates the intensity, horizontal (two theta for diffraction),
        vertical (azimuth for diffraction) and converted (d-spacing for
        diffraction) arrays from the Detector and the data.

        If data is provided this makes the intensity array otherwise
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
            raise ValueError("No diffreaction file or settings have been given.")

        if mask == None:
            if settings.calibration_mask:  # in settings.items():
                mask = settings.calibration_mask

        if settings.reduce_by is not None:
            self.reduce_by = settings.reduce_by
            
        if self.detector == None:
            self.get_detector(settings=settings)

        # get the intensities (without mask)
        self.intensity = self.import_image(diff_file, reduce_by=False)
        self.tth = ma.array(self.detector.get_horizontal())
        self.azm = ma.array(self.detector.get_vertical())
        # self.dspace = self._get_d_space()

        # get and apply mask
        mask_array = self.get_mask(mask, self.intensity)
        self.mask_apply(mask_array, debug=debug)

        if self.reduce_by is not None:
            self.intensity = self._reduce_array(self.intensity)
            self.tth = self._reduce_array(self.tth)
            if (re.findall("azimuth", self.calibration["y_label"].lower())
                or 
                re.findall("theta", self.calibration["x_label"].lower())
                or
                self.azm_end == self.azm_start + 360):
                # the data would appear to be diffraction data. 
                self.azm = self._reduce_array(self.azm, polar=True)
            else:
                # data is some other sort of data.
                self.azm = self._reduce_array(self.azm, polar=False)

        # get new azm_blocks from detector.
        self.azm_blocks = self.detector.azm_blocks
        
        self.azm_start = np.floor(np.min(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        self.azm_end   =  np.ceil(np.max(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        self.tth_start = np.min(self.tth.flatten())
        self.tth_end   = np.max(self.tth.flatten())
        
        self.Dispersionlabel = self.detector.calibration["x_label"]
        self.DispersionUnits = self.detector.calibration["x_unit"]
        self.Azimuthlabel = self.detector.calibration["y_label"]
        self.AzimuthUnits = self.detector.calibration["y_unit"]
        
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
        required_list = []

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

    # add common functions
    _get_d_space = _AngleDispersive_common._get_d_space
    conversion = _AngleDispersive_common.conversion
    bins = _AngleDispersive_common.bins
    set_limits = _AngleDispersive_common.set_limits
    test_azims = _AngleDispersive_common.test_azims
    GetDataType = _AngleDispersive_common.GetDataType
    duplicate_without_detector = _AngleDispersive_common.duplicate_without_detector
    _reduce_array = _AngleDispersive_common._reduce_array
    """
    FIXME: add more flxibility to conversion    
    XYFunctions does not have to be X-ray diffraction but it could be. To pass a
    empty conversion throgh the function set DataClass.conversion_factor=False.
    But other conversions might be necessary in fiture and so the function may 
    need further generalisationby by accepting a lambda function. 
    """

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
    # this function is added because it requires access to self:
    dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks


class OrthogonalDetector:
    def __init__(
        self,
        calibration=None,
        diffraction_data=None,
        mask=None,
        max_shape=None,
        debug=False,
    ):
        self.calib = calibration
        self.max_shape = max_shape

        # initiate a bunch of defaults. 
        self.calibration = {}

        self.calibration["x_dim"] = 0
        self.calibration["x"] = [0,1]
        self.calibration["x_start"] = np.nan
        self.calibration["x_end"] = np.nan
        self.calibration["x_scale"] = "linear"
        self.calibration["y"] = [0,1]
        self.calibration["y_start"] = np.nan
        self.calibration["y_end"] = np.nan
        self.calibration["y_scale"] = "linear"
        self.calibration["rotation"] = 0  # in degrees

        self.calibration["x_unit"] = "num"
        self.calibration["x_label"] = "pixels"
        self.calibration["y_unit"] = "num"
        self.calibration["y_label"] = "pixels"

        #set a default spacing
        self.azm_blocks = 100
        

        self.calibration["conversion_constant"] = 1

        # if given a calibration then load it.
        if calibration:
            self.load_calibration(
                calibration=calibration, diffraction_data=diffraction_data
            )

    def load_calibration(self, calibration=None, diffraction_data=None):
        """
        Populates the calibration

        Parameters
        ----------
        calibration : dictionary, οπτιοναλ
            Dictionary containing the calibration parameters for the detector.
            Otherwise returns a blank detector

        Raises
        ------
        ValueError
            Unrecognised key in the dictionary.

        Returns
        -------
        None.

        """

        # fill in the calibration if it exists.
        if calibration:
            for key, value in calibration.items():
                if key == "max_shape":
                    pass
                elif key in list(self.calibration.keys()):
                    self.calibration[key] = value
                else:
                    error_str = (
                        "Key '" + key + "' in not recognised as a calibration parameter"
                    )
                    raise ValueError(error_str)

        if "max_shape" in list(self.calibration.keys()):
            self.max_shape = self.calibration["max_shape"]
        elif diffraction_data:
            if (
                os.path.splitext(diffraction_data)[1] == ".txt"
                or os.path.splitext(diffraction_data)[1] == ".csv"
            ):
                self.max_shape = csv_to_image(diffraction_data).shape
            else:
                with Image.open(diffraction_data) as img:
                    self.max_shape = img.size
        if self.calibration["x_dim"] != 0:
            self.max_shape = np.flip(self.max_shape)

        # either "x" or x_start and x_end are allowed. Fill in the other values.
        if self.calibration["x"] == None:
            self.calibration["x"] = (
                self.calibration["x_start"],
                (self.calibration["x_end"] - self.calibration["x_start"])
                / (self.max_shape[self.calibration["x_dim"]] - 1),
            )
        if self.calibration["x_start"] == np.nan:
            self.calibration["x_start"] = self.calibration["x"][0]
            self.calibration["x_end"] = self.calibration["x"][0] + self.calibration[
                "x"
            ][1] * (self.max_shape[self.calibration["x_dim"]] - 1)
        # ditto for y.
        if self.calibration["x_dim"] == 0:
            y_dim = 1
        else:
            y_dim = 0
        
        if "y" in calibration:
            self.calibration["y"] = calibration["y"]
        elif "y_start" in calibration:
            self.calibration["y"][0] = self.calibration["y_start"]
            self.calibration["y"][1] = (self.calibration["y_end"] - self.calibration["y_start"]) / (self.max_shape[y_dim] - 1),

        # print(self.calibration["y_label"].lower(), "azimuth", self.calibration["y_label"].lower() == "azimuth")
        if self.calibration["y_label"].lower() == "azimuth":
            self.azm_blocks = 45
            if "y_unit" not in calibration:
                self.calibration["y_unit"] = "deg"

        if "theta" in self.calibration["x_label"].lower():
            if "$" not in self.calibration["x_label"].lower():
                self.calibration["x_label"] = re.sub(r'\$\\theta\$|\\theta|theta', "$\theta$", self.calibration["x_label"].lower())
                # force string to be a raw string. 
                self.calibration["x_label"] = self.calibration["x_label"].encode('unicode_escape').decode()
            if "x_unit" not in calibration:
                self.calibration["x_unit"] = "deg"
            
        self.calibration_check

    def calibration_check(self):
        # FIXME (SAH, June 2024) This should have a validation proess in here.
        pass

    def get_xy(self):
        """
        Gets the x (i.e. uncalibrated) pixel positions for the data.

        Returns
        -------
        array
            DESCRIPTION.
        array
            DESCRIPTION.
        """

        if self.calibration["x_dim"] == 1:
            nx, ny = self.max_shape
        else:
            nx, ny = self.max_shape

        if self.calibration["x_scale"] == "linear":
            x = np.linspace(0, nx - 1, nx)
        else:
            x = np.logspace(0, nx - 1, nx)
        if self.calibration["y_scale"] == "linear":
            y = np.linspace(0, ny - 1, ny)
        else:
            y = np.logspace(0, ny - 1, ny)
        x_array, y_array = np.meshgrid(y, x)

        if self.calibration["rotation"] != 0:
            c, s = (
                np.cos(self.calibration["rotation"]),
                np.sin(self.calibration["rotation"]),
            )
            R = np.array(((c, -s), (s, c)))

            tmp = np.vstack((x_array, x_array))
            tmp = R @ tmp
            x_array = np.reshape(tmp[:, 0], x_array.shape)
            y_array = np.reshape(tmp[:, 1], y_array.shape)

        return x_array, y_array

    def get_horizontal(self):
        """
        Gets the horizontal (i.e. 'x') values of the data.
        :param mask:
        :return:
        """
        x_array, _ = self.get_xy()
        calibrated_x = np.ones(x_array.shape) * self.calibration["x"][0]
        for i in range(len(self.calibration["x"]) - 1):
            calibrated_x += x_array * self.calibration["x"][i + 1] * (i + 1)

        return calibrated_x

    def get_vertical(self, mask=None):
        """
        Gets the vertical (i.e. 'y') values of the data.
        :param mask:
        :return:
        """
        _, y_array = self.get_xy()
        calibrated_y = np.ones(y_array.shape) * self.calibration["y"][0]
        for i in range(len(self.calibration["y"]) - 1):
            calibrated_y += y_array * self.calibration["y"][i + 1] * (i + 1)
        return calibrated_y
