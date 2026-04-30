#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

__all__ = ["ESRFlvpDetector"]

import glob
import json
import sys
import os
import re
from pathlib import Path
from copy import copy, deepcopy
from importlib.metadata import version

import fabio
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyFAI
from numpy import cos as cos  # used in self.calibration["trans_function"]
from numpy import pi as pi  # used in self.calibration["trans_function"]
from numpy import sin as sin  # used in self.calibration["trans_function"]
from packaging.version import Version

# Logic to support multiple PyFAI versions
if Version(version("pyFAI")).major >= 2025:
    from pyFAI.integrator.azimuthal import AzimuthalIntegrator
else:
    from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from pyFAI.detectors._common import Detector
from pyFAI.goniometer import MultiGeometry

import cpf # need to import whole package to avoind trying to import part of incompletely iniated method (cpf.settings.issettings for _get_metadata)
from cpf.input_types._AngleDispersive_common import _AngleDispersive_common
from cpf.input_types._metadata_common import _metadata_common
from cpf.input_types._Masks import _masks
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive
from cpf.util.logging import get_logger
from cpf import h5_functions

logger = get_logger("cpf.input_types.ESRFlvpFunctions")


"""
25th April 2024

The funtion names need some reordering/rejigging I think.

Current useage in xrd fit pattern:
    data.filldata -- loads everything. Makes detector, imports data, adds mask, returns everything.
    data.import_image -- replaces the intensity information without touching anything else (at least it should)
        sets:
            self.get_calibration
            self.get_detector
            self.intensity
            self.tth
            self.azm
            self.dspace
            self.azm_start
            self.azm_end

    data.import_image -- replaces all the intensity data with thatof the new image
Also exists:


    data.get_calibration -- loads calibration parameters.
        sets self.calibration
        sets self.conversion_constant
    data.get_detector --
        sets data.detector - which is a detector class. Requires a calibration
    data.get_masked_calibration -- reads the calibration image and applies a mask to it.
        returns masked intensity array
    data.get_mask -- read mask file
        sets data.original_mask
    data.set_mask
    data.mask_resore.

My problem is that the function names are not clear. And that the functionality overlaps somewhat.

I think that we should have a simplifed structure.
Something along the following lines:
    data = settings.data_class
        initalises a blank class.

    data.set_calibration(calibration from settings class)
        sets data.calibration

    data.set_detector(calibraiton = ??)
        calls set_calibration if not set.
        sets data.detector -- which is an empty detector class

    data.load_data(fname, mask)
        fname optional -- if not called from subfile name
        mask optional -- either use set maks or apply a new one.
        cannot be set without a detector specificed.

    data.set_mask
        applies mask to data

"""
"""
13th June 2024.

This has now been rewritten using new python code from Wilson C. at the ESRF.
The new code allows the construction of a multigeometry detector without having to use
any of the calibration functions. This has simplified the class structure again
and we can live for now with the method names.

The methods that are needed and built from code unique to the ESRF data sets are:
    [the order is the order the methods appear in the file; numbers are the order should appear in.]
    -  1 - __init__
    -  3 - get_calibration(self, file_name=None, settings=None, debug=False)
    -   - _get_pos(self, frame, unit="radians")
    -   - _get_sorted_files(self, file_string, debug=False)
    -  4 - get_detector(self, settings=None, calibration_file=None, calibration_data=None, debug=False)
    -  6 - import_image(self, image_name=None, settings=None, mask=None, dtype=None, debug=False)
    -  7 - fill_data(self, diff_file=None, settings=None, mask=None, debug=False)
    -  2 - get_requirements(self, parameter_settings=None)
    -  5 - detector_check(calibration_data, settings=None)
    -  8 - get_masked_intensity(self, mask, im_ints)
    -  9 - set_mask(self, range_bounds=[-np.inf, np.inf], i_max=np.inf, i_min=-np.inf, mask=None    )
    - 10 - mask_restore(self)

    - get_masked_calibration(self, calibration_data, debug, calibration_mask=None)   --- uncalled. can possibly be removed.

methods that can be removed:
    - set_detector(self, settings=None, file_name=None, debug=False)
    - fill_calibration(self, diff_file, settings=None, debug=False)

methods that should be the same as Dioptas
     - _get_d_space(self, mask=None)     --- might not be needed.

Methods that are replicated with the Dioptas functions are:
    - duplicate
    - conversion(self, tth_in, azm=None, reverse=False)
    - bins(self, orders_class, cascade=False)
    - equalObs(self, x, nbin)                         --- [should be _equalObs]
    - test_azims(self, steps = 360)                   --- is this used? can it be removed?
    - set_limits(self, range_bounds=[-np.inf, np.inf], azm_bounds=[-np.inf, np.inf])



"""

""" 7th June 2025
The 360 files from a single spin of the detector have now been wrapped in a h5 file. 
A continuous scan at the same angle is also wrapped in a h5 file. 
Both need different use cases for the processing of them. 

The h5 file from continuous scan at a constant angle requires:
    - a poni calibration file
    - each image treating as a separate step. 
    - h5 data commands to allow sifting through the scans in time.
    
The h5 file from a spin of the detector needs:
    - a json calibration file
    - the whole scan integrating into a single image/step.
    
How do we distinguish between them?
- The spin has a json calibration and will be fed no h5 instricutions. 
- the time scan at constant angle will need to be def the h5 time series instructions for where to start and stop. 

"""
"""
10th June 2025 
Ok so both of these scan types can appear in the same h5 file. 

therefore need to be able to distinguish between a constant angle scan and a detector spin scan. 

No not true -- a constant angle scan can bt treated as a Dioptas data type, saving only rotation scans for the ESRFlvp data type. 
Good

But still need to determine between them.
"""



"""
mask functions

Need functions that can:
    a. make mask from settings (store as original mask)
    b. make mask of smaller range (azm and tth masks)
    c. apply mask to all arrays
    d. restore original mask
    e. remove all masks


Remove mask functions from these class files and put in separate common file.

"""


class ESRFlvpDetector:
    # For data from LVP at ESRF collected using spinning detector. Crichton et al 2023.
    #  functions to change
    #   -   Load/import data
    #   -   Load mask
    #   -   Load Calibration

    def __init__(self, settings_class=None):
        """
        :param settings_class:
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
        self.azm_start = 0
        self.azm_end = 360
        self.tth_start = None
        self.tth_end = None

        # Moving image plate that be integrated into single image: so contonuous data.
        # basically the data can be treated as having a single claibration.
        self.Dispersion = "Angle"
        self.continuous_azm = True

        self.Dispersionlabel = r"2$\theta$"
        self.DispersionUnits = r"$^\circ$"
        self.Azimuthlabel = r"Azimuth"
        self.AzimuthUnits = r"$^\circ$"
        self.Observationslabel = r"Intensity"
        self.ObservationsUnits = r"counts"
        
        self._default_h5_datakey  = '/*.1/measurement/p900kw/'
        self._default_h5_azimuths = '/*.1/measurement/azim/'
        self._default_h5_iterate = [{"from": 0, "to": 0, "step": 1, "label":["pos"], "do":"iterate"},
                           {"do":"combine", 
                 "from": 0, 
                 "to": -1, 
                 "step": 1, 
                 "using":"position",
                 "label": ['/*.1/measurement/azim/'],
                 # "pos": '/*.1/measurement/azim/',
                 "dim": 0}]
        
        self.mask_default = {"threshold": [1, np.inf]}
                
        self.azm_blocks = 2
        # default blocks are 2 degrees incase using only a single detector position
        # if the detector is being spun then the blocks are changed to a larger value.

        self.reduce_by = None
        
        self._default_metadata_labels_hdf5  = {"time_label": '/*.1/measurement/epoch_trig/', # time stamps in ESRF edf file.
                                      "exposure_label": '/*.1/measurement/timer_period/', # exposure times
                                      }        
        self._default_metadata_labels_edf = {"time_label": "time_of_day", # time stamps in ESRF edf file.
                                      "exposure_label": "acq_expo_time", # exposure times
                                      }

        self.calibration = None
        self.conversion_constant = None
        self.detector = None

        if settings_class:
            if settings_class.calibration_parameters != None:
                self.get_calibration(settings=settings_class)
            if self.calibration:
                self.detector = self.get_detector(settings=settings_class)

    def duplicate(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf], with_detector=True, as_masked=None):
        """
        Makes an independent copy of a ESRFlvpDetector Instance.

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
        new : ESRFlvpDetector Instance.
            Copy of ESRFlvpDetector with independedent data values

        """

        #validate the ranges
        range_bounds, azi_bounds = self.check_bounds(range_bounds, azi_bounds)
        
        if with_detector:
            new = copy(self)
        else:
            # cannot deepcopy ESRFlvpDetector instance. So make a new one and fill.
            new = ESRFlvpDetector()
            new.conversion_constant = self.conversion_constant
            new.azm_start = self.azm_start
            new.azm_end = self.azm_end
            new.Dispersionlabel = self.Dispersionlabel
            new.DispersionUnits = self.DispersionUnits
            new.Azimuthlabel = self.Azimuthlabel
            new.AzimuthUnits = self.AzimuthUnits
            new.Observationslabel = self.Observationslabel
            new.ObservationsUnits = self.ObservationsUnits
            new.azm_blocks = self.azm_blocks

        # set new range.
        new.tth_start = range_bounds[0]
        new.tth_end = range_bounds[1]

        # restrict the data. 
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
            if self.dspace is not None:
                new.dspace = deepcopy(self.dspace[local_mask])

        if "x" in dir(self):
            if self.x is not None:
                new.x = deepcopy(self.x[local_mask])
        if "y" in dir(self):
            if self.y is not None:
                new.y = deepcopy(self.y[local_mask])
        if "z" in dir(self):
            if self.z is not None:
                new.z = deepcopy(self.z[local_mask])

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

    def get_calibration(self, file_name=None, settings=None, debug=False):
        """
        Opens the file containing the calibration data and updates
        ESRFlvpClass.calibration and
        ESRFlvpClass.conversion_constant

        Either file_name or settings class are required.

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

        with open(parms_file, "r") as f:
            self.calibration = json.load(f)

        self.conversion_constant = self.calibration["wavelength"] * 1e10  # in angstroms

    def _get_pos(self, frame, unit="radians"):
        """
        Takes the name of a single date file name and extracts the orientation
        of the detector from the filename. The returned angle is within the maximum
        and minimum allowed azimuths (azm_start and azm_end).

        Parameters
        ----------
        frame : string
            Name of the data file.
        unit : string, optional
            Units to return orientation in (either radians or degrees). The default is "radians".

        Returns
        -------
        pos : number
            Orientation of the detector from the file name.
        """
        pos = float(frame.split(".")[-2].split("_")[-1]) + 0.5
        if pos > self.azm_end:
            pos = pos - 360
        elif pos <= self.azm_start:
            pos = pos + 360
        if unit.find("deg") == -1:
            pos = np.deg2rad(pos)
        return pos

    def _get_sorted_files(self, file_string, reduce_by=None, debug=False):
        """
        Sort the glob string for the files and returns a sorted list. The sorting
        is perfomed using the position of the detectors. These postions are
        returned also (in degrees).

        Parameters
        ----------
        file_string : string
            String containing the glob needed to list all the files in the data collection.

        Returns
        -------
        files_list : list
            List of files returned by the glob, ordered by ascending rotation position.
        positions : array
            Array of detector positions in degrees.

        """
        # load the list of files
        # print("file_string", file_string)
        # print(file_string)
        if isinstance(file_string, list) and os.path.splitext(os.path.basename(file_string[0]))[1] == ".h5":
            #define where data locations are in the initaition of the class.
            
            #file string is a list of format 
            # [ file name, h5 kiy, list of frames wanteed, label]
            
            files_list = h5_functions.get_image_keys_new(str(file_string[0]), self.h5_datakey, self.h5_iterate)
            azm_list = h5_functions.get_image_keys_new(str(file_string[0]), self.h5_azimuths, self.h5_iterate)
            # print("for get images", [str(file_string[0])] + azm_list[0])
            positions = h5_functions.get_images([str(file_string[0])] + azm_list[0])
            
            # print("reduce_by", reduce_by, self.reduce_by)
            if reduce_by is not None or self.reduce_by is not None:
                if reduce_by is False:
                    # used to allow the full data image to be read by data_fill as part of reading the calibrations
                    pass
                elif reduce_by is not None and reduce_by != 1:
                    if reduce_by < 1:
                        reduce_by = 1/reduce_by
                    keep = np.int_(np.linspace(0,len(file_string[2]),int(np.floor((len(file_string[2]))/reduce_by)), endpoint=False))
                    files_list = [file_string[2][i] for i in keep]
                    positions = positions[keep]
                    
                elif self.reduce_by is not None and self.reduce_by != 1:
                    if self.reduce_by < 1:
                        self.reduce_by = 1/self.reduce_by
                    keep = np.int_(np.linspace(0,len(file_string[2]),int(np.floor((len(file_string[2]))/self.reduce_by)), endpoint=False))
                    files_list = [file_string[2][i] for i in keep]
                    positions = positions[keep]
                else:
                    # print("ran through here")
                    pass            

            
            # print("files_list", files_list)
            # print("azm_list", azm_list)
            # print("positions", positions)
            
            """
            # reduce the size of the data set (if called for) by skipping over images.
            # this will only work for none h5 ESRF data sets.
            # FIXME: need to adjuct reduction when get h5 data sets -- excpt this function will not be called for h5 data.
            if reduce_by is not None or self.reduce_by is not None:
                if reduce_by is False:
                    # used to allow the full data image to be read by data_fill as part of reading the calibrations
                    pass
                elif reduce_by is not None and reduce_by != 1:
                    if reduce_by < 1:
                        reduce_by = 1/reduce_by
                    keep = np.int_(np.linspace(0,len(files_list),int(np.floor((len(files_list))/reduce_by)), endpoint=False))
                    files_list = [files_list[i] for i in keep]
                    positions = [positions[i] for i in keep]
                    
                elif self.reduce_by is not None and self.reduce_by != 1:
                    if self.reduce_by < 1:
                        self.reduce_by = 1/self.reduce_by
                    keep = np.int_(np.linspace(0,len(files_list),int(np.floor((len(files_list))/self.reduce_by)), endpoint=False))
                    files_list = [files_list[i] for i in keep]
                    positions = [positions[i] for i in keep]
                else:
                    pass            
            """
            
            
            
            
            # positions do not need to be adjusted becuase they are correct from the h5 file.
            # it is the separate (tiff/edf) image file positions that need to be moved relative to the name
        else:
            files_list = glob.glob(str(file_string))
            positions = []
            for i in range(len(files_list)):
                positions.append(self._get_pos(files_list[i], unit="degrees"))
    
            order = np.array(np.argsort(positions))
            positions = [positions[i] for i in order]
            files_list = [files_list[i] for i in order]
            if not files_list:
                raise ValueError("No image files are found")
                
                # the data files are absent. Issue a major warming and assume 
                # that there should be 360 files
                # logger.warning("The detector files are missing. Assume there are 360 files and proceed.")
                # positions = np.linspace(0, 360, 361)
                # files_list = list(positions)
    
            # reduce the size of the data set (if called for) by skipping over images.
            # this will only work for none h5 ESRF data sets.
            # FIXME: need to adjuct reduction when get h5 data sets -- excpt this function will not be called for h5 data.
            if reduce_by is not None or self.reduce_by is not None:
                if reduce_by is False:
                    # used to allow the full data image to be read by data_fill as part of reading the calibrations
                    pass
                elif reduce_by is not None and reduce_by != 1:
                    if reduce_by < 1:
                        reduce_by = 1/reduce_by
                    keep = np.int_(np.linspace(0,len(files_list),int(np.floor((len(files_list))/reduce_by)), endpoint=False))
                    files_list = [files_list[i] for i in keep]
                    positions = [positions[i] for i in keep]
                    
                elif self.reduce_by is not None and self.reduce_by != 1:
                    if self.reduce_by < 1:
                        self.reduce_by = 1/self.reduce_by
                    keep = np.int_(np.linspace(0,len(files_list),int(np.floor((len(files_list))/self.reduce_by)), endpoint=False))
                    files_list = [files_list[i] for i in keep]
                    positions = [positions[i] for i in keep]
                else:
                    pass
            # positions = np.deg2rad(positions)


        return files_list, positions

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
        # Parts of this code are copied from the juypiter notebooks associated with Chriton et al., 2023.

        if self.calibration == None:
            self.get_calibration(
                settings=settings, file_name=calibration_file, debug=debug
            )

        # set the file to use for the shape of the data.
        if diffraction_data is not None:
            calib_frames = diffraction_data
        elif settings.calibration_data is not None:
            calib_frames = settings.calibration_data
        elif settings.image_list is not None:
            # if settings.image_number == 1:
            #     # there is only 1 file and so take the entire list
            #     calib_frames = settings.image_list
            # else:
            # stop
            calib_frames = settings.image_list[0]
        else:
            raise ValueError("There is no data to construct the detector from.")

        # print("calib_frames", calib_frames)
        
        # copy h5 settings in to self if they exist. 
        # used to make sure the data class has these properties
        if "h5_datakey" in dir(settings):
            self.h5_datakey = settings.h5_datakey
        else:
            self.h5_datakey = self._default_h5_datakey
        if "h5_iterate" in dir(settings):
            self.h5_iterate = settings.h5_iterate
        else:
            self.h5_iterate = self._default_h5_iterate
        if "h5_azimuths" in dir(settings):
            self.h5_azimuths = settings.h5_azimuths
        else:
            self.h5_azimuths = self._default_h5_azimuths

        
        # load the list of files
        # if os.path.splitext(os.path.basename(calib_frames))[1] == ".h5":
        #     #define where data locations are in the initaition of the class.
        #     imgs_ = h5_functions.get_image_keys_new(calib_frames, self.h5_data, self.h5_data_iterate)
        #     positions = h5_functions.get_images([[calib_frames, [self.h5_azimuths,list(range(len(imgs_)))]]])
        #     #positions needs to be adjusted to match with the non h5 file types
        #     # positions -= 0.5
        # else:
        imgs_, positions = self._get_sorted_files(calib_frames, debug=False)
        positions = np.deg2rad(positions)
        frames = int(len(positions))
        
        # make list of AzimuthalIntegrator objects for all detector postions
        ais = []
        for i in range(frames):
            pos = positions[i]

            # these lines are required because they are called by the calibration trans_function.
            rot_x = self.calibration["param"][
                self.calibration["param_names"].index("rot_x")
            ]
            rot_y = self.calibration["param"][
                self.calibration["param_names"].index("rot_y")
            ]
            rot1 = self.calibration["param"][
                self.calibration["param_names"].index("rot1")
            ]
            rot2 = self.calibration["param"][
                self.calibration["param_names"].index("rot2")
            ]
            poni1 = self.calibration["param"][
                self.calibration["param_names"].index("poni1")
            ]
            poni2 = self.calibration["param"][
                self.calibration["param_names"].index("poni2")
            ]
            dist = self.calibration["param"][
                self.calibration["param_names"].index("dist")
            ]

            rot1_expr = eval(self.calibration["trans_function"]["rot1_expr"])
            rot2_expr = eval(self.calibration["trans_function"]["rot2_expr"])
            poni1_expr = eval(self.calibration["trans_function"]["poni1_expr"])
            poni2_expr = eval(self.calibration["trans_function"]["poni2_expr"])
            dist_expr = eval(self.calibration["trans_function"]["dist_expr"])

            # make pyFAI AzimuthalIntegrator object for each detector position and append.
            # edited from ESRP code - which edited AzimuthalIntegrator properties and is comparatively very slow.
            # makeing a new AzimuthalIntegrator each time is 100s-1000s of times faster.
            my_ai = AzimuthalIntegrator(
                detector=self.calibration["detector"],
                wavelength=self.calibration["wavelength"],
                dist=dist_expr,
                poni1=poni1_expr,
                poni2=poni2_expr,
                rot1=rot1_expr,
                rot2=rot2_expr,
                rot3=pos,
            )
            ais.append(my_ai)

        # create the multigeometry detector object.
        self.detector = MultiGeometry(
            ais,
            unit="2th_deg",
            radial_range=(1, 13),
            azimuth_range=(self.azm_start, self.azm_end),
        )

        logger.moreinfo(" ".join(map(str, ["Detector is: %s" % self.detector])))
        if logger.is_below_level(level="DEBUG"):
            # plot all the positions of the AzimuthalIntegrators.
            p_angles = []
            p_dist = []
            p_poni1 = []
            p_poni2 = []
            p_rot1 = []
            p_rot2 = []
            p_rot3 = []

            fig, ax = plt.subplots(6, figsize=(9, 9))

            p_angles = np.array(
                [np.rad2deg(ais[frame].rot3) for frame in range(frames)]
            )
            p_dist = np.array([ais[frame].dist for frame in range(frames)])
            p_poni1 = np.array([ais[frame].poni1 for frame in range(frames)])
            p_poni2 = np.array([ais[frame].poni2 for frame in range(frames)])
            p_rot1 = np.array([ais[frame].rot1 for frame in range(frames)])
            p_rot2 = np.array([ais[frame].rot2 for frame in range(frames)])
            p_rot3 = np.array(np.rad2deg([ais[frame].rot3 for frame in range(frames)]))

            ax[0].plot(p_angles, p_dist, marker=".", ls="--")
            ax[1].plot(p_angles, p_poni1, marker=".", ls="--")
            ax[2].plot(p_angles, p_poni2, marker=".", ls="--")
            ax[3].plot(p_angles, p_rot1, marker=".", ls="--")
            ax[4].plot(p_angles, p_rot2, marker=".", ls="--")
            ax[5].plot(p_angles, p_rot3, marker=".", ls="--")

            ax[0].set_ylabel("Distance (m)")
            ax[1].set_ylabel("Poni1 (m)")
            ax[2].set_ylabel("Poni2 (m)")
            ax[3].set_ylabel("Rot1 (rad)")
            ax[4].set_ylabel("Rot2 (rad)")
            ax[5].set_ylabel("Rot3 (°)")

            ax[-1].set_xlabel("Rot3 (°)")
            plt.tight_layout()



    def _read_frames(self, frames, dtype, reduce_by=None, return_metadata=False):
        """
        Reads iamge frames and their metadata from edf images.

        Parameters
        ----------
        frames : TYPE
            DESCRIPTION.
        dtype : TYPE
            DESCRIPTION.
        reduce_by : TYPE, optional
            DESCRIPTION. The default is None.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        
        imagedata = []
        md_tmp = []
        for frame in frames:
            with fabio.open(frame) as f:
                imagedata.append(self._reduce_array(np.array(f.data, dtype=dtype), reduce_by=reduce_by))
                f_without_data = f
                f_without_data.data = None
                md_tmp.append(f_without_data)
        imagedata = np.flipud(imagedata)
        if return_metadata:
            return imagedata, md_tmp
        else:
            return imagedata
        # return ma.array(
        #     [np.flipud(
        #         self._reduce_array(np.array(fabio.open(f).data, dtype=dtype), reduce_by=reduce_by)
        #         ) for f in frames]
        # )


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

        # read image
        if isinstance(image_name, list):
            # then it is a h5 file containing data from a the spin of the detector.
                        
            if self.intensity is None: 
                # Convert the input data from integer to float because the lmfit model values
                # inherits integer properties from the data.
                #
                # Allow option for the data type to be set.
                if dtype == None:
                    if self.intensity.size > 2:
                        # self.intensity has been set before. Inherit the dtype.
                        dtype = self.intensity.dtype
                    else:
                        tmp_list = image_name
                        tmp_list[3] = list(tmp_list[3][0])
                        tmp_image = ma.array(h5_functions.get_images(tmp_list))
                        dtype = self.GetDataType(tmp_image[0], minimumPrecision=False)                
            else:
                # inherit the data type from previosuly.
                dtype = self.intensity.dtype
            
            # print("image_name", image_name)
            self.intensity.data[:] = self._reduce_array(h5_functions.get_images(image_name).astype(dtype), keep_FirstDim=False)
            # flip along axis 2 because elements are otherwise upside down.
            for i in range(self.intensity.shape[0]):
                self.intensity.data[i,:,:] = np.flipud(self.intensity.data[i,:,:])
            
            self._set_metadata(None, settings=settings)
            
        elif os.path.splitext(os.path.basename(image_name))[1] == ".h5":
            # then it is a h5 file containing data from a the spin of the detector.
            image_list = h5_functions.get_image_keys_new(str(image_name), self.h5_data, self.h5_data_return)
                        
            if self.intensity is None: 
                # Convert the input data from integer to float because the lmfit model values
                # inherits integer properties from the data.
                #
                # Allow option for the data type to be set.
                if dtype == None:
                    if self.intensity.size > 2:
                        # self.intensity has been set before. Inherit the dtype.
                        dtype = self.intensity.dtype
                    else:
                        tmp_image = ma.array(h5_functions.get_images([[image_name, [self.h5_data,0]]]))
                        dtype = self.GetDataType(tmp_image[0], minimumPrecision=False)                
            else:
                # inherit the data type from previosuly.
                dtype = self.intensity.dtype
            
            self.intensity.data[:] = self._reduce_array(h5_functions.get_images([[image_name, image_list[0]]]).astype(dtype), keep_FirstDim=True)
            # flip along axis 2 because elements are otherwise upside down.
            for i in range(self.intensity.shape[0]):
                self.intensity.data[i,:,:] = np.flipud(self.intensity.data[i,:,:])
                
            self._set_metadata(None, settings=settings)
            
        else:
            # get ordered list of separate tiff, edf, etc. images
            # reduce the size of the data while listing (if called for)
            frames, detectorangles = self._get_sorted_files(image_name, reduce_by=reduce_by, debug=debug)
            self.detector_check(detectorangles=detectorangles)
            
            if logger.is_below_level(level="MOREINFO"):
                import time
                st = time.time()
    
            if self.intensity is None:
                # Convert the input data from integer to float because the lmfit model values
                # inherits integer properties from the data.
                #
                # Allow option for the data type to be set.
                if dtype == None:
                    if self.intensity.size > 2:
                        # self.intensity has been set before. Inherit the dtype.
                        dtype = self.intensity.dtype
                    else:
                        tmp_image = ma.array(fabio.open(frames[0]).data)
                        dtype = self.GetDataType(tmp_image[0], minimumPrecision=False)
    
                self.intensity, metadata_tmp = self._read_frames(frames, dtype, reduce_by, return_metadata=True)
                
                #make a full size mask and then reduce it if necessary. 
                frame_mask = self._reduce_array(
                    self.get_mask(mask, self._read_frames(frames, dtype, reduce_by=False)), 
                    keep_FirstDim=True
                )
                # resize mask to be the same as the intensity
                self.intensity = ma.array(
                    self.intensity, mask=frame_mask
                )
            else:
                # inherit the data type from previosuly.
                dtype_tmp = self.intensity.dtype
                
                # does not need keep_FirstDim=True because iterating over the frame list which 
                # is already reduced
                self.intensity.data[:], metadata_tmp = self._read_frames(frames, dtype_tmp, reduce_by, return_metadata=True)
            
            self._set_metadata(metadata_tmp, settings=settings)
                
            # 13th June 2024 - Note on flipud: the flipud command is included to invert the short axis of the detector intensity.
            # If I flip the data then the 'spots' in the reconstructed data are spot like, rather than incoherent
            # intensity diffraction peaks.
            # this is possibly due to confusion between clockwise and anti-clockwise data reading between pyFAI
            # ESRF and LVP beamline multidetector objects.
            # I (SAH) do not believe this is the same feature as the Dioptas and Fit2D
            # adjustment required in the Dioptas class.
            if logger.is_below_level(level="MOREINFO"):
                logger.moreinfo(
                    " ".join(map(str, [f"Image import took {time.time()-st} seconds"]))
                )
                # print("image import took", time.time()-st, "seconds")

            # if np.max(angles) - np.min(angles) >= 45:
            #     self.azm_blocks = 45


    def _det_check():
        
        
        if len(np.unique(detectorangles)) == 1:
            raise ValueError("The detector angles are all the same. this is the")


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
        # parts of this methos are copied from ESRF/ID7 (Wilson Crichton's) code.

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

        if mask == None and settings is not None:
            if settings.calibration_mask is not None:  # in settings.items():
                mask = settings.calibration_mask
            else:
                mask = self.mask_default
            
        if settings.reduce_by is not None:
            self.reduce_by = settings.reduce_by
            # self.h5_azimuths_return[1]["step"] = self.reduce_by
            # self.h5_data_iterate[1]["step"] = self.reduce_by
            # self.h5_data_return[1]["step"] = self.reduce_by

        if settings.metadata_labels is not None:
            self.metadata_labels = settings.metadata_labels
        if (isinstance(diff_file, list) or 
              os.path.splitext(os.path.basename(diff_file))[1] == ".h5"):
            self._default_metadata_labels = self._default_metadata_labels_hdf5
            self.metadata_labels = self._default_metadata_labels_hdf5
        else:
            self._default_metadata_labels = self._default_metadata_labels_edf
            self.metadata_labels = self._default_metadata_labels_edf
        # need both metadata_labels and _default_metadata_labels for _get_metadata to work.
            
        if self.detector == None:
            # if reduce_by or self.reduce_by then a reduced list of image files is returned
            self.get_detector(settings=settings)
            
        # set default for data type.
        # Fabio defaults to float64 if nothing is set.
        # Float32 or float16 take up much less memoary than flost64.
        # N.B. Float 16 does not have sufficient precision to be used.
        # Float32 is the minimum
        array_dtype = np.float32

        # get ordered list of images
        frames, detectorangles = self._get_sorted_files(diff_file, reduce_by=self.reduce_by, debug=debug)
        detectorangles = np.deg2rad(detectorangles)
        self.detector_check(calibration_data=diff_file, detectorangles=detectorangles)
        
        # use import_image to fill in the intensity data.
        # No reduction because if reduced then the calibration is nonsence.
        if isinstance(diff_file, list):
            df = diff_file[0]
        else:
            df = diff_file
        if np.size(self.intensity) <= 1:
            if diff_file != None and os.path.splitext(os.path.basename(str(df)))[1] != ".h5":
                # sets self.intensity
                self.import_image(diff_file, settings=settings, mask=mask, dtype=array_dtype, reduce_by=self.reduce_by)
                
            else:
                # empty array
                self.intensity = self._reduce_array(ma.zeros(
                    [
                        len(detectorangles),
                        self.detector.ais[0].detector.shape[0],
                        self.detector.ais[0].detector.shape[1],
                    ],
                    dtype=array_dtype,
                ), keep_FirstDim=True)

                self.import_image(diff_file, settings=settings, mask=mask, dtype=array_dtype, reduce_by=self.reduce_by)
                
                
        # create emmpty arrays
        self.tth = ma.zeros(self.intensity.shape,
            dtype=array_dtype,
        )
        self.azm = ma.zeros(self.intensity.shape,
            dtype=array_dtype,
        )
        if make_zyx:
            # the x, y, z, arrays are not usually needed.
            # but if created can fill the memory.
            # switched incase ever needed/wanted
            self.x = ma.zeros(self.intensity.shape,
                dtype=array_dtype,
            )
            self.y = ma.zeros(self.intensity.shape,
                dtype=array_dtype,
            )
            self.z = ma.zeros(self.intensity.shape,
                dtype=array_dtype,
            )
        logger.debug(
            " ".join(
                map(
                    str,
                    [
                        "self.tth.shape, self.azm.shape (%i, %i, %i)" % self.tth.shape,
                        self.azm.shape,
                    ],
                )
            )
        )

        # fill the arrays
        # print("frames", frames)
        # print(len(self.detector.ais))
        for i in range(len(self.detector.ais)):
            # print(self._reduce_array(np.rad2deg(self.detector.ais[i].twoThetaArray())).shape)
            self.tth[i, :, :] = self._reduce_array(np.rad2deg(self.detector.ais[i].twoThetaArray()))
            self.azm[i, :, :] = self._reduce_array(np.rad2deg(self.detector.ais[i].chiArray()), polar=True, keep_FirstDim=False)
            if make_zyx:
                zyx = self.detector.ais[i].calc_pos_zyx()
                self.z[i, :, :] = self._reduce_array(zyx[0])
                self.y[i, :, :] = self._reduce_array(zyx[1])
                self.x[i, :, :] = self._reduce_array(zyx[2])
        # correct azimuths from pyfai to ID06 refererence frame
        self.azm = (-self.azm+180)
        #force self.azm back to be within azm_start -- azm_end
        self.azm[self.azm < self.azm_start] = self.azm[self.azm < self.azm_start] + 360
        self.azm[self.azm >= self.azm_end]  = self.azm[self.azm >= self.azm_end] - 360
        
        logger.debug(" ".join(map(str, ["Detector is: %s" % self.detector])))
        if logger.is_below_level(level="DEBUG"):
            # plot all the positions of the AzimuthalIntegrators.
            p_angles = []
            p_Chi = []

            fig, ax = plt.subplots(1, figsize=(9, 9))

            frames = int(len(self.detector.ais))

            p_angles = np.array(
                [(self.detector.ais[frame].rot3) for frame in range(frames)]
            )
            p_Chi = np.array(
                (
                    [
                        np.mean(self.detector.ais[frame].chiArray())
                        for frame in range(frames)
                    ]
                )
            )

            ax.plot(p_angles, p_Chi, marker=".", ls="--")
            ax.set_ylabel("Mean ChiArray from azimuthal integrators (°)")
            ax.set_xlabel("Rot3, from azimuthal integrators  (°)")
            plt.title("rot3 vs Chi for EXRF lvp Multigeometry")

        # F IX ME: (June 2024) i dont know that the d-space array is needed.
        # Check for calls and if this is the only one then remove it
        # self.dspace    = self._get_d_space()

        # add masks to arrays
        # FIX ME: should we apply the mask as the arrays are populated rather than here?
        mask_array = self.get_mask(mask, self.intensity)
        self.mask_apply(mask_array, debug=debug)

        self.azm_start = (
            np.floor(np.min(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )
        self.azm_end = (
            np.ceil(np.max(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )
        
        if np.max(self.azm_end) - np.min(self.azm_start) >= 90:
            self.azm_blocks = 45
        elif np.max(self.azm_end) - np.min(self.azm_start) >= 45:
            self.azm_blocks = 22.5
        
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

        For ESRFlvp detectors the default metadata_dictionary is either a 
        a fabio.open(file).header dictionary or (for hdf5 files) a dictionary containing 
        the file name to be read using the hdf5 metadata keys
        
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
        elif (not image_obj and settings) or cpf.settings.is_settings(image_obj):
            # when calling hdf5 file there is no image_obj to send (= None) and settings is
            # provided instead. 
            
            # here we set pointers to the things needed when the metadata is read.
            # assuming that it is not wise (or possible) to list all the possible hdf5
            # keys which could be read as metadata. 
            metadata_dictionary = {}
            metadata_dictionary["image"] = settings.subfit_filename
            metadata_dictionary["note"] = "It is not reasonable to load all hdf5 keys into a dictionary as metadata. Instead carry file name and use keys"
            metadata_dictionary["h5_datakey"] = settings.h5_datakey
            # metadata_dictionary = self._get_file_created_modified(metadata_dictionary, metadata_dictionary["image"][0])
        else:
            metadata_dictionary = {}
            for obj in image_obj:
                # obj is a fabio image instance with the data removed. 
                for j in obj.header:
                    if j not in metadata_dictionary:
                        metadata_dictionary[j] = []
                    try:
                        metadata_dictionary[j].append(float(obj.header.get(j, None)))
                    except:
                        metadata_dictionary[j].append(obj.header.get(j, None))
                
        if settings and "metadata_read_func" in settings.__dict__:
            metadata_dictionary.update(settings.metadata_read_func(settings, image_obj=image_obj))            
        # add the file creation and modifications time
        metadata_dictionary.update(self._get_file_created_modified(metadata_dictionary["image"][0]))
        self.metadata = metadata_dictionary


    def get_requirements(self, parameter_settings=None):
        """
        #Returns the parameters required for this detector class.
        #:return: String parameters
        """

        # FIX ME: this should be called as part of the settings initiation. But it is not.
        # FIX ME: the different bits of this function are not checked and unlikely ot be correct.
        # FIX ME: this function should return the lists of parameters and the checking
        # done in the settings class file.

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
            "Calib_type",  # Must be 'ESRFlvp'
            # "Calib_detector"      # not required -- this is a unique deteor type for ESRF LVP.
            # 'Calib_data',        # Removed because this is not strictly required.
            "Calib_param",  # *.json file the holds the multigeometry calibration
            # 'Calib_pixels',      # this is only needed for GSAS-II and should be read from image file. FIX
            "datafile_directory",  # yes
            "datafile_Basename",  # here this is the file name but needs a * for the detector steps
            "datafile_Ending",  # may contain the * depending where the collection number is specified.
            "datafile_NumDigit",
            "AziBins",  # required based on detector type
            "fit_orders",
            # 'Output_type',		   # should be optional
            # 'Output_NumAziWrite',  # should be optional
            # 'Output_directory']	   # should be optional
        ]

        optional_list = [
            "Calib_mask",  # a mask file is not strictly required.
            "datafile_StartNum",  # now optionally replaced by datafile_Files
            "datafile_EndNum",  # now optionally replaced by datafile_Files
            "datafile_Files",  # optionl replacement for start and end.
        ]

        # Check required against inputs if given
        if parameter_settings is not None:
            # properties of the data files.
            all_present = 1
            for par in parameter_settings:
                if par in required_list:
                    logger.info(" ".join(map(str, [("Got: ", par)])))
                else:
                    logger.warning(
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

    @staticmethod
    def detector_check(calibration_data=None, settings=None, detectorangles=None, stop_on_error=False):
        """
        Function to check if the data is a compound data set from the spinning 
        detector. 
        
        It is is not it either: raises an error or issues an error string. 
        
        
        
        Get detector information
        :param settings:
        :param calibration_data:
        :return: detector:
        """
        if detectorangles is not None:
            if len(np.unique(detectorangles)) == 1:
                err_str = ""
                if stop_on_error:
                    raise ValueError(err_str)
                else:
                    return err_str
        else:
            #need to work this out.
            pass
        
        
        
            
            # FIX ME: SAH, June 2024: I dont know if we need this function or if it is ever called.
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
                detector = Detector(
                    pixel1=sz, pixel2=sz, splineFile=None, max_shape=im_all.shape
                )
            # FIX ME: check the detector type is valid.
            return detector


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
    plot_integrated = _Plot_AngleDispersive.plot_integrated
    plot_fitted = _Plot_AngleDispersive.plot_fitted
    plot_collected = _Plot_AngleDispersive.plot_collected
    plot_calibrated = _Plot_AngleDispersive.plot_calibrated
    plot_integrated = _Plot_AngleDispersive.plot_integrated
    what_plot_type = _Plot_AngleDispersive.what_plot_type

    # this function is added because it requires access to self:
    dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks
