#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["ESRFlvpDetector"]

import sys
import json
import glob
import re
from copy import deepcopy, copy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import fabio
import pyFAI
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from pyFAI.goniometer import MultiGeometry
from numpy import cos as cos #used in self.calibration["trans_function"]
from numpy import sin as sin #used in self.calibration["trans_function"]
from numpy import pi as pi #used in self.calibration["trans_function"]
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive
from cpf.input_types._AngleDispersive_common import _AngleDispersive_common
from cpf.input_types._Masks import _masks
from cpf.XRD_FitPattern import logger
import cpf.logger_functions as lg


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

class ESRFlvpDetector():
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
        self.azm_start = -180
        self.azm_end   =  180
        self.tth_start = None
        self.tth_end   = None
        # Moving image plate that be integrated into single image: so contonuous data.
        # basically the data can be treated as having a single claibration. 
        self.DispersionType = "AngleDispersive"
        self.Dispersionlabel = r"2$\theta$"
        self.DispersionUnits = r"$^\circ$"
        self.continuous_azm = True

        self.azm_blocks = 2
        # default blocks are 2 degrees incase using only a single detector position
        # if the detector is being spun then the blocks are changed to a larger value.
        
        self.calibration = None
        self.conversion_constant = None
        self.detector = None
        
        if settings_class:
            self.get_calibration(settings=settings_class)
        if self.calibration:
            self.detector = self.get_detector(settings=settings_class)


    def duplicate(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf]):
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
        
        new = copy(self)
        
        local_mask = np.where(
            (self.tth >= range_bounds[0]) & 
            (self.tth <= range_bounds[1]) &
            (self.azm >= azi_bounds[0]) & 
            (self.azm <= azi_bounds[1])
        )
        
        new.intensity = deepcopy(self.intensity[local_mask])
        new.tth       = deepcopy(self.tth[local_mask])
        new.azm       = deepcopy(self.azm[local_mask])
        if "dspace" in dir(self):
            new.dspace = deepcopy(self.dspace[local_mask])
            
        #set nee range.
        new.tth_start = range_bounds[0]
        new.tth_end   = range_bounds[1]
        
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

        self.conversion_constant = self.calibration['wavelength'] * 1e10 #in angstroms



    def _get_pos(self, frame, unit="radians"):
        """
        Takes the name of a single date file name and extracts the orientation 
        of the detector from the name. The returned angle is within the maximum 
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
        pos = float(frame.split(".")[-2].split("_")[-1])+0.5
        if pos>self.azm_end:
            pos=pos-360
        elif pos <= self.azm_start:
            pos = pos+360
        if unit.find('deg') == -1:
            pos = np.deg2rad(pos)
        return pos
    


    def _get_sorted_files(self, file_string, debug=False):
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
        #load the list of files
        files_list = glob.glob(file_string)
        
        positions = []
        for i in range(len(files_list)):
            positions.append(self._get_pos(files_list[i], unit="degrees"))
            
        order = np.array(np.argsort(positions)) 
        positions = [positions[i] for i in order]
        files_list = [files_list[i] for i in order]
        
        if not files_list:
            raise ValueError("No image files are found")
        
        return files_list, positions
        
        
    
    def get_detector(self, settings=None, calibration_file=None, diffraction_data=None, debug=False):
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
            self.get_calibration(settings=settings, file_name=calibration_file, debug=debug)
        
        # set the file to use for the shape of the data.     
        if diffraction_data is not None:
            calib_frames = diffraction_data
        elif settings.calibration_data is not None:
            calib_frames = settings.calibration_data
        elif settings.image_list is not None:
            calib_frames = settings.image_list[0]
        else:
            raise ValueError("There is no data to construct the detector from.")
    
        #load the list of files
        imgs_, positions = self._get_sorted_files(calib_frames, debug=False)
        positions = np.deg2rad(positions)
        frames = int(len(imgs_))

        # make list of AzimuthalIntegrator objects for all detector postions
        ais = []
        for i in range(frames):
            pos = positions[i]

            #these lines are required because they are called by the calibration trans_function.
            rot_x = self.calibration["param"][self.calibration['param_names'].index("rot_x")]
            rot_y = self.calibration["param"][self.calibration['param_names'].index("rot_y")]
            rot1 = (self.calibration["param"][self.calibration['param_names'].index("rot1")])
            rot2 = (self.calibration["param"][self.calibration['param_names'].index("rot2")])
            poni1 = (self.calibration["param"][self.calibration['param_names'].index("poni1")])
            poni2 = (self.calibration["param"][self.calibration['param_names'].index("poni2")])
            dist = self.calibration["param"][self.calibration['param_names'].index("dist")]

            rot1_expr = eval(self.calibration["trans_function"]['rot1_expr']) 
            rot2_expr = eval(self.calibration["trans_function"]['rot2_expr'])
            poni1_expr = eval(self.calibration["trans_function"]['poni1_expr'])
            poni2_expr = eval(self.calibration["trans_function"]['poni2_expr'])
            dist_expr = eval(self.calibration["trans_function"]['dist_expr'])
            
            # make pyFAI AzimuthalIntegrator object for each detector position and append.
            # edited from ESRP code - which edited AzimuthalIntegrator properties and is comparatively very slow.
            # makeing a new AzimuthalIntegrator each time is 100s-1000s of times faster. 
            my_ai = AzimuthalIntegrator(detector=self.calibration["detector"], wavelength=self.calibration["wavelength"],
                        dist = dist_expr,
                        poni1 = poni1_expr,
                        poni2 = poni2_expr,
                        rot1 = rot1_expr,
                        rot2 = rot2_expr,
                        rot3 = pos)
            ais.append(my_ai)    
            
        # create the multigeometry detector object.
        self.detector = MultiGeometry(ais, unit="2th_deg", radial_range=(1, 13), azimuth_range=(self.azm_start, self.azm_end)) 
        

        logger.moreinfo(" ".join(map(str, ["Detector is: %s" % self.detector] )) )
        if lg.make_logger_output(level="DEBUG"):            
            
            #plot all the positions of the AzimuthalIntegrators. 
            p_angles = []
            p_dist = []
            p_poni1=[]
            p_poni2 = []
            p_rot1 = []
            p_rot2 = []
            p_rot3 = []
        
            fig,ax = plt.subplots(6, figsize=(9,9))
            
            p_angles = np.array([np.rad2deg(ais[frame].rot3) for frame in range(frames)])
            p_dist   = np.array([ais[frame].dist for frame in range(frames)])
            p_poni1  = np.array([ais[frame].poni1 for frame in range(frames)])
            p_poni2  = np.array([ais[frame].poni2 for frame in range(frames)])
            p_rot1   = np.array([ais[frame].rot1 for frame in range(frames)])
            p_rot2   = np.array([ais[frame].rot2 for frame in range(frames)])
            p_rot3   = np.array(np.rad2deg([ais[frame].rot3 for frame in range(frames)]))
        
            ax[0].plot(p_angles, p_dist, marker='.', ls='--')
            ax[1].plot(p_angles, p_poni1, marker='.', ls='--')
            ax[2].plot(p_angles, p_poni2, marker='.', ls='--')
            ax[3].plot(p_angles, p_rot1, marker='.', ls='--')
            ax[4].plot(p_angles, p_rot2, marker='.', ls='--')
            ax[5].plot(p_angles, p_rot3, marker='.', ls='--')
        
            ax[0].set_ylabel("Distance (m)")
            ax[1].set_ylabel("Poni1 (m)")
            ax[2].set_ylabel("Poni2 (m)")
            ax[3].set_ylabel("Rot1 (rad)")
            ax[4].set_ylabel("Rot2 (rad)")
            ax[5].set_ylabel("Rot3 (°)")
        
            ax[-1].set_xlabel("Rot3 (°)")
            plt.tight_layout()
            
            

    # @staticmethod
    def import_image(self, image_name=None, settings=None, mask=None, dtype=None, debug=False):
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
    
        #check inputs
        if image_name == None and settings.subpattern == None:
            raise ValueError("Settings are given but no subpattern is set.")
        
        if self.detector == None:
            self.get_detector(settings)

        if image_name == None:
            #load the data for the chosen subpattern.
            image_name = settings.subfit_filename
            
        #get ordered list of images
        frames, angles = self._get_sorted_files(image_name, debug=debug)     
        
        
        # if mask is False and ma.is_masked(self.intensity) == True:
        #     mask =self.intensity.mask
        if lg.make_logger_output(level="moreinfo"):  
            import time
            st = time.time()
            
        if self.intensity is None:
            
            # Convert the input data from integer to float because the lmfit model values
            # inherits integer properties from the data.
            #
            # Allow option for the data type to be set.
            if dtype==None:
                if self.intensity.size > 2:
                    # self.intensity has been set before. Inherit the dtype.
                    dtype = self.intensity.dtype                
                else:
                    tmp_image = ma.array(fabio.open(frames[0]).data)
                    dtype = self.GetDataType(tmp_image[0], minimumPrecision=False)
                
            self.intensity = ma.array([np.flipud(fabio.open(f).data) for f in frames], dtype=dtype)
            self.intensity = ma.array(self.intensity, mask=self.get_mask(mask, self.intensity))
        else:
            #inherit the data type from previosuly.
            dtype_tmp = self.intensity.dtype
            self.intensity.data[:] = ma.array([np.flipud(fabio.open(f).data) for f in frames], dtype = dtype_tmp)
            
        # 13th June 2024 - Note on flipud: the flipud command is included to invert the short axis of the detector intensity. 
        # If I flip the data then the 'spots' in the reconstructed data are spot like, rather than incoherent 
        #intensity diffraction peaks. 
        # this is possibly due to confusion between clockwise and anti-clockwise data reading between pyFAI 
        # ESRF and LVP beamline multidetector objects. 
        # I (SAH) do not believe this is the same feature as the Dioptas and Fit2D 
        # adjustment required in the Dioptas class. 
        if lg.make_logger_output(level="moreinfo"):  
            logger.moreinfo(" ".join(map(str, [f"Image import took {time.time()-st} seconds"] )) )
            # print("image import took", time.time()-st, "seconds")
            
        if max(angles) - min(angles) >= 45:
            self.azm_blocks = 45
        
        # apply mask to the intensity array
        # print("mask", mask)
        # print(ma.is_masked(self.intensity))
        # try:
        #     print("mask 1", self.intensity.mask[0][0])
        # except:
        #     pass
        
        # if mask is not False:
        #     self.intensity = ma.array(self.intensity, mask=self.get_mask(mask, self.intensity))
        
        # if mask is False:
        #     return self.intensity
        # elif mask is not False:
        #     # apply given mask
        #     self.intensity = ma.array(self.intensity, mask=self.fill_mask(mask, self.intensity))
        #     return ma.array(im, mask=mask)
        # else:
        #     #apply mask from intensities
        #     # self.intensity = ma.array(im, mask=self.intensity.mask)
        #     return self.intensity

        # if mask is None and ma.is_masked(self.intensity) == False:
        #     self.intensity = ma.array(im)
        #     return ma.array(im)
        # elif mask is not None:
        #     # apply given mask
        #     self.intensity = ma.array(im, mask=self.fill_mask(mask, im))
        #     return ma.array(im, mask=mask)
        # else:
        #     #apply mask from intensities
        #     self.intensity = ma.array(im, mask=self.intensity.mask)
        #     return ma.array(im)



    def fill_data(
        self, diff_file=None, settings=None, mask=None, make_zyx=False, 
        debug=False
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
        
        #check inputs
        if diff_file != None:
            pass
        elif settings != None: # and diff_file == None
            #load the data for the chosen subpattern.
            if settings.subfit_filename != None:
                diff_file = settings.subfit_filename
            else:
                raise ValueError("Settings are given but no subpattern is set.")
        else:
            raise ValueError("No diffraction file or settings have been given.")
        
        if mask==None and settings is not None:
            if settings.calibration_mask is not None:# in settings.items():
                mask = settings.calibration_mask
        
        if self.detector == None:
            self.get_detector(settings=settings)

        # set default for data type. 
        # Fabio defaults to float64 if nothing is set.
        # Float32 or float16 take up much less memoary than flost64.
        # N.B. Float 16 does not have sufficient precision to be used.
        # Float32 is the minimum
        array_dtype = np.float32

        #get ordered list of images
        frames, detectorangles = self._get_sorted_files(diff_file, debug=debug)  
        detectorangles = np.deg2rad(detectorangles)
        
        # use import_image to fill in the intensity data. 
        if np.size(self.intensity) <= 1:
            if diff_file != None:
                #sets self.intensity
                #self.import_image(diff_file, mask=mask, dtype=array_dtype)
                self.import_image(diff_file, dtype=array_dtype)
            else:
                # empty array
                self.intensity = ma.zeros([len(frames), self.detector.ais[0].detector.shape[0], self.detector.ais[0].detector.shape[1]], dtype=array_dtype)
        
                
        #create emmpty arrays
        self.tth = ma.zeros([len(frames), self.detector.ais[0].detector.shape[0], self.detector.ais[0].detector.shape[1]], dtype=array_dtype)
        self.azm = ma.zeros([len(frames), self.detector.ais[0].detector.shape[0], self.detector.ais[0].detector.shape[1]], dtype=array_dtype)
        if make_zyx:
            # the x, y, z, arrays are not usually needed. 
            # but if created can fill the memory.
            # switched incase ever needed/wanted
            self.x = ma.zeros([len(frames), self.detector.ais[0].detector.shape[0], self.detector.ais[0].detector.shape[1]], dtype=array_dtype)
            self.y = ma.zeros([len(frames), self.detector.ais[0].detector.shape[0], self.detector.ais[0].detector.shape[1]], dtype=array_dtype)
            self.z = ma.zeros([len(frames), self.detector.ais[0].detector.shape[0], self.detector.ais[0].detector.shape[1]], dtype=array_dtype)
        logger.debug(" ".join(map(str, ["self.tth.shape, self.azm.shape (%i, %i, %i)" % self.tth.shape, self.azm.shape] )) )
        
        # fill the arrays
        for i in range(len(self.detector.ais)):
            self.tth[i,:,:] = np.rad2deg(self.detector.ais[i].twoThetaArray())
            self.azm[i,:,:] = np.rad2deg(self.detector.ais[i].chiArray())
            if make_zyx:
                zyx = self.detector.ais[i].calc_pos_zyx() 
                self.z[i,:,:] = zyx[0]
                self.y[i,:,:] = zyx[1] 
                self.x[i,:,:] = zyx[2] 

        logger.debug(" ".join(map(str, ["Detector is: %s" % self.detector] )) )
        if lg.make_logger_output(level="DEBUG"):  
            
            #plot all the positions of the AzimuthalIntegrators. 
            p_angles = []
            p_Chi = []
        
            fig,ax = plt.subplots(1, figsize=(9,9))
            
            frames = int(len(self.detector.ais))
            
            p_angles = np.array([np.rad2deg(self.detector.ais[frame].rot3) for frame in range(frames)])
            p_Chi   = np.array(np.rad2deg([np.mean(self.detector.ais[frame].chiArray()) for frame in range(frames)]))
                
            ax.plot(p_angles, p_Chi, marker='.', ls='--')
            ax.set_ylabel("Mean ChiArray from azimuthal integrators (°)")
            ax.set_xlabel("Rot3, from azimuthal integrators  (°)")
            plt.title("rot3 vs Chi for EXRF lvp Multigeometry")
            
        # from sys import getsizeof
        
        #F IX ME: (June 2024) i dont know that the d-space array is needed. 
        # Check for calls and if this is the only one then remove it
        # self.dspace    = self._get_d_space()
        
        # add masks to arrays
        #FIX ME: should we apply the mask as the arrays are populated rather than here?
        mask_array = self.get_mask(mask, self.intensity)
        
        self.mask_apply(mask_array, debug=debug)
        #FIX ME: should we apply the mask as the arrays are populated rather than 
        # in a separate function?
               
        self.azm_start = np.floor(np.min(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        self.azm_end   =  np.ceil(np.max(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        self.tth_start = np.min(self.tth.flatten())
        self.tth_end   = np.max(self.tth.flatten())
        
        
  
        
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
            "Calib_type",          # Must be 'ESRFlvp'
            #"Calib_detector"      # not required -- this is a unique deteor type for ESRF LVP. 
            # 'Calib_data',        # Removed because this is not strictly required.
            "Calib_param",         # *.json file the holds the multigeometry calibration
            # 'Calib_pixels',      # this is only needed for GSAS-II and should be read from image file. FIX
            
            "datafile_directory",  # yes
            "datafile_Basename",   # here this is the file name but needs a * for the detector steps
            "datafile_Ending",     # may contain the * depending where the collection number is specified.
            "datafile_NumDigit",
            
            "AziBins",  # required based on detector type
            "fit_orders",
            # 'Output_type',		   # should be optional
            # 'Output_NumAziWrite',  # should be optional
            # 'Output_directory']	   # should be optional
        ]
        
        optional_list = [
            'Calib_mask',         # a mask file is not strictly required.
            
            'datafile_StartNum',  # now optionally replaced by datafile_Files
            'datafile_EndNum',    # now optionally replaced by datafile_Files
            "datafile_Files",     # optionl replacement for start and end. 
            
            ]

        # Check required against inputs if given
        if parameter_settings is not None:
            # properties of the data files.
            all_present = 1
            for par in parameter_settings:
                if par in required_list:
                    logger.info(" ".join(map(str, [("Got: ", par)] )) )
                else:
                    logger.warning(" ".join(map(str, [("The settings file requires a parameter called  '", par, "'")] )) )
                    all_present = 0
            if all_present == 0:
                sys.exit(
                    "The highlighted settings are missing from the input file. Fitting cannot proceed until they "
                    "are all present."
                )
        return required_list        
    
    

    @staticmethod
    def detector_check(calibration_data, settings=None):
        """
        Get detector information
        :param settings:
        :param calibration_data:
        :return: detector:
        """
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
            detector = pyFAI.detectors.Detector(
                pixel1=sz, pixel2=sz, splineFile=None, max_shape=im_all.shape
            )
        # FIX ME: check the detector type is valid.
        return detector


# add common function. 
ESRFlvpDetector._get_d_space = _AngleDispersive_common._get_d_space
ESRFlvpDetector.conversion   = _AngleDispersive_common.conversion
ESRFlvpDetector.bins         = _AngleDispersive_common.bins
ESRFlvpDetector.set_limits   = _AngleDispersive_common.set_limits
ESRFlvpDetector.test_azims   = _AngleDispersive_common.test_azims

#add masking functions to detetor class.
ESRFlvpDetector.get_mask     = _masks.get_mask
ESRFlvpDetector.set_mask     = _masks.set_mask
ESRFlvpDetector.mask_apply   = _masks.mask_apply
ESRFlvpDetector.mask_restore = _masks.mask_restore
ESRFlvpDetector.mask_remove  = _masks.mask_remove

#these methods are all called from _Plot_AngleDispersive as they are shared with other detector types.
#Each of these methods remains here because they are called by higher-level functions: 
ESRFlvpDetector.plot_masked      = _Plot_AngleDispersive.plot_masked
ESRFlvpDetector.plot_integrated  = _Plot_AngleDispersive.plot_integrated
ESRFlvpDetector.plot_fitted      = _Plot_AngleDispersive.plot_fitted
ESRFlvpDetector.plot_collected   = _Plot_AngleDispersive.plot_collected
ESRFlvpDetector.plot_calibrated  = _Plot_AngleDispersive.plot_calibrated
ESRFlvpDetector.plot_integrated  = _Plot_AngleDispersive.plot_integrated
#this function is added because it requires access to self:
ESRFlvpDetector.dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks    