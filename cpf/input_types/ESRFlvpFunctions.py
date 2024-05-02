#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["ESRFlvpDetector"]

import sys
import json
import glob
from copy import deepcopy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import fabio
import pyFAI
from pyFAI.goniometer import GoniometerRefinement, GeometryTransformation
from PIL import Image, ImageDraw
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive


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


class ESRFlvpDetector():
    # For data from LVP at ESRF collected using spinning detector. Crichton et al 2023. 
    #  functions to change
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

        self.x = None
        self.y = None
        self.intensity = None
        self.tth = None
        self.azm = None
        self.dspace = None
        self.azm_start = -180
        self.azm_end   =  180
        self.tth_start = None
        self.tth_end   = None
        # Moving image plate that be integrated into single image: so contonuous data.
        # basically the data can be treated as having a single claibration. 
        self.DispersionType = "AngleDispersive"
        self.continuous_azm = True

        self.azm_blocks = 2.5
        
        self.calibration = None
        self.conversion_constant = None
        self.detector = None
        
        if settings_class:
            self.load_calibration(settings=settings_class)
        
        if self.calibration:
            self.detector = self.load_detector(settings=settings_class)


    def duplicate(self):
        """
        Makes an independent copy of a ESRFlvpDetector Instance
                
        Parameters
        ----------
        None.

        Returns
        -------
        ESRFlvpDetector Instance.

        """
        new = deepcopy(self)
        return new

    
    
    def load_calibration(self, file_name=None, settings=None, debug=False):
        """
        Here for capmatibility with Dioptas functions until I make sure that everything is consistent and works
        """
        self.get_calibration(file_name=file_name, settings=settings)
        
        
        
    def get_calibration(self, file_name=None, settings=None):
        """
        Opens the file with the calibration data in it and updates 
        ESRFlvpClass.Calibration.calibration and 
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



    def load_detector(self, file_name=None, settings=None):
        """
        Here for capmatibility with Dioptas functions until I make sure that everything is consistent and works
        """
        self.get_detector(file_name=file_name, settings=settings)
        
    
    def get_detector(self, settings=None, file_name=None, debug=False):
        """
        Takes the detector information from the settings class, or the 
        calibration file and creates a detector-type instance which converts 
        the x,y,theta of the diffraction data pixels into two theta vs. 
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
        # this is a replaceent function for get_detector and get_calibration used in DioptasFunctions.py. 
        # Parts of this code are copied from the juypiter notebooks associated with Chriton et al., 2023. 
        
        if self.calibration == None:
            self.load_calibration(settings=settings, file_name=file_name, debug=debug)
        
        
        def get_pos(frame):
            return np.deg2rad(float(frame.split(".")[-2].split("_")[-1]) + 0.5)
    
        if settings != None:
            parms_file = settings.calibration_parameters
        else:
            parms_file = file_name
            
        # with open(parms_file, "r") as f:
        #     self.calibration = json.load(f)
            
        # self.conversion_constant = self.calibration["wavelength"]
        
        goniotrans = GeometryTransformation(param_names = ["dist", 
                                                   "poni1",
                                                   "poni2",
                                                   "rot1", 
                                                   "rot2", 
                                                   "rot_x", "rot_y"], #rotation centre 
                                    dist_expr= "dist", 
                                    poni1_expr= "poni1 + ((rot_x*cos(2*pi-(pos))) - rot_y*sin(2*pi-(pos)))",
                                    poni2_expr= "poni2 + ((rot_x*sin(2*pi-(pos))) + rot_y*cos(2*pi-(pos)))",
                                    rot1_expr= "(rot1*cos(2*pi-(pos)) - rot2*sin(2*pi-(pos)))", 
                                    rot2_expr= "(rot1*sin(2*pi-(pos)) + rot2*cos(2*pi-(pos)))",
                                    rot3_expr= "pos")          
        param = {"dist":2.150,
                 "poni1": 0,
                 "poni2":0,
                 "rot1": 0,
                 "rot2": 0,
                 "rot_x": 0,
                 "rot_y": 0
                }
        bounds = {"dist": (1.8,5.0),
                  "poni1": (-0.1, 0.1),
                  "poni2": (0.2, 1.2),
                  "rot1": (-0.1, 0.1),
                  "rot2": (-0.1, 0.1),
                  "rot_x":(-0.1,0.1),
                  "rot_y":(-0.1,0.1)
                 }
        gonioref = GoniometerRefinement(param, #self.calibration.get("param", []), #initial guess
                                        bounds=bounds,
                                        pos_function=get_pos,
                                        trans_function=goniotrans,
                                        #detector=ai.detector, 
                                        #wavelength=ai.wavelength
                                        )

        # needs the json file from ESRF ID7 LVP beamline
        # which is a goniometer valibration V2 file. 
        # gonioref.sload requires a file name. Passing a previous data array (in any format) is not possible
        self.detector = gonioref.sload(parms_file)
        
        

    # @staticmethod
    def import_image(self, image_name=None, settings=None, mask=None, debug=False):
        """
        Import the data image
        If there is no new mask applies the curent mask to the image data. 
        :param mask:
        :param image_name: Name of file
        :param debug: extra output for debugging - plot
        :return: intensity array
        """
        # FIX ME: why is mask in here? should it inherit from previous data if it exists?

        # FIX ME: nxs files need to be checked. But they should be read like h5 files.

    
        #check inputs
        if image_name == None and settings.subpattern == None:
            raise ValueError("Settings are given but no subpattern is set.")
        
        if self.detector == None:
            self.load_detector(settings)

        if image_name == None:
            #load the data for the chosen subpattern.
            image_name = settings.subfit_filename
            
        #copied from Wilson Crichton's code.
        frames = glob.glob(image_name)
        frames.sort()
        #print(len(frames))
        #print(frames[:10])
        angles = [int(frames[f].split(".")[-2].split("_")[-1]) for f in range(len(frames))]
        
        print("here)")
        im = [fabio.open(f).data for f in frames]
        if max(angles) - min(angles) >= 45:
            self.azm_blocks = 45
        print(len(im))
        
        # Convert the input data from integer to float because otherwise the lmfit model
        # inherits Integer properties somewhere along the line.
        
        # here this seems to make arrays that are too large and crash the computer. 
        #im = np.array(im, dtype="f")
        
        
        if mask == None and ma.is_masked(self.intensity) == False:
            self.intensity = ma.array(im)
            print("exit load1")
            return ma.array(im)
        elif mask is not None:
            
            self.intensity = self.get_masked_intensity(mask, im)
            
            #self.intenstiy = [...(self.intensity[f,:,:] for f in range(np.shape(asdf.intensity)[0])]
            #self.intensity = ma.array(im, mask=mask)
            print("exit load2")
            return ma.array(im, mask=mask)
        else:
            self.intensity = ma.array(im, mask=self.intensity.mask)
            
            print("exit load3")
            return ma.array(im)



    def fill_data(
        self, diff_file=None, settings=None, mask=None, debug=False
    ):
        """
        reads all the edf data images in and applies the calibration to them. 
        
        Parameters
        ----------
        diff_file : string, optional
            Filename of the diffraction data files to import.
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
        
        #check inputs
        if diff_file == None and settings == None:
            raise ValueError("No diffreaction file or settings have been given.")
        elif settings != None and settings.subpattern == None:
            raise ValueError("Settings are given but no subpattern is set.")
        if diff_file == None:
            #load the data for the chosen subpattern.
            diff_file = settings.subfit_filename
        
        if self.detector == None:
            self.set_detector(settings)

        import time
        
        st = time.time()

        # use import_image to fill in the intensity data. 
        # (if needed)
        if np.size(self.intensity) <= 1:
            #<=1 replaces self.intesntiy== None as the test. because there is a any/all truth error with the old test
            self.import_image(diff_file, mask=mask)
        
        print(time.time()-st)
        
        #copied from Wilson Crichton's code.
        frames = glob.glob(diff_file)
        frames.sort()
        #print(len(frames))
        #print(frames[:10])
        detectorangles = np.deg2rad([int(frames[f].split(".")[-2].split("_")[-1]) for f in range(len(frames))])
        #print(angles)
        #images = [fabio.open(f).data for f in frames]
        #print(len(images2))
        
        print(time.time()-st)
        
        # NOTE: self.data was called multigeo in ESRF/LVP code. 
        # NOTE: multigeo.ais (or self.data.ais) is a list of azimuthalIntregator objects for each orientation. 
        # this is what the DioptasFunctions 
        filled_detector = self.detector.get_mg(detectorangles)
        #FIX ME: these should be read from the data not fixed.
        filled_detector.radial_range=(0, 13)    #change me to zero centre of res2 image for P2R conversion
        filled_detector.azimuth_range=(self.azm_start, self.azm_end)
        
        
        print(time.time()-st)
        
        #print(multigeo)
        
        #FIX ME: this is probably not needed but it might be a good way of plotting the data rather than using a scatter plot with too many points!!!        
        #self.integrated = self.data.integrate2d(images, npt_rad=5000, npt_azim=360, method=("no", "histogram", "cython"))
        



        #get chi array for 1st detector position. 
        #print(data.ais[0].chiArray())
        
        # z,y,x can be determined from the detector using 
        #self.z, self.y, self.x = filled_detector.ais[0].calc_pos_zyx()
        #FIX ME: needs to be looped over to get all positions.
        
        # FIX ME: Do I collapse the data into a single 1D array or keep it as a 3D array? 
        # Not that it matters outside of this class -- the data should be comparable in all cases.         
        self.tth = []
        self.azm = []
        #self.intensity = []
        #self.x = [] #commented to save memory
        #self.y = [] #commented to save memory
        self.dspace = []
        
        for i in range(len(filled_detector.ais)):
            
            #print(i, "tth", np.shape(filled_detector.ais[i].twoThetaArray()))
            #print(i, "azm", np.shape(filled_detector.ais[i].chiArray()))
            self.tth.append(filled_detector.ais[i].twoThetaArray())
            self.azm.append(filled_detector.ais[i].chiArray()) 
            #self.intensity.append(images[i])
            
            # zyx = filled_detector.ais[i].calc_pos_zyx() #commented to save memory
            # self.z.append(zyx[0]) #commented because not needed, even in plots
            # self.y.append(zyx[1]) #commented to save memory
            # self.x.append(zyx[2]) #commented to save memory
        
        print(time.time()-st)
        
        # from sys import getsizeof
        # print(getsizeof(self.tth))
        
        self.tth       = ma.array(np.rad2deg(self.tth), mask=self.intensity.mask)
        self.azm       = ma.array(np.rad2deg(self.azm), mask=self.intensity.mask)
        #self.intensity = np.array(self.intensity) #set in load_image
        #self.z         = np.array(self.z) #commented because not needed, even in plots
        #self.y         = np.array(self.y) #commented to save memory
        #self.x         = np.array(self.x) #commented to save memory
         
        
        print(time.time()-st)
        
        # print(getsizeof(self.tth))
        
        print("now here")
        self.dspace    = self.get_d_space()
        print("done")
               
        self.azm_start = np.around(self.azm.min() / self.azm_blocks) * self.azm_blocks
        self.azm_end = np.around(self.azm.max() / self.azm_blocks) * self.azm_blocks
        
        print(self.azm.min(), self.azm.max(), self.azm_start, self.azm_end)

        
        
    def get_d_space(self, mask=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        Get as twotheta and then convert into d-spacing
        :param mask:
        :return:
        """
        if mask is not None:
            return ma.array(
                self.conversion(self.tth, reverse=False),
                mask=mask,
            )
        else:
            return ma.array(
                self.conversion(self.tth, reverse=False)
            )
      
  

    def conversion(self, tth_in, azm=None, reverse=False):
        """
        Convert two theta values into d-spacing.
        azm is needed to enable compatibility with the energy dispersive detectors
        :param tth_in:
        :param reverse:
        :param azm:
        :return:
        """
        wavelength = self.conversion_constant
        if not reverse:
            # convert tth to d_spacing
            dspc_out = wavelength / 2 / np.sin(np.radians(tth_in / 2))
        else:
            # convert d-spacing to tth.
            # N.B. this is the reverse function so that labels tth and d_spacing are not correct.
            # print(tth_in)
            dspc_out = 2 * np.degrees(np.arcsin(wavelength / 2 / tth_in))
        return dspc_out
        
  
        
    def get_requirements(self, parameter_settings=None):
        """
        #Returns the parameters required for this detector class.
        #:return: String parameters
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
            'Calib_mask',        # a mask file is not strictly required.
            
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
                    print("Got: ", par)
                else:
                    print("The settings file requires a parameter called  '", par, "'")
                    all_present = 0
            if all_present == 0:
                sys.exit(
                    "The highlighted settings are missing from the input file. Fitting cannot proceed until they "
                    "are all present."
                )
        return required_list

    

    def set_detector(self, settings=None, file_name=None, debug=False):
        """
        Takes the detector information from the settings class, or the 
        calibration file and creates a detector-type instance which converts 
        the x,y,theta of the diffraction data pixels into two theta vs. 
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
        print("this is now a dummy function")
        pass
        """
        # this is a replaceent function for get_detector and get_calibration used in DioptasFunctions.py. 
        
        # Parts of this code are copied from the juypiter notebooks associated with Chriton et al., 2023. 
        
        def get_pos(frame):
            return np.deg2rad(float(frame.split(".")[-2].split("_")[-1]) + 0.5)
    
        if settings != None:
            parms_file = settings.calibration_parameters
        else:
            parms_file = file_name
            
        with open(parms_file, "r") as f:
            self.calibration = json.load(f)
            
        self.conversion_constant = self.calibration["wavelength"]
        
        goniotrans = GeometryTransformation(param_names = ["dist", 
                                                   "poni1",
                                                   "poni2",
                                                   "rot1", 
                                                   "rot2", 
                                                   "rot_x", "rot_y"], #rotation centre 
                                    dist_expr= "dist", 
                                    poni1_expr= "poni1 + ((rot_x*cos(2*pi-(pos))) - rot_y*sin(2*pi-(pos)))",
                                    poni2_expr= "poni2 + ((rot_x*sin(2*pi-(pos))) + rot_y*cos(2*pi-(pos)))",
                                    rot1_expr= "(rot1*cos(2*pi-(pos)) - rot2*sin(2*pi-(pos)))", 
                                    rot2_expr= "(rot1*sin(2*pi-(pos)) + rot2*cos(2*pi-(pos)))",
                                    rot3_expr= "pos")          
        param = {"dist":2.150,
                 "poni1": 0,
                 "poni2":0,
                 "rot1": 0,
                 "rot2": 0,
                 "rot_x": 0,
                 "rot_y": 0
                }
        bounds = {"dist": (1.8,5.0),
                  "poni1": (-0.1, 0.1),
                  "poni2": (0.2, 1.2),
                  "rot1": (-0.1, 0.1),
                  "rot2": (-0.1, 0.1),
                  "rot_x":(-0.1,0.1),
                  "rot_y":(-0.1,0.1)
                 }
        gonioref = GoniometerRefinement(param, #self.calibration.get("param", []), #initial guess
                                        bounds=bounds,
                                        pos_function=get_pos,
                                        trans_function=goniotrans,
                                        #detector=ai.detector, 
                                        #wavelength=ai.wavelength
                                        )

        # needs the json file from ESRF ID7 LVP beamline
        # which is a goniometer valibration V2 file. 
        self.detector = gonioref.sload(parms_file)
        
        print(self.detector)
        
        print(self.detector.to_dict())
        
        print("END")
        
        """
        
        
        
        
        
        
    def fill_calibration(self, diff_file, settings=None, debug=False):
        pass
        
    

    
    def get_masked_calibration(self, calibration_data, debug, calibration_mask=None):
        """
        load calibration data/image -- to get intensities for mask
        :param calibration_data:
        :param debug:
        :param calibration_mask:
        :return:
        """

        im = self.import_image(calibration_data)
        # im=im_all.data
        # im = np.array(im)[::-1]

        # im = im_a.data # self.ImportImage(calibration_data, debug=debug)
        intens = ma.array(im.data)
        # create mask from mask file if present. If not make all values valid
        if calibration_mask:
            intens = self.get_mask(calibration_mask, intens)
        else:
            im_mask = np.zeros_like(intens)
            intens = ma.array(intens, mask=im_mask)
            intens = ma.masked_outside(intens, 0, np.inf)
        return intens
    
    
    

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





    def get_masked_intensity(self, mask, im_ints):
        """
        ESFRlvp mask is currently applied to each file (detector position) as it is imported.
        This mask can either be am image of the same size as the detector 
        or a dictionary with exclusions and limits in it. 
        
        Parameters
        ----------
        mask : dictionary or array
            mask for the intensity data. It is either a dictionary or an image. 
        im_ints : array
            intensity array of the detector

        Returns
        -------
        im_ints : array
            masked image array.

        """
        # FIX ME: this is the same mask function as exists in XYFunctions. And should be the same as DioptasFunctions. 
        # it could/shold be removed from here. 
        if isinstance(mask, dict):
                        
            im_mask = np.zeros(np.array(im_ints)[0].shape, 'bool')
            if "polygon" in mask:
                polygons = mask['polygon']
                for i in polygons:
                    img = Image.new('L', im_mask.shape, 0)
                    ImageDraw.Draw(img).polygon(i, outline=1, fill=1)
                    im_mask = im_mask +ma.array(img)
            im_mask = np.repeat(im_mask[np.newaxis, :, :], np.array(im_ints).shape[0], axis=0)

            if "threshold" in mask:
                threshold = mask['threshold']
                im_mask = np.asarray(im_mask) | ma.masked_outside(np.array(im_ints),threshold[0],threshold[1]).mask
            #mask invalid values
            mask2 = ma.masked_invalid(np.array(im_ints)).mask
            #combine masks
            im_mask = np.asarray(im_mask) | np.asarray(mask2)
            im_ints = ma.array(im_ints, mask=im_mask)
            
        else:
            #the mask is an image.
            im_mask = np.array(Image.open(mask))
            im_ints = ma.array(im_ints, mask=im_mask)    
            im_ints = ma.masked_less(im_ints, 0)
        
        self.original_mask = im_mask
        return im_ints
    
    
    
    # def get_mask(self, msk_file, im_ints):
    #     """
    #     Dioptas mask is compressed Tiff image.
    #     Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py
    #     :param msk_file:
    #     :param im_ints:
    #     :return:
    #     """
    #     im_mask = np.array(Image.open(msk_file))
    #     # im_mask = np.array(im_mask)[::-1]
    #     im_ints = ma.array(im_ints, mask=im_mask)

    #     im_ints = ma.masked_less(im_ints, 0)

    #     self.original_mask = im_mask
    #     # TALK TO SIMON ABOUT WHAT TO DO WITH THIS
    #     # SAH!! July 2021: we need a mask function that makes sure all the arrays have the same mask applied to them.
    #     #     == and possibly one that returns just the mask.

    #     """
    #     if debug:
    #         # N.B. The plot is a pig with even 1000x1000 pixel images and takes a long time to render.
    #         if ImTTH.size > 100000:
    #             print(' Have patience. The mask plot will appear but it can take its time to render.')
    #         fig_1 = plt.figure()
    #         ax1 = fig_1.add_subplot(1, 3, 1)
    #         ax1.scatter(ImTTH, ImAzi, s=1, c=(ImInts.data), edgecolors='none', cmap=plt.cm.jet, vmin=0,
    #                     vmax=np.percentile(ImInts.flatten(), 98))
    #         ax1.set_title('All data')

    #         ax2 = fig_1.add_subplot(1, 3, 2)
    #         ax2.scatter(ImTTH, ImAzi, s=1, c=im_mask, edgecolors='none', cmap='Greys')
    #         ax2.set_title('Mask')
    #         ax2.set_xlim(ax1.get_xlim())

    #         ax3 = fig_1.add_subplot(1, 3, 3)
    #         ax3.scatter(ImTTH, ImAzi, s=1, c=(ImInts), edgecolors='none', cmap=plt.cm.jet, vmin=0,
    #                     vmax=np.percentile(ImInts.flatten(), 98))
    #         ax3.set_title('Masked data')
    #         ax3.set_xlim(ax1.get_xlim())
    #         # ax2.colorbar()
    #         plt.show()

    #         plt.close()
    #     """

        # # FIX ME: need to validate size of images vs. detector name. Otherwise, the mask can be the wrong size
        # # det_size = pyFAI.detectors.ALL_DETECTORS['picam_v1'].MAX_SHAPE
        # # FIX ME : this could probably all be done by using the detector class in fabio.
        # return im_ints

    

    def bins(self, orders_class, cascade=False):
        """
        Determine bins to use in initial fitting.
        Assign each data to a chunk corresponding to its azimuth value
        Returns array with indices for each bin and array of bin centroids
        :param orders_class:
        :return chunks:
        :return bin_mean_azi:
        """

        # determine how to divide the data into bins and how many.
        if cascade:
            bt = orders_class.cascade_bin_type
            if bt == 1:
                b_num = orders_class.cascade_number_bins
            else:
                b_num = orders_class.cascade_per_bin
        else:
            bt = orders_class.fit_bin_type
            if bt == 1:
                b_num = orders_class.fit_number_bins
            else:
                b_num = orders_class.fit_per_bin

        # make the bins
        if bt == 0:
            # split the data into bins with an approximately constant number of data.
            # uses b_num to determine bin size
            num_bins = int(np.round(len(self.azm[self.azm.mask == False]) / b_num))
            bin_boundaries = self.equalObs(
                np.sort(self.azm[self.azm.mask == False]), num_bins
            )
        elif bt == 1:
            # split the data into a fixed number of bins
            # uses b_num to determine bin size
            lims = np.array(
                [
                    np.min(self.azm[self.azm.mask == False]),
                    np.max(self.azm[self.azm.mask == False] + 0.01),
                ]
            )
            lims = np.around(lims / 45) * 45
            bin_boundaries = np.linspace(lims[0], lims[1], num=b_num + 1)

        else:
            # issue an error
            err_str = "The bin type is not recognised. Check input file."
            raise ValueError(err_str)

        if 0:  # for debugging
            print(bin_boundaries)
            # print(np.sort(bin_boundaries))
            # create histogram with equal-frequency bins
            n, bins, patches = plt.hist(
                self.azm[self.azm.mask == False], bin_boundaries, edgecolor="black"
            )
            plt.show()
            # display bin boundaries and frequency per bin
            print("bins and occupancy", bins, n)
            print("expected number of data per bin", orders_class.cascade_per_bin)
            print("total data", np.sum(n))

        # fit the data to the bins
        chunks = []
        bin_mean_azi = []
        temp_azimuth = self.azm.flatten()
        for i in range(len(bin_boundaries) - 1):
            start = bin_boundaries[i]
            end = bin_boundaries[i + 1]
            azi_chunk = np.where((temp_azimuth > start) & (temp_azimuth <= end))
            chunks.append(azi_chunk)
            bin_mean_azi.append(np.mean(temp_azimuth[azi_chunk]))

        return chunks, bin_mean_azi

    def equalObs(self, x, nbin):
        """
        get equally populated bins for data set.
        copied from: https://www.statology.org/equal-frequency-binning-python/ on 26th May 2022.

        Parameters
        ----------
        x : TYPE
            data to disperse.
        nbin : TYPE
            number of bins.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        nlen = len(x)
        x = np.sort(x)
        return np.interp(np.linspace(0, nlen, nbin + 1), np.arange(nlen), np.sort(x))

    def test_azims(self, steps = 360):
        """
        Returns equally spaced set of aximuths within possible range.

        Parameters
        ----------
        steps : Int, optional
            DESCRIPTION. The default is 360.

        Returns
        -------
        array
            list of possible azimuths.

        """
        return np.linspace(self.azm_start,self.azm_end, steps+1)



    def set_limits(self, range_bounds=[-np.inf, np.inf]):
        """
        Set limits to data in two theta
        :param range_bounds:
        :param i_max:
        :param i_min:
        :return:
        """
        local_mask = np.where(
            (self.tth >= range_bounds[0]) & (self.tth <= range_bounds[1])
        )
        self.intensity = self.intensity[local_mask]
        self.tth = self.tth[local_mask]
        self.azm = self.azm[local_mask]
        self.dspace = self.dspace[local_mask]

    # def set_azimuth_limits(self, azi_bounds=[-np.inf, np.inf], azi_flatten=True):
    #     """
    #     Set limits to data in azimuth
    #     :param range_bounds:
    #     :param i_max:
    #     :param i_min:
    #     :return:
    #     """
    #     local_mask = np.where((self.azm >= azi_bounds[0]) & (self.azm <= azi_bounds[1]))
    #     self.intensity = self.intensity[local_mask]
    #     self.tth = self.tth[local_mask]
    #     self.azm = self.azm[local_mask]
    #     self.dspace = self.dspace[local_mask]

    #     if azi_flatten:
    #         self.azm = np.mean(self.azm)

    def set_mask(
        self, range_bounds=[-np.inf, np.inf], i_max=np.inf, i_min=-np.inf, mask=None
    ):
        """
        Set limits to data
        :param range_bounds:
        :param i_max:
        :param i_min:
        :return:
        """
        local_mask = ma.getmask(
            ma.masked_outside(self.tth, range_bounds[0], range_bounds[1])
        )
        local_mask2 = ma.getmask(ma.masked_outside(self.intensity, i_min, i_max))
        combined_mask = np.ma.mask_or(ma.getmask(self.intensity), local_mask)
        combined_mask = np.ma.mask_or(combined_mask, local_mask2)
        NoneType = type(None)
        if not isinstance(mask, NoneType):  # or mask.all() != None:
            combined_mask = np.ma.mask_or(combined_mask, mask)
        self.intensity.mask = combined_mask
        self.tth.mask = combined_mask
        self.azm.mask = combined_mask
        self.dspace.mask = combined_mask

    def mask_restore(self):
        """
        Restores the loaded mask.
        If the image data is still the same size as the original.
        """

        print("Restore original mask.")
        self.intensity.mask = self.original_mask
        self.tth.mask = self.original_mask
        self.azm.mask = self.original_mask
        self.dspace.mask = self.original_mask



#these methods are all called from _Plot_AngleDispersive as they are shared with other detector types.
#Each of these methods remains here because they are called by higher-level functions: 
ESRFlvpDetector.plot_masked      = _Plot_AngleDispersive.plot_masked
ESRFlvpDetector.plot_fitted      = _Plot_AngleDispersive.plot_fitted
ESRFlvpDetector.plot_collected   = _Plot_AngleDispersive.plot_collected
ESRFlvpDetector.plot_calibrated  = _Plot_AngleDispersive.plot_calibrated
#this function is added because it requires access to self:
ESRFlvpDetector.dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks

    
    