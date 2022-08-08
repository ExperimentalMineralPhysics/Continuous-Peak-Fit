#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["DioptasDetector"]

import os
import sys
from copy import deepcopy
import fabio
import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyFAI
from PIL import Image
from matplotlib import gridspec
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from matplotlib import cm, colors


class DioptasDetector:
    # For Dioptas functions to change
    #   -   Load/import data
    #   -   Load mask
    #   -   Load Calibration

    def __init__(self, calibration_parameters=None, fit_parameters=None):
        """
        :param calibration_parameters:
        """
        # self.parameters = self.get_calibration(calibration_parameters)
        self.requirements = self.get_requirements(fit_parameters)
        self.detector = None
        self.intensity = None
        self.tth = None
        self.azm = None
        self.dspace = None
        self.calibration = None
        
        # image plate data to all azimuths are possible -- i.e. continuous   
        self.continuous_azm = True

        if calibration_parameters:
            self.calibration = self.get_calibration(calibration_parameters)

    def duplicate(self):
        """
        Makes a "deep" copy of an Dioptas Instance, using copy.deepcopy()
        """
        
        new = deepcopy(self)
        
        return new

    def fill_data(self, diff_file, settings=None, debug=None):#, calibration_mask=None):
        """
        :param settings: -- now a class.
        :param diff_file:
        :param debug:
        :param calibration_mask: -- now removed because calibration mask file is contained in settings class
        """
        # self.detector = self.detector_check(calibration_data, detector)
        #        print(settings.Calib_param)
        self.calibration = self.get_calibration(settings.calibration_parameters)
        #self.calibration = self.get_calibration(settings.Calib_param)
        self.get_detector(diff_file, settings)

        self.intensity = self.get_masked_calibration(diff_file, debug, calibration_mask=settings.calibration_mask)
        self.tth = self.get_two_theta(mask=self.intensity.mask)
        self.azm = self.get_azimuth(mask=self.intensity.mask)
        self.dspace = self.get_d_space(mask=self.intensity.mask)



    def get_detector(self, diff_file, settings=None):
        """
        :param diff_file:
        :param settings:
        :return:
        """

        #if hasattr(settings, "Calib_detector"):
        if settings.calibration_detector is not None:
            self.detector = pyFAI.detector_factory(settings.calibration_detector)
        else:
            # if settings is None or detector == 'unknown' or detector == 'other' or detector == 'blank':
            im_all = fabio.open(diff_file)
            if "pixelsize1" in self.calibration:
                sz1 = self.calibration["pixelsize1"]
                sz2 = self.calibration["pixelsize2"]
            elif settings.calibration_pixel_size is not None:
                sz1 = settings.calibration_pixel_size * 1e-6
                sz2 = settings.calibration_pixel_size * 1e-6
            else:
                err_str = (
                    "The pixel size or the detector are is not defined. Add to inputs."
                )
                raise ValueError(err_str)

            self.detector = pyFAI.detectors.Detector(
                pixel1=sz1, pixel2=sz2, splineFile=None, max_shape=im_all.shape
            )
        # FIX ME: check the detector type is valid.
        # self.detector = self.detector_check(calibration_data, detector=None)
        # return detector

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
            #sz = calibration_data.Calib_pixels  # Pixel_size
            sz = calibration_data.calibration_pixel_size  # Pixel_size
            if sz > 1:
                sz = sz * 1e-6
            detector = pyFAI.detectors.Detector(
                pixel1=sz, pixel2=sz, splineFile=None, max_shape=im_all.shape
            )
        # FIX ME: check the detector type is valid.
        return detector

    def get_masked_calibration(self, calibration_data, debug, calibration_mask=None):
        """
        load calibration data/image -- to get intensities for mask
        :param calibration_data:
        :param debug:
        :param calibration_mask:
        :return:
        """
        im_all = fabio.open(calibration_data)
        im=im_all.data
        im = np.array(im)[::-1]
        
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

    #@staticmethod
    def import_image(self,image_name, mask=None, debug=False):
        """
        Import the data image
        :param mask:
        :param image_name: Name of file
        :param debug: extra output for debugging - plot
        :return: intensity array
        """
        # FIX ME: why is mask in here? should it inherit from previous data if it exists?
        # im = Image.open(ImageName) ##always tiff?- no
        filename, file_extension = os.path.splitext(image_name)
        if file_extension == ".nxs":
            im_all = h5py.File(image_name, "r")
            # FIX ME: This assumes that there is only one image in the nxs file.
            # If there are more, then it will only read 1.
            im = np.array(im_all["/entry1/instrument/detector/data"])
        else:
            im_all = fabio.open(image_name)
            im = im_all.data
        
        # FIX ME: Copied from MED_functions. 
        # FIX ME: Convert the input data from integer to float because otherwise the model inherits Integer properties somewhere along the line.
        # FIX ME: I dont know if this needed here but given everything else is the same I have copied it from MED_functions.
        im = np.array(im, dtype="f")
        
        # Dioptas flips the images to match the orientations in Plot2D.
        # Therefor implemented here to be consistent with Dioptas.
        im = np.array(im)[::-1]
        
        self.intensity = ma.array(im, mask=self.intensity.mask)
                
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)            
            self.plot_collected(fig_plot=fig, axis_plot=ax)
            plt.title(os.path.split(image_name)[1])
            plt.show()
            plt.close()
            
        
        if mask is not None:
            return ma.array(im, mask=mask)
        else:
            return ma.array(im)


    def conversion(self, tth_in, azm=None, reverse=False):
        """
        Convert two theta values into d-spacing.
        azm is needed to enable compatibility with the energy dispersive detectors
        :param tth_in:
        :param reverse:
        :param azm:
        :return:
        """
        
        wavelength = self.calibration["wavelength"] * 1e10
        if not reverse:
            # convert tth to d_spacing
            dspc_out = wavelength / 2 / np.sin(np.radians(tth_in / 2))
        else:
            # convert d-spacing to tth.
            # N.B. this is the reverse function so that labels tth and d_spacing are not correct.
            # print(tth_in)
            dspc_out = 2 * np.degrees(np.arcsin(wavelength / 2 / tth_in))
        return dspc_out


    def get_mask(self, msk_file, im_ints):
        """
        Dioptas mask is compressed Tiff image.
        Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py
        :param msk_file:
        :param im_ints:
        :return:
        """
        im_mask = np.array(Image.open(msk_file))
        #im_mask = np.array(im_mask)[::-1]
        im_ints = ma.array(im_ints, mask=im_mask)

        im_ints = ma.masked_less(im_ints, 0)

        self.original_mask = im_mask
        # TALK TO SIMON ABOUT WHAT TO DO WITH THIS
        # SAH!! July 2021: we need a mask function that makes sure all the arrays have the same mask applied to them.
        #     == and possibly one that returns just the mask.

        """
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
            ax2.scatter(ImTTH, ImAzi, s=1, c=im_mask, edgecolors='none', cmap='Greys')
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
        """

        # FIX ME: need to validate size of images vs. detector name. Otherwise, the mask can be the wrong size
        # det_size = pyFAI.detectors.ALL_DETECTORS['picam_v1'].MAX_SHAPE
        # FIX ME : this could probably all be done by using the detector class in fabio.
        return im_ints

    def get_calibration(self, file_nam):
        """
        Parse inputs file, create type specific inputs
        :rtype: object
        :param file_nam:
        :return:
        """
        parms_file = open(file_nam, "r")
        file_lines = parms_file.readlines()
        parms_dict = {}
        for item in file_lines:
            new_parms = item.strip("\n").split(":", 1)
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
                    if parm.startswith("["):
                        list_vals = parm.strip("[").strip("]").split(",")
                        new_list = []
                        for val in list_vals:
                            new_value = None
                            try:
                                new_value = int(val)
                                new_list.append(new_value)
                            except ValueError:
                                try:
                                    new_value = float(val)
                                    new_list.append(new_value)
                                except ValueError:
                                    new_list.append(
                                        val.replace("'", "")
                                            .replace(" ", "")
                                            .replace('"', "")
                                    )
                        parms_dict[new_parms[0]] = new_list
                    elif parm.startswith("{"):
                        # print parm
                        list_vals = parm.strip("{").strip("}").split(",")
                        new_dict = {}
                        for key_val in list_vals:
                            # print key_val
                            new_key = (
                                key_val.split(":")[0]
                                    .replace("'", "")
                                    .replace(" ", "")
                                    .replace('"', "")
                            )
                            val = key_val.split(":")[1]
                            new_value = None
                            try:
                                new_value = int(val)
                                new_dict[str(new_key)] = new_value
                            except ValueError:
                                try:
                                    new_value = float(val)
                                    new_dict[str(new_key)] = new_value
                                except ValueError:
                                    new_dict[str(new_key)] = (
                                        val.replace("'", "")
                                            .replace(" ", "")
                                            .replace('"', "")
                                    )
                        parms_dict[new_parms[0]] = new_dict
                    elif not parm:
                        parms_dict[new_parms[0]] = ""

                    else:
                        parms_dict[new_parms[0]] = str(parm)

        # process 'detector_config' into pixel sizes if needed.
        if "Detector_config" in parms_dict:
            # print(parms_dict['Detector_config'])
            parms_dict["pixelsize1"] = parms_dict["Detector_config"]["pixel1"]
            parms_dict["pixelsize2"] = parms_dict["Detector_config"]["pixel2"]
        # get wavelengths in Angstrom
        parms_dict["conversion_constant"] = parms_dict["Wavelength"] * 1e10
        # force all dictionary labels to be lower case -- makes s
        parms_dict = {k.lower(): v for k, v in parms_dict.items()}
        # report type
        parms_dict["DispersionType"] = "AngleDispersive"
        
        # print("parms_dict", parms_dict)
        return parms_dict


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
            if bt==1:
                b_num = orders_class.cascade_number_bins
            else:
                b_num = orders_class.cascade_per_bin
        else:
            bt    = orders_class.fit_bin_type
            if bt==1:
                b_num = orders_class.fit_number_bins
            else:
                b_num = orders_class.fit_per_bin
        
        #make the bins
        if bt == 0:
            # split the data into bins with an approximately constant number of data.      
            # uses b_num to determine bin size
            num_bins = int(np.round(len(self.azm[self.azm.mask == False])/b_num))
            bin_boundaries = self.equalObs(np.sort(self.azm[self.azm.mask == False]), num_bins)
        elif bt==1:
            # split the data into a fixed number of bins     
            # uses b_num to determine bin size            
            lims = np.array([np.min(self.azm[self.azm.mask == False]), np.max(self.azm[self.azm.mask == False]+0.01)])
            lims = np.around(lims / 45) * 45
            bin_boundaries = np.linspace(lims[0],lims[1],num=b_num+1)
            
        else:
            # issue an error
            err_str = (
                "The bin type is not recognised. Check input file."
            )
            raise ValueError(err_str)
            
        if 0: #for debugging
            print(bin_boundaries)
            # print(np.sort(bin_boundaries))
            #create histogram with equal-frequency bins 
            n, bins, patches = plt.hist(self.azm[self.azm.mask == False], bin_boundaries, edgecolor='black')
            plt.show()
            #display bin boundaries and frequency per bin 
            print("bins and occupancy", bins, n)
            print("expected number of data per bin", orders_class.cascade_per_bin)
            print("total data", np.sum(n))
            
        #fit the data to the bins
        chunks = []
        bin_mean_azi = []
        temp_azimuth = self.azm.flatten()
        for i in range(len(bin_boundaries)-1):
            start = bin_boundaries[i]
            end = bin_boundaries[i+1]
            azi_chunk = np.where((temp_azimuth > start) & (temp_azimuth <= end))
            chunks.append(azi_chunk)
            bin_mean_azi.append(np.mean(temp_azimuth[azi_chunk]))
            
        return chunks, bin_mean_azi


    def equalObs(self,x, nbin):
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
        return np.interp(np.linspace(0, nlen, nbin + 1),
                         np.arange(nlen),
                         np.sort(x))
        
    def get_azimuthal_integration(self):
        """

        :return:
        """
        ai = AzimuthalIntegrator(
            self.calibration["distance"],
            self.calibration["poni1"],
            self.calibration["poni2"],
            self.calibration["rot1"],
            self.calibration["rot2"],
            self.calibration["rot3"],
            pixel1=self.calibration["pixelsize1"],
            pixel2=self.calibration["pixelsize2"],
            detector=self.detector,
            wavelength=self.calibration["wavelength"],
        )
        return ai


    def get_two_theta(self, mask=None):
        """
        Give 2-theta value for detector x,y position; calibration info in data
        :param mask:
        :return:
        """

        ai = self.get_azimuthal_integration()
        if mask is not None:
            return ma.array(np.degrees(ai.twoThetaArray()), mask=mask)
        else:
            return ma.array(np.degrees(ai.twoThetaArray()))

    def get_azimuth(self, mask=None):
        """
        Give azimuth value for detector x,y position; calibration info in data
        :param mask:
        :return:
        """
        ai = self.get_azimuthal_integration()
        if mask is not None:
            return ma.array(np.degrees(ai.chiArray()), mask=mask)
        else:
            return ma.array(np.degrees(ai.chiArray()))

    def get_d_space(self, mask=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        Get as twotheta and then convert into d-spacing
        :param mask:
        :return:
        """
        ai = self.get_azimuthal_integration()
        if mask is not None:
            return ma.array(
                self.calibration["wavelength"]
                * 1e10
                / 2
                / np.sin(ai.twoThetaArray() / 2),
                mask=mask,
            )
        else:
            return ma.array(
                self.calibration["wavelength"]
                * 1e10
                / 2
                / np.sin(ai.twoThetaArray() / 2)
            )



    def set_limits(self, range_bounds=[-np.inf, np.inf]):
        """
        Set limits to data in two theta
        :param range_bounds:
        :param i_max:
        :param i_min:
        :return:
        """    
        local_mask = np.where((self.tth >= range_bounds[0]) & (self.tth <= range_bounds[1]))
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
            

    def set_mask(self, range_bounds=[-np.inf, np.inf], i_max=np.inf, i_min=-np.inf, mask=None):
        """
        Set limits to data 
        :param range_bounds:
        :param i_max:
        :param i_min:
        :return:
        """     
        local_mask = ma.getmask(ma.masked_outside(self.tth, range_bounds[0], range_bounds[1]))
        local_mask2 = ma.getmask(ma.masked_outside(self.intensity, i_min, i_max))
        combined_mask = np.ma.mask_or(ma.getmask(self.intensity), local_mask)
        combined_mask = np.ma.mask_or(combined_mask, local_mask2)
        NoneType = type(None)
        if not isinstance(mask, NoneType):# or mask.all() != None:
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
        
    def dispersion_ticks(self, disp_ticks=None, unique=10):
        """
        Returns the labels for the dispersion axis/colour bars.
        
        :param disp_ticks: -- unused for maintained for compatibility with MED functions
        :param unique:     
        :return new_tick_positions:
        """    
        
        if len(np.unique(self.azm)) >= unique:
            disp_lims = np.array(
                [np.min(self.azm.flatten()), np.max(self.azm.flatten())]
            )
            disp_lims = np.around(disp_lims / 180) * 180
            disp_ticks = list(range(int(disp_lims[0]),int(disp_lims[1]+1),45))
            
        return disp_ticks
    

    # def full_plot(self,
    #         plot_type="data",
    #         masked=True,
    #         im_intensity_2=None,
    #         name="",
    # ):
    #     """
    #     :param im_dispersion:
    #     :param im_azimuths:
    #     :param im_intensity:
    #     :param plot_type:
    #     :param masked:
    #     :param im_intensity_2:
    #     :param name:
    #     :return:
    #     """
    #     # Plot data.
    #     # possibilities:
    #     #   1. just plotting the data or model - with or without mask: label 'data'
    #     #   2. plot of data, model and differences: label 'model'
    #     #   3. plot of all data, mask and masked data. label 'mask'
    #     to_plot = []
    #     if plot_type == "data":
    #         x_plots = 1
    #         y_plots = 1
    #         to_plot.append(im_intensity)
    #         plot_mask = [masked]
    #         plot_title = ["Data"]
    #         plot_cmap = [plt.cm.jet]
    #         spec = gridspec.GridSpec(
    #             ncols=x_plots,
    #             nrows=y_plots,
    #             width_ratios=[1],
    #             wspace=0.5,
    #             hspace=0.5,
    #             height_ratios=[1],
    #         )
    #     elif plot_type == "model":
    #         x_plots = 3
    #         y_plots = 1
    #         to_plot.append(im_intensity)
    #         to_plot.append(im_intensity_2)
    #         to_plot.append(im_intensity - im_intensity_2)
    #         plot_mask = [masked, masked, masked]
    #         plot_title = ["Data", "Model", "Residuals"]
    #         plot_cmap = ["jet", "jet", "jet"]
    #     elif plot_type == "mask":
    #         x_plots = 3
    #         y_plots = 2
    #         to_plot.append(im_intensity)
    #         to_plot.append(np.array(ma.getmaskarray(im_intensity), dtype="uint8") + 1)
    #         to_plot.append(im_intensity)
    #         plot_mask = [False, False, True]
    #         plot_title = ["All Data", "Mask", "Masked Data"]
    #         plot_cmap = "jet", "Greys", "jet"
    #         spec = gridspec.GridSpec(
    #             ncols=x_plots,
    #             nrows=y_plots,
    #             width_ratios=[1, 1, 1],
    #             wspace=0.5,
    #             hspace=0.5,
    #             height_ratios=[2, 1],
    #         )
    #     else:
    #         print("Plotting set up incorrect. Continuing without plotting.")
    #         return

    #     y_lims = np.array(
    #         [np.min(im_azimuths.flatten()), np.max(im_azimuths.flatten())]
    #     )
    #     y_lims = np.around(y_lims / 180) * 180

    #     # N.B. The plot is a pig with even 1000x1000 pixel images and takes a long time to render.
    #     if im_intensity.size > 100000:
    #         print(
    #             " Have patience. The plot(s) will appear but it can take its time to render."
    #         )

    #     fig_1 = plt.figure()
    #     for i in range(x_plots):
    #         ax1 = fig_1.add_subplot(spec[i])
    #         if plot_mask[i]:
    #             ax1.scatter(
    #                 im_dispersion,
    #                 im_azimuths,
    #                 s=1,
    #                 c=(to_plot[i]),
    #                 edgecolors="none",
    #                 cmap=plot_cmap[i],
    #                 vmin=0,
    #             )  # plot all data including masked data.
    #         else:
    #             ax1.scatter(
    #                 im_dispersion.data,
    #                 im_azimuths.data,
    #                 s=1,
    #                 c=to_plot[i].data,
    #                 edgecolors="none",
    #                 cmap=plot_cmap[i],
    #                 vmin=0,
    #             )
    #         ax1.set_title(plot_title[i])
    #         ax1.set_ylim(y_lims)
    #         locs, labels = plt.xticks()
    #         plt.setp(labels, rotation=90)
    #         if i == 0:
    #             ax1.set_ylabel(r"Azimuth ($^\circ$)")
    #         ax1.set_xlabel(r"$2\theta$ ($^\circ$)")

    #     if y_plots > 1:
    #         for i in range(x_plots):
    #             ax1 = fig_1.add_subplot(spec[i + x_plots])
    #             if plot_type == "mask" and i == 1:
    #                 # plot cdf of the intensities.
    #                 # sort the data in ascending order
    #                 x1 = np.sort(im_intensity.data)
    #                 x2 = np.sort(im_intensity)

    #                 # get the cdf values of y
    #                 y1 = np.arange(np.size(x1)) / float(np.size(x1))
    #                 y2 = np.arange(np.size(x2)) / float(ma.count(x2))

    #                 # ax1 = fig_1.add_subplot(1, 1, 1)
    #                 ax1.plot(
    #                     x1,
    #                     y1,
    #                 )
    #                 ax1.plot(
    #                     x2,
    #                     y2,
    #                 )
    #                 # ax1.legend()
    #                 ax1.set_title("CDF of the intensities")

    #             else:
    #                 if plot_mask[i]:
    #                     ax1.scatter(
    #                         im_dispersion,
    #                         to_plot[i],
    #                         s=1,
    #                         c=im_azimuths,
    #                         edgecolors="none",
    #                         cmap=plot_cmap[i],
    #                         vmin=0,
    #                     )  # plot all data including masked data.
    #                 else:
    #                     ax1.scatter(
    #                         im_dispersion.data,
    #                         to_plot[i].data,
    #                         s=1,
    #                         c=im_azimuths.data,
    #                         edgecolors="none",
    #                         cmap=plot_cmap[i],
    #                         vmin=0,
    #                     )
    #                 # ax1.set_title(plot_title[i])
    #                 # ax1.set_ylim(y_lims)
    #                 locs, labels = plt.xticks()
    #                 plt.setp(labels, rotation=90)
    #                 if i == 0:
    #                     ax1.set_ylabel(r"Intensity (a.u.)")
    #                 ax1.set_xlabel(r"$2\theta$ ($^\circ$)")

    #     plt.suptitle(name + "; masking")
    #     return fig_1

    # @staticmethod
    # def debug_plot(title_str='Unspecified', two_theta=None, azimuth=None, intensity=None, show=0):
    #     """

    #     :param title_str:
    #     :param two_theta:
    #     :param azimuth:
    #     :param intensity:
    #     :param show:
    #     :return:
    #     """
    #     '''
    #     fig = plt.figure()
    #     if fit_settings.Calib_type == "Med":
    #         # plt.scatter(twotheta, azimuth, s=intens/10, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0,
    #         # vmax = 1000)  #Better for energy dispersive data.  #FIX ME: need variable for vmax.
    #         for x in range(9):
    #             plt.plot(
    #                 twotheta[x], intens[x] + 100 * x
    #             )  # Better for energy dispersive data. #FIX ME: need variable
    #             # for vmax.
    #         plt.xlabel("Energy (keV)")
    #         plt.ylabel("Intensity")
    #     else:
    #     '''
    #     # plt.scatter(twotheta, azimuth, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 1000)
    #     # Better for monochromatic data.
    #     # FIX ME: need variable for vmax.
    #     plt.scatter(
    #         two_theta,
    #         azimuth,
    #         s=4,
    #         c=intensity,
    #         edgecolors="none",
    #         cmap=plt.cm.jet,
    #         vmin=0,
    #         vmax=4000,  # FIX ME: generally applicable?
    #     )
    #     # Better for monochromatic data.
    #     # FIX ME: need variable for vmax.
    #     plt.xlabel(r"2$\theta$")
    #     plt.ylabel("Azimuth")
    #     plt.colorbar()
    #     plt.title(title_str)
    #     plt.draw()
    #     if show:
    #         plt.show()
    #     plt.close()




    def plot_masked(self, fig_plot=None):
        """
        Plot all the information needed to mask the data well. 
        :param fig:
        :return:
        """
        
        x_plots = 3
        y_plots = 2
        spec = gridspec.GridSpec(
            ncols=x_plots,
            nrows=y_plots,
            width_ratios=[1, 1, 1],
            wspace=0.5,
            hspace=0.5,
            height_ratios=[2, 1],
        )
                
        ax1 = fig_plot.add_subplot(spec[0])
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax1, show="unmasked_intensity", x_axis="default", limits=[0, 100])
        ax1.set_title("All Data")
        ax2 = fig_plot.add_subplot(spec[1])
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax2, show="mask", x_axis="default", limits=[0, 100])
        ax2.set_title("Mask")
        ax3 = fig_plot.add_subplot(spec[2])
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax3, show="intensity", x_axis="default", limits=[0, 100])
        ax3.set_title("Masked Data")
        
        
        ax4 = fig_plot.add_subplot(spec[3])
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax4, show="unmasked_intensity", x_axis="default", y_axis="intensity", limits=[0, 100])
        
        
        ax5 = fig_plot.add_subplot(spec[4])
        # plot cdf of the intensities.
        # sort the data in ascending order
        x1 = np.sort(self.intensity.data)
        x2 = np.sort(self.intensity)

        # get the cdf values of y
        y1 = np.arange(np.size(x1)) / float(np.size(x1))
        y2 = np.arange(np.size(x2)) / float(ma.count(x2))

        # ax1 = fig_1.add_subplot(1, 1, 1)
        ax5.plot(
            x1,
            y1,
        )
        ax5.plot(
            x2,
            y2,
        )
        ax5.set_title("CDF of the intensities")
        
        ax6 = fig_plot.add_subplot(spec[5])
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax6, show="intensity", x_axis="default", y_axis="intensity", limits=[0, 100])
           
        
        
    def plot_fitted(self, fig_plot=None, model=None, fit_centroid=None ):
        """
        add data to axes. 
        :param ax:
        :param show:
        :return:
        """
        
        #match max and min of colour scales
        limits={'max': np.max([np.max(self.intensity), np.max(model)]),
                'min': np.min([np.min(self.intensity), np.min(model)])}
        
        #plot data
        ax1 = fig_plot.add_subplot(1,3,1)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax1, show="intensity", x_axis="default", limits=limits, colourmap = "magma_r")
        ax1.set_title("Data")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        #plot model
        ax2 = fig_plot.add_subplot(1,3,2)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax2, data=model, limits=limits, colourmap = "magma_r")
        if fit_centroid is not None: 
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax2.set_title("Model")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
            
        #plot residuals
        ax3 = fig_plot.add_subplot(1,3,3)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax3, data=self.intensity-model, limits=[0, 100], colourmap = "residuals-blanaced")
        if fit_centroid is not None: 
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax3.set_title("Residuals")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        
        # tidy layout
        plt.tight_layout()
        
        
    @staticmethod
    def residuals_colour_scheme(maximum_value, minimum_value, **kwargs):
        # @staticmethod is needed or get an error saying:
        # "TypeError: residuals_colour_scheme() takes 2 positional arguments but 3 were given"
    
        # create custom colormap for residuals
        # ---------------
        #Need: a colour map that is white at 0 and the colours are equally scaled on each side. So it will match the intensity in black and white. Also one this is truncated so dont have lots of unused colour bar.
        #This can't be done using DivergingNorm(vcenter=0) or CenteredNorm(vcenter=0) so make new colourmap.
        #
        #create a colour map that truncates seismic so balanced around 0. 
        # It is not perfect because the 0 point insn't necessarily perfectly white but it is close enough (I think). 
        n_entries = 256
        all_colours = cm.seismic(np.arange(n_entries))
        
        if np.abs(maximum_value)> np.abs(minimum_value):
            n_cut = np.int(((2*maximum_value-(maximum_value-np.abs(minimum_value)))/(2*maximum_value))*n_entries)
            keep = n_entries-n_cut
            all_colours = all_colours[keep:]
        else:
            n_cut = np.int(((2*np.abs(minimum_value)-(maximum_value-np.abs(minimum_value)))/(2*np.abs(minimum_value)))*n_entries)
            keep = n_entries-n_cut
            all_colours = all_colours[:keep]
        all_colours = colors.ListedColormap(all_colours, name='myColorMap', N=all_colours.shape[0])
    
        return all_colours



    def plot_collected(self, fig_plot=None, axis_plot=None, show='intensity', limits=[0, 99.9] ):
        """
        add data to axes. 
        :param ax:
        :param show:
        :return:
        """
        
        IMax = np.percentile(self.intensity, limits[1])
        IMin = np.percentile(self.intensity, limits[0])
        #if IMin < 0:
        #    IMin = 0
            
        if limits[0] > 0 and limits[1] < 100: 
            cb_extend="both"
        elif limits[1] < 100:
            cb_extend="max"
        elif limits[0]>0:
            cb_extend="min"
        else:
            cb_extend="neither"        
                
        the_plot = axis_plot.imshow(self.intensity, vmin=IMin, vmax=IMax,
            cmap=plt.cm.jet)    
        axis_plot.set_xlabel("x")
        axis_plot.set_ylabel("y")
        axis_plot.invert_yaxis()
        
        fig_plot.colorbar(mappable=the_plot, extend=cb_extend)
        


    def plot_calibrated(self, fig_plot=None, axis_plot=None, show="intensity", x_axis="default", y_axis="default", data=None, limits=[0, 99.9], colourmap = "jet"):
        """
        add data to axes. 
        :param ax:
        :param show:
        :return:
        """
        rstr=False
        if self.intensity.size > 50000:
            print(
                " Have patience. The plot(s) will appear but it can take its time to render."
            )
            rstr = True
            
        if x_axis == "default":
            plot_x = self.tth
        else:
            plot_x = self.tth
        plot_y = self.azm  
                    
        if y_axis == "intensity":
            #plot y rather than azimuth on the y axis 
            plot_y = self.intensity
            #organise colour scale as azimuth
            plot_i = self.azm
            label_y = "Intensity (a.u.)"
            y_ticks = False
        else: #if y_axis is "default" or "azimuth"
            plot_y = self.azm
            plot_i = self.intensity
            label_y = "Azimuth (deg)"
            
            y_lims = np.array(
                [np.min(plot_y.flatten()), np.max(plot_y.flatten())]
            )
            y_lims = np.around(y_lims / 180) * 180
            axis_plot.set_ylim(y_lims)
            # y_ticks = list(range(int(y_lims[0]),int(y_lims[1]+1),45))

            y_ticks = self.dispersion_ticks()
            
        # organise the data to plot
        if data is not None:
            plot_i = data
        elif show == "unmasked_intensity":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_i = plot_i.data
        elif show == "mask":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_i = np.array(ma.getmaskarray(self.intensity), dtype="uint8") + 1
            colourmap="Greys"
        else: #if show == "intensity"
            plot_i = plot_i
            
        #set colour map
        if colourmap == "residuals-blanaced":
            colourmap = self.residuals_colour_scheme(np.max(plot_i.flatten()), np.min(plot_i.flatten())) 
        else:
            colourmap = colourmap
            
        #set colour bar and colour maps.   
        if colourmap=="Greys":
            IMax=2.01
            IMin=0
            cb_extend="neither"
        elif isinstance(limits, dict):
            IMax = limits['max']
            IMin = limits['min']
            cb_extend="neither"
        else:
            if limits[1] == 100:
                IMax = np.max(plot_i)
            else:
                IMax = np.percentile(plot_i, limits[1])
            if limits[0] == 0:
                IMin = np.min(plot_i)
            else:
                IMin = np.percentile(plot_i, limits[0])
            if limits[0] > 0 and limits[1] < 100: 
                cb_extend="both"
            elif limits[1] < 100:
                cb_extend="max"
            elif limits[0]>0:
                cb_extend="min"
            else:
                cb_extend="neither"
        
        #set axis limits
        x_lims = [np.min(plot_x.flatten()), np.max(plot_x.flatten())]
        axis_plot.set_xlim(x_lims)
        if y_ticks:
            axis_plot.set_yticks(y_ticks)
        
        the_plot = axis_plot.scatter(
            plot_x,
            plot_y,
            s=1,
            c=plot_i,
            edgecolors="none",
            cmap=colourmap,
            vmin=IMin,
            vmax=IMax,
            rasterized=rstr 
        )    
        axis_plot.set_xlabel(r"2$\theta$ (deg)")
        axis_plot.set_ylabel(label_y)
        
        
        fig_plot.colorbar(mappable=the_plot, extend=cb_extend)

