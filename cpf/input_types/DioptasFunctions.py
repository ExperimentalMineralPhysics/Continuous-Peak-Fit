#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["DioptasDetector"]


import sys
from copy import deepcopy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import fabio
import pyFAI
from pyFAI.io import ponifile
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
from PIL import Image
from cpf.input_types._Plot_AngleDispersive import _Plot_AngleDispersive
from cpf import IO_functions
import cpf.h5_functions as h5_functions



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
        
        self.azm_blocks = 45
        
        # Single image plate detector so contonuous data.
        self.continuous_azm = True
        self.DispersionType = "AngleDispersive"

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



    def fill_data(
        self, diff_file, settings=None, debug=None
    ):  # , calibration_mask=None):
        """
        :param settings: -- now a class.
        :param diff_file:
        :param debug:
        :param calibration_mask: -- now removed because calibration mask file is contained in settings class
        """

        # self.detector = self.detector_check(calibration_data, detector)
        #        print(settings.Calib_param)
        self.get_calibration(settings.calibration_parameters)
        self.get_detector(diff_file, settings)

        self.intensity = self.get_masked_calibration(
            diff_file, debug, calibration_mask=settings.calibration_mask
        )
        self.tth = self.get_two_theta(mask=self.intensity.mask)
        self.azm = self.get_azimuth(mask=self.intensity.mask)
        self.dspace = self.get_d_space(mask=self.intensity.mask)
        
        blocks = 180
        self.azm_start = np.around(np.min(self.azm.flatten()) / blocks) * blocks
        self.azm_end = np.around(np.max(self.azm.flatten()) / blocks) * blocks

    def get_detector(self, diff_file=None, settings=None):
        """
        :param diff_file:
        :param settings:
        :return:
        """
        # if self.calibration.detector:
        #     #have read the poni calibration and the detector is present.
        #     #this seems to be the updated Dioptas/pyFAI way of storing the pixel size [i.e. poni version 2]???
        #     self.detector = self.calibration.detector
        # el
        if settings.calibration_detector is not None:
            #if we tell it the type of detector.
            self.detector = pyFAI.detector_factory(settings.calibration_detector)
        else:
            #we dont know anything. so recreate a detector. 
                
            #open the file to get the shape of the data.     
            if diff_file is not None: 
                im_all = fabio.open(diff_file)
            elif settings.calibration_data is not None:
                im_all = fabio.open(settings.calibration_data)
            else:
                im_all = fabio.open(settings.image_list[0])
                

            if self.calibration.detector:
                 #have read the poni calibration and the detector is present.
                 #this seems to be the updated Dioptas/pyFAI way of storing the pixel size [i.e. poni version 2]???
                 #self.detector = self.calibration.detector
                 sz1 = self.calibration.detector.pixel1
                 sz2 = self.calibration.detector.pixel2
                 
            elif "pixelsize1" in self.calibration:
                #old Dioptas/pyFAI way of storing the pixel size???
                sz1 = self.calibration["pixelsize1"]
                sz2 = self.calibration["pixelsize2"]
            elif settings.calibration_pixel_size is not None:
                #get the pixel size from the input file
                sz1 = settings.calibration_pixel_size * 1e-6
                sz2 = settings.calibration_pixel_size * 1e-6
            else:
                #this won't work
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
            # sz = calibration_data.Calib_pixels  # Pixel_size
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

    # @staticmethod
    def import_image(self, image_name, mask=None, debug=False):
        """
        Import the data image
        :param mask:
        :param image_name: Name of file
        :param debug: extra output for debugging - plot
        :return: intensity array
        """
        # FIX ME: why is mask in here? should it inherit from previous data if it exists?

        # FIX ME: nxs files need to be checked. But they should be read like h5 files.

        # im = Image.open(ImageName) ##always tiff?- no
        # filename, file_extension = os.path.splitext(image_name)
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

        # FIX ME: Copied from MED_functions.
        # FIX ME: Convert the input data from integer to float because otherwise the model inherits Integer properties somewhere along the line.
        # FIX ME: I dont know if this needed here but given everything else is the same I have copied it from MED_functions.
        im = np.array(im, dtype="f")

        # Dioptas flips the images to match the orientations in Plot2D (should this be Fit2D?).
        # Therefore implemented here to be consistent with Dioptas.
        # FIX ME: what to do about this?
        im = np.array(im)[::-1]

        if debug:
            print("min+max:", np.min(im), np.max(im))
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            self.plot_collected(fig_plot=fig, axis_plot=ax)
            # plt.title(os.path.split(image_name)[1])
            plt.title(IO_functions.title_file_names(image_name=image_name))
            plt.show()
            plt.close()

        if mask == None and ma.is_masked(self.intensity) == False:
            self.intensity = ma.array(im)
            return ma.array(im)
        elif mask is not None:
            self.intensity = ma.array(im, mask=mask)
            return ma.array(im, mask=mask)
        else:
            self.intensity = ma.array(im, mask=self.intensity.mask)
            return ma.array(im)

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

        wavelength = self.calibration.wavelength * 1e10
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
        # im_mask = np.array(im_mask)[::-1]
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



    def get_calibration(self, file_name=None, settings=None):
        """
        Opens the file with the calibration data in it and updates 
        DioptasCalibration.calibration and 
        DioptasCalibration.conversion_constant
        
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
        self.conversion_constant = pf.wavelength * 1e10



    def bins(self, orders_class, cascade=False):
        """
        Assign data into azimuthal bins.
        This is primarily called as part of the initial fitting.
        Assign each data to a chunk corresponding to its azimuth value
        Returns array with indices for each bin and array of bin centroids
        
        
        Parameters
        ----------
        orders_class : settings class, optional
            cpf settings class.
        cascade : True/False, optional
            Switch that determines type of outputs returned.
            The default is False.

        Returns
        -------
        chunks : list
            list of arrays. Each element of the array is an array of data indicies that
            are part of the bin
        bin_mean_azi : array
            list containing mean azimuth of data in each bin.
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

        # assign the data to the bins
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


    def get_azimuthal_integration(self):
        """

        :return:
        """

        ai = AzimuthalIntegrator(
            dist=self.calibration.dist,
            poni1=self.calibration.poni1,
            poni2=self.calibration.poni2,
            rot1=self.calibration.rot1,
            rot2=self.calibration.rot2,
            rot3=self.calibration.rot3,
            detector=self.detector,
            wavelength=self.calibration.wavelength,
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
                self.calibration.wavelength * 1e10 / 2 / np.sin(ai.twoThetaArray() / 2),
                mask=mask,
            )
        else:
            return ma.array(
                self.calibration.wavelength * 1e10 / 2 / np.sin(ai.twoThetaArray() / 2)
            )

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
DioptasDetector.plot_masked      = _Plot_AngleDispersive.plot_masked
DioptasDetector.plot_fitted      = _Plot_AngleDispersive.plot_fitted
DioptasDetector.plot_collected   = _Plot_AngleDispersive.plot_collected
DioptasDetector.plot_calibrated  = _Plot_AngleDispersive.plot_calibrated
#this function is added because it requires access to self:
DioptasDetector.dispersion_ticks = _Plot_AngleDispersive._dispersion_ticks

