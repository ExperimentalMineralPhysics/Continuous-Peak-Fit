#!/usr/bin/env python
# -*- coding: utf-8 -*-


# XYFunctions.
#
# For simple x, y data the does not require any complex conversions. 
# The peaks are asumed to run vertically in the image (Y direction)
# If the peaks run the other way then use the YXFunctions. 
#
# The 'Calibration', such as it is, are simple finctions of x or y. 
#
# Functions to change for another imput type:
#   -   get_calibration
#   -   get_detector 
#   -   Load/import data
#   -   Load mask
#   -   Load Calibration



__all__ = ["XYDetector"]

import os
import sys
from copy import deepcopy
from matplotlib import image
#import h5py
import cpf.h5_functions as h5_functions
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
#import pyFAI
#from PIL import Image
from matplotlib import gridspec
#from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
#from pyFAI.io import ponifile
from matplotlib import cm, colors
from cpf import IO_functions
import pickle
import cpf.lmfit_model as llm

class XYDetector:

    def __init__(self, settings_class=None):
        """
        :param calibration_parameters:
        """

        self.requirements = self.get_requirements()

        self.calibration = None
        self.detector = None

        self.x = None
        self.y = None
        self.azm_start = None
        self.azm_end   = None
        self.tth_start = None
        self.tth_end   = None
        self.intensity = None
        self.tth = None
        self.azm = None
        self.dspace = None
        self.continuous_azm = True
        
        #set the start and end limits for the data
        self.start = 0
        self.end   = np.inf

        if settings_class:
            self.get_calibration(settings=settings_class)
        if self.calibration:
            self.detector = self.get_detector(settings=settings_class)


    def duplicate(self):
        """
        Makes a "deep" copy of an XY Instance, using copy.deepcopy()
        :return: deepcopy of self.
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
        self.get_calibration(settings=settings)
        self.get_detector(diff_file, settings)
        
        self.intensity = self.get_masked_calibration(
            diff_file, debug, calibration_mask=settings.calibration_mask
        )
        self.tth = self.get_two_theta(mask=self.intensity.mask)
        self.azm = self.get_azimuth(mask=self.intensity.mask)
        self.dspace = self.get_d_space(mask=self.intensity.mask)
        if self.azm_start is None:
            self.azm_start = np.min(self.azm)
        if self.azm_end is None:
            self.azm_end   = np.max(self.azm)

        
    def get_detector(self, diff_file=None, settings=None):
        """
        :param diff_file:
        :param settings:
        """

        # the detector has nothing to contain for a single x, y image with no conversions. 
        # but it might be referenced by other bits and bobs so has to exist as something. 
        self.detector = None


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
        
        intens = ma.array(im.data)
        # create mask from mask file if present. If not make all values valid
        if calibration_mask:
            intens = self.get_mask(calibration_mask, intens)
        else:
            im_mask = np.zeros_like(intens)
            intens = ma.array(intens, mask=im_mask)
        return intens


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
        if isinstance(image_name, list):
            # then it is a h5 type file
            im = h5_functions.get_images(image_name)
        elif os.path.splitext(image_name)[1] == ".txt":
            # then it is a text file, assume there are less than 20 header rows. 
            done = 0
            n = 0
            while done == 0:
                try:
                    im = np.loadtxt(image_name, skiprows=n)
                    done = 1
                except:
                    done = 0
                    n += 1
                    if n>20:
                        # issue an error
                        err_str = "There seems to be more than 20 header rows in the data file. Is this correct?"
                        raise ValueError(err_str)
        else:
            im = image.imread(image_name)
            im = im.data

        #self.intensity = ma.array(im, dtype="f")
        
        if 0:#debug:
            print(self.intensity)
            print("min+max:", np.min(self.intensity), np.max(self.intensity))
            print("min+max:", np.nanmin(self.intensity), np.nanmax(self.intensity))
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
        Convert between two theta and d-spacing. 
        This is not used for XY detector class.
        So Value in = value out
        :param tth_in:
        :param reverse:
        :param azm:
        :return:
        """

        return tth_in

    def get_mask(self, msk_file, im_ints):
        """
        Dioptas mask is compressed Tiff image.
        Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py
        :param msk_file:
        :param im_ints:
        :return:
        """
        from PIL import Image, ImageDraw
        if isinstance(msk_file, dict):
                        
            im_mask = np.zeros(self.intensity.shape, 'bool')
            
            polygons = msk_file['polygon']
            for i in polygons:
                img = Image.new('L', self.intensity.shape, 0)
                ImageDraw.Draw(img).polygon(i, outline=1, fill=1)
                im_mask = im_mask +ma.array(img)
            threshold = msk_file['threshold']
            im_mask = np.asarray(im_mask) | ma.masked_outside(self.intensity,int(threshold[0]),threshold[1]).mask
            #mask invalid values
            mask2 = ma.masked_invalid(self.intensity).mask
            #combine masks
            im_mask = np.asarray(im_mask) | np.asarray(mask2)
            im_ints = ma.array(im_ints, mask=im_mask)
            plt.imshow(im_mask)
            
        else:
            #the mask is an image.
            im_mask = np.array(Image.open(msk_file))
            im_ints = ma.array(im_ints, mask=im_mask)    
            im_ints = ma.masked_less(im_ints, 0)
        
        self.original_mask = im_mask
        return im_ints

    def f(x=None, m=0, c=1):
        print(type(x), type(m), type(c))
        return (x * m) + c
    
    def get_calibration(self, file_name=None, settings=None):
        """
        Parse inputs file, create type specific inputs
        :rtype: object
        :param file_nam:
        :return:
        """
        #make empty calibration
        self.calibration = {
                    "x_dim":   0, 
                    "x":       [0,1], 
                    "x_start": None, # begining of the 'tth' axis
                    "x_end":   None, # end of the 'tth' axis
                    "y":       [0,1], 
                    "y_start": None, # begining of the 'azimith' axis
                    "y_end":   None, # end of the 'azimith' axis
                    "conversion_constant": 1
                    }

        # parms_file is file
        if settings != None:
            calib_temp = settings.calibration_parameters
        else:
            parms_file = file_name
            with open(parms_file, 'wb') as f:
                calib_temp = pickle.dump("calibration", f)

        #pass calib_temp into the calibration
        for key, val in calib_temp.items():
            self.calibration[key] = val
            
        # if parms_file != None and isinstance(parms_file, dict):
        #     # here parms_file is not a file name but a dictionary of parameters.
        #     calib_temp = parms_file
        # else:
        #     # parms_file is file
        #     if settings != None:
        #         parms_file = settings.calibration_parameters
        #     else:
        #         parms_file = file_name
                
        #     with open(parms_file, 'wb') as out_file:
        #         calib_temp = pickle.dump(parameters, out_file)


        # Currently no conversion from X,Y to aximuth and two theta (or whatever the units are)
        #pf = ponifile.PoniFile()
        #pf.read_from_file(parms_file)

        #self.calibration = pf
        #self.conversion_constant = pf.wavelength * 1e10

        # if parms_file != None and isinstance(parms_file, dict):
        #     # here parms_file is not a file name but a dictionary of parameters.
        #     calib_temp = parms_file
        #     if "x_dim" in calib_temp:
        #         self
            
            
            
        # elif parms_file != None and parms_file != "":
            
        #     with open(parms_file, 'wb') as out_file:
        #         self.calibration = pickle.dump(parameters, out_file)
        # else:
        #     # default empty calibration
        #     self.calibration = {"x_dim": 0, "x": [0,1], "y":[0,1], "y_start":0, "y_end":1}
        
        self.conversion_constant = self.calibration["conversion_constant"]
        self.azm_start = self.calibration["y_start"]
        self.azm_end = self.calibration["y_end"]
        

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
            if orders_class.azi_limits == "tight":
                pass
            elif orders_class.azi_limits == "full":
                lims = [0, 360]
            else: #orders_class.azi_bounds == "default":
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

    def test_azims(self, steps = 1000):
        """
        Returns equally spaced set of aximuths within possible range.

        Parameters
        ----------
        steps : Int, optional
            DESCRIPTION. The default is 1000.

        Returns
        -------
        array
            list of possible azimuths.

        """
        return np.linspace(self.azm_start,self.azm_end, steps+1)
    
    def get_xy(self, mask=None):
        """
        Gets the twotheta (i.e. calibrated) values of the data. 
        :param mask:
        :return:
        """

        nx, ny = self.intensity.shape
        x = np.linspace(0, nx-1, nx)
        y = np.linspace(0, ny-1, ny)
        xv, yv = np.meshgrid(x, y)

        if self.calibration["x_dim"]==0:
            return xv, yv
        else:
            return ma.array(yv) * self.calibration["x"][1] + self.calibration["x"][0]
        

    def get_two_theta(self, mask=None):
        """
        Gets the twotheta (i.e. calibrated) values of the data. 
        :param mask:
        :return:
        """

        nx, ny = self.intensity.shape
        x = np.linspace(0, nx-1, nx)
        y = np.linspace(0, ny-1, ny)
        xv, yv = np.meshgrid(x, y)

        if self.calibration["x_dim"]==0:
            return ma.array(xv) * self.calibration["x"][1] + self.calibration["x"][0]
        else:
            return ma.array(yv) * self.calibration["x"][1] + self.calibration["x"][0]


    def get_azimuth(self, mask=None):
        """
        Give azimuth value for detector x,y position; calibration info in data
        :param mask:
        :return:
        """
        
        nx, ny = self.intensity.shape
        x = np.linspace(0, nx-1, nx)
        y = np.linspace(0, ny-1, ny)
        xv, yv = np.meshgrid(x, y)

        if self.calibration["x_dim"]==0:
            return ma.array(yv) * self.calibration["y"][1] + self.calibration["y"][0]
        else:
            return ma.array(xv) * self.calibration["y"][1] + self.calibration["y"][0]



    def get_d_space(self, mask=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        Get as twotheta and then convert into d-spacing
        :param mask:
        :return:
        """
        
        return self.get_two_theta()

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

    def dispersion_ticks(self, disp_ticks=None, unique=10):
        """
        Returns the labels for the dispersion axis/colour bars.

        :param disp_ticks: -- unused for maintained for compatibility with MED functions
        :param unique:
        :return new_tick_positions:
        """
        if len(np.unique(self.azm)) >= unique:
            disp_lims = np.array(
                [self.azm_start, self.azm_end]
            )
            if disp_lims[1]-disp_lims[0] == 360:
                disp_lims = np.around(disp_lims / 180) * 180
                disp_ticks = list(range(int(disp_lims[0]), int(disp_lims[1] + 1), 45))

        return disp_ticks

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
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax1,
            show="unmasked_intensity",
            x_axis="default",
            limits=[0, 100],
        )
        ax1.set_title("All Data")
        ax2 = fig_plot.add_subplot(spec[1])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax2,
            show="mask",
            x_axis="default",
            limits=[0, 100],
        )
        ax2.set_title("Mask")
        ax3 = fig_plot.add_subplot(spec[2])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax3,
            show="intensity",
            x_axis="default",
            limits=[0, 100],
        )
        ax3.set_title("Masked Data")

        ax4 = fig_plot.add_subplot(spec[3])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax4,
            show="unmasked_intensity",
            x_axis="default",
            y_axis="intensity",
            limits=[0, 100],
        )

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
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax6,
            show="intensity",
            x_axis="default",
            y_axis="intensity",
            limits=[0, 100],
        )

    def plot_fitted(self, fig_plot=None, model=None, fit_centroid=None):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """

        # match max and min of colour scales
        limits = {
            "max": np.max([np.max(self.intensity), np.max(model)]),
            "min": np.min([np.min(self.intensity), np.min(model)]),
        }

        # plot data
        ax1 = fig_plot.add_subplot(1, 3, 1)
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax1,
            show="intensity",
            x_axis="default",
            limits=limits,
            colourmap="magma_r",
        )
        ax1.set_title("Data")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        # plot model
        ax2 = fig_plot.add_subplot(1, 3, 2)
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax2,
            data=model,
            limits=limits,
            colourmap="magma_r",
        )
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax2.set_title("Model")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)

        # plot residuals
        ax3 = fig_plot.add_subplot(1, 3, 3)
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax3,
            data=self.intensity - model,
            limits=[0, 100],
            colourmap="residuals-blanaced",
        )
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
        # Need: a colour map that is white at 0 and the colours are equally scaled on each side. So it will match the intensity in black and white. Also one this is truncated so dont have lots of unused colour bar.
        # This can't be done using DivergingNorm(vcenter=0) or CenteredNorm(vcenter=0) so make new colourmap.
        #
        # create a colour map that truncates seismic so balanced around 0.
        # It is not perfect because the 0 point insn't necessarily perfectly white but it is close enough (I think).
        n_entries = 256
        all_colours = cm.seismic(np.arange(n_entries))

        if np.abs(maximum_value) > np.abs(minimum_value):
            n_cut = np.int(
                (
                    (2 * maximum_value - (maximum_value - np.abs(minimum_value)))
                    / (2 * maximum_value)
                )
                * n_entries
            )
            keep = n_entries - n_cut
            all_colours = all_colours[keep:]
        else:
            n_cut = np.int(
                (
                    (
                        2 * np.abs(minimum_value)
                        - (maximum_value - np.abs(minimum_value))
                    )
                    / (2 * np.abs(minimum_value))
                )
                * n_entries
            )
            keep = n_entries - n_cut
            all_colours = all_colours[:keep]
        all_colours = colors.ListedColormap(
            all_colours, name="myColorMap", N=all_colours.shape[0]
        )

        return all_colours

    def plot_collected(
        self, fig_plot=None, axis_plot=None, show="intensity", limits=[0, 99.9]
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """
        IMax = np.nanpercentile(self.intensity.compressed(), limits[1])
        IMin = np.nanpercentile(self.intensity.compressed(), limits[0])
        
        if limits[0] > 0 and limits[1] < 100:
            cb_extend = "both"
        elif limits[1] < 100:
            cb_extend = "max"
        elif limits[0] > 0:
            cb_extend = "min"
        else:
            cb_extend = "neither"

        the_plot = axis_plot.imshow(
            self.intensity, vmin=IMin, vmax=IMax, cmap=plt.cm.jet
        )
        axis_plot.set_xlabel("x")
        axis_plot.set_ylabel("y")
        axis_plot.invert_yaxis()

        fig_plot.colorbar(mappable=the_plot, extend=cb_extend)

    def plot_calibrated(
        self,
        fig_plot=None,
        axis_plot=None,
        show="intensity",
        x_axis="default",
        y_axis="default",
        data=None,
        limits=[0, 99.9],
        colourmap="jet",
        rastered=False,
        point_scale=3,
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """

        if self.intensity.size > 50000:
            print(
                " Have patience. The plot(s) will appear but it can take its time to render."
            )
            rastered = True
            # print(type(axis_plot))
            # axis_plot = raster_axes.RasterAxes(axes=axis_plot)
            # print(type(axis_plot))

        if x_axis == "default":
            plot_x = self.tth
        else:
            plot_x = self.tth
        plot_y = self.azm

        if y_axis == "intensity":
            # plot y rather than azimuth on the y axis
            plot_y = self.intensity
            # organise colour scale as azimuth
            plot_i = self.azm
            label_y = "Intensity (a.u.)"
            y_ticks = False
        else:  # if y_axis is "default" or "azimuth"
            plot_y = self.azm
            plot_i = self.intensity
            label_y = self.calibration["y_label"]

            #y_lims = np.array([np.min(plot_y.flatten()), np.max(plot_y.flatten())])
            #y_lims = np.around(y_lims / 180) * 180
            y_lims = [self.azm_start, self.azm_end]
            axis_plot.set_ylim(y_lims)
            # y_ticks = list(range(int(y_lims[0]),int(y_lims[1]+1),45))
            
            y_ticks = False
            if y_lims[1]-y_lims[0] == 360:
                y_ticks = list(range(int(y_lims[0]), int(y_lims[1] + 1), 45))


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
            colourmap = "Greys"
        else:  # if show == "intensity"
            plot_i = plot_i

        # set colour map
        if colourmap == "residuals-blanaced":
            colourmap = self.residuals_colour_scheme(
                np.max(plot_i.flatten()), np.min(plot_i.flatten())
            )
        else:
            colourmap = colourmap

        # set colour bar and colour maps.
        if colourmap == "Greys":
            IMax = 2.01
            IMin = 0
            cb_extend = "neither"
        elif isinstance(limits, dict):
            IMax = limits["max"]
            IMin = limits["min"]
            cb_extend = "neither"
        else:
            if limits[1] == 100:
                IMax = np.nanmax(plot_i)
            else:
                IMax = np.nanpercentile(plot_i.compressed(), limits[1])
            if limits[0] == 0:
                IMin = np.nanmin(plot_i)
            else:
                IMin = np.nanpercentile(plot_i.compressed(), limits[0])
            if limits[0] > 0 and limits[1] < 100:
                cb_extend = "both"
            elif limits[1] < 100:
                cb_extend = "max"
            elif limits[0] > 0:
                cb_extend = "min"
            else:
                cb_extend = "neither"

        # set axis limits
        x_lims = [np.min(plot_x.flatten()), np.max(plot_x.flatten())]
        axis_plot.set_xlim(x_lims)
        if y_ticks:
            axis_plot.set_yticks(y_ticks)

        # the_plot = axis_plot.imshow(plot_i, 
        #                             cmap=colourmap
        #                             vmin=IMin,
        #                             vmax=IMax,)

        the_plot = axis_plot.scatter(
            plot_x,
            plot_y,
            s=1,
            c=plot_i,
            edgecolors="none",
            cmap=colourmap,
            vmin=IMin,
            vmax=IMax,
            rasterized=rastered,
        )
        axis_plot.set_xlabel(self.calibration["x_label"])
        axis_plot.set_ylabel(label_y)

        fig_plot.colorbar(mappable=the_plot, extend=cb_extend)
