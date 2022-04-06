#!/usr/bin/env python
# -*- coding: utf-8 -*-

# __all__ = ['Requirements', 'ImportImage', 'DetectorCheck', 'Conversion', 'GetMask', 'GetCalibration', 'GetTth',
#            'GetAzm', 'GetDsp']

import sys
import os
from copy import deepcopy, copy
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from cpf import Med


# For MED files from Energy dispersive diffraction patterns at X17B2, 6-BMB and others.

# For ease of use the *.med files are imported using Mark Rivers' Med.py and Mca.py functions.
# These were downloaded from ??? in June 2019.
# and edited to remove unnecessary dependencies e.g. 'Numeric' and 'LinearAlgebra' packages and functions that had
# further dependencies.
# I don't know how much of this is required but for now I am just adding the packages.

# FIX ME: Need a dependency check.


class MedDetector:
    def __init__(self, calibration_parameters=None, fit_parameters=None):
        """
        :param calibration_parameters:
        """
        self.parameters = None
        if calibration_parameters is not None:
            self.parameters = self.get_calibration(calibration_parameters)
        
        self.requirements = self.get_requirements(fit_parameters)
        self.detector = None
        self.intensity = None
        self.tth = None
        self.azm = None
        self.dspace = None
        self.calibration = None

    def __deepcopy__(self,a):
        """
        Makes a "deep" copy of an Mca instance, using copy.deepcopy() on all of
        the attributes of the Mca instance. All of the attribute will point to
        new objects.
        """
        
        new = MedDetector()
        new.parameters =copy(self.parameters)
        new.calibration = copy(self.calibration)
        
        new.intensity = deepcopy(self.intensity)
        new.tth = deepcopy(self.tth)
        new.azm  = deepcopy(self.azm)
        new.dspace = deepcopy(self.dspace) 
        
        return new

    def fill_data(self, diff_file, settings=None, debug=None, calibration_mask=None):
        """
        :param diff_file:
        :param settings:
        :param debug:
        :param calibration_mask:
        :return:
        """
        self.calibration = self.get_calibration(settings.Calib_param)
        self.get_detector(diff_file, settings)

        self.intensity = self.get_masked_calibration(diff_file, debug, calibration_mask)
        self.tth = self.get_two_theta(mask=self.intensity.mask)
        self.azm = self.get_azimuth(mask=self.intensity.mask)
        self.dspace = self.get_d_space(mask=self.intensity.mask)

    def get_detector(self, calibration_data, detector=None):
        self.detector = self.detector_check(calibration_data, detector)

    @staticmethod
    def detector_check(image_name, detector=None):
        """
        :param image_name:
        :param detector:
        :return:
        """
        # FIX ME: we should try importing the image and reading the detector type from it.
        if detector is None:
            sys.exit(
                "\n\nDioptas requires a detector type.\nTo list the possible detector types in the command line type:"
                "\n   import pyFAI\n   pyFAI.detectors.Detector.registry\n\n"
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
        im = self.import_image(calibration_data, debug=debug)
        intens = ma.array(im)
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
            "Calib_detector"
            # 'Calib_data',        # Removed because this is not strictly required.
            "Calib_param",  # Now added to calibration function file.
            # 'Calib_pixels',      # this is only needed for GSAS-II and should be read from image file. FIX
            # ME: add to GSAS-II inputfile.
            # 'Calib_mask',        # a mask file is not strictly required.
            "datafile_directory",
            "datafile_Basename",
            "datafile_Ending",
            # 'datafile_StartNum',  # now optionally replaced by datafile_Files
            # 'datafile_EndNum',    # now optionally replaced by datafile_Files
            "datafile_NumDigit",
            # 'AziBins',            # required based on detector type
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


    def import_image(self,image_name, mask=None, debug=False):
        """
        :param image_name:
        :param debug:
        :return:
        """
        # FIX ME: why is mask in here? should it inherit from previous data if it exists?
        filename, file_extension = os.path.splitext(image_name)
        if file_extension == ".med":
            dat = Med.Med()
            dat.read_file(image_name)
            im_all = []
            for x in range(dat.n_detectors):
                im_all.append(dat.mcas[x].get_data())
            if mask is not None:
                im_all = ma.array(im_all, mask=mask)
            else:
                im_all = ma.array(im_all)
        else:
            sys.exit("Unknown file type. A *.med file is required")
            
        if mask is not None and ma.is_masked(self.intensity):
            self.intensity = ma.array(im_all, mask=self.intensity.mask)
        else:
            self.intensity = np.array(im_all)
                
        return im_all

    def conversion(self, e_in, calib, reverse=0, azm=None):
        """
        :param e_in:
        :param calib:
        :param reverse:
        :param azm:
        :return:
        """
        # convert energy into d-spacing.
        if azm is not None:
            e_in = np.array([[azm, e_in]])
            # FIX ME: this is a horrible way to make the conversion work, but it is needed for
            # processing the guesses of MED detectors.
        else:
            e_in = np.array([[0, e_in]])

        # print('E_in',E_in, 'calib,', calib)
        # print('type E_in', type(E_in))
        # print(E_in[:,1])
        # Convert Energy (keV) to wavelength (Angstroms)
        # E = hc/lambda;  h = 4.135667662(25)×10−15 eV s; c = 299792458 m/s
        wavelength = 4.135667662e-15 * (299792458 * 1e10) / (e_in[:, 1] * 1000)

        # determine which detector calibration goes with each energy.
        sort_idx = np.array(calib["azimuths"]).argsort()
        de = sort_idx[np.searchsorted(calib["azimuths"], e_in[:, 0], sorter=sort_idx)]

        if not reverse:
            # convert energy to d-spacing
            dspc_out = []
            for i in range(len(de)):
                dspc_out.append(
                    wavelength[i]
                    / 2
                    / np.sin(
                        np.radians(calib["calibs"].mcas[i].calibration.two_theta / 2)
                    )
                )
            dspc_out = np.squeeze(np.array(dspc_out))
        else:
            # convert d-spacing to energy.
            # N.B. this is the reverse function so that labels tth and d-spacing are not correct.
            dspc_out = []
            for i in range(len(de)):
                dspc_out.append(
                    wavelength[i]
                    / 2
                    / np.sin(
                        np.radians(calib["calibs"].mcas[i].calibration.two_theta / 2)
                    )
                )
                print("dspc_out", dspc_out)
            dspc_out = np.array(dspc_out)

        return dspc_out

    def get_mask(self, mask, im_ints=None, im_two_theta=None, im_azimuth=None, image_y=None, image_x=None, debug=False):
        """
        :param mask:
        :param im_ints:
        :param im_two_theta:
        :param im_azimuth:
        :param image_y:
        :param image_x:
        :param debug:
        :return:
        """
        # Masks for 10 element detector are either for full detector or for an energy range.
        # Currently, this just works by removing complete detectors.
        im_mask = ma.zeros(im_ints.shape)
        for x in range(len(mask)):
            im_mask[mask[x] - 1] = 1

        # FIX ME: DMF update to use class plot function!
        if debug:
            # Plot mask.
            # This is left in here for debugging.
            fig = plt.figure()
            ax = fig.add_subplot(1, 2, 1)
            plt.subplot(121)
            plt.scatter(im_two_theta, im_azimuth, s=4, c=im_ints, edgecolors="none", cmap=plt.cm.jet)
            ax.set_ylim([0, 360])
            ax = fig.add_subplot(1, 2, 2)
            plt.subplot(122)
            plt.scatter(
                ma.array(im_two_theta, mask=im_mask),
                ma.array(im_azimuth, mask=im_mask),
                s=4,
                c=im_ints,
                edgecolors="none",
                cmap=plt.cm.jet,
            )
            ax.set_ylim([0, 360])
            plt.colorbar()
            plt.show()
            plt.close()

        im_ints = ma.array(im_ints, mask=im_mask)
        return im_ints

    def get_calibration(self, file_name):
        """
        Parse inputs file, create type specific inputs.
        :param file_name:
        :return:
        """
        parms_dict = {}
        dat = Med.Med(file=file_name)
        parms_dict["DispersionType"] = "EnergyDispersive"
        parms_dict["calibs"] = dat
        if len(dat.mcas) == 10:
            az = np.linspace(0, 180, 9)
            az = np.append(az, 270)
            az[0] = az[0] + 2
            az[8] = az[8] - 2
            parms_dict["azimuths"] = az
        else:
            sys.exit(
                "\n\nDetector is unknown. Edit function to add number of detectors and associated azimuths.\n\n"
            )

        # 2 theta angle in degrees
        # parms_dict['conversion_constant'] = 6.5
        all_angle = []
        for x in range(parms_dict["calibs"].n_detectors):
            all_angle.append(parms_dict["calibs"].mcas[x].calibration.two_theta)
        parms_dict["conversion_constant"] = all_angle
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
        temp_azimuth = azimu.flatten()
        for i in range(len(bin_vals)):
            azi_chunk = np.where((temp_azimuth == bin_vals[i]))
            chunks.append(azi_chunk)
        return chunks, bin_vals

    # FIX ME: DMF missing function from dioptas get_azimuthal_integration -- is it needed here?

    def get_two_theta(self, mask=None, pix=None, det=None):
        """
        Called get two theta, but it is really getting the energy.
        :param mask:
        :param pix:
        :param det:
        :return:
        """
        channels = np.arange(self.calibration["calibs"].mcas[0].data.shape[0])
        en = []
        for x in range(self.calibration["calibs"].n_detectors):
            en.append(self.calibration["calibs"].mcas[x].channel_to_energy(channels))
        if mask is not None:
            en = ma.array(en, mask=mask)
        else:
            en = ma.array(en)
        return en

    def get_azimuth(self, mask=None, pix=None, det=None):
        """
        Give azimuth value for detector x,y position; have to assume geometry of the detectors
        :param mask:
        :param pix:
        :param det:
        :return:
        """
        l = self.calibration["calibs"].mcas[0].data.shape
        azi = np.tile(self.calibration["azimuths"], [l[0], 1]).T
        if mask is not None:
            azi = ma.array(azi, mask=mask)
        else:
            azi = ma.array(azi, mask=None)
        return azi

    def get_d_space(self, mask=None, pix=None, det=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        :param mask:
        :param pix:
        :param det:
        :return:
        """
        channels = np.arange(self.calibration["calibs"].mcas[0].data.shape[0])
        dspc = []
        for x in range(self.calibration["calibs"].n_detectors):
            dspc.append(self.calibration["calibs"].mcas[x].channel_to_d(channels))
        if mask is not None:
            d_space = ma.array(dspc, mask=mask)
        else:
            d_space = ma.array(dspc)
        return d_space


    def set_limits(self, range_bounds=[-np.inf, np.inf], i_max=np.inf, i_min=-np.inf):
        """
        Set limits to data 
        :param range_bounds:
        :param i_max:
        :param i_min:
        :return:
        """
                
        local_mask = np.array((self.tth <= range_bounds[0]) | (self.tth >= range_bounds[1]) | (self.intensity > i_max) | (self.intensity < i_min) )
        masks = (local_mask) + (self.intensity.mask)
   
        self.intensity = ma.masked_array(self.intensity,mask = masks)
        self.tth = ma.masked_array(self.tth,mask=self.intensity.mask)
        self.azm = ma.masked_array(self.azm,mask=self.intensity.mask)
        self.dspace = ma.masked_array(self.dspace,mask=self.intensity.mask)


    # @staticmethod
    # def residuals_colour_scheme(maximum_value, minimum_value, **kwargs):
    
    #     # create custom colormap for residuals
    #     # ---------------
    #     #Need: a colour map that is white at 0 and the colours are equally scaled on each side. So it will match the intensity in black and white. Also one this is truncated so dont have lots of unused colour bar.
    #     #This can't be done using DivergingNorm(vcenter=0) or CenteredNorm(vcenter=0) so make new colourmap.
    #     #
    #     #create a colour map that truncates seismic so balanced around 0. 
    #     # It is not perfect because the 0 point insn't necessarily perfectly white but it is close enough (I think). 
    #     n_entries = 256
    #     all_colours = cm.seismic(np.arange(n_entries))
        
    #     if np.abs(maximum_value)> np.abs(minimum_value):
    #         n_cut = np.int(((2*maximum_value-(maximum_value-np.abs(minimum_value)))/(2*maximum_value))*n_entries)
    #         keep = n_entries-n_cut
    #         all_colours = all_colours[keep:]
    #     else:
    #         n_cut = np.int(((2*np.abs(minimum_value)-(maximum_value-np.abs(minimum_value)))/(2*np.abs(minimum_value)))*n_entries)
    #         keep = n_entries-n_cut
    #         all_colours = all_colours[:keep]
    #     all_colours = colors.ListedColormap(all_colours, name='myColorMap', N=all_colours.shape[0])
    
    #     return all_colours
    
    

    def plot_masked(self, fig_plot=None):
        """
        Plot all the information needed to mask the data well. 
        :param fig:
        :return:
        """
        
        ax1 = fig_plot.add_subplot(111)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax1, x_axis="default", limits=[0, 100])
        
        
        
    def plot_fitted(self, fig_plot=None, model=None, fit_centroid=None ):
        """
        add data to axes. 
        :param ax:
        :param show:
        :return:
        """
        #reshape model to match self
        model2 = deepcopy(self.intensity)
        p = 0
        for i in range(len(self.intensity)):
            for j in range(len(self.intensity[0])):
                if not ma.is_masked(self.intensity[i][j]):
                    model2[i][j] = model[p]
                    p=p+1
        
        
        #match max and min of colour scales
        limits={'max': np.max([np.max(self.intensity), np.max(model)]),
                'min': np.min([np.min(self.intensity), np.min(model)])}
        
        #plot data
        ax1 = fig_plot.add_subplot(1,3,1)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax1, show="intensity", x_axis="default", limits=limits)
        ax1.set_title("Data")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        #plot model
        ax2 = fig_plot.add_subplot(1,3,2)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax2, data=model2, limits=limits)
        if fit_centroid is not None: 
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax2.set_title("Model")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
            
        #plot residuals
        ax3 = fig_plot.add_subplot(1,3,3)
        self.plot_calibrated(fig_plot=fig_plot, axis_plot=ax3, data=self.intensity-model2, limits=[0, 100])
        if fit_centroid is not None: 
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax3.set_title("Residuals")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        
        # tidy layout
        plt.tight_layout()

    def plot_collected(self, fig_plot=None, axis_plot=None, show='intensity', x_axis="default", limits=[0, 100]):
        """
        add data to axes in form collected in. 
        :param fig_plot:
        :param axis_plot:
        :param show:
        :param x_axis:
        :param limits:
        :return:
        """
        
        max_counts = np.max(self.intensity)
        spacing = max_counts/10
        
        for x in range(len(self.intensity)):
            axis_plot.plot(list(range(len(self.intensity[x]))), self.intensity[x] + spacing * x)  
        axis_plot.set_xlabel("Channel number")
        axis_plot.set_ylabel("Counts")


    def plot_calibrated(self, fig_plot=None, axis_plot=None, show='intensity', x_axis="default", y_axis="default", data=None, limits=[0, 99.9], colourmap = "jet"):
        """
        add data to axes in form collected in. 
        :param fig_plot:
        :param axis_plot:
        :param show:
        :param x_axis:
        :param limits:
        :return:
        """
        
        max_counts = np.max(self.intensity)
        spacing = max_counts/10/20
        
        if x_axis == "default":
            plot_x = self.tth
        else:
            plot_x = self.tth
        label_x = "Energy (keV)"
        
        
        if data is not None:
            plot_i = data
        else:
            plot_i = self.intensity
        
        
        if y_axis == "intensity":
            #y is unmodified intensity; colour of data by azimuth
            plot_y = plot_i
            plot_c = self.azm
            label_y = "Intensity (a.u.)"
        elif y_axis == "azimuth":
            #y is unmodified azimuth, intensity is size of dots amd colour scale 
            plot_y = self.azm
            plot_c = plot_i/np.max(self.intensity)*200
            plot_s = plot_i/np.max(self.intensity)*200
            label_y = "Azimuth (deg)"
        else: #if y_axis is default
            #y is intensity distributed by azimuth; colour is azimuth 
            plot_y = plot_i + self.azm * spacing# + dispersion by azimuth
            plot_c = self.azm
            label_y = "Counts"
            
        # sort the data
        # here assume x isenergy, y is azzimuth, z is intensity or counts.
        if show == "unmasked_intensity":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_c = plot_c
        elif show == "mask":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_c = np.array(ma.getmaskarray(self.intensity), dtype="uint8") + 1
            colourmap="Greys"
        else: 
            plot_c = plot_c
        
        #set colour map
        if colourmap == "residuals-blanaced":
            colourmap = self.residuals_colour_scheme(np.max(plot_c.flatten()), np.min(plot_c.flatten())) 
        else:
            colourmap = colourmap
        

        if 'plot_s' in locals():
                    
            the_plot = axis_plot.scatter(plot_x, plot_y, s=plot_s, c=plot_c, cmap=colourmap, alpha=.2)
            fig_plot.colorbar(mappable=the_plot).set_label(label="Counts")
            
        else:
                
            # Set colour map range - to be in 90 degree blocks
            blocks = 90
            colour_start = np.floor((np.min(self.azm)-.1)/blocks)*blocks
            colour_end   = np.ceil((np.max(self.azm))/blocks)*blocks
            normalize = colors.Normalize(vmin=colour_start, vmax=colour_end)
            c_map = cm.get_cmap(name=colourmap)
            
            for x in range(len(plot_y)-1,-1,-1):
                 if not plot_c.mask.all():
                     colour = c_map(normalize(np.mean(plot_c[x])))
                     the_plot=axis_plot.plot(plot_x[x], plot_y[x], color=colour)
            
            # Colorbar
            s_map = cm.ScalarMappable(norm=normalize, cmap=c_map)
            # label colour bar with azimuths if there are less than 10
            if len(self.intensity) <= 10:
                ticks = np.unique(plot_c)
            else:
                ticks = None
            cbar = fig_plot.colorbar(s_map, ticks=ticks)
            cbar.set_label("Azimuth")
            
        axis_plot.set_xlabel(label_x)
        axis_plot.set_ylabel(label_y)