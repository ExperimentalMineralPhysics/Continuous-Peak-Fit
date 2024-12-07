#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["MedDetector"]


import os
import re
import sys
from copy import copy, deepcopy

import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
from matplotlib import cm, colors, gridspec

from cpf.input_types import Med, med_detectors
from cpf.input_types._AngleDispersive_common import _AngleDispersive_common
from cpf.input_types._Masks import _masks
from cpf.logger_functions import CPFLogger

logger = CPFLogger("cpf.input_types.MedFunctions")


class MedDetector:
    """
    Data class for Energy dispersive diffraction.

    The defines the class for energy dispersive diffraction detectors and is
    currently implemented for the *.med files used by GSECARS and the APS,
    specifically APS beamline 6-BMB.

    The data reading uses Mark Rivers' Med.py and Mca.py functions, which were
    downloaded from Github in June 2019. These were then edited to remove
    unnecessary dependencies e.g. the 'Numeric' and 'LinearAlgebra' packages
    and other functions that had further dependencies.
    """

    def __init__(self, settings_class=None):
        """
        :param calibration_parameters:
        """

        self.requirements = self.get_requirements()

        self.calibration = None
        self.detector = None

        self.plot_orientation = "horizontal"

        self.intensity = None
        self.tth = None
        self.azm = None

        self.dspace = None
        self.x = None
        self.y = None
        self.azm_start = None
        self.azm_end = None
        self.tth_start = None
        self.tth_end = None
        self.DispersionType = "EnergyDispersive"
        # separate detectors around the ring so not continuous
        self.continuous_azm = False

        self.azm_blocks = 45

        self.calibration = None
        self.conversion_constant = None
        self.detector = None

        if settings_class:
            self.calibration = self.get_calibration(settings=settings_class)
        if self.calibration:
            self.detector = self.get_detector(settings=settings_class)

    def duplicate(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf]):
        """
        Makes an independent copy of a MedDetector Instance.

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
        new : MedDetector Instance.
            Copy of MedDetector with independedent data values

        """

        new = copy(self)

        local_mask = np.where(
            (self.tth >= range_bounds[0])
            & (self.tth <= range_bounds[1])
            & (self.azm >= azi_bounds[0])
            & (self.azm <= azi_bounds[1])
        )

        new.intensity = deepcopy(self.intensity[local_mask])
        new.tth = deepcopy(self.tth[local_mask])
        new.azm = deepcopy(self.azm[local_mask])

        new.intensity = new.intensity[new.intensity.mask == False]
        new.tth = new.tth[new.tth.mask == False]
        new.azm = new.azm[new.azm.mask == False]

        if "dspace" in dir(new):
            if self.dspace is not None:
                new.dspace = deepcopy(self.dspace[local_mask])
                new.dspace = new.dspace[new.dspace.mask == False]

        if "x" in dir(new):
            if self.x is not None:
                new.x = deepcopy(self.x.squeeze()[local_mask])
                new.x = new.x[new.x.mask == False]
        if "y" in dir(new):
            if self.y is not None:
                new.y = deepcopy(self.y.squeeze()[local_mask])
                new.y = new.y[new.y.mask == False]
        if "z" in dir(new):
            if self.z is not None:
                new.z = deepcopy(self.z.squeeze()[local_mask])
                new.z = new.z[new.z.mask == False]

        return new

    def get_calibration(self, file_name=None, settings=None, debug=False):
        """
        Opens the *.med file with the calibration data in it and updates
        MedClass.calibration and MedClass.conversion_constant

        Either file_name or settings should be set.

        Parameters
        ----------
        file_name : string, optional
            Filename of the calibration file. In this case a *.med file.
            The default is None.
        settings : settings class, optional
            cpf settings class containing the calibration file name.
            The default is None.

        Returns
        -------
        None.

        """
        parms_dict = {}

        if settings != None:
            dat = Med.Med(file=settings.calibration_parameters)
        else:
            dat = Med.Med(file=file_name)

        # save this block of code incase it is needed for future studies/datasets
        # this is a blank calibratiom from a dead detector. Make something more reasonable.
        # if i in range(len(dat)):
        #     if (#dat.mcas[i].calibration.two_theta == 6.5
        #         and dat.mcas[i].calibration.offset = 0.0
        #         and dat.mcas[i].calibration.slope = 0.0
        #         and dat.mcas[i].calibration.quad = 0.0
        #         )
        #     dat.mcas[i].calibration.two_theta == 6.5
        #     dat.mcas[i].calibration.offset = 0.1
        #     dat.mcas[i].calibration.slope = 0.06
        #     dat.mcas[i].calibration.quad = 4E-8

        if settings.calibration_detector == "10-element":
            det_name = "W2010_10element"
        else:
            det_name = settings.calibration_detector

        # get the azimuthal postions of the detectors.
        if hasattr(med_detectors, det_name):
            # the detector is in the default list.
            # so read it.
            parms_dict["azimuths"] = getattr(med_detectors, det_name)()
        else:
            # the detector positions are in a file.
            a, ext = os.path.splitext(det_name)
            if len(ext) == 0:
                det_name += ".det"

            try:
                with open(det_name, "r") as f:
                    data = []
                    for i in range(
                        50
                    ):  # read only 50 rows here. more than there are detectors
                        for line in f:
                            if re.match("^#", line):
                                data.append(line)
                start_col = max(enumerate(data))[0]
                det = pd.read_csv(
                    det_name, sep=",", skiprows=start_col
                ).values.tolist()  # use your actual delimiter.

                det = [item for sublist in det for item in sublist]

                parms_dict["azimuths"] = det
            except:
                sys.exit("\n\n Energy dispersive detector is unknown.\n\n")

        # get and store the calibration parameters
        two_theta = []
        offset = []
        slope = []
        quad = []
        for x in range(dat.n_detectors):
            two_theta.append(dat.mcas[x].calibration.two_theta)
            offset.append(dat.mcas[x].calibration.offset)
            slope.append(dat.mcas[x].calibration.slope)
            quad.append(dat.mcas[x].calibration.quad)
        parms_dict["conversion_constant"] = two_theta
        parms_dict["two_theta"] = two_theta
        parms_dict["offset"] = offset
        parms_dict["slope"] = slope
        parms_dict["quad"] = quad

        parms_dict["mcas"] = dat

        self.azm_start = (
            np.floor((np.min(np.array(parms_dict["azimuths"])) / self.azm_blocks))
            * self.azm_blocks
        )
        self.azm_end = (
            np.ceil((np.max(np.array(parms_dict["azimuths"])) / self.azm_blocks))
            * self.azm_blocks
        )

        self.calibration = parms_dict

    def get_detector(
        self, settings=None, calibration_file=None, diffraction_data=None, debug=False
    ):
        """
        Takes the detector information from the settings class, or the
        calibration *.med file and creates a detector-type instance which converts
        the x,y,azimith of the diffraction data pixels into two-theta vs.
        azimuth.

        Either settings or the file_name are required.

        Parameters
        ----------
        file_name : string, optional
            Filename of the calibration file. In this case a *.med file.
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

        # cannot load calibration into detector simply. Much simpler to re-read the files.
        if hasattr(settings, "calibration_data"):
            self.detector = Med.Med(file=settings.calibration_parameters)
        elif calibration_file != None:
            self.detector = Med.Med(file=calibration_file)

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

        # check inputs
        if image_name == None and settings.subpattern == None:
            raise ValueError("Settings are given but no subpattern is set.")

        if self.detector == None:
            self.get_detector(settings)

        if image_name == None:
            # load the data for the chosen subpattern.
            image_name = settings.subfit_filename

        _, file_extension = os.path.splitext(image_name)
        if file_extension == ".med":
            self.detector.read_file(image_name)
            im_all = []
            for x in range(self.detector.n_detectors):
                im_all.append(self.detector.mcas[x].get_data())
        else:
            sys.exit("Unknown file type. A *.med file is required")

        # Convert the input data from integer to float because the lmfit model values
        # inherits integer properties from the data.
        #
        # Allow option for the data type to be set.
        if dtype == None:
            if self.intensity is None and np.size(self.intensity) > 2:
                # self.intensity has been set before. Inherit the dtype.
                dtype = self.intensity.dtype
            elif "int" in im_all[0].dtype.name:
                # the type is either int or uint - convert to float
                # using same bit precision
                precision = re.findall("\d+", im_all[0].dtype.name)[0]
                dtype = np.dtype("float" + precision)
        im_all = ma.array(im_all, dtype=dtype)

        # apply mask to the intensity array
        if mask == None and ma.is_masked(self.intensity) == False:
            self.intensity = ma.array(im_all)
            return ma.array(im_all)
        elif mask is not None:
            # apply given mask
            self.intensity = ma.array(im_all, mask=self.fill_mask(mask, im_all))
            return ma.array(im_all, mask=mask)
        else:
            # apply mask from intensities
            self.intensity = ma.array(im_all, mask=self.intensity.mask)
            return ma.array(im_all)

    def fill_data(
        self, diff_file=None, settings=None, mask=None, make_zyx=False, debug=False
    ):
        """
        Initiates the data arrays.
        Creates the intensity, two theta, azimuth and d-spacing arrays from the
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
            Switch to make x,y,z arrays for the pixels. These arrays are not
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

        self.tth = self._get_two_theta()
        self.azm = self._get_azimuth()
        self.dspace = self._get_d_space()
        if make_zyx:
            raise ValueError("'make_zyx' is not implemented for MedDetector. ")
            #  FIXME: (SAH 19th June 2024) implement this.

        # get and apply mask
        # catch old formtting of the mask
        if isinstance(mask, list):
            mask = {"detector": mask}
        mask_array = self.get_mask(mask, debug=debug)
        self.mask_apply(mask_array, debug=debug)

        self.azm_start = (
            np.around(np.min(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )
        self.azm_end = (
            np.around(np.max(self.azm.flatten()) / self.azm_blocks) * self.azm_blocks
        )

    @staticmethod
    def detector_check(image_name, settings=None):
        """
        :param image_name:
        :param detector:
        :return:
        """
        # FIX ME: we should try importing the image and reading the detector type from it.
        if settings is None:
            sys.exit(
                "\n\nDioptas requires a detector type.\nTo list the possible detector types in the command line type:"
                "\n   import pyFAI\n   pyFAI.detectors.Detector.registry\n\n"
            )
        # FIX ME: check the detector type is valid.
        return None  # detector

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

    def _get_two_theta(self, mask=None, pix=None, det=None):
        """
        Called get two theta, but it is really getting the energy.
        :param mask:
        :param pix:
        :param det:
        :return:
        """
        channels = np.arange(self.calibration["mcas"].mcas[0].data.shape[0])
        en = []
        for x in range(self.calibration["mcas"].n_detectors):
            en.append(self.calibration["mcas"].mcas[x].channel_to_energy(channels))
        if mask is not None:
            en = ma.array(en, mask=mask)
        else:
            en = ma.array(en)
        return en

    def _get_azimuth(self, mask=None, pix=None, det=None):
        """
        Give azimuth value for detector x,y position; have to assume geometry of the detectors
        :param mask:
        :param pix:
        :param det:
        :return:
        """
        l = self.calibration["mcas"].mcas[0].data.shape
        azi = np.tile(self.calibration["azimuths"], [l[0], 1]).T
        if mask is not None:
            azi = ma.array(azi, mask=mask)
        else:
            azi = ma.array(azi, mask=None)
        return azi

    # =============================================================================
    #  start of functions defined as 'common' for angle dispersive classes.
    #  but the functions here are not common with the angle dispersive versions.
    # =============================================================================

    def _get_d_space(self, mask=None, pix=None, det=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        :param mask:
        :param pix:
        :param det:
        :return:
        """
        channels = np.arange(self.calibration["mcas"].mcas[0].data.shape[0])
        dspc = []
        for x in range(self.calibration["mcas"].n_detectors):
            dspc.append(self.calibration["mcas"].mcas[x].channel_to_d(channels))
        if mask is not None:
            d_space = ma.array(dspc, mask=mask)
        else:
            d_space = ma.array(dspc)
        return d_space

    def conversion(self, e_in, azm=None, reverse=False):
        """
        Convert energy values into d-spacing.
        azm is needed to enable compatibility with the energy dispersive detectors
        :param e_in:
        :param azm:
        :param reverse:
        :return:
        """
        # convert energy into d-spacing.
        if azm is not None:
            e_in = np.array([[azm, e_in]])
            # FIX ME: this is a horrible way to make the conversion work, but it is needed for
            # processing the guesses of MED detectors.
        else:
            e_in = np.array([[0, e_in]], dtype=object)

        # Convert Energy (keV) to wavelength (Angstroms)
        # E = hc/lambda;  h = 4.135667662(25)×10−15 eV s; c = 299792458 m/s
        wavelength = 4.135667662e-15 * (299792458 * 1e10) / (e_in[:, 1] * 1000)

        # determine which detector calibration goes with each energy.
        sort_idx = np.array(self.calibration["azimuths"]).argsort()
        de = sort_idx[
            np.searchsorted(self.calibration["azimuths"], e_in[:, 0], sorter=sort_idx)
        ]

        if not reverse:
            # convert energy to d-spacing
            dspc_out = []
            for i in range(len(de)):
                dspc_out.append(
                    wavelength[i]
                    / 2
                    / np.sin(
                        np.radians(
                            self.calibration["mcas"].mcas[i].calibration.two_theta / 2
                        )
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
                        np.radians(
                            self.calibration["mcas"].mcas[i].calibration.two_theta / 2
                        )
                    )
                )
            dspc_out = np.array(dspc_out)

        return dspc_out

    def bins(self, orders_class, **kwargs):
        """
        Determine bins to use in initial fitting.
        Assign each data to a chunk corresponding to its azimuth value
        Returns array with indices for each bin and array of bin centroids
        :param azimu:
        :param orders:
        :param **kwargs: - to ensure compatibility
        :return:
        """
        bin_mean_azi = np.unique(self.azm.data)
        chunks = []
        # azichunks = []
        temp_azimuth = self.azm.flatten()
        for i in range(len(bin_mean_azi)):
            azi_chunk = np.where((temp_azimuth == bin_mean_azi[i]))
            chunks.append(azi_chunk)

        return chunks, bin_mean_azi

    def test_azims(self, steps=None):
        """
        Returns unique azimuths in the data set.

        Parameters
        ----------
        steps : None, optional
            DESCRIPTION. Does nothing in this case. Maintained for compatibility

        Returns
        -------
        array
            list of possible azimuths.

        """
        return np.unique(self.azm)

    # =============================================================================
    # Start of plotting functions.
    # Not common with the angle dispersive plotting
    # =============================================================================

    def _dispersion_ticks(self, disp_ticks=None, unique=10, disp_lims=None):
        """
        Returns the labels for the dispersion axis/colour bars.

        :param disp_ticks:
        :param unique:
        :param disp_lims: not used maintained for compatibility with angle dispersive
        :return new_tick_positions:
        """
        if len(np.unique(self.azm)) <= unique:
            disp_ticks = np.unique(np.unique(self.azm))

        return disp_ticks

    def _plot_spacing(self, gap=200, Imax=None):
        """
        Sets the spacing for the detectors in the plots.

        Parameters
        ----------
        gap : number, optional
            Desired spacing for the azimuths. The default is 200.

        Returns
        -------
        number
            Spacing to use for the plotting of azimuths.

        """
        if Imax:
            return Imax / gap
        else:
            return np.max(self.intensity) / gap

    def plot_masked(self, fig_plot=None):
        """
        Plot all the information needed to mask the data well.
        :param fig:
        :return:
        """

        ax1 = fig_plot.add_subplot(111)
        self.plot_calibrated(
            fig_plot=fig_plot, axis_plot=ax1, x_axis="default", limits=[0, 100]
        )

    def plot_fitted(
        self,
        fig_plot=None,
        model=None,
        fit_centroid=None,
        plot_ColourRange=None,
        **kwargs,
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """
        # reshape model to match self
        model2 = deepcopy(self.intensity)
        p = 0

        for i in range(len(self.intensity)):
            # for j in range(len(self.intensity[0])):
            if not ma.is_masked(self.intensity[i]):
                model2[i] = model.flatten()[p]
                p = p + 1

        # match max and min of colour scales
        # FIXME: this is not used but should be used to set to colour ranges of the data -- by azimuth.
        # limits = {"max": self.azm_start, "min": self.azm_end}

        if plot_ColourRange:
            spcng = (self._plot_spacing(gap=200, Imax=plot_ColourRange["max"]),)
        else:
            spcng = self._plot_spacing(gap=200)

        # make axes
        gs = gridspec.GridSpec(1, 3, wspace=0.0)
        ax = []
        for i in range(3):
            ax.append(fig_plot.add_subplot(gs[i]))

        # plot data
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[0],
            show="intensity",
            x_axis="default",
            limits=deepcopy(plot_ColourRange),
            colourbar=False,
        )
        ax[0].set_title("Data")
        # plot model
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                ax[1].plot(
                    np.squeeze(fit_centroid[1][i]),
                    spcng * np.unique(self.azm),
                    ".k--",
                    linewidth=0.5,
                )
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[1],
            data=model2,
            colourbar=False,
            limits=deepcopy(plot_ColourRange),
        )
        ax[1].set_title("Model")
        # plot residuals
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                ax[2].plot(
                    np.squeeze(fit_centroid[1][i]),
                    spcng * np.unique(self.azm),
                    ".k--",
                    linewidth=0.5,
                )
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[2],
            data=self.intensity - model2,
            limits=deepcopy(plot_ColourRange),
            colourbar=True,
        )
        ax[2].set_title("Residuals")

        # organise the axes and labelling.
        bottom0, top0 = ax[0].get_ylim()
        bottom1, top1 = ax[1].get_ylim()
        bottom2, top2 = ax[2].get_ylim()
        top_max = np.max([top0, top1, top2])
        bottom_min = np.min([bottom0, bottom1, bottom2])
        ax[0].set_ylim(top=top_max, bottom=bottom_min)
        ax[1].set_ylim(top=top_max, bottom=bottom_min)
        ax[2].set_ylim(top=top_max, bottom=bottom_min)
        ax[2].set_ylim(top=top_max, bottom=bottom_min)
        ax[1].set(yticklabels=[])
        ax[2].set(yticklabels=[])
        ax[1].set(ylabel=None)
        ax[2].set(ylabel=None)
        ax[1].set(ylabel=None)

        # # Colorbar
        # # set the colour bar underneath the axes.
        # blocks = 90
        # colour_start = np.floor((np.min(self.azm)-.1)/blocks)*blocks
        # colour_end   = np.ceil((np.max(self.azm))/blocks)*blocks
        # normalize = colors.Normalize(vmin=colour_start, vmax=colour_end)
        # c_map = cm.get_cmap(name="jet")
        # s_map = cm.ScalarMappable(norm=normalize, cmap=c_map)
        # # label colour bar with azimuths if there are less than 10
        # if len(self.intensity) <= 10:
        #     ticks = np.unique(self.azm)
        # else:
        #     ticks = None
        # cbar = fig_plot.colorbar(ax[1], ax=[ax[0], ax[1], ax[2]], s_map, ticks=ticks, orientation="horizontal")
        # cbar.set_label("Azimuth")

        # tidy layout
        plt.tight_layout()

    def plot_collected(
        self,
        fig_plot=None,
        axis_plot=None,
        show="intensity",
        x_axis="default",
        limits=[0, 100],
        spacing=None,
    ):
        """
        add data to axes in form collected in.
        :param fig_plot:
        :param axis_plot:
        :param show:
        :param x_axis:
        :param limits:
        :return:
        """

        if spacing == None:
            spacing = self._plot_spacing(gap=10)

        for x in range(len(self.intensity)):
            axis_plot.plot(
                list(range(len(self.intensity[x]))), self.intensity[x] + spacing * x
            )
        axis_plot.set_xlabel("Channel number")
        axis_plot.set_ylabel("Counts")

    def plot_calibrated(
        self,
        fig_plot=None,
        axis_plot=None,
        show="intensity",
        x_axis="default",
        y_axis="default",
        spacing=None,
        data=None,
        limits=None,
        colourmap="jet",
        colourbar=True,
        debug=False,
    ):
        """
        add data to axes in form collected in.
        :param fig_plot:
        :param axis_plot:
        :param show:
        :param x_axis:
        :param limits:
        :return:
        """

        if spacing == None:
            if isinstance(limits, dict):
                spacing = self._plot_spacing(gap=200, Imax=limits["max"])
                lims = [limits["rmin"], limits["max"] + np.max(self.azm) * spacing]
                limits["max"] = limits["max"] + np.max(self.azm) * spacing
                # spacing = self._plot_spacing(gap=200, Imax=limits["max"])
                #
            else:
                spacing = self._plot_spacing(gap=200)
                lims = None

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
            # y is unmodified intensity; colour of data by azimuth
            plot_y = plot_i
            plot_c = self.azm
            label_y = "Intensity (a.u.)"
        elif y_axis == "azimuth":
            # y is unmodified azimuth, intensity is size of dots amd colour scale
            plot_y = self.azm
            plot_c = plot_i / np.max(self.intensity) * 200
            plot_s = plot_i / np.max(self.intensity) * 200
            label_y = "Azimuth (deg)"
        else:  # if y_axis is default
            # y is intensity distributed by azimuth; colour is azimuth
            plot_y = plot_i + self.azm * spacing
            plot_c = self.azm
            label_y = "Counts"

            # organise to plot 0 count lines.
            plot_y0 = []
            plot_x0 = []
            plot_c0 = []
            for i in range(len(np.unique(self.azm))):
                if ma.MaskedArray.all(plot_x[self.azm == np.unique(self.azm)[i]]):
                    plot_x0.append(
                        [
                            plot_x[self.azm == np.unique(self.azm)[i]][0],
                            plot_x[self.azm == np.unique(self.azm)[i]][-1],
                        ]
                    )
                    plot_y0.append(
                        [
                            self.azm[self.azm == np.unique(self.azm)[i]][0] * spacing,
                            self.azm[self.azm == np.unique(self.azm)[i]][-1] * spacing,
                        ]
                    )
                    plot_c0.append(np.mean(plot_c[self.azm == np.unique(self.azm)[i]]))
                else:
                    plot_y0.append([np.nan, np.nan])
                    plot_x0.append([np.nan, np.nan])
            plot_x0 = np.array(plot_x0)
            plot_y0 = np.array(plot_y0)

        # organise the data to plot
        # here assume x isenergy, y is azzimuth, z is intensity or counts.
        if show == "unmasked_intensity":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_c = plot_c
        elif show == "mask":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_c = np.array(ma.getmaskarray(self.intensity), dtype="uint8") + 1
            colourmap = "Greys"
        else:
            plot_c = plot_c

        # set colour map
        if colourmap == "residuals-blanaced":
            colourmap = self.residuals_colour_scheme(
                np.max(plot_c.flatten()), np.min(plot_c.flatten())
            )
        else:
            colourmap = colourmap

        if "plot_s" in locals():
            the_plot = axis_plot.scatter(
                plot_x, plot_y, s=plot_s, c=plot_c, cmap=colourmap, alpha=0.2
            )
            fig_plot.colorbar(mappable=the_plot).set_label(label="Counts")
        else:
            # Set colour map range - to be in 90 degree blocks using azm_start and azm_end
            normalize = colors.Normalize(vmin=self.azm_start, vmax=self.azm_end)
            c_map = cm.get_cmap(name=colourmap)
            for i in range(len(np.unique(self.azm)) - 1, -1, -1):
                colour = c_map(
                    normalize(np.mean(plot_c[self.azm == np.unique(self.azm)[i]]))
                )
                if y_axis == "default":
                    the_plot = axis_plot.plot(
                        plot_x0[i], plot_y0[i], "--", linewidth=0.5, color=colour
                    )
                the_plot = axis_plot.plot(
                    plot_x[self.azm == np.unique(self.azm)[i]],
                    plot_y[self.azm == np.unique(self.azm)[i]],
                    color=colour,
                )

            # Colorbar
            if colourbar == "horizontal":
                orientation = "horizontal"
            else:
                orientation = "vertical"
            s_map = cm.ScalarMappable(norm=normalize, cmap=c_map)
            # label colour bar with unique azimuths if there are less than 10
            # set colour bar labels with unique azimuths (if there are less than 'unique' azimuths - see function for value of unique).
            ticks = self._dispersion_ticks()
            if colourbar is True:
                cbar = plt.colorbar(
                    s_map, ticks=ticks, orientation=orientation, ax=axis_plot
                )
                cbar.set_label("Azimuth")
            else:
                cbar = []

        axis_plot.set_xlabel(label_x)
        axis_plot.set_ylabel(label_y)
        if lims:
            axis_plot.set_ylim(lims)
        elif isinstance(limits, dict):
            axis_plot.set_ylim([limits["min"], limits["max"]])

        return the_plot, cbar

    # add common functions.
    set_limits = _AngleDispersive_common.set_limits

    # add masking functions to detetor class.
    get_mask = _masks.get_mask
    set_mask = _masks.set_mask
    mask_apply = _masks.mask_apply
    mask_restore = _masks.mask_restore
    mask_remove = _masks.mask_remove
