#!/usr/bin/env python

__all__ = ["execute", "write_output"]

import json
import os
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import multiprocessing as mp
import importlib.util
import cpf.DioptasFunctions as DioptasFunctions
import cpf.GSASIIFunctions as GSASIIFunctions
import cpf.MedFunctions as MedFunctions
from cpf.XRD_FitSubpattern import fit_sub_pattern
from cpf.Cosmics import cosmicsimage
from cpf.IO_functions import json_numpy_serializer, file_list, make_outfile_name, peak_string

np.set_printoptions(threshold=sys.maxsize)
# FIX ME: Need to add complexity here.
# Need option to check required functions are present if adding new output format.
# Also need to do something similar with data class


def register_default_formats() -> object:
    """
    Load all available output modules
    :return:
    """
    # FIX ME: We could add extra checks here to make sure the required functions exist in each case.
    # FIX ME: Really these should live in their own sub-folder
    output_list = [
        "CoefficientTable",
        "DifferentialStrain",
        "MultiFit",
        "Polydefix",
        "PolydefixED",
    ]
    new_module = {}
    for output_module in output_list:
        module = __import__("cpf.Write" + output_module, fromlist=[None])
        new_module[output_module] = module
    return new_module


# Load potential output formats
output_methods_modules = register_default_formats()


def get_output_options(output_type):
    """
    Check if input is string or list of strings
    :param output_type: string or list of strings
    :return: list of strings
    """
    output_mod_type = []
    if isinstance(output_type, str):
        output_mod_type.append(output_type)
    else:
        output_mod_type = output_type
    return output_mod_type


def detector_factory(calibration_param, calibration_type, fit_settings=None):
    """
    Factory function to provide appropriate class for data dependent on type.
    Currently, supported options are Dioptas, GSASII and Med.
    :rtype: object
    :param fit_settings:
    :param calibration_param:
    :param calibration_type:
    :return:
    """
    if calibration_type == "Dioptas":
        detector_class = DioptasFunctions.DioptasDetector
        return detector_class(calibration_param, fit_settings)
    if calibration_type == "GSASII":
        detector_class = GSASIIFunctions.GSASIIDetector
        return detector_class(calibration_param, fit_settings)
    if calibration_type == "Med":
        detector_class = MedFunctions.MedDetector
        return detector_class(calibration_param, fit_settings)
    else:
        raise ValueError("Unrecognized calibration type.")


def initiate(
    setting_file=None,
    inputs=None,
    out_type=None,
    report=False,
    **kwargs
):
    """
    Run checks on input files, initiate data class and check output options
    :param report:
    :param out_type:
    :param initiate_data:
    :param setting_file:
    :param inputs:
    :return fit_parameters:
    :return fit_settings:
    """
    # Check we can find and load settings file
    if settings_file:
        try:
            print("\nThe name of the settings file is: ", setting_file)
            if not str.endswith(setting_file, ".py"):
                setting_file = setting_file + ".py"
            spec = importlib.util.spec_from_file_location(
                "settings_module", setting_file
            )
            fit_settings = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(fit_settings)
            fit_parameters = dir(fit_settings)
            fit_settings.inputfile = setting_file

        except FileNotFoundError:
            raise FileNotFoundError
    elif inputs:
        raise NotImplementedError
    else:
        raise ValueError("No inputs provided")

    # Data directory and file validation: check they exist.
    if os.path.exists(fit_settings.datafile_directory) is False:
        raise ImportError(
            "The data directory " + fit_settings.datafile_directory + " does not exist."
        )
    else:
        print("The data directory exists.")
    diff_files, n_diff_files = file_list(fit_parameters, fit_settings)
    for j in range(n_diff_files):
        if os.path.isfile(diff_files[j]) is False:
            raise ImportError(
                "The data file "
                + diff_files[j]
                + " does not exist. Check the input file."
            )
    print("All the data files exist.")

    # Confirm that an output directory is listed.
    if "Output_directory" not in fit_parameters:
        fit_settings.Output_directory = "."
        fit_parameters.append("Output_directory")
    if not os.path.isdir(fit_settings.Output_directory):
        raise ValueError("Output directory not found")

    # Load the detector class here to access relevant functions and check required parameters are present
    if "Calib_type" not in fit_parameters:
        raise ValueError(
            "There is no 'Calib_type' in the settings. The fitting cannot proceed until a recognised "
            "calibration type is present."
        )
    else:
        data_class = detector_factory(fit_settings.Calib_param, fit_settings.Calib_type)

    # FIX ME: This doesn't seem to be used, if it should be this needs moving to class structure.
    alternatives_list = [[["datafile_StartNum", "datafile_EndNum"], ["datafile_Files"]]]
    optional_list = ["datafile_Step"]  # FIX ME: This doesn't seem to be used
    possible = [[[], []] * len(alternatives_list)]
    for x in range(len(alternatives_list)):
        for y in range(2):
            for z in range(len(alternatives_list[x][y])):
                if alternatives_list[x][y][z] in fit_parameters:
                    possible[x][y].append(1)
                else:
                    possible[x][y].append(0)
    # exit if all parameters are not present

    # check the peak fitting options in the input file are not illicit.
    missing = []
    extras = []
    fit_orders = fit_settings.fit_orders

    for i in range(len(fit_orders)):
        # FIX ME: we should check for incorrect or unneeded options
        required = ["background", "peak", "range"]
        possible = ["PeakPositionSelection", "Imax"]
        peak_required = ["d-space", "width", "height", "profile"]

        # check range and background
        if "range" not in fit_orders[i]:
            missing.append("fit_orders " + str(i) + " is missing a " "'range'")
        elif (
                not isinstance(fit_orders[i]["range"], list)
                and len(fit_orders[i]["range"]) != 2
        ):
            missing.append(
                "fit_orders " + str(i) + " has an incorrectly formatted "+"'range'")

        if "background" not in fit_orders[i]:
            missing.append("fit_orders " + str(i) + " is missing a "+"'background'")
        elif not isinstance(fit_orders[i]["background"], list):
            missing.append(
                "fit_orders " + str(i) + " has an incorrectly formatted"+"'background'")

        # check peaks
        if "peak" not in fit_orders[i]:
            missing.append("fit_orders" + str(i) + "has no "+"'peak'")
        else:
            for j in range(len(fit_orders[i]["peak"])):
                for k in range(len(peak_required)):
                    if not peak_required[k] in fit_orders[i]["peak"][j]:
                        missing.append(
                            "fit_orders "
                            + str(i)
                            + ", peak "
                            + str(j)
                            + " has no "
                            + peak_required[k]
                        )
                    elif not isinstance(
                        fit_orders[i]["peak"][j][peak_required[k]], list
                    ) and not isinstance(
                        fit_orders[i]["peak"][j][peak_required[k]], int
                    ):
                        missing.append(
                            "fit_orders "
                            + str(i)
                            + ", peak "
                            + str(j)
                            + " incorrectly formatted "
                            "" + peak_required[k] + " "
                            " "
                        )

        # if PeakPositionSelection - check within range
        if "PeakPositionSelection" in fit_orders[i]:
            # how many peaks in list
            tthguesses = np.array(fit_orders[i]["PeakPositionSelection"])

            # FIX ME: use of j here outside of loop - need to check meaning!
            if max(tthguesses[:, 0]) > len(fit_orders[i]["peak"][j]):
                missing.append(
                    "fit_orders "
                    + str(i)
                    + ": PeakPositionSelection contains too many peaks"
                )

            # if positions outside of range.
            if np.min(tthguesses[:, 2]) < fit_orders[i]["range"][0][0]:
                missing.append(
                    "fit_orders "
                    + str(i)
                    + ": PeakPositionSelection has at least one two theta value that is too small"
                )
            if np.max(tthguesses[:, 2]) > fit_orders[i]["range"][0][1]:
                missing.append(
                    "fit_orders "
                    + str(i)
                    + ": PeakPositionSelection has at least one two theta value that is too large"
                )
        else:
            if len(fit_orders[i]["peak"]) > 1:
                missing.append(
                    "fit_orders " + str(i) + ": There are more peaks than listed in "
                    "PeakPositionSelection"
                    " than peaks listed"
                )

        # list unrecognised entries

    if missing:
        print("\nMissing Values:")
        for i in range(len(missing)):
            print(missing[i])
        if not report:
            raise ValueError("The problems listed above will prevent the data fitting.")
        else:
            print(
                "The problems listed above will prevent the data fitting and need to be rectified before execution"
            )
    else:
        print("Fit_orders appears to be correct")

    # FIX ME: do we need to check fit bounds?

    # Load output function(s) and check output options are present and valid
    if out_type is not None:
        output_mod = get_output_options(out_type)
    elif "Output_type" in fit_parameters:
        output_mod = get_output_options(fit_settings.Output_type)
    else:
        raise ValueError(
            "No output type. Add 'Output_type' to input file or specify 'out_type' in command."
        )

    # Check output format exists
    for mod in output_mod:
        if mod not in output_methods_modules.keys():
            raise ImportError(
                "The 'Output_type' "
                + mod
                + " is not recognised; the file '"
                + mod
                + "' does not "
                "exist. Check if "
                "the calibration "
                "type exists."
            )

    return fit_settings, fit_parameters, data_class


def set_range(
    setting_file=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    propagate=True,
    iterations=1,
    track=False,
    parallel=True,
    subpattern="all",
    **kwargs
):
    """
    :param setting_file:
    :param inputs:
    :param debug:
    :param refine:
    :param save_all:
    :param propagate:
    :param iterations:
    :param track:
    :param parallel:
    :param subpattern:
    :param kwargs:
    :return:
    """
    fit_settings, fit_parameters, new_data = initiate(
        setting_file, inputs=inputs, report=True
    )

    if "datafile_Files" not in fit_parameters:
        fit_settings.datafile_Files = fit_settings.datafile_Files[0]
    # search over the first file only
    # strip file list to first file
    else:
        fit_settings.datafile_EndNum = fit_settings.datafile_StartNum

    # restrict to sub-patterns listed
    if subpattern == "all":
        sub_pats = list(range(0, len(fit_settings.fit_orders)))
    elif isinstance(subpattern, list):
        sub_pats = subpattern
    else:
        sub_pats = [int(x) for x in str(subpattern)]

    # make new order search list
    orders_tmp = []
    for i in range(len(sub_pats)):
        j = sub_pats[i]
        orders_tmp.append(fit_settings.fit_orders[j])
    fit_settings.fit_orders = orders_tmp

    execute(inputs=new_data, debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
            parallel=parallel, mode="set-range", report=True, fit_settings=fit_settings, fit_parameters=fit_parameters)


def initial_peak_position(
    setting_file=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    propagate=True,
    iterations=1,
    track=False,
    parallel=True,
    sub_pattern="all",
    **kwargs
):
    """
    :param setting_file:
    :param inputs:
    :param debug:
    :param refine:
    :param save_all:
    :param propagate:
    :param iterations:
    :param track:
    :param parallel:
    :param sub_pattern:
    :param kwargs:
    :return:
    """

    # see https://matplotlib.org/stable/users/event_handling.html for how to make work
    fit_settings, fit_parameters, new_data = initiate(
        setting_file, inputs=inputs, report=True
    )

    # search over the first file only
    # strip file list to first file
    if "datafile_Files" not in fit_parameters:
        fit_settings.datafile_EndNum = fit_settings.datafile_StartNum
    else:
        fit_settings.datafile_Files = fit_settings.datafile_Files[0]

    # restrict to sub-patterns listed
    if sub_pattern == "all":
        sub_pats = list(range(0, len(fit_settings.fit_orders)))
    elif isinstance(sub_pattern, list):
        sub_pats = sub_pattern
    else:
        sub_pats = [int(x) for x in str(sub_pattern)]

    # make new order search list
    orders_tmp = []
    for i in range(len(sub_pats)):
        j = sub_pats[i]
        orders_tmp.append(fit_settings.fit_orders[j])
    fit_settings.fit_orders = orders_tmp

    execute(fit_settings=fit_settings, fit_parameters=fit_parameters, inputs=new_data, debug=debug, refine=refine,
            save_all=save_all, propagate=propagate, iterations=iterations, parallel=parallel, mode="set-guess")


class PointBuilder:
    # copied from https://matplotlib.org/stable/users/event_handling.html on 21 July 2021
    def __init__(self, points, fig):
        self.points = points
        self.lines = []
        self.ax = fig
        self.xs = list(points.get_xdata())
        self.ys = list(points.get_ydata())
        self.ks = list([])
        self.cid1 = points.figure.canvas.mpl_connect("button_press_event", self)
        self.cid2 = points.figure.canvas.mpl_connect("key_press_event", self)

    def __call__(self, event):
        # print('click', event)
        if event.inaxes != self.points.axes:
            return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)

        try:
            new_event = event.button
        except:
            new_event = event.key
        # replace mouse clicks with numbers and make string numbers in to numbers.
        # left mouse --> 1
        # right mouse --> 2
        # number (0-9) --> number 0-9 as number not string
        # other characters (a-z, etc...) --> ASCII equivalent.
        if isinstance(new_event, str):  # keyboard button press
            try:
                new_event = int(new_event)
            except:
                new_event = ord(new_event)  # replace letter with its ASCII value
        elif new_event == 1:  # left mouse button
            new_event = 1
        elif new_event == 3:  # right mouse button
            new_event = 2
        self.ks.append(new_event)

        order = np.argsort(np.array(self.ys))
        self.xs = [self.xs[i] for i in order]
        self.ys = [self.ys[i] for i in order]
        self.ks = [self.ks[i] for i in order]

        sets = list(set(self.ks))
        self.lines = []
        for i in range(len(sets)):
            index = []
            j = 0
            while j < len(self.ks):
                if sets[i] == self.ks[j]:
                    index.append(j)
                j += 1
            index = np.array(index, dtype="int_")
            self.lines = self.ax.scatter(
                [self.xs[j] for j in index], [self.ys[j] for j in index]
            )
            # self.ax.plot(px,py, linewidth=2, marker='o')

        self.lines.figure.canvas.draw()

    def array(self):
        # make array for input file -- and sort it

        # replace mouse clicks with numbers and make string numbers in to numbers.
        # left mouse --> 1
        # right mouse --> 2
        # number (0-9) --> number 0-9 as number not string
        # other characters (a-z, etc...) --> ASCII equivalent.
        for i in range(len(self.ks)):
            if isinstance(self.ks[i], str):  # keyboard button press
                try:
                    self.ks[i] = int(self.ks[i])
                except:
                    self.ks[i] = ord(self.ks[i])  # replace letter with its ASCII value
            elif self.ks[i] == 1:  # left mouse button
                self.ks[i] = 1
            elif self.ks[i] == 3:  # right mouse button
                self.ks[i] = 2

        self.array = []
        for i in range(len(self.xs)):
            self.array.append([self.ks[i], self.ys[i], self.xs[i]])
        self.array.sort()

        return self.array


def order_search(
    setting_file=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    propagate=True,
    iterations=1,
    track=False,
    parallel=True,
    search_parameter="height",
    search_over=[0, 20],
    sub_pattern="all",
    search_peak=0,
    search_series=["fourier", "spline"],
    **kwargs
):
    """
    :param search_series:
    :param search_peak:
    :param sub_pattern:
    :param search_over:
    :param track:
    :param search_parameter:
    :param parallel:
    :param propagate:
    :param save_all:
    :param setting_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """
    fit_settings, fit_parameters, new_data = initiate(setting_file, inputs=inputs)

    # force it to write the required output type.
    fit_settings.Output_type = "DifferentialStrain"

    # search over the first file only
    # strip file list to first file
    if "datafile_Files" not in fit_parameters:
        fit_settings.datafile_EndNum = fit_settings.datafile_StartNum
    else:
        fit_settings.datafile_Files = fit_settings.datafile_Files[0]

    # restrict to sub-patterns listed
    if sub_pattern == "all":
        sub_patterns = list(range(0, len(fit_settings.fit_orders)))
    elif isinstance(sub_pattern, list):
        sub_patterns = sub_pattern
    else:
        sub_patterns = [int(x) for x in str(sub_pattern)]
    # make new order search list
    if isinstance(search_over, list) and len(search_over) == 2:
        search = list(range(search_over[0], search_over[1]))
    else:
        search = [int(x) for x in str(search_over)]

    ordered_search = []
    for i in range(len(sub_patterns)):
        for j in range(len(search_series)):
            tmp_order = fit_settings.fit_orders[sub_patterns[i]]
            for k in range(len(search)):
                orders_s = deepcopy(tmp_order)
                if search_parameter != "background":
                    orders_s["peak"][search_peak][search_parameter] = search[k]
                    orders_s["peak"][search_peak][
                        search_parameter + "-type"
                    ] = search_series[j]
                else:
                    orders_s["background"][search_peak] = search[k]
                orders_s["note"] = (
                    search_parameter
                    + "="
                    + str(search[k])
                    + "; type="
                    + search_series[j]
                )
                ordered_search.append(orders_s)
    fit_settings.fit_orders = ordered_search

    execute(fit_settings=fit_settings, fit_parameters=fit_parameters, inputs=new_data, debug=debug, refine=refine,
            save_all=save_all, propagate=propagate, iterations=iterations, parallel=parallel)

    # write a differential strain output file
    # FIX ME: using debug as a switch to get all info in file.
    # should probably use another switch
    write_output(setting_file, fit_settings=fit_settings, fit_parameters=fit_parameters, debug=True,
                 outtype="DifferentialStrain")


def write_output(
    setting_file=None,
    fit_settings=None,
    fit_parameters=None,
    parms_dict=None,
    out_type=None,
    det=None,
    use_bounds=False,
    differential_only=False,
    debug=False,
    **kwargs
):
    """
    :param debug:
    :param fit_parameters:
    :param use_bounds:
    :param differential_only:
    :param setting_file:
    :param fit_settings:
    :param parms_dict:
    :param out_type:
    :param det:
    :return:
    """
    # Fail gracefully
    if (
        parms_dict is None
        and setting_file is None
        and fit_settings is None
        and fit_parameters is None
    ):
        raise ValueError(
            "Either the settings file or the parameter dictionary need to be specified."
        )
    # If no params_dict then initiate. Check all the output functions are present and valid.
    if parms_dict is None and fit_settings is None and fit_parameters is None:
        fit_settings, fit_parameters, new_data = initiate(setting_file, outtype=out_type)
        parms_dict = new_data.get_calibration(os.path.abspath(fit_settings.Calib_param))
    elif parms_dict is None:
        _, _, new_data = initiate(setting_file, outtype=out_type)
        parms_dict = new_data.get_calibration(os.path.abspath(fit_settings.Calib_param))
    if out_type is not None:
        output_mod = get_output_options(out_type)
    elif "Output_type" in fit_parameters:
        output_mod = get_output_options(fit_settings.Output_type)
    else:
        raise ValueError(
            "No output type. Add 'Output_type' to input file or specify 'out_type' in command."
        )

    for mod in output_mod:
        print("\nWrite output file(s) using", mod)
        wr = output_methods_modules[mod]
        wr.WriteOutput(
            fit_settings, parms_dict, differential_only=differential_only, debug=debug
        )


def execute(
    setting_file=None,
    fit_settings=None,
    fit_parameters=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    propagate=True,
    iterations=1,
    track=False,
    parallel=True,
    mode="fit",
    report=False,
    **kwargs
):
    """
    :param fit_parameters: 
    :param fit_settings: 
    :param setting_file:
    :param parallel:
    :param report:
    :param mode:
    :param track:
    :param propagate:
    :param save_all:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """
    if fit_settings is None:
        fit_settings, fit_parameters, new_data = initiate(setting_file, report=report)
    else:
        new_data = inputs

    # Define locally required names
    temporary_data_file = make_outfile_name(
        "PreviousFit_JSON",
        directory=fit_settings.Output_directory,
        extension=".dat",
        overwrite=True,
    )
    # temporary_data_file = make_outfile_name('PreviousFit_JSON', extension='dat', overwrite=True)

    # Get list of diffraction patterns #
    diff_files, n_diff_files = file_list(fit_parameters, fit_settings)

    # Initiate data in data class with appropriate detector class
    if "Calib_mask" in fit_parameters:
        use_mask = fit_settings.Calib_mask
    else:
        use_mask = None

    if "Calib_data" in fit_parameters:
        new_data.fill_data(
            os.path.abspath(fit_settings.Calib_data),
            settings=fit_settings,
            debug=debug,
            calibration_mask=use_mask,
        )
    else:
        new_data.fill_data(
            os.path.abspath(diff_files[0]),
            settings=fit_settings,
            debug=debug,
            calibration_mask=use_mask,
        )

    # Get calibration parameter file
    parms_dict = new_data.parameters

    # Get calibration for image set.
    # Get array for azimuth, two theta and d-space the same size as the image.
    azimuth = new_data.azm
    twotheta = new_data.tth
    dspace = new_data.dspace
    original_mask = np.copy(new_data.azm.mask)
    # print(azimuth, twotheta, dspace)

    # plot calibration file
    if (
        0
    ):  # debug and 'Calib_data' in FitParameters:
        # FIX ME: SAH - I have just removed this so that the debug runs without the calibration file.
        new_data.debug_plot(two_theta=twotheta, azimuth=azimuth, intensity=new_data.intensity)

    # if parallel processing start the pool
    if parallel is True:
        p = mp.Pool(processes=mp.cpu_count())

    # Process the diffraction patterns #
    for j in range(n_diff_files):
        print("Process ", diff_files[j])

        # Get diffraction pattern to process.
        intens = new_data.import_image(diff_files[j], mask=azimuth.mask, debug=debug)
        # Replace this with call through to class structure
        # FIX ME: the mask call is needed here to pass the mask in. but should it inherit the mask from the preceding
        # data?
        # DMF FIX ME: No test data seems to use the image_prepare setting. What should be the equivalent when it's not
        # used?

        # FIX ME: replace mask with image_prepare mask.
        if "Image_prepare" in fit_parameters:
            if "cosmics" in fit_settings.Image_prepare:
                print("Remove Cosmics")
                # set defaults
                gain = 2.2
                sigclip = 1.5  # 3 is dioptas default
                objlim = 3.0  # 3 is dioptas default
                sigfrac = 0.3
                if bool(fit_settings.Image_prepare["cosmics"]):
                    if "gain" in fit_settings.Image_prepare["cosmics"]:
                        gain = fit_settings.Image_prepare["cosmics"]["gain"]
                    if "sigclip" in fit_settings.Image_prepare["cosmics"]:
                        sigclip = fit_settings.Image_prepare["cosmics"]["sigclip"]
                    if "objlim" in fit_settings.Image_prepare["cosmics"]:
                        objlim = fit_settings.Image_prepare["cosmics"]["objlim"]
                    if "sigfrac" in fit_settings.Image_prepare["cosmics"]:
                        sigfrac = fit_settings.Image_prepare["cosmics"]["sigfrac"]
                    # FIX ME: use argparse or someway of passing any argument into cosmics.

                test = cosmicsimage(
                    intens, sigclip=sigclip, objlim=objlim, gain=gain, sigfrac=sigfrac
                )
                num = 2
                for i in range(num):
                    test.lacosmiciteration(True)
                    # print(asdf["nnew"], asdf["niter"])
                    test.clean()
                    msk = np.logical_or(original_mask, np.array(test.mask, dtype="bool"))
                    # need to keep original mask and keep referencing to it because we overwrite the mask in azimuth
                    # below.
                intens = ma.array(intens, mask=msk)
                azimuth = ma.array(azimuth, mask=intens.mask)
                twotheta = ma.array(twotheta, mask=intens.mask)
                dspace = ma.array(dspace, mask=intens.mask)
                # FIX ME: applying a mask to the data should mask all the arrays. This should be a class function.
                # FIX ME: I think these calls should all be class applications rather than as it is now. -- but it seems
                # to work

        # plot input file
        if debug:
            new_data.debug_plot(os.path.basename(diff_files[j]), twotheta, azimuth, intens, show=1)
            '''
            fig = plt.figure()
            if fit_settings.Calib_type == "Med":
                # Better for energy dispersive data.  #FIX ME: need variable for vmax.
                for x in range(9):
                    plt.plot(
                        twotheta[x], intens[x] + 100 * x
                    )  # Better for energy dispersive data.  #FIX ME: need variable for vmax.
                plt.xlabel("Energy (keV)")
                plt.ylabel("Intensity")
            else:
                # FIX ME: need variable for vmax.
                plt.scatter(
                    twotheta,
                    azimuth,
                    s=4,
                    c=intens,
                    edgecolors="none",
                    cmap=plt.cm.jet,
                    vmin=0,
                    vmax=4000,
                )  # FIX ME: need variable for vmax.
                plt.xlabel(r"2$\theta$")
                plt.ylabel("Azimuth")
                plt.colorbar()
            plt.title(os.path.basename(diff_files[j]))
            plt.draw()
            plt.close()
            '''

        # Get previous fit (if it exists and is required)
        if os.path.isfile(temporary_data_file) and propagate is True and mode == "fit":
            # Read JSON data from file
            print("Loading previous fit results from %s" % temporary_data_file)
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)

        # Switch to save the first fit in each sequence.
        if j == 0 or save_all is True:
            save_figs = 1
        else:
            save_figs = 0

        fit_orders = fit_settings.fit_orders
        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        # Get number of sub-patterns to be fitted
        n_sub_patterns = len(fit_orders)
        fitted_param = []
        lmfit_models = []
        parallel_pile = []

        for i in range(n_sub_patterns):
            tth_range = fit_orders[i]["range"][0]
            if hasattr(fit_settings, "fit_orders"):  # FitSettings.fit_orders:
                orders = fit_orders[i]
            if hasattr(fit_settings, "AziBins"):
                orders["AziBins"] = fit_settings.AziBins
            if "fit_bounds" in fit_parameters:
                bounds = fit_settings.fit_bounds
            else:
                bounds = None
            if "previous_fit" in locals() and mode == "fit":
                params = previous_fit[i]
            else:
                params = []

            # Track the position of the peak centroid
            # FIX ME; This is crude - the range doesn't change width. so can't account for massive change in stress.
            # But does it need to?
            if track is True and "previous_fit" in locals():
                print("tth_range", tth_range)
                mid = []
                for k in range(len(params["peak"])):
                    mid.append(params["peak"][k]["d-space"][0])
                # Replace this with call through to class structure?
                cent = new_data.conversion(
                    np.sum(mid) / (k + 1), parms_dict, reverse=1, azm=None
                )
                tth_range[0] = tth_range[0] - ((tth_range[1] + tth_range[0]) / 2) + cent
                tth_range[1] = tth_range[1] - ((tth_range[1] + tth_range[0]) / 2) + cent
                print("move by:", ((tth_range[1] + tth_range[0]) / 2) + cent)
                print("tth_range", tth_range)

                # The PeakPositionSelections are only used if the fits are not being propagated
                if "PeakPositionSelection" in orders:
                    for k in range(len(orders["PeakPositionSelection"])):
                        orders["PeakPositionSelection"][k][2] = (
                            orders["PeakPositionSelection"][k][2]
                            - ((tth_range[1] + tth_range[0]) / 2)
                            + cent
                        )

            # DISCUSS WITH SIMON Can replace with sending new_data through to subpattern call and add functions to
            # class to find sub-arrays - will depend on parallelisation though....

            # Get subsection of data to pass
            sub_pattern = np.where((twotheta >= tth_range[0]) & (twotheta <= tth_range[1]))
            twotheta_sub = twotheta[sub_pattern]
            d_spacing_sub = dspace[sub_pattern]
            azimuth_sub = azimuth[sub_pattern]
            intens_sub = intens[sub_pattern]

            # Mask the subpattern by intensity if called for
            if "imax" in orders or "imin" in orders:
                imax = np.inf
                imin = 0
                if "imax" in orders:
                    imax = int(orders["imax"])
                if "imin" in orders:
                    imin = int(orders["imin"])
                intens_sub = ma.masked_outside(intens_sub, imin, imax)
                azimuth_sub = ma.array(azimuth_sub, mask=intens_sub.mask)
                twotheta_sub = ma.array(twotheta_sub, mask=intens_sub.mask)
                d_spacing_sub = ma.array(d_spacing_sub, mask=intens_sub.mask)

            if mode == "setrange":
                # FIX ME: this plots the input data. perhaps it should have its own switch rather than being 
                # subservient to Debug.
                # It is not a debug it is a setup thing.
                # plot the data and the mask.
                fig_1 = new_data.full_plot(twotheta_sub,
                                           azimuth_sub,
                                           intens_sub,
                                           plot_type="mask",
                                           name=peak_string(orders)
                                           )
                filename = make_outfile_name(
                    diff_files[j],
                    directory=fit_settings.Output_directory,
                    additional_text="mask",
                    orders=orders,
                    extension=".png",
                    overwrite=True,
                )

                fig_1.savefig(filename)

                if debug:
                    plt.draw()
                plt.close()

                step = 1000

            elif mode == "setguess":

                fig_1 = new_data.full_plot(twotheta_sub,
                                           azimuth_sub,
                                           intens_sub,
                                           plot_type="data",
                                           name=peak_string(orders))
                (points,) = fig_1.get_axes()[0].plot([], [], )
                point_builder = PointBuilder(points, fig_1.get_axes()[0])
                plt.show(block=True)

                selection_arr = point_builder.array()
                # print(selection_arr)
                print("")
                print("Selected points, %s:" % peak_string(orders))
                print("[")
                for k in range(len(selection_arr)):
                    print(json.dumps(selection_arr[k]) + ",")
                print("]")
                print("")
                step = 1000
            else:
                if parallel is True:  # setup parallel version
                    kwargs = {
                        "save_fit": save_figs,
                        "debug": debug,
                        "refine": refine,
                        "iterations": iterations,
                        "f_name": diff_files[j],
                        "fdir": fit_settings.Output_directory,
                    }
                    args = (
                        [twotheta_sub, d_spacing_sub, parms_dict],
                        azimuth_sub,
                        intens_sub,
                        new_data,
                        orders,
                        params,
                        bounds,
                    )
                    parallel_pile.append((args, kwargs))

                else:  # non-parallel version
                    tmp = fit_sub_pattern(
                        [twotheta_sub, d_spacing_sub, parms_dict],
                        azimuth_sub,
                        intens_sub,
                        new_data,
                        orders,
                        params,
                        bounds,
                        save_fit=save_figs,
                        debug=debug,
                        refine=refine,
                        iterations=iterations,
                        f_name=diff_files[j],
                        fdir=fit_settings.Output_directory,
                    )
                    fitted_param.append(tmp[0])
                    lmfit_models.append(tmp[1])

        # write output files
        if mode == "fit":
            if parallel is True:
                tmp = p.map(parallel_processing, parallel_pile)
                for i in range(n_sub_patterns):
                    fitted_param.append(tmp[i][0])
                    lmfit_models.append(tmp[i][1])

            # store the fit parameters' information as a JSON file.
            filename = make_outfile_name(
                diff_files[j],
                directory=fit_settings.Output_directory,
                extension=".json",
                overwrite=True,
            )
            with open(filename, "w") as TempFile:
                # Write a JSON string into the file.
                json.dump(
                    fitted_param,
                    TempFile,
                    sort_keys=True,
                    indent=2,
                    default=json_numpy_serializer,
                )

            # if propagating the fits write them to a temporary file
            if propagate:
                # print json_string
                with open(temporary_data_file, "w") as TempFile:
                    # Write a JSON string into the file.
                    json.dump(
                        fitted_param,
                        TempFile,
                        sort_keys=True,
                        indent=2,
                        default=json_numpy_serializer,
                    )

    if mode == "fit":
        # Write the output files.
        write_output(settings_file, fit_settings=fit_settings, fit_parameters=fit_parameters, parms_dict=parms_dict,
                     outtype=fit_settings.Output_type)

    if parallel is True:
        p.close()


def parallel_processing(p):
    a, kw = p
    return fit_sub_pattern(*a, **kw)


if __name__ == "__main__":
    # Load settings fit settings file.
    sys.path.append(os.getcwd())
    settings_file = sys.argv[1]
    print(settings_file)
    execute(inputs=None, debug=False, refine=True, save_all=False, iterations=1, setting_file=settings_file, parallel=False)
