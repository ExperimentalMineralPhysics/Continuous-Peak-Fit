#!/usr/bin/env python

__all__ = ["execute", "write_output"]

import json
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import pathos.multiprocessing as mp
from cpf.settings import settings
from cpf.XRD_FitSubpattern import fit_sub_pattern
from cpf.Cosmics import image_preprocess as cosmicsimage_preprocess
from cpf.IO_functions import (
    json_numpy_serializer,
    make_outfile_name,
    peak_string,
    any_terms_null,
    title_file_names,
)

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
        module = __import__(
            "cpf.output_formatters.Write" + output_module, fromlist=[None]
        )
        new_module[output_module] = module
    return new_module


# Load potential output formats
output_methods_modules = register_default_formats()


def initiate(setting_file=None, inputs=None, out_type=None, report=False, **kwargs):
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

    # Fail gracefully
    if setting_file is None:
        raise ValueError(
            "Either the settings file or the parameter dictionary need to be specified."
        )
    # If no params_dict then initiate. Check all the output functions are present and valid.

    settings_for_fit = settings()
    settings_for_fit.populate(settings_file=setting_file, report=report)

    return settings_for_fit


def set_range(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    # propagate=True,
    iterations=1,
    # track=False,
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

    if setting_class is None:
        settings_for_fit = initiate(setting_file, inputs=inputs, report=True)
    else:
        settings_for_fit = setting_class

    # search over the first file only
    # restrict file list to first file
    settings_for_fit.set_data_files(keep=0)

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    execute(
        setting_class=settings_for_fit,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        mode="set-range",
        report=True,
    )


def initial_peak_position(
    setting_file=None,
    setting_class=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    # propagate=True,
    iterations=1,
    # track=False,
    parallel=True,
    subpattern="all",
    **kwargs
):
    """
    Calls interactive graph to set the inital peak postion guesses. 
    
    The event handler code is copied from: 
    https://matplotlib.org/stable/users/event_handling.html for how to make work
    
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

    if setting_class is None:
        settings_for_fit = initiate(setting_file, inputs=inputs, report=True)
    else:
        settings_for_fit = setting_class

    # search over the first file only
    settings_for_fit.set_data_files(keep=0)

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    print("\n'initial_peak_position' needs an interactive matplotlib figure.")
    print(
        "If you are using sypder with inline figures, call '%matplotlib qt', then rerun the script"
    )
    print("To restore the inline plotting afterwards call '%matplotlib inline'")
    print("To move to the next peak selection close the window.\n")

    execute(
        setting_class=settings_for_fit,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        mode="set-guess",
        report=True,
    )


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
    setting_class= None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    # propagate=True,
    iterations=1,
    # track=False,
    parallel=True,
    search_parameter="height",
    search_over=[0, 20],
    subpattern="all",
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

    if setting_class is None:
        settings_for_fit = initiate(setting_file, inputs=inputs, report=True)
    else:
        settings_for_fit = setting_class

    # force it to write the required output type.
    settings_for_fit.set_output_types(out_type_list="DifferentialStrain")

    # search over the first file only
    settings_for_fit.set_data_files(keep=0)

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    # set search orders
    settings_for_fit.set_order_search(
        search_parameter=search_parameter,
        search_over=search_over,
        subpatterns=subpattern,
        search_peak=search_peak,
        search_series=search_series,
    )
    
    settings_for_fit.file_label = "search="+search_parameter+"_subpattern="+str(subpattern)+"_peak="+str(search_peak)

    execute(
        setting_class=settings_for_fit,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        report=True,
    )

    # write a differential strain output file
    # FIX ME: using debug as a switch to get all info in file.
    # should probably use another switch
    write_output(
        setting_file=setting_file,
        setting_class=settings_for_fit,
        debug=True,
        out_type="DifferentialStrain",
    )


def write_output(
    setting_file=None,
    setting_class=None,
    # fit_settings=None,
    # fit_parameters=None,
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

    if setting_class is None:
        setting_class = initiate(setting_file, report=True, out_type=out_type)

    if out_type is not None:
        print("Changing output_type to " + out_type)
        setting_class.set_output_types(out_type_list=out_type)

    if setting_class.output_types is None:
        print("Thre are no output types. Add 'Output_type' to input file or specify 'out_type' in command.")
    else:
        for mod in setting_class.output_types:
            print("\nWrite output file(s) using", mod)
            wr = output_methods_modules[mod]
            wr.WriteOutput(
                setting_class=setting_class,
                setting_file=setting_file,
                differential_only=differential_only,
                debug=debug,
            )


def execute(
    setting_file=None,
    setting_class=None,
    # fit_settings=None,
    # fit_parameters=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    # propagate=True, #moved this option to settings file
    iterations=1,
    # track=False,  #moved this option to settings file
    parallel=True,
    mode="fit",
    report=False,
    intensity_threshold = 0,
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

    if setting_class is None:
        settings_for_fit = initiate(setting_file, inputs=inputs, report=True)
    else:
        settings_for_fit = setting_class
    new_data = settings_for_fit.data_class

    # Define locally required names
    temporary_data_file = make_outfile_name(
        "PreviousFit_JSON",
        directory=settings_for_fit.output_directory,
        extension=".dat",
        overwrite=True,
    )

    if settings_for_fit.calibration_data:
        data_to_fill = os.path.abspath(settings_for_fit.calibration_data)
    else:
        # data_to_fill = os.path.abspath(settings_for_fit.datafile_list[0])
        data_to_fill = settings_for_fit.image_list[0]
    print(data_to_fill)
    new_data.fill_data(
        data_to_fill,
        settings=settings_for_fit,
        debug=debug,
    )

    # Get calibration parameter file
    parms_dict = new_data.calibration
    # FIXME this should be removable.

    # plot calibration file
    if debug and settings_for_fit.calibration_data is not None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        new_data.plot_collected(fig_plot=fig, axis_plot=ax)
        plt.title("Calibration data")
        plt.show()
        plt.close()

    # if parallel processing start the pool
    if parallel is True:
        # p = mp.Pool(processes=mp.cpu_count())
        p = mp.ProcessPool(nodes=mp.cpu_count())
        # p = mp.Pool()

        # Since we may have already closed the pool, try to restart it
        try:
            p.restart()
        except AssertionError:
            pass

    # Process the diffraction patterns #
    for j in range(settings_for_fit.image_number):

        # print("Process ", settings_for_fit.image_list[j])
        print("Process ", title_file_names(image_name=settings_for_fit.image_list[j]))

        # Get diffraction pattern to process.
        new_data.import_image(settings_for_fit.image_list[j], debug=debug)

        if settings_for_fit.datafile_preprocess is not None:
            # needed because image preprocessing adds to the mask and is different for each image.
            new_data.mask_restore()
            if "cosmics" in settings_for_fit.datafile_preprocess:
                new_data = cosmicsimage_preprocess(new_data, settings_for_fit)
        else:
            # nothing is done here.
            pass

        # plot input file
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity")
            # plt.title(os.path.basename(settings_for_fit.datafile_list[j]))
            plt.title(title_file_names(settings_for_fit=settings_for_fit, num=j))
            plt.show()
            plt.close()

        # Get previous fit (if it exists and is required)
        if (
            os.path.isfile(temporary_data_file)
            and settings_for_fit.fit_propagate is True
            and mode == "fit"
        ):
            # Read JSON data from file
            print("Loading previous fit results from %s" % temporary_data_file)
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)
                
                # if the previous_fit is not the same size as fit_orders the inout file must have been changed. 
                # so discard the previous fit and start again.
                if len(previous_fit) != len(settings_for_fit.fit_orders):
                    del previous_fit

        # Switch to save the first fit in each sequence.
        if j == 0 or save_all is True:
            save_figs = 1
        else:
            save_figs = 0

        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        fitted_param = []
        lmfit_models = []
        parallel_pile = []

        for i in range(len(settings_for_fit.fit_orders)):

            # get settings for current subpattern
            settings_for_fit.set_subpattern(j, i)

            if "previous_fit" in locals() and mode == "fit":
                params = previous_fit[i]
            else:
                params = []

            # Track the position of the peak centroid
            # FIXME: This is crude - the range doesn't change width. so can't account for massive change in stress.
            # But does it need to?
            tth_range = settings_for_fit.subfit_orders["range"]
            if settings_for_fit.fit_track is True and "previous_fit" in locals():
                clean = any_terms_null(params, val_to_find=None)
                if clean == 0:
                    # the previous fit has problems so discard it
                    print(
                        "Propagated fit has problems so not sensible to track the centre of the fit."
                    )
                    params = []
                else:
                    print("old subpattern range", tth_range)
                    mid = []
                    for k in range(len(params["peak"])):
                        mid.append(params["peak"][k]["d-space"][0])
                        # FIXME: replace with caluculation of mean d-spacing.

                    # Replace this with call through to class structure?
                    # print(mid)
                    cent = new_data.conversion(np.mean(mid), reverse=True)
                    # print("cent", cent)
                    # print("mean tth_range", np.mean(tth_range))
                    # print((tth_range[1] + tth_range[0]) / 2)
                    move_by = cent - np.mean(tth_range)
                    print("move by", move_by)
                    tth_range = tth_range + move_by
                    # tth_range[0] = tth_range[0] - ((tth_range[1] + tth_range[0]) / 2) + cent
                    # tth_range[1] = tth_range[1] - ((tth_range[1] + tth_range[0]) / 2) + cent
                    # print("move by:", ((tth_range[1] + tth_range[0]) / 2) - cent)
                    print("new subpattern range", tth_range)

                    # copy new positions back into settings_for_fit
                    settings_for_fit.fit_orders[i]["range"] = (
                        settings_for_fit.fit_orders[i]["range"] + move_by
                    )
                    print(settings_for_fit.fit_orders[i]["range"])

                    # The PeakPositionSelections are only used if the fits are not being propagated
                    if "PeakPositionSelection" in settings_for_fit.fit_orders[i]:
                        for k in range(
                            len(settings_for_fit.fit_orders[i]["PeakPositionSelection"])
                        ):
                            settings_for_fit.fit_orders[i]["PeakPositionSelection"][k][
                                2
                            ] = (
                                settings_for_fit.fit_orders[i]["PeakPositionSelection"][
                                    k
                                ][2]
                                + move_by
                                # - ((tth_range[1] + tth_range[0]) / 2)
                                # + cent
                            )

                    # re-get settings for current subpattern
                    settings_for_fit.set_subpattern(j, i)

            sub_data = new_data.duplicate()
            sub_data.set_limits(range_bounds=tth_range)

            # Mask the subpattern by intensity if called for
            if (
                "imax" in settings_for_fit.subfit_orders
                or "imin" in settings_for_fit.subfit_orders
            ):
                sub_data = SpotProcess(sub_data, settings_for_fit)
            if mode == "set-range":

                fig_1 = plt.figure()
                sub_data.plot_masked(fig_plot=fig_1)
                plt.suptitle(peak_string(settings_for_fit.subfit_orders) + "; masking")

                # filename = make_outfile_name(
                #     settings_for_fit.datafile_list[j],
                #     directory=settings_for_fit.output_directory,
                #     additional_text="mask",
                #     orders=settings_for_fit.subfit_orders,
                #     extension=".png",
                #     overwrite=True,
                # )
                filename = make_outfile_name(
                    settings_for_fit.image_list[j],
                    directory=settings_for_fit.output_directory,
                    additional_text="mask",
                    orders=settings_for_fit.subfit_orders,
                    extension=".png",
                    overwrite=True,
                )

                fig_1.savefig(filename)

                # if debug:
                plt.show()
                plt.close()

            elif mode == "set-guess":

                fig_1 = plt.figure()
                ax = fig_1.add_subplot(1, 1, 1)
                ax_o1 = plt.subplot(111)
                sub_data.plot_calibrated(
                    fig_plot=fig_1, axis_plot=ax, y_axis="azimuth", limits=[0, 100]
                )
                plt.title(peak_string(settings_for_fit.subfit_orders))

                (points,) = fig_1.get_axes()[0].plot(
                    [],
                    [],
                )
                point_builder = PointBuilder(points, fig_1.get_axes()[0])
                plt.show(block=True)

                selection_arr = point_builder.array()
                # print(selection_arr)
                print("")
                print(
                    "Selected points, %s:" % peak_string(settings_for_fit.subfit_orders)
                )
                print("[")
                for k in range(len(selection_arr)):
                    print(json.dumps(selection_arr[k]) + ",")
                print("]")
                print("")

            else:
                if parallel is True:  # setup parallel version
                    kwargs = {
                        "previous_params": params,
                        "save_fit": save_figs,
                        "debug": debug,
                        "refine": refine,
                        "iterations": iterations,
                        "intensity_threshold":intensity_threshold,
                    }
                    arg = (sub_data, settings_for_fit.duplicate())
                    parallel_pile.append((arg, kwargs))

                else:  # non-parallel version
                    tmp = fit_sub_pattern(
                        sub_data,
                        settings_for_fit,  # added
                        params,
                        save_fit=save_figs,
                        debug=debug,
                        refine=refine,
                        iterations=iterations,
                        intensity_threshold=intensity_threshold,
                    )
                    fitted_param.append(tmp[0])
                    lmfit_models.append(tmp[1])

        # write output files
        if mode == "fit":
            if parallel is True:
                tmp = p.map(parallel_processing, parallel_pile)
                for i in range(len(settings_for_fit.fit_orders)):
                    fitted_param.append(tmp[i][0])
                    lmfit_models.append(tmp[i][1])

            # store the fit parameters' information as a JSON file.
            filename = make_outfile_name(
                settings_for_fit.subfit_filename,
                directory=settings_for_fit.output_directory,
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
            if settings_for_fit.fit_propagate:
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
        write_output(
            setting_file=setting_file, setting_class=settings_for_fit, debug=True
        )

    if parallel is True:
        p.clear()


def parallel_processing(p):
    a, kw = p
    return fit_sub_pattern(*a, **kw)


if __name__ == "__main__":
    # Load settings fit settings file.
    sys.path.append(os.getcwd())
    settings_file = sys.argv[1]
    print(settings_file)
    execute(
        inputs=None,
        debug=False,
        refine=True,
        save_all=False,
        iterations=1,
        setting_file=settings_file,
        parallel=False,
    )
