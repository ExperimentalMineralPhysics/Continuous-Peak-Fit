#!/usr/bin/env python

from __future__ import annotations  # Enables additional Python features

__all__ = ["execute", "write_output"]

import json
import logging
import sys
from importlib import import_module
from os import cpu_count
from pathlib import Path
from types import ModuleType
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
import proglog
from pathos.pools import ParallelPool

from cpf import output_formatters
from cpf.BrightSpots import SpotProcess
from cpf.data_preprocess import remove_cosmics as cosmicsimage_preprocess
from cpf.IO_functions import (
    any_terms_null,
    json_numpy_serializer,
    make_outfile_name,
    peak_string,
    title_file_names,
)
from cpf.logging import CPFLogger
from cpf.settings import Settings
from cpf.XRD_FitSubpattern import fit_sub_pattern

np.set_printoptions(threshold=sys.maxsize)
# FIX ME: Need to add complexity here.
# Need option to check required functions are present if adding new output format.
# Also need to do something similar with data class


logger = CPFLogger("cpf.XRD_FitPattern")


def register_default_formats() -> dict[str, ModuleType]:
    """
    Load all available output modules
    :return:
    """
    # FIX ME: We could add extra checks here to make sure the required functions exist in each case.
    output_list = output_formatters.module_list
    new_module = {}
    for output_module in output_list:
        module: ModuleType = import_module(f"cpf.output_formatters.{output_module}")
        new_module[output_module[5:]] = module
    return new_module


# Load potential output formats
output_methods_modules = register_default_formats()


def initialise_logger(
    settings_file: Optional[str | Path] = None,
    report: str | bool = False,
):
    # Start the logger with the desired outputs
    logger.handlers.clear()
    format = "%(asctime)s [%(levelname)s] %(message)s"
    formatter = logging.Formatter(format)
    if isinstance(report, (str)):
        # Fail gracefully
        if settings_file is None:
            raise ValueError(
                "Settings file needs to be specified in order to create a log file."
            )

        # Set the logging level and log file name
        level = report.upper()
        log_name = make_outfile_name(settings_file, extension=".log", overwrite=True)

        # Create and add the file handler
        fh = logging.FileHandler(log_name, mode="a", encoding="utf-8")
        fh.setFormatter(formatter)
        fh.setLevel(level)
        logger.addHandler(fh)
    else:
        level = "INFO"

    # Create the stream handler
    sh = logging.StreamHandler()
    sh.setFormatter(formatter)
    sh.setLevel(level)
    logger.addHandler(sh)


def initiate(
    settings_file: Optional[Union[str, Path]] = None,
    inputs=None,
    out_type=None,
    report: str | bool = False,
    **kwargs,
):
    """
    Run checks on input files, initiate data class and check output options

    :param report:
    :param out_type:
    :param initiate_data:
    :param settings_file:
    :param inputs:
    :return fit_parameters:
    :return fit_settings:
    """

    # Fail gracefully
    if settings_file is None:
        raise ValueError(
            "Either the settings file or the parameter dictionary need to be specified."
        )
    # Convert to Path object
    if isinstance(settings_file, str):
        try:
            settings_file = Path(settings_file)
        except Exception as error:
            raise error

    # make a header in the log file so that we know where the processing starts
    logger.info("")
    logger.info("=================================================================")
    logger.info("")
    logger.info("Starting data proceesing using: %s.py" % settings_file)
    logger.info("")
    logger.info("=================================================================")
    logger.info("")

    # If no params_dict then initiate. Check all the output functions are present and valid.
    settings_for_fit = Settings()
    settings_for_fit.populate(settings_file=settings_file, report=report)

    return settings_for_fit


def view(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug=False,
    refine=True,
    save_all=False,
    # propagate=True,
    iterations=1,
    # track=False,
    parallel=True,
    pattern="all",
    subpattern="all",
    report=False,
):
    """
    :param settings_file:
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

    if settings_class is None:
        settings_for_fit = initiate(settings_file, inputs=inputs, report=report)
    else:
        settings_for_fit = settings_class

    # view the listed file only
    if pattern != "all":
        # restrict file list to first file
        settings_for_fit.set_data_files(keep=pattern)

    execute(
        settings_class=settings_for_fit,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        mode="view",
        report=True,
    )


def set_range(
    settings_file: Optional[Union[str, Path]] = None,
    settings_class: Optional[Settings] = None,
    inputs=None,
    debug: bool = False,
    refine: bool = True,
    save_all: bool = False,
    # propagate: bool = True,
    iterations: int = 1,
    # track: bool = False,
    parallel: bool = True,
    subpattern: str = "all",
    report: bool = False,
):
    """
    :param settings_file:
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

    settings_for_fit: Settings = (
        initiate(settings_file, inputs=inputs, report=report)
        if settings_class is None
        else settings_class
    )

    # search over the first file only
    # restrict file list to first file
    settings_for_fit.set_data_files(keep=0)

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    execute(
        settings_class=settings_for_fit,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        mode="set-range",
        report=True,
    )


def initial_peak_position(
    settings_file: Optional[Union[str, Path]] = None,
    settings_class: Optional[Settings] = None,
    inputs=None,
    debug: bool = False,
    refine: bool = True,
    save_all: bool = False,
    # propagate: bool = True,
    iterations: int = 1,
    # track: bool = False,
    parallel: bool = True,
    subpattern: str = "all",
    report: bool = False,
):
    """
    Calls interactive graph to set the inital peak postion guesses.

    The event handler code is copied from:
    https://matplotlib.org/stable/users/event_handling.html for how to make work

    :param settings_file:
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

    settings_for_fit: Settings = (
        initiate(settings_file, inputs=inputs, report=report)
        if settings_class is None
        else settings_class
    )

    # search over the first file only
    settings_for_fit.set_data_files(keep=0)

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    logger.info("\n'initial_peak_position' needs an interactive matplotlib figure.")
    logger.info(
        " ".join(
            map(
                str,
                [
                    (
                        "If you are using sypder with inline figures, call '%matplotlib qt', then rerun the script"
                    )
                ],
            )
        )
    )
    logger.info(
        " ".join(
            map(
                str,
                [
                    (
                        "To restore the inline plotting afterwards call '%matplotlib inline'"
                    )
                ],
            )
        )
    )
    logger.info(
        " ".join(map(str, [("To move to the next peak selection close the window.\n")]))
    )

    execute(
        settings_class=settings_for_fit,
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
        # logger.info(" ".join(map(str, [('click', event)])))
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
    settings_file: Optional[Union[str, Path]] = None,
    settings_class: Optional[Settings] = None,
    inputs=None,
    debug: bool = False,
    refine: bool = True,
    save_all: bool = False,
    # propagate: bool = True,
    iterations: int = 1,
    # track: bool = False,
    parallel: bool = True,
    search_parameter: str = "height",
    search_over: list[int] = [0, 20],
    subpattern: str = "all",
    search_peak: int = 0,
    search_series: list[str] = ["fourier", "spline"],
    report: bool = False,
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
    :param settings_file:
    :param inputs:
    :param debug:
    :param refine:
    :param iterations:
    :return:
    """

    settings_for_fit: Settings = (
        initiate(settings_file, inputs=inputs, report=report)
        if settings_class is None
        else settings_class
    )

    # force it to write the required output type.
    settings_for_fit.set_output_types(out_type_list=["DifferentialStrain"])

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

    settings_for_fit.file_label = (
        "search="
        + search_parameter
        + "_subpattern="
        + str(subpattern)
        + "_peak="
        + str(search_peak)
    )

    execute(
        settings_class=settings_for_fit,
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
        settings_file=settings_file,
        settings_class=settings_for_fit,
        debug=True,
        out_type="DifferentialStrain",
    )


def write_output(
    settings_file: Optional[Union[str, Path]] = None,
    settings_class: Optional[Settings] = None,
    # fit_settings=None,
    # fit_parameters=None,
    parms_dict=None,
    out_type=None,
    det=None,
    use_bounds: bool = False,
    differential_only: bool = False,
    debug: bool = False,
    report: bool = False,
):
    """
    :param debug:
    :param fit_parameters:
    :param use_bounds:
    :param differential_only:
    :param settings_file:
    :param fit_settings:
    :param parms_dict:
    :param out_type:
    :param det:
    :return:
    """

    settings_class = (
        initiate(settings_file, report=report, out_type=out_type)
        if settings_class is None
        else settings_class
    )

    if out_type is not None:
        logger.moreinfo(  # type: ignore
            " ".join(
                map(
                    str,
                    [
                        (
                            "Output_type was provided as an option; will use %s "
                            % out_type
                        )
                    ],
                )
            )
        )
        settings_class.set_output_types(out_type_list=out_type)

    if settings_class.output_types is None:
        logger.warning(
            " ".join(
                map(
                    str,
                    [
                        (
                            "There are no output types. Add 'Output_type' to input file or specify 'out_type' in command."
                        )
                    ],
                )
            )
        )
    else:
        for mod in settings_class.output_types:
            logger.info(" ".join(map(str, [("Writing output file(s) using %s" % mod)])))
            wr = output_methods_modules[mod]
            wr.WriteOutput(
                settings_class=settings_class,
                settings_file=settings_file,
                differential_only=differential_only,
                debug=debug,
            )


def execute(
    settings_file: Optional[Union[str, Path]] = None,
    settings_class=None,
    # fit_settings=None,
    # fit_parameters=None,
    inputs=None,
    debug: bool = False,
    refine: bool = True,
    save_all: bool = False,
    # propagate: bool = True, #moved this option to settings file
    iterations: int = 1,
    # track: bool = False,  #moved this option to settings file
    parallel: bool = True,
    mode: str = "fit",
    report: bool = False,
    fit_method: str = "leastsq",
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
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

    initialise_logger(settings_file=settings_file, report=report)

    if settings_class is None:
        settings_for_fit = initiate(settings_file, inputs=inputs, report=report)
    else:
        settings_for_fit = settings_class
    new_data = settings_for_fit.data_class

    # Define locally required names
    temporary_data_file = make_outfile_name(
        "PreviousFit_JSON",
        directory=settings_for_fit.output_directory,
        extension=".dat",
        overwrite=True,
    )

    if settings_for_fit.calibration_data:
        data_to_fill = Path(settings_for_fit.calibration_data).resolve()
    else:
        data_to_fill = settings_for_fit.image_list[0]

    new_data.fill_data(
        data_to_fill,
        settings=settings_for_fit,
        debug=debug,
    )

    # Get calibration parameter file
    parms_dict = new_data.calibration
    # FIXME this should be removable.

    # plot calibration file
    if (
        logger.is_below_level(level="DEBUG")
        and settings_for_fit.calibration_data is not None
    ):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        new_data.plot_collected(fig_plot=fig, axis_plot=ax)
        plt.title("Calibration data")
        plt.show()
        plt.close()

    # if parallel processing start the pool
    if parallel is True:
        pool = ParallelPool(nodes=cpu_count())

        # Since we may have already closed the pool, try to restart it
        try:
            pool.restart()
        except AssertionError:
            pass

    # Process the diffraction patterns
    # for j in range(settings_for_fit.image_number):
    progress = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger
    for j in progress.iter_bar(iteration=range(settings_for_fit.image_number)):
        logger.info(
            " ".join(
                map(
                    str,
                    [
                        (
                            "Processing %s"
                            % title_file_names(
                                image_name=settings_for_fit.image_list[j]
                            )
                        )
                    ],
                )
            )
        )

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
        if logger.is_below_level(level="DEBUG") or mode == "view":
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity")
            plt.title(title_file_names(settings_for_fit=settings_for_fit, num=j))
            plt.show()
            if mode == "view":
                sys.exit("Plotted data.")
            else:
                plt.close()

        # Get previous fit (if it exists and is required)
        if (
            Path(temporary_data_file).is_file()
            and settings_for_fit.fit_propagate is True
            and mode == "fit"
        ):
            # Read JSON data from file
            logger.moreinfo(  # type: ignore
                " ".join(
                    map(
                        str,
                        [
                            (
                                "Loading previous fit results from %s"
                                % temporary_data_file
                            )
                        ],
                    )
                )
            )
            with open(temporary_data_file) as json_data:
                previous_fit = json.load(json_data)

                # if the previous_fit is not the same size as fit_orders the inout file must have been changed.
                # so discard the previous fit and start again.
                if len(previous_fit) != len(settings_for_fit.fit_orders):
                    del previous_fit

        # Switch to save the first fit in each sequence.
        save_figs = True if (j == 0 or save_all is True) else False

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
            tth_range = np.array(settings_for_fit.subfit_orders["range"])
            if settings_for_fit.fit_track is True and "previous_fit" in locals():
                clean = any_terms_null(params, val_to_find=None)
                if clean == 0:
                    # the previous fit has problems so discard it
                    logger.moreinfo(  # type: ignore
                        " ".join(
                            map(
                                str,
                                [
                                    (
                                        "Tracking peak centre but propagated fit has problems. Not sensible to track the centre of the fit for this step."
                                    )
                                ],
                            )
                        )
                    )
                    params = []
                else:
                    mid = []
                    for k in range(len(params["peak"])):
                        mid.append(params["peak"][k]["d-space"][0])
                        # FIXME: replace with caluculation of mean d-spacing.

                    cent = new_data.conversion(np.mean(mid), reverse=True)

                    move_by = cent - np.mean(tth_range)
                    move_by = move_by[
                        0
                    ]  # this is needed to turn move_by from array to float

                    # update tth_range and settings
                    tth_range = tth_range + move_by
                    settings_for_fit.fit_orders[i]["range"] = (
                        settings_for_fit.fit_orders[i]["range"] + move_by
                    )

                    logger.moreinfo(  # type: ignore
                        " ".join(
                            map(
                                str,
                                [
                                    (
                                        f"Move range for fitting. \n Initial range: [{tth_range[0]:4.2f},{tth_range[1]:4.2f}]; will be moved by {move_by:4.2f}; the new range is [{tth_range[0]+move_by:4.2f},{tth_range[1]+move_by:4.2f}]"
                                    )
                                ],
                            )
                        )
                    )

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
                            )

                    # re-get settings for current subpattern
                    settings_for_fit.set_subpattern(j, i)

            sub_data = new_data.duplicate(range_bounds=tth_range)
            # sub_data.set_limits(range_bounds=tth_range)

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
                # plt.close()

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

                # report the points selected.
                # set to critical to ensure that they are printed -- because HAS to be done.
                logger.critical(
                    " ".join(
                        map(
                            str,
                            [
                                (
                                    "Selected points for %s peak(s): ["
                                    % peak_string(settings_for_fit.subfit_orders)
                                )
                            ],
                        )
                    )
                )
                for k in range(len(selection_arr)):
                    logger.critical(
                        " ".join(map(str, [(json.dumps(selection_arr[k]) + ",")]))
                    )
                logger.critical(" ".join(map(str, [("]")])))

            else:
                if parallel is True:  # setup parallel version
                    kwargs = {
                        "previous_params": params,
                        "save_fit": save_figs,
                        "debug": debug,
                        "refine": refine,
                        "iterations": iterations,
                        "min_data_intensity": settings_for_fit.fit_min_data_intensity,
                        "min_peak_intensity": settings_for_fit.fit_min_peak_intensity,
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
                        min_data_intensity=settings_for_fit.fit_min_data_intensity,
                        min_peak_intensity=settings_for_fit.fit_min_peak_intensity,
                        fit_method=fit_method,
                    )
                    fitted_param.append(tmp[0])
                    lmfit_models.append(tmp[1])

        # write output files
        if mode == "fit" or mode == "search":
            if parallel is True:
                tmp = pool.map(parallel_processing, parallel_pile)
                for i in range(len(settings_for_fit.fit_orders)):
                    fitted_param.append(tmp[i][0])
                    lmfit_models.append(tmp[i][1])

            # store the fit parameters' information as a JSON file.
            if mode == "search":
                additional_text = settings_for_fit.file_label
            else:
                additional_text = None

            filename = make_outfile_name(
                settings_for_fit.subfit_filename,
                directory=settings_for_fit.output_directory,
                additional_text=additional_text,
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
            settings_file=settings_file, settings_class=settings_for_fit, debug=debug
        )

    if parallel is True:
        pool.clear()


def parallel_processing(p):
    a, kw = p
    return fit_sub_pattern(*a, **kw)


if __name__ == "__main__":
    # Load settings fit settings file.
    sys.path.append(str(Path().cwd()))
    settings_file = sys.argv[1]
    logger.info(" ".join(map(str, [(settings_file)])))
    # Safely exit the program
    for handler in logging.getLogger().handlers:
        handler.flush()
    # sys.exit()
    execute(
        inputs=None,
        debug=False,
        refine=True,
        save_all=False,
        iterations=1,
        settings_file=settings_file,
        parallel=False,
    )
