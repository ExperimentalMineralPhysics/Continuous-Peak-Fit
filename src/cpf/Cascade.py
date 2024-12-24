#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 14:16:58 2021

@author: simon


This files makes the cascade plots of azimuthal intensity vis time.
It needs to do a few things.

1. Call XRD_FitPattern to make the azimual distribution files.
    a. by fitting a peak to the data
    b. max
    c. mean of the 90th+ percentiles.

    To be consistent with the philosophy of the code, the bins for the aximuths should
    be determined by a combination of width and the number of data contained. By
    keeping the number of data constant we will have more points at large two theta
    than at small two theta.
    - this will need a switch in the input file somewhere. (Azibins needs another call
      e.g. AziNum?)

2. Plot the outputs as a cascade plot.

3. Determine the number of peaks in each steip and plot n vs. time.
    Peaks need to be determinable by:
        a. above or below average (e.g. He 2000)
        b. peak finding algorithm (various)


# Another algorithm that could be tried is:
https://scikit-image.org/docs/stable/auto_examples/features_detection/plot_blob.html
https://scikit-image.org/docs/stable/api/skimage.feature.html#skimage.feature.blob_log
"""

from __future__ import annotations

__all__ = ["initiate", "set_range", "execute"]

import json
import logging
from os import cpu_count
from pathlib import Path
from typing import Literal, Optional

import matplotlib.colors as colours
import matplotlib.pyplot as plt
import numpy as np
import proglog
from pathos.pools import ParallelPool

# import plotly.graph_objects as go
# import pandas as pd
from scipy.signal import find_peaks

import cpf.IO_functions as IO
import cpf.XRD_FitPattern as XRD_FitPattern
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

logger = CPFLogger("cpf.Cascade")


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


def initiate(*args, **kwargs):
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

    # pass everything through to XRD_FitPattern.initiate
    settings_for_fit = XRD_FitPattern.initiate(*args, **kwargs)

    return settings_for_fit


def set_range(*args, **kwargs):
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

    # pass everything through to XRD_FitPattern.set_range
    XRD_FitPattern.set_range(*args, **kwargs)


def execute(
    settings_file: Optional[Path] = None,
    settings_class: Optional[Settings] = None,
    inputs=None,
    debug: bool = False,
    save_all: bool = False,
    parallel: bool = True,
    subpattern: str = "all",
    mode: str = "cascade",
    report: bool = False,
    show_plots: bool = False,
    **kwargs,
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
        settings_for_fit = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_for_fit = settings_class
    new_data = settings_for_fit.data_class

    # Define locally required names
    temporary_data_file = make_outfile_name(
        "PreviousChunkFit_JSON",
        directory=settings_for_fit.output_directory,
        extension=".dat",
        overwrite=True,
    )

    data_to_fill: Path
    if settings_for_fit.calibration_data:
        data_to_fill = settings_for_fit.calibration_data.resolve()
    else:
        data_to_fill = settings_for_fit.image_list[0]
        if not isinstance(data_to_fill, Path):
            data_to_fill = Path(data_to_fill)
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
        pool = ParallelPool(nodes=cpu_count())

        # Since we may have already closed the pool, try to restart it
        try:
            pool.restart()
        except AssertionError:
            pass

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

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
                            "Process %s"
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
        if debug:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity")
            plt.title(title_file_names(image_name=settings_for_fit.image_list[j]))
            plt.show()
            plt.close()

        # Get previous fit (if it exists and is required)
        if (
            Path(temporary_data_file).is_file()
            and settings_for_fit.fit_propagate is True
        ):  # and mode == "fit":
            # Read JSON data from file
            logger.info(
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

        # Switch to save the first fit in each sequence.
        if j == 0 or save_all is True:
            save_figs = 1
        else:
            save_figs = 0

        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        all_fitted_chunks = []
        all_chunk_positions = []
        parallel_pile = []

        for i in range(len(settings_for_fit.fit_orders)):
            # get settings for current subpattern
            settings_for_fit.set_subpattern(j, i)

            if "previous_fit" in locals():  # and mode == "fit":
                params = previous_fit
            else:
                params = []

            # Track the position of the peak centroid
            # FIXME: This is crude - the range doesn't change width. so can't account
            # for massive change in stress.
            # But does it need to?
            tth_range = settings_for_fit.subfit_orders["range"]
            if settings_for_fit.cascade_track is True and "previous_fit" in locals():
                clean = any_terms_null(params, val_to_find=None)
                if clean == 0:
                    # the previous fit has problems so discard it
                    logger.info(
                        " ".join(
                            map(
                                str,
                                [
                                    (
                                        "Propagated fit has problems so not sesible to track the centre of the fit."
                                    )
                                ],
                            )
                        )
                    )
                else:
                    logger.info(
                        " ".join(map(str, [("old subpattern range", tth_range)]))
                    )
                    mid = []
                    for k in range(len(params["peak"])):
                        mid.append(params["peak"][k]["d-space"][0])
                        # FIXME: replace with caluculation of mean d-spacing.

                    # Replace this with call through to class structure?
                    # logger.info(" ".join(map(str, [(mid)])))
                    cent = new_data.conversion(np.mean(mid), reverse=True)
                    # logger.info(" ".join(map(str, [("cent", cent)])))
                    # logger.info(" ".join(map(str, [("mean tth_range", np.mean(tth_range))])))
                    # logger.info(" ".join(map(str, [((tth_range[1] + tth_range[0]) / 2)])))
                    move_by = cent - np.mean(tth_range)
                    logger.info(" ".join(map(str, [("move by", move_by)])))
                    tth_range = tth_range + move_by
                    # tth_range[0] = tth_range[0] - ((tth_range[1] + tth_range[0]) / 2) + cent
                    # tth_range[1] = tth_range[1] - ((tth_range[1] + tth_range[0]) / 2) + cent
                    # logger.info(" ".join(map(str, [("move by:", ((tth_range[1] + tth_range[0]) / 2) - cent)])))
                    logger.info(
                        " ".join(map(str, [("new subpattern range", tth_range)]))
                    )

                    # copy new positions back into settings_for_fit
                    settings_for_fit.fit_orders[i]["range"] = (
                        settings_for_fit.fit_orders[i]["range"] + move_by
                    )
                    logger.info(
                        " ".join(map(str, [(settings_for_fit.fit_orders[i]["range"])]))
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

                filename = make_outfile_name(
                    settings_for_fit.datafile_list[j],
                    directory=settings_for_fit.output_directory,
                    additional_text="mask",
                    orders=settings_for_fit.subfit_orders,
                    extension=".png",
                    overwrite=True,
                )

                fig_1.savefig(filename)

                if debug:
                    plt.show()
                plt.close()

            else:
                if parallel is True:  # setup parallel version
                    kwargs = {
                        "save_fit": save_figs,
                        "debug": debug,
                        "mode": mode,
                        "histogram_type": settings_for_fit.cascade_histogram_type,
                        "histogram_bins": settings_for_fit.cascade_histogram_bins,
                    }
                    arg = (sub_data, settings_for_fit.duplicate())
                    parallel_pile.append((arg, kwargs))

                else:  # non-parallel version
                    tmp = fit_sub_pattern(
                        sub_data,
                        settings_for_fit,  # added
                        None,  # do not pass params they are not needed.
                        save_fit=save_figs,
                        debug=debug,
                        mode=mode,
                        histogram_type=settings_for_fit.cascade_histogram_type,
                        histogram_bins=settings_for_fit.cascade_histogram_bins,
                    )
                    all_fitted_chunks.append(tmp[0])
                    all_chunk_positions.append(tmp[1])

        # write output files
        if mode != "set-range":
            if parallel is True:
                tmp = pool.map(parallel_processing, parallel_pile)
                for i in range(len(settings_for_fit.fit_orders)):
                    # fitted_param.append(tmp[i][0])
                    # lmfit_models.append(tmp[i][1])
                    all_fitted_chunks.append(tmp[i][0])
                    all_chunk_positions.append(tmp[i][1])

                # tmp = fit_sub_pattern(
                #     sub_data,
                #     settings_for_fit,  # added
                #     None,  # do not pass params they are not needed.
                #     save_fit=save_figs,
                #     debug=debug,
                #     mode=mode,
                #     histogram_type = settings_for_fit.cascade_histogram_type,
                #     histogram_bins = settings_for_fit.cascade_histogram_bins,
                # )
                # all_fitted_chunks.append(tmp[0])
                # all_chunk_positions.append(tmp[1])

            # paste the fits to an output file.
            # store the chunk fits' information as a JSON file.
            filename = make_outfile_name(
                settings_for_fit.subfit_filename,
                directory=settings_for_fit.output_directory,
                additional_text="chunks",
                # orders=settings_for_fit.subfit_orders,
                extension=".json",
                overwrite=True,
            )
            with open(filename, "w") as TempFile:
                # Write a JSON string into the file.
                json.dump(
                    (all_fitted_chunks, all_chunk_positions),
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
                        (all_fitted_chunks, all_chunk_positions),
                        TempFile,
                        sort_keys=True,
                        indent=2,
                        default=json_numpy_serializer,
                    )

    if parallel is True:
        pool.close()

    # plot the fits
    plot_cascade_chunks(
        inputs=settings_for_fit, report=report, subpattern="all", **kwargs
    )

    plot_peak_count(inputs=settings_for_fit, report=report, subpattern="all", **kwargs)


def parallel_processing(p):
    a, kw = p
    return fit_sub_pattern(*a, **kw)


def read_saved_chunks(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug=False,
    report=False,
    **kwargs,
):
    """


    Parameters
    ----------
    settings_file : TYPE, optional
        DESCRIPTION. The default is None.
    settings_class : TYPE, optional
        DESCRIPTION. The default is None.
    inputs : TYPE, optional
        DESCRIPTION. The default is None.
    debug : TYPE, optional
        DESCRIPTION. The default is False.
    report : TYPE, optional
        DESCRIPTION. The default is False.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    all_azis : TYPE
        DESCRIPTION.
    all_heights : TYPE
        DESCRIPTION.

    """

    if inputs:
        settings_class = inputs
    elif settings_class is None:
        settings_class = initiate(settings_file, report=True)
    else:
        settings_class = settings_class

    all_azis = []
    all_fits = []

    print("Reading chunk files")
    logger = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger

    print("Reading chunk files")
    logger = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger

    print("Reading chunk files")
    logger = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger

    # for f in range(settings_class.image_number):
    for f in logger.iter_bar(iteration=range(settings_class.image_number)):
        settings_class.set_subpattern(f, 0)
        filename = IO.make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            additional_text="chunks",
            extension=".json",
            overwrite=True,
        )
        if Path(filename).is_file():
            # Read JSON data from file
            with open(filename) as json_data:
                fit = json.load(json_data)
                all_azis.append(fit[1])
                all_fits.append(fit[0])

    print("Finished reading chunk files")
    return all_azis, all_fits


def get_chunks_range(chunks, series="h"):
    # make length of g
    min_all = np.inf * np.ones(len(chunks[0]))
    max_all = -np.inf * np.ones(len(chunks[0]))

    for e in range(len(chunks)):
        for g in range(len(chunks[e])):
            for i in range(len(chunks[e][g][series])):
                v_min = np.min(chunks[e][g][series][i])
                v_max = np.max(chunks[e][g][series][i])
                min_all[g] = np.min([min_all[g], v_min])
                max_all[g] = np.max([max_all[g], v_max])

    return min_all, max_all


def plot_cascade_chunks(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug=False,
    report=False,
    plot_type="timeseries",
    subpattern="all",
    scale="linear",
    azi_range="all",
    vmax=np.inf,
    vmin=0,
    show_plots: bool = False,
    **kwargs,
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    if inputs:
        settings_class = inputs
    elif settings_class is None:
        settings_class = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_class = settings_class

    # get file times
    modified_time_s = np.array(
        [Path(file).stat().st_mtime for file in settings_class.image_list]
    )
    modified_time_s -= float(modified_time_s[0])

    y_label_str = r"Time (s)"
    # use file numbers if all times are the same
    if len(np.unique(modified_time_s)) == 1:
        modified_time_s = list(range(settings_class.image_number))
        y_label_str = r"Image in sequence"

    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_plots = len(settings_class.fit_orders)
    else:
        num_plots = len(subpattern)

    all_azis, all_data = read_saved_chunks(
        inputs=settings_class, debug=debug, report=report, subpattern=subpattern
    )
    min_all, max_all = get_chunks_range(all_data, series="h")

    for j in range(num_plots):
        # loop over the number of sets of peaks fit for (i.e. len(settings_class.fit_orders))

        settings_class.set_subpattern(0, j)

        # set plot limits
        if vmax == np.inf:
            vmax = max_all[j]
        # set colour bar ends
        if vmax < max_all[j] and vmin > min_all[j]:
            cb_extend = "both"
        elif vmax < max_all[j]:
            cb_extend = "max"
        elif vmin > min_all[j]:
            cb_extend = "min"
        else:
            cb_extend = "neither"
        # set colour scale normalisation
        norm = None
        if scale == "sqrt":
            norm = colours.PowerNorm(gamma=0.5)
        elif scale == "log":
            norm = colours.LogNorm(vmin=vmin)
        elif scale == "linear":
            norm = None
        for k in range(len(settings_class.subfit_orders["peak"])):
            if plot_type == "timeseries":
                # loop over the number of peaks in each fit_orders
                print(
                    "Making cascade plot for "
                    + IO.peak_string(settings_class.fit_orders[j], peak=k)
                )
                if all_data[0][j]["h"][k]:
                    # if there is some data in the array plot it.
                    fig, ax = plt.subplots()
                    logger = proglog.default_bar_logger(
                        "bar"
                    )  # shorthand to generate a bar logger
                    for i in logger.iter_bar(
                        iteration=range(settings_class.image_number)
                    ):
                        plt.scatter(
                            all_data[i][j]["chunks"],
                            modified_time_s[i]
                            * np.ones(np.shape(all_data[i][j]["chunks"])),
                            s=0.05,
                            c=(all_data[i][j]["h"][k]),
                            # vmax= vmax,
                            # vmin= vmin,
                            cmap="YlOrBr",
                            norm=norm,
                        )

                    # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                    pk: int | Literal["all"] = k
                    for l in range(len(settings_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = IO.peak_string(settings_class.fit_orders[j], peak=pk)
                    plt.title(ttlstr)
                    plt.xlabel(r"Azimuth (deg)")
                    plt.ylabel(y_label_str)

                    if azi_range == "all":
                        azi_range = [
                            np.min(all_data[i][j]["chunks"]),
                            np.max(all_data[i][j]["chunks"]),
                        ]
                    x_ticks = settings_class.data_class.dispersion_ticks(
                        disp_lims=azi_range
                    )
                    ax.set_xticks(x_ticks)

                    cb = plt.colorbar(extend=cb_extend)
                    # cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")

                    # Save the figure
                    filename = IO.make_outfile_name(
                        settings_class.datafile_basename,
                        directory=settings_class.output_directory,
                        additional_text="CascadePlot",
                        orders=settings_class.subfit_orders,
                        peak=pk,
                        extension=".png",
                        overwrite=True,
                    )
                    fig.savefig(filename, transparent=True, bbox_inches="tight")
                    if show_plots is True:
                        plt.show()
                    else:
                        plt.close()

            elif plot_type == "map":
                leng = 500

                pass
            elif plot_type == "map_video":
                pass
            elif plot_type == "data_cube":
                # plot data cube.
                # interactive 3D graph. (i.e. not spyder inline)
                # x, z are position in space
                # y is azimuth
                # plot points for positions with peaks greater than vmin
                # coloured by azimuth
                # size is size of symbol (and limited by vmin, vmax and normalisation)
                """
                print("Making data cube plot for "+ IO.peak_string(settings_class.fit_orders[j], peak=k))
                if all_data[0][j]["h"][k]:
                    # if there is some data in the array plot it.
                    fig, ax = plt.subplots()
                    logger = proglog.default_bar_logger('bar')  # shorthand to generate a bar logger
                    for i in logger.iter_bar(iteration=range(settings_class.image_number)):
                        points_plot = all_data[i][j]["h"][k] >= vmin
                        x = all_data[i][j]["chunks"][points_plot]
                        x_label = r"Azimuth (deg)"
                        z = modified_time_s[i] * np.ones(np.shape(all_data[i][j]["chunks"]))[points_plot]
                        z_label = y_label_str
                        y = all_data[i][j]["h"][k][points_plot]
                        y_label =
                        plt.scatter3(
                            x, y, z,
                            s=1,
                            c=(all_data[i][j]["h"][k]),
                            vmax= vmax,
                            vmin= vmin,
                            cmap="YlOrBr",
                            norm=norm
                        )
                    # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                    pk = k
                    for l in range(len(settings_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = IO.peak_string(settings_class.fit_orders[j], peak=pk)
                    plt.title(ttlstr)
                    plt.xlabel(x_label)
                    plt.ylabel(y_label)
                    cb = plt.colorbar(extend=cb_extend)
                    cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")
                    plt.show()
                    # save the figure
                    filename = IO.make_outfile_name(
                        settings_class.datafile_basename,
                        directory=settings_class.output_directory,
                        additional_text="CascadePlot",
                        orders=settings_class.subfit_orders,
                        peak=pk,
                        extension=".png",
                        overwrite=True,
                    )
                    fig.savefig(filename, transparent=True, bbox_inches="tight")
                print("\n")
                """
                pass

            elif plot_type == "data_cube_video":
                pass
            else:
                raise ValueError("Plot type is not recognised.")


def peak_count(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug: bool = False,
    report: bool = False,
    prominence: int = 15,
    subpattern: str = "all",
    show_plots: bool = False,
    **kwargs,
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    if inputs:
        settings_class = inputs
    elif settings_class is None:
        settings_class = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_class = settings_class

    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_orders = len(settings_class.fit_orders)
    else:
        num_orders = len(subpattern)

    all_azis, all_data = read_saved_chunks(
        inputs=settings_class, debug=debug, report=report, subpattern=subpattern
    )

    # get the number of peaks
    # all_peaks, all_properties, count = peak_count(settings_class=settings_class, prominence=prominence)

    all_peaks: list[list] = []
    all_peakAzis: list = []
    all_properties: list[list] = []
    count = []
    peak_labels = []
    for j in range(num_orders):
        # loop over the number of sets of peaks fit for (i.e. len(settings_class.fit_orders))

        settings_class.set_subpattern(0, j)

        for k in range(len(settings_class.subfit_orders["peak"])):
            # loop over the number of peaks in each fit_orders

            all_peaks_tmp = []
            all_properties_tmp = []
            count_tmp = []
            if all_data[0][j]["h"][k]:
                # if there is some data in the array plot it.
                for i in range(settings_class.image_number):
                    peaks, properties = find_peaks(
                        all_data[i][j]["h"][k], prominence=prominence, width=0
                    )
                    properties["PeakAzis"] = np.array(all_azis[i][j])[peaks]
                    all_peaks_tmp.append(peaks)
                    # all_peakAzis.append(np.array(all_azis[i][j])[peaks])
                    all_properties_tmp.append(properties)
                    count_tmp.append(len(peaks))

                    if 1 and i == 0:
                        fig = plt.figure()
                        plt.plot(all_data[i][j]["h"][k])
                        plt.plot(peaks, np.array(all_data[i][j]["h"][k])[peaks], "x")

                        # plt.plot(np.zeros_like(x), "--", color="gray")

                        if show_plots is True:
                            plt.show()
                        else:
                            plt.close()

                all_peaks.append(all_peaks_tmp)
                all_properties.append(all_properties_tmp)
                count.append(count_tmp)
                # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                pk: int | Literal["all"] = k
                for l in range(len(settings_class.subfit_orders["peak"])):
                    if not all_data[0][j]["h"][l]:
                        pk = "all"
                # make the figure title
                ttlstr = IO.peak_string(settings_class.subfit_orders, peak=pk)
                peak_labels.append(ttlstr)

    return all_peaks, all_properties, count, peak_labels


def plot_peak_count(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug: bool = False,
    report: bool = False,
    prominence: int = 1,
    subpattern: str = "all",
    rotate: bool = False,
    show_plots: bool = False,
    **kwargs,
):
    """
    :param fit_parameters:
    :param fit_settings:
    :param settings_file:
    :param parallel:
    :param report:
    :param mode:
    :param save_all:
    :param inputs:
    :param debug:
    :param subpattern:
    :return:
    """

    if inputs:
        settings_class = inputs
    elif settings_class is None:
        settings_class = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_class = settings_class

    # get file times
    modified_time_s = np.array(
        [Path(file).stat().st_mtime for file in settings_class.image_list]
    )
    modified_time_s -= float(modified_time_s[0])

    y_label_str = r"Time (s)"
    # use file numbers if all times are the same
    if len(np.unique(modified_time_s)) == 1:
        modified_time_s = list(range(settings_class.image_number))
        y_label_str = r"Image in sequence"

    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)

    if subpattern == "all":
        num_orders = len(settings_class.fit_orders)
    else:
        num_orders = len(subpattern)

    # all_azis, all_data = read_saved_chunks(
    #    inputs=settings_class, debug=debug, report=report, subpattern=subpattern
    # )

    # get the number of peaks
    all_peaks, all_properties, count, titles = peak_count(
        settings_class=settings_class, prominence=prominence
    )

    fig, ax = plt.subplots()
    if not rotate:
        for j in range(len(all_peaks)):
            plt.plot(modified_time_s, count[j], ".-", label=titles[j])
        plt.xlabel(y_label_str)
        plt.ylabel(r"Number peaks")
    else:
        for j in range(len(all_peaks)):
            plt.plot(count[j], modified_time_s, ".-", label=titles[j])
        plt.xlabel(r"Number peaks")
        plt.ylabel(y_label_str)
    plt.legend()

    # Save the plot
    filename = IO.make_outfile_name(
        "PeakCountTime",
        directory=settings_class.output_directory,
        additional_text="prominence" + str(prominence),
        extension=".png",
        overwrite=True,
    )
    fig.savefig(filename, transparent=True)

    if show_plots is True:
        plt.show()
    else:
        plt.close()


"""
ths is probably a better way of finding peaks using ImageD111
Copied from ImageD11/test/peaksearchtiftest/scriptedpeaksearch.py
https://github.com/FABLE-3DXRD/ImageD11/blob/194d2fda453ee3e259e67651c3fcc1442dc3014b/test/peaksearchtiftest/scriptedpeaksearch.py
24th March 2023
"""
# import fabio, numpy as np
import glob

import fabio
import numpy.ma as ma
from ImageD11.blobcorrector import correctorclass, perfect
from ImageD11.columnfile import columnfile
from ImageD11.labelimage import labelimage
from ImageD11.peaksearcher import peaksearch


def execute2(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug=False,
    save_all=False,
    parallel=True,
    subpattern="all",
    mode="cascade",
    report=False,
    threshold=[1, 10, 100, 1000, 10000],
    **kwargs,
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
        settings_for_fit = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_for_fit = settings_class
    new_data = settings_for_fit.data_class

    if settings_for_fit.calibration_data:
        data_to_fill = settings_for_fit.calibration_data.resolve()
    else:
        data_to_fill = settings_for_fit.image_list[0].resolve()

    new_data.fill_data(
        data_to_fill,
        settings=settings_for_fit,
        debug=debug,
    )

    # plot calibration file
    if debug and settings_for_fit.calibration_data is not None:
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        new_data.plot_collected(fig_plot=fig, axis_plot=ax)
        plt.title("Calibration data")
        plt.show()
        plt.close()

    if not isinstance(threshold, list):
        threshold = [threshold]

    # if parallel processing start the pool
    # if parallel is True:
    #     p = mp.Pool(processes=mp.cpu_count())
    # p = mp.Pool()

    # restrict to sub-patterns listed
    settings_for_fit.set_subpatterns(subpatterns=subpattern)

    # Process the diffraction patterns #
    for j in range(settings_for_fit.datafile_number):
        print("Process ", settings_for_fit.datafile_list[j])
        # Get diffraction pattern to process.
        new_data.import_image(settings_for_fit.datafile_list[j], debug=debug)

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
            plt.title(settings_for_fit.datafile_list[j].stem)
            plt.show()
            plt.close()

        # give somes minimal options
        corrector = perfect()  # no spatial disortion
        dims = new_data.intensity.shape
        label_ims = {
            t: labelimage(
                shape=dims,
                fileout=make_outfile_name(
                    settings_for_fit.datafile_list[j],
                    directory=settings_for_fit.output_directory,
                )
                + "_t%d.flt" % (t),
                sptfile=make_outfile_name(
                    settings_for_fit.datafile_list[j],
                    directory=settings_for_fit.output_directory,
                )
                + "_t%d.spt" % (t),
                spatial=corrector,
            )
            for t in threshold
        }

        # for filename, omega in lines:

        # open image with Fabio, required by peaksearch.
        # can't do this as a data class function because it cant pickle a Fabio instance.
        frame = fabio.open(settings_for_fit.datafile_list[j])
        frame.data = new_data.intensity

        frame.header["Omega"] = 0
        frame.data = frame.data.astype(np.float32)
        # corrections like dark/flat/normalise would be added here
        peaksearch(
            settings_for_fit.datafile_list[j].stem,
            frame,
            corrector,
            threshold,
            label_ims,
        )
        for t in threshold:
            label_ims[t].finalise()

    # if debug:


def load_flts(settings_file=None, settings_class=None, inputs=None, **kwargs):
    """
    loads flt files for the selected subpattern.
    Makes an array of the peaks, where x,y are the centres and i is the size of the peak.
    Adds this array to the data class.

    Parameters
    ----------
    settings_file : TYPE, optional
        DESCRIPTION. The default is None.
    settings_class : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """

    if settings_class is None:
        settings_for_fit = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_for_fit = settings_class

    # if nothing has been set assume the first pattern
    if settings_for_fit.subfit_filename == None:
        settings_for_fit.set_subpattern(0, 0)

    # find all the *.flt files for this image
    fnam = make_outfile_name(
        settings_for_fit.subfit_filename, directory=settings_for_fit.output_directory
    )
    fnams = glob.glob(fnam + "*.flt")
    print(fnams)
    obj = []
    # obj = columnfile()
    for i in range(len(fnams)):
        try:
            obj.append(columnfile(fnams[i]))

            plt.scatter(-obj[i].dety, obj[i].detz, 1, c=(obj[i].sum_intensity))
            plt.colorbar()
            plt.title(fnams[i])
            plt.show()

            plt.scatter((obj[i].Number_of_pixels), (obj[i].sum_intensity))
            plt.title(fnams[i])
            plt.show()

        except:
            obj.append([])
    return obj


def make_im_from_flts(
    settings_file=None,
    settings_class=None,
    data_class=None,
    inputs=None,
    debug=False,
    **kwargs,
):
    if settings_class is None:
        settings_for_fit = initiate(settings_file, inputs=inputs, report=True)
    else:
        settings_for_fit = settings_class

    # if nothing has been set assume the first pattern
    if settings_for_fit.subfit_filename == None:
        settings_for_fit.set_subpattern(0, 0)

    if data_class == None:
        data_class = settings_for_fit.data_class
        data_class.fill_data(
            settings_for_fit.subfit_filename, settings=settings_for_fit
        )

    peaks_im = ma.zeros(data_class.intensity.shape)

    pks = load_flts(
        settings_file=settings_file,
        settings_class=settings_for_fit,
        data_class=data_class,
    )

    for i in range(len(pks)):
        try:
            x = -pks[i].dety
            y = pks[i].detz
            z = pks[i].sum_intensity

            for j in range(len(x)):
                peaks_im[int(y[j]), int(x[j])] = z[j]
        except:
            pass

    # peaks_im = peaks_im(mask=data_class.intensity.mask)
    peaks_im = ma.masked_where(ma.getmask(data_class.intensity), peaks_im)

    if debug:
        fig = plt.figure()
        plt.imshow(peaks_im, vmax=30)
        plt.colorbar()
        plt.show()

        fig = plt.figure()
        plt.imshow(np.log10(data_class.intensity))
        plt.colorbar()
        plt.show()

        fig = plt.figure()
        plt.scatter(data_class.tth, data_class.azm, s=0.1, c=peaks_im, vmax=30)
        plt.colorbar()
        plt.show()

    # print(pks[0].get_bigarray)
    # print(dir(pks[0]))

    data_class.peaks_image = peaks_im

    return data_class


"""
# ============================
# NOTES:
# These are a set of functions that attempted to find the spots in the raw data and count them.
# They were either very slow or didn't work.
# They are here to remind me to implement them at some time in the future because they are a 'better' way of making the measurements.
# ============================

def PeakCount(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):


    FitSettings, FitParameters, new_data = XRDFit.initiate(settings_file, inputs=inputs, report=True)

    # Get list of diffraction patterns #
    diff_files, n_diff_files = IO.FileList(FitParameters, FitSettings)

    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60

    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]

    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])

    FitSettings.fit_orders = orders_tmp


    FitSettings_files = copy.deepcopy(FitSettings.datafile_Files)
    num_peaks_all = []
    visible = []

    for h in range(n_diff_files):

        #logger.info(" ".join(map(str, [(FitParameters)])))
        #logger.info(" ".join(map(str, [(diff_files)])))

        #FitSettings_tmp = FitSettings
        FitSettings_tmp = FitSettings
        FitSettings_tmp.datafile_Files = [FitSettings_files[h]]
        #FitSettings_tmp.datafile_Files[0] = diff_files[h]

        #logger.info(" ".join(map(str, [(FitSettings_tmp.datafile_Files)])))
        #logger.info(" ".join(map(str, [(FitSettings_tmp)])))
        #logger.info(" ".join(map(str, [(FitParameters)])))

        outfname = IO.make_outfile_name(diff_files[h], directory=FitSettings_tmp.Output_directory, extension='csv', additional_text='peaks')
        logger.info(" ".join(map(str, [(outfname)])))

        if not os.path.exists(outfname):
            #This currently gets all the peaks in the image, not just those in the regios of interest
            intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data,
                    debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                    parallel=parallel,
                    mode='ReturnImage', report=True)

            # initialize with default parameters. The "denoise" parameter can be of use in your case
            # import 2D example dataset

            plt.imshow(twotheta)
            plt.colorbar()

            plt.imshow(np.array(twotheta))
            plt.colorbar()

            img = ma.array(intens)#, mask=False)
            plt.imshow(img)

            img.filled(0)
            img = ma.array(img, mask = None)
            plt.imshow(img)

            fp = findpeaks(limit=1, whitelist=['peak'], scale=False, togray=False, denoise=None)

            # make the fit
            fp.fit(img)
            # Make plot
            #fp.plot()
            #fp.plot_persistence()

            fp.results['persistence']
            logger.info(" ".join(map(str, [(fp.results['persistence'])])))
            logger.info(" ".join(map(str, [(type(fp.results['persistence']))])))

            peaks = fp.results
            #logger.info(" ".join(map(str, [(type(peaks))])))
            #logger.info(" ".join(map(str, [(peaks)])))
            #logger.info(" ".join(map(str, [(type(peaks))])))

            fp.results['persistence'].to_csv(outfname)
        else:

            peaks=pd.read_csv(outfname)
            #logger.info(" ".join(map(str, [(fp)])))
            #peaks = fp.results
            if h==0:
                #This currently gets all the peaks in the image, not just those in the regios of interest
                intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data,
                    debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                    parallel=parallel,
                    mode='ReturnImage', report=True)

        # get peaks in each diffraction peak window

        num_peaks = []
        for i in range(len(FitSettings.fit_orders)):

            tthRange = FitSettings.fit_orders[i]['range'][0]

            msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])


            # find and list all peaks in unmasked area.
            #FIXME: how does this interact with the region specific masks?
            pks=[]
            for j in range(len(peaks)):
                x = peaks['x'].loc[j]
                y = peaks['y'].loc[j]
                if twotheta.mask[y][x]==False:
                    if ((twotheta[y][x] > tthRange[0]) & (twotheta[y][x] < tthRange[1]) & (twotheta.mask[y][x] == False)):
                        pks.append(peaks.loc[j])
            num_peaks.append(len(pks))

            #approximate fraction of ring that is masked.
            #number of pixels in two theta range
            num_pix = ((np.array(twotheta) > tthRange[0]) & (np.array(twotheta) < tthRange[1])).sum()

            #number of unmasked pixles in two theta range
            num_vis_pix = ((tthRange[0] < twotheta) & (twotheta < tthRange[1])).sum()


            visible.append(num_vis_pix/num_pix)

            # # find and list all peaks in unmasked area.
            # #FIXME: how does this interact with the region specific masks?
            # peaks = []
            # for j in range(len(fp.results['persistence'])):
            #     x = fp.results['persistence']['x'].loc[j]
            #     y = fp.results['persistence']['y'].loc[j]
            #     if twotheta.mask[y][x]==False:
            #         if ((twotheta[y][x] > tthRange[0]) & (twotheta[y][x] < tthRange[1]) & (twotheta.mask[y][x] == False)):
            #             peaks.append(fp.results['persistence'].loc[j])
            # num_peaks.append(len(peaks))

            # #approximate fraction of ring that is masked.
            # #number of pixels in two theta range
            # num_pix = ((np.array(twotheta) > tthRange[0]) & (np.array(twotheta) < tthRange[1])).sum()

            # #number of unmasked pixles in two theta range
            # num_vis_pix = ((tthRange[0] < twotheta) & (twotheta < tthRange[1])).sum()


            # visible.append(num_vis_pix/num_pix)
        num_peaks_all.append(num_peaks)
        #logger.info(" ".join(map(str, [(num_peaks)])))
        #logger.info(" ".join(map(str, [(visible)])))

    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60
    num_peaks_all=np.array(num_peaks_all)
    logger.info(" ".join(map(str, [(num_peaks_all)])))
    logger.info(" ".join(map(str, [(modified_time)])))
    if 1:
        fig,ax = plt.subplots()
        for i in range(len(FitSettings.fit_orders)):
            logger.info(" ".join(map(str, [(num_peaks_all[:,i])])))
            plt.plot(modified_time, num_peaks_all[:,i], '.-', label=IO.peak_string(FitSettings.fit_orders[i]))

        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Intensity Maxima Count')
        plt.legend()
        plt.show()

        fig.savefig('IntensityMaximaTime2.png')
        fig.savefig('IntensityMaximaTime2.pdf')
        fig.savefig('IntensityMaximaTime2.eps')

    #plt.scatter([*range(len(num_peaks_all))], num_peaks_all)
    #plt.scatter

        # #fraction that is unmasked (as approximation for how much of ring is visible.)

        # tthRange = FitSettings.fit_orders[i]['range'][0]

        # msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])

        # img = ma.array(intens, mask=msk.mask)

        # plt.imshow(img)

        # img.filled(0)
        # img = ma.array(img, mask = None)
        # plt.imshow(img)

        # fp = findpeaks(limit=2, whitelist=['peak'])

        # # make the fit
        # fp.fit(img)
        # # Make plot

        # fp.plot()

        # fp.plot_persistence()

        # fp.results['persistence']

        # logger.info(" ".join(map(str, [(fp.results)])))










    for i in range(len(FitSettings.fit_orders)):


        tthRange = FitSettings.fit_orders[i]['range'][0]

        msk = ma.masked_outside(twotheta, tthRange[0], tthRange[1])

        img = ma.array(intens, mask=msk.mask)

        plt.imshow(img)

        img.filled(0)
        img = ma.array(img, mask = None)
        plt.imshow(img)

        fp = findpeaks(limit=2, whitelist=['peak'])

        # make the fit
        fp.fit(img)
        # Make plot

        fp.plot()

        fp.plot_persistence()

        fp.results['persistence']

        logger.info(" ".join(map(str, [(fp.results)])))







    if 1:
        fig,ax = plt.subplots()

        for i in range(len(FitSettings.fit_orders)):
            plt.plot(modified_time, peak_count[:][i], '.-', label=IO.peak_string(FitSettings.fit_orders[i]))

        plt.xlabel(r'Time (min)')
        plt.ylabel(r'Intensity Maxima Count')
        plt.legend()
        plt.show()

        fig.savefig('IntensityMaximaTime.png')
        fig.savefig('IntensityMaximaTime.pdf')
        fig.savefig('IntensityMaximaTime.eps')



def WatershedCount(settings_file=None, inputs=None, debug=False, refine=True, save_all=False, propagate=True, iterations=1,
            track=False, parallel=True, subpattern='all', **kwargs):


    FitSettings, FitParameters, new_data = XRDFit.initiate(settings_file, inputs=inputs, report=True)

    # Get list of diffraction patterns #
    diff_files, n_diff_files = IO.FileList(FitParameters, FitSettings)

    #get file times
    modified_time = []
    for i in range(n_diff_files):
        modified_time.append(os.path.getmtime(diff_files[i]))
    modified_time = np.array(modified_time)
    modified_time = (modified_time - modified_time[0])/60

    # restrict to subpatterns listed
    if subpattern=='all':
        subpats = list(range(0, len(FitSettings.fit_orders)))
    elif isinstance(subpattern,list):
        subpats = subpattern
    else:
        subpats = [int(x) for x in str(subpattern)]

    # make new order search list
    orders_tmp = []
    for i in range(len(subpats)):
        j = subpats[i]
        orders_tmp.append(FitSettings.fit_orders[j])

    FitSettings.fit_orders = orders_tmp


    FitSettings_files = copy.deepcopy(FitSettings.datafile_Files)
    num_peaks_all = []
    visible = []

    for h in range(n_diff_files):

        #This currently gets all the peaks in the image, not just those in the regios of interest
        intens, twotheta, azimu = XRDFit.execute(FitSettings=FitSettings_tmp, FitParameters=FitParameters, inputs=new_data,
                debug=debug, refine=refine, save_all=save_all, propagate=propagate, iterations=iterations,
                parallel=parallel,
                mode='ReturnImage', report=True)

        im = img_as_float(intens)

        # image_max is the dilation of im with a 20*20 structuring element
        # It is used within peak_local_max function
        image_max = ndi.maximum_filter(im, size=20, mode='constant')

        # Comparison between image_max and im to find the coordinates of local maxima
        coordinates = peak_local_max(im, min_distance=20)

        # display results
        fig, axes = plt.subplots(1, 3, figsize=(8, 3), sharex=True, sharey=True)
        ax = axes.ravel()
        ax[0].imshow(im, cmap=plt.cm.gray)
        ax[0].axis('off')
        ax[0].set_title('Original')

        ax[1].imshow(image_max, cmap=plt.cm.gray)
        ax[1].axis('off')
        ax[1].set_title('Maximum filter')

        ax[2].imshow(im, cmap=plt.cm.gray)
        ax[2].autoscale(False)
        ax[2].plot(coordinates[:, 1], coordinates[:, 0], 'r.')
        ax[2].axis('off')
        ax[2].set_title('Peak local max')

        fig.tight_layout()

        plt.show()





"""
