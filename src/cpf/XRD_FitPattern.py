#!/usr/bin/env python

from __future__ import annotations  # Enables additional Python features

__all__ = ["execute", "write_output"]

import json
import logging
import sys
import inspect
from importlib import import_module
from os import cpu_count
from pathlib import Path
from types import ModuleType
from typing import Literal, Optional

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
from cpf.output_formatters.fits_io import WriteFits, ReadFits_to_list
from cpf.series_functions import get_series_mean
from cpf.settings import Settings, is_settings, get_settings
from cpf.util.logging import get_logger, set_global_log_level
from cpf.XRD_FitSubpattern import fit_sub_pattern

np.set_printoptions(threshold=sys.maxsize)
# FIX ME: Need to add complexity here.
# Need option to check required functions are present if adding new output format.
# Also need to do something similar with data class


logger = get_logger("cpf.XRD_FitPattern")



__doc__ = "This is the main function for fitting the diffraction (or dispersed) peak data. "



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


def initiate(
    settings: Optional[str | Path | dict | Settings()] = None,
    inputs=None,
    out_type=None,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs,
):
    """
    Run checks on input files, initiate data class and check output options

    :param settings:
    :param report:
    :param out_type:
    :param initiate_data:
    :param inputs:
    :return fit_parameters:
    :return fit_settings:
    """

    # Set the logger level for this run
    set_global_log_level(report)

    # Add a file handler to this logger
    if isinstance(settings, dict):
        running_name = settings.get("run_name", "cpf_logging_file")
        setting_type = "dictionary"
    elif isinstance(settings, str) or isinstance(settings, Path):
        running_name = settings
        setting_type = "file" #assume string deontes file name
    elif is_settings(settings): #isinstance(settings, type(Settings()))
        setting_type = "class"
        if settings.settings_file is not None:
            #there is a file name or run label in the settings class
            running_name = settings.settings_file
        else:
            running_name = "cpf settings class"     
    else: #unknown
        setting_type = "Unknown"
        running_name = "cpf_run"     
    log_file = log_file = make_outfile_name(
        base_filename=running_name, extension=".log", overwrite=True
    )
    logger.add_file_handler(log_file)

    # make a header in the log file so that we know where the processing starts
    # It is the first pass through the method if: 
    # 1. if settings is not a settings class then it has to be new.
    # 2. if it is a settings class and is not validated then is likely to be new
    # 3. but it could be a validated settings class and also new executeion. 
    #       but in normal use this will not be the case. Therefore can ignore for now.
    #       and if the settings class has been validated then it will appear in the log 
    #       file as an initiation or validation.
    # all calls to the methods in here pass around a settings class once it is made and 
    # so is_settings == True.
    if (not is_settings(settings) or 
        (is_settings(settings) and settings.is_valid()==False)
        ): 
        logger.info("")
        logger.info("=================================================================")
        logger.info("")
        logger.info(f"Starting data proceesing using settings {setting_type}{' from' if setting_type != 'file' else ''}: {running_name}")    
        logger.info("")
        logger.info("=================================================================")
        logger.info("")
        
    settings_class = get_settings(settings, **kwargs)
    
    return settings_class



def view(
    settings: [str | Path | dict | Settings()],
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
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs
):
    """
    :param settings:
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
    
    settings_class = initiate(settings, report=report, **kwargs)
    # make a note in the logger.
    # suppress output if called by another module.
    for i in range(len(inspect.stack())-1,-1,-1):
        if inspect.stack()[i].function == "<module>":          
            base_call = inspect.stack()[i-1].function
    if base_call == "view":
        logger.info("")
        logger.info(f"Running: XRD_FitPattern.view with settings: {settings_class.settings_file}")
        logger.info("")

    # view the listed file only
    if pattern != "all":
        # restrict file list to first file
        settings_class.set_data_files(keep=pattern)

    write_output(settings_class, out_type="CollectionMovie", **kwargs)
    
    # write_output(settings_file=settings_file, out_type="RangesMovie")

    # execute(
    #     settings_class=settings_class,
    #     debug=debug,
    #     refine=refine,
    #     save_all=save_all,
    #     iterations=iterations,
    #     parallel=parallel,
    #     mode="view",
    #     report=True,
    # )


def set_range(
    settings: [str | Path | dict | Settings()],
    inputs=None,
    debug: bool = False,
    refine: bool = True,
    save_all: bool = False,
    # propagate: bool = True,
    iterations: int = 1,
    # track: bool = False,
    parallel: bool = True,
    subpattern: str = "all",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs
):
    """
    :param settings:
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

    settings_class = initiate(settings, report=report, **kwargs)
    # make a note in the logger.
    # suppress output if called by another module.
    for i in range(len(inspect.stack())-1,-1,-1):
        if inspect.stack()[i].function == "<module>":          
            base_call = inspect.stack()[i-1].function
    if base_call == "set_range":
        logger.info("")
        logger.info(f"Running: XRD_FitPattern.set_range with settings: {settings_class.settings_file}")
        logger.info("")
    
    # search over the first file only
    # restrict file list to first file
    settings_class.set_data_files(keep=0)
    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)
    # circulment revalidating the settings class
    settings_class._unmodified_self = settings_class._validation_copy()
    
    execute(
        settings_class,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        mode="set-range",
        report=report,
    )


def initial_peak_position(
    settings: [str | Path | dict | Settings()],
    inputs=None,
    debug: bool = False,
    refine: bool = True,
    save_all: bool = False,
    # propagate: bool = True,
    iterations: int = 1,
    # track: bool = False,
    parallel: bool = True,
    subpattern: str = "all",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs
):
    """
    Calls interactive graph to set the inital peak postion guesses.

    The event handler code is copied from:
    https://matplotlib.org/stable/users/event_handling.html for how to make work

    :param settings:
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
    
    settings_class = initiate(settings, report=report, **kwargs)
    # make a note in the logger.
    # suppress output if called by another module.
    for i in range(len(inspect.stack())-1,-1,-1):
        if inspect.stack()[i].function == "<module>":          
            base_call = inspect.stack()[i-1].function
    if base_call == "initial_peak_position":
        logger.info("")
        logger.info(f"Running: XRD_FitPattern.initial_peak_position with settings: {settings_class.settings_file}")
        logger.info("")
    
    # search over the first file only
    settings_class.set_data_files(keep=0)
    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)
    # circulment revalidating the settings class
    settings_class._unmodified_self = settings_class._validation_copy()

    logger.info("\n'initial_peak_position' needs an interactive matplotlib figure.")
    logger.info("If you are using sypder with inline figures, call '%matplotlib qt', then rerun the script")
    logger.info("To restore the inline plotting afterwards call '%matplotlib inline'")
    logger.info("To move to the next peak selection close the window.\n")

    execute(
        settings_class,
        debug=debug,
        refine=refine,
        save_all=save_all,
        iterations=iterations,
        parallel=parallel,
        mode="set-guess",
        report=report,
        **kwargs
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
    settings: [str | Path | dict | Settings()],
    inputs=None,
    refine: bool = True,
    save_all: bool = False,
    parallel: bool = False,
    search_image: [int | list | str] = 0,
    search_parameter: str = "height",
    search_over: list[int] = [0, 20],
    subpattern: str = "all",
    search_peak: int = "all",
    search_series: list[str] = ["fourier", "spline"],
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs
):
    """
    Searches for the best order to use for 'search_parameter', where 'search_over'
    is one of the model parameters (e.g. 'background', 'width', etc). 
    The data is fit with orders in the range defined by 'search_over'.
    The resultant fits are plotted by cpf.output_formatters.WriteOrderSearchFigures
    
    Makes a json file with all the fits that is named:
        *settings_file*__search=*search_parameter*_subpattern=*subpattern*_peak=*search_peak*

    The function will technically work in parallel but due to memory limitations 
    and the way the code is structured, it is set by default to run this function in 
    series. 
    For the same reason, althogh the code can run 'all' the peaks at once, it is 
    split up and loops over each peak separately.


    Parameters
    ----------
    settings : *.py file, string, Path or cpf Settings
        text file containing all the fitting parameters. The default is None.
    inputs : TYPE, optional
        DESCRIPTION. The default is None.
    refine : bool, optional
        DESCRIPTION. The default is True.
    save_all : bool, optional
        DESCRIPTION. The default is False.
    parallel : bool, optional
        Process the data in parallel? Set to false because parallel fills the memory with data.
        USE WITH CAUTION. The default is False.
    search_image : str, int, optional
        Which image to use in the order search. Can be integer image in list, "mid" or -1 for last image.
        The default is 0.
    search_parameter : str, optional
        DESCRIPTION. The default is "height".
    search_over : list[int], optional
        DESCRIPTION. The default is [0, 20].
    subpattern : str, optional
        DESCRIPTION. The default is "all".
    search_peak : int, optional
        DESCRIPTION. The default is "all".
    search_series : list[str], optional
        DESCRIPTION. The default is ["fourier", "spline"].
    report : Literal[        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"    ], optional
        DESCRIPTION. The default is "INFO".

    Returns
    -------
    None.

    """
    
    settings_class = initiate(settings, report=report, **kwargs)
    # make a note in the logger.
    # suppress output if called by another module.
    for i in range(len(inspect.stack())-1,-1,-1):
        if inspect.stack()[i].function == "<module>":          
            base_call = inspect.stack()[i-1].function
    if base_call == "order_search":
        logger.info("")
        logger.info(f"Running: XRD_FitPattern.order_search with settings: {settings_class.settings_file}")
        logger.info("")

    # search over the first file only
    settings_class.set_data_files(keep=search_image)
    settings_class.fit_propagate = False
    
    # loop over the peaks in turn unless forced
    if subpattern =="force all":
        subpattern = ["all"]
    elif subpattern =="all":
        subpattern = list(range(len(settings_class.fit_orders)))
    elif not isinstance(subpattern, list):
        subpattern = [subpattern]

    for i in range(len(subpattern)):

        logger.info(f"Performing order_search for peak {i}")
        
        # set search orders and execute
        settings_class.set_order_search(
            search_parameter=search_parameter,
            search_over=search_over,
            subpatterns=subpattern[i],
            search_peak=search_peak,
            search_series=search_series,
        )
        settings_class.fit_propagate = False
        settings_class.file_label = (
            "search="
            + search_parameter
            + "_subpattern="
            + str(subpattern[i])
            + "_peak="
            + str(search_peak)
        )
        # circulment revalidating the settings class
        settings_class._unmodified_self = settings_class._validation_copy()
        
        execute(
            settings_class,
            refine=refine,
            save_all=save_all,
            mode="search",
            parallel=parallel,
            report=report,
        )
    
        # call WriteOrderSearchFigures to make the figures.
        write_output(
            settings_class,
            debug=True,
            out_type="OrderSearchFigures",
        )
        
        write_output(
            settings_class,
            debug=True,
            out_type="OrderSearchMovie",
        )

        settings_class.unset_order_search()
    logger.info("Order searches are completed.")

def write_output(
    settings,
    parms_dict=None,
    out_type=None,
    det=None,
    use_bounds: bool = False,
    differential_only: bool = False,
    debug: bool = False,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs,
):
    """
    
    :param settings : *.py file, string, Path or cpf Settings
    :param debug:
    :param fit_parameters:
    :param use_bounds:
    :param differential_only:
    :param fit_settings:
    :param parms_dict:
    :param out_type:
    :param det:
    :return:
    """
    
    settings_class = initiate(settings, report=report, **kwargs)
    # make a note in the logger.
    # suppress output if called by another module.
    for i in range(len(inspect.stack())-1,-1,-1):
        if inspect.stack()[i].function == "<module>":          
            base_call = inspect.stack()[i-1].function
    if base_call == "write_output":
        logger.info("")
        logger.info(f"Running: XRD_FitPattern.write_output with settings: {settings_class.settings_file}")
        logger.info("")

    if out_type is not None:
        logger.moreinfo(f"Output_type was provided as an option; will use {out_type}")
        settings_class.set_output_types(out_type_list=out_type)

    if settings_class.output_types is None:
        logger.warning("There are no output types. Add 'Output_type' to input file or specify 'out_type' in command.")
    else:
        for mod in settings_class.output_types:
            logger.info(" ".join(map(str, [("Writing output file(s) using %s" % mod)])))
            wr = output_methods_modules[mod]
            wr.WriteOutput(
                settings_class,
                differential_only=differential_only,
                debug=debug,
                **kwargs,
            )


def execute(
    settings: [str | Path | dict | Settings()],
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
    resume: bool = False,
    mode: str = "fit",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    fit_method: str = "leastsq",
    **kwargs,
):
    """
    :param settings : *.py file, string, Path or cpf Settings
    :param fit_parameters:
    :param fit_settings:
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
    
    # if not is_settings(settings):
    settings_class = initiate(settings, report=report, **kwargs)
    # make note in logger
    # suppress output if called by another module.
    for i in range(len(inspect.stack())-1,-1,-1):
        if inspect.stack()[i].function == "<module>":          
            base_call = inspect.stack()[i-1].function
    if base_call == "execute":
        logger.info("")
        logger.info(f"Running: XRD_FitPattern.execute with settings: {settings_class.settings_file}")
        logger.info("")
    
    #get data from settings class
    new_data = settings_class.data_class

    # Define locally required names
    temporary_data_file = make_outfile_name(
        "PreviousFit_JSON",
        directory=settings_class.output_directory,
        extension=".dat",
        overwrite=True,
    )

    as_masked = kwargs.pop('as_masked', False)
    if (mode == "set-range" or mode == "view"):
        as_masked = True
    else:
        as_masked = as_masked
        if as_masked == True:
            logger.warning("'as_masked'==True changes the fit for some masked datasts. I dont know why. Check fits with and without this setting")
        
    if settings_class.calibration_data:
        data_to_fill = Path(settings_class.calibration_data).resolve()
    else:
        data_to_fill = settings_class.image_list[0]

    new_data.fill_data(
        data_to_fill,
        settings=settings_class,
        debug=debug,
    )

    # Get calibration parameter file
    parms_dict = new_data.calibration
    # FIXME this should be removable.

    # plot calibration file
    if (
        logger.is_below_level(level="DEBUG")
        and settings_class.calibration_data is not None
    ):
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        new_data.plot_collected(fig_plot=fig, axis_plot=ax)
        plt.title("Calibration data")
        plt.show()
        plt.close()

    # if parallel processing start the pool
    if parallel is True:
        nodes = int(np.min([cpu_count(), len(settings_class.fit_orders)]))
        pool = ParallelPool(nodes=nodes)
        # Since we may have already closed the pool, try to restart it
        try:
            pool.restart()
        except AssertionError:
            pass

    # Process the diffraction patterns
    # for j in range(settings_class.image_number):
    progress = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger
    for j in progress.iter_bar(image=range(settings_class.image_number)):
        logger.info(f"Processing {title_file_names(image_name=settings_class.image_list[j])}")

        settings_class.set_subpattern(j, 0)

        # Get diffraction pattern to process.
        new_data.import_image(settings=settings_class, debug=debug)

        # get json file name for outputs.
        if mode == "search":
            additional_text = settings_class.file_label
        else:
            additional_text = None
        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            additional_text=additional_text,
            extension=".json",
            overwrite=True,
        )
        
        # if the output file already exists and resume is true then skip
        # this iteration        
        if resume == True and Path(filename).is_file():
            logger.info(f"  {title_file_names(image_name=settings_class.image_list[j])} has already been processed -- skipping")
            continue
        # else do the process.

        if ((isinstance(settings_class.datafile_preprocess, dict) ) or #settings_class.datafile_preprocess is not None or 
            (isinstance(settings_class.calibration_mask, dict) and "threshold" in settings_class.calibration_mask)
            ):
            # needed because image preprocessing adds to the mask and is different for each image.
            # set intenstiy threshold and/or cosmics for each frame. 
            # only applies these because all other mask functions are static.
            new_data.mask_restore()
            if (isinstance(settings_class.datafile_preprocess, dict) and "cosmics" in settings_class.datafile_preprocess):
                new_data = cosmicsimage_preprocess(new_data, settings_class)
            if (isinstance(settings_class.calibration_mask, dict) and "threshold" in settings_class.calibration_mask):
                # set intenstiy threshold for each frame. 
                # only applies to threshold because all other mask functions are static
                # (cannot change between frames)
                new_data.set_mask(intensity_bounds = settings_class.calibration_mask["threshold"])
        else:
            # nothing is done here.
            pass

        # plot input file
        if logger.is_below_level(level="DEBUG") or mode == "view":
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity")
            plt.title(title_file_names(settings_for_fit=settings_class, num=j))
            plt.show()
            if mode == "view":
                filename = make_outfile_name(
                    settings_class.image_list[j],
                    directory=settings_class.output_directory,
                    extension=".png",
                    overwrite=True,
                )

                fig.savefig(filename)

            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax_o1 = plt.subplot(111)
            new_data.plot_integrated(fig_plot=fig, axis_plot=ax, show="intensity")
            # plt.title(os.path.basename(settings_class.datafile_list[j]))
            plt.title(title_file_names(settings_for_fit=settings_class, num=j))
            plt.show()
            if mode == "view":
                filename = make_outfile_name(
                    settings_class.image_list[j],
                    directory=settings_class.output_directory,
                    additional_text = "integrated",
                    extension=".png",
                    overwrite=True,
                )

                fig.savefig(filename)

                print("Plotted data.")
            else:
                plt.close()

        # Get previous fit (if it exists and is required)
        if (
            Path(temporary_data_file).is_file()
            and settings_class.fit_propagate is True
            and mode == "fit"
            and j != 0  # not the first data in series.
        ):
            # Read JSON data from file
            logger.moreinfo(f"Loading previous fit results from {temporary_data_file}.")
            previous_fit, _ = ReadFits_to_list(temporary_data_file, replace=False)
            # if the previous_fit is not the same size as fit_orders the inout file must have been changed.
            # so discard the previous fit and start again.
            if len(previous_fit) != len(settings_class.fit_orders):
                del previous_fit

        # Switch to save the first fit in each sequence.
        save_figs = False#True if (j == 0 or save_all is True) else False

        # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        fitted_param = []
        lmfit_models = []
        parallel_pile = []

        for i in range(len(settings_class.fit_orders)):
            
            if parallel == False:
                serial_string = f"Fitting range {i+1}/{len(settings_class.fit_orders)}"
                logger.info(serial_string)

            # get settings for current subpattern
            settings_class.set_subpattern(j, i)

            if "previous_fit" in locals() and mode == "fit":
                params = previous_fit[i]
                params.pop("correlation_coeffs", None)
            else:
                params = None

            # Track the position of the peak centroid
            # FIXME: This is crude - the range doesn't change width. so can't account for massive change in stress.
            # But does it need to?
            tth_range = np.array(settings_class.subfit_orders["range"])
            if settings_class.fit_track is True and "previous_fit" in locals():
                null_terms = any_terms_null(params, val_to_find=None)
                if null_terms == True:
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
                    # if tacking things start from positon of previous fit positions
                    tth_range = previous_fit[i]["range"][0]
                    mid = []
                    for k in range(len(params["peak"])):
                        mid.append(get_series_mean(params['peak'][k], "d-space"))

                    if mid == None or mid == 0:
                        # the previous fits failed in some way.
                        move_by = 0
                    else:
                        cent = new_data.conversion(np.mean(mid), reverse=True)
                        move_by = cent - np.mean(tth_range)

                    # update tth_range and settings
                    tth_range = tth_range + move_by
                    settings_class.fit_orders[i]["range"] = (
                        settings_class.fit_orders[i]["range"] + move_by
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
                    if "PeakPositionSelection" in settings_class.fit_orders[i]:
                        for k in range(
                            len(settings_class.fit_orders[i]["PeakPositionSelection"])
                        ):
                            settings_class.fit_orders[i]["PeakPositionSelection"][k][
                                2
                            ] = (
                                settings_class.fit_orders[i]["PeakPositionSelection"][
                                    k
                                ][2]
                                + move_by
                            )

                    # re-get settings for current subpattern
                    settings_class.set_subpattern(j, i)

            sub_data = new_data.duplicate_without_detector(range_bounds=tth_range, as_masked=as_masked)

            # Mask the subpattern by intensity if called for
            if (
                "imax" in settings_class.subfit_orders
                or "imin" in settings_class.subfit_orders
            ):
                sub_data = SpotProcess(sub_data, settings_class, as_masked=as_masked)

            if mode == "set-range":
                fig_1 = plt.figure()
                sub_data.plot_masked(fig_plot=fig_1, **kwargs)
                plt.suptitle(peak_string(settings_class.subfit_orders) + "; masking")

                filename = make_outfile_name(
                    settings_class.image_list[j],
                    directory=settings_class.output_directory,
                    additional_text="mask",
                    orders=settings_class.subfit_orders,
                    extension=".png",
                    overwrite=True,
                )

                fig_1.savefig(filename)

                # if debug:
                plt.show()
                plt.close()

            elif mode == "view":
                fig = plt.figure()
                ax = fig.add_subplot(1, 1, 1)
                ax_o1 = plt.subplot(111)
                sub_data.plot_calibrated(
                    fig_plot=fig, axis_plot=ax, show="intensity", **kwargs #rastered="scatter"
                )
                plt.suptitle(
                    peak_string(settings_class.subfit_orders) + "; calibrated"
                )

                filename = make_outfile_name(
                    settings_class.image_list[j],
                    directory=settings_class.output_directory,
                    additional_text="range",
                    orders=settings_class.subfit_orders,
                    extension=".png",
                    overwrite=True,
                )

                fig.savefig(filename)

                plt.show()
                if mode == "view":
                    print("Plotted data.")
                else:
                    plt.close()

            elif mode == "set-guess":
                fig_1 = plt.figure()
                ax = fig_1.add_subplot(1, 1, 1)
                ax_o1 = plt.subplot(111)
                sub_data.plot_calibrated(
                    fig_plot=fig_1, axis_plot=ax, y_axis="azimuth", limits=[0, 100],
                    **kwargs
                )
                plt.title(peak_string(settings_class.subfit_orders))

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
                                    % peak_string(settings_class.subfit_orders)
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
                        "min_data_intensity": settings_class.fit_min_data_intensity,
                        "min_peak_intensity": settings_class.fit_min_peak_intensity,
                    }
                    arg = (sub_data, settings_class.duplicate_without_dataclass())
                    parallel_pile.append((arg, kwargs))

                else:  # non-parallel version
                    tmp = fit_sub_pattern(
                        sub_data,
                        settings_class.duplicate_without_dataclass(),  # added
                        params,
                        save_fit=save_figs,
                        debug=debug,
                        refine=refine,
                        iterations=iterations,
                        min_data_intensity=settings_class.fit_min_data_intensity,
                        min_peak_intensity=settings_class.fit_min_peak_intensity,
                        fit_method=fit_method,
                        **kwargs
                    )
                    fitted_param.append(tmp)

        # write output files
        if mode == "fit" or mode == "search":
            if parallel is True:
                tmp = pool.map(parallel_processing, parallel_pile)
                for i in range(len(settings_class.fit_orders)):
                    fitted_param.append(tmp[i])
            
            # store the fit parameters' information as a JSON file.
            WriteFits(settings_class, fitted_param, data_class=new_data, mode=mode)

            # if propagating the fits write them to a temporary file
            if settings_class.fit_propagate:
                WriteFits(settings_class, fitted_param, filename_to_write=temporary_data_file)

    if mode == "fit":
        # Write the output files.
        write_output(
            settings_class, debug=debug
        )

    if parallel is True:
        pool.clear()

def parallel_processing(p):
    a, kw = p
    return fit_sub_pattern(*a, **kw)

if __name__ == "__main__":
    # Load settings fit settings file.
    sys.path.append(str(Path().cwd()))
    settings = Path(sys.argv[1])
    logger.info(" ".join(map(str, [(settings)])))
    # Safely exit the program
    for handler in logging.getLogger().handlers:
        handler.flush()
    # sys.exit()
    execute(
        settings,
        inputs=None,
        debug=False,
        refine=True,
        save_all=False,
        iterations=1,
        parallel=False,
    )
