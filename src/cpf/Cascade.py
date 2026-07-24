




"""

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

"""
SAH notes:
    
cascade calls XRD_FitPAttern to get bits to run. 

Execute runs the code and returns a list of 'spots' (however defined) and associated positions (e.g. location and intensity).     

write a panda data frame with columns:
    filename,
    (metadata)
    possible peak, (is this peak in the two theta range of a defined peak?)
    tth,
    tth_err,
    tth_width,
    tth_width_err,
    azm,
    azm_err,
    axm_width,
    azm_width_err,
    intensity, (height of maximum abve background as determined by algorithm)
    intensity_err,
    total counts, 
    total_counts_err,
    + extra columns depending on method used.



"""


# from __future__ import annotations

__all__ = ["initiate", "set_range", "execute"]

import json
import logging
import sys
import inspect
from importlib import import_module
from os import cpu_count
from pathlib import Path
from typing import Literal, Optional
# import pandas as pd

# import matplotlib.colors as colours
import matplotlib.pyplot as plt
import numpy as np
import proglog
from pathos.pools import ParallelPool
# from scipy.signal import find_peaks

# from cpf.XRD_FitPattern import initiate, view, set_range#, write_output
from cpf.XRD_FitPattern import initiate as xrd_fp_initiate
from cpf.XRD_FitPattern import view as xrd_fp_view
from cpf.XRD_FitPattern import set_range as xrd_fp_set_range


# from cpf.BrightSpots import SpotProcess
from cpf.data_preprocess import remove_cosmics as cosmicsimage_preprocess
from cpf.settings import Settings
from cpf.util.io import (
    # has_value,
    make_outfile_name,
    # numpy_to_json,
    # peak_string,
    title_file_names,
)
from cpf.output_formatters.fits_io import WriteFits
from cpf.util.logging import get_logger
# from cpf.XRD_FitSubpattern import fit_sub_pattern
from types import ModuleType
from cpf import spot_methods
from cpf import spot_outputs

logger = get_logger("cpf.Cascade")



def register_default_spot_finders() -> dict[str, ModuleType]:
    """
    Load all available spot finding modules
    These are files in the output_formatters folder that have a name that is 
    of the form `Write*.pn', where * is the name used to call the output formatter

    Returns
    -------
    dict[str, ModuleType]
        doctionary of possuble spot finding modules
    """
    # FIX ME: We could add extra checks here to make sure the required functions exist in each case.
    output_list = spot_methods.module_list
    new_module = {}
    for output_module in output_list:
        module: ModuleType = import_module(f"cpf.spot_methods.{output_module}")
        new_module[output_module.replace("spots", "").strip("_")] = module
    return new_module

# Load potential output formats
spot_finding_methods_modules = register_default_spot_finders()


def register_default_spot_output_formats() -> dict[str, ModuleType]:
    """
    Load all available output modules. 
    These are files in the pattern_processing_methods folder that have a name that is 
    of the form `Write*.pn', where * is the name used to call the output formatter

    Returns
    -------
    dict[str, ModuleType]
        doctionary of possuble output modules.
    """
    # FIX ME: We could add extra checks here to make sure the required functions exist in each case.
    output_list = spot_outputs.module_list
    new_module = {}
    for output_module in output_list:
        module: ModuleType = import_module(f"cpf.spot_outputs.{output_module}")
        new_module[output_module.replace("Write", "")] = module
    return new_module

# # Load potential output formats
spot_output_methods_modules = register_default_spot_output_formats()



def initiate(*args, **kwargs
) -> Settings():
    """
    Takes input and creates a Settings class object, which is used 
    to run the fitting processes.
        
    This method is a wrapper for cpf.XRD_FitPattern.initiate()
    
    Parameters
    ----------
    settings : Optional[str | Path | dict | Settings()]
        Pointer to information needed for settings class. Can be of the form:        
        string -- filename of python formatted file 
        Path -- path for python formatted file  
        dict -- dictionary of settings
        Settings() -- cpf Settings class 
    report : Literal[ "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"    ], optional
        Logger level for how much information to write to the log files.
        The default is "INFO".
    **kwargs : key, value pairs
        key, value arguments arguments. Passed though the method to cpf.settings.Settings() but not used here. 
        See documention of these for kwarg use.

    Returns
    -------
    settings_class : Settings()
        A class containing the processing parameters for continuous peak fit.
    """
    
    # pass everything through to XRD_FitPattern.initiate
    settings_class = xrd_fp_initiate(*args, **kwargs)
    return settings_class


def view(*args, **kwargs):
    """
    Plot the data using the calubration.
    
    Makes a movie file of the data using output "CollectionMovie".  
    
    This method is a wrapper for cpf.XRD_FitPattern.view()
    
    Parameters
    ----------
    settings : Optional[str | Path | dict | Settings()]
        Pointer to information needed for settings class. Can be of the form:        
        string -- filename of python formatted file 
        Path -- path for python formatted file  
        dict -- dictionary of settings
        Settings() -- cpf Settings class 
    pattern : int, str, optional
        Which patterns in the series to plot. 
        Either "all", "mid", or integer list of patterns to keep.
        The default is "all".
    report : Literal[ "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR" ], optional
        Logger level for how much information to write to the log files.
        The default is "INFO".
    **kwargs : key, value pairs
        key, value arguments arguments. Passed though the method to to cpf.XRD_Fitpattern.initiate() and 
        cpf.XRD_Fitpattern.execute(). See documention of these for kwarg use.

    Returns
    -------
    None.

    """
    
    # pass everything through to XRD_FitPattern.view
    xrd_fp_view(*args, **kwargs)
    
    
    
def set_range(*args, **kwargs):
    """
    Plot data in the elected ranges using the calubration and how the masks are applied to the 
    data. Subgigues produced:
        - unmasked data 
            - azimuth vs two theta, colours scale = intensity
            - intensity vs two theta, colours scale = azimuth
        - mask
        - masked data
            - azimuth vs two theta, colours scale = intensity
            - intensity vs two theta, colours scale = azimuth
        - Cumulative distribution function of unmasked and masked intensities
    
    Makes a separate image file for each of the elected ranges.  

    This method is a wrapper for cpf.XRD_FitPattern.set_range()
    
    Parameters
    ----------
    settings : Optional[str | Path | dict | Settings()]
        Pointer to information needed for settings class. Can be of the form:        
        string -- filename of python formatted file 
        Path -- path for python formatted file  
        dict -- dictionary of settings
        Settings() -- cpf Settings class 
    pattern : int, str, optional
        Which patterns in the series to plot. 
        Either "all", "mid", or integer list of patterns to keep.
        The default is 0.
    subpattern : int, optional
        Which subpattern in the series to plot. 
        The default is "all".
    report : Literal[ "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR" ], optional
        Logger level for how much information to write to the log files.
        The default is "INFO".
    **kwargs : key, value pairs
        key, value arguments arguments. Passed though the method to cpf.XRD_Fitpattern.initiate() and 
        cpf.XRD_Fitpattern.execute(). See documention of these for kwarg use.
 
    """

    # pass everything through to XRD_FitPattern.set_range
    xrd_fp_set_range(*args, **kwargs)


def write_output(
    settings,
    out_type: [str | list] = None,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs,
):
    """
    Write output types for either list of out_type or the list in the settings.
    
    Parameters
    ----------
    settings : Optional[str | Path | dict | Settings()]
       Pointer to information needed for settings class. Can be of the form:        
       string -- filename of python formatted file 
       Path -- path for python formatted file  
       dict -- dictionary of settings
       Settings() -- cpf Settings class 
    out_type : [str, list], optional
        String or list of outputs to be processed. If absent then outputs from settings are used. 
        The default is None (and outputs from settings are used).
    report : Literal[        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"    ], optional
        Logger level for how much information to write to the log files.
        The default is "INFO".
    **kwargs : key, value pairs
        key, value arguments arguments. Passed though the method to cpf.XRD_Fitpattern.initiate() and 
        called output methods. See documention of these for kwarg use.

    Returns
    -------
    None.

    """
    
    settings_class = initiate(settings, report=report, **kwargs)
    # make a note in the logger.
    # suppress output if called by another module.
    for i in range(len(inspect.stack()) - 1, -1, -1):
        if inspect.stack()[i].function == "<module>":
            base_call = inspect.stack()[i - 1].function
    if base_call == "write_output":
        logger.info("")
        logger.info(
            f"Running: Cascade.write_output with settings: {settings_class.settings_file}"
        )
        logger.info("")

    if out_type is not None:
        logger.moreinfo(f"Output_type was provided as an option; will use {out_type}")
        settings_class.set_output_types(out_type_list=out_type)

    if settings_class.output_types is None:
        logger.warning(
            "There are no output types. Add 'Output_type' to input file or specify 'out_type' in command."
        )
    else:
        for mod in settings_class.output_types:
            if mod in spot_output_methods_modules:
                logger.info(" ".join(map(str, [("Writing output file(s) using %s" % mod)])))
                wr = spot_output_methods_modules[mod]
                wr.WriteOutput(
                    settings_class,
                    **kwargs,
                )


def execute(
    settings,#: [str | Path | dict | Settings()],
    subpattern: str = "all",
    parallel: bool = True,
    resume: bool = False,
    mode: str = "cascade",
    # show_plots: bool = False,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    spot_finding_method = "scipy_peakfind",
    **kwargs,
):
    """
    Runs the spot finding algorithms on the diffraction data. As a method for constraining grain-size. 
    
    All parameters the affect the fitting are determined from the settings class (formed by settings and cpf.XRD_FitPattern.initiate())
    The arguments input directly into here only affect how the code is run here.
    For example, if the code is run in parallel or not. Other kwargs are passed though but do not affect the running of the code.    
    
    Parameters
    ----------
    settings : Optional[str | Path | dict | Settings()]
        Pointer to information needed for settings class. Can be of the form:        
        string -- filename of python formatted file 
        Path -- path for python formatted file  
        dict -- dictionary of settings
        Settings() -- cpf Settings class 
    subpattern : int, optional
        Which subpattern in the series to plot. 
        The default is "all".
    parallel : bool, optional
        Turns parallel processing on (if True) or off (if False). The default is True.
    resume : bool, optional
        Resume processing the fits from last completed (if True) or 
        from the begining (if False). The default is False.
    mode : str
        How to process the provided data, changes path through the method.
        Options are:
        - view -- plots data. mode set by XRD_FitPattern.view()
        - set-range -- plots data for ranges. mode set by XRD_FitPattern.set_range()
        - set_guess -- used for setting inital position guesses for peaks. mode set by XRD_FitPattern.initial_peak_position()
        - search -- used for searching over orders and series types. mode set by XRD_FitPattern.order_search()
        - fit -- fits the data.             
        The default is "fit".
    report : Literal[ "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"    ], optional
        Logger level for how much information to write to the log files.
        The default is "INFO".
    **kwargs : key, value pairs
        key, value arguments arguments. Passed though to called methods. Not used here.

    """
    
    settings_class = initiate(settings, report=report, **kwargs)
    # make note in logger
    # suppress output if called by another module.
    for i in range(len(inspect.stack()) - 1, -1, -1):
        if inspect.stack()[i].function == "<module>":
            base_call = inspect.stack()[i - 1].function
    if base_call == "execute":
        logger.info("")
        logger.info(
            f"Running: Cascade.execute with settings: {settings_class.settings_file}"
        )
        logger.info("")
    
    # parse values that affect the fitting from settings (if they are set). 
    if "fit_options" not in settings_class.__dict__:
        settings_class.fit_options = {}
    as_masked = settings_class.fit_options.get('as_masked', False)
    if (mode == "set-range" or mode == "view"):
        as_masked = True
    elif (mode == "fit" or mode == "search"):
        as_masked = as_masked
        if as_masked == True:
            logger.warning("'as_masked'==True changes the fit for some masked datasts. I dont know why. Check fits with and without this setting")
    else:
        logger.critical(f"Unknown mode '{mode}'.")
        
    # get spot finding method
    wr = spot_finding_methods_modules[spot_finding_method]
        
    #get data from settings class
    new_data = settings_class.data_class

    # Define locally required names
    # temporary_data_file = make_outfile_name(
    #     "PreviousChunkFit_JSON",
    #     directory=settings_class.output_directory,
    #     extension=".dat",
    #     overwrite=True,
    # )

    data_to_fill: Path
    if settings_class.calibration_data:
        data_to_fill = settings_class.calibration_data.resolve()
    else:
        data_to_fill = settings_class.image_list[0]
        if not isinstance(data_to_fill, Path):
            data_to_fill = Path(data_to_fill)

    new_data.fill_data(
        data_to_fill,
        settings=settings_class,
        # report=report
    )

    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)
    
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
    else:
        pool = None

    # all_spots = pd.DataFrame()
    
    # Process the diffraction patterns
    # for j in range(settings_class.image_number):
    progress = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger
    for j in progress.iter_bar(iteration=range(settings_class.image_number)):
        logger.info(
            f"Processing {title_file_names(image_name=settings_class.image_list[j])}"
        )
        
        # Get diffraction pattern to process.
        settings_class.set_subpattern(j, 0)
        new_data.import_image(settings=settings_class)#, debug=debug)
        # get metadata
        metadata = new_data.get_metadata(settings_class=settings_class)
        
        if (
            isinstance(settings_class.datafile_preprocess, dict)
            or (
                isinstance(settings_class.calibration_mask, dict)
                and "threshold" in settings_class.calibration_mask
            )
        ):
            # needed because image preprocessing adds to the mask and is different for each image.
            # set intenstiy threshold and/or cosmics for each frame.
            # only applies these because all other mask functions are static.
            new_data.mask_restore()
            if (
                isinstance(settings_class.datafile_preprocess, dict)
                and "cosmics" in settings_class.datafile_preprocess
            ):
                new_data = cosmicsimage_preprocess(new_data, settings_class)
            if (
                isinstance(settings_class.calibration_mask, dict)
                and "threshold" in settings_class.calibration_mask
            ):
                # set intenstiy threshold for each frame.
                # only applies to threshold because all other mask functions are static
                # (cannot change between frames)
                new_data.set_mask(
                    intensity_bounds=settings_class.calibration_mask["threshold"]
                )
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
                    additional_text="integrated",
                    extension=".png",
                    overwrite=True,
                )

                fig.savefig(filename)

                print("Plotted data.")
            else:
                plt.close()

        """# Get previous fit (if it exists and is required)
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
                del previous_fit"""

        # Switch to save the first fit in each sequence.
        save_figs = False # True if (j == 0 or save_all is True) else False

        # call image processing function and do.
        logger.info(f"Finding spots in patterns using {spot_finding_method}")
        
        #perform the actual function here. 
        all_spots = wr.spot_find(
            settings_class,
            new_data,
            parallel_pool = pool,
            **kwargs,
        )
        
        # write output files
        # store the fit parameters' information as a JSON file.
        filename_to_write = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            additional_text="spots",#settings_class.file_label,
            extension=".json",
            overwrite=True,
        )
        WriteFits(settings_class, all_spots, filename_to_write=filename_to_write, metadata=metadata, mode=mode)

    if mode == "cascade":
        # Write the output files.
        write_output(settings_class, **kwargs)

    if parallel is True:
        pool.clear()


    # # plot the fits
    # wr.write_grains_list(
    #     inputs=settings_class, report=report, subpattern="all", **kwargs
    # )

    # plot_peak_count(inputs=settings_class, report=report, subpattern="all", **kwargs)

        # # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
        # all_fitted_chunks = []
        # all_chunk_positions = []
        # parallel_pile = []

        # for i in range(len(settings_class.fit_orders)):
        #     if parallel == False:
        #         serial_string = f"Fitting range {i+1}/{len(settings_class.fit_orders)}"
        #         logger.info(serial_string)

        #     # get settings for current subpattern
        #     settings_class.set_subpattern(j, i)

        #     """if "previous_fit" in locals() and mode == "fit":
        #         params = previous_fit[i]
        #         params.pop("correlation_coeffs", None)
        #     else:
        #         params = None

        #     # Track the position of the peak centroid
        #     # FIXME: This is crude - the range doesn't change width. so can't account
        #     # for massive change in stress.
        #     # But does it need to?
        #     tth_range = np.array(settings_class.subfit_orders["range"])
        #     if settings_class.fit_track is True and "previous_fit" in locals():
        #         null_terms = has_value(params, val=None)
        #         if null_terms == True:
        #             # the previous fit has problems so discard it
        #             logger.moreinfo(  # type: ignore
        #                 "Tracking peak centre but propagated fit has problems. Not sensible to track the centre of the fit for this step."
        #             )
        #             params = []
        #         else:
        #             # if tacking things start from positon of previous fit positions
        #             tth_range = previous_fit[i]["range"][0]
        #             mid = []
        #             for k in range(len(params["peak"])):
        #                 mid.append(get_series_mean(params["peak"][k], "d-space"))
            
        #             if mid == None or mid == 0:
        #                 # the previous fits failed in some way.
        #                 move_by = 0
        #             else:
        #                 cent = new_data.conversion(np.mean(mid), reverse=True)
        #                 move_by = cent - np.mean(tth_range)
            
        #             # update tth_range and settings
        #             tth_range = tth_range + move_by
        #             settings_class.fit_orders[i]["range"] = (
        #                 settings_class.fit_orders[i]["range"] + move_by
        #             )
            
        #             logger.moreinfo(
        #                 f"Move range for fitting. \n Initial range: [{tth_range[0]:4.2f},{tth_range[1]:4.2f}]; will be moved by {move_by:4.2f}; the new range is [{tth_range[0]+move_by:4.2f},{tth_range[1]+move_by:4.2f}]"
        #             )
            
        #             # The PeakPositionSelections are only used if the fits are not being propagated
        #             if "PeakPositionSelection" in settings_class.fit_orders[i]:
        #                 for k in range(
        #                     len(settings_class.fit_orders[i]["PeakPositionSelection"])
        #                 ):
        #                     settings_class.fit_orders[i]["PeakPositionSelection"][k][
        #                         2
        #                     ] = (
        #                         settings_class.fit_orders[i]["PeakPositionSelection"][
        #                             k
        #                         ][2]
        #                         + move_by
        #                     )
            
        #             # re-get settings for current subpattern
        #             settings_class.set_subpattern(j, i)"""
        #     tth_range = settings_class.subfit_orders["range"]

        #     sub_data = new_data.duplicate_without_detector(
        #         range_bounds=tth_range, as_masked=as_masked
        #     )

        #     # Mask the subpattern by intensity if called for
        #     if (
        #         "imax" in settings_class.subfit_orders
        #         or "imin" in settings_class.subfit_orders
        #     ):
        #         sub_data = SpotProcess(sub_data, settings_class, as_masked=as_masked)

        #     if mode == "set-range":
        #         fig_1 = plt.figure()
        #         sub_data.plot_masked(fig_plot=fig_1, **kwargs)
        #         plt.suptitle(peak_string(settings_class.subfit_orders) + "; masking")

        #         filename = make_outfile_name(
        #             settings_class.image_list[j],
        #             directory=settings_class.output_directory,
        #             additional_text="mask",
        #             orders=settings_class.subfit_orders,
        #             extension=".png",
        #             overwrite=True,
        #         )

        #         fig_1.savefig(filename)

        #         # if debug:
        #         plt.show()
        #         plt.close()

        #     elif mode == "view":
        #         fig = plt.figure()
        #         ax = fig.add_subplot(1, 1, 1)
        #         ax_o1 = plt.subplot(111)
        #         sub_data.plot_calibrated(
        #             fig_plot=fig,
        #             axis_plot=ax,
        #             show="intensity",
        #             **kwargs,  # rastered="scatter"
        #         )
        #         plt.suptitle(peak_string(settings_class.subfit_orders) + "; calibrated")

        #         filename = make_outfile_name(
        #             settings_class.image_list[j],
        #             directory=settings_class.output_directory,
        #             additional_text="range",
        #             orders=settings_class.subfit_orders,
        #             extension=".png",
        #             overwrite=True,
        #         )

        #         fig.savefig(filename)

        #         plt.show()
        #         if mode == "view":
        #             print("Plotted data.")
        #         else:
        #             plt.close()
            
        #     else: 
        #         if parallel is True:  # setup parallel version
        #             kwargs = {
        #                 "save_fit": save_figs,
        #                 "mode": mode,
        #                 # "histogram_type": settings_class.cascade_histogram_type,
        #                 # "histogram_bins": settings_class.cascade_histogram_bins,
        #                 **kwargs
        #             }
        #             arg = (sub_data, settings_class.duplicate())
        #             parallel_pile.append((arg, kwargs))

        #         else:  # non-parallel version
        #             tmp = fit_sub_pattern(
        #                 sub_data,
        #                 settings_class,  # added
        #                 None,  # do not pass params they are not needed.
        #                 save_fit=save_figs,
        #                 mode=mode,
        #                 # histogram_type=settings_class.cascade_histogram_type,
        #                 # histogram_bins=settings_class.cascade_histogram_bins,
        #                 **kwargs
        #             )
        #             all_fitted_chunks.append(tmp[0])
        #             all_chunk_positions.append(tmp[1])

        # # write output files
        # if mode != "set-range":
        #     if parallel is True:
        #         tmp = pool.map(parallel_processing, parallel_pile)
        #         for i in range(len(settings_class.fit_orders)):
        #             # fitted_param.append(tmp[i][0])
        #             # lmfit_models.append(tmp[i][1])
        #             all_fitted_chunks.append(tmp[i][0])
        #             all_chunk_positions.append(tmp[i][1])

        #         # tmp = fit_sub_pattern(
        #         #     sub_data,
        #         #     settings_class,  # added
        #         #     None,  # do not pass params they are not needed.
        #         #     save_fit=save_figs,
        #         #     debug=debug,
        #         #     mode=mode,
        #         #     histogram_type = settings_class.cascade_histogram_type,
        #         #     histogram_bins = settings_class.cascade_histogram_bins,
        #         # )
        #         # all_fitted_chunks.append(tmp[0])
        #         # all_chunk_positions.append(tmp[1])

        #     # paste the fits to an output file.
        #     # store the chunk fits' information as a JSON file.
        #     filename = make_outfile_name(
        #         settings_class.subfit_filename,
        #         directory=settings_class.output_directory,
        #         additional_text="chunks",
        #         # orders=settings_class.subfit_orders,
        #         extension=".json",
        #         overwrite=True,
        #     )
        #     with open(filename, "w") as TempFile:
        #         # Write a JSON string into the file.
        #         json.dump(
        #             (all_fitted_chunks, all_chunk_positions),
        #             TempFile,
        #             sort_keys=True,
        #             indent=2,
        #             default=numpy_to_json,
        #         )

        #     # if propagating the fits write them to a temporary file
        #     if settings_class.fit_propagate:
        #         # print json_string
        #         with open(temporary_data_file, "w") as TempFile:
        #             # Write a JSON string into the file.
        #             json.dump(
        #                 (all_fitted_chunks, all_chunk_positions),
        #                 TempFile,
        #                 sort_keys=True,
        #                 indent=2,
        #                 default=numpy_to_json,
        #             )

    # if parallel is True:
    #     pool.close()

    # # plot the fits
    # plot_cascade_chunks(
    #     inputs=settings_class, report=report, subpattern="all", **kwargs
    # )

    # plot_peak_count(inputs=settings_class, report=report, subpattern="all", **kwargs)
