
import os as _os
from importlib import import_module as _import_module
from types import ModuleType as _ModuleType

import pandas as pd
import matplotlib.pyplot as plt
import proglog
import numpy as np
from cpf.util.io import (
    has_value,
    make_outfile_name,
    numpy_to_json,
    peak_string,
    title_file_names,
    peak_hkl
)
from pathlib import Path

from scipy.signal import find_peaks
from typing import Literal, Optional
import json
from cpf.XRD_FitSubpattern import fit_sub_pattern
from cpf.fitsubpattern_chunks import fit_chunks
from cpf.BrightSpots import SpotProcess

from cpf.util.logging import get_logger

logger = get_logger("cpf.pattern_processing_methods.spots_scipy_peakfind")



def method_defaults():
    """List non-universally required parameters for this image processing method."""
    # required parameters -- ones that cannot be guessed at in advance or given a pre-existing value
    # optional parameters -- values that can be guessed at in advance

    required_params = [
        #'apparently none!
    ]
    optional_params = {
        "waterfall_type": "range",
        "prominence": 15,
        "width": 0
    }

    return required_params, optional_params



def spot_find(
        settings_class,
        data_class,
        mode="cascade",
        resume = False,
        parallel_pool = None,
        subpattern = "all",
        # parallel,
        # pool,
        # as_masked,
        report: Literal[
            "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
        ] = "INFO",
        **kwargs
        ):
    """
    Calculate spot positions using scipy peak find. 
    
    This algorithm fits peaks to chunks of the data around 

    Parameters
    ----------
    settings_class : TYPE
        DESCRIPTION.
    data_class : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    #parse inputs 
    if "spot_find_settings" in settings_class.__dict__:
        prominence = settings_class.spot_find_settings.get("prominence", method_defaults()[1]["prominence"])
        width = settings_class.spot_find_settings.get("width", method_defaults()[1]["width"])
    else:
        prominence = method_defaults()[1]["prominence"]
        width = method_defaults()[1]["width"]
    as_masked = settings_class.fit_options.get('as_masked', False)
    # the followshing should be in 'settings_class.spot_find_settings' but this is not implemented yet.
    waterfall_type = settings_class.fit_options.get("waterfall_type", method_defaults()[1]["waterfall_type"])
        
    
    # restrict to patterns listed
    # settings_class.set_data_files(keep=pattern)
    # -- not needed bacuse pattern is already set in cascade.execute
    # restrict to sub-patterns listed
    settings_class.set_subpatterns(subpatterns=subpattern)

    num_orders = len(settings_class.fit_orders)
    
    if not resume:
        all_fitted_chunks, all_chunk_positions = waterfall_intensities(
                settings_class=settings_class,
                data_class=data_class,
                mode=mode,
                parallel_pool=parallel_pool,
                **kwargs
                )
        
        # store the chunk fits' information as a JSON file.
        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            additional_text="chunks",
            # orders=settings_class.subfit_orders,
            extension=".json",
            overwrite=True,
        )
        print("filename", filename)
        write_intermediate_files(filename, (all_fitted_chunks, all_chunk_positions))
    else:
        all_fitted_chunks, all_chunk_positions = read_intermediate_files(
            settings_class,
            report = report,
            **kwargs,
        )

    # spots_via_scipy(
    #     settings_class,
    #     **kwargs)
    
    # all_azis, all_data = read_intermediate_files(
    #     settings_class, subpattern=subpattern
    # )
    
    #turn output from spots_via_scipy() into a panda DataFrame
    # set up data array
    # headers_use =  ["image_position",
    #                "DataFile",
    #                'phase', 
    #                'peak',   
    #                ]
    # headers_use += settings_class.metadata
    headers_use = ["2theta",
                    "2theta_err",  
                    '2theta extent',
                    '2theta extent_err',
                    'azimuth', 
                    'azimuth_err',
                    'azimuth_extent',
                    'azimuth_extent_err',
                    'intensity',
                    'intensity_err'
                   ]

    all_azis, all_data = all_chunk_positions, all_fitted_chunks

    spot_array = []
    # for ptrn in range(len(settings_class.image_list)):
    for ordr in range(num_orders):
        settings_class.set_subpattern(settings_class.subfit_filename_position, ordr)
        for pk in range(len(settings_class.subfit_orders["peak"])):
            
            spots, properties = find_peaks(
                all_data[ordr]["h"][pk], prominence=prominence, width=width
            )
                
            for pos in range(len(spots)):
                spts = spots[pos]
                # spot_array_row = [
                #     _os.path.basename(settings_class.subfit_filename), 
                #     settings_class.subfit_orders['peak'][pk]['phase'],
                #     peak_hkl(settings_class.subfit_orders)[pk]]
                # if settings_class.metadata:
                #     spot_array_row += list(data_class.get_metadata(settings_class).values())
                if len(all_data[ordr]["d"][pk]) != 0:
                    spot_array_row = [
                        settings_class.data_class.conversion(np.array(all_data[ordr]["d"][pk])[spts]), #tth
                        settings_class.data_class.conversion(np.array(all_data[ordr]["d_err"][pk])[spts]), #tth error
                        np.array(all_data[ordr]["w"][pk])[spts], #width in tth
                        np.array(all_data[ordr]["w_err"][pk])[spts], #width in tth error
                        np.array(all_azis[ordr])[spts],#azimuth values
                        np.array(all_azis[ordr])[spts]*0 ,#azimuth error values
                        properties['widths'][pos],  # azimith extent
                        properties['widths'][pos]*0,
                        np.array(all_data[ordr]["h"][pk])[spts], # height
                        np.array(all_data[ordr]["h_err"][pk])[spts], #height error
                    ]
                else:
                    spot_array_row = [
                        np.mean(settings_class.fit_orders[ordr]["range"]), #middle of tth range
                        settings_class.fit_orders[ordr]["range"][1]-settings_class.fit_orders[ordr]["range"][0], #tth error: width of range
                        0,
                        0,
                        np.array(all_azis[ordr])[spts],#azimuth values
                        np.array(all_azis[ordr])[spts]*0 ,#azimuth error values
                        properties['widths'][pos], # azimith extent
                        properties['widths'][pos]*0,                        
                        np.array(all_data[ordr]["h"][pk])[spts], # height
                        0#np.array(all_data[ptrn][ordr]["h_err"][pk])[spts], #height error
                    ]

                spot_array.append(spot_array_row)
    
    spot_table = pd.DataFrame(spot_array, columns=headers_use)
    
    # write_grains_list(filename, save=True)
    
    # plot_grains(filename, save=True)
        
    # return spot_array, headers_use
    return spot_table
    
    
def waterfall_intensities(
        settings_class,
        data_class,
        mode="cascade",
        parallel_pool = None,
        **kwargs
        ):
    """
    Calculates azimuthal intensities for the fit_ranges defined in settings.
    The number and type of bins used to do this are defined in settings.
    
    FIXME: document names of bin options properly.
    
    Possible values for waterfall_type are:
        'cascade',
        'range'
        'maxima'
        '98thPercentile'
        [not implemented but should be when call correct function]


    Parameters
    ----------
    settings_class : TYPE
        DESCRIPTION.
    data_class : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    
    # parse inputs. 
    as_masked = settings_class.fit_options.get('as_masked', False)
    # the followshing should be in 'settings_class.spot_find_settings' but this is not implemented yet.
    waterfall_type = settings_class.fit_options.get("waterfall_type", method_defaults()[1]["waterfall_type"])
    
    if parallel_pool:
        parallel = True
        pool = parallel_pool
    else:
        parallel = False
    
    #get metadata
    metadata = data_class.get_metadata(settings_class)
    
    # Switch to save the first fit in each sequence.
    save_figs = False # True if (j == 0 or save_all is True) else False

    # Pass each sub-pattern to Fit_Subpattern for fitting in turn.
    all_fitted_chunks = []
    all_chunk_positions = []
    parallel_pile = []
    for i in range(len(settings_class.fit_orders)):
        if parallel == False:
            serial_string = f"Fitting range {i+1}/{len(settings_class.fit_orders)}"
            logger.info(serial_string)

        # get settings for current subpattern
        settings_class.set_subpattern(settings_class.subfit_filename_position, i)
        
        print(settings_class.subfit_filename_position, i, settings_class.subfit_filename)


        """if "previous_fit" in locals() and mode == "fit":
            params = previous_fit[i]
            params.pop("correlation_coeffs", None)
        else:
            params = None

        # Track the position of the peak centroid
        # FIXME: This is crude - the range doesn't change width. so can't account
        # for massive change in stress.
        # But does it need to?
        tth_range = np.array(settings_class.subfit_orders["range"])
        if settings_class.fit_track is True and "previous_fit" in locals():
            null_terms = has_value(params, val=None)
            if null_terms == True:
                # the previous fit has problems so discard it
                logger.moreinfo(  # type: ignore
                    "Tracking peak centre but propagated fit has problems. Not sensible to track the centre of the fit for this step."
                )
                params = []
            else:
                # if tacking things start from positon of previous fit positions
                tth_range = previous_fit[i]["range"][0]
                mid = []
                for k in range(len(params["peak"])):
                    mid.append(get_series_mean(params["peak"][k], "d-space"))
        
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
        
                logger.moreinfo(
                    f"Move range for fitting. \n Initial range: [{tth_range[0]:4.2f},{tth_range[1]:4.2f}]; will be moved by {move_by:4.2f}; the new range is [{tth_range[0]+move_by:4.2f},{tth_range[1]+move_by:4.2f}]"
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
                settings_class.set_subpattern(j, i)"""
        tth_range = settings_class.subfit_orders["range"]

        sub_data = data_class.duplicate_without_detector(
            range_bounds=tth_range, as_masked=as_masked
        )
        # Mask the subpattern by intensity if called for
        if (
            "imax" in settings_class.subfit_orders
            or "imin" in settings_class.subfit_orders
        ):
            sub_data = SpotProcess(sub_data, settings_class, as_masked=as_masked)

        if parallel is True:  # setup parallel version
            kwargs = {
                "save_fit": save_figs,
                "mode": waterfall_type,
                # "histogram_type": settings_class.cascade_histogram_type,
                # "histogram_bins": settings_class.cascade_histogram_bins,
                **kwargs
            }
            arg = (sub_data, settings_class.duplicate_without_dataclass())
            parallel_pile.append((arg, kwargs))

        else:  # non-parallel version
        
            # *** SHOULD USE fit_chunks DIRECTLY HERE ***
            # *** THERE IS NO NEED TO CALL fit_sub_pattern ***
            tmp = fit_chunks(
                sub_data,
                settings_class,
                mode=waterfall_type,
                # histogram_type=histogram_type,
                # histogram_bins=histogram_bins,
                # fit_method=fit_method,
            )
            """
            tmp = fit_sub_pattern(
                sub_data,
                settings_class.duplicate_without_dataclass(),  # added
                None,  # do not pass params they are not needed.
                save_fit=save_figs,
                mode=mode,
                # histogram_type=settings_class.cascade_histogram_type,
                # histogram_bins=settings_class.cascade_histogram_bins,
                **kwargs
            )
            """
            all_fitted_chunks.append(tmp[0])
            all_chunk_positions.append(tmp[1])

    # write output files
    if parallel is True:
        tmp = pool.map(parallel_processing, parallel_pile)
        for i in range(len(settings_class.fit_orders)):
            all_fitted_chunks.append(tmp[i][0])
            all_chunk_positions.append(tmp[i][1])
    
    return all_fitted_chunks, all_chunk_positions


def parallel_processing(p):
    a, kw = p
    # return fit_sub_pattern(*a, **kw)
    return fit_chunks(*a, **kw)




def write_intermediate_files(filename, fits_to_write, file_type="lists"):
    """
    Write files for the the chunk fits. 

    These can either be 
        "lists" -- compound lists and dictionaries of parameters, or 
        "tabulated" -- the same data arranged in tables.     
    
    Parameters
    ----------
    filename : str
        Filename to write.
    fits_to_write : 
        The chunk fits to write to the files.
    file_type : str
        either "lists" or "tabulated". Default is "lists"

    Returns
    -------
    None.

    """
    
    all_fitted_chunks = fits_to_write[0]
    all_chunk_positions = fits_to_write[1]
        
    if file_type.lower()=="lists":
        with open(filename, "w") as TempFile:
            # Write a JSON string into the file.
            json.dump(
                (all_fitted_chunks, all_chunk_positions),
                TempFile,
                sort_keys=True,
                indent=2,
                default=numpy_to_json,
            )
    elif file_type.lower()=="tabulated":
        #write everything as data frame --> csv
        raise NotImplementedError("saving chunks as tabulated data is not implemented")
        # this doesnt actually work. Needs fixing. 
        
        import pandas as pd 
        
        for fits in all_fitted_chunks:
           
            for peaks in range(len(fits['d'])):
            # fits = all_fitted_chunks[i]
                
                fits['peak'][peaks] = [fits['peak'][peaks]] * len(fits['d'][peaks])

            fits_pd = pd.DataFrame.from_dict(fits)
        [fits_to_write[0][3]['peak']]* 4
        pd.DataFrame.from_dict(fits_to_write[0][1])
    else:
        raise ValueError(f"The file_type {file_type} is not recognised or implemented.")


def read_intermediate_files(
    settings_class,
    file_type = "lists",
    pattern = "all",
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
    **kwargs,
):
    """
    Reads the intermediate chunk files.


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

    # if inputs:
    #     settings_class = inputs
    # elif settings_class is None:
    #     settings_class = initiate(settings_file, report=report)
    # else:
    #     settings_class = settings_class

    all_azis = []
    all_fits = []


    if file_type.lower()=="lists":

        print("Reading chunk files")
        lgr = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger
    
        # for f in range(setting_class.image_number):
        for f in lgr.iter_bar(iteration=range(settings_class.image_number)):
            settings_class.set_subpattern(f, 0)
            filename = make_outfile_name(
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
    
    elif file_type.lower()=="tabulated":
        #write everything as data frame --> csv
        raise NotImplementedError("reading chunks as tabulated data is not implemented")
        # this doesnt actually work. Needs fixing. 
    else:
        raise ValueError(f"The file_type {file_type} is not recognised or implemented.")

    
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



# # def write_grains_list(
# def spots_via_scipy(
#     settings_class,
#     pattern: [str | list] = "all",
#     subpattern: [str | list] = "all",
#     report: Literal[
#         "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
#     ] = "INFO",
    
#     ):
#     """
#     Find spots in the azimuthal waterfall data and return a list of positions. 
    
#     Returns dataframe or table/list of lists/something with columns:
#         filename,
#         (metadata)
#         possible peak, (is this peak in the two theta range of a peak)
#         tth,
#         tth_err,
#         tth_width,
#         tth_width_err,
#         azm,
#         azm_err,
#         axm_width,
#         azm_width_err,
#         intensity, (height of maximum abve background as determined by algorithm)
#         intensity_err,
#         total counts, 
#         total_counts_err,
#         + extra columns depending on method used.

#     Parameters
#     ----------
#     settings : TYPE
#         DESCRIPTION.
#      : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     """

#     #parse inputs 
#     if "spot_find_settings" in settings_class.__dict__:
#         prominence = settings_class.spot_find_settings.get("prominence", method_defaults()[1]["prominence"])
#         width = settings_class.spot_find_settings.get("width", method_defaults()[1]["width"])
#     else:
#         prominence = method_defaults()[1]["prominence"]
#         width = method_defaults()[1]["width"]
        

#     # restrict to patterns listed
#     # settings_class.set_data_files(keep=pattern)
#     # -- not needed bacuse pattern is already set in cascade.execute
#     # restrict to sub-patterns listed
#     settings_class.set_subpatterns(subpatterns=subpattern)

#     num_orders = len(settings_class.fit_orders)

#     all_azis, all_data = read_intermediate_files(
#         settings_class, subpattern=subpattern
#     )

#     all_peaks: list[list] = []
#     all_peakAzis: list = []
#     all_properties: list[list] = []
#     count = []
#     peak_labels = []
    
#     for ptrn in range(len(settings_class.image_list)):
        
#         for ordr in range(num_orders):
#             # loop over the number of sets of peaks fit for (i.e. len(settings_class.fit_orders))
    
#             settings_class.set_subpattern(ptrn, ordr)
    
#             for pk in range(len(settings_class.subfit_orders["peak"])):
#                 # loop over the number of peaks in each fit_orders
    
#                 all_peaks_tmp = []
#                 all_properties_tmp = []
#                 count_tmp = []
#                 if all_data[0][ordr]["h"][pk]:
#                     # if there is some data in the array plot it.
#                     # for i in range(settings_class.image_number):
                        
#                     spots, properties = find_peaks(
#                         all_data[ptrn][ordr]["h"][pk], prominence=prominence, width=width
#                     )
#                     properties["PeakAzis"] = np.array(all_azis[ptrn][ordr])[spots]
#                     # properties["PeakTths"] = np.array(all_azis[i][ordr])[spots]
#                     all_peaks_tmp.append(spots)
#                     all_peakAzis.append(np.array(all_azis[ptrn][ordr])[spots])
#                     all_properties_tmp.append(properties)
#                     count_tmp.append(len(spots))
#                     if 1 and ptrn == 0:
#                         fig = plt.figure()
#                         plt.plot(all_azis[ptrn][ordr], all_data[ptrn][ordr]["h"][pk])
#                         plt.plot(np.array(all_azis[ptrn][ordr])[spots], np.array(all_data[ptrn][ordr]["h"][pk])[spots], "x")

#                         # plt.plot(np.zeros_like(x), "--", color="gray")
#                         show_plots = True
#                         if show_plots is True:
#                             plt.show()
#                         else:
#                             plt.close()

#                 all_peaks.append(all_peaks_tmp)
#                 all_properties.append(all_properties_tmp)
#                 count.append(count_tmp)
#                 # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
#                 # pk: int | Literal["all"] = k
#                 for l in range(len(settings_class.subfit_orders["peak"])):
#                     if not all_data[0][ordr]["h"][l]:
#                         pk = "all"
#                 # make the figure title
#                 ttlstr = peak_string(settings_class.subfit_orders, peak=pk)
#                 peak_labels.append(ttlstr)

#     # return all_peaks, all_properties, count, peak_labels
    
#     #turn output from spots_via_scipy() into a panda DataFrame
#     # set up data array
#     headers_use =  ["image_position",
#                    "DataFile",
#                    'phase', 
#                    'peak',   
#                    ]
#     headers_use += settings_class.metadata
#     headers_use += ["2theta",
#                     "2theta_err",  
#                     '2theta extent',
#                     '2theta extent_err',
#                     'azimuth', 
#                     'azimuth_err',
#                     'azimuth_extent',
#                     'azimuth_extent_err',
#                     'intensity'
#                     'intensity_err'
#                    ]

    
#     spot_array = []
#     for ptrn in range(len(settings_class.image_list)):
#         for ordr in range(num_orders):
#             settings_class.set_subpattern(ptrn, ordr)
#             for pk in range(len(settings_class.subfit_orders["peak"])):
                
#                 spots, properties = find_peaks(
#                     all_data[ptrn][ordr]["h"][pk], prominence=prominence, width=width
#                 )
                    
#                 for spts in spots:
#                     spot_array_row = [
#                         _os.path.basename(settings_class.subfit_filename), 
#                         settings_class.subfit_orders['peak'][pk]['phase'],
#                         peak_hkl(settings_class.subfit_orders)[pk]]
                    
#                     for md in settings_class.metadata:
#                         spot_array_row += [
                            
#                             ]
                    
#                     if len(all_data[ptrn][ordr]["d"][pk]) != 0:
#                         spot_array_row += [
#                             settings_class.data_class.conversion(np.array(all_data[ptrn][ordr]["d"][pk])[spts]), #tth
#                             settings_class.data_class.conversion(np.array(all_data[ptrn][ordr]["d_err"][pk])[spts]), #tth error
#                             np.array(all_data[ptrn][ordr]["w"][pk])[spts], #width in tth
#                             np.array(all_data[ptrn][ordr]["w_err"][pk])[spts], #width in tth error
#                             np.array(all_azis[ptrn][ordr])[spts],#azimuth values
#                             np.array(all_azis[ptrn][ordr])[spts]*0 ,#azimuth error values
#                             np.array(all_data[ptrn][ordr]["h"][pk])[spts], # height
#                             np.array(all_data[ptrn][ordr]["h_err"][pk])[spts], #height error
#                         ]
#                     else:
#                         spot_array_row += [
#                             np.mean(settings_class.fit_orders[ordr]["range"]), #middle of tth range
#                             settings_class.fit_orders[ordr]["range"][1]-settings_class.fit_orders[ordr]["range"][0], #tth error: width of range
#                             0,
#                             0,
#                             np.array(all_azis[ptrn][ordr])[spts],#azimuth values
#                             np.array(all_azis[ptrn][ordr])[spts]*0 ,#azimuth error values
#                             np.array(all_data[ptrn][ordr]["h"][pk])[spts], # height
#                             0#np.array(all_data[ptrn][ordr]["h_err"][pk])[spts], #height error
#                         ]

#                     spot_array.append(spot_array_row)
        
#     import pandas as pd
#     spot_table = pd.DataFrame(spot_array, columns=headers_use)
#     stop
#     return spot_table
        
# %% 
# =============================================================================
# HERE AFTER ARE FUNCIOTNS THAT SHOULD BE MOEVED INTO OUTPUTS
# =============================================================================
        
def peak_count(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug: bool = False,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
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
        settings_class = initiate(settings_file, inputs=inputs, report=report)
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
                ttlstr = peak_string(settings_class.subfit_orders, peak=pk)
                peak_labels.append(ttlstr)

    return all_peaks, all_properties, count, peak_labels


def plot_peak_count(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug: bool = False,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
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
        settings_class = initiate(settings_file, inputs=inputs, report=report)
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
    filename = make_outfile_name(
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
        
        


def plot_cascade_chunks(
    settings_file=None,
    settings_class=None,
    inputs=None,
    debug=False,
    report: Literal[
        "DEBUG", "EFFUSIVE", "MOREINFO", "INFO", "WARNING", "ERROR"
    ] = "INFO",
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
        settings_class = initiate(settings_file, inputs=inputs, report=report)
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
        if "min_all" in locals():
            vmin = min_all[j]
        if "max_all" in locals():
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
                    + peak_string(settings_class.fit_orders[j], peak=k)
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
                        if 1:  # len(all_data[i][j]["h"]) == 1:
                            tms = modified_time_s[j]
                        else:
                            tms = modified_time_s[i] * np.ones(
                                np.shape(all_data[i][j]["chunks"])
                            )

                        plt.scatter(
                            all_data[i][j]["chunks"],
                            modified_time_s[i]
                            * np.ones(np.shape(all_data[i][j]["chunks"])),
                            s=10,  # s=.05,
                            c=(all_data[i][j]["h"][k]),
                            vmax=vmax,
                            vmin=vmin,
                            # cmap="YlOrBr",
                            norm=norm,
                        )

                    # determine the label for the figure -- if there is data in the other peaks then just label as single peak otherwise it is all the peaks
                    pk: int | Literal["all"] = k
                    for l in range(len(settings_class.subfit_orders["peak"])):
                        if not all_data[0][j]["h"][l]:
                            pk = "all"
                    # make the figure title
                    ttlstr = peak_string(settings_class.fit_orders[j], peak=pk)
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
                    # x.set_xticks(x_ticks)
                    cb = plt.colorbar(extend=cb_extend)
                    # cb.set_label(r"Log$_{10}$(Intensity)")
                    cb.set_label(r"Intensity")

                    # Save the figure
                    filename = make_outfile_name(
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

