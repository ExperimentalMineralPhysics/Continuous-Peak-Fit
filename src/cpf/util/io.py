#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 05:44:13 2021
@author: simon
"""

# Enables modern type hinting features on older Python versions
from __future__ import annotations

import os
import re
from copy import deepcopy
from pathlib import Path
from typing import Any, Literal, TypeVar, overload

import numpy as np
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)

import cpf.peak_functions as pf
from cpf.util.logging import get_logger

logger = get_logger("cpf.util.io")

T = TypeVar("T")


# Needed for JSON to save fitted parameters.
# Copied from https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array#24375113
# on 13th November 2018
def numpy_to_json(o):
    """
    Serialize numpy types for json
    Parameters:
        o (object): any python object which fails to be serialized by json
    Example:
        >>> import json
        >>> a = np.array([1, 2, 3])
        >>> json.dumps(a, default=numpy_to_json)
    """
    numpy_types = (
        np.bool_,
        # np.bytes_, -- python `bytes` class is not json serializable
        # np.complex64,  -- python `complex` class is not json serializable
        # np.complex128,  -- python `complex` class is not json serializable
        # np.complex256,  -- special handling below
        # np.datetime64,  -- python `datetime.datetime` class is not json serializable
        np.float16,
        np.float32,
        np.float64,
        # np.float128,  -- special handling below
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        # np.object_  -- should already be evaluated as python native
        np.str_,
        np.timedelta64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.void,
    )

    if isinstance(o, np.ndarray):
        return o.tolist()
    elif isinstance(o, numpy_types):
        return o.item()
    elif isinstance(o, np.float128):
        return o.astype(np.float64).item()
    # elif isinstance(o, np.complex256): -- no python native for np.complex256
    #     return o.astype(np.complex128).item() -- python `complex` class is not json serializable
    else:
        raise TypeError(
            "{} of type {} is not JSON serializable".format(repr(o), type(o))
        )


def file_list(fit_parameters):
    """
    From the Settings make a list of all the data files.
    This function is called by the output writing scripts to make sure the file names are called consistently.

    if datafile_StartNum and datafile_EndNum are in the settings a list of file names is made. If step is present
    it is also used. datafile_StartNum is always present in the file list and if step > 0 then is the fisrt file in the
    list. If step <0 then datafile_StartNum is the last file in the list.
    e.g.
    1. datafile_StartNum = 1
       datafile_EndNum   = 10
       step = 1
       gives file list: 1,2,3,4,5,6,7,8,9,10

    2. datafile_StartNum = 1
       datafile_EndNum   = 10
       step = 4
       gives file list: 1,5,9

    3. datafile_StartNum = 1
       datafile_EndNum   = 10
       step = -4
       gives file list: 9,5,1

    4. datafile_StartNum = 10
       datafile_EndNum   = 1
       step = 4
       gives file list: 10,6,2

    5. datafile_StartNum = 10
       datafile_EndNum   = 1
       step = -4
       gives file list: 2,6,10


    :param fit_parameters:
    :param fit_settings:
    :return:
    """
    # Define step
    if "datafile_Step" in list(fit_parameters):
        step = fit_parameters["datafile_Step"]
    else:
        step = 1

    if "datafile_NumDigit" not in list(fit_parameters):
        fit_parameters["datafile_NumDigit"] = 1

    # Diffraction patterns -- make list of files
    diff_files = []
    if "datafile_Files" not in list(fit_parameters) and "datafile_StartNum" not in list(
        fit_parameters
    ):
        # There is only a single file because nothing else is defined
        n_diff_files = 1
        diff_files.append(
            os.path.abspath(
                fit_parameters.get("datafile_directory", ".")
                + os.sep
                + fit_parameters.get("datafile_Basename", "")
                + fit_parameters.get("datafile_Ending", "")
            )
        )
        # diff_files.append(
        #     os.path.abspath(
        #         getattr(fit_settings, 'datafile_directory', '.')
        #         + os.sep
        #         + getattr(fit_settings, 'datafile_Basename', '')
        #         + getattr(fit_settings, 'datafile_Ending', '')
        #     )
        # )
    elif "datafile_Files" not in list(fit_parameters):
        n_diff_files = int(
            np.floor(
                np.abs(
                    fit_parameters["datafile_EndNum"]
                    - fit_parameters["datafile_StartNum"]
                )
                / np.abs(step)
            )
            + 1
        )
        for j in range(n_diff_files):
            # Make list of diffraction pattern names and no. of pattern
            if fit_parameters["datafile_EndNum"] >= fit_parameters["datafile_StartNum"]:
                n = str(fit_parameters["datafile_StartNum"] + (j * np.abs(step))).zfill(
                    fit_parameters["datafile_NumDigit"]
                )
            else:
                n = str(
                    fit_parameters["datafile_StartNum"] + (j * -np.abs(step))
                ).zfill(fit_parameters["datafile_NumDigit"])
            # Append diffraction pattern name and directory
            diff_files.append(
                os.path.abspath(
                    fit_parameters.get("datafile_directory", ".")
                    + os.sep
                    + fit_parameters.get("datafile_Basename", "")
                    + n
                    + fit_parameters.get("datafile_Ending", "")
                )
            )
        if step < 0:
            diff_files = diff_files[::-1]

    elif "datafile_Files" in list(fit_parameters):
        n_diff_files = int(
            np.round(len(fit_parameters["datafile_Files"]) / np.abs(step))
        )
        for j in range(n_diff_files):
            # Make list of diffraction pattern names and no. of pattern
            if step < 0:
                n = str(fit_parameters["datafile_Files"][j * step - 1]).zfill(
                    fit_parameters["datafile_NumDigit"]
                )
            else:
                n = str(fit_parameters["datafile_Files"][j * step]).zfill(
                    fit_parameters["datafile_NumDigit"]
                )
            # Append diffraction pattern name and directory
            diff_files.append(
                os.path.abspath(
                    fit_parameters.get("datafile_directory", ".")
                    + os.sep
                    + fit_parameters.get("datafile_Basename", "")
                    + n
                    + fit_parameters.get("datafile_Ending", "")
                )
            )
    else:
        n_diff_files = int(len(fit_parameters["datafile_Files"]) / step + 1)
    return diff_files, n_diff_files


def image_list(fit_parameters, files_only=False):
    """
    From the Settings make a list of all the data images to be processed.
    If the images are h5 type files a list of files is made first then the list is expanded for the images in the h5 files.

    #FIXME: This function is called by the output writing scripts to make sure the file names are called consistently.


    Parameters
    ----------
    fit_parameters : TYPE
        DESCRIPTION.
    fit_settings : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # Local import to avoid circular errors
    import cpf.h5_functions as h5_functions

    # make the file list
    diff_files, n_diff_files = file_list(fit_parameters)

    if files_only != True and "*" not in diff_files[0]:
        # iterate for h5 files.
        image_list = []
        if "h5_datakey" in fit_parameters:
            # if new h5 format is present then call it.
            for i in range(n_diff_files):
                h5_list = h5_functions.get_image_keys_new(
                    diff_files[i],
                    h5key_data=fit_parameters["h5_datakey"],
                    h5_iterate=fit_parameters["h5_iterate"],
                )
                for j in range(len(h5_list)):
                    tmp = [diff_files[i]]
                    tmp.extend(h5_list[j])
                    image_list.append(tmp)

        elif "h5_key_list" in fit_parameters:
            # if the input contains the old hdf5 file instircutions make the new
            # dictionary based format and call that
            h5datakey, h5iterations = h5_functions.update_key_structure(fit_parameters)
            for i in range(n_diff_files):
                h5_list = h5_functions.get_image_keys_new(
                    diff_files[i],
                    h5key_data=h5datakey,
                    h5_iterate=h5iterations,
                )
                for j in range(len(h5_list)):
                    tmp = [diff_files[i]]
                    tmp.extend(h5_list[j])
                    image_list.append(tmp)

        else:
            image_list = diff_files

    else:
        image_list = diff_files

    n_images = len(image_list)

    return diff_files, n_diff_files, image_list, n_images


def get_file_indices(
    param_dict: dict | None = None,
    start_num: int | None = None,
    end_num: int | None = None,
    step: int = 1,
    files: list[int] | None = None,
    keys: list[int] | None = None,
):
    """
    Converts parameters into a list of numbers corresponding to the file numbers
    or HDF5 keys of the data to be processed.

    The parameters can come from a Settings class, a dictionary, or the 'start',
    'stop', 'step', 'files', and 'keys' parameters directly.

    If 'param_dict' is populated, the 'start', 'stop' and 'step' values should be
    derived entirely from it. All other parameters will be ignored.

    Other supported combinations of parameters include:
    - 'files'
    - 'files', 'step'
    - 'files', 'start_num', 'end_num'
    - 'files', 'start_num', 'end_num', 'step'
    - 'keys'
    - 'keys', 'step'
    - 'keys', 'start_num', 'end_num'
    - 'keys', 'start_num', 'end_num', 'step'
    - 'start_num', 'end_num'
    - 'start_num', 'end_num', 'step'

    'keys' and 'files' cannot both be provided simultaneously.

    The hierarchy works as follow:
    1. If files (or keys) is given and start_num, end_num, step are not given
    the contents of files/keys is returned.
    2. If files (or keys) is given and start_num, end_num or step are also given
    files/keys are trimmed to values between start_num and end_num using step.
    - if only step and step = -1, the order of files/keys is reversed.
    3. If only start_num, end_num +/- step are given a list of numbers is
    made from these values.
    - start_num is present in the file list if step > 0
    - end_num is present in the list if step < 0

    Examples
    --------
    A.  files = [1,3,4,6,7,9,10]
        gives: [1,3,4,6,7,9,10]

    B.  files = [1,3,4,6,7,9,10]
        start_num = 2
        end_num   = 9
        gives: [3,4,6,7,9]

    C.  files = [1,3,4,6,7,9,10]
        start_num = 2
        end_num   = 9
        step = -1
        gives: [9,7,6,4,3]

    D.  start_num = 1
        end_num   = 10
        step = 1
        gives: 1,2,3,4,5,6,7,8,9,10

    E.  start_num = 1
        end_num   = 10
        step = 4
        gives: 1,5,9

    F.  start_num = 1
        end_num   = 10
        step = -4
        gives: 9,5,1

    G.  start_num = 10
        end_num   = 1
        step = 4
        gives: 10,6,2

    H.  start_num = 10
        end_num   = 1
        step = -4
        gives: 2,6,10


    Parameters
    ----------
    param_dict : dict[str, int | None] | None
        The dictionary containing parameters on the start and end numbers of the files
        or HDF5 keys to be analysed, along with the step interval to use.
        The default value is None.
    start_num : int | None
        Starting value for the file index/number. The default is None.
    end_num : int | None
        End value for the file index/number. The default is None.
    step : int
        step value for the indices. The default is 1.
    files : list[int] | None
        List of file indices. The default is None.
    keys : list[int] | None
        List of key indices for hdf5 file. The default is None.

    Returns
    -------
    indices_list : list
        List of numbers in the file names or the HDF5 keys to process.
    num_indices : int
        Number of entries in indices_list.
    """

    indices_list: list[int] = []
    match (
        param_dict is not None,
        files is not None,
        keys is not None,
    ):
        case (True, True, True) | (True, True, False) | (True, False, True):
            # 'param_dict' and 'files'/'keys' cannot be provided together
            raise ValueError(
                "'param_dict' cannot be provided together with 'files' and/or 'keys'."
            )
        case (False, True, True):
            # Raise error if 'files' and 'keys' are both present,
            raise ValueError("'files' and 'keys' cannot both be provided.")
        case (False, True, False):
            if files is None:
                raise ValueError("'files'' was not provided")
            indices_list = files
        case (False, False, True):
            if keys is None:
                raise ValueError("'keys' was not provided")
            indices_list = keys
        case (True, False, False):
            if param_dict is None:
                raise ValueError("'param_dict' was not provided")
            # Get 'start', 'stop', and 'step' from 'param_dict'
            if isinstance(param_dict, dict):
                # if is a dictionary assume from hdf5 files.
                # dictionary from hdf5 functions
                start_num = param_dict.get("from", None)
                end_num = param_dict.get("to", None)
                step = param_dict.get("step", 1)
            else:
                start_num = param_dict.get("datafile_StartNum", None)
                end_num = param_dict.get("datafile_EndNum", None)
                step = param_dict.get("datafile_Step", 1)
            # Error if 'start_num' and 'end_num' could not be determined
            if not (start_num is not None and end_num is not None):
                raise ValueError("Could not construct indices list from 'param_dict'.")
            start_num, end_num, step = map(int, (start_num, end_num, step))

    # Process differently depending on the absence/presence of 'start_num',
    # 'end_num', 'indices_list'
    match (start_num is not None, end_num is not None, len(indices_list) > 0):
        case (True, False, False) | (False, True, False) | (False, False, False):
            # Raise error if 'start_num' and/or 'end_num' are missing when 'indices_list'
            # has been provided
            raise ValueError(
                "Both 'start_num' and 'end_num' must be set if neither 'files' or "
                "'keys' were provided."
            )
        case (False, _, True):
            # Reverse the list if step is negative
            if step < 0:
                indices_list = list(reversed(indices_list))
            step = abs(step)
            # Use 'step' to return every nth item in the list
            indices_subset = [
                index for i, index in enumerate(indices_list) if i % step == 0
            ]
            return indices_subset, len(indices_subset)
        case (True, _, True) | (True, True, False):
            # Assure type checker that 'start_num' is an int by this point
            if start_num is None:
                raise ValueError("'start_num' was not set")

            # Do not allow negative values for 'start_num' and 'end_num'
            if start_num < 0 or (end_num is not None and end_num < 0):
                raise ValueError("'start_num' and 'end_num' cannot be negative values")

            # If 'step' is negative, reverse the indices list after creation
            reverse_list = step < 0

            if indices_list:
                # Create a default 'end_num' if none was provided
                end_num = end_num if end_num is not None else indices_list[-1]

                # Find the indices nearest to the specified 'start_num' and 'end_num'
                idx_start: int | None = None
                idx_end: int | None = None
                if start_num > end_num:
                    indices_list = list(reversed(indices_list))
                    for i, num in enumerate(indices_list):
                        if num <= start_num and idx_start is None:
                            idx_start = i
                        if num == end_num and idx_end is None:
                            idx_end = i
                        elif num < end_num and idx_end is None:
                            idx_end = i - 1
                else:
                    for i, num in enumerate(indices_list):
                        if num >= start_num and idx_start is None:
                            idx_start = i
                        if num == end_num and idx_end is None:
                            idx_end = i
                        elif num > end_num and idx_end is None:
                            idx_end = i - 1
                if not (idx_start is not None and idx_end is not None):
                    raise ValueError(
                        "Could not determine start and end indices using the start "
                        "and end numbers provided"
                    )
                positions = range(idx_start, idx_end + 1, abs(step))
                indices_list = [indices_list[i] for i in positions]
                if reverse_list:
                    indices_list = list(reversed(indices_list))
                return indices_list, len(indices_list)
            else:
                # 'end_num' will be set at this point
                if end_num is None:
                    raise ValueError("'end_num' was not set")

                # Generate number list that includes 'end_num' itself
                if start_num < end_num:
                    end_num += 1
                    step = abs(step)
                else:
                    end_num -= 1
                    step = -abs(step)
                indices_list = list(range(start_num, end_num, step))
                if reverse_list:
                    indices_list = list(reversed(indices_list))
                return indices_list, len(indices_list)
    # This should never trigger, but is here to ensure a list is always returned
    logger.warning(
        "Unexpected combination of parameters detected, returning empty list"
    )
    return [], 0


@overload
def has_value(
    obj: dict,
    val: str | int | float | None = None,
    path: str = "",
    is_present: bool = False,
) -> bool: ...


@overload
def has_value(
    obj: list,
    val: str | int | float | None = None,
    path: str = "",
    is_present: bool = False,
) -> bool: ...


@overload
def has_value(
    obj: pd.DataFrame,
    val: str | int | float | None = None,
    path: str = "",
    is_present: bool = False,
) -> bool: ...


def has_value(
    obj: Any,
    val: Any = None,
    path: str = "",
    is_present: bool = False,
):
    """
    This function recursively searches through lists, dictionaries, and
    Pandas DataFrames for the specified value ("None" by default) and
    returns True if any are detected.

    Parameters
    ----------
    obj : dict, list, pd.DataFrame
        Nested dictionary/list of parameters or Pandas DataFrame to inspect.
    val : str, int, float, None
        Value or string to find in the dictionary. The default is None.
    path : str
        Path taken through the dictionary/list. The default is "".
    is_present : bool
        Boolean for if 'val' are in dictionary. Used for iterating through nested structures.
        The default is False.

    Returns
    -------
    is_present : bool
        True - if any instances of 'val' have been found in the dicionary
        False - if 'val' is not in dictionary.

    """
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021
    if isinstance(obj, dict):
        for key, value in obj.items():
            if is_present := has_value(value, val, path + f"['{key}']", is_present):
                break  # Stop once one is found
    elif isinstance(obj, list):
        for key, value in enumerate(obj):
            if is_present := has_value(value, val, path + f"[{key}]", is_present):
                break  # Stop once one is found
    elif isinstance(obj, pd.DataFrame):
        # Look in DataFrame for values
        is_present = (
            bool((obj.isna() if pd.isna(val) else obj.eq(val)).any().any())
            or is_present
        )
    elif obj == val:
        is_present = True
        logger.moreinfo(" ".join(map(str, [(f"Value {val} found at {path}")])))
        # could be verbose if verbose logger.
    return is_present


@overload
def replace_value(
    obj: dict,
    old: str | int | float | None = None,
    new: str | int | float | None = 0,
    path: str = "",
) -> dict: ...


@overload
def replace_value(
    obj: list,
    old: str | int | float | None = None,
    new: str | int | float | None = 0,
    path: str = "",
) -> list: ...


@overload
def replace_value(
    obj: pd.DataFrame,
    old: str | int | float | None = None,
    new: str | int | float | None = 0,
    path: str = "",
) -> pd.DataFrame: ...


def replace_value(
    obj: T,
    old: Any = None,
    new: Any = 0,
    path: str = "",
) -> T:
    """
    Recursively replaces the specified old value in a dictionary, list, or Pandas
    DataFrame with the desired new value, and returns the updated object.

    Parameters
    ----------
    obj : dict, list, pd.DataFrame
        Nested dictionary/list or Pandas DataFrame to inspect.
    old : str, int, float, None
        Value or string to find in the dictionary. The default is None.
    new : str, int, float, None
        Value or string to use as replacement in the dictionary.
        The default is 0.
    path : str
        Path taken through the dictionary/list. The default is "".

    Returns
    -------
    obj :  dict, list, pd.DataFrame
        Nested dictionary or list of parameters.

    """
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021

    if isinstance(obj, (dict, list)):
        for key, value in obj.items() if isinstance(obj, dict) else enumerate(obj):
            obj[key] = replace_value(deepcopy(value), old, new, path + f"['{key}']")
    elif isinstance(obj, pd.DataFrame):
        # replace contents of panda data frame
        old = np.nan if old is None else old
        obj = obj.replace(old, new).infer_objects(copy=False)
    elif obj == old:  # and old is not None:
        obj = new
        logger.moreinfo(" ".join(map(str, [(f"Value {old} found at {path}")])))
    return obj


def has_huge_errors(obj: dict, min_ratio: int | float = 3, is_huge: bool = False):
    """
    This function accepts a nested dictionary and list as argument and iterates over
    all values of nested dictionaries and lists.

    Huge errors are flagged if:
        1. value_err/value >= min_ratio
        2. abs(value)-value_err >= 0  (i.e. not within error of 0)

    Parameters
    ----------
    obj : dict
        Nested dictionary of parameters to inspect.
    min_ratio : float
        Minimum ratio for how big large errors are before being flagged. The default is 3.
    is_huge : bool
        Boolean for if large errors are found. Used for iterating through nested structures.
        The default is False.

    Returns
    -------
    is_huge : bool
        True - if large errors have been found in the dicionary
        False - if no large errors are present.

    """
    for k in range(len(obj["background"])):
        for j in range(len(obj["background"][k])):
            if (
                obj["background"][k][j] != 0
                and obj["background"][k][j] != None
                and obj["background_err"][k][j] != 0
                and obj["background_err"][k][j] != None
                and obj["background_err"][k][j] / obj["background"][k][j] >= min_ratio
            ):
                is_huge = True
                err_rat = obj["background_err"][k][j] / obj["background"][k][j]
                logger.moreinfo(
                    f"Huge errors found in background {k}, {j}: "
                    f"value= {obj['background'][k][j]: 3.2e}; "
                    f"error={obj['background_err'][k][j]: 3.2e}; "
                    f"fractional error = {err_rat: 5.1f}"
                )
    comp_list, comp_names = pf.peak_components(include_profile=True)
    for k in range(len(obj["peak"])):
        for cp in range(len(comp_list)):
            comp = comp_names[cp]
            for j in range(len(obj["peak"][k][comp])):
                # check if
                # - there is a value and an error
                # - the error is less than "min_ratio" * error
                # - the value is not within error of 0
                #       [this is a sanity check to the preceeding check -- the ratio of error/value tends to inifinty as the
                #         value becomes very small.]
                #       [without this there is lots of discarding the previous fit when the profile values are close to 0]
                if (
                    obj["peak"][k][comp][j] != 0
                    and obj["peak"][k][comp][j] != None
                    and obj["peak"][k][comp + "_err"][j] != 0
                    and obj["peak"][k][comp + "_err"][j] != None
                    and obj["peak"][k][comp + "_err"][j] / obj["peak"][k][comp][j]
                    >= min_ratio
                    and np.abs(obj["peak"][k][comp][j])
                    - obj["peak"][k][comp + "_err"][j]
                    >= 0
                ):
                    is_huge = True
                    err_rat = obj["peak"][k][comp + "_err"][j] / obj["peak"][k][comp][j]
                    logger.moreinfo(
                        f"Huge error found in peak {k}, {comp} {j}: "
                        f"value= {obj['peak'][k][comp][j]: 3.2e}; "
                        f"error={obj['peak'][k][comp+'_err'][j]: 3.2e}; "
                        f"fractional error = {err_rat: 5.1f}"
                    )
    return is_huge


def peak_string(
    orders: dict[str, Any],
    peak: int | list[int] | Literal["all"] = "all",
    fname=False,
):
    """
    Parameters
    ----------
    orders: dict[str, list[dict[str, int | str]]]
        A nested dictionary containing information about the peaks to be fitted and
        how to go about fitting them.
        The "peak" key contains a list of dictionaries of the individual peaks to be
        fitted. The dictionaries in turn should contain the 'hkl' key, which holds
        the Miller indices in the form of a string or an integer array, as well as
        the 'phase' key, which contains information about which material this peak
        belongs to.
    peak: int | list[int] | str | Literal['all']
        Indices of the peak descriptors to use to generate the peak string with.
        Takes an integer, a list of integers, or 'all'. The default is 'all'.
    fname: bool
        Toggle whether to construct a file name-compatible or human-readable string.
        If true, peak indices will be placed in parentheses (e.g. Peak (110)), whereas
        they will be written as 'Peak-110_' if it's set to False).
        The default is False.

    Returns
    -------
    p_str: str
        A string listing the peaks processed by the function. They take the format
        "Peak (110) & Peak (220) & ..." if 'fname' is set to False, and are returned
        as "Peak-110_Peak-220_..." if 'fname' is True.
    """
    # Construct list of indices to parse
    if peak == "all":
        peaks = list(range(len(orders["peak"])))
    # If a list of ints is provided
    elif isinstance(peak, list) and all(isinstance(x, int) for x in peak):
        peaks = peak
    # If an int was provided
    elif isinstance(peak, str) or isinstance(peak, int) or np.issubdtype(peak, np.integer):
        peaks = [peak]
    # Raise a TypeError otherwise
    else:
        raise TypeError(
            f"'peak' received an unsupported value: {peak} "
            "It must be 'all', an integer, or a list of integers. "
        )

    # Construct peak string
    p_str = ""
    for n, x in enumerate(peaks):
        # Determine name of peak from 'phase' key; fall back to using 'Peak'
        if "phase" in orders["peak"][x]:
            p_str = p_str + str(peak_phase(orders, peak=x)[0])
        else:
            p_str = p_str + "Peak"
        # Encase indices in brackets if it's human-readable
        if fname is False:
            p_str = p_str + " ("
        else:
            p_str = p_str + "-"
        # Use Miller indices if 'hkl' key is present; fall back to auto-incrementing
        if "hkl" in orders["peak"][x]:
            p_str = p_str + str(peak_hkl(orders, peak=x, as_string=True)[0])
        else:
            p_str = p_str + str(x + 1)
        # Close the bracket if it's human-readable
        if fname is False:
            p_str = p_str + ")"
        # Separate peaks using either '&' if human-readable or '_'
        if len(peaks) > 1 and n < len(peaks) - 1:
            if fname is False:
                p_str = p_str + " & "
            else:
                p_str = p_str + "_"
    return p_str


def peak_hkl(
    orders: dict[str, Any],
    peak: int | list[int] | Literal["all"] = "all",
    as_string: bool = True,
) -> list[str | list[int]]:
    """
    Takes the values stored under the 'hkl' keys in the desired peak descriptors (which
    are stored as a list of dictionaries under the 'peak' key in 'fit_orders') and
    returns them as a list of Miller indices. If 'as_string' is True, the indices are
    returned as filename-friendly strings. Otherwise, they are integer arrays.

    Possible data types stored in 'hkl':
    * string -- hkls are stored as a string. Generally used for indices with a leading
    0 or with negative indices. E.g. "-110", "001"
    * int --  hkls are all positive, with no leading 0. E.g. 100, 200
    * list -- hkls are stored as a list of ints. E.g. [-1, 1, 0]

    Possible outputs:
    * string -- Contents of the returned list are strings. E.g. ["-110", "220"]
    * list -- Contents of the list are integer arrays. E.g. [[1, 1, 0], ...]

    Parameters
    ----------
    orders: dict[str, list[dict[str, int | str]]]
        A nested dictionary containing information about the peaks to be fitted and
        how to go about fitting them.
        The "peak" key contains a list of dictionaries of the individual peaks to be
        fitted. The dictionaries in turn should contain the 'hkl' key, which holds
        the Miller indices in the form of a string or an integer array.
    peak: int | list[int] | Literal['all']
        Indices of the peak descriptors to extract the hkl peaks from. Takes an
        integer, a list of integers, or 'all'. The default is 'all'.
    as_string: bool
        Toggle whether to return the Miller indices as a string or a list of integers.

    Returns
    -------
    out: list[str | list[int]]
       The list of peaks selected, stored as either a list of strings or a list of
       integer arrays.
    """
    # Construct list of peaks to parse
    if peak == "all":
        peaks = list(range(len(orders["peak"])))
    # Accept Python and NumPy ints and convert to Python int
    elif isinstance(peak, (int, np.integer)):
        peaks = [int(peak)]
    elif isinstance(peak, list) and all(isinstance(p, (int, np.integer)) for p in peak):
        peaks = [int(p) for p in peak]
    else:
        raise TypeError(
            f"'peak' received an unsupported value: {peak} "
            "It must be 'all', an integer, or a list of integers. "
        )

    out: list[str | list[int]] = []
    for x in peaks:
        # Handle 0 and if no 'hkl' key is found
        if "hkl" in orders["peak"][x]:
            hkl = orders["peak"][x]["hkl"]
            if hkl == "0" or hkl == 0:
                hkl = "000"
        else:
            hkl = "000"
        # Store peaks as strings
        if not isinstance(hkl, list) and as_string == True:
            # String in, string out
            out.append(str(hkl))
        elif isinstance(hkl, list) and as_string == True:
            # List in, string out
            hkl_string = ""
            for y in range(len(hkl)):
                hkl_string = hkl_string + str(hkl[y])
            out.append(hkl_string)
        # Store peaks as lists of integers
        elif isinstance(hkl, list) and as_string == False:
            # List in, list out
            out.append(hkl)
        elif not isinstance(hkl, list) and as_string == False:
            # String in, list out
            pos = 0
            hkl = str(hkl)
            if hkl[0] == "-":
                h = hkl[pos : pos + 2]
                pos = pos + 2
            else:
                h = hkl[pos : pos + 1]
                pos = pos + 1
            if hkl[pos] == "-":
                k = hkl[pos : pos + 2]
                pos = pos + 2
            else:
                k = hkl[pos : pos + 1]
                pos = pos + 1
            if hkl[pos] == "-":
                l = hkl[pos : pos + 2]
                pos = pos + 2
            else:
                l = hkl[pos : pos + 1]
                pos = pos + 1
            hkl_list = [int(h), int(k), int(l)]
            if len(hkl) > pos:
                if hkl[pos] == "-":
                    m = hkl[pos : pos + 2]
                    pos = pos + 2
                else:
                    m = hkl[pos : pos + 1]
                    pos = pos + 1
                hkl_list.append(int(m))
            out.append(hkl_list)
    return out


def peak_phase(
    orders: dict[str, Any],
    peak: int | list[int] | Literal["all"] = "all",
):
    """
    Extracts and returns the values stored under the 'phase' keys in the desired peak
    descriptors (which are stored as a list of dictionaries under the 'peak' key in
    'fit_orders') and returns them as a list of strings. If the 'phase' key is not
    present, returns 'Unknown' as the phase name instead.

    Parameters
    ----------
    orders: dict[str, list[dict[str, int | str]]]
        A nested dictionary containing information about the peaks to be fitted and
        how to go about fitting them.
        The "peak" key contains a list of dictionaries of the individual peaks to be
        fitted. The dictionaries in turn should contain the 'phase' key, which holds
        the name of the phase/material this peak belongs to.
    peak: int | list[int] | Literal['all']
        The indices of the peak descriptors to extract the phase name from. Takes an
        integer, a list of integers, or 'all'. The default is 'all'.

    Returns
    -------
    out: list[str]
       The list of phase names for the selected peaks.
    """
    # Convert input into a list of integers
    if peak == "all":
        peaks = list(range(len(orders["peak"])))
    # Accept NumPy and Python ints and convert to Python ints
    elif isinstance(peak, (int, np.integer)):
        peaks = [int(peak)]
    elif isinstance(peak, list) and all(isinstance(p, (int, np.integer)) for p in peak):
        peaks = [int(p) for p in peak]
    else:
        raise TypeError(
            f"'peak' received an unsupported value: {peak} "
            "It must be 'all', an integer, or a list of integers. "
        )
    # Load names of specified phases, with a fallback to 'Unknown'
    out: list[str] = []
    for x in peaks:
        if "phase" in orders["peak"][x]:
            out.append(orders["peak"][x]["phase"])
        else:
            out.append("Unknown")
    return out


def title_file_names(settings_for_fit=None, num=0, image_name=None, string=True):
    if string == True:
        joint = ";  "
    else:
        joint = "_"

    if image_name == None:
        image_name = settings_for_fit.image_list[num]

    if isinstance(image_name, list):
        t_f_str = os.path.split(image_name[0])[1]
        t_f_str, _ = os.path.splitext(t_f_str)
        t_f_str += joint
        t_f_str += str(image_name[-1])
    else:
        _, t_f_str = os.path.split(image_name)
        t_f_str, _ = os.path.splitext(t_f_str)

    return t_f_str


def make_outfile_name(
    base_filename: list | Path | str,
    directory=None,
    additional_text=None,
    extension=None,
    orders=None,
    peak="all",
    overwrite=True,
):
    """
    Make file names for output files.
    Needs to account for different file name formats
    :param base_filename:
    :param directory:
    :param additional_text:
    :param extension:
    :param orders:
    :param overwrite:
    :return:
    """

    # If the file type is h5 the name arrives as a list of bits
    if isinstance(base_filename, list):
        # Get the file name and file extension
        logger.debug(f"Received h5 file type with base filename {base_filename}")
        filename = title_file_names(image_name=base_filename, string=False)
        ending = Path(base_filename[0]).suffix[1:]  # Get file extension and remove dot
    elif isinstance(base_filename, (str, Path)):
        # Convert to a Path object and get the file name and file extension
        logger.debug(f"Received file with base filename {base_filename}")
        base_filename = (
            Path(base_filename) if isinstance(base_filename, str) else base_filename
        )
        ending = base_filename.suffix[1:]  # Remove leading dot from file extension
        filename = base_filename.stem  # Take last level of the file path
    else:
        err_str = "Unexpected type passed in for base_filename parameter"
        logger.error(err_str)
        raise ValueError(err_str)

    # Strip directory if it is in the name and there is a new directory
    if directory or directory == "":
        filename = Path(filename).stem

    # Remove trailing "_" from base filename if present
    filename = filename[0:-1] if filename.endswith("_") else filename

    # assume here that the only illicit symbol in the file name is *, from
    # ESRF compound detector. Replace with text. If any other illicit characters
    # are present it is going to look very strange.
    filename = licit_filename(filename, replacement="multi")

    # check the iteration is not in the file ending.
    # another way to do this would be to check in the file ending in the input file is empty.
    try:
        int(ending)
        in_ending = True  # the file ending is convertable to a number -- therefore assume the extension
        # contains the iterating parameter.
    except:
        # FIX ME: DMF this exception needs handling properly!
        in_ending = False
    if in_ending:
        filename = filename + "_" + ending

    if orders:  # add phase and hkl to file name
        filename = filename + "__" + peak_string(orders, fname=True, peak=peak)
    if additional_text:  # add additional text -- i.e. notes from orders
        filename = filename + "__" + additional_text
    if orders and "note" in orders:  # add additional text from note in orders.
        filename = (
            filename + "__" + "".join(i for i in orders["note"] if i not in r"\/:;*?<>|")
        )
    filename = filename.strip("_")
    if directory:
        filename = os.path.join(directory, filename)

    if extension and extension[0] == ".":
        extension = extension[1:]

    if overwrite is False:
        i = 0
        if os.path.exists("{}.{}".format(filename, extension)):
            i += 1
        while os.path.exists("{}_{:d}.{}".format(filename, i, extension)):
            i += 1
        if i != 0:
            filename = "{}_{:d}".format(filename, i)
    if extension:
        filename = "{}.{}".format(filename, extension)

    return filename


def lmfit_fix_int_data_type(fname):
    """
    fixes problem with lmfit save/load model.
    lmfit load model cannot read int32 data with nulls in it.
    if replace 'int32' with 'float32' it will read.
    """

    obj_read = open(fname, "r")
    txt_content = obj_read.read()
    obj_read.close()

    txt_content = txt_content.replace("uint", "float")
    txt_content = txt_content.replace("int", "float")
    logger.effusive(" ".join(map(str, [("    Rewriting", fname)])))

    obj_read = open(fname, "w")
    obj_read.write(txt_content)
    obj_read.close()


def number_to_string(number, replace=".", withthis="pt"):
    """
    Turns a number into a string and then replaces the decimal place with a "pt".
    """
    number = str(number)
    number = re.sub(r"\.+", "pt", number)
    return number


def licit_filename(fname, replacement="==", exclude_dir=True):
    """
    This makes sure that a file name generated from a string is licit.
    It replaces the illegal characters [<>:|?*] with a replacement character
    If exclude_dir==False it also replaces / and \
    It also replaces all '.' with 'pt' -- assuming any occurance is a number.

    after https://gist.github.com/AaronLaw/a936bebfbbd691fc954252444767e6de -- Find NTFS illegal characters in black list and rename filename.
    """
    blacklist = r"[<>:|?*]"
    # the file name might include a directory link...
    # if exclude_dir == False:
    #     blacklist += r"/\\"

    fname = re.sub(blacklist, replacement, fname)

    # the file name might include a directory link...
    if exclude_dir == False:
        fname = fname.replace("/", replacement)
        fname = fname.replace("\\", replacement)

    fname = number_to_string(fname)

    return fname


def figure_suptitle_space(figure, topmargin=1):
    """increase figure size to make topmargin (in inches) space for
    titles, without changing the axes sizes.
    after: https://stackoverflow.com/questions/55767312/how-to-position-suptitle#55768955

            Acutally now does this by compresssing the axes away from the top of the figure.
    """

    axes = figure.axes
    pos = []
    for i in range(len(axes)):
        pos.append(axes[i].get_position().bounds)
    w, h = figure.get_size_inches()
    figh = h - topmargin  # - (1-s.y1)*h
    for i in range(len(axes)):
        al = pos[i][0]
        ab = pos[i][1] / h * figh
        aw = pos[i][2]
        ah = pos[i][3] / h * figh
        axes[i].set_position((al, ab, aw, ah))
