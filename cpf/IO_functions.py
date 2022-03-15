#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 05:44:13 2021
@author: simon
"""

import numpy as np
import os


# Needed for JSON to save fitted parameters.
# Copied from https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array#24375113
# on 13th November 2018
def json_numpy_serializer(o):
    """
    Serialize numpy types for json
    Parameters:
        o (object): any python object which fails to be serialized by json
    Example:
        >>> import json
        >>> a = np.array([1, 2, 3])
        >>> json.dumps(a, default=json_numpy_serializer)
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


def file_list(fit_parameters, fit_settings):
    """
    From the Settings make a list of all the data files.
    This function is called by the output writing scripts to make sure the file names are called consistently.
    :param fit_parameters:
    :param fit_settings:
    :return:
    """
    # Define step
    if ('datafile_Step' in fit_parameters):
        step = fit_settings.datafile_Step
    elif ('datafile_Files' in fit_parameters):
        step = 1
    else: #if just datafile_StartNum and datafile_EndNum
        if fit_settings.datafile_EndNum >= fit_settings.datafile_StartNum:
            step = 1
        else:
            step = -1
    
    if "datafile_NumDigit" not in fit_parameters:
        fit_settings.datafile_NumDigit = 1

    # Diffraction patterns -- make list of files
    diff_files = []
    if "datafile_Files" not in fit_parameters:
        n_diff_files = int(
            np.floor(
                (fit_settings.datafile_EndNum - fit_settings.datafile_StartNum) / step + 1
            )
        )
        for j in range(n_diff_files):
            # Make list of diffraction pattern names and no. of pattern
            n = str(fit_settings.datafile_StartNum + j * step).zfill(
                fit_settings.datafile_NumDigit
            )
            # Append diffraction pattern name and directory
            diff_files.append(
                os.path.abspath(
                    fit_settings.datafile_directory
                    + os.sep
                    + fit_settings.datafile_Basename
                    + n
                    + fit_settings.datafile_Ending
                )
            )
    elif "datafile_Files" in fit_parameters:
        n_diff_files = int(len(fit_settings.datafile_Files) / step)
        for j in range(n_diff_files):
            # Make list of diffraction pattern names and no. of pattern
            n = str(fit_settings.datafile_Files[j * step]).zfill(
                fit_settings.datafile_NumDigit
            )
            # Append diffraction pattern name and directory
            diff_files.append(
                os.path.abspath(
                    fit_settings.datafile_directory
                    + os.sep
                    + fit_settings.datafile_Basename
                    + n
                    + fit_settings.datafile_Ending
                )
            )
    else:
        n_diff_files = int(len(fit_settings.datafile_Files) / step + 1)
    return diff_files, n_diff_files


def any_terms_null(obj_to_inspect, val_to_find=None, index_path="", clean=None):
    """
    This function accepts a nested dictionary and list as argument
    and iterates over all values of nested dictionaries and lists.
    If any of the values are "Null" it returns 0
    :param obj_to_inspect:
    :param val_to_find:
    :param index_path:
    :param clean:
    :return:
    """
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021
    if clean is None:
        clean = 1

    if isinstance(obj_to_inspect, dict):
        for key, value in obj_to_inspect.items():
            clean = any_terms_null(
                value, val_to_find, index_path + f"['{key}']", clean=clean
            )

    if isinstance(obj_to_inspect, list):
        for key, value in enumerate(obj_to_inspect):
            clean = any_terms_null(
                value, val_to_find, index_path + f"[{key}]", clean=clean
            )

    if obj_to_inspect == val_to_find:
        clean = 0
        print(f"Value {val_to_find} found at {index_path}")

    return clean


def replace_null_terms(obj_to_inspect, val_to_find=None, index_path="", clean=None):
    """
    This function accepts a nested dictionary and list as argument
    and iterates over all values of nested dictionaries and lists.
    If any of the values are "Null" it returns 0
    :param obj_to_inspect:
    :param val_to_find:
    :param index_path:
    :param clean:
    :return:
    """
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021

    if isinstance(obj_to_inspect, dict):
        for key, value in obj_to_inspect.items():
            obj_to_inspect[key] = replace_null_terms(
                value, val_to_find, index_path + f"['{key}']", clean=clean
            )

    if isinstance(obj_to_inspect, list):
        for key, value in enumerate(obj_to_inspect):
            obj_to_inspect[key] = replace_null_terms(
                value, val_to_find, index_path + f"[{key}]", clean=clean
            )

    if obj_to_inspect == val_to_find:
        obj_to_inspect = np.nan
        print(f"Value {val_to_find} found at {index_path}")

    return obj_to_inspect


def peak_string(orders, fname=False, peak="all"):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param fname: switch indicting if the string is to be a file name or not.
    :param peak: switch to add restricted peaks to the output string
    :return: string listing peak names
    """
    if peak == "all":
        peek = list(range(len(orders["peak"])))
    elif not isinstance(peak, list):
        peek = [int(x) for x in str(peak)]
    else:
        peek = peak
    p_str = ""
    for x in peek:
        if "phase" in orders["peak"][x]:
            p_str = p_str + orders["peak"][x]["phase"]
        else:
            p_str = p_str + "Peak"
        if fname is False:
            p_str = p_str + " ("
        else:
            p_str = p_str + "-"
        if "hkl" in orders["peak"][x]:
            p_str = p_str + str(orders["peak"][x]["hkl"])
        else:
            p_str = p_str + str(x + 1)
        if fname is False:
            p_str = p_str + ")"

        if peak == "all" and x < len(orders["peak"]) - 1 and len(orders["peak"]) > 1:
            if fname is False:
                p_str = p_str + " & "
            else:
                p_str = p_str + "_"
    return p_str


def make_outfile_name(
        base_filename,
        directory=None,
        additional_text=None,
        extension=None,
        orders=None,
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

    # strip directory if it is in the name and there is a new directory
    if directory or directory == "":
        base_filename = os.path.basename(base_filename)

    if base_filename:  # catch in case the string does not exist.
        filename, ending = os.path.splitext(base_filename)
        ending = ending[1:]
    else:
        filename = "Fit2Peak"
        ending = ""

    if filename[-1:] == "_":  # if the base file name ends in an '_' remove it.
        filename = filename[0:-1]

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
        filename = filename + "__" + peak_string(orders, fname=True)
    if additional_text:  # add additional text -- i.e. notes from orders
        filename = filename + "__" + additional_text
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
    print("    Rewriting", fname)

    obj_read = open(fname, "w")
    obj_read.write(txt_content)
    obj_read.close()
