#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 05:44:13 2021
@author: simon
"""

import numpy as np
import os
import re
import cpf.h5_functions as h5_functions
import cpf.peak_functions as pf
from copy import deepcopy

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


def image_list(fit_parameters, fit_settings):
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

    # make the file list
    diff_files, n_diff_files = file_list(fit_parameters, fit_settings)

    # iterate for h5 files.
    image_list = []
    if "h5_key_list" in fit_parameters:

        # FIX ME: all this code should be moved to settings and validation.
        h5_key_list = fit_settings.h5_key_list

        if "h5_key_names" in fit_parameters:
            h5_key_names = fit_settings.h5_key_names
        else:
            h5_key_names = []

        if "h5_key_start" in fit_parameters:
            h5_key_start = fit_settings.h5_key_start
        else:
            h5_key_start = 0
        # h5_key_start = fit_settings.h5_key_start

        if "h5_key_end" in fit_parameters:
            h5_key_end = fit_settings.h5_key_end
        else:
            h5_key_end = -1
        #  h5_key_end   = fit_settings.h5_key_end
        if isinstance(h5_key_end, int):
            h5_key_end = [h5_key_end]

        if "h5_key_step" in fit_parameters:
            h5_key_step = fit_settings.h5_key_step
        else:
            h5_key_step = 1
        # h5_key_step  = fit_settings.h5_key_step

        if "h5_data" in fit_parameters:
            h5_data = fit_settings.h5_data
        else:
            h5_data = "iterate"
        # h5_data      = fit_settings.h5_data

        for i in range(n_diff_files):
            h5_list = h5_functions.get_image_keys(
                diff_files[i],
                h5_key_list,
                h5_key_names,
                key_start=h5_key_start,
                key_end=deepcopy(h5_key_end),
                key_step=h5_key_step,
                bottom_level=h5_data,
            )
            # N.B. deepcopying of h5_key_end is needed otherwise it is reset for subsequent h5 files.

            for j in range(len(h5_list)):
                image_list.append([diff_files[i], h5_list[j]])

    else:
        image_list = diff_files

    n_images = len(image_list)

    return diff_files, n_diff_files, image_list, n_images


def file_list(fit_parameters, fit_settings):
    """
    From the Settings make a list of all the data files.
    This function is called by the output writing scripts to make sure the file names are called consistently.
    :param fit_parameters:
    :param fit_settings:
    :return:
    """
    # Define step
    if "datafile_Step" in fit_parameters:
        step = fit_settings.datafile_Step
    elif "datafile_Files" in fit_parameters:
        step = 1
    elif "datafile_StartNum" in fit_parameters and "datafile_EndNum" in fit_parameters:
        if fit_settings.datafile_EndNum >= fit_settings.datafile_StartNum:
            step = 1
        else:
            step = -1
    else: # it must be a single file with no numberical compoent.
        step = None

    if "datafile_NumDigit" not in fit_parameters:
        fit_settings.datafile_NumDigit = 1

    # Diffraction patterns -- make list of files
    diff_files = []
    if step == None:
        # There is only a single file because nothing else is defined
        n_diff_files = 1
        diff_files.append(
            os.path.abspath(
                fit_settings.datafile_directory
                + os.sep
                + fit_settings.datafile_Basename
                + fit_settings.datafile_Ending
            )
        )
    elif "datafile_Files" not in fit_parameters:
        n_diff_files = int(
            np.floor(
                (fit_settings.datafile_EndNum - fit_settings.datafile_StartNum) / step
                + 1
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


def replace_null_terms(obj_to_inspect, val_to_find=None, index_path="", clean=None, replace_with=0):
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
        obj_to_inspect = replace_with
        print(f"Value {val_to_find} found at {index_path}")

    return obj_to_inspect


def any_errors_huge(obj_to_inspect, large_errors=3, clean=None):
    """
    This function accepts a nested dictionary and list as argument
    and iterates over all values of nested dictionaries and lists.
    If any of the error values are more than scale times the fitted value it 
    flags the errors as huge 
    :param obj_to_inspect:
    :param large_errors:
    :param clean:
    :return:
    """
    if clean is None:
        clean = 1
    for k in range(len(obj_to_inspect["background"])):
        for j in range(len(obj_to_inspect["background"][k])):
            if (obj_to_inspect["background"][k][j] != 0 and 
                obj_to_inspect["background_err"][k][j] != 0 and 
                obj_to_inspect["background_err"][k][j] != None and 
                obj_to_inspect["background_err"][k][j]/obj_to_inspect["background"][k][j] >= large_errors):
                clean = 0
                err_rat = obj_to_inspect["background_err"][k][j] / obj_to_inspect["background"][k][j]
                print(f"Huge errors found in background {k}, {j}: value= {obj_to_inspect['background'][k][j]: 3.2e}; error={obj_to_inspect['background_err'][k][j]: 3.2e}; fractional error = {err_rat: 5.1f}")
        
    comp_list, comp_names = pf.peak_components(include_profile=True)
    for k in range(len(obj_to_inspect["peak"])):
        for cp in range(len(comp_list)):
            comp = comp_names[cp]
            for j in range(len(obj_to_inspect["peak"][k][comp])):
                if (obj_to_inspect["peak"][k][comp+"_err"][j] != 0 and 
                    obj_to_inspect["peak"][k][comp+"_err"][j] != None and 
                    obj_to_inspect["peak"][k][comp+"_err"][j]/obj_to_inspect["peak"][k][comp][j] >= large_errors):
                    clean = 0
                    err_rat = obj_to_inspect["peak"][k][comp+"_err"][j]/obj_to_inspect["peak"][k][comp][j]
                    print(f"Huge error found in peak {k}, {comp} {j}: value= {obj_to_inspect['peak'][k][comp][j]: 3.2e}; error={obj_to_inspect['peak'][k][comp+'_err'][j]: 3.2e}; fractional error = {err_rat: 5.1f}")

    return clean

        
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
            # p_str = p_str + orders["peak"][x]["phase"]
            p_str = p_str + peak_phase(orders, peak=x)[0]
        else:
            p_str = p_str + "Peak"
        if fname is False:
            p_str = p_str + " ("
        else:
            p_str = p_str + "-"
        if "hkl" in orders["peak"][x]:
            # p_str = p_str + str(orders["peak"][x]["hkl"])
            p_str = p_str + peak_hkl(orders, peak=x, string=True)[0]
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


def peak_hkl(orders, peak="all", string=True):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param string: switch indicting if the output is a string or numeric list of hkl.
    :param peak: switch to add restricted peaks to the output string
    :return: string or list of peak hkl.
    """
    # possible hkl input formats:
    #   string -- hkls as a string, generally used for hkls with leading 0 or negative hkl indicy
    #   int -- hkls all positive and no leading 0
    #   list -- e.g. [h,k,l]. New format.
    # desired outputs:
    #   string -- hkls as a string using String==True
    #   list -- e.g. [h,k,l]. New format. using String==False

    if peak == "all":
        peek = list(range(len(orders["peak"])))
    elif not isinstance(peak, list):
        peek = [int(x) for x in str(peak)]
    else:
        peek = peak

    out = []
    for x in peek:
        if "hkl" in orders["peak"][x]:
            hkl = orders["peak"][x]["hkl"]
            if hkl == "0" or hkl == 0:
                hkl = "000"
        else:
            hkl = "000"

        if not isinstance(hkl, list) and string == True:
            # string in, string out
            out.append(str(hkl))
        elif isinstance(hkl, list) and string == True:
            # convert numbers to string
            out_tmp = ""
            for y in range(len(hkl)):
                out_tmp = out_tmp + str(hkl[y])
            out.append(out_tmp)
        elif isinstance(hkl, list) and string == False:
            # list in, list out
            out.append(hkl)
        elif not isinstance(hkl, list) and string == False:
            # convert string to array
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
            out_tmp = [h, k, l]
            if len(hkl) > pos:
                if hkl[pos] == "-":
                    m = hkl[pos : pos + 2]
                    pos = pos + 2
                else:
                    m = hkl[pos : pos + 1]
                    pos = pos + 1
                out_tmp.append(m)
            out.append(out_tmp)

    return out


def peak_phase(orders, peak="all"):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param string: switch indicting if the output is a string or numeric list of hkl.
    :param peak: switch to add restricted peaks to the output string
    :return: string or list of peak hkl.
    """
    # possible hkl input formats:
    #   string -- hkls as a string, generally used for hkls with leading 0 or negative hkl indicy
    #   int -- hkls all positive and no leading 0
    #   list -- e.g. [h,k,l]. New format.
    # desired outputs:
    #   string -- hkls as a string using String==True
    #   list -- e.g. [h,k,l]. New format. using String==False

    if peak == "all":
        peek = list(range(len(orders["peak"])))
    elif not isinstance(peak, list):
        peek = [int(x) for x in str(peak)]
    else:
        peek = peak

    out = []
    for x in peek:
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

    # print(image_name)

    if isinstance(image_name, list):
        t_f_str = os.path.split(image_name[0])[1]
        t_f_str, _ = os.path.splitext(t_f_str)
        t_f_str += joint
        t_f_str += image_name[1][2]
        # t_f_str += joint
        # t_f_str += str(image_name[1][1])

    else:
        _, t_f_str = os.path.split(image_name)
        t_f_str, _ = os.path.splitext(t_f_str)

    return t_f_str


def make_outfile_name(
    base_filename,
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

    # if the file type is h5 the name arrives as a list of bits
    if isinstance(base_filename, list):
        # root_name = base_filename[0]
        filename = title_file_names(image_name=base_filename, string=False)

        # get the file extension
        _, ending = os.path.splitext(base_filename[0])
        ending = ending[1:]
    else:
        filename, ending = os.path.splitext(base_filename)
        ending = ending[1:]

    # strip directory if it is in the name and there is a new directory
    if directory or directory == "":
        # base_filename = os.path.basename(base_filename)
        filename = os.path.basename(filename)

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
        filename = filename + "__" + peak_string(orders, fname=True, peak=peak)
    if additional_text:  # add additional text -- i.e. notes from orders
        filename = filename + "__" + additional_text
    if orders and "note" in orders:  # add additional text from note in orders.
        filename = (
            filename + "__" + "".join(i for i in orders["note"] if i not in "\/:;*?<>|")
        )
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


def number_to_string(number, replace=".", withthis="pt"):
    """
    Turns a number into a string and then replaces the decimal place with a "pt".
    """
    number = str(number)
    number = re.sub(r"\.+", "pt", number)
    return number


def licit_filename(fname, replacement="=="):
    """
    This makes sure that a file name generated from a string is licit.
    It replaces the illegal characters [<>:/\|?*] with a replacement character
    It also replaces all '.' with 'pt' -- assuming any occurance is a number.

    after https://gist.github.com/AaronLaw/a936bebfbbd691fc954252444767e6de -- Find NTFS illegal characters in black list and rename filename.
    """
    blacklist = r"[<>:/\|?*]"

    fname = re.sub(blacklist, replacement, fname)

    fname = number_to_string(fname)

    return fname
