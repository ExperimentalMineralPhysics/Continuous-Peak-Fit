#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
from copy import deepcopy

import cv2  # opencv-python
import h5py
import hdf5plugin  # have to import it to read the files. not sure why
import matplotlib.pyplot as plt
import numpy as np

import cpf.XRD_FitPattern as fp
from cpf.IO_functions import (
    StartStopFilesToList,
    licit_filename,
    make_outfile_name,
    title_file_names,
)
from cpf.util.logging import get_logger

logger = get_logger("cpf.h5_functions")


# copied from internet somewhere August 2022
def allkeys(obj):
    # =================
    # FIXME: this function is used for the original hdf5 processing formats.
    # FIXME: it has been replaced with the new dictionary based format.
    # FIXME: the old format is now parsed through 'update_key_structure' and this functions should not be called.
    # FIXME: will be removed in a furture update.
    # =================
    "Recursively find all keys in an h5py.Group."
    keys = (obj.name,)
    if isinstance(obj, h5py.Group):
        for key, value in obj.items():
            if isinstance(value, h5py.Group):
                keys = keys + allkeys(value)
            else:
                keys = keys + (value.name,)
    return keys


# copied from internet somewhere August 2022
def getkeys(obj, key_list, level):
    # =================
    # FIXME: this function is used for the original hdf5 processing formats.
    # FIXME: it has been replaced with the new dictionary based format.
    # FIXME: the old format is now parsed through 'update_key_structure' and this functions should not be called.
    # FIXME: will be removed in a furture update.
    # =================
    keys = list(obj[key_list[0]].keys())
    keys = (obj.name,)
    if isinstance(obj, h5py.Group):
        for key, value in obj.items():
            keys = keys + (value.name,)

            # if isinstance(value, h5py.Group):
            #     keys = keys + allkeys(value)
            # else:

    return keys


def image_key_validate(h5key_list, key_start, key_end, key_step):
    """
    Validates the lengths and structure of the h5 file key values.
    Does not validate if all the keys exist.


    Parameters
    ----------
    h5key_list : TYPE
        DESCRIPTION.
    key_start : TYPE, optional
        DESCRIPTION. The default is key_start.
    key_end : TYPE, optional
        DESCRIPTION. The default is key_end.
    key_step : TYPE, optional
        DESCRIPTION. The default is key_step.

    Returns
    -------
    None.

    """

    # =================
    # FIXME: this function is used for the original hdf5 processing formats.
    # FIXME: it has been replaced with the new dictionary based format.
    # FIXME: the old format is now parsed through 'update_key_structure' and this functions should not be called.
    # FIXME: will be removed in a furture update.
    # =================

    fail = 0

    number_indicies = len(h5key_list)

    if len(key_start) != (1 | number_indicies):
        logger.error(
            " ".join(map(str, [("There are the wrong number of start positions")]))
        )
        fail = +1

    if len(key_end) != (1 | number_indicies):
        logger.error(
            " ".join(map(str, [("There are the wrong number of end positions")]))
        )
        fail = +1

    if len(key_step) != (1 | number_indicies):
        logger.error(
            " ".join(map(str, [("There are the wrong number of step lengths")]))
        )
        fail = +1

    if fail != 0:
        err_str = "The h5 settings do not balance. They need fixing."
        logger.error(" ".join(map(str, [(err_str)])))
        raise ValueError(err_str)

    return fail


def unique_labels(labels, i=3, number_data=None):
    """
    Takes the label array and turns it into strings as sensible lengths.
    The array return will be of unique values if the input array was unique,

    Parameters
    ----------
    labels : TYPE
        the list of labels to be rounded/cleaned up.
    i : number
        The number of decimal places to start the rounding at.

    Returns
    -------
    Labels rounded

    """
    number_labels = np.size(labels)
    number_unique = len(np.unique(labels))

    if number_labels == number_unique and (
        number_data == number_labels or number_labels == 1
    ):
        done = 0
        while done == 0:
            try:  # catch the unique labels not being numbers
                labels_rounded = np.round(labels, i)
                if len(np.unique(labels_rounded)) == number_labels:
                    done = 1
                    labels = np.round(labels, i + 1)
                    labels = np.array(labels)
                else:
                    i += 1
            except:
                done = 1
                pass

    else:
        err_str = "There are insufficient unique labels for the h5 file levels."
        logger.error(" ".join(map(str, [(err_str)])))
        raise ValueError(err_str)

    return labels


def get_image_keys(
    datafile,
    h5key_list,
    h5key_names,
    key_start=0,
    key_end=-1,
    key_step=1,
    key_route="",
    bottom_level="iterate",
    index=0,
    key_str="",
):
    # =================
    # FIXME: this function is used for the original hdf5 processing formats.
    # FIXME: it has been replaced with the new dictionary based format.
    # FIXME: the old format is now parsed through 'update_key_structure' and this functions should not be called.
    # FIXME: will be removed in a furture update.
    # =================

    # validate the inputs
    # image_key_validate(h5key_list, key_start=key_start, key_end=key_end, key_step=key_step)

    if isinstance(datafile, str):
        datafile = h5py.File(datafile)

    keys = []
    sep1 = "_"
    sep2 = "="
    if len(h5key_list) == 1:
        # we are at the end of the list and need to find out how many values are here, or what the values are
        key = key_route + "/" + h5key_list[0]

        # make sure the key exists
        if not key in datafile.keys():
            err_str = "The key, '%s' does not exist in '%s'" % (key, key_route)
            logger.warning(" ".join(map(str, [(err_str)])))
        else:
            number_data = datafile[key].shape[index]

            if key_start[0] < 0:
                key_start[0] = number_data + key_start[0]
            if key_end[0] < 0:
                key_end[0] = number_data + key_end[0]

            # if the order is reverse switch the values
            # FIXME: this only works for the input as set needs to be updatedto work for all logics.
            # FIXME: Should use cpf.IO_functions/file_list work for this.
            if key_step[0] < 0:
                key_start[0], key_end[0] = key_end[0], key_start[0]
                key_end[0] -= 1
                # key_start[0] -= 1
            else:
                key_end[0] += 1

            if bottom_level == "sum":
                index_values = list([*range(key_start[0], key_end[0], key_step[0])])

                # get the labels.
                # FIXME: this doesnt work if the last enty in the h5key_names list is blank.
                # FIXME: Should use a function to cover all options and combine with iterate call.
                labels = key_route + datafile[key_route + "/" + h5key_names[0]][()]
                # FIXME! this is the line to use if the last entry in the list in empty
                # labels = str(key_str)

                return [[key, index_values, labels]]

            elif bottom_level == "iterate":
                # iterate over the size of the array in the h5 group.

                # make list of key labels
                key_strings_new = get_image_key_strings(
                    datafile,
                    key_route=key_route,
                    key_names=h5key_names[0],
                    key_measure=h5key_list[0],
                    key_start=0,
                    key_end=-1,
                    key_step=1,
                    index=index,
                )
                out = []

                if number_data != len(key_strings_new):
                    err_str = (
                        "The number of data (%i) does not match the number of keys (%i). "
                        "It is not possible to continue until this error is corrected."
                        % (number_data, len(key_strings_new))
                    )
                    raise ValueError(err_str)

                # get the labels:
                lbl_temp = []
                for i in [*range(key_start[0], key_end[0], key_step[0])]:
                    lbl_temp.append([key, i, key_str + sep1 + key_strings_new[i]])

                return lbl_temp

            else:
                err_str = "This h5 process is not recognised."
                raise ValueError(err_str)
    else:
        # step inside the key_list

        # make current key
        key_str = key_str + "/" + h5key_list[0]

        # make list of key labels
        key_strings = get_image_key_strings(
            datafile,
            key_route=key_str,
            key_names=h5key_names[0],
            key_measure=h5key_list[0],
            key_start=0,
            key_end=-1,
            key_step=1,
            index=0,
        )

        # make list of values
        if isinstance(h5key_names[0], str) and len(h5key_names[0]) == 0:
            key_lst = list(datafile[h5key_list[0]].keys())
        else:
            key_lst = list(datafile[h5key_list[0]].keys())

        # cut to what we want ot iterate over
        # if we are using all the data make sure we run to the end.
        if key_start[0] < 0:
            key_start[0] = len(key_lst) + key_start[0]
        if key_end[0] < 0:
            key_end[0] = len(key_lst) + key_end[0]

        # if we are looking for more images than there are cut the list
        # key_end[0] = np.min(number_data, key_end[0])

        # if the order is reverse switch the values
        if key_step[0] < 0:
            key_start[0], key_end[0] = key_end[0], key_start[0]
            key_end[0] -= 1
            # key_start[0] -= 1
        else:
            key_end[0] += 1
        # populate the list
        for i in [*range(key_start[0], key_end[0], key_step[0])]:
            # use the index to pass either the index of the array or the upper level group index
            if bottom_level == "label":
                if index == 0:
                    index = str(key_lst[i])
                else:
                    index = key_lst[i]

            keys.extend(
                get_image_keys(
                    datafile,
                    h5key_list[1:],
                    h5key_names[1:],
                    key_start=key_start[1:],
                    key_end=key_end[1:],
                    key_step=key_step[1:],
                    key_route=key_lst[i],
                    bottom_level=bottom_level,
                    index=index,
                    key_str=str(key_strings[i]),
                )
            )

    datafile.close()
    return keys


def get_image_key_strings(
    datafile,
    key_route="",
    key_names=[""],
    key_measure=[""],
    key_start=0,
    key_end=-1,
    key_step=1,
    key_str="",
    index=0,
    sep1="_",
    sep2="=",
):
    """
    Makes the key labels for use in the filenames and so forth.
    Trasposes the key_names into strings, numberical values or reads from other groups depending on the content.
    The key names is a list of lists or strings.
    There are three different actions depending on the string:
    '' - empty string -- adds the numberical position of the key as a label
    '/' - slash -- adds the name of the key as a label
    '...' - key name -- adds keyname=keyvalue as a label.
    Multiple keynames can be concatiated within a list and all will be added to the label.

    Parameters
    ----------
    data_file : file handle
        Filehandle for h5 data file.
    key_route : TYPE, optional
        DESCRIPTION. The default is "".
    key_names : TYPE, optional
        DESCRIPTION. The default is [""].
    key_start : TYPE, optional
        DESCRIPTION. The default is 0.
    key_end : TYPE, optional
        DESCRIPTION. The default is -1.
    key_step : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    key_str : list
        List of labels for the images in the data group.

    """

    # =================
    # FIXME: this function is used for the original hdf5 processing formats.
    # FIXME: it has been replaced with the new dictionary based format.
    # FIXME: the old format is now parsed through 'update_key_structure' and this functions should not be called.
    # FIXME: will be removed in a furture update.
    # =================

    if not isinstance(key_names, list):
        key_names = [key_names]
    number_data_tmp1 = []
    number_data_tmp = []
    for i in range(len(key_names)):
        if (
            isinstance(datafile[key_route + "/" + key_names[i]], h5py.Group)
            and key_names[i] == "/"
        ):
            logger.debug(" ".join(map(str, [(i, "/")])))
            number_data_tmp.append(
                len(list(datafile[key_route + "/" + key_names[i]].keys()))
            )

        elif (
            isinstance(datafile[key_route + "/" + key_names[i]], h5py.Group)
            and key_names[i] == ""
        ):
            try:
                if key_measure == "/":
                    number_data_tmp.append(
                        len(list(datafile[key_route + "/" + key_names[i]].keys()))
                    )
                else:
                    attr = datafile[key_route + "/" + key_measure]
                    tmp = datafile.get(key_route + "/" + key_measure)
                    number_data_tmp.append(tmp.shape[index])

            except:
                logger.debug(" ".join(map(str, [(i, "value")])))
                number_data_tmp.append(0)

        else:
            number_data_tmp.append(datafile[key_route + "/" + key_names[i]].size)

    number_data = np.max(number_data_tmp)

    if isinstance(key_names, str):
        key_names = [key_names]

    # get the labels
    labels = []
    for i in range(len(key_names)):
        if key_names[i] == "":
            # if empty then list numbers.
            labels_temp = np.array(range(number_data))
            # FIXME: We are zero counting the images and the indecies. It might be better to 1 count them.
            # if so to 1 count the indicies we add 1 to the prewvious line.
            # the counting over the arrays needs to be done in get_image_keys
            labels_temp = unique_labels(labels_temp, number_data=number_data)
        elif key_names[i] == "/":
            # if "/" then list names of subgroups
            labels_temp = list(datafile[key_route + "/" + key_names[i]].keys())
        else:
            # it is a key and so list the key contents.
            try:
                # number_data = datafile[key_route + "/" + key_names[i]].shape[index]
                labels_temp = datafile[key_route + "/" + key_names[i]][()]
            except:
                # catch incase 1 iterable index is the length of the data and the other is only 1.
                labels_temp = datafile[key_route + "/" + key_names[i]][()]
            labels_temp = unique_labels(labels_temp, number_data=labels_temp.size)
        labels.append(labels_temp)

    # if we are using all the data make sure we run to the end.
    if key_end == -1:
        key_end = number_data - 1
    # else:
    #     key_end += 1
    # if we are looking for more images than there are cut the list
    # key_end[0] = np.min(number_data, key_end[0])

    # if the order is reverse switch the values
    # FIXME: this only works for the input as set needs to be updatedto work for all logics.
    # FIXME: Should use cpf.IO_functions/file_list work for this.
    if key_step <= -1:
        key_start, key_end = key_end, key_start
        key_end -= 1
        key_start -= 1

    # make the list of labels
    out = []
    for i in [*range(key_start, key_end + 1, key_step)]:
        lbl_str = ""
        for j in range(len(labels)):
            if isinstance(key_names[j], str) and len(key_names[j]) == 0:
                pass
            else:
                lbl_str = lbl_str + os.path.basename(
                    os.path.normpath(key_route + "/" + key_names[j])
                )

            if (
                len(key_names[j]) == 0 or key_names[j] == "/"
            ):  # then there is no name to paste into the file name string
                if np.size(labels[j]) == 1:
                    lbl_str = lbl_str + str(labels[j])
                else:
                    lbl_str = lbl_str + str(labels[j][i])
            elif labels[j].size == 1:
                lbl_str = lbl_str + sep2 + str(labels[j])
            else:
                if np.size(labels[j]) == 1:
                    # if labels are a single item they can come wrapped in lists.
                    # not sure why but this is here to capture it.
                    labels = [item for items in labels for item in items]

                    lbl_str = lbl_str + sep2 + str(labels[j])
                else:
                    lbl_str = lbl_str + sep2 + str(labels[j][i])
            if j != len(key_names) - 1:  # np.size(labels[j]):
                lbl_str = lbl_str + sep1

        # IO.licit_filename(lbl_str)
        if len(key_str) != 0:
            lbl_str = key_str + sep1 + lbl_str

        lbl_str = licit_filename(lbl_str)

        out.append(lbl_str)

    return out


# %%


def DefaultProcessDictionary(types=False):
    """
    Definition of default process dictionary for iterating over hdf5 files.

    Either it returns the default dictionary, if types==False:
       {"do":   "iterate",
        "from":  0,
        "to":    -1,
        "step":  1,
        "using": "position",
        "label": ["", "*"]}
    or it returns a dictionary with the expected formats or values:
        {"do":    {"type": str,
                   "values": ["sum", "iterate"]},
        "from":  {"type": (int, float, np.ndarray)},
        "to":    {"type": (int, float, np.ndarray)},
        "step":  {"type": (int)},
        "dim":   {"type": (int)},
        "using": {"type": str,
                  "values": ["value", "position"]},
        "label": {"type": (list, str)},
        }

    Returns
    -------
    dict
        Default dictionary for h5 key descriptions.

    """
    if types == False:
        return {
            "do": "iterate",
            "from": 0,
            "to": -1,
            "step": 1,
            "dim": 0,
            "using": "position",
            "label": ["", "*"],
        }
    else:
        return {
            "do": {"type": str, "values": ["sum", "iterate"]},
            "from": {"type": (int, float, np.ndarray)},
            "to": {"type": (int, float, np.ndarray)},
            "step": {"type": (int)},
            # "list":   {"type": list}, # this is possible but not in the detault set.
            "dim": {"type": (int)},
            "using": {"type": str, "values": ["value", "position"]},
            "label": {"type": (list, str)},
        }


def update_key_structure(fit_parameters, fit_settings):
    """
    For backwards compatability. Can be removed in future release.

    Convert the original hdf5 key structure into the new (improved?) structure.

    1. Key for data to be iterated over:
        "h5_key_list" (old)     ==> "h5datakey" (new)
        ['/', 'measurement/p3'] ==> '/*/measurement/p3'

    2. Indicies to iterate over and labels are now a list of dictionaries
        "h5_key_names" & 'h5_key_start' & 'h5_key_end' & 'h5_key_step' &
        'h5_data' ==> 'h5iterations'

        h5_key_names = [["/",''],
                        ["", 'instrument/positioners/dz1','instrument/positioners/dy']]
        h5_key_start = [0, 0]
        h5_key_end   = [42,-1]
        h5_key_step  = [-1,1]
        h5_data      = 'iterate'
            ==>
        h5iterations = [{"from": 0,
                       "to": 42,
                       "step": -1,
                       "using":"value",
                       "label":["*", 'pos']}
                      {"do":"iterate",
                       "from": 0,
                       "to": -1,
                       "step": 1,
                       "using":"position",
                       "label":['*','/*/instrument/positioners/dz1',
                                '/*/instrument/positioners/dy']}]

    There is no provision in the old system for a list in keys, than may or may not then be trucated.

    Parameters
    ----------
    fit_parameters : cpf parameeter list

    fit_settings : cpf settings class


    Returns
    -------
    h5datakey : str
        Data key for hdf5 file's data
    h5iterations : list
        list of dictionaries for how to loop over the hdf5 file.


    """
    # convert old "h5_key_list" into new "h5datakey".
    h5datakey = ""
    for i in range(len(fit_settings.h5_key_list)):
        if fit_settings.h5_key_list[i][0] != "/":
            h5datakey += "/"
        h5datakey += fit_settings.h5_key_list[i]
        if i < len(fit_settings.h5_key_list) and i != len(fit_settings.h5_key_list) - 1:
            h5datakey += "*"

    # convert everything else into a dictionary.
    h5iterations = []
    lowerKey = ""
    for i in range(len(fit_settings.h5_key_list)):
        # make 1 dictionary for each level.
        loop = DefaultProcessDictionary()

        loop["from"] = fit_settings.h5_key_start[i]
        loop["to"] = fit_settings.h5_key_end[i]
        loop["step"] = fit_settings.h5_key_step[i]

        if i < len(fit_settings.h5_key_list):
            loop["using"] = "value"
        else:
            loop["using"] = "position"

        if isinstance(fit_settings.h5_key_names, str):
            fit_settings.h5_key_names == [fit_settings.h5_key_names]
        loop["label"] = []
        for j in range(len(fit_settings.h5_key_names[i])):
            if fit_settings.h5_key_names[i][j] == "":
                loop["label"].append("pos")
            elif fit_settings.h5_key_names[i][j] == "/":
                loop["label"].append("value")
            else:
                loop["label"].append(lowerKey + fit_settings.h5_key_names[i][j])
        lowerKey += "*/"
        h5iterations.append(loop)

    return h5datakey, h5iterations


def image_key_validate_new(
    fit_parameters=None, fit_settings=None, h5_iterate=None, end_if_errors=False
):
    """
    Validates the lengths and structure of h5_iterate, which is the list/dictionary
    of h5 file key values and labels.

    It does not though validate if the keys exist is the hdf5 files.

    This script calls the detauls structures defined in 'DefaultProcessDictionary'

    Parameters
    ----------
    fit_parameters : continuous peak fit parameter list, optional
    fit_settings : continuous peak fit settings class, optional
        continuous peak fit settings class with h5 description in it.
    h5_iterate : dictionary
        containing information used to iterare over the h5 files,
    end_if_errors : bool, optional
        Switch to stop execution and spit out all an error. The default is False.

    Returns
    -------
    errors : list
        A list of all the problems with the h5 descriptions.
    h5_iterate : list
        An updated version of the input list/dictionary.
    """

    # get expected data types
    expected_entries = DefaultProcessDictionary(types=True)
    expected_keys = list(DefaultProcessDictionary().keys())

    errors = []
    warnings = []
    if h5_iterate == None:
        try:
            h5_iterate = fit_settings.h5_iterate
        except:
            errors.append("h5_iterate is not in the settings")

    if isinstance(h5_iterate, dict) == True:
        h5_iterate = [h5_iterate]

    if "h5_iterate" != None:
        num_iter = len(h5_iterate)
        if fit_settings:
            number_indicies = len(
                [a.start() for a in list(re.finditer("\*", fit_settings.h5_datakey))]
            )
            if not (num_iter == number_indicies or num_iter == number_indicies - 1):
                errors.append(
                    "The number of iterations is not the same as, or 1 less than, the number of indicies in the data"
                )
        for i in range(num_iter):
            dict_keys = list(h5_iterate[i].keys())
            # check keys that are present
            for j in range(len(dict_keys)):
                if dict_keys[j] in expected_entries:
                    if expected_entries[dict_keys[j]]["type"] == str:
                        if (
                            not h5_iterate[i][dict_keys[j]]
                            in expected_entries[dict_keys[j]]["values"]
                        ):
                            errors.append(
                                f"h5_iterate[{i}]: the '{dict_keys[j]}' value must be one of {expected_entries[dict_keys[j]]['values']}"
                            )
                        else:
                            pass
                    else:
                        if not isinstance(
                            h5_iterate[i][dict_keys[j]],
                            expected_entries[dict_keys[j]]["type"],
                        ):
                            errors.append(
                                f"h5_iterate[{i}]: the '{dict_keys[j]}' value must be a {expected_entries[dict_keys[j]]['type']}"
                            )
                        else:
                            pass

            # check for missing keys
            for j in range(len(expected_keys)):
                if expected_keys[j] not in dict_keys:
                    warnings.append(
                        f"h5_iterate[{i}]: does not contain '{expected_keys[j]}'; a default value has been added"
                    )
                    h5_iterate[i][expected_keys[j]] = DefaultProcessDictionary()[
                        expected_keys[j]
                    ]

            # check for extra keys
            extra = [x for x in dict_keys if x not in expected_entries]
            if len(extra) != 0:
                errors.append(
                    f"h5_iterate[{i}]: there are unrecognised entries: {extra}. These will not be used."
                )

    for i in range(len(errors)):
        logger.error(" ".join(map(str, [(errors[i])])))
    for i in range(len(warnings)):
        logger.warning(" ".join(map(str, [(warnings[i])])))
    if end_if_errors == True and len(errors) > 0:
        err_str = "The h5 settings do not balance. They need fixing."
        logger.error(" ".join(map(str, [(err_str)])))
        raise ValueError(err_str)

    if errors == []:
        errors = False

    return errors, h5_iterate


def h5_tree(val, keyup="", list_of_keys=[], recurse=True):
    """
    Recursively list all the keys in the h5 file.

    Can be slow especially if the files are big.

    edited after: https://stackoverflow.com/questions/61133916/is-there-in-python-a-single-function-that-shows-the-full-structure-of-a-hdf5-fi

    Parameters
    ----------
    val : file or h5py.Group
        The object to in iterated over now.
    keyup : string, optional
        strong for keys at levels above the one being iterated over now.
        The default is "".
    list_of_keys : list, optional
        List of all the keys found in the file so far. The default is [].
    recurse : bool, optional
        switch to just look at current level or all levels below. The default is True.

    Returns
    -------
    list_of_keys : list
        list of keys in hdf5 file.
    """
    for key, val in val.items():
        # if group iterate over it, otherwise it will be a dataset
        if type(val) == h5py._hl.group.Group:
            if recurse == True:
                h5_tree(val, keyup=keyup + "/" + key)
            else:
                list_of_keys.append(keyup + "/" + key)
        else:
            list_of_keys.append(keyup + "/" + key)
    return list_of_keys


def get_image_keys_new(datafile, h5key_data, h5_iterate, sep1="_", sep2="="):
    """
    Makes a list of all the images addressed in the h5file 'datafile', by 'h5_iterate'.
    This is used by continuous peak fit as the list of images to iterate over.

    The list contains a unique label for each image; the label is either extracted from
    the h5 file using h5_iterate["label"] or listed as a number.

    The function returns a list of images, each element of which is:
        ['h5/file/key', [number list for images], 'uniquielabelforimage']

    the expected format of 'h5_iterate' is defined in DefaultProcessDictionary


    Parameters
    ----------
    datafile : str or open file pointer
        Datafile to be worked on
    h5key_data : str
        key in h5 file that contains the data.
    h5_iterate : list or dictionary
        list of doctionaryies containing the information on iteratiing over the keys.
    sep1 : str, optional
        seperator to go between labels. The default is "_".
    sep2 : str, optional
        seperator to go between key name and value in labels. The default is "=".

    Raises
    ------
    ValueError
        DESCRIPTION.
    NotImplementedError
        DESCRIPTION.

    Returns
    -------
    out : list
        List of the labels, keys and data number.

    """

    # set some defaults for regualr expreassions
    regexp_num = "([-+]?[0-9]*\.[0-9]+|[-+]?[0-9]+)"
    regexp_end = "$"

    # validate the inputs
    _, h5_iterate = image_key_validate_new(h5_iterate=h5_iterate, end_if_errors=True)
    if isinstance(h5key_data, list):
        h5key_data = h5key_data[0]
    elif not isinstance(h5key_data, str):
        raise ValueError("h5_data must be a single string")
    if os.path.isfile(datafile) == False:
        raise ValueError("The hdf5 file does not exist")

    # open file if needed
    if isinstance(datafile, str):
        df = h5py.File(datafile)
    else:
        df = datafile

    # get number of searches (i.e. number *) in datakey
    loops = [a.start() for a in list(re.finditer("\*", h5key_data))]
    dividers = [a.start() for a in list(re.finditer("/", h5key_data))]
    if len(loops) > 1:
        raise NotImplementedError()

    # find all keys that match the wildcards.
    # this is necessary incase we are looping over all valid entries without
    # directly knowing the key names
    keylist = []
    for i in range(len(loops)):
        # FIXME this loop should be an iterative loop. so thatit can iterate
        # over as many keys as required.

        # find part of key before wildcard.
        h5key_data_sub = h5key_data[: max(j for j in dividers if j < loops[i]) + 1]
        # get all matching keys
        lst = h5_tree(df[h5key_data_sub], list_of_keys=[], recurse=False)
        # FIXME: the h5tree function call needs list_of_keys=[] otherwise during
        # testing when the function is run it appends to the previous list of keys.

        ## get list of values that matches the wildcard.
        # find part of key before separation after wildcard.
        # cannot assume that the wildcard is the last part of the keyname.
        h5key_data_sub = h5key_data[: min(j for j in dividers if j > loops[i])]
        # replace * with regular expresion search for number
        h5key_data_sub = re.sub(r"\*", regexp_num, h5key_data_sub)
        # match lst to subkey.
        lst = list(filter(re.compile(h5key_data_sub).match, lst))

        ## get keys that match that have subkeys matching h5key_data
        # cannot assume that all the keys have the same subdata.

        # capture numerical values in keys.
        vals = []
        for i in range(len(lst)):
            vals_tmp = re.search(h5key_data_sub, lst[i])
            # get values from regexpression
            vals.append(list(vals_tmp.groups()))

        # check if the key exists.
        for j in range(len(vals)):
            # replace * with number using regular expresion
            h5key_data_sub = re.sub(r"\*", vals[j][0], h5key_data)
            try:
                # if the data key does not exist it will skip over it.
                df[h5key_data_sub]
                keylist.append(h5key_data_sub)
            except:
                pass

    if 0:  # for debugging
        # get all keys in h5 file [CAUTION could be slow if files are large]
        key_list_all = h5_tree(df)

        # setup h5key_data to be used.
        # replace * with number using regular expresion
        h5key_data = re.sub(r"\*", regexp_num, h5key_data)
        # make sure ends with a $
        if h5key_data[-1] != "$":
            h5key_data = h5key_data + regexp_end

        # get list of matching keys
        keylist = list(filter(re.compile(h5key_data).match, key_list_all))

    # capture numerical values in keys.
    vals = []
    for i in range(len(keylist)):
        vals_tmp = re.search(re.sub(r"\*", regexp_num, h5key_data), keylist[i])
        # get values from regexpression
        vals.append(list(vals_tmp.groups()))

    # sort values and keys into order.
    vals_num = np.array(vals).astype(float)
    o = np.argsort(vals_num, axis=0)[:, 0]
    vals_num = vals_num[o]
    keylist = [keylist[i] for i in o]
    vals = [vals[i] for i in o]

    # cut list to size.
    present = []
    for i in range(len(loops)):
        # FIXME. need to make h5 iterate as long as ther are loops+1
        itera = h5_iterate[i]

        if itera["using"] == "value":
            # make list of values to keep from parameter dictionary

            if itera["from"] == 0:
                itera["from"] = float(vals[0][0])
            if itera["to"] == -1:
                itera["to"] = float(vals[-1][0])

            keep, _ = StartStopFilesToList(paramDict=itera)
            present = np.atleast_1d(
                np.array(
                    [
                        x
                        for x in [np.where(np.isclose(i, vals_num))[0] for i in keep]
                        if len(x) > 0
                    ]
                ).squeeze()
            )
            keylist = [keylist[x] for x in present]
            vals = [vals[x] for x in present]

        elif itera["using"] == "position":
            # make list of positions to keep from parameter dictionary

            if itera["to"] == -1:
                itera["to"] = int(len(vals)) - 1

            keep, _ = StartStopFilesToList(paramDict=itera)
            keylist = [keylist[x] for x in keep]
            vals = [vals[x] for x in keep]

        else:
            raise ValueError("unrecognased type")

    # get labels.
    labels = []
    itera = h5_iterate[0]
    for i in range(len(keylist)):
        # FIXME. need to make h5 iterate as long as ther are loops+1
        label_tmp = []
        if not isinstance(h5_iterate[0]["label"], list):
            h5_iterate[0]["label"] = list(h5_iterate[0]["label"])
        # make the labels.
        label_tmp = get_labels(
            df, itera["label"], 1, i, vals[i], sep1=sep1, sep2=sep2, key=keylist[i]
        )
        labels.append(label_tmp)

    if len(keylist) == 0:
        raise ValueError(
            "No keys have been found. Check that the indicies iterating over are in the hdf5 file."
        )

    # iterate over the keylist and expand with labels.
    out = []
    for i in range(len(keylist)):
        # the key list should contain everything that is not the bottom level
        # of the h5 keys. this loop only works on the bottom levels with the data

        # look at the last of the iteration keys.
        # deepcopy so can check for "to": -1 every time round.
        itera = deepcopy(h5_iterate[-1])

        # make sure the key exists
        if not keylist[i] in df.keys():
            err_str = "The key, '%s' does not exist in '%s'" % (keylist[i], df)
            logger.warning(" ".join(map(str, [(err_str)])))
        else:
            number_data = df[keylist[i]].shape[itera["dim"]]
            if itera["do"] == "sum":
                # index_values = list([*range(itera["start"], itera["stop"], itera["step"])])
                if itera["to"] == -1:
                    itera["to"] = number_data - 1
                index_values, _ = StartStopFilesToList(paramDict=itera)
                # get the labels -- only need labels from layers above because summing the data.
                lbls = labels[i]
                # lbls = get_labels(df, itera["label"],  number_data, j, vals[i], sep1=sep1, sep2=sep2, key=labels[i])
                out.append([keylist[i], index_values, lbls])

            elif itera["do"] == "iterate":
                # iterate over the size of the array in the h5 group.
                if itera["to"] == -1:
                    itera["to"] = number_data - 1
                for j in np.arange(itera["from"], itera["to"] + 1, itera["step"]):
                    index_values = j
                    lbls = labels[i]
                    additional_label = get_labels(
                        df,
                        itera["label"],
                        number_data,
                        j,
                        vals[i],
                        sep1=sep1,
                        sep2=sep2,
                        key=labels[i],
                    )
                    # print(additional_label)
                    if len(additional_label) != "":
                        lbls += f"{sep1}{additional_label}"

                    # if itera["using"] == "value" and :
                    out.append([keylist[i], index_values, lbls])

                # for j in range(number_data):
                #     index_values = j

                #     lbls = labels[i]
                #     additional_label = get_labels(df, itera["label"],  number_data, j, vals[i], sep1=sep1, sep2=sep2, key=labels[i])
                #     # print(additional_label)
                #     if len(additional_label) != "":
                #         lbls += f"{sep1}{additional_label}"

                #     # if itera["using"] == "value" and :
                #     out.append([keylist[i], index_values, lbls])

            else:
                err_str = "This h5 process is not recognised."
                raise ValueError(err_str)

    df.close()
    return out


def get_labels(
    datafile, label_list, num_entries, seq, vals, sep1="_", sep2="=", key=None
):
    """
    Makes a label string using keys and values listed in "label_list".

    The possible entries in the label_list are:
        'pos': returns the position in the loop
        'value': returns the numeric value of the key's address
        'key/in/h5_file': returns value in key, labelled with the end key.
                            (i.e. h5_file=0.1234)
    Each entry is separated in the returned string by 'sep1' (default='_')
    If the label is from a h5key then the keyname and value are separated by
        'sep2' (default='=')
    The returned string is checked for illicit filename characters using
        cpf.IO_functions.licit_filename.

    Parameters
    ----------
    datafile : h5file pointer
        H5 file for which the labels are being made.
    label_list : list
        List containing keys of required labels.
    seq : List,
        DESCRIPTION. The default is None.
    vals : List
        DESCRIPTION. The default is None.
    sep1 : string, optional
        seperator to go between different labels. The default is "_".
    sep2 : string, optional
        seperator to go between label name and value. The default is "=".

    Returns
    -------
    label : string
        String of labels.

    """
    label = ""
    if not isinstance(label_list, list):
        label_list = list(label_list)

    for j in range(len(label_list)):
        if j != 0:
            label += sep1
        lb = label_list[j]
        if lb.lower() == "pos":
            # looking for position in the sequence
            label += str(seq)
        elif lb.lower() == "value":
            # looking for position in the data file
            if isinstance(vals, list):
                label += str(vals[0])
            else:
                label += str(vals)
        elif "key" in lb.lower():
            # add key value to the label
            label += str(key)
        else:
            # look up key value
            lb = re.sub(r"\*", vals[0], lb)
            try:
                ky = datafile[lb][()]

                # ky = unique_labels(ky)

                # if number_data != len(ky):
                #     err_str = (
                #         "The number of data (%i) does not match the number of keys (%i). "
                #         "It is not possible to continue until this error is corrected."
                #         % (number_data, len(key_strings_new)))
                #     raise ValueError(err_str)
                if isinstance(ky, list):
                    ky = ky[seq]
                elif isinstance(ky, np.ndarray):
                    ky = ky[seq]
                elif isinstance(ky, bytes):
                    ky = "".join(map(chr, ky))
                else:  # is string, or something else.
                    pass
                div = [a.start() for a in list(re.finditer("/", lb))]
                if lb[-1] == "/":
                    label += licit_filename(f"{lb[div[-2]+1:div[-1]]}{sep2}{ky}")
                else:
                    label += licit_filename(f"{lb[div[-1]+1:]}{sep2}{ky}")
            except:
                err_str = "The key, '%s' does not exist in '%s, for dataframe %s'" % (
                    lb,
                    datafile,
                    seq,
                )
                logger.warning(" ".join(map(str, [(err_str)])))
                raise ValueError(err_str)

    return label


# %%


def get_images(
    image_list=None,
    settings_file=None,
    settings_class=None,
    image_num=None,
    debug=False,
):
    """
    Export images from h5 files for use with Dioptas.

    Parameters
    ----------
    image_list : LIST
        This is a list that contains the file names and h5 location of the images.
    settings_file : TYPE
        DESCRIPTION.
    export_type : TYPE, optional
        DESCRIPTION. The default is "tiff".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    data : numpy array
        The image(s) data array.

    """

    if image_list != None:
        if len(image_list) == 2 and len(image_list[1]) == 3:
            # then the list is for a single file and needs to be wrapped in another list
            image_list = [image_list]
    elif settings_class != None:
        settings_for_fit = settings_class
        image_list = settings_for_fit.image_list
    elif settings_file != None:
        settings_for_fit = fp.initiate(settings_file)
        image_list = settings_for_fit.image_list
    else:
        raise ValueError(
            "There are no settings or image_list. The fitting cannot proceed until a recognised "
            "settings file or class is present."
        )

    if image_num == None:
        image_num = list(range(len(image_list)))
    elif isinstance(image_num, int):
        image_num = [image_num]
    # image_num could also be a list -- in which case leave it alone.

    # get image
    for i in range(len(image_num)):
        n = image_num[i]
        datafile = h5py.File(image_list[n][0], "r")
        datakey = image_list[n][1][0]
        data_position_in_key = image_list[n][1][1]
        data_tmp = np.array(datafile[datakey][data_position_in_key])

        if len(image_num) == 1:
            data = data_tmp
        else:
            axis = 2
            if i == 0:
                data = data_tmp
                data = np.expand_dims(data, axis=axis)
            else:
                data = np.append(data, np.expand_dims(data_tmp, axis=axis), axis=axis)

    return data


def plot_images(
    settings_file=None,
    settings_class=None,
    image_num=None,
    export_type="tiff",
    debug=False,
):
    """
    Export images from h5 files for use with Dioptas.

    Parameters
    ----------
    settings_file : TYPE
        DESCRIPTION.
    export_type : TYPE, optional
        DESCRIPTION. The default is "tiff".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    # FIX ME: there should be an option to override the settings file and just import a subset of the images -- with more precision than just an image_num

    # FIX ME: the plotting of the imags should be done through the data class.

    if settings_class != None:
        settings_for_fit = settings_class
    elif settings_file != None:
        settings_for_fit = fp.initiate(settings_file)
    else:
        raise ValueError(
            "There are no settings. The fitting cannot proceed until a recognised "
            "settings file or class is present."
        )

    if image_num == None:
        image_num = list(range(settings_for_fit.image_number))

    for i in range(len(image_num)):
        n = image_num[i]
        im_data = get_images(settings_class=settings_for_fit, image_num=n)
        plt.imshow(np.log10(im_data))
        # plt.imshow(np.log10(im_data[0,:,:]))
        # plt.title(IO.make_outfile_name(os.path.split(settings_for_fit.image_list[n][0])[1],additional_text=settings_for_fit.image_list[n][1][2]))
        plt.title(title_file_names(settings_for_fit=settings_for_fit, num=i))
        plt.show()


def save_images(
    settings_file=None,
    settings_class=None,
    image_num=None,
    export_type="tif",
    debug=False,
):
    """
    Export images from h5 files for use with Dioptas.

    Parameters
    ----------
    settings_file : TYPE
        DESCRIPTION.
    export_type : TYPE, optional
        DESCRIPTION. The default is "tiff".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    # FIX ME: ther should be an option to override the settings file and just import a subset of the images -- with more precision than just an image_num

    if settings_class != None:
        settings_for_fit = settings_class
    elif settings_file != None:
        settings_for_fit = fp.initiate(settings_file)
    else:
        raise ValueError(
            "There are no settings. The fitting cannot proceed until a recognised "
            "settings file or class is present."
        )

    if image_num == None:
        image_num = list(range(settings_for_fit.image_number))

    for i in range(len(image_num)):
        n = image_num[i]
        im_data = get_images(settings_class=settings_for_fit, image_num=n)

        fnam = make_outfile_name(
            settings_for_fit.image_list[n][0],
            additional_text=settings_for_fit.image_list[n][1][2],
            directory=settings_for_fit.output_directory,
            extension=export_type,
        )

        cv2.imwrite(fnam, im_data[0, :, :])
        logger.moreinfo(" ".join(map(str, [("Writing:", fnam)])))
