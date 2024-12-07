#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 08:02:51 2024

@author: j033794sh
"""

import re

import numpy as np


def ParseOrders(seriesOrder, param_str, component, position):
    """
    Takes the constraints in the order setttings and makes a constraint expression
    to apply in lmfit.


    Parameters
    ----------
    orders : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    expression = None
    value = 0  # np.nan
    vary = True
    if len(SeriesConstraints(seriesOrder)) == 0:  # isinstance(seriesOrder, list):
        # this is the default and the code already knows what to do
        expression = None
        vary = True
        value = np.nan

    else:
        # get list of constraints.
        Series_Constraints = SeriesConstraints(seriesOrder)
        positionsFixed = []  # empty inace
        # parse list of constrainsts for current parameter
        for i in range(len(Series_Constraints)):
            if "fixed" in Series_Constraints[i].lower():
                # fix the parameter to a value
                vary = False
                # convert value(s) to tuple
                value = tuple(
                    map(
                        float,
                        (
                            Series_Constraints[i]
                            .removeprefix("fixed")
                            .removeprefix("(")
                            .removesuffix(")")
                            .split(",")
                        ),
                    )
                )
                expression = None

            elif "equal" in Series_Constraints[i]:
                if "equal" in Series_Constraints[i]:
                    relation = ""  # '='
                elif "greater" in Series_Constraints[i]:
                    relation = ">"  # this iss not work.
                elif "less" in Series_Constraints[i]:
                    relation = "<"  # this iss not work.
                # else:
                #     relation = "" # '='

                parts = (
                    Series_Constraints[i]
                    .replace("greater", "equal")
                    .replace("less", "equal")
                    .split("equal")
                )

                # from first bit of parts determine what part of series are constrained
                if "all" in parts[0]:
                    positionsFixed = list(range(BiggestValue(seriesOrder) * 2 + 1))

                elif "order" in parts[0]:
                    # split the string up assuming the orders are split either by comma
                    # or by english suffixies for positions. e.g. 1st2nd3rd,4th
                    # the complete list of suffixes is "st", "nd", "rd" and "th".
                    split_orders = list(
                        map(
                            int,
                            parts[0]
                            .replace("order", "")
                            .replace("st", ",")
                            .replace("nd", ",")
                            .replace("rd", ",")
                            .replace("th", ",")
                            .removesuffix(",")
                            .split(","),
                        )
                    )
                    positionsFixed = []
                    for j in range(len(split_orders)):
                        if split_orders[j] == 0:
                            positionsFixed.append(0)
                        else:
                            positionsFixed.append(split_orders[j] * 2 - 1)
                            positionsFixed.append(split_orders[j] * 2)
                else:  # no order and possible list of positions
                    # split on ','
                    split_positions = parts[0].removesuffix(",").split(",")
                    # find and replce "-" with list of values
                    for j in range(len(split_positions)):
                        if "-" in split_positions[j]:
                            split_positions[j] = split_positions[j].split("-")
                            split_positions[j] = list(
                                range(split_positions[j][0], split_positions[j][1], 1)
                            )
                    # make sure the list is flat.
                    positionsFixed = [
                        item for items in split_positions for item in items
                    ]

                # from second  bit of parts determine what these are constrained to.
                parts[1] = list(
                    map(int, (parts[1].replace("peaks", "").replace("peak", "")))
                )
                if len(parts[1]) != 1:
                    # if there is more than 1 value here then it is not clear what value to constrain
                    # this value to.
                    # raise ValueError("The peak to fix this to is not unique.")

                    # fix to all
                    positionsFixed = list(range(100))

            elif (
                "delta" in Series_Constraints[i]
                or "greater" in Series_Constraints[i]
                or "less" in Series_Constraints[i]
            ):
                raise NotImplementedError(
                    "non equal constrains on parameters are not implemented"
                )

            else:
                # raise an error.
                raise ValueError("There is no recognised relationship in this string. ")

        if position in positionsFixed:
            # make param_str point to preak this one is fixed to.
            param_str_old = re.sub(r"\d+$", "", param_str) + str(parts[1][0])
            expression = relation + param_str_old + "_" + component + str(position)
    return expression, vary, value


def BiggestValue(input_list):
    """
    Returns the biggest numberical value from list with mixed contents.
    Ignores nested lists

    Parameters
    ----------
    input_list : list
        List containing mixtures of data types.

    Returns
    -------
    BiggestValue : number
        Largest value in the base of the list.
    """

    return np.max(SeriesValues(input_list))


def SeriesValues(input_list):
    """
    Returns only the numerical value from list with mixed contents.
    Ignores nested lists

    Parameters
    ----------
    input_list : list
        List containing mixtures of data types.

    Returns
    -------
    BiggestValue : list
        list of values in input list
    """

    tmp_arr = []
    if isinstance(input_list, list):
        for i in range(len(input_list)):
            try:
                input_list[i] + 0
                tmp_arr.append(input_list[i])
            except:
                pass

    else:  # assume integer
        tmp_arr = input_list

    if not isinstance(tmp_arr, int) and len(tmp_arr) == 1:
        tmp_arr = tmp_arr[0]
    return tmp_arr


def SeriesConstraints(input_list):
    """
    Returns a list of the constraints for the parameter.

    Parameters
    ----------
    input_list : list
        Order parmeter list or integer.

    Returns
    -------
    SeriesConstraints : list
        List of the constraints for the series.

    """
    SeriesConstraints = []
    if isinstance(input_list, list):
        for i in range(len(input_list)):
            try:
                input_list[i] + 0
            except:
                SeriesConstraints.append(input_list[i].lower())
    return SeriesConstraints
