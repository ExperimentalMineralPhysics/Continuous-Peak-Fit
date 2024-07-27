#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import logging

# from cpf.XRD_FitPattern import logger

# =============================================================================
# Logger levels
# =============================================================================

# logging levels and what they record
# DEBUG   10   Detailed information, typically of interest only when diagnosing problems.
# INFO    20   Confirmation that things are working as expected.
#              Usually at this level the logging output is so low level that it’s not useful to users who are not familiar with the software’s internals.
# WARNING 30   An indication that something unexpected happened, or indicative of some problem in the near future (e.g. ‘disk space low’). The software is still working as expected.
# ERROR    40 	Due to a more serious problem, the software has not been able to perform some function.
# CRITICAL 50 	A serious error, indicating that the program itself may be unable to continue running.


# Configure the logging module
logging.basicConfig(
    level=logging.CRITICAL,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        # logging.FileHandler('CPF.log', mode='a', encoding=None, delay=False),  # Log to a file
        logging.StreamHandler()  # Log to stdout
    ],
)
# add custom logging levels.
levelNum = 17
levelName = "MOREINFO"


def logForMI(self, message, *args, **kwargs):
    if self.isEnabledFor(levelNum):
        self._log(levelNum, message, args, **kwargs)


def logToRoot(message, *args, **kwargs):
    logging.log(levelNum, message, *args, **kwargs)


logging.addLevelName(levelNum, levelName)
setattr(logging, levelName, levelNum)
setattr(logging.getLoggerClass(), levelName.lower(), logForMI)
setattr(logging, levelName.lower(), logToRoot)


levelName2 = "EFFUSIVE"
levelNum2 = 13


def logForE(self, message, *args, **kwargs):
    if self.isEnabledFor(levelNum2):
        self._log(levelNum2, message, args, **kwargs)


def logToRoot(message, *args, **kwargs):
    logging.log(levelNum2, message, *args, **kwargs)


logging.addLevelName(levelNum2, levelName2)
setattr(logging, levelName2, levelNum2)
setattr(logging.getLoggerClass(), levelName2.lower(), logForE)
setattr(logging, levelName2.lower(), logToRoot)

# Create a logger instance
logger = logging.getLogger(__name__)

"""
Logging levels and what they need to record
DEBUG       10  Detailed information, typically of interest only when diagnosing
                problems.
INFO        20  Confirmation that things are working as expected.
                Usually at this level the logging output is so low level that
                it’s not useful to users who are not familiar with the
                software’s internals.
WARNING     30  An indication that something unexpected happened, or indicative of
                some problem in the near future (e.g. ‘disk space low’). The
                software is still working as expected.
ERROR       40 	Due to a more serious problem, the software has not been able to
                perform some function.
CRITICAL    50 	A serious error, indicating that the program itself may be unable
                to continue running.

"""


def make_logger_output(level="DEBUG"):
    """
    Returns true/flase if logger is above or below 'level'

    Parameters
    ----------
    level : string, optional
        Logging module level label. The default is "DEBUG".

    Returns
    -------
    bool.

    """
    return logger.getEffectiveLevel() <= logger_level(level=level)


def logger_level(level="DEBUG"):
    """
    Get numberical value for logging level.

    Parameters
    ----------
    level : string, optional
        Debug level strings. The default is "DEBUG".

    Returns
    -------
    integer
        Logger level as integer.

    """
    return getattr(logging, level)


def log_pretty_print(
    lmtif_obj,
    oneline=False,
    colwidth=8,
    precision=4,
    fmt="g",
    columns=["value", "min", "max", "stderr", "vary", "expr", "brute_step"],
):
    """Pretty-print of parameters data.
    Edited after lmfit parmaeters pretty_print to return text as list.

    Parameters
    ----------
    oneline : bool, optional
        If True prints a one-line parameters representation [False]
    colwidth : int, optional
        Column width for all columns specified in `columns` [8]
    precision : int, optional
        Number of digits to be printed after floating point [4]
    fmt : {'g', 'e', 'f'}, optional
        Single-character numeric formatter. Valid values are: `'g'`
        floating point and exponential (default), `'e'` exponential,
        or `'f'` floating point.
    columns : :obj:`list` of :obj:`str`, optional
        List of :class:`Parameter` attribute names to print (default
        is to show all attributes).

    Returns
    --------
    output : list
        List of lines of text normally returned to stdout.

    """

    """
    pretty_print.
    https://github.com/lmfit/lmfit-py/blob/master/lmfit/parameter.py#L295

    License:
    BSD-3

    Copyright 2022 Matthew Newville, The University of Chicago
                   Renee Otten, Brandeis University
                   Till Stensitzki, Freie Universitat Berlin
                   A. R. J. Nelson, Australian Nuclear Science and Technology Organisation
                   Antonino Ingargiola, University of California, Los Angeles
                   Daniel B. Allen, Johns Hopkins University
                   Michal Rawlik, Eidgenossische Technische Hochschule, Zurich

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:

      1. Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

      2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

      3. Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from this
      software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
    LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    POSSIBILITY OF SUCH DAMAGE.
    """

    output = []
    if oneline:
        output.append(lmtif_obj.pretty_repr(oneline=oneline))
        return output

    name_len = max(len(s) for s in lmtif_obj)
    allcols = ["name"] + columns
    title = "{:{name_len}} " + len(columns) * " {:>{n}}"
    output.append(title.format(*allcols, name_len=name_len, n=colwidth).title())
    numstyle = "{%s:>{n}.{p}{f}}"  # format for numeric columns
    otherstyles = dict(
        name="{name:<{name_len}} ",
        stderr="{stderr!s:>{n}}",
        vary="{vary!s:>{n}}",
        expr="{expr!s:>{n}}",
        brute_step="{brute_step!s:>{n}}",
    )
    line = " ".join(otherstyles.get(k, numstyle % k) for k in allcols)
    for name, values in sorted(lmtif_obj.items()):
        pvalues = {k: getattr(values, k) for k in columns}
        pvalues["name"] = name
        # stderr is a special case: it is either numeric or None (i.e. str)
        if "stderr" in columns and pvalues["stderr"] is not None:
            pvalues["stderr"] = (numstyle % "").format(
                pvalues["stderr"], n=colwidth, p=precision, f=fmt
            )
        elif "brute_step" in columns and pvalues["brute_step"] is not None:
            pvalues["brute_step"] = (numstyle % "").format(
                pvalues["brute_step"], n=colwidth, p=precision, f=fmt
            )
        output.append(
            line.format(name_len=name_len, n=colwidth, p=precision, f=fmt, **pvalues)
        )
    output.append("")

    return output


def pretty_print_to_logger(lmtif_obj, space=True, level="DEBUG", **kwargs):
    """
    Wrapper for lmtif_obj.pretty_print() replacement that writes the
    pretty_print() output to the python logger.

    Parameters
    ----------
    lmtif_obj : lmfit.Parameters object
        lmfit object to be sent to the logger.
    space : bool, optional
        Create blank lines before and after the pretty_print(). The default is True.
    level : logging level string, optional
        Determines what logging level the output is written. The default is "DEBUG".
    **kwargs : TYPE
        Kwargs for pretty_print.

    """
    text = log_pretty_print(lmtif_obj, kwargs)
    if space == True:
        logger.log(getattr(logging, level.upper()), " ".join(map(str, [("")])))
    for i in range(len(text)):
        logger.log(getattr(logging, level.upper()), " ".join(map(str, [text[i]])))
        # equivalent to logger.LEVEL(msg)
        # getattr(logging, level) converts level to a number
    if space == True:
        logger.log(getattr(logging, level.upper()), " ".join(map(str, [("")])))
