from __future__ import annotations

import logging
import os
from typing import Literal, Optional, cast

import lmfit


class CPFLogger(logging.Logger):
    """
    This class contains the additional attributes that will be grafted onto whichever
    logger class the Python session will be run with.
    """

    """
    Define logging levels here
    """

    NOTSET = logging.NOTSET
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL
    EFFUSIVE = 13
    MOREINFO = 17

    levels = {
        "DEBUG": DEBUG,
        "INFO": INFO,
        "WARNING": WARNING,
        "ERROR": ERROR,
        "CRITICAL": CRITICAL,
        "EFFUSIVE": EFFUSIVE,
        "MOREINFO": MOREINFO,
    }

    def debug(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.DEBUG):
            self._log(self.DEBUG, msg, args, **kwargs)

    def info(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.INFO):
            self._log(self.INFO, msg, args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.WARNING):
            self._log(self.WARNING, msg, args, **kwargs)

    def error(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.ERROR):
            self._log(self.ERROR, msg, args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.CRITICAL):
            self._log(self.CRITICAL, msg, args, **kwargs)

    def moreinfo(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.MOREINFO):
            self._log(self.MOREINFO, msg, args, **kwargs)

    def effusive(self, msg, *args, **kwargs):
        self.check_level()
        if self.isEnabledFor(self.EFFUSIVE):
            self._log(self.EFFUSIVE, msg, args, **kwargs)

    """
    Helper functions and attributes
    """
    formatter: Optional[logging.Formatter] = None

    def set_formatter(
        self,
        format: Optional[str] = "%(asctime)s [%(levelname)s] %(message)s",
    ):
        """
        Set the formatter to be used by the logger.
        """
        if format is not None:
            self.formatter = logging.Formatter(format)

    def add_stream_handler(
        self,
    ):
        # Remove all existing StreamHandlers
        for handler in self.handlers[:]:
            if isinstance(handler, logging.StreamHandler):
                self.removeHandler(handler)
                # print(f"Removed stream handler {handler.name}")  # DEBUG MESSAGE
                handler.close()
        new_handler = logging.StreamHandler()
        new_handler.setFormatter(self.formatter)
        new_handler.setLevel(logging.NOTSET)  # Let the logger decide what to log
        self.addHandler(new_handler)

    def add_file_handler(
        self,
        file_name: str,
    ):
        # Remove all existing FileHandlers
        for handler in self.handlers[:]:  # Iterate over a copy to avoid mutation issues
            if isinstance(handler, logging.FileHandler):
                self.removeHandler(handler)
                # print(f"Removed file handler {handler.name}")  # DEBUG MESSAGE
                handler.close()
        # Create and attach a new FileHandler
        new_handler = logging.FileHandler(file_name, mode="a", encoding="utf-8")
        new_handler.setFormatter(self.formatter)
        new_handler.setLevel(logging.NOTSET)  # Let the logger decide what to log
        self.addHandler(new_handler)

    def check_level(self):
        """
        Checks to see if the current logger level has been changed from its previous one.
        If it has, sets the logger level to the new value. This is run before every logger
        call to ensure the logger is dynamically updated.
        """
        current_level = self.levels[os.environ["CPF_LOG_LEVEL"]]
        if self.level != current_level:
            self.level = current_level
            self.setLevel(self.level)
            # print(f"Logger {self.name} is now {self.level}")  # DEBUG MESSAGE

    def is_below_level(self, level: str = "DEBUG") -> bool:
        """
        Checks to see if the logger is operating above or below the target level,
        returning False and True respectively.
        """

        return self.getEffectiveLevel() <= self.levels[level]

    def log_lmfit_obj(
        self,
        lmfit_obj: lmfit.Parameters,
        space: bool = True,
        level: str = "DEBUG",
        **kwargs,
    ):
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

        def pretty_print_lmfit_obj(
            lmfit_obj: lmfit.Parameters,
            one_line: bool = False,
            col_width: int = 8,
            precision: int = 4,
            fmt: str = "g",
            columns: list[str] = [
                "value",
                "min",
                "max",
                "stderr",
                "vary",
                "expr",
                "brute_step",
            ],
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
            if one_line is True:
                output.append(lmfit_obj.pretty_repr(oneline=one_line))
                return output

            name_len = max(len(s) for s in lmfit_obj)
            allcols = ["name"] + columns
            title = "{:{name_len}} " + len(columns) * " {:>{n}}"
            output.append(
                title.format(*allcols, name_len=name_len, n=col_width).title()
            )
            numstyle = "{%s:>{n}.{p}{f}}"  # format for numeric columns
            otherstyles = dict(
                name="{name:<{name_len}} ",
                stderr="{stderr!s:>{n}}",
                vary="{vary!s:>{n}}",
                expr="{expr!s:>{n}}",
                brute_step="{brute_step!s:>{n}}",
            )
            line = " ".join(otherstyles.get(k, numstyle % k) for k in allcols)
            for name, values in sorted(lmfit_obj.items()):
                pvalues = {k: getattr(values, k) for k in columns}
                pvalues["name"] = name
                # stderr is a special case: it is either numeric or None (i.e. str)
                if "stderr" in columns and pvalues["stderr"] is not None:
                    pvalues["stderr"] = (numstyle % "").format(
                        pvalues["stderr"], n=col_width, p=precision, f=fmt
                    )
                elif "brute_step" in columns and pvalues["brute_step"] is not None:
                    pvalues["brute_step"] = (numstyle % "").format(
                        pvalues["brute_step"], n=col_width, p=precision, f=fmt
                    )
                output.append(
                    line.format(
                        name_len=name_len, n=col_width, p=precision, f=fmt, **pvalues
                    )
                )
            output.append("")

            return output

        """
        Start of 'log_lmfit_obj' function
        """

        text = pretty_print_lmfit_obj(lmfit_obj, **kwargs)

        if space is True:  # Add space before message
            self.log(getattr(CPFLogger, level.upper()), " ".join(map(str, [("")])))
        for i in range(len(text)):  # Print contents
            self.log(getattr(CPFLogger, level.upper()), " ".join(map(str, [text[i]])))
        if space is True:  # Add space after message
            self.log(getattr(CPFLogger, level.upper()), " ".join(map(str, [("")])))


def set_global_log_level(
    level: str,
):
    """
    This function sets an environment variable that stores the desired level of the
    logger in memory for as long as the current Python session is live, and can be
    updated in real-time when used as part of a Jupyter Notebook session. With this,
    it should be possible to switch log levels in between runs without having to
    restart the Python kernel.
    """
    level = level.upper()  # Standardise as uppercase
    if level not in CPFLogger.levels.keys():
        raise ValueError("Unrecognised logger level specified")
    os.environ["CPF_LOG_LEVEL"] = level
    print(f"Global logger level has been set to {level}")  # DEBUG


def configure_logger():
    """
    Comprehensively adds the attributes of CPFLogger class defined above to whatever
    logger object is currently being configured. When in a Python console or Jupyter
    Notebook, a Logger instance will be instantiated by default, so if another logger
    of our own is set up, we will get duplicate logs. By dynamically grafting the
    CPFLogger's attributes onto the existing Logger object, this avoids that issue.

    This method was inspired by the answers to Stack Overflow post
    http://stackoverflow.com/q/2183233/2988730, especially
    http://stackoverflow.com/a/13638084/2988730
    """

    # Set the global log level to INFO for now
    print("Performing initial loggger setup")
    set_global_log_level("INFO")

    # Add the additional attributes from the CPFLogger to the logger class being used
    # Overwrite what is already present
    for attribute, value in CPFLogger.__dict__.items():
        # Skip modification of built-in methods (they start with "__")
        if attribute.startswith("__"):
            continue
        setattr(logging.getLoggerClass(), attribute, value)
        # print(f"Added {attribute} to {logging.getLoggerClass()}")  # DEBUG message

    # Add the new levels from the CPFLogger into the main logger (overwrite existing)
    for levelName, levelNum in CPFLogger.levels.items():
        methodName = levelName.lower()

        def logToRoot(message, *args, **kwargs):
            logging.log(levelNum, message, *args, **kwargs)

        # Add additional log levels to the logging module
        logging.addLevelName(levelNum, levelName)
        setattr(logging, levelName, levelNum)
        setattr(logging, methodName, logToRoot)
        # print(f"Added logging level {levelName} to {logging.getLoggerClass()}")  # DEBUG

    # Load and configure the root logger used by this Python environment
    # It will have CPFLogger's attributes now
    logger = cast(CPFLogger, logging.getLogger())
    logger.set_formatter()
    logger.add_stream_handler()
    logger.setLevel(os.environ["CPF_LOG_LEVEL"])


# Run the function to configure the logging environment
configure_logger()


def get_logger(name: Optional[str] = None) -> CPFLogger:
    """
    Creates an instance of the logger class used in the current Python environment, with
    the logger taking on whatever name was specified. If no name was provided, the root
    logger is returned instead.

    In order to satisfy type checkers, this function re-casts the logger class used by
    the Python session as the CPFLogger.

    Normally, each module should have its own unique logger instance, but when working
    with Python consoles or Jupyter Notebooks, which spawn their own logger instances,
    this will lead to duplicate logs, so for now, the root logger is reused in each file.
    """
    logger = cast(CPFLogger, logging.getLogger())
    return logger
