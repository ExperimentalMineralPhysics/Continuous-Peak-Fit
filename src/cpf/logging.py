from __future__ import annotations

import logging
from typing import Optional

import lmfit


class CPFLogger(logging.Logger):
    def __init__(
        self,
        name: str,
        level: str | int = "INFO",
        format: Optional[str] = "%(asctime)s [%(levelname)s] %(message)s",
        log_dir: Optional[str] = None,
    ):
        # Inherit the properties of the Logger class this is derived from
        super().__init__(name, level)

        # Clear all pre-existing handlers
        self.handlers.clear()

        # Set the default logging level
        self.setLevel(level)

        # Set the format to use
        formatter = logging.Formatter(format) if isinstance(format, str) else None

        # Set up the stream handler (logs to the console)
        sh = logging.StreamHandler()
        sh.setFormatter(formatter)
        sh.setLevel(level)
        self.addHandler(sh)

        # Set up file handler if a file path is given
        if log_dir is not None:
            fh = logging.FileHandler(
                filename=log_dir,
                encoding="utf-8",
                mode="a",  # Append to existing log file
            )
            fh.setFormatter(formatter)
            fh.setLevel(level)
            self.addHandler(fh)

    """
    Define logging levels here
    """

    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL
    EFFUSIVE = 13
    MOREINFO = 17

    def debug(self, msg, *args, **kwargs):
        levelNum = self.DEBUG
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    def info(self, msg, *args, **kwargs):
        levelNum = self.INFO
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    def warning(self, msg, *args, **kwargs):
        levelNum = self.WARNING
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    def error(self, msg, *args, **kwargs):
        levelNum = self.ERROR
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    def critical(self, msg, *args, **kwargs):
        levelNum = self.CRITICAL
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    def moreinfo(self, msg, *args, **kwargs):
        levelNum = self.MOREINFO
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    def effusive(self, msg, *args, **kwargs):
        levelNum = self.EFFUSIVE
        if self.isEnabledFor(levelNum):
            self._log(levelNum, msg, args, **kwargs)

    """
    Helper functions
    """

    def is_below_level(self, level: str = "DEBUG") -> bool:
        """
        Checks to see if the logger is operating above or below the target level,
        returning False and True respectively.
        """

        return self.getEffectiveLevel() <= getattr(self, level)

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


# Basic test to check that logger properties are correct
def run_test(level: str | int = "INFO"):
    # Assign the logger to a variable
    logger = CPFLogger("CPFLogger")
    logger.setLevel(level)
    print(f"Logger currently set to level {logger.getEffectiveLevel()}")

    if len(logger.handlers) > 0:
        print("This logger instance has the following handlers:")
    for h in logger.handlers:
        print(h)

    # Print it at different levels
    logger.debug("Printing at the DEBUG level")
    logger.moreinfo("Printing at the MOREINFO level")
    logger.effusive("Printing at the EFFUSIVE level")
    logger.info("Printing at the INFO level")
    logger.warning("Printing at the WARNING level")
    logger.error("Printing at the ERROR level")
    logger.critical("Printing at the CRITICAL level")

    print(f"Current effective level: {logger.getEffectiveLevel()}")
    print(f"DEBUG attribute: {getattr(CPFLogger, 'DEBUG')}")
    print(f"EFFUSIVE attribute: {getattr(CPFLogger, 'EFFUSIVE')}")
    print(f"MOREINFO attribute: {getattr(CPFLogger, 'MOREINFO')}")
    print(f"INFO attribute: {getattr(CPFLogger, 'INFO')}")
    print(f"WARNING attribute: {getattr(CPFLogger, 'WARNING')}")
    print(f"ERROR attribute: {getattr(CPFLogger, 'ERROR')}")
    print(f"CRITICAL attribute: {getattr(CPFLogger, 'CRITICAL')}")


if __name__ == "__main__":
    """
    A simple way of running the test
    """
    run_test(level="INFO")
