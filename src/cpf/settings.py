#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations

__all__ = ["settings", "get_output_options", "detector_factory"]

import importlib.util
import json
from copy import copy, deepcopy
from pathlib import Path
from typing import Any, Literal, Optional

import numpy as np

import cpf.input_types as input_types
import cpf.output_formatters as output_formatters
from cpf.IO_functions import (
    image_list,
    json_numpy_serializer,
)

# , get_output_options, detector_factory, register_default_formats
from cpf.logging import CPFLogger
from cpf.series_functions import (
    coefficient_type_as_number,
    coefficient_type_as_string,
    coefficient_types,
    get_number_coeff,
)

# import logging

# logger = logging.getLogger(__name__)
logger = CPFLogger("cpf.settings")

# from cpf.XRD_FitPattern import logger


class Settings:
    """
    Settings class definitions.
    The settings class is contains all the variables/informtion needed to execute
    continuous peak fit.
    """

    """
    # Need to end up with:
        self.inputfile -- name of the input file.

        self.datafilelist -- list of file names to process

        self.datafile... ---

        self.orders -- dirctionary of the orders for fitting. needs to be ediatable. (e.g. add more, remove some)
        self.limits -- dictionary of the limits for the fitted functions.

        self.outputs -- list of output processes to run

    """
    print(logger.handlers)

    def __init__(
        self,
        settings_file: Optional[Path] = None,
        out_type=None,
        report: bool = False,
        debug: bool = False,
        mode: Literal["fit"] = "fit",
    ):
        """
        Initialise the cpf settings class.

        Parameters
        ----------
        option : settings_file
            *.py file containing the fit settings

        option : out_type
            Output type as list to override the settings in the file.

        option : report
            Not implemented.

        option : debug
            Verbose outputs to find errors.

        Notes
        -----
        Each required value is initialised as a blank instance.
        These are then populated by the .populate(...) function if there is
        a settings_file.

        Required settings are callable as direct functions.
        Optional settings for the outputs are sored in a dictionary.

        """

        self.datafile_list: list[str | Path] = []
        self.datafile_number: int = 0
        self.image_list: list[str | Path] = []
        self.image_number: int = 0
        self.datafile_directory = Path(".")

        self.datafile_preprocess = None

        self.file_label: Optional[str] = None

        # calibration type: dioptas etc.
        self.calibration_type: Literal["Dioptas", ""] = ""
        # file on which the calibration was done
        self.calibration_data: Optional[Path] = None
        # mask file for data
        self.calibration_mask: Optional[Path] = None
        # file with the calibration in it.
        self.calibration_parameters: Optional[Path] = None
        # FIXME: these are optional and should probalably be burried in an optional dictionary.
        self.calibration_detector: Literal["Pilatus1M", ""] = ""
        self.calibration_pixel_size = None

        self.fit_bin_type: Optional[int] = None
        self.fit_per_bin: Optional[int] = None
        self.fit_number_bins: Optional[int] = None
        self.fit_orders: list[dict[str, Any]] = []
        self.fit_bounds: dict[str, list[Any]] = {
            "background": ["0.95*min", "1.05*max"],
            "d-space": ["min", "max"],
            "height": [0, "1.05*max"],
            "profile": [0, 1],
            "width": ["range/(ndata)", "range/2"],
        }
        self.fit_min_data_intensity = 0
        self.fit_min_peak_intensity = "0.25*std"

        self.fit_track: bool = False
        self.fit_propagate: bool = True

        self.cascade_bin_type: Optional[int] = (
            0  # set default type - number data per bin
        )
        self.cascade_per_bin: Optional[int] = 50  # set default value
        self.cascade_number_bins: Optional[int] = 900  # set default value
        self.cascade_track: bool = False
        self.cascade_histogram_type: Literal["data"] = "data"
        self.cascade_histogram_bins: Optional[int] = None

        # output requirements
        self.output_directory: Optional[Path] = Path(".")
        self.output_types: list[str] = []
        self.output_settings: dict = {}

        # initiate the subpattern settings.
        # set to save diggging through self.fit_orders and carring values around
        # not set until called by 'set_subpatterns'.
        self.subfit_file_position = None
        self.subfit_filename = None
        self.subfit_order_position = None
        self.subfit_orders = None

        self.settings_file = settings_file
        # read the settings file given
        if self.settings_file is not None:
            self.populate(report=report)

        # #set the fitting defaults to carry around
        # self.refine=True,
        # self.save_all=False,
        # self.propagate=True,
        # self.iterations=1,
        # self.track=False,
        # self.parallel=True,
        # self.mode="fit",
        # self.report=False,

    def duplicate(self):
        """
        Makes a copy of an settings instance. But make deep copies of the sub_fit settings.
        This is for the parallel processing that needs an immutable object to work.
        """
        new = copy(self)
        new.subfit_file_position = deepcopy(self.subfit_file_position)
        new.subfit_filename = deepcopy(self.subfit_filename)
        new.subfit_order_position = deepcopy(self.subfit_order_position)
        new.subfit_orders = deepcopy(self.subfit_orders)
        return new

    def populate(
        self,
        settings_file: Optional[str | Path] = None,
        out_type=None,
        report=False,
        debug=False,
    ):
        """
        Fills the settings class from the settings file.

        Parameters
        ----------
        option : settings_file
            *.py file containing the fit settings

        option : out_type
            Output type as list to override the settings in the file.

        option : report
            Not implemented.

        option : debug
            Verbose outputs to find errors.

        """

        print("settings populate", logger.handlers)

        # Fail gracefully
        if settings_file is None:
            raise ValueError("The settings file needs to be specified.")

        # Convert to a Path object
        if isinstance(settings_file, str):
            try:
                settings_file = Path(settings_file)
            except Exception as error:
                raise error
        self.settings_file = settings_file
        if not self.settings_file.suffix == ".py":
            self.settings_file = self.settings_file.with_suffix(".py")

        self.check_files_exist(self.settings_file)
        self.read_settings_file()

        # override the files output settings.
        if not out_type is None:
            self.set_output_types(out_type_list=out_type)

    def reset(self):
        """
        set the values in the settings class back to those in the settings file.
        """
        self.populate()

    def read_settings_file(self):
        """
        adds values to the class from the settings file.
        Fails with a list of missing parameters if not complete.
        """

        # store all the settings from file in a mocule class.
        module_name = self.settings_file.stem
        spec = importlib.util.spec_from_file_location(module_name, self.settings_file)
        self.settings_from_file = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.settings_from_file)
        # then sort them in a useful way...

        ##all_settings_from_file = dir(self.settings_from_file)#

        # add data directory and data files
        self.datafile_directory = self.settings_from_file.datafile_directory
        if isinstance(self.datafile_directory, str):  # Convert to Path object
            self.datafile_directory = Path(self.datafile_directory)

        (
            self.datafile_list,
            self.datafile_number,
            self.image_list,
            self.image_number,
        ) = image_list(dir(self.settings_from_file), self.settings_from_file)
        # Convert datafile list entries to Path objects, if they exist
        if len(self.datafile_list) > 0:
            try:
                self.datafile_list = [
                    Path(file) for file in self.datafile_list if isinstance(file, str)
                ]
            except Exception as error:
                raise error

        # FIXME: datafile_base name should probably go because it is not a required variable it is only used in writing the outputs.
        if "datafile_Basename" in dir(self.settings_from_file):
            self.datafile_basename: str = self.settings_from_file.datafile_Basename
        if "datafile_Ending" in dir(self.settings_from_file):
            self.datafile_ending: str = self.settings_from_file.datafile_Ending

        # add output directory if listed.
        # change if listed among the inputs
        if "Output_directory" in dir(self.settings_from_file):
            self.output_directory = self.settings_from_file.Output_directory
            if isinstance(self.output_directory, str):
                self.output_directory = Path(self.output_directory)

        # Load the detector class here to access relevant functions and check required parameters are present
        if "Calib_type" in dir(self.settings_from_file):
            self.calibration_type = self.settings_from_file.Calib_type
        if "Calib_param" in dir(self.settings_from_file):
            self.calibration_parameters = self.settings_from_file.Calib_param
            if isinstance(self.calibration_parameters, str):
                self.calibration_parameters = Path(self.calibration_parameters)

        if "Calib_data" in dir(self.settings_from_file):
            self.calibration_data = self.settings_from_file.Calib_data
            if isinstance(self.calibration_data, str):
                self.calibration_data = Path(self.calibration_data)
        if "Calib_mask" in dir(self.settings_from_file):
            self.calibration_mask = self.settings_from_file.Calib_mask
            if isinstance(self.calibration_mask, str):
                self.calibration_mask = Path(self.calibration_mask)
        if "Calib_detector" in dir(self.settings_from_file):
            self.calibration_detector = self.settings_from_file.Calib_detector
        if "Calib_pixels" in dir(self.settings_from_file):
            self.calibration_pixel_size: int = self.settings_from_file.Calib_pixels

        # load the data class.
        self.data_class = detector_factory(fit_settings=self)

        if "Image_prepare" in dir(self.settings_from_file):
            logger.warning(
                "'Image_prepare' is depreciated nomenclature. Has been replased by 'image_preprocess'"
            )
            self.settings_from_file.image_preprocess = (
                self.settings_from_file.Image_prepare
            )

        if "image_preprocess" in dir(self.settings_from_file):
            self.datafile_preprocess = self.settings_from_file.Image_prepare

        #     # FIX ME: This doesn't seem to be used, if it should be this needs moving to class structure.
        #     alternatives_list = [[["datafile_StartNum", "datafile_EndNum"], ["datafile_Files"]]]
        #     optional_list = ["datafile_Step"]  # FIX ME: This doesn't seem to be used
        #     possible = [[[], []] * len(alternatives_list)]
        #     for x in range(len(alternatives_list)):
        #         for y in range(2):
        #             for z in range(len(alternatives_list[x][y])):
        #                 if alternatives_list[x][y][z] in fit_parameters:
        #                     possible[x][y].append(1)
        #                 else:
        #                     possible[x][y].append(0)
        #     # exit if all parameters are not present

        # organise the cascade properties
        if "cascade_number_bins" in dir(self.settings_from_file):
            self.cascade_number_bins = self.settings_from_file.cascade_number_bins
            self.cascade_bin_type = 1
        if "cascade_per_bin" in dir(self.settings_from_file):
            self.cascade_per_bin = self.settings_from_file.cascade_per_bin
            self.cascade_bin_type = 0
        if "cascade_bin_type" in dir(self.settings_from_file):
            self.cascade_bin_type = self.settings_from_file.cascade_bin_type
        if "cascade_historgram_type" in dir(self.settings_from_file):
            self.cascade_historgram_type = (
                self.settings_from_file.cascade_historgram_type
            )
        if "cascade_historgram_bins" in dir(self.settings_from_file):
            self.cascade_historgram_bins = (
                self.settings_from_file.cascade_historgram_bins
            )

        # organise the fits
        if "fit_orders" in dir(self.settings_from_file):
            self.fit_orders = self.settings_from_file.fit_orders
        if "fit_bounds" in dir(self.settings_from_file):
            self.fit_bounds = self.settings_from_file.fit_bounds
        if "fit_track" in dir(self.settings_from_file):
            self.fit_track = self.settings_from_file.fit_track
        if "fit_propagate" in dir(self.settings_from_file):
            self.fit_propagate = self.settings_from_file.fit_propagate
        if "fit_min_data_intensity" in dir(self.settings_from_file):
            self.fit_min_data_intensity = self.settings_from_file.fit_min_data_intensity
        if "fit_min_peak_intensity" in dir(self.settings_from_file):
            self.fit_min_peak_intensity = self.settings_from_file.fit_min_peak_intensity

        if "AziDataPerBin" in dir(self.settings_from_file):
            self.fit_per_bin = self.settings_from_file.AziDataPerBin
            self.fit_bin_type = 0
        elif "AziBins" in dir(self.settings_from_file):
            self.fit_number_bins = self.settings_from_file.AziBins
            self.fit_bin_type = 1
        if "AziBinType" in dir(self.settings_from_file):
            self.fit_bin_type = self.settings_from_file.AziBinType

        if "Output_type" in dir(self.settings_from_file):
            self.set_output_types(out_type_list=self.settings_from_file.Output_type)

        self.validate_settings_file()
        # FIXME: it needs to fail if everything is not present as needed and report what is missing

    def validate_settings_file(self):
        """
        Does all the validation of the settings file.
        Fails with missing parameters if not complete.
        """

        # check data files and directory
        if self.datafile_basename != None or self.datafile_directory != None:
            self.validate_datafiles()
        # FIX ME: check all the images exist.

        # check output directoy
        if self.output_directory != None:
            self.check_directory_exists(self.output_directory)

        # check calibration
        if self.calibration_type == None:
            raise ValueError(
                "There is no 'Calib_type' in the settings. The fitting cannot proceed until a recognised "
                "calibration type is present."
            )
        if self.calibration_parameters == None:
            raise ValueError(
                "There is no 'Calib_param' in the settings. The fitting cannot proceed until recognised "
                "calibration parameters are present."
            )

        # validate fit_orders and bounds
        if self.fit_orders:
            self.validate_fit_orders()
        if self.fit_bounds:
            self.validate_fit_bounds()

        # validate output types
        if self.output_types:
            self.validate_output_types()

    def check_files_exist(
        self,
        files_to_check: list[Path] | Path,
        write: bool = False,
    ):
        """
        Check if a file exists. If not issue an error
        """

        # Validate input
        if isinstance(files_to_check, (list, Path)) is False:
            raise TypeError("Input is of an unsupported type")

        if isinstance(files_to_check, Path):
            files_to_check = [files_to_check]

        for j in range(len(files_to_check)):
            # q = glob.iglob(files_to_check[j])
            # if not q:#glob.glob(files_to_check[j]):
            if not Path(".").glob(str(files_to_check[j])):
                # use glob.glob for a file search to account for compund detectors of ESRFlvp detectors
                raise ImportError(
                    f"The file {files_to_check[j]} is not found but is required."
                )
            else:
                if write == True:
                    logger.info(" ".join(map(str, [(f"{files_to_check[j]} exists.")])))

    def check_directory_exists(
        self,
        directory: Path,
        write: bool = False,
    ):
        """
        Check if a directory exists. If not issue an error
        """
        if directory.exists() is False:
            raise FileNotFoundError(
                f"The directory {directory.name!r} is not found but is required."
            )
        else:
            if write == True:
                logger.info(" ".join(map(str, [(f"{directory.name!r} exists.")])))

    def validate_datafiles(self):
        """
        Checks if the input directory and files all exist.
        """
        self.check_directory_exists(self.datafile_directory, write=True)
        self.check_files_exist(self.datafile_list, write=False)
        logger.info(" ".join(map(str, [("All the data files exist.")])))

    def validate_fit_orders(
        self,
        report: bool = False,
        peak=None,
        orders=None,
    ):
        """
        check that the orders of the fit contain all the needed parameters
        """

        if not self.fit_orders and orders is None:
            raise ValueError("There are no fit orders.")

        if not self.fit_orders:
            raise ValueError("No values for 'fit_orders' were provided")

        # enable orders as the variable to check over -- then it is possible to validate orders that are not in the class.
        if orders is None:
            orders = self.fit_orders
            order_str = "fit_orders"
        else:
            order_str = "orders"

        # allow just one peak-set to be validated.
        if peak == None:
            validate = list(range(len(orders)))
        else:
            validate = peak

        # check the peak fitting options in the input file are not illicit.
        missing: list = []
        extras: list = []

        for i in validate:
            # FIX ME: we should check for incorrect or unneeded options
            required = ["background", "peak", "range"]
            possible = ["PeakPositionSelection", "imax", "imin"]
            comp_list = ["d-space", "width", "height", "profile"]
            comp_modifications = ["fixed", "type"]

            # check range
            if "range" not in orders[i]:
                missing.append(order_str + " " + str(i) + " is missing a " "'range'")
            else:
                # check to see if range is list in list, e.g. [[16.0, 16.1]]. If so extract it to signle list.
                # this is old notation and now depreciated
                if (
                    isinstance(orders[i]["range"], list)
                    and len(orders[i]["range"]) == 1
                    and len(orders[i]["range"][0]) == 2
                ):
                    logger.warning(
                        " ".join(
                            map(
                                str,
                                [
                                    (
                                        "  input fix: subpattern "
                                        + str(i)
                                        + ": range is a now a simple list."
                                    )
                                ],
                            )
                        )
                    )
                    self.fit_orders[i]["range"] = self.fit_orders[i]["range"][0]
                # check if range is valid
                if (
                    not isinstance(orders[i]["range"], list)
                    or len(orders[i]["range"]) != 2
                ):
                    missing.append(
                        order_str
                        + "["
                        + str(i)
                        + "] has an incorrectly formatted "
                        + "'range'"
                    )
                if orders[i]["range"][0] >= orders[i]["range"][1]:
                    missing.append(
                        order_str
                        + "["
                        + str(i)
                        + "]['range'] has values in wrong order"
                    )

            # check background
            if "background" not in orders[i]:
                missing.append(
                    order_str + " " + str(i) + " is missing a " + "'background'"
                )
            elif not isinstance(orders[i]["background"], list):
                missing.append(
                    order_str
                    + " "
                    + str(i)
                    + " has an incorrectly formatted"
                    + "'background'"
                )
            else:
                # replace old component extansion naming convention.
                # bascially -- check if it is old nomlencature (*-*) and replace (with *_*)
                for l in range(len(comp_modifications)):
                    if "background" + "-" + comp_modifications[l] in self.fit_orders[i]:
                        logger.warning(
                            " ".join(
                                map(
                                    str,
                                    [
                                        (
                                            "  input fix: subpattern %s: 'background-%s' replaced with 'background_%s'"
                                            % (
                                                str(i),
                                                comp_modifications[l],
                                                comp_modifications[l],
                                            )
                                        )
                                    ],
                                )
                            )
                        )
                        self.fit_orders[i][
                            "background" + "_" + comp_modifications[l]
                        ] = self.fit_orders[i].pop(
                            "background" + "-" + comp_modifications[l]
                        )
                if "background_type" in self.fit_orders[i]:
                    status = self.validate_order_type(
                        self.fit_orders[i]["background_type"],
                        peak_set=i,
                        peak=0,
                        component="background",
                        report=report,
                    )
                    if isinstance(status, str):
                        missing.append(status)
                else:
                    status = 0

            # check peaks
            if "peak" not in orders[i]:
                missing.append(order_str + " " + str(i) + "has no " + "'peak'")
            else:
                for j in range(len(orders[i]["peak"])):
                    for k in range(len(comp_list)):
                        if not comp_list[k] in orders[i]["peak"][j]:
                            missing.append(
                                order_str
                                + " "
                                + str(i)
                                + ", peak "
                                + str(j)
                                + " has no "
                                + comp_list[k]
                            )
                        elif not isinstance(
                            orders[i]["peak"][j][comp_list[k]], list
                        ) and not isinstance(orders[i]["peak"][j][comp_list[k]], int):
                            missing.append(
                                order_str
                                + " "
                                + str(i)
                                + ", peak "
                                + str(j)
                                + " incorrectly formatted "
                                "" + comp_list[k] + " "
                                " "
                            )

                        # replace old component extansion naming convention.
                        # bascially -- check if it is old nomlencature (*-*) and replace (with *_*)
                        for l in range(len(comp_modifications)):
                            if (
                                comp_list[k] + "-" + comp_modifications[l]
                                in self.fit_orders[i]["peak"][j]
                            ):
                                logger.warning(
                                    " ".join(
                                        map(
                                            str,
                                            [
                                                (
                                                    "  input fix: subpattern %s %s: '%s-%s' replaced with '%s_%s'"
                                                    % (
                                                        str(i),
                                                        str(j),
                                                        comp_list[k],
                                                        comp_modifications[l],
                                                        comp_list[k],
                                                        comp_modifications[l],
                                                    )
                                                )
                                            ],
                                        )
                                    )
                                )
                                self.fit_orders[i]["peak"][j][
                                    comp_list[k] + "_" + comp_modifications[l]
                                ] = self.fit_orders[i]["peak"][j].pop(
                                    comp_list[k] + "-" + comp_modifications[l]
                                )

                        # validate component types
                        if comp_list[k] + "_type" in self.fit_orders[i]["peak"][j]:
                            status = self.validate_order_type(
                                self.fit_orders[i]["peak"][j][comp_list[k] + "_type"],
                                peak_set=i,
                                peak=j,
                                component=comp_list[k],
                                report=report,
                            )
                            if isinstance(status, str):
                                missing.append(status)
                        else:
                            status = 0
                        # validate fixed components
                        if comp_list[k] + "_fixed" in self.fit_orders[i]["peak"][j]:
                            if not isinstance(status, str):
                                mssng = self.validate_order_fixed(
                                    peak_set=i,
                                    peak=j,
                                    component=comp_list[k],
                                    report=report,
                                )
                                for m in range(len(mssng)):
                                    missing.append(mssng[m])

                        # check if phase and hkl are present. If neither then add them to make a label.
                        if (
                            not "hkl" in self.fit_orders[i]["peak"][j]
                            and not "phase" in self.fit_orders[i]["peak"][j]
                        ):
                            self.fit_orders[i]["peak"][j]["phase"] = "Region"
                            self.fit_orders[i]["peak"][j]["hkl"] = i + 1

                if "PeakPositionSelection" in self.fit_orders[i]:
                    mssng = self.validate_position_selection(peak_set=i, report=report)
                    for m in range(len(mssng)):
                        missing.append(mssng[m])

        missing = [x for x in missing if x != []]

        # report missing bits and bobs
        if len(missing) > 0:
            logger.error(
                " ".join(map(str, [("Input problems that will prevent execution: ")]))
            )
            for i in range(len(missing)):
                logger.error(" ".join(map(str, [(missing[i])])))
            logger.error(
                " ".join(
                    map(
                        str,
                        [
                            (
                                "The problems listed above will prevent the data fitting and need to be rectified before execution"
                            )
                        ],
                    )
                )
            )
            if not report:
                raise ValueError(
                    "The problems listed above will prevent the data fitting."
                )
        else:
            logger.info(" ".join(map(str, [("fit_orders appears to be correct")])))

    def validate_order_type(
        self,
        comp_type,
        peak_set,
        peak,
        component,
        report: bool = False,
    ):
        """
        Checks that a set component type is valid -- i.e. it exists in PeakFunctions.
        """
        # rtrn = coefficient_type_as_number(comp_type, return_error=0)

        if isinstance(comp_type, str):
            if comp_type in coefficient_types():
                status = 0
            else:
                status = (
                    "subpattern "
                    + str(peak_set)
                    + ", peak"
                    + str(peak)
                    + " "
                    + component
                    + ": the series type is not recognised"
                )
        elif isinstance(comp_type, int):
            try:
                coefficient_type_as_string(comp_type)
                status = 0
            except:
                status = (
                    "subpattern "
                    + str(peak_set)
                    + ", peak"
                    + str(peak)
                    + " "
                    + component
                    + ": the series type is not recognised"
                )

        return status

    def validate_order_fixed(
        self,
        peak_set,
        peak,
        component,
        report: bool = False,
    ):
        """
        Checks that a fixed component set is the same size of the orders that govern it.

        FIXME: This should possibly override one of the settings but I am not sure which.
        """

        if self.fit_orders is None:
            raise ValueError("No values were provided for 'fit_orders'")

        # get coefficient type and number of coefficients expected
        if component + "_type" in self.fit_orders[peak_set]["peak"][peak]:
            comp_type = coefficient_type_as_number(component, return_error=1)
        else:
            comp_type = 0

        num_coeffs = get_number_coeff(
            self.fit_orders[peak_set], component, peak=0, azimuths=self.data_class.azm
        )

        out = []

        if comp_type == 5:
            # indepdentent values for each azimuth!
            # this cannot be validated without referece to the data.
            # FIXME: use the data_class when loaded to check this.
            out.append(
                "    subpattern "
                + str(peak_set)
                + ", peak "
                + str(peak)
                + ": "
                + component
                + "_fixed could be not validated because order type is independent"
            )

        else:
            # make sure the fixed profile is a list
            # the only reason for it not to be a list is if it is a single value.
            if not isinstance(
                self.fit_orders[peak_set]["peak"][peak][component + "_fixed"], list
            ):
                self.fit_orders[peak_set]["peak"][peak][component + "_fixed"] = [
                    self.fit_orders[peak_set]["peak"][peak][component + "_fixed"]
                ]
                logger.warning(
                    " ".join(
                        map(
                            str,
                            [
                                (
                                    "    subpattern %s, peak %s: %s_fixed changed to a list"
                                    % (str(peak_set), str(peak), component)
                                )
                            ],
                        )
                    )
                )

            # validate
            if not num_coeffs == len(
                self.fit_orders[peak_set]["peak"][peak][component + "_fixed"]
            ):
                out.append(
                    "subpattern "
                    + str(peak_set)
                    + ", peak"
                    + str(peak)
                    + " "
                    + component
                    + "_fixed: The order does not match that of the fixed component. "
                )
        return out

    def validate_position_selection(
        self,
        peak_set: int | list[int] = 0,
        report: bool = False,
    ):
        """
        Checks that the multiple peak position selections have the right number of parts.
        """
        if self.fit_orders is None:
            raise ValueError("No values were provided for 'fit_orders'")

        if isinstance(peak_set, int):
            peak_set = [peak_set]

        miss = []
        for i in range(len(peak_set)):
            # if PeakPositionSelection - there might be more than one peak
            if "PeakPositionSelection" in self.fit_orders[i]:
                # how many peaks in list
                tthguesses = np.array(self.fit_orders[i]["PeakPositionSelection"])

                # if there are multiple peaks
                # FIX ME: use of j here outside of loop - need to check meaning!
                # ANSWER: it was there in error I think.
                if max(tthguesses[:, 0]) > len(self.fit_orders[i]["peak"]):
                    miss.append(
                        "fit_orders "
                        + str(i)
                        + ": PeakPositionSelection contains too many peaks"
                    )
                elif max(tthguesses[:, 0]) < len(self.fit_orders[i]["peak"]):
                    miss.append(
                        "fit_orders "
                        + str(i)
                        + ": PeakPositionSelection contains too few peaks"
                    )

                # if positions outside of range.
                if np.min(tthguesses[:, 2]) < self.fit_orders[i]["range"][0]:
                    miss.append(
                        "fit_orders "
                        + str(i)
                        + ": PeakPositionSelection has at least one dispersion value (2theta or energy) that is too small"
                    )
                if np.max(tthguesses[:, 2]) > self.fit_orders[i]["range"][1]:
                    miss.append(
                        "fit_orders "
                        + str(i)
                        + ": PeakPositionSelection has at least one dispersion value (2theta or energy) that is too large"
                    )

                for j in range(len(self.fit_orders[i]["peak"])):
                    # check that the number of selections is enough to make the expanded series from.
                    # FIXME: implement
                    pass

            else:  # there is no PeakPositionSelection. can only be 1 peak.
                if len(self.fit_orders[i]["peak"]) > 1:
                    miss.append(
                        "fit_orders " + str(i) + ": There are multiple peaks but no "
                        "PeakPositionSelection"
                        " listed"
                    )

        return miss

    def validate_fit_bounds(
        self,
        report: bool = False,
    ):
        """
        check the peak fitting bounds in the input file are valid.
        """

        if not self.fit_bounds:
            raise ValueError("There are no fit bounds.")

        required = ("background", "d-space", "height", "profile", "width")

        # check all are present
        missing = []
        for cp in range(len(required)):
            comp = required[cp]
            if comp not in self.fit_bounds:
                missing.append("fit_bounds is missing a " + comp)
            elif (
                not isinstance(self.fit_bounds[comp], list)
                and len(self.fit_bounds[comp]) != 2  # type: ignore  # Ignoring hard-to-fix type checking for now
            ):
                missing.append("fit_bounds has an incorrectly formatted " + comp)

        # list unrecognised entries
        # entries = self.fit_bounds.keys()
        extras = [x for x in self.fit_bounds.keys() if x not in required]

        # report findings
        incorrect = False
        if missing:
            logger.error(" ".join(map(str, [("")])))
            logger.error(" ".join(map(str, [("Missing Values:")])))
            for i in range(len(missing)):
                logger.error(" ".join(map(str, [(missing[i])])))
            incorrect = True
        if extras:
            logger.warning(" ".join(map(str, [("")])))
            logger.warning(" ".join(map(str, [("Extra Values:")])))
            for i in range(len(extras)):
                logger.warning(" ".join(map(str, [(extras[i])])))
            incorrect = True
        if incorrect:
            err_str = "The problems listed above will prevent the data fitting and need to be rectified before execution."
            if not report:
                logger.error(" ".join(map(str, [(err_str)])))
                raise ValueError(err_str)
            else:
                logger.error(" ".join(map(str, [(err_str)])))
        else:
            print("here")
            logger.info(" ".join(map(str, [("fit_bounds appears to be correct")])))

    def set_output_types(
        self,
        out_type_list: list[str] = [],
        report: bool = False,
    ):
        """
        set output types in settings class given a list of output types
        """
        if out_type_list:
            self.output_types = get_output_options(out_type_list)

    def validate_output_types(self, report=False):
        """
        validate output types
        """
        # Check output format exists
        for mod in self.output_types:
            if "Write" + mod not in output_formatters.module_list:
                raise ImportError(
                    "The 'Output_type' "
                    + mod
                    + " is not recognised; the file '"
                    + mod
                    + "' does not "
                    "exist. Check if "
                    "the output "
                    "type exists."
                )

        missing = []
        for i in range(len(self.output_types)):
            wr = getattr(output_formatters, "Write" + self.output_types[i])
            required, optional = wr.Requirements()
            for j in range(len(required)):
                try:
                    self.output_settings[required[j]] = getattr(
                        self.settings_from_file, required[j]
                    )
                except:
                    missing.append(
                        "The output "
                        + self.output_types[i]
                        + " requires the setting "
                        + required[j]
                    )
            for j in range(len(optional)):
                try:
                    self.output_settings[optional[j]] = getattr(
                        self.settings_from_file, optional[j]
                    )
                except:
                    missing.append(
                        "The output '"
                        + self.output_types[i]
                        + "' is missing the optional setting '"
                        + optional[j]
                        + "'"
                    )

        if missing:
            logger.warning(" ".join(map(str, [("Missing output settings:")])))
            for i in range(len(missing)):
                logger.warning(" ".join(map(str, [(missing[i])])))
            logger.warning(
                " ".join(
                    map(
                        str,
                        [
                            (
                                "The issues listed above may prevent outputs being written correctly"
                            )
                        ],
                    )
                )
            )
        else:
            logger.info(
                " ".join(map(str, [("The output settings appear to be in order")]))
            )

    def set_data_files(
        self,
        start=0,
        end=None,
        keep=None,
    ):
        """
        Cut the number of data files.

        Parameters
        ----------
        keep : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        if keep is not None:
            start = keep
            end = keep + 1
        elif end is None:
            end = len(self.datafile_list)

        # self.datafile_list = self.datafile_list[start:end]
        # self.datafile_number = len(self.datafile_list)
        self.image_list = self.image_list[start:end]
        self.image_number = len(self.image_list)

    def set_subpatterns(self, subpatterns="all"):
        """
        restrict to sub-patterns listed
        """
        if subpatterns == "all":
            sub_pats = list(range(0, len(self.fit_orders)))
        elif isinstance(subpatterns, list):
            sub_pats = subpatterns
        else:
            sub_pats = [int(x) for x in str(subpatterns)]
        # make new order search list
        orders_tmp = []
        for i in range(len(sub_pats)):
            j = sub_pats[i]
            orders_tmp.append(self.fit_orders[j])
        self.fit_orders = orders_tmp

    def set_order_search(
        self,
        search_parameter="height",
        search_over=[0, 20],
        subpatterns="all",
        search_peak=0,
        search_series=["fourier", "spline"],
    ):
        """
        set a range of orders to fit.
        This is used when determining what is best orders to use for fit.
        """
        if subpatterns == "all":
            subpatterns = list(range(0, len(self.fit_orders)))
        elif isinstance(subpatterns, list):
            subpatterns = subpatterns
        else:
            subpatterns = [int(x) for x in str(subpatterns)]

        # make new order search list
        if isinstance(search_over, list) and len(search_over) == 2:
            search = list(range(search_over[0], search_over[1]))
        else:
            search = [int(x) for x in str(search_over)]

        orders_search = []
        for i in range(len(subpatterns)):
            for j in range(len(search_series)):
                tmp_order = self.fit_orders[subpatterns[i]]
                for k in range(len(search)):
                    orders_s = deepcopy(tmp_order)
                    if search_parameter != "background":
                        orders_s["peak"][search_peak][search_parameter] = search[k]
                        orders_s["peak"][search_peak][search_parameter + "_type"] = (
                            search_series[j]
                        )
                    else:
                        orders_s["background"][search_peak] = search[k]
                    if len(tmp_order) > 1:
                        intro_string = "peak=" + str(search_peak) + "_"
                    else:
                        intro_string = ""
                    orders_s["note"] = (
                        search_parameter
                        + "="
                        + str(search[k])
                        + "_type="
                        + search_series[j]
                    )
                    orders_search.append(orders_s)
        self.fit_orders = orders_search

    def set_subpattern(self, file_number, number_subpattern):
        """
        Set the parameters for the subpattern to be fit as immediately accesible.
        It makes for shorter calls in XRD_Fit_Subpatten
        """

        self.subfit_file_position = number_subpattern
        # self.subfit_filename = self.datafile_list[file_number]
        self.subfit_filename = self.image_list[file_number]
        self.subfit_order_position = number_subpattern
        self.subfit_orders = self.fit_orders[number_subpattern]

    def save_settings(
        self, filename: str = "settings.json", filepath: Path = Path(".")
    ):
        """
        Saves the settings class to file.

        Parameters
        ----------
        filename : String
            filename to save the file as.
        filepath : string
            Filepath to save the file in.

        Returns
        -------
        None.

        """
        logger.warning(
            " ".join(
                map(
                    str,
                    [
                        (
                            "Caution: save_settings writes a temporary file with no content"
                        )
                    ],
                )
            )
        )

        fnam = filepath / filename
        with open(fnam, "w") as TempFile:
            # Write a JSON string into the file.
            json.dump(
                "This is a temporary file with no content",
                TempFile,
                sort_get_image_key_strings=False,
                indent=2,
                default=json_numpy_serializer,
            )
        logger.info(" ".join(map(str, [("Done writing", filename)])))


def get_output_options(output_type: list[str]):
    """
    Check if input is string or list of strings
    :param output_type: string or list of strings
    :return: list of strings
    """
    output_mod_type = []
    if isinstance(output_type, str):
        output_mod_type.append(output_type)
    else:
        output_mod_type = output_type
    return output_mod_type


def detector_factory(fit_settings: Settings):
    """
    Factory function to provide appropriate class for data dependent on type.
    *should* support any option that is named *Functions and contains *Detector as class.
    :rtype: object
    :param fit_settings:
    :param calibration_type:
    :param calibration_param:
    :return:
    """

    def_func = fit_settings.calibration_type + "Functions"
    def_def = fit_settings.calibration_type + "Detector"

    if def_func in input_types.module_list:
        detector = getattr(input_types, def_func)
        detector_class = getattr(detector, def_def)
        return detector_class(settings_class=fit_settings)
    else:
        raise ValueError("Unrecognized calibration type.")


if __name__ == "__main__":
    Settings
