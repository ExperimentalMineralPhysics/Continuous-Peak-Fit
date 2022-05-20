#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["Settings"]

import json
import os
import sys
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import multiprocessing as mp
import importlib.util
import cpf.DioptasFunctions as DioptasFunctions
import cpf.GSASIIFunctions as GSASIIFunctions
import cpf.MedFunctions as MedFunctions
import cpf
from cpf.XRD_FitSubpattern import fit_sub_pattern
from cpf.Cosmics import cosmicsimage
from cpf.IO_functions import json_numpy_serializer, file_list, make_outfile_name, peak_string
from cpf.PeakFunctions import coefficient_type_as_number, get_number_coeff


def register_default_formats() -> object:
    """
    Load all available output modules
    :return:
    """
    # FIX ME: We could add extra checks here to make sure the required functions exist in each case.
    # FIX ME: Really these should live in their own sub-folder
    output_list = [
        "CoefficientTable",
        "DifferentialStrain",
        "MultiFit",
        "Polydefix",
        "PolydefixED",
    ]
    new_module = {}
    for output_module in output_list:
        module = __import__("cpf.Write" + output_module, fromlist=[None])
        new_module[output_module] = module
    return new_module


# Load potential output formats
output_methods_modules = register_default_formats()


class Settings:
    """Settings class definitions.
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
    
    def __init__(self, 
                 settings_file=None, 
                 out_type=None,
                 report=False):
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
    
        Notes
        -----
        Each required value is initialised as a blank instance. 
        These are then populated by the .populate(...) function if there is 
        a settings_file.
        
        Required settings are callable as direct functions. 
        Optional settings for the outputs are sored in a dictionary.
        
        """
        
        self.datafile_list = None
        self.datafile_number = 0
        self.datafile_directory = "."
        
        self.datafile_preprocess = None
        
        self.fit_bin_type = None
        self.fit_per_bin = None
        self.fit_number_bins = None
        self.fit_orders = None
        self.fit_bounds = {
              "background": ['min', 'max'],
              "d-space":    ['min', 'max'],
              "height":     [ 0,    '2*max'],
              "profile":    [ 0,     1],
              "width":      [ 'range/(ndata)',  'range/2'],
              }
        
        self.subfit_file_position = None
        self.subfit_filename = None
        self.subfit_order_position = None       
        self.subfit_orders = None
        
        self.output_directory = "."
        self.output_types = None
        self.output_settings = dict()
                
        self.calibration_type = None
        self.calibration_data = None
        self.calibration_mask = None
        self.calibration_parameters = None
        self.calibration_detector = None
        self.calibration_pixel_size = None
        
        # initiate the subpattern settings.
        # not set until called by 'set_subpatterns'. 
        self.subfit_file_position  = None     
        self.subfit_filename       = None
        self.subfit_order_position = None       
        self.subfit_orders         = None
        
        
        self.settings_file = settings_file
        #read the settings file given
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
     
    # def __getstate__(self):
    #     state = self
    #     # Don't pickle baz
    #     del builtins
    #     return state

    # def __setstate__(self, state):
    #     self.__dict__.update(state)
    #     # Add baz back since it doesn't exist in the pickle
    #     self.baz = 0
        
    
    def populate(self, 
        settings_file=None,
        inputs=None,
        out_type=None,
        report=False,
        ):
        """
        Fills the settings class from the settings file. 
        :param report:
        :param out_type:
        :param initiate_data:
        :param setting_file:
        :param inputs:
        :return fit_parameters:
        :return fit_settings:
        """
        
        # Fail gracefully
        if (
            settings_file is None
        ):
            raise ValueError(
                "Either the settings file or the parameter dictionary need to be specified."
            )
            
        
        if settings_file is not None:
            self.settings_file = settings_file
            
            if not str.endswith(self.settings_file, ".py"):
                self.settings_file = self.settings_file + ".py"
        
        
        self.check_files_exist(self.settings_file)
        self.read_settings_file()
        
        # override the files output settings.
        if not out_type is None:
            self.set_output_types(out_type_list = out_type)
            
            
    
    def reset(self):
        # set the values in the settings class back to those in the settings file. 
        self.populate()
        
                         
    def read_settings_file(self):
        #reads the settings file and adds values to the class.
        
        spec = importlib.util.spec_from_file_location(
                    "settings_module", self.settings_file
                )
        # store all the settings from file in the class. then sort them in a useful way. 
        self.settings_from_file = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.settings_from_file)
        
        all_settings_from_file = dir(self.settings_from_file)
        
        #add and check data directory        
        self.datafile_directory = self.settings_from_file.datafile_directory
        self.check_directory_exists(self.datafile_directory, write=True)
        
        #add and check data files
        self.datafile_list, self.datafile_number = file_list(dir(self.settings_from_file), self.settings_from_file)
        self.check_files_exist(self.datafile_list, write=False)
        print("All the data files exist.")
        
        # FIXME: datafile_base name should probably go because it is not a required variable it is only used in writing the outputs.
        if "datafile_Basename" in dir(self.settings_from_file):
            self.datafile_basename  = self.settings_from_file.datafile_Basename
        
        #add output directory if listed. 
        # change if listed among the inputs
        if "Output_directory" in dir(self.settings_from_file):
            self.output_directory = self.settings_from_file.Output_directory
            self.check_directory_exists(self.output_directory)
            
        # Load the detector class here to access relevant functions and check required parameters are present
        if "Calib_type" not in dir(self.settings_from_file):
            raise ValueError(
                "There is no 'Calib_type' in the settings. The fitting cannot proceed until a recognised "
                "calibration type is present."
            )
        else:
            self.calibration_type = self.settings_from_file.Calib_type
        if "Calib_param" not in dir(self.settings_from_file):
            raise ValueError(
                "There is no 'Calib_param' in the settings. The fitting cannot proceed until recognised "
                "calibration parameters are present."
            )
        else:
            self.calibration_parameters = self.settings_from_file.Calib_param
        if "Calib_data" in dir(self.settings_from_file):
            self.calibration_data = self.settings_from_file.Calib_data
        if "Calib_mask" in dir(self.settings_from_file):
            self.calibration_mask = self.settings_from_file.Calib_mask
        if "Calib_detector" in dir(self.settings_from_file):
            self.calibration_detector = self.settings_from_file.Calib_detector
        if "Calib_pixels" in dir(self.settings_from_file):
            self.calibration_pixel_size = self.settings_from_file.Calib_pixels
            
        # load the data class.
        self.data_class = detector_factory(self.calibration_type, self.calibration_parameters)
        
        
        if "Image_prepare" in dir(self.settings_from_file):
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
    

            
            
        self.fit_orders = self.settings_from_file.fit_orders
        self.validate_fit_orders()

        if "fit_bounds" in dir(self.settings_from_file):
            self.fit_bounds = self.settings_from_file.fit_bounds
            self.validate_fit_bounds()
        
        if "AziDataPerBin" in dir(self.settings_from_file):
            self.fit_per_bin = self.AziDataPerBin
            self.fit_bin_type = 0
        elif "AziBins" in dir(self.settings_from_file):
            self.fit_number_bins = self.settings_from_file.AziBins
            self.fit_bin_type = 1
        if "AziBinType" in dir(self.settings_from_file):
            self.fit_bin_type = self.AziBinType        
    
        if "Output_type" in dir(self.settings_from_file):
            #self.output_types = get_output_options(fit_settings.Output_type)
            self.set_output_types(out_type_list = self.settings_from_file.Output_type)
        
    
    def check_files_exist(self,files_to_check, write=False):
        # Check if a file exists. If not issue an error
        
        if isinstance(files_to_check, str):
            files_to_check = [files_to_check]    
        
        for j in range(len(files_to_check)):
            if os.path.isfile(files_to_check[j]) is False:
                raise ImportError(
                    "The file "
                    + files_to_check[j]
                    + " is not found but is required."
                )
            else:
                if write==True:
                    print(files_to_check[j] + " is avaliable.")
                    
        
    def check_directory_exists(self, directory, write=False):
        # Check if a directory exists. If not issue an error
        if os.path.exists(directory) is False:
            raise ImportError(
                "The directory " + directory + " is not found but is required."
            )
        else:
            if write==True:
                print(directory + " is avaliable.")

    
    
    def validate_fit_orders(self, report=False, peak=None, orders=None):
        
        if not self.fit_orders and not orders:
            raise ValueError("There are no fit orders.")
            
        # enable orders as the variable to check over -- then it is possible to validate orders that are not in the class. 
        if not orders:
            orders = self.fit_orders
            order_str = "fit_orders"
        else:
            order_str = "orders"
        # check the peak fitting options in the input file are not illicit.
        missing = []
        extras = []
    
        # allow just one peak-set to be validated.
        if peak == None:
            validate = list(range(len(orders)))
        else:
            validate = peak
    
        
        for i in validate:
            # FIX ME: we should check for incorrect or unneeded options
            required = ["background", "peak", "range"]
            possible = ["PeakPositionSelection", "Imax", "Imin"]
            comp_list = ["d-space", "width", "height", "profile"]
            comp_modifications= ["fixed", "type"]
    
    
            # check range 
            if "range" not in orders[i]:
                missing.append(order_str + " " + str(i) + " is missing a " "'range'")
                
            #check to see if range is list in list, e.g. [[16.0, 16.1]]. If so extract it to signle list.
            #this is old notation and now depreciated
            if (isinstance(orders[i]["range"], list) 
                    and len(orders[i]["range"]) == 1
                    and len(orders[i]["range"][0]) == 2):
                self.fit_orders[i]["range"] = self.fit_orders[i]["range"][0] 
            # check if  range is valid
            if (
                    not isinstance(orders[i]["range"], list)
                    or len(orders[i]["range"]) != 2
            ):
                missing.append(
                    order_str + "[" + str(i) + "] has an incorrectly formatted "+"'range'")
            
            if (orders[i]["range"][0] > orders[i]["range"][1]):
                missing.append(
                    order_str + "[" + str(i) + "]['range'] has values in wrong order")
            
            # check background
            if "background" not in orders[i]:
                missing.append(order_str + " " + str(i) + " is missing a "+"'background'")
            elif not isinstance(orders[i]["background"], list):
                missing.append(
                    order_str + " " + str(i) + " has an incorrectly formatted"+"'background'")
    
            # check peaks
            if "peak" not in orders[i]:
                missing.append(order_str + " " + str(i) + "has no "+"'peak'")
            else:
                for j in range(len(orders[i]["peak"])):
                    for k in range(len(comp_list)):
                        if not comp_list[k] in orders[i]["peak"][j]:
                            missing.append(
                                order_str + " "
                                + str(i)
                                + ", peak "
                                + str(j)
                                + " has no "
                                + comp_list[k]
                            )
                        elif not isinstance(
                            orders[i]["peak"][j][comp_list[k]], list
                        ) and not isinstance(
                            orders[i]["peak"][j][comp_list[k]], int
                        ):
                            missing.append(
                                order_str + " "
                                + str(i)
                                + ", peak "
                                + str(j)
                                + " incorrectly formatted "
                                "" + comp_list[k] + " "
                                " "
                            )
                            
                        #validate component type and length if fixed
                        if comp_list[k]+"_type" in self.fit_orders[i]["peak"][j]:
                            status = self.validate_order_type(self.fit_orders[i]["peak"][j][comp_list[k]+"_type"], self.fit_orders[i]["peak"][j][comp_list[k]])
                            if isinstance(status,str):
                                missing.append(status)
     
                        if "PeakPositionSelection" in self.fit_orders[i]:
                            missing.append(self.validate_position_selection(peak_set=i, report=False))
            
        #report missing bits and bobs
        if len(missing) > 0:
            for i in range(len(missing)):
                print(missing[i])
                
            if not report:
                raise ValueError("The problems listed above will prevent the data fitting.")
            else:
                print(
                    "The problems listed above will prevent the data fitting and need to be rectified before execution"
                )
        else:
            print("fit_orders appears to be correct")
                        
                            
                            
    def validate_order_type(comp_type):
        """
        Checks that a set component type is valid -- i.e. it exists in PeakFunctions.
        """
        
        rtrn = coefficient_type_as_number(comp_type, return_error=0)
        
        if isinstance(rtrn,str):
            status = comp_type + ": " + rtrn
        else:
            status = status
            
        return status
        
        
    # def validate_order_fixed(self, comp_fixed, comp_orders, comp_type=None):
    #     """
    #     Checks that a fixed component set is the same size of the orders that govern it. 
    #     FIXME: This should probably override one of the settings butI am not sure which.
    #     """
        
    #     # get component type as a number
    #     if isinstance(comp_type, str):
    #         coefficient_type_as_number(comp_type, return_error=1)
    #     elif isinstance(comp_type, None):
    #         comp_type = 0
            
            
    #     if 
        
        
    #     get_number_coeff(orders, comp, peak=0, azimuths=self.azm)
        
        
        
    #     comp_list = ["d-space", "width", "height", "profile"]
                
    #     for cp in range(len(comp_list)):
    #         comp = comp_list[cp]
            
    #         if comp_names[cp] + "_fixed" in self.fit_orders["peak"][peak]:
    #             fixed = 1
    #             vals = orders["peak"][j][
    #                 comp_names[cp] + "_fixed"
    #             ]  # values if fixed.
    #     if "profile_fixed" in orders["peak"][y]:
    #         if (
    #             not pf.fourier_order(orders["peak"][y]["profile_fixed"])
    #                 == orders["peak"][y]["profile"]
    #         ):
    #             print(
    #                 "Peak "
    #                 + str(y)
    #                 + ": The order of the profile does not match that of the fixed profile. "
    #                 "Changing profile to match fixed profile."
    #             )
    #             orders["peak"][y]["profile"] = pf.fourier_order(
    #                 orders["peak"][y]["profile_fixed"]
    #             )
    #         # make sure the fixed profile is a list
    #         if not isinstance(orders["peak"][y]["profile_fixed"], list):
    #             orders["peak"][y]["profile_fixed"] = [
    #                 orders["peak"][y]["profile_fixed"]
    #             ]
    # return orders, backgnd

    def validate_position_selection(self, peak_set=0, report=False):
    
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
                if np.min(tthguesses[:, 2]) < self.fit_orders[i]["range"][0][0]:
                    miss.append(
                        "fit_orders "
                        + str(i)
                        + ": PeakPositionSelection has at least one dispersion value (2theta or energy) that is too small"
                    )
                if np.max(tthguesses[:, 2]) > self.fit_orders[i]["range"][0][1]:
                    miss.append(
                        "fit_orders "
                        + str(i)
                        + ": PeakPositionSelection has at least one dispersion value (2theta or energy) that is too large"
                    )
            else: #there is no PeakPositionSelection. can only be 1 peak.
                if len(self.fit_orders[i]["peak"]) > 1:
                    miss.append(
                        "fit_orders " + str(i) + ": There are more peaks than listed in "
                        "PeakPositionSelection"
                        " than peaks listed"
                    )
        return miss
    
    
    def validate_fit_bounds(self, report=False):
        # check the peak fitting bounds in the input file are not illicit.
        if not self.fit_bounds:
            raise ValueError("There are no fit bounds.")
        
        required = ["background", "d-space", "height", "profile", "width"]

        # check all are present
        missing = []    
        for cp in range(len(required)):
            comp = required[cp]
            if comp not in self.fit_bounds:
                missing.append("fit_bounds is missing a " + comp)
            elif (
                not isinstance(self.fit_bounds[comp], list)
                and len(self.fit_bounds[comp]) != 2
                ):
                missing.append(
                    "fit_bounds has an incorrectly formatted " + comp)

        # list unrecognised entries
        # entries = self.fit_bounds.keys()
        extras = [x for x in self.fit_bounds.keys() if x not in required]
        
        #report findings
        incorrect=False
        if missing:
            print("\nMissing Values:")
            for i in range(len(missing)):
                print(missing[i])
            incorrect=True
        if extras:
            print("\nExtra Values:")
            for i in range(len(extras)):
                print(extras[i])
            incorrect=True
        if incorrect:
            if not report:
                raise ValueError("The problems listed above will prevent the data fitting.")
            else:
                print(
                    "The problems listed above will prevent the data fitting and need to be rectified before execution"
                )
        else:
            print("fit_bounds appears to be correct")
            
                     

        
    def set_output_types(self, out_type_list=None, report=False):
    
    
        # Load output function(s) and check output options are present and valid
        if out_type_list is not None:
            self.output_types = get_output_options(out_type_list)
            
        self.validate_output_types(self.output_types)
        

            
            
        
    def validate_output_types(self, report=False):
    
        # Check output format exists
        for mod in self.output_types:
            if mod not in output_methods_modules.keys():
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
            wr = output_methods_modules[self.output_types[i]]
            #print(self.output_types[i])
            required, optional = wr.Requirements()
            for j in range(len(required)):
                try:
                    self.output_settings[required[j]] = getattr(self.settings_from_file,required[j])
                except:
                    missing.append("The output " +self.output_types[i] + " requires the setting " + required[j])
            for j in range(len(optional)):
                try:
                    self.output_settings[optional[j]] = getattr(self.settings_from_file,optional[j])
                except:
                    missing.append("The output '" +self.output_types[i] + "' is missing the optional setting '" + optional[j] + "'") 
            
        if missing:
            print("\nMissing output settings:")
            for i in range(len(missing)):
                print(missing[i])
            print("The issues listed above may prevent outputs being written correctly")
        else:
            print("The output settings appear to be in order")  
                
                       
    def set_data_files(self, start=0, end=None, keep=None):
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
            end=keep+1
        elif end is None:
            end=len(self.datafile_list)
        
        self.datafile_list = self.datafile_list[start:end]
        self.datafile_number = len(self.datafile_list)   
        

    def set_subpatterns(self, subpatterns="all"):

        # restrict to sub-patterns listed
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


    def set_order_search(self, 
                        search_parameter="height", 
                        search_over=[0, 20],
                        subpatterns="all",
                        search_peak=0,
                        search_series=["fourier", "spline"]):
        
        
        # restrict to sub-patterns listed
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
                        orders_s["peak"][search_peak][
                            search_parameter + "-type"
                        ] = search_series[j]
                    else:
                        orders_s["background"][search_peak] = search[k]
                    if len(tmp_order) > 1:
                        intro_string = "peak="+str(search_peak)+"_"
                    else:
                        intro_string = ''
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
        Set the parameters for the subpattern to be fit as immetiately accesible. 
        It makes for shorter calls in XRD_Fit_Subpatten
        """
        
        self.subfit_file_position = number_subpattern   
        self.subfit_filename = self.datafile_list[file_number]
        self.subfit_order_position = number_subpattern        
        self.subfit_orders = self.fit_orders[number_subpattern]
    
def get_output_options(output_type):
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



def detector_factory(calibration_type, calibration_param, fit_settings=None):
    """
    Factory function to provide appropriate class for data dependent on type.
    Currently, supported options are Dioptas, GSASII and Med.
    :rtype: object
    :param fit_settings:
    :param calibration_type:
    :param calibration_param:
    :return:
    """
    if calibration_type == "Dioptas":
        detector_class = DioptasFunctions.DioptasDetector
        return detector_class(calibration_param, fit_settings)
    if calibration_type == "GSASII":
        detector_class = GSASIIFunctions.GSASIIDetector
        return detector_class(calibration_param, fit_settings)
    if calibration_type == "Med":
        detector_class = MedFunctions.MedDetector
        return detector_class(calibration_param, fit_settings)
    else:
        raise ValueError("Unrecognized calibration type.")