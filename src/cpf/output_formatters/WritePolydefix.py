__all__ = ["Requirements", "WriteOutput"]

import json
import os
import glob
from pathlib import Path

import numpy as np

import cpf.output_formatters.WriteMultiFit as WriteMultiFit
from cpf.output_formatters.crystallographic_operations import plane_indices_4_to_3
from cpf.settings import get_settings
from cpf.IO_functions import make_outfile_name, peak_hkl
from cpf.output_formatters.ReadFits import ReadFits_to_list, ReadFits_to_dataframe
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WritePolydefix")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = { # dictionary
        # "time": "time_label", # time stamps for the diffraction paterns
        # "temperature": "*tc1_calcs.I", # which thermocouple to call.       
        "Phase": True,  # the phase we are interested in -- if True then guesses most common phase
        "ElasticProperties": True, # default is to use the phase name of the material. If more than 1 material need wild cards to match phase names. 
        "differential_only": False,
        # "which_thermocouple": 1,  # which thermocoule to include from 6BMB/X17B2 collection system. default to 1. #FIX ME: this needs to be included
        ###"Output_directory",  # if no direcrtory is specified write to current directory.
        # "Output_TemperaturePower",
        # "Output_tc",
        # "SampleGeometry": "3d",  # changes the strain tensor calucaltion from 2d to 3d. This determines how the cetroid and differnetial strain of the dpsaice are extracted from the fourier series.
        # "SampleDeformation": "compression",  # changes calculation between 'compression' and 'extension'.
    }
    # OptionalParams = [
    #     "Output_ElasticProperties",  # FIX ME: this needs to be included
    #     "tc",  # which thermocoule to include from 6BMB/X17B2 collection system. default to 1. #FIX ME: this needs to be included
    #     # "Output_directory",  # if no direcrtory is specified write to current directory. -- not set by default in settings class.
    #     "phase",  # the phase we are interested in
    #     "datafile_StartNum",  # start num and end num are needed for Polydefix which assumes continuous data files.
    #     "datafile_EndNum",
    #     "datafile_NumDigit",
    #     # FIXME: we should be able to re-number the data files so that processing them is possible with polydexif.
    #     # FIXME: these are addressed in the code as .
    # ]

    # append the requirements from WriteMultiFit to the lists because PolyDefix requires WriteMultiFit.
    r, o = WriteMultiFit.Requirements()
    RequiredParams = RequiredParams + r
    OptionalParams = OptionalParams | o

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, differential_only=False, **kwargs):
def WriteOutput(
    settings,
    differential_only=False,
    debug=False,
    **kwargs,
):
    """
    Writes *.exp files required by polydefix. 
    Calls WriteMltiFit to create *.fit files needed by polydefix
    
    Polydefix: Merkel and Hilairet (2015) http://dx.doi.org/10.1107/S1600576715010390.
    
    N.B. this is a different file than that required by polydefix for energy dispersive diffraction.
    
    Parameters
    ----------
    settings_class : cpf settings class
        Settings class used for fitting the data.
    differential_only : bool, optional
        Only include the differential part of the strain (cos^2 and sin^2 parts of the Fourier series). 
        Ignore the offset (cos and sin parts of the Fourier series).
        The default is False.
        
    Returns
    -------
    None.

    """
    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters
    Phase             = settings_class.output_settings.get("Phase", Requirements()[1]["Phase"])
    ElasticProperties = settings_class.output_settings.get("ElasticProperties", Requirements()[1]["ElasticProperties"])
    differential_only = settings_class.output_settings.get("differential_only", Requirements()[1]["differential_only"])
    #override with kwargs
    Phase             = kwargs.get("Phase", Phase)
    ElasticProperties = kwargs.get("ElasticProperties", ElasticProperties)
    differential_only = kwargs.get("differential_only", differential_only)
    
        
    # write *.fit files
    WriteMultiFit.WriteOutput(
        settings_class, differential_only=differential_only, debug=debug
    )

    # get the fits
    fits, _ = ReadFits_to_list(settings=settings_class)
    fitsDF = ReadFits_to_dataframe(settings=settings_class)
    
    #parse Phase and ElasticProperties
    if Phase is True:
       phases = fitsDF["phase"].unique()
       num_occurences = []
       for i in phases:
           num_occurences.append(fitsDF["phase"].str.count(i).sum())
       Phase = [phases[num_occurences.index(max(num_occurences))]] 
    elif isinstance(Phase, str):
        Phase = [Phase]
    else:
        #phase is a list
        pass
    
    
    base = settings_class.datafile_basename
    if base is None:
        logger.info(
            " ".join(map(str, [("No base filename, using input filename instead.")]))
        )
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    if differential_only is not False:
        base = base + "_DiffOnly"

    # write all the files
    for i in range(len(Phase)):
        fnam = base
        add_txt = None
        
        if len(Phase) > 1:
            add_txt = Phase[i]
        else:
            add_txt = None # because Files = 1
        out_file = make_outfile_name(
            fnam,
            directory=settings_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".exp",
            overwrite=True,
            additional_text=add_txt
        )
        text_file = open(out_file, "w")
        logger.info(" ".join(map(str, [("Writing %s" % out_file)])))

        # headers. set file version to be 1.
        text_file.write("# Experiment analysis file. to be used with Polydefix\n")
        text_file.write("# For more information: http://merkel.zoneo.net/Polydefix/\n")
        text_file.write(
            "# File Created by WritePolydefix function in ContinuousPeakFit\n"
        )
        text_file.write("# For more information: http://www.github.com/me/something\n")

        text_file.write("# File version\n")
        text_file.write("     2\n")

        # Write data properties
        text_file.write("# Directory with FIT files\n")
        # needs absolute path for the data files and the string has to end with a '/' e.g. Users/me/data/BCC1_2GPa_10s_e/
        text_file.write(
            "     %s/\n"
            % str(Path((os.getcwd())).resolve() / str(settings_class.output_directory))
        )
        text_file.write("# Basename for FIT files\n")
        # if last symbol in the file name is '_' then we need to strip it from the name
        text_file.write(
            "     %s\n" % settings_class.datafile_basename.strip("_").strip(".")
        )
        # if "datafile_startnum" in self.settings_from_input:
        #     self.datafile_startnum  = self.settings_from_input["datafile_StartNum"]
        #     self.datafile_endnum    = self.settings_from_input["datafile_EndNum"]
        #     self.datafile_numdigits = self.settings_from_input[datafile_NumDigit"]
        if (
            settings_class.settings_from_input["datafile_StartNum"]
            > settings_class.settings_from_input["datafile_EndNum"]
        ):
            logger.info(" ".join(map(str, [("start>end")])))
            strt = settings_class.settings_from_input["datafile_EndNum"]
            eend = settings_class.settings_from_input["datafile_StartNum"]
        else:
            strt = settings_class.settings_from_input["datafile_StartNum"]
            eend = settings_class.settings_from_input["datafile_EndNum"]
        text_file.write("# First index for FIT files\n")
        text_file.write("     %i\n" % strt)
        text_file.write("# Last index for FIT files\n")
        text_file.write("     %i\n" % eend)
        text_file.write("# Number of digits for FIT files\n")
        text_file.write(
            "     %i\n" % settings_class.settings_from_input["datafile_NumDigit"]
        )
        text_file.write("# Wavelength\n")
        # text_file.write("     %8.7g\n" % settings_class.data_class.calibration["conversion_constant"])
        text_file.write("     %8.7g\n" % settings_class.data_class.conversion_constant)
        text_file.write("# Fit offset for maximum stress 1 for yes, 0 for no\n")
        text_file.write("     %i\n" % 1)
        text_file.write("# Starting offset value, in degrees\n")
        text_file.write("     %8.4g\n" % 10.0000)
        text_file.write("# Fit beam center; 1 for yes, 0 for no\n")
        text_file.write("     %i\n" % 0)
        text_file.write("# Material properties set (1/0)\n")
        text_file.write("     %i\n" % 1)
        text_file.write("# Peaks properties set (1/0)\n")
        text_file.write("     %i\n" % 1)

        # write number of peaks and hkls.
        numpeaks = 0
        for x in range(len(settings_class.fit_orders)):
            numpeaks = numpeaks + len(settings_class.fit_orders[x]["peak"])
        text_file.write("# Number of peaks\n")
        text_file.write("     %i\n" % numpeaks)
        text_file.write("# Peaks info (use, h, k, l)\n")
        for x in range(len(settings_class.fit_orders)):
            for y in range(len(settings_class.fit_orders[x]["peak"])):
                
                # FIXME: use this line below as a shortening for all the x and y pointers
                settings_class.set_subpattern(i, x)
                
                use = 1
                if Phase[i] != settings_class.fit_orders[x]["peak"][y]["phase"]:
                    use = 0
                
                # check if the d-spacing fits are NaN or not. if NaN switch off.
                if type(fits[i][x]["peak"][y]["d-space"][0]) == type(None) or np.isnan(
                    fits[i][x]["peak"][y]["d-space"][0]
                ):
                    use = 0

                if "hkl" in settings_class.fit_orders[x]["peak"][y]:
                    hkl = str(settings_class.fit_orders[x]["peak"][y]["hkl"])
                    if hkl == "0" or hkl == 0:
                        hkl = "000"
                else:
                    hkl = "000"
                    use = 0
                    
                hkl = peak_hkl(settings_class.fit_orders[x], peak=y, string=False)[0]
                h,k,l = hkl
                
                text_file.write(" %5i    %s    %s    %s\n" % (use, h, k, l))

        # material properties
        text_file.write("# Material properties\n")
        
        if ElasticProperties:
            if ElasticProperties is True or ElasticProperties == Phase[i]:
                #use phase name to get elastic properties
                fname = glob.glob(f"*{Phase[i]}*")
            elif isinstance(ElasticProperties, str) and "*" in ElasticProperties:
                fname = glob.glob(ElasticProperties)
            else:
                fname = ElasticProperties
            if isinstance(fname, list):
                if len(fname) == 1:
                    fname= fname[0]
                else:
                    raise ValueError("More than 1 property file has been identified.")
            fid = open(fname, "r")     
            # pipe ealstic properties to the output file.
            text_file.write(fid.read())   
            
        elif not Phase:
            # left empty on purpose
            text_file.write("")

        else:
            # FIX ME: here I am just writing the elastic properties of BCC iron with no regard for the structure of the file or the avaliable data.
            #         This should really be fixed.
            text_file.write("# Name:  (N.B. this is the default phase)\n")
            text_file.write("BCC-Fe\n")
            text_file.write("# Symmetry\n")
            text_file.write("cubic  \n")
            text_file.write("# EOS stuff (v0, k0, k'0)\n")
            text_file.write("      23.5530      159.900      6.52000\n")
            text_file.write("# Elastic model\n")
            text_file.write("       1\n")
            text_file.write(
                "# Parameters for isotropic elastic model (k0 k1 k2, g0, g1, g2)\n"
            )
            text_file.write("      0.00000      0.00000      0.00000\n")
            text_file.write("      0.00000      0.00000      0.00000\n")
            text_file.write("# Parameters for anisotropic elastic model (Cij)\n")
            text_file.write(
                "      223.000      127.000      127.000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      127.000      223.000      127.000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      127.000      127.000      223.000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      122.000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      122.000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      122.000\n"
            )
            text_file.write("# Parameters for anisotropic elastic model (dCij/dp)\n")
            text_file.write(
                "      6.12245      4.08163      4.08163      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      4.08163      6.12245      4.08163      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      4.08163      4.08163      6.12245      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      2.44898      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      2.44898      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      2.44898\n"
            )
            text_file.write("# Parameters for anisotropic elastic model (d2Cij/dp2)\n")
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            )
            text_file.write(
                "     0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n"
            )

        text_file.close()
