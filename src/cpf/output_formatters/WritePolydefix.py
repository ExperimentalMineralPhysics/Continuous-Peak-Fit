__all__ = ["Requirements", "WriteOutput"]

import json
import os
from pathlib import Path

import numpy as np

import cpf.output_formatters.WriteMultiFit as WriteMultiFit
from cpf.output_formatters.crystallographic_operations import indicies4to3
from cpf.IO_functions import make_outfile_name
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WritePolydefix")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        "Output_ElasticProperties",  # FIX ME: this needs to be included
        "tc",  # which thermocoule to include from 6BMB/X17B2 collection system. default to 1. #FIX ME: this needs to be included
        # "Output_directory",  # if no direcrtory is specified write to current directory. -- not set by default in settings class.
        "phase",  # the phase we are interested in
        "datafile_StartNum",  # start num and end num are needed for Polydefix which assumes continuous data files.
        "datafile_EndNum",
        "datafile_NumDigit",
        # FIXME: we should be able to re-number the data files so that processing them is possible with polydexif.
        # FIXME: these are addressed in the code as .
    ]

    # append the requirements from WriteMultiFit to the lists because PolyDefix requires WriteMultiFit.
    r, o = WriteMultiFit.Requirements()
    RequiredParams = RequiredParams + r
    OptionalParams = OptionalParams + o

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, differential_only=False, **kwargs):
def WriteOutput(
    settings_class=None,
    settings_file=None,
    differential_only=False,
    debug=False,
    **kwargs,
):
    # writes *.fit files produced by multifit (as input for polydefix)
    # writes *.exp files required by polydefix.
    # N.B. this is a different file than that required by polydefix for energy dispersive diffraction.

    if settings_class is None and settings_file is None:
        raise ValueError(
            "bummer Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        from cpf.XRD_FitPattern import initiate

        settings_class = initiate(settings_file)

    # Write fit files
    # WriteMultiFit.WriteOutput(
    #     FitSettings, parms_dict, differential_only=differential_only
    # )
    WriteMultiFit.WriteOutput(
        settings_class=settings_class, differential_only=False, debug=debug
    )

    # FitParameters = dir(FitSettings)

    # Parse the required inputs.
    # base_file_name = FitSettings.datafile_Basename

    # diffraction patterns
    # diff_files, n_diff_files = IO.file_list(FitParameters, FitSettings)

    # base, ext = os.path.splitext(os.path.split(FitSettings.datafile_Basename)[1])
    base = settings_class.datafile_basename
    if base is None:
        logger.info(
            " ".join(map(str, [("No base filename, using input filename instead.")]))
        )
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    # if not base:
    #     logger.info(" ".join(map(str, [("No base filename, using input filename instead.")])))
    #     base = os.path.splitext(os.path.split(FitSettings.inputfile)[1])[0]
    if differential_only is not False:
        base = base + "_DiffOnly"
    # out_file = IO.make_outfile_name(base, directory=FitSettings.Output_directory, extension='.exp', overwrite=True)

    # #Parse the required inputs.
    # base_file_name = FitSettings.datafile_Basename
    # diff_files, n_diff_files = FileList(FitParameters, FitSettings)
    # if 'Output_directory' in FitParameters:
    #     out_dir = FitSettings.Output_directory
    # else:
    #     out_dir = './'

    # # create output file name from passed name
    # path, filename = os.path.split(base_file_name)
    # base, ext = os.path.splitext(filename)
    # if base[-1:] == '_': #if the base file name ends in an '_' remove it.
    #     base = base[0:-1]

    # if differential_only is not False:
    #     base = base+'_DiffOnly'
    # out_file = out_dir + base + '.exp'

    # get how many files to write -- one for each phase or set of elastic properties
    if "Output_ElasticProperties" in settings_class.output_settings:
        if isinstance(settings_class.output_settings["Output_ElasticProperties"], list):
            Files = len(settings_class.output_settings["Output_ElasticProperties"])
        else:
            Files = 1
            # force the elastic properties to be a list
            settings_class.output_settings["Output_ElasticProperties"] = [
                settings_class.output_settings["Output_ElasticProperties"]
            ]
    elif "phase" in settings_class.output_settings:
        if isinstance(settings_class.output_settings["phase"], list):
            Files = len(settings_class.output_settings["phase"])
        else:
            Files = 1
            # force the phase to be a list
            settings_class.output_settings["phase"] = [
                settings_class.output_settings["phase"]
            ]
    else:
        Files = 1

    # if "Output_ElasticProperties" in FitParameters:
    #     if isinstance(FitSettings.Output_ElasticProperties, list):
    #         Files = len(FitSettings.Output_ElasticProperties)
    #     else:
    #         Files = 1
    #         FitSettings.Output_ElasticProperties = [
    #             FitSettings.Output_ElasticProperties
    #         ]
    # elif "phase" in FitParameters:
    #     if isinstance(FitSettings.phase, list):
    #         Files = len(FitSettings.phase)
    #     else:
    #         Files = 1
    #         FitSettings.phase = [FitSettings.phase]
    # else:
    #     Files = 1

    for i in range(Files):
        fnam = base
        if Files > 1:
            # if "Output_ElasticProperties" in FitParameters:
            if "Output_ElasticProperties" in settings_class.output_settings:
                fnam = (
                    base
                    + "_"
                    + os.path.splitext(
                        os.path.split(
                            settings_class.output_settings["Output_ElasticProperties"][
                                i
                            ]
                        )[1]
                    )[0]
                )
            # elif "phase" in FitParameters:
            elif "phase" in settings_class.output_settings:
                fnam = (
                    base
                    + "_"
                    + os.path.splitext(
                        os.path.split(settings_class.output_settings["phase"][i])[1]
                    )[0]
                )

        out_file = make_outfile_name(
            fnam,
            directory=settings_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".exp",
            overwrite=True,
        )
        text_file = open(out_file, "w")
        logger.info(" ".join(map(str, [("Writing %s" % out_file)])))

        # headers. set file version to be 1.
        text_file.write("# Experiment analysis file. to be used with Polydefix\n")
        text_file.write("# For more information: http://merkel.zoneo.net/Polydefix/\n")
        text_file.write(
            "# File Created by WritePolydefixED function in ContinuousPeakFit\n"
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
        # if "datafile_startnum" in self.settings_from_file:
        #     self.datafile_startnum  = self.settings_from_file.datafile_StartNum
        #     self.datafile_endnum    = self.settings_from_file.datafile_EndNum
        #     self.datafile_numdigits = self.settings_from_file.datafile_NumDigit
        if (
            settings_class.settings_from_file.datafile_StartNum
            > settings_class.settings_from_file.datafile_EndNum
        ):
            logger.info(" ".join(map(str, [("start>end")])))
            strt = settings_class.settings_from_file.datafile_EndNum
            eend = settings_class.settings_from_file.datafile_StartNum
        else:
            strt = settings_class.settings_from_file.datafile_StartNum
            eend = settings_class.settings_from_file.datafile_EndNum
        text_file.write("# First index for FIT files\n")
        text_file.write("     %i\n" % strt)
        text_file.write("# Last index for FIT files\n")
        text_file.write("     %i\n" % eend)
        text_file.write("# Number of digits for FIT files\n")
        text_file.write(
            "     %i\n" % settings_class.settings_from_file.datafile_NumDigit
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

                if "hkl" in settings_class.fit_orders[x]["peak"][y]:
                    hkl = str(settings_class.fit_orders[x]["peak"][y]["hkl"])

                    if hkl == "0" or hkl == 0:
                        hkl = "000"

                    # check if the d-spacing fits are NaN or not. if NaN switch off (0) otherwise use (1).
                    # Got first file name
                    # filename = os.path.splitext(os.path.basename(diff_files[0]))[0]
                    # filename = filename+'.json'

                    filename = make_outfile_name(
                        settings_class.subfit_filename,  # diff_files[z],
                        directory=settings_class.output_directory,  # directory=FitSettings.Output_directory,
                        # diff_files[0],
                        # directory=FitSettings.Output_directory,
                        extension=".json",
                        overwrite=True,
                    )  # overwrite =false to get the file name without incrlemeting it.

                    # Read JSON data from file
                    with open(filename) as json_data:
                        fit = json.load(json_data)
                    # check if the d-spacing fits are NaN or not. if NaN switch off.
                    if type(fit[x]["peak"][y]["d-space"][0]) == type(None) or np.isnan(
                        fit[x]["peak"][y]["d-space"][0]
                    ):
                        use = 0
                    else:
                        use = 1

                    # check if the phase name matches the material (if exists)
                    # if "phase" in FitParameters:
                    if "phase" in settings_class.output_settings:
                        if (
                            # FitSettings.phase[i]
                            # == FitSettings.fit_orders[x]["peak"][y]["phase"]
                            settings_class.output_settings["phase"]
                            == settings_class.fit_orders[x]["peak"][y]["phase"]
                        ):
                            use = use
                        else:
                            use = 0
                    elif (
                        "Output_ElasticProperties" in settings_class.output_settings
                    ):  # elif "Output_ElasticProperties" in FitParameters:
                        # phase_str = FitSettings.Output_ElasticProperties[i]
                        # phase_str = settings_class.output_settings["Output_ElasticProperties"][i]
                        phase_str = os.path.splitext(
                            os.path.split(
                                settings_class.output_settings[
                                    "Output_ElasticProperties"
                                ][i]
                            )[1]
                        )[0]
                        strt = phase_str.find("_")
                        phase_str = phase_str[strt + 1 :]
                        if (
                            phase_str
                            == settings_class.fit_orders[x]["peak"][y]["phase"]
                        ):
                            use = use
                        else:
                            use = 0
                    else:
                        use = 0
                else:
                    hkl = "000"
                    use = 0
                pos = 0
                
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
                if len(hkl) > pos:
                    if hkl[pos] == "-":
                        m = hkl[pos : pos + 2]
                        pos = pos + 2
                    else:
                        m = hkl[pos : pos + 1]
                        pos = pos + 1
                    HKL = indicies4to3(np.array([h,k,l,m], dtype=int))
                    h, k, l = HKL
                text_file.write(" %5i    %s    %s    %s\n" % (use, h, k, l))

        # material properties
        text_file.write("# Material properties\n")
        if "Material" in settings_class.output_settings:
            raise ValueError(
                "''Material'' is depreciated as an option. Change your input file to have a ''Phase'' instead."
            )

        elif "Output_ElasticProperties" in settings_class.output_settings:
            fid = open(
                settings_class.output_settings["Output_ElasticProperties"][i], "r"
            )

            # pipe ealstic properties to the output file.
            text_file.write(fid.read())

        elif "phase" in settings_class.output_settings:
            try:
                fid = open(
                    (
                        settings_class.output_settings["phase"][i]
                        + "_elastic_properties.txt"
                    ),
                    "r",
                )

                # pipe ealstic properties to the output file.
                text_file.write(fid.read())
            except:
                try:
                    fid = open(
                        (
                            "Properties_"
                            + settings_class.output_settings["phase"][i]
                            + ".txt"
                        ),
                        "r",
                    )
                    # pipe ealstic properties to the output file.
                    text_file.write(fid.read())
                except:
                    raise ValueError(
                        "The files "
                        + settings_class.output_settings["phase"][i]
                        + "_elastic_properties.txt and Properties_"
                        + settings_class.output_settings["phase"][i]
                        + " are not found. Edit input file and try again."
                    )

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
