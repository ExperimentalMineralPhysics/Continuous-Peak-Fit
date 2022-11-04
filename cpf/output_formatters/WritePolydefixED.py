__all__ = ["Requirements", "WriteOutput"]


import os

import numpy as np

# import cpf.PeakFunctions as ff
import cpf.series_functions as sf
import cpf.IO_functions as IO
import json
import re
import datetime


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        "ElasticProperties",  # FIX ME: this needs to be included
        "tc",  # which thermocoule to include from 6BMB/X17B2 collection system. default to 1. #FIX ME: this needs to be included
        ###"Output_directory",  # if no direcrtory is specified write to current directory.
        "phase",  # the phase we are interested in
        "Output_TemperaturePower",
        "Output_tc",
    ]

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(
    setting_class=None,
    setting_file=None,
    differential_only=False,
    debug=False,
    **kwargs
):
    # writes *.exp files required by polydefixED.
    # N.B. this is a different file than that required by polydefix for monochromatic diffraction.
    # This file contains a list of all the diffraction information. Hence it has to be written after the fitting as a single operation.

    if setting_class is None and setting_file is None:
        raise ValueError(
            "bummer Either the settings file or the setting class need to be specified."
        )
    elif setting_class is None:
        import cpf.XRD_FitPattern.initiate as initiate

        setting_class = initiate(setting_file)

    # FitParameters = dir(FitSettings)

    # Parse the required inputs.
    # base_file_name = FitSettings.datafile_Basename

    # get diffration file names and number.
    # diff_files, n_diff_files = IO.file_list(FitParameters, FitSettings)

    # base, ext = os.path.splitext(os.path.split(FitSettings.datafile_Basename)[1])
    # if not base:
    #     print("No base filename, using input filename instead.")
    #     base = os.path.splitext(os.path.split(FitSettings.inputfile)[1])[0]

    base = setting_class.datafile_basename
    if base is None:
        print("No base filename, using input filename instead.")
        base = os.path.splitext(os.path.split(setting_class.settings_file)[1])[0]
    if differential_only is not False:
        base = base + "_DiffOnly"
    # out_file = IO.make_outfile_name(
    #     base, directory=FitSettings.Output_directory, extension=".exp", overwrite=True
    # )
    out_file = IO.make_outfile_name(
        base,
        directory=setting_class.output_directory,  # directory=FitSettings.Output_directory,
        extension=".exp",
        overwrite=True,
    )

    # if 'Output_directory' in FitParameters:
    #     out_dir = FitSettings.Output_directory
    # else:
    #     out_dir = './'

    # # create output file name from passed name
    # path, filename = os.path.split(base_file_name)
    # base, ext = os.path.splitext(filename)
    # if base[-1:] == '_': #if the base file name ends in an '_' remove it.
    #     base = base[0:-1]
    # out_file = out_dir + base + '.exp'

    text_file = open(out_file, "w")
    print("Writing", out_file)

    # headers. set file version to be 1.
    text_file.write("# Experiment analysis file. to be used with PolydefixED\n")
    text_file.write("# For more information: http://merkel.zoneo.net/Polydefix/\n")
    text_file.write(
        "# File Created by WritePolydefixED function in ContinuousPeakFit\n"
    )
    text_file.write("# For more information: http://www.github.com/me/something\n")
    text_file.write("# File version\n")
    text_file.write("     1\n")

    # write two theta angle.
    # Writes two theta from the first detector.
    # FIX ME: should this be the average of all the two theat angles?
    text_file.write("# 2 theta angle (in degrees)\n")
    # text_file.write("%12.5f\n" % parms_dict["calibs"].mcas[0].calibration.two_theta)
    text_file.write(
        "%12.5f\n"
        % setting_class.data_class.calibration["calibs"].mcas[0].calibration.two_theta
    )

    # Write detector properties.
    text_file.write("# Number of detector positions\n")
    text_file.write("%6i\n" % len(setting_class.data_class.calibration["azimuths"]))
    text_file.write("# Angles for detector positions: Number. Use (1/0). Angle. \n")

    az_used = []
    for x in range(len(setting_class.data_class.calibration["azimuths"])):

        # determine if detector masked
        if setting_class.calibration_mask is not None:
            if x + 1 in setting_class.calibration_mask:
                use = 0
            else:
                use = 1
                az_used.append(setting_class.data_class.calibration["azimuths"][x])

        # print detector number, if it is being used, then angle, lastly just a number.
        text_file.write(
            "%6i  %i  %7.3f  %i \n"
            % (x + 1, use, setting_class.data_class.calibration["azimuths"][x], 1)
        )

        # FIX ME: if 10 element detector then turn detector 10 off.

    # Write number of peaks
    num_subpatterns = len(setting_class.fit_orders)
    num_peaks = 0
    for x in range(num_subpatterns):
        num_peaks = num_peaks + len(setting_class.fit_orders[x]["peak"])
    text_file.write("# Number of peaks\n")
    text_file.write("%6i \n" % num_peaks)

    # Write peak properties
    text_file.write("# Peak information \n")
    text_file.write("# Number. use (1/0). h. k. l\n")
    peak = 1
    for x in range(num_subpatterns):
        for y in range(len(setting_class.fit_orders[x]["peak"])):
            if "hkl" in setting_class.fit_orders[x]["peak"][y]:
                hkl = str(setting_class.fit_orders[x]["peak"][y]["hkl"])

                # check if the d-spacing fits are NaN or not. if NaN switch off (0) otherwise use (1).
                # Got first file name
                # filename = os.path.splitext(os.path.basename(diff_files[0]))[0]
                # filename = filename+'.json'

                filename = IO.make_outfile_name(
                    setting_class.datafile_list[0],
                    directory=setting_class.output_directory,
                    extension=".json",
                    overwrite=True,
                )
                # Read JSON data from file
                with open(filename) as json_data:
                    fit = json.load(json_data)
                # check if the d-spacing fits are NaN or not. if NaN switch off.
                if np.isnan(fit[x]["peak"][y]["d-space"][0]):
                    use = 0
                else:
                    use = 1
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
            text_file.write(" %5i %5i    %s    %s    %s\n" % (peak, use, h, k, l))
            peak = peak + 1

    # Write fetting behaviour
    text_file.write("# Fit offset for maximum stress 1 for yes. 0 for no\n")
    text_file.write("     1\n")
    text_file.write("# Starting offset value. in degrees\n")
    text_file.write("%6i \n" % setting_class.data_class.calibration["azimuths"][0])

    # Write material properties
    # FIX ME: need to be able to import material properties and write them to the file without messing about.
    text_file.write("# Material properties set (1/0)\n")
    text_file.write("     1\n")
    text_file.write("# Material properties\n")
    text_file.write("# Version\n")
    text_file.write("2       \n")

    if "Material" in setting_class.output_settings:

        raise ValueError(
            "''Material'' is depreciated as an option. Change your input file to have a ''Phase'' instead."
        )

    elif "phase" in setting_class.output_settings:
        # left empty on purpose
        text_file.write("")

    elif "ElasticProperties" in setting_class.output_settings:

        fid = open(setting_class.output_settings["ElasticProperties"], "r")

        # pipe ealstic properties to the output file.
        text_file.write(fid.read())

    else:

        # FIX ME: here I am just writing the elastic properties of olivine with no regard for the structure of the file or the avaliable data.
        #         This should really be fixed.
        text_file.write("# Name\n")
        text_file.write("Olivine\n")
        text_file.write("# Symmetry\n")
        text_file.write("ortho  \n")
        text_file.write("# EOS stuff (v0. k0. dk0/P. dK0/dT) (Li Li2006)\n")
        text_file.write("# Thermal expansion coefficients (a. b. c)\n")
        text_file.write("# 2.69000e-05 2.12000e-08 0.00000\n")
        text_file.write("# EOS stuff (v0. k0. dk0/P. dK0/dT)       \n")
        text_file.write(" 292.58 128.800 4.2000 0.00000         \n")
        text_file.write(
            "# Thermal expansion coefficients (a. b. c) (ISAAK ET AL.. 1989)   \n"
        )
        text_file.write(" 3.40e-05 0.7100e-08 0.00000          \n")
        text_file.write("# Elastic model           \n")
        text_file.write("1             \n")
        text_file.write(
            "# Parameters for isotropic elastic model (g0. g1. g2. g3. g4)   \n"
        )
        text_file.write(
            "# G = G0 + G1*P + G2*P*P + G3*(T-300) + G4*(T-300)*(T-300)  \n"
        )
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000        \n")
        text_file.write(
            "# Parameters for anisotropic elastic model (Cij) (ISAAK ET AL.. 1989)  \n"
        )
        text_file.write(" 323.700 66.4000 71.6000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 66.4000 197.600 67.2300 0.00000 0.00000 0.00000       \n")
        text_file.write(" 71.6000 67.2300 235.100 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 64.6200 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 78.0500 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 79.0400       \n")
        text_file.write(
            "# Parameters for anisotropic elastic model (dCij/dp) (ZHA. 1996)     \n"
        )
        text_file.write(" 7.98000 4.74000 4.48000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 4.74000 6.37000 3.76000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 4.48000 3.76000 6.38000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 2.17000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 1.64000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 2.31000       \n")
        text_file.write(
            "# Parameters for anisotropic elastic model (d2Cij/dp2)       \n"
        )
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(
            "# Parameters for anisotropic elastic model (dCij/dT) (ISAAK ET AL.. 1989)   \n"
        )
        text_file.write(
            " -0.0340000 -0.0105000 -0.00940000 0.00000 0.00000 0.00000       \n"
        )
        text_file.write(
            " -0.0105000 -0.0285000 -0.00510000 0.00000 0.00000 0.00000       \n"
        )
        text_file.write(
            " -0.00940000 -0.00510000 -0.0286000 0.00000 0.00000 0.00000       \n"
        )
        text_file.write(" 0.00000 0.00000 0.00000 -0.0128000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 -0.0130000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 -0.0157000       \n")
        text_file.write(
            "# Parameters for anisotropic elastic model (d2Cij/dT2)       \n"
        )
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000      \n")

    # Write details of experiment.

    # write number of steps
    text_file.write("# Number of time steps in the experiment\n")
    text_file.write("%i \n" % setting_class.datafile_number)
    text_file.write("# Information on time step       \n")
    text_file.write(
        "# Step number.  Step name.         Step time.    Step temperature (K).    taux de def E\n"
    )
    for x in range(setting_class.datafile_number):
        temp = np.array(
            0.0
        )  # pre-set temperature to 0 incase no value is found in the data file.
        with open(setting_class.datafile_list[x], "r") as myfile:
            for line in myfile:  # For each line, stored as line,

                # get date.
                if "DATE:" in line:
                    ftime = datetime.datetime.strptime(
                        line, "DATE:  %b %d, %Y %H:%M:%S.%f "
                    )
                    if x == 0:
                        t_o = ftime
                    t = ftime - t_o  # time relative to first diffraction

                # get temperature
                if (
                    "Output_TemperaturePower" in setting_class.output_settings
                ):  # get temperture from furnace power
                    st = "LVP_furnace_calcs.D"
                    if st in line:
                        pwr = float(re.findall(r'"(.*?)"', line)[0])
                        for i in range(
                            len(
                                setting_class.output_settings["Output_TemperaturePower"]
                            )
                        ):
                            temp = (
                                temp
                                + setting_class.output_settings[
                                    "Output_TemperaturePower"
                                ][i]
                                * pwr**i
                            )
                        temp = temp + 273
                else:  # get temperature from thermocouple reading
                    if not "Output_tc" in setting_class.output_settings:
                        st = "LVP_tc1_calcs.I"
                    else:
                        st = (
                            "LVP_tc"
                            + setting_class.output_settings["Output_tc"]
                            + "_calcs.I"
                        )
                    if st in line:
                        temp = re.findall(r'"(.*?)"', line)
                        temp = float(temp[0]) + 273.0

        # strip directory and ending from filename
        diff_name = os.path.splitext(os.path.basename(setting_class.datafile_list[x]))[
            0
        ]

        text_file.write(
            "%8i        %s  %8.8g            %6.1f\n"
            % (x + 1, "{:<18}".format(diff_name), t.seconds, temp)
        )

    # Write the experimental data fits.
    text_file.write("# Experimental data \n")
    text_file.write(
        "# Peak number.  h.  k.  l.  d-spacing.  intensity.  detector number. step number E t \n"
    )
    for z in range(setting_class.datafile_number):

        setting_class.set_subpattern(z, 0)

        filename = IO.make_outfile_name(
            setting_class.subfit_filename,
            directory=setting_class.output_directory,
            extension=".json",
            overwrite=True,  # overwrite = True to get the file name without incrlemeting it.
        )

        # Read JSON data from file
        with open(filename) as json_data:
            fit = json.load(json_data)

        peak = 1
        for x in range(num_subpatterns):
            for y in range(len(setting_class.fit_orders[x]["peak"])):

                # Write the following structure to the data file.
                # [peak number     H     K     L     d-spacing     Intensity     detector number     step number   ]  repeat for the number of
                # ...                                                                                ]  peaks*number detectors*number patterns
                # [peak number     H     K     L     d-spacing     Intensity     detector number     step number   ]  in data set.

                if "hkl" in setting_class.fit_orders[x]["peak"][y]:
                    hkl = str(setting_class.fit_orders[x]["peak"][y]["hkl"])
                    use = 1
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

                if "symmetry" in setting_class.fit_orders[x]["peak"][y]:
                    sym = setting_class.fit_orders[x]["peak"][y]["symmetry"]
                else:
                    sym = 1

                az = setting_class.data_class.calibration["azimuths"]
                coef_type = sf.params_get_type(fit[x], "d", peak=y)
                peak_d = sf.coefficient_expand(
                    np.array(az), fit[x]["peak"][y]["d-space"], coeff_type=coef_type
                )
                # peak_tth = 2.*np.degrees(np.arcsin(wavelength/2/peak_d))
                coef_type = sf.params_get_type(fit[x], "h", peak=y)
                peak_i = sf.coefficient_expand(
                    np.array(az_used) * sym,
                    fit[x]["peak"][y]["height"],
                    coeff_type=coef_type,
                )
                n = -1
                for w in range(
                    len(setting_class.data_class.calibration["calibs"].mcas)
                ):
                    # determine if detector masked
                    if setting_class.calibration_mask:
                        if w + 1 in setting_class.calibration_mask:
                            use = 0
                        else:
                            use = 1
                            n = n + 1
                    else:
                        use = 1
                        n = n + 1

                    if use == 1:  # write output for unmasked detectors
                        text_file.write(
                            "%8i        %s   %s   %s   %8.4f  %10.3f       %3i           %3i\n"
                            % (peak, h, k, l, peak_d[w], peak_i[n], w + 1, z + 1)
                        )
                    else:  # if the detetor is masked fill with zero heights
                        text_file.write(
                            "%8i        %s   %s   %s   %8.4f  %10.3f       %3i           %3i\n"
                            % (peak, h, k, l, peak_d[w], 0, w + 1, z + 1)
                        )

                peak = peak + 1

    text_file.close()
