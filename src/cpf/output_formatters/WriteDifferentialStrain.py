__all__ = ["Requirements", "WriteOutput"]


# import cpf.PeakFunctions as pf
import json
import os

import numpy as np
from lmfit.model import load_modelresult

import cpf.lmfit_model as lmm
import cpf.output_formatters.convert_fit_to_crystallographic as cfc
from cpf.IO_functions import (
    lmfit_fix_int_data_type,
    make_outfile_name,
    peak_string,
    replace_null_terms,
)
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteDifferentialStrain")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        "SampleGeometry"  # changes the strain tensor calucaltion from 2d to 3d. This determines how the cetroid and differnetial strain of the dpsaice are extracted from the fourier series.
        "SampleDeformation"  # changes calculation between 'compression' and 'extension'.
    ]

    return RequiredParams, OptionalParams


def WriteOutput(settings_class=None, settings_file=None, debug=True, **kwargs):
    """

    writes some of the fitted coeficients to a table. With a focus on the differential strain coefficents
    Parameters
    ----------
    settings_class : TYPE
        DESCRIPTION.
    parms_dict : TYPE
        DESCRIPTION.
    debug : TYPE, optional
        DESCRIPTION. The default is True.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # version 1.1 has elasped time and reduced chi squared added to the output table.
    # version 1.2 has corrected the orientation calculation. The reported angle is now that of the maximum compression
    # version 2  accounts for 2d and 3d experimntal gemoetries, compression and extension. The crystallographic values are now calculated in a different file/function.
    f_version = 2

    if settings_class is None and settings_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        from cpf.XRD_FitPattern import initiate

        settings_class = initiate(settings_file)

    # Parse optional parameters
    SampleGeometry = "3D".lower()
    if "SampleGeometry" in settings_class.output_settings:
        SampleGeometry = settings_class.output_settings["SampleGeometry"].lower()
    SampleDeformation = "compression".lower()
    if "SampleDeformation" in settings_class.output_settings:
        SampleDeformation = settings_class.output_settings["SampleDeformation"].lower()

    base = settings_class.datafile_basename

    # if not base:
    if base is None or len(base) == 0:
        logger.info(
            " ".join(("No base filename, trying ending without extension instead."))
        )
        base = settings_class.datafile_ending

    if base is None:
        logger.info(" ".join(("No base filename, using input filename instead.")))
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    out_file = make_outfile_name(
        base,
        directory=settings_class.output_directory,
        extension=".dat",
        overwrite=True,
        additional_text=settings_class.file_label,
    )

    text_file = open(out_file, "w")
    logger.info(" ".join(["Writing %s" % out_file]))
    text_file.write(
        "# Summary of fits produced by continuous_peak_fit for input file: %s.\n"
        % settings_class.settings_file  # FitSettings.inputfile
    )
    text_file.write(
        f"# Sample Geometry = {SampleGeometry}; Sample Deformation = {SampleDeformation}\n"
    )
    text_file.write("# \n")
    text_file.write(
        "# For more information: https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit\n"
    )
    text_file.write("# File version: %i \n" % f_version)
    text_file.write("# \n")

    # write header
    width_col = 12
    dp = 5  # used later or the number of decial places in the numbers -- but set here so all settings are in the same place.
    width_fnam = 25
    width_hkl = 15
    text_file.write(("# {0:<" + str(width_fnam - 2) + "}").format("Data File" + ","))
    text_file.write(("{0:<" + str(width_hkl) + "}").format("Peak" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d_mean" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d_mean_err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d2cos" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d2cos_err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d2sin" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d2sin_err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("corr coef" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("diff strain" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("diff s err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("orientation" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("orient err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d_max" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("d_min" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("mean h" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("h_err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("mean w" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("w_err" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("mean p" + ","))
    text_file.write(("{0:>" + str(width_col) + "}").format("p0_err" + ","))
    if debug:
        text_file.write(("{0:>" + str(width_col) + "}").format("Time taken" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Chunk time" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Sum Resid^2" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Status" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Func eval" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Num vars" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Num data" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Deg Freedom" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("ChiSq" + ","))
        text_file.write(("{0:>" + str(width_col) + "}").format("Red. ChiSq" + ","))
        text_file.write(
            ("{0:<" + str(width_col) + "}").format("Akaike Information Criterion" + ",")
        )
        text_file.write(
            ("{0:<" + str(width_col) + "}").format(
                "Bayesian Information Criterion" + ","
            )
        )
    text_file.write("\n")

    all_fits = []
    for z in range(settings_class.image_number):
        settings_class.set_subpattern(z, 0)

        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            extension=".json",
            overwrite=True,
        )

        if os.path.isfile(filename):
            # Read JSON data from file
            with open(filename) as json_data:
                fit = json.load(json_data)

            fit = replace_null_terms(fit)

            all_fits.append(fit)

            # calculate the required parameters.
            num_subpatterns = len(fit)
            for y in range(num_subpatterns):
                out_name = make_outfile_name(
                    settings_class.subfit_filename,
                    directory="",
                    overwrite=True,
                )
                # logger.info(" ".join(('  Incorporating ' + subfilename)))
                logger.info(
                    " ".join(
                        ["  Incorporating: %s,%s" % (out_name, peak_string(fit[y]))]
                    )
                )
                # try reading an lmfit object file.
                savfilename = make_outfile_name(
                    settings_class.subfit_filename,
                    directory=settings_class.output_directory,
                    orders=fit[y],
                    extension=".sav",
                    overwrite=True,
                )
                if os.path.isfile(savfilename):
                    try:
                        gmodel = load_modelresult(
                            savfilename, funcdefs={"peaks_model": lmm.peaks_model}
                        )
                    except:
                        lmfit_fix_int_data_type(savfilename)
                        try:
                            gmodel = load_modelresult(
                                savfilename, funcdefs={"peaks_model": lmm.peaks_model}
                            )
                        except:
                            # raise FileNotFoundError
                            logger.info(
                                " ".join(
                                    map(
                                        str,
                                        [
                                            (
                                                "    Can't open file with the correlation coefficients in..."
                                            )
                                        ],
                                    )
                                )
                            )
                    # FIX ME: this will only work for one peak. Needs fixing if more then one peak in subpattern
                    try:
                        corr = gmodel.params["peak_0_d3"].correl["peak_0_d4"]
                    except:
                        corr = np.nan

                elif (
                    "correlation_coeffs" in fit[y]
                ):  # try reading correlation coefficient from fit file.
                    try:
                        corr_all = json.loads(fit[0]["correlation_coeffs"])
                        corr = corr_all["peak_0_d3"]["peak_0_d4"]
                    except:
                        corr = np.nan
                else:  # no correlation coefficient.
                    corr = np.nan

                width_fnam = np.max((width_fnam, len(out_name)))

                for x in range(len(fit[y]["peak"])):
                    out_peak = []

                    if fit[y]["peak"][x]["d-space_type"] != "fourier":
                        raise ValueError(
                            "This output type is not setup to process non-Fourier peak centroids."
                        )

                    # peak
                    out_peak = peak_string(fit[y], peak=[x])
                    width_hkl = np.max((width_hkl, len(out_peak)))

                    # %%% get converted values.
                    crystallographic_values = cfc.fourier_to_crystallographic(
                        fit,
                        SampleGeometry=SampleGeometry,
                        SampleDeformation=SampleDeformation,
                        subpattern=y,
                        peak=x,
                    )

                    try:
                        # centroid
                        out_d0 = crystallographic_values["dp"]
                        out_d0err = crystallographic_values["dp_err"]

                        # differential components
                        out_dcos2 = fit[y]["peak"][x]["d-space"][3]
                        out_dsin2 = fit[y]["peak"][x]["d-space"][4]
                        out_dcos2err = fit[y]["peak"][x]["d-space_err"][3]
                        out_dsin2err = fit[y]["peak"][x]["d-space_err"][4]

                        out_dd = crystallographic_values["differential"]
                        out_dderr = crystallographic_values["differential_err"]
                        out_ang = crystallographic_values["orientation"]
                        out_angerr = crystallographic_values["orientation_err"]

                        out_dcorr = (
                            corr  # gmodel.params['peak_0_d3'].correl['peak_0_d4']
                        )

                        # differential max
                        out_dmax = crystallographic_values["d_max"]
                        out_dmaxerr = crystallographic_values["d_max_err"]
                        # differential min
                        out_dmin = crystallographic_values["d_min"]
                        out_dminerr = crystallographic_values["d_min_err"]

                        # centre position (not used)
                        # FIX ME: This is not correct at the time of writing. See function called above.
                        out_x0 = crystallographic_values["x0"]
                        out_y0 = crystallographic_values["y0"]
                    except:
                        out_dcos2 = np.nan
                        out_dsin2 = np.nan
                        out_dcos2err = np.nan
                        out_dsin2err = np.nan
                        out_dd = np.nan
                        out_dderr = np.nan
                        out_ang = np.nan
                        out_angerr = np.nan

                    # height mean
                    if fit[y]["peak"][x]["height_type"] == "fourier":
                        out_h0 = fit[y]["peak"][x]["height"][0]
                        out_h0err = fit[y]["peak"][x]["height_err"][0]
                    else:
                        tot = np.sum(fit[y]["peak"][x]["height"])
                        errsum = np.sqrt(
                            np.sum(np.array(fit[y]["peak"][x]["height_err"]) ** 2)
                        )
                        num = np.shape(fit[y]["peak"][x]["height"])
                        out_h0 = float(tot / num)
                        out_h0err = float(errsum / num)
                    if out_h0 is None:  # catch  'null' as an error
                        out_h0 = np.nan
                    if out_h0err is None:  # catch  'null' as an error
                        out_h0err = np.nan

                    # width mean
                    if fit[y]["peak"][x]["width_type"] == "fourier":
                        out_w0 = fit[y]["peak"][x]["width"][0]
                        out_w0err = fit[y]["peak"][x]["width_err"][0]
                    else:
                        tot = np.sum(fit[y]["peak"][x]["width"])
                        errsum = np.sqrt(
                            np.sum(np.array(fit[y]["peak"][x]["width_err"]) ** 2)
                        )
                        num = np.shape(fit[y]["peak"][x]["width"])
                        out_w0 = float(tot / num)
                        out_w0err = float(errsum / num)
                    if out_w0 is None:  # catch  'null' as an error
                        out_w0 = np.nan
                    if out_w0err is None:  # catch  'null' as an error
                        out_w0err = np.nan

                    # profile mean
                    if fit[y]["peak"][x]["profile_type"] == "fourier":
                        out_p0 = fit[y]["peak"][x]["profile"][0]
                        out_p0err = fit[y]["peak"][x]["profile_err"][0]
                    else:
                        tot = np.sum(fit[y]["peak"][x]["profile"])
                        errsum = np.sqrt(
                            np.sum(np.array(fit[y]["peak"][x]["profile_err"]) ** 2)
                        )
                        num = np.shape(fit[y]["peak"][x]["profile"])
                        out_p0 = float(tot / num)
                        out_p0err = float(errsum / num)
                    if out_p0 is None:  # catch  'null' as an error
                        out_p0 = np.nan
                    if out_p0err is None:  # catch  'null' as an error
                        out_p0err = np.nan

                    # %% write numbers to file

                    text_file.write(
                        ("{0:<" + str(width_fnam) + "}").format(out_name + ",")
                    )
                    text_file.write(
                        ("{0:<" + str(width_hkl) + "}").format(out_peak + ",")
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_d0
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_d0err
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dcos2
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dcos2err
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dsin2
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dsin2err
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dcorr
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dd
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dderr
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_ang
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_angerr
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dmax
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_dmin
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_h0
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_h0err
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_w0
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_w0err
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_p0
                        )
                    )
                    text_file.write(
                        ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                            out_p0err
                        )
                    )
                    if debug:
                        # include properties from the lmfit output that were passed with the fits.
                        text_file.write(
                            ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                                fit[y]["FitProperties"]["time-elapsed"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                                fit[y]["FitProperties"]["chunks-time"]
                            )
                        )
                        text_file.write(
                            (
                                "{0:" + str(width_col - 1) + "." + str(dp - 1) + "e},"
                            ).format(fit[y]["FitProperties"]["sum-residuals-squared"])
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + ".0f},").format(
                                fit[y]["FitProperties"]["status"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + ".0f},").format(
                                fit[y]["FitProperties"]["function-evaluations"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + ".0f},").format(
                                fit[y]["FitProperties"]["n-variables"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + ".0f},").format(
                                fit[y]["FitProperties"]["n-data"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + ".0f},").format(
                                fit[y]["FitProperties"]["degree-of-freedom"]
                            )
                        )
                        text_file.write(
                            (
                                "{0:" + str(width_col - 1) + "." + str(dp - 1) + "e},"
                            ).format(fit[y]["FitProperties"]["ChiSq"])
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                                fit[y]["FitProperties"]["RedChiSq"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                                fit[y]["FitProperties"]["aic"]
                            )
                        )
                        text_file.write(
                            ("{0:" + str(width_col - 1) + "." + str(dp) + "f},").format(
                                fit[y]["FitProperties"]["bic"]
                            )
                        )
                    if "note" in fit[y]:
                        text_file.write(
                            ("{0:<" + str(width_hkl) + "}").format(fit[y]["note"])
                        )

                    text_file.write("\n")
        else:
            logger.info(" ".join(("%s does not exist on the path" % filename)))

    text_file.close()
