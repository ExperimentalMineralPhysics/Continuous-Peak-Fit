__all__ = ["Requirements", "WriteOutput"]


import os

import numpy as np
import cpf.PeakFunctions as ff
import json
import lmfit
from lmfit.model import load_modelresult

# import cpf.XRD_FitPattern as XRD_FP
import cpf.IO_functions as IO


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        "Output_directory"  # if no direcrtory is specified write to current directory.
    ]

    return RequiredParams


def WriteOutput(FitSettings, parms_dict, **kwargs):
    # writes some of the fitted coeficients to a table.
    # More general version of WriteDifferentialStrain.
    # Focused on d-spacing coefficents but generalisible to all coefficients.

    coefs_vals_write = "d-space"

    FitParameters = dir(FitSettings)

    # Parse the required inputs.
    # base_file_name = FitSettings.datafile_Basename

    # diffraction patterns
    diff_files, n_diff_files = IO.file_list(FitParameters, FitSettings)
    # if 'Output_directory' in FitParameters:
    #     out_dir = FitSettings.Output_directory
    # else:
    #     out_dir = './'

    # # create output file name from passed name
    # path, filename = os.path.split(base_file_name)
    # base, ext = os.path.splitext(filename)
    # if not base:
    #     print("No base filename, using input filename instead.")
    #     base =  os.path.splitext(FitSettings.inputfile)[0]

    # if base[-1:] == '_': #if the base file name ends in an '_' remove it.
    #     base = base[0:-1]

    # out_file = out_dir + base + '.dat'

    base, ext = os.path.splitext(os.path.split(FitSettings.datafile_Basename)[1])
    if not base:
        print("No base filename, using input filename instead.")
        base = os.path.splitext(os.path.split(FitSettings.inputfile)[1])[0]
    out_file = IO.make_outfile_name(
        base, directory=FitSettings.Output_directory, extension=".dat", overwrite=True
    )

    # get number of coefficients.
    num_subpatterns = len(FitSettings.fit_orders)
    max_coefs = 0
    for y in range(num_subpatterns):
        for x in range(len(FitSettings.fit_orders[y]["peak"])):
            max_coefs = np.max(
                [max_coefs, FitSettings.fit_orders[y]["peak"][x]["d-space"]]
            )

    text_file = open(out_file, "w")
    print("Writing", out_file)

    text_file.write(
        "# Summary of fits produced by continuous_peak_fit for input file: %s.\n"
        % FitSettings.inputfile
    )
    text_file.write("# For more information: http://www.github.com/me/something\n")
    text_file.write("# File version: %i \n" % 1)
    text_file.write("# \n")

    # write header
    col_width = 15
    text_file.write(("# {0:<" + str(col_width + 5) + "}").format("Data File" + ","))
    text_file.write(("{0:<" + str(col_width) + "}").format("Peak" + ","))
    # coefficient header list
    text_file.write(("{0:>" + str(col_width) + "}").format("d0" + ","))
    for w in range(max_coefs):
        text_file.write(
            ("{0:>" + str(col_width) + "}").format("d" + str(w + 1) + "cos,")
        )
        text_file.write(
            ("{0:>" + str(col_width) + "}").format("d" + str(w + 1) + "sin,")
        )
    # coefficient errors header list
    text_file.write(("{0:>" + str(col_width) + "}").format("d0" + " error,"))
    for w in range(max_coefs):
        text_file.write(
            ("{0:>" + str(col_width) + "}").format("d" + str(w + 1) + "cos error,")
        )
        text_file.write(
            ("{0:>" + str(col_width) + "}").format("d" + str(w + 1) + "sin error,")
        )

    text_file.write("\n")

    for z in range(n_diff_files):

        filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        filename = filename + ".json"

        # Read JSON data from file
        with open(filename) as json_data:
            fit = json.load(json_data)

        num_subpatterns = len(FitSettings.fit_orders)
        for y in range(num_subpatterns):

            orders = FitSettings.fit_orders[y]
            """
            subfilename = os.path.splitext(os.path.basename(diff_files[z]))[0] + '_'
            
            for x in range(len(orders['peak'])):
                if 'phase' in orders['peak'][x]:
                    subfilename = subfilename + orders['peak'][x]['phase']
                else:
                    subfilename = subfilename + "Peak"
                if 'hkl' in orders['peak'][x]:
                    subfilename = subfilename + str(orders['peak'][x]['hkl'])
                else:
                    subfilename = subfilename + x
                        
                if x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
                    subfilename = subfilename + '_'
                        
            print('  Incorporating ' + subfilename)
            try:
                gmodel = load_modelresult(subfilename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
            except:
                IO.lmfit_fix_int_data_type(subfilename+'.sav')
                try:
                    gmodel = load_modelresult(subfilename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
                except:
                    raise FileNotFoundError     
        
            # FIX ME: this will only work for one peak. Needs fixing if more then one peak in subpattern
            corr = gmodel.params['peak_0_d3'].correl['peak_0_d4']
            #ci, trace = lmfit.conf_interval(mini, gmodel, sigmas=[1, 2], trace=True)
            #lmfit.printfuncs.report_ci(ci)
            """

            # get parameters to write to output file.
            # file name
            out_name = os.path.splitext(os.path.basename(diff_files[z]))[0]

            for x in range(len(orders["peak"])):

                out_peak = []

                # peak
                if "phase" in orders["peak"][x]:
                    out_peak = orders["peak"][x]["phase"]
                else:
                    out_peak = out_peak + "Peak"
                if "hkl" in orders["peak"][x]:
                    out_peak = out_peak + " (" + str(orders["peak"][x]["hkl"]) + ")"
                else:
                    out_peak = out_peak + " " + str(x)
                """
                # d0
                out_d0       = fit[y]['peak'][x]['d-space'][0]
                out_d0err    = fit[y]['peak'][x]['d-space_err'][0]
                #differential coefficients, errors and covarience
                out_dcos2    = fit[y]['peak'][x]['d-space'][3]
                out_dsin2    = fit[y]['peak'][x]['d-space'][4]
                out_dcos2err = fit[y]['peak'][x]['d-space_err'][3]
                out_dsin2err = fit[y]['peak'][x]['d-space_err'][4]
                out_dcorr    = gmodel.params['peak_0_d3'].correl['peak_0_d4']
                
                #differential strain
                # a.cos(2t) + b.sin(2t) = (a^2+b^2)^(1/2) . cos(2t - atan(b/a))
                # out_dd  = (a^2+b^2)^(1/2)
                out_dd    = (fit[y]['peak'][x]['d-space'][3]**2 + fit[y]['peak'][x]['d-space'][4]**2)**(1/2)
                #out_dderr= [(2.a.da.)^2 + (2.b.db)^2]^(1/2)]^(1/2)
                out_dderr = ((2 * fit[y]['peak'][x]['d-space'][3] * fit[y]['peak'][x]['d-space_err'][3])**2 + 
                             (2 * fit[y]['peak'][x]['d-space'][4] * fit[y]['peak'][x]['d-space_err'][4])**2 ) ** (1/4)
                
                #angle
                # out_ang  = atan(b/a)
                out_ang    = (np.arctan(fit[y]['peak'][x]['d-space'][4] / fit[y]['peak'][x]['d-space'][3]))
                # d (atan(c))/dc = 1/(c^2+1). c = b/a. dc = c.((da/a)^2 + (db/b)^2)^(1/2)
                # out_angerr = dc. 
                out_angerr = 1 / ((fit[y]['peak'][x]['d-space'][4] / fit[y]['peak'][x]['d-space'][3])**2 + 1) * (np.abs(out_ang) * ((fit[y]['peak'][x]['d-space_err'][3]/fit[y]['peak'][x]['d-space'][3])**2 + (fit[y]['peak'][x]['d-space_err'][4]/fit[y]['peak'][x]['d-space'][4])**2)**(1/2))
                # correction to make angle correct (otherwise out by pi)
                if fit[y]['peak'][x]['d-space'][3] < 0:
                    out_ang = out_ang + np.pi
                # force angles to be between -pi and pi
                if out_ang > np.pi/2:
                    out_ang = out_ang - np.pi
                elif out_ang < -np.pi/2:
                    out_ang = out_ang + np.pi
                # convert angles to degrees.
                out_ang    = np.rad2deg(out_ang)
                out_angerr = np.rad2deg(out_angerr)
                
                
                #differential max
                out_dmax    = out_dd * np.cos(2*np.pi/4 - np.deg2rad(out_ang)) + out_d0
                #differential min 
                out_dmin    = out_dd * np.cos(2*3*np.pi/4 - np.deg2rad(out_ang)) + out_d0
                
                #height mean
                out_h0    = fit[y]['peak'][x]['height'][0]
                out_h0err = fit[y]['peak'][x]['height_err'][0]
                #width mean
                out_w0    = fit[y]['peak'][x]['width'][0]
                out_w0err = fit[y]['peak'][x]['width_err'][0]
                #profile mean
                out_p0    = fit[y]['peak'][x]['profile'][0]
                out_p0err = fit[y]['peak'][x]['profile_err'][0]
                """

                # write numbers to file
                dp = col_width - 7
                text_file.write(
                    ("{0:<" + str(col_width + 7) + "}").format(out_name + ",")
                )
                text_file.write(("{0:<" + str(col_width) + "}").format(out_peak + ","))

                # coefficients
                text_file.write(
                    ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                        fit[y]["peak"][x]["d-space"][0]
                    )
                )
                for w in range(max_coefs):
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x]["d-space"][2 * w + 1]
                        )
                    )
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x]["d-space"][2 * w + 2]
                        )
                    )
                # coefficient errors
                text_file.write(
                    ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                        fit[y]["peak"][x]["d-space_err"][0]
                    )
                )
                for w in range(max_coefs):
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x]["d-space_err"][2 * w + 1]
                        )
                    )
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x]["d-space_err"][2 * w + 2]
                        )
                    )

                text_file.write("\n")

    # col_width = 10
    # text_file.write(("# {0:<"+str(col_width+5)+"}").format("Data File"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("Peak"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("d0"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("d0_err"+","))

    text_file.close()
