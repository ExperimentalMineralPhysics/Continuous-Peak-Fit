__all__ = ["Requirements", "WriteOutput"]


import os
import numpy as np
import json
from itertools import product
# import cpf.XRD_FitPattern as XRD_FP
import cpf.IO_functions as IO
import cpf.peak_functions as pf


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        ##"Output_directory"  # if no direcrtory is specified write to current directory.
        "coefs_vals_write"  # -- pick which set of coefficients to write
    ]

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(setting_class=None, setting_file=None, debug=False, **kwargs):
    # writes some of the fitted coeficients to a table.
    # More general version of WriteDifferentialStrain.
    # Focused on d-spacing coefficents but generalisible to all coefficients.

    if setting_class is None and setting_file is None:
        raise ValueError(
            "bummer Either the settings file or the setting class need to be specified."
        )
    elif setting_class is None:
        import cpf.XRD_FitPattern.initiate as initiate

        setting_class = initiate(setting_file)

    #get what to write
    if 'coefs_write' in kwargs:
        coefs_vals_write = kwargs.get('coefs_write')
    elif "coefs_vals_write" in setting_class.output_settings:
        coefs_vals_write = setting_class.output_settings["coefs_vals_write"]
    else:
        coefs_vals_write = "all"    


    #make filename for output
    base = setting_class.datafile_basename
    if base is None:
        print("No base filename, using input filename instead.")
        base = os.path.splitext(os.path.split(setting_class.settings_file)[1])[0]
    # base, ext = os.path.splitext(os.path.split(FitSettings.datafile_Basename)[1])
    # if not base:
    #     print("No base filename, using input filename instead.")
    #     base = os.path.splitext(os.path.split(FitSettings.inputfile)[1])[0]
    out_file = IO.make_outfile_name(
        base, directory=setting_class.output_directory, extension=".dat", overwrite=True, additional_text="all_coefficients",
    )

    # get number of coefficients.
    num_subpatterns = len(setting_class.fit_orders)
    peak_properties = pf.peak_components(full=True)
    if coefs_vals_write == "all":
        coefs_vals_write = peak_properties
    max_coef = {}
    for w in range(len(coefs_vals_write[1])):
        ind = coefs_vals_write[1][w]
        if ind != "background":
            max_coef[coefs_vals_write[1][w]] = 0
        else:
            max_coef[coefs_vals_write[1][w]] = [0]
    max_coefs = 0
    for y in range(num_subpatterns):
        for x in range(len(setting_class.fit_orders[y]["peak"])):

            for w in range(len(coefs_vals_write[1])):
                ind = coefs_vals_write[1][w]
                if ind != "background":
                    max_coef[ind] = np.max(
                        [max_coef[ind], 2*setting_class.fit_orders[y]["peak"][x][ind]+1]
                    )
                else:
                    for v in range(len(setting_class.fit_orders[y][ind])):
                        if v > len(max_coef[ind])-1:
                            max_coef[ind].append(0)
                        max_coef[ind][v] = np.max(
                            [max_coef[ind][v], 2*setting_class.fit_orders[y][ind][v]+1]
                        )
                
            #     max_coef[coefs_vals_write[1][w]] = 0
            # print(coefs_vals_write[1][y])
            # max_coef[coefs_vals_write[1][y]] = 0
            
            
            # max_coefs = np.max(
            #     [max_coefs, setting_class.fit_orders[y]["peak"][x]["d-space"]]
            # )
            

    text_file = open(out_file, "w")
    print("Writing", out_file)

    text_file.write(
        "# Summary of fits produced by continuous_peak_fit for input file: %s.\n"
        % setting_class.settings_file
    )
    text_file.write("# For more information: http://www.github.com/me/something\n")
    text_file.write("# File version: %i \n" % 2)
    text_file.write("# \n")

    # write header
    col_width = 17
    dp = col_width - 7
    text_file.write(("# {0:<" + str(col_width + 5) + "}").format("Data File" + ","))
    text_file.write(("{0:<" + str(col_width) + "}").format("Peak" + ","))
    
    #parmeter header list
    
    for w in range(len(max_coef)):
        ind = coefs_vals_write[1][w]
        if ind != "background" and ind != "symmetry":
            for v in range(max_coef[ind]):
                text_file.write(
                    ("{0:>" + str(col_width) + "}").format(ind + "_" + str(v)+ ",")
                )
                text_file.write(
                    ("{0:>" + str(col_width) + "}").format(ind + "_" + str(v)+ " err,")
                )
        elif ind == "symmetry":
            text_file.write(
                ("{0:>" + str(col_width) + "}").format(ind + ",")
            )
        else:
            for u in range(len(setting_class.fit_orders[y][ind])):
                
                for v in range(max_coef[ind][u]):
                        text_file.write(
                            ("{0:>" + str(col_width) + "}").format(ind + str(u) + "_" + str(v) + ",")
                        )
                        text_file.write(
                            ("{0:>" + str(col_width) + "}").format(ind + str(u) + "_" + str(v) + " err,")
                        )

    text_file.write("\n")




    # read all the data. 
    fits=[]
    for z in range(setting_class.image_number):
        setting_class.set_subpattern(z, 0)

        filename = IO.make_outfile_name(
            # diff_files[z],
            # directory=FitSettings.Output_directory,
            setting_class.subfit_filename,  # diff_files[z],
            directory=setting_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".json",
            overwrite=True,
        )  # overwrite =false to get the file name without incrlemeting it.

        # filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        # filename = filename + ".json"

        # Read JSON data from file

        with open(filename) as json_data:
            fits.append(json.load(json_data))

    #print(fits[0][1])    
    #stop
    ordering_of_output = "peak"
    
    #make lists of the parameters to iterate over
    images = list(range(setting_class.image_number))
    subpatterns = list(range(num_subpatterns))
    max_peaks = 1
    for i in range(len(setting_class.fit_orders)):
        max_peaks = np.max([max_peaks, len(setting_class.fit_orders[i]["peak"])])

    
    lists = images, subpatterns, list(range(max_peaks))
    lists = np.array(list(product(*lists)))
    
    #sort lists
    if ordering_of_output == "peak":
        lists=lists[np.lexsort((lists[:, 2], ))]
        lists=lists[np.lexsort((lists[:, 1], ))]
    else:
        lists=lists[np.lexsort((lists[:, 0], ))]
    
    #print the data.
    for z in range(len(lists)):
    
        setting_class.set_subpattern(lists[z,0],lists[z,1])
        data_to_write = fits[lists[z,0]][lists[z,1]]
        
        #print(data_to_write)
        #print(len(data_to_write))
        
        
        if len(data_to_write["peak"]) > lists[z,2]:
            
            
             
            out_name = IO.make_outfile_name(
                setting_class.subfit_filename,  # diff_files[z],
                directory="",
                extension=".json",
                overwrite=True,
            )
            # print(data_to_write["peak"])
            # print(len(data_to_write["peak"]))
            # print(len(data_to_write["peak"]) , lists[z,2])
            out_peak = []
    
            # peak
            print(setting_class.fit_orders[lists[z,1]])
            print(setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]])
            if "phase" in setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]:
                out_peak = setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]["phase"]
            else:
                out_peak = out_peak + "Peak"
            if "hkl" in setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]:
                out_peak = out_peak + " (" + str(setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]["hkl"]) + ")"
            else:
                out_peak = out_peak + " " + str([lists[z,2]])
            # write numbers to file
            dp = col_width - 7
            text_file.write(
                ("{0:<" + str(col_width + 7) + "}").format(out_name + ",")
            )
            text_file.write(("{0:<" + str(col_width) + "}").format(out_peak + ","))
            
            
            print(z)
            print(len(data_to_write["peak"]) , lists[z,2])
             #     print(lists[z,:])
        
        
            for w in range(len(coefs_vals_write[1])):
                ind = coefs_vals_write[1][w]
                ind_err = ind+"_err"
                
                if ind != "background" and ind != "symmetry":   
                    print(data_to_write["peak"], [lists[z,2]], ind)
                    for v in range(len(data_to_write["peak"][lists[z,2]][ind])):
                                   
                        text_file.write(
                            ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                                data_to_write["peak"][lists[z,2]][ind][v]
                            )
                        )
                        try:
                            text_file.write(
                                ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                                    data_to_write["peak"][lists[z,2]][ind_err][v]
                                )
                            )
                        except:
                            text_file.write(("{0:<" + str(col_width) + "}").format("None" + ","))
                            
                                
                elif ind == "symmetry":
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            data_to_write["peak"][lists[z,2]][ind]
                        )
                    )
                else:
                    for u in range(len(data_to_write[ind])):
                        print(data_to_write[ind])
                        for v in range(len(data_to_write[ind][u])):
                            
                                text_file.write(
                                    ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                                        data_to_write[ind][u][v]
                                    )
                                )
                                text_file.write(
                                    ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                                        data_to_write[ind][u][v]
                                    )
                                )
        
        
            text_file.write("\n")
    text_file.close()
        
        # print(lists[z,:])
        # print(len(fits[lists[z,0]][lists[z,1]]["peak"]), lists[z,2])
        
        # if len(fits[lists[z,0]][lists[z,1]]["peak"])
        # for w in range(len(max_coef)):
        #     ind = coefs_vals_write[1][w]
        #     if ind != "background" and ind != "symmetry":
        #         for v in range(max_coef[ind]):
                    
        #             text_file.write(
        #                 ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
        #                     fits[lists[z,0]][y]["peak"][x][coefs_vals_write][2 * w + 1]
        #                 )
        #             )
                    
                    
        #             text_file.write(
        #                 ("{0:>" + str(col_width) + "}").format(ind + "_" + str(v)+ ",")
        #             )
        #             text_file.write(
        #                 ("{0:>" + str(col_width) + "}").format(ind + "_" + str(v)+ " err,")
        #             )
        #     elif ind == "symmetry":
        #         text_file.write(
        #             ("{0:>" + str(col_width) + "}").format(ind + ",")
        #         )
        #     else:
        #         for u in range(len(setting_class.fit_orders[y][ind])):
                    
        #             for v in range(max_coef[ind][u]):
        #                     text_file.write(
        #                         ("{0:>" + str(col_width) + "}").format(ind + str(u) + "_" + str(v) + ",")
        #                     )
        #                     text_file.write(
        #                         ("{0:>" + str(col_width) + "}").format(ind + str(u) + "_" + str(v) + " err,")
        #                     )
        
        
        
        

    stop

           # get parameters to write to output file.
           # strip directory and ending from filename
           # out_name = os.path.splitext(os.path.basename(diff_files[z]))[0]
           # out_name = os.path.splitext(os.path.basename(setting_class.datafile_list[z]))[0]

           # setting_class.set_subpattern(z,0)

           # out_name = IO.make_outfile_name(
           #     setting_class.subfit_filename,  # diff_files[z],
           #     directory=setting_class.output_directory,
           #     extension=".json",
           #     overwrite=True,
           # )




    
    stop    
    
    
    
    
    lists = images, subpatterns
    a = product(*lists)        
    print(a,b)
    stop
    print(images)
    print(subpatterns)
    
    #sort order to write over
    if ordering_of_output == "peak":
        lists = subpatterns, images
    else:
        lists = images, subpatterns

    print(product(*lists))
    for x,y in product(*lists):
        
        print(x,y)

        setting_class.set_subpattern(z, 0)
        
        
        
    

    print(arr)
    stop









    for z in range(setting_class.image_number):

        filename = IO.make_outfile_name(
            # diff_files[z],
            # directory=FitSettings.Output_directory,
            setting_class.subfit_filename,  # diff_files[z],
            directory=setting_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".json",
            overwrite=True,
        )  # overwrite =false to get the file name without incrlemeting it.

        # filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        # filename = filename + ".json"

        # Read JSON data from file
        with open(filename) as json_data:
            fit = json.load(json_data)

        for y in range(num_subpatterns):

            orders = setting_class.fit_orders[y]
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
            # strip directory and ending from filename
            # out_name = os.path.splitext(os.path.basename(diff_files[z]))[0]
            # out_name = os.path.splitext(os.path.basename(setting_class.datafile_list[z]))[0]

            setting_class.set_subpattern(z,0)

            out_name = IO.make_outfile_name(
                setting_class.subfit_filename,  # diff_files[z],
                directory=setting_class.output_directory,
                extension=".json",
                overwrite=True,
            )

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
                print(y,x,coefs_vals_write)
                print(fit[y]["peak"][x][coefs_vals_write][0])
                text_file.write(
                    ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                        fit[y]["peak"][x][coefs_vals_write][0]
                    )
                )
                for w in range(max_coefs):
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x][coefs_vals_write][2 * w + 1]
                        )
                    )
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x][coefs_vals_write][2 * w + 2]
                        )
                    )
                # coefficient errors
                text_file.write(
                    ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                        fit[y]["peak"][x][coefs_vals_write + "_err"][0]
                    )
                )
                for w in range(max_coefs):
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x][coefs_vals_write + "_err"][2 * w + 1]
                        )
                    )
                    text_file.write(
                        ("{0:" + str(col_width - 1) + "." + str(dp) + "f},").format(
                            fit[y]["peak"][x][coefs_vals_write + "_err"][2 * w + 2]
                        )
                    )

                text_file.write("\n")

    # col_width = 10
    # text_file.write(("# {0:<"+str(col_width+5)+"}").format("Data File"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("Peak"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("d0"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("d0_err"+","))

    text_file.close()
