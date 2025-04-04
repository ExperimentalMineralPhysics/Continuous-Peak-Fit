__all__ = ["Requirements", "WriteOutput"]


import os
import numpy as np
import json
from itertools import product
# import cpf.XRD_FitPattern as XRD_FP
import cpf.IO_functions as IO
import cpf.peak_functions as pf
from cpf.XRD_FitPattern import logger
import pandas as pd


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
def WriteOutput(setting_class=None, setting_file=None, debug=False, writefile=True, *args, **kwargs):
    # writes some of the fitted coeficients to a table.
    # More general version of WriteDifferentialStrain.
    # Focused on d-spacing coefficents but generalisible to all coefficients.

    #define defaults
    dp = 6 #how many decimal points to write out
    col_width = 15 # default column width for csv file.
    
    ## output file version
    # 1: ?
    # 2: ?
    # 3: rewritten as panda data frame - to force columns to line up.   
    version = 3
        
    
    ordering_of_output = "peak"

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
    if isinstance(coefs_vals_write, str):
        coefs_vals_write = [coefs_vals_write]
    if coefs_vals_write == ["all"]:
        peak_properties = pf.peak_components(full=True)
        coefs_vals_write = peak_properties[1]

    #make filename for output
    base = setting_class.datafile_basename
    if base is None:
        logger.info(" ".join(map(str, [("No base filename, using input filename instead.")])))
        base = os.path.splitext(os.path.split(setting_class.settings_file)[1])[0]
    out_file = IO.make_outfile_name(
        base, directory=setting_class.output_directory, extension=".dat", overwrite=True, additional_text="all_coefficients",
    )

    # get number of coefficients.
    num_subpatterns = len(setting_class.fit_orders)
    max_coef = {}
    for w in range(len(coefs_vals_write)):
        ind = coefs_vals_write[w]
        if ind != "background":
            max_coef[coefs_vals_write[w]] = 0
        else:
            max_coef[coefs_vals_write[w]] = [0]
    max_coefs = 0
    for y in range(num_subpatterns):
        for x in range(len(setting_class.fit_orders[y]["peak"])):

            for w in range(len(coefs_vals_write)):
                ind = coefs_vals_write[w]
                if ind != "background" and ind != "symmetry":
                    max_coef[ind] = np.max(
                        [max_coef[ind], np.max(2*setting_class.fit_orders[y]["peak"][x][ind])+1]
                    )
                elif ind == "symmetry":
                    # symmetry might not be present in the fitting files.
                    try:
                        max_coef[ind] = np.max(
                            [max_coef[ind], np.max(2*setting_class.fit_orders[y]["peak"][x][ind])+1]
                        )
                    except:
                        pass                    
                else:
                    for v in range(len(setting_class.fit_orders[y][ind])):
                        if v > len(max_coef[ind])-1:
                            max_coef[ind].append(0)
                        max_coef[ind][v] = np.max(
                            [max_coef[ind][v], 2*setting_class.fit_orders[y][ind][v]+1]
                        )
                
    # store headers
    headers = []
    headers.append("DataFile")
    headers.append("Peak")
    
    #parmeter header list
    for w in range(len(max_coef)):
        ind = coefs_vals_write[w]
        if ind != "background" and ind != "symmetry":
            for v in range(max_coef[ind]):
                headers.append(ind + str(v))
                headers.append(ind + str(v)+ "err")
        elif ind == "symmetry":
            headers.append(ind)
        else:
            for u in range(len(setting_class.fit_orders[y][ind])):
                
                for v in range(max_coef[ind][u]):
                        headers.append(ind + str(u) + "_" + str(v))
                        headers.append(ind + str(u) + "_" + str(v)+ "err")

    # read all the data. 
    fits=[]
    for z in range(setting_class.image_number):
        setting_class.set_subpattern(z, 0)

        filename = IO.make_outfile_name(
            setting_class.subfit_filename,  # diff_files[z],
            directory=setting_class.output_directory,  # directory=FitSettings.Output_directory,
            extension=".json",
            overwrite=True,
        )  # overwrite =false to get the file name without incrlemeting it.

        # Read JSON data from file
        with open(filename) as json_data:
            fits.append(json.load(json_data))
    
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
    
    
    RowsList = []
    #print the data.
    for z in range(len(lists)):
    
        RowLst = {}   
    
        setting_class.set_subpattern(lists[z,0],lists[z,1])
        data_to_write = fits[lists[z,0]][lists[z,1]]
        
        
        if len(data_to_write["peak"]) > lists[z,2]:
            
            out_name = IO.make_outfile_name(
                setting_class.subfit_filename,
                directory="",
                extension=".json",
                overwrite=True,
            )
            
            out_peak = []
    
            # peak
            if "phase" in setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]:
                out_peak = setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]["phase"]
            else:
                out_peak = out_peak + "Peak"
            if "hkl" in setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]:
                out_peak = out_peak + " (" + str(setting_class.fit_orders[lists[z,1]]["peak"][lists[z,2]]["hkl"]) + ")"
            else:
                out_peak = out_peak + " " + str([lists[z,2]])
            RowLst["DataFile"] = out_name
            RowLst["Peak"] = out_peak
            
            for w in range(len(coefs_vals_write)):
                ind = coefs_vals_write[w]
                ind_err = ind+"_err"
                
                if ind != "background" and ind != "symmetry":   
                    for v in range(len(data_to_write["peak"][lists[z,2]][ind])):
                        if data_to_write["peak"][lists[z,2]][ind][v] is None:  # catch  'null' as an error
                            data_to_write["peak"][lists[z,2]][ind][v] = np.nan
                        if data_to_write["peak"][lists[z,2]][ind_err][v] is None:  # catch  'null' as an error
                            data_to_write["peak"][lists[z,2]][ind_err][v] = np.nan   
                        RowLst[ind+str(v)] = data_to_write["peak"][lists[z,2]][ind][v]
                        try:
                            RowLst[ind+str(v)+"err"] = data_to_write["peak"][lists[z,2]][ind_err][v]
                        except:
                            RowLst[ind+str(v)+"err"] = "None"
                
                        
                        
                                
                elif ind == "symmetry":
                    try:
                        if data_to_write["peak"][lists[z,2]][ind] is None:  # catch  'null' as an error
                            data_to_write["peak"][lists[z,2]][ind] = np.nan
                        RowLst[ind] = data_to_write["peak"][lists[z,2]][ind]
                    except:
                        pass
                    
                else: #background
                    for u in range(len(data_to_write[ind])):
                        for v in range(len(data_to_write[ind][u])):
                            if data_to_write[ind][u][v] is None:  # catch  'null' as an error
                                data_to_write[ind][u][v] = np.nan
                            if data_to_write[ind_err][u][v] is None:  # catch  'null' as an error
                                data_to_write[ind_err][u][v] = np.nan 
                            RowLst[ind+str(u)+"_"+str(v)] = data_to_write[ind][u][v]
                            try:
                                RowLst[ind+str(u)+"_"+str(v)+"err"] = data_to_write[ind_err][u][v]
                            except:
                                RowLst[ind+str(u)+"_"+str(v)+"err"] = "None"
        
            RowsList.append(RowLst)
            # text_file.write("\n")
   
    # make data frame using headers - so columns are in good order.
    df = pd.DataFrame(RowsList,columns=headers)
    # keep unmodificed version of the dataframe to output
    out = df    
    
    if writefile == True:
        ## format dateframe for writing to file neatly. 
        #make strings in DateFile and Peak columns all the same length
        len_datafile = np.max(df["DataFile"].str.len())
        len_peaks    = np.max(df["Peak"].str.len())
        df_tmp = df["DataFile"].str.pad(np.max([len_datafile, col_width]),side='left',fillchar=' ')
        df["DataFile"] = df_tmp
        df_tmp = df["Peak"].str.pad(np.max([len_peaks, col_width]),side='left',fillchar=' ')
        df["Peak"] = df_tmp
        
        # rename the columns so that the headers are the same width as the columns
        class NewClass(object): pass
        columns = NewClass()
        for i in range(len(headers)):
            if headers[i] == "DataFile":
                setattr(columns, headers[i], headers[i].rjust(np.max([len_datafile, col_width])))
            elif headers[i] == "Peak":
                setattr(columns, headers[i], headers[i].rjust(np.max([len_peaks, col_width])))
            else:
                setattr(columns, headers[i], headers[i].rjust(col_width))
        columns = columns.__dict__
        df.rename(columns=columns, inplace=True)
            
        # write data frame to csv file
        with open(out_file, 'w') as f:
            logger.info(" ".join(map(str, [("Writing %s" % out_file)])))
    
            f.write(
                "# continuous_peak_fit : Table of all coefficients from fits for input file: %s.\n"
                % setting_class.settings_file
            )
            f.write("# For more information: http://www.github.com/me/something\n")
            f.write("# File version: %i \n" % version)
            f.write("# \n")
            
            df.to_csv(f, index=False, header=True,
                      na_rep = "".ljust(col_width),
                      float_format = lambda x: f"{x:{str(col_width)}d}" if isinstance(x, (int, np.integer)) else (f"{x:{str(col_width)}.{str(dp)}f}" if np.abs(x) > 0.1 else f"{x:{str(col_width)}.{str(dp)}e}"))
    
    return out