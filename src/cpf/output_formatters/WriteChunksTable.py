__all__ = ["Requirements", "WriteOutput"]

import os
import cpf.IO_functions as IO
from cpf.XRD_FitPattern import logger
from cpf.Cascade import read_saved_chunks
import pandas as pd

def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        "coefs_vals_write"  # -- pick which set of coefficients to write
        "write_order"
    ]

    return RequiredParams, OptionalParams


# def WriteOutput(FitSettings, parms_dict, **kwargs):
def WriteOutput(setting_class=None, setting_file=None, debug=False, *args, **kwargs):
    # writes some of the fitted coeficients to a table.
    # More general version of WriteDifferentialStrain.
    # Focused on d-spacing coefficents but generalisible to all coefficients.

    if setting_class is None and setting_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
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
        
    write_order = ["peak", "file"]
    if "write_order" in kwargs:
        write_order = kwargs.get('write_order')

    #make filename for output
    base = setting_class.datafile_basename
    if base is None:
        logger.info(" ".join(map(str, [("No base filename, using input filename instead.")])))
        base = os.path.splitext(os.path.split(setting_class.settings_file)[1])[0]
    out_file = IO.make_outfile_name(
        base, directory=setting_class.output_directory, extension=".dat", overwrite=True, additional_text="chunkfits",
    )

    # get values from chunk files    
    all_azis, all_fits = read_saved_chunks(setting_class=setting_class)
    
    #integrate the data and add to 

    # create the data table.
    azm = []
    tbl = []
    tbl_errs = []
    for i in range(setting_class.image_number):
        for j in range(len(setting_class.fit_orders)):
            for k in range(len(setting_class.fit_orders[j]["peak"])):
                azm.append([
                    IO.make_outfile_name(setting_class.image_list[i],directory="",overwrite=True), 
                    IO.peak_string(setting_class.fit_orders[j], peak=k), 
                    all_azis[i][j],
                    # all_fits[i][j]["d"][k], 
                    # all_fits[i][j]["h"][k], 
                    # all_fits[i][j]["w"][k], 
                    # all_fits[i][j]["p"][k] 
                    ])
                tbl.append([
                    IO.make_outfile_name(setting_class.image_list[i],directory="",overwrite=True), 
                    IO.peak_string(setting_class.fit_orders[j], peak=k), 
                    #all_azis[i][j], 
                    all_fits[i][j]["d"][k], 
                    all_fits[i][j]["h"][k], 
                    all_fits[i][j]["w"][k], 
                    all_fits[i][j]["p"][k] 
                    ])
                tbl_errs.append([
                    IO.make_outfile_name(setting_class.image_list[i],directory="",overwrite=True), 
                    IO.peak_string(setting_class.fit_orders[j], peak=k), 
                    #all_azis[i][j], 
                    all_fits[i][j]["d_err"][k], 
                    all_fits[i][j]["h_err"][k], 
                    all_fits[i][j]["w_err"][k], 
                    all_fits[i][j]["p_err"][k] 
                    ])

    # get number of coefficients.
    # num_subpatterns = len(setting_class.fit_orders)
    # peak_properties = pf.peak_components(full=True)
    # if coefs_vals_write == "all":
    #     coefs_vals_write = peak_properties
    # max_coef = {}
    # for w in range(len(coefs_vals_write[1])):
    #     ind = coefs_vals_write[1][w]
    #     if ind != "background":
    #         max_coef[coefs_vals_write[1][w]] = 0
    #     else:
    #         max_coef[coefs_vals_write[1][w]] = [0]
    # max_coefs = 0
    # for y in range(num_subpatterns):
    #     for x in range(len(setting_class.fit_orders[y]["peak"])):

    #         for w in range(len(coefs_vals_write[1])):
    #             ind = coefs_vals_write[1][w]
    #             if ind != "background":
    #                 max_coef[ind] = np.max(
    #                     # [max_coef[ind], np.max(setting_class.fit_orders[y]["peak"][x][ind])]
    #                     [max_coef[ind], np.max(2*setting_class.fit_orders[y]["peak"][x][ind])+1]
    #                 )
    #             else:
    #                 for v in range(len(setting_class.fit_orders[y][ind])):
    #                     if v > len(max_coef[ind])-1:
    #                         max_coef[ind].append(0)
    #                     max_coef[ind][v] = np.max(
    #                         [max_coef[ind][v], 2*setting_class.fit_orders[y][ind][v]+1]
    #                     )

    # make a pnada data frame
    azm = pd.DataFrame(columns=["file", "peak", "azimuths"], data=azm)
    tbl = pd.DataFrame(columns=["file", "peak", "d-space", "height","width", "profile"], data=tbl)
    tbl_errs = pd.DataFrame(columns=["file", "peak", "d-space_err", "height_err","width_err", "profile_err"], data=tbl_errs)
    # sort the data table 
    azm = azm.sort_values(write_order)
    tbl = tbl.sort_values(write_order)
    tbl_errs = tbl_errs.sort_values(write_order)
    
    #write the table to file
    logger.info(" ".join(map(str, [("Writing %s" % out_file)])))
    with open(out_file, 'w') as text_file:
        text_file.write(
            "# Summary of chunk fits produced by continuous_peak_fit.Cascade for input file: %s.\n"
            % setting_class.settings_file
        )
        text_file.write("# For more information: http://www.github.com/me/something\n")
        text_file.write("# File version: %i \n" % 1)
        text_file.write("# The columns are separated by ';' while values for each parameter are separated by ','. \n")
        text_file.write("# \n")
        azm.to_csv(text_file, index=False, sep=";")
        tbl.to_csv(text_file, index=False, sep=";")
        tbl_errs.to_csv(text_file, index=False, sep=";")