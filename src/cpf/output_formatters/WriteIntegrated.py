__all__ = ["Requirements", "WriteOutput"]


import os
import numpy as np
from cpf.IO_functions import make_outfile_name, title_file_names
from cpf.util.logging import get_logger
from cpf.histograms import histogram1d
from cpf.data_preprocess import remove_cosmics as cosmicsimage_preprocess

logger = get_logger("cpf.output_types.WriteIntegrated")

def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = ["bin_n", "histogram_type"]

    return RequiredParams, OptionalParams


def WriteOutput(settings_class=None, setting_file=None, **kwargs):
    """
    Writes integrated data to text file. The default settings write a file in which 
    the data is binned into data of constant number of pixels contributing. The
    withs of each bin therefore change with two theta angle for rectangular area 
    detectors. 
    
    N.B. this output requires the data files to be present to work.

    Parameters
    ----------
    FitSettings : TYPE
        DESCRIPTION.
    parms_dict : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    histogram_type = kwargs.get("histogram_type", "data")
    bin_n          = kwargs.get("bin_n", None)
    # bin_data       = kwargs.get("bin_n", 200)

    if np.array(bin_n).size != 1:
        raise NotImplementedError("Two dimensional historgrams are not implemented yet.")
    
    # set default values if dasta in bin is not given.
    if bin_n == None:
        if histogram_type == "data":
            bin_n = 400
        elif histogram_type == "width":
            bin_n = 2000
    
    if settings_class is None and setting_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        from cpf.XRD_FitPattern import initiate

        settings_class = initiate(setting_file)

    # make the base file name
    if setting_file:
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    else:
        base = settings_class.datafile_basename
    if base is None or len(base) == 0:
        logger.info(
            " ".join(
                map(
                    str,
                    [("No base filename, trying ending without extension instead.")],
                )
            )
        )
        base = settings_class.datafile_ending
    if base is None:
        logger.info(
            " ".join(map(str, [("No base filename, using input filename instead.")]))
        )
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]

    # make the data class.
    data_to_fill = settings_class.image_list[0]
    new_data = settings_class.data_class
    new_data.fill_data(
        data_to_fill,
        settings=settings_class,
    )
    
    # integrate the data
    for z in range(settings_class.image_number):

        # Get diffraction pattern to process.
        new_data.import_image(settings_class.image_list[z])
        settings_class.set_subpattern(z, 0)

        if settings_class.datafile_preprocess is not None:
            # needed because image preprocessing adds to the mask and is different for each image.
            new_data.mask_restore()
            if "cosmics" in settings_class.datafile_preprocess:
                new_data = cosmicsimage_preprocess(new_data, settings_class)
        else:
            # nothing is done here.
            pass

        p, i, a = histogram1d(new_data.tth, new_data.intensity, bin_n=bin_n, histogram_type=histogram_type)

        # GSAS-ii cannot take negative intensities. 
        # Force positive intensities here incase darkfield applied to images incorrectly.
        moveby = np.nanmin(i)
        if np.nanmin(i) < 0:
            i = i - moveby
        else: 
            moveby = 0


        filename = make_outfile_name(
            settings_class.subfit_filename,
            directory=settings_class.output_directory,
            extension=".xy",
            overwrite=True,
        )
        
        
        #make a header
        header = f"# Integrated data from: {title_file_names(settings_class)}. \n#\n"
        header += "# Integrated using Continuous Peak Fit \n"
        header += "#     https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit \n"
        header += f"# from inputfile: \n#     {settings_class.settings_file} \n#\n"
        header += f"# Historgram bin type: '{histogram_type}'\n"
        if histogram_type == "data":
            header += "\n# The bins are all made with approximately equal numbers of unmasked pixels\n"
        else:
            header += "# The bins all have an equal width in two theta.\n"
        if 'moveby' in locals() and moveby != 0:
            header += f"# Data has been adjusted. Subtraction of {moveby} from intrgated values.\n"
        header += "#\n#\n"
            
        try:
            header += new_data.detector_description()
            header += "\n#"
        except:
            header += ""
        header += '\n# ' + "2th_deg" + '\t I'

        with open(filename, "w") as xy:
            xy.write(header)
            xy.write("\n")
            for y in range(len(p)):
                xy.write((f"{p[y]}, {i[y]}\n"))
