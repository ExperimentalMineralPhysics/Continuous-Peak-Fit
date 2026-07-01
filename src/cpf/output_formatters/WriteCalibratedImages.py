__all__ = ["Requirements", "WriteOutput"]

import json
import numpy as np
import os
import proglog
from  cpf.settings import get_settings
from cpf.util.io import make_outfile_name
from cpf.util.logging import get_logger
from cpf.histograms import histogram2d, calculate_n_bins
from fabio.TiffIO import TiffIO
from cpf.util.io import numpy_to_json

logger = get_logger("cpf.output_types.WriteCalibratedImages")

def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        "tth_bins": "auto",
        "azimuth_bins": "auto",
    }

    return RequiredParams, OptionalParams


def WriteOutput(settings, debug=False, **kwargs):
    """
    Writes a *.mov file of raw data.

    N.B. this output requires the data files to be present to work.

    Parameters
    ----------
    settings : [str | Path | dict | Settings()]
        Class containing all variables and options needed for the fitting, or 
        dictionary of all the settings or 
        string or path to a file with the settings in.
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

    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters
    tth_bins     = settings_class.output_settings.get("tth_bins", Requirements()[1]["tth_bins"])
    azimuth_bins = settings_class.output_settings.get("azimuth_bins", Requirements()[1]["azimuth_bins"])
    #override with kwargs
    tth_bins        = kwargs.get("tth_bins", tth_bins)
    azimuth_bins = kwargs.get("azimuth_bins", azimuth_bins)

    # make the data class.
    data_to_fill = settings_class.image_list[0]
    data_class = settings_class.data_class
    data_class.fill_data(
        data_to_fill,
        settings=settings_class,
        debug=debug,
    )

    # parse auto bin sizes
    if azimuth_bins == "auto":
        azimuth_bins = 360
    if tth_bins == "auto":
        tth_bins = 1000#calculate_n_bins(data_class.tth, num_azi_bins=azimuth_bins)

    progress = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger
    for z in progress.iter_bar(image=range(settings_class.image_number)):
    # for z in range(settings_class.image_number):
        # read data file
        data_class.import_image(settings_class.image_list[z])
        ortho_data = histogram2d(data_class.intensity, 
                                 data_class.tth, 
                                 data_class.azm, 
                                 x_bins=tth_bins, 
                                 y_bins=azimuth_bins,
                                 azm_bounds = [data_class.azm_start, data_class.azm_end])
        
        intensities = np.nan_to_num(ortho_data[0].T, copy=True, nan=-1, posinf=-1, neginf=-1)

        settings_class.set_subpattern(z, 0)
        
        # make the base file name
        base = os.path.splitext(os.path.split( settings_class.subfit_filename )[1])[0]
        
        out_file = make_outfile_name(
            base,
            directory=settings_class.output_directory,
            extension="tif",
            overwrite=True,
            additional_text="rebinned"
        )
        logger.info(" ".join(map(str, [("Writing %s" % out_file)])))
        tif = TiffIO(out_file, mode='w')
        tif.writeImage(intensities)
        
        
    # make calibration parameters and save
    tth_centers = (ortho_data[1][:-1] + ortho_data[1][1:]) / 2
    azm_centers = (ortho_data[2][:-1] + ortho_data[2][1:]) / 2
    
    x_0 = tth_centers[0]
    x_slope = (tth_centers[-1]-tth_centers[0])/(len(tth_centers)-1)
    y_0 = azm_centers[0]
    y_slope = (azm_centers[-1]-azm_centers[0])/(len(azm_centers)-1)
    
    out = {}
    out['calibration'] = {"x_dim": 0, 
     "x": [x_0, x_slope], 
     "x_start": np.nanmin(ortho_data[1]), 
     "x_end": np.nanmax(ortho_data[1]), 
     "x_label": r"2theta",  
     "x_unit": r"$^\circ$", 
     "y":[y_0, y_slope], 
     "y_start": np.nanmin(ortho_data[2]), 
     "y_end": np.nanmax(ortho_data[2]), 
     "y_label": "Azimuth",
     "y_unit": r"$^\circ$", 
     "max_shape": ortho_data[0].shape}
    out['tth_bin_boundaries'] = ortho_data[1]
    out['azimuth_bin_boundaries'] = ortho_data[2]
        
    calibs_out_file = make_outfile_name(
        settings_class.subfit_filename,
        directory=settings_class.output_directory,
        extension="json",
        overwrite=True,
        additional_text="calibration"
    )
    with open(calibs_out_file, "w") as fout:
        # Write a JSON string into the file.
        json.dump(
            out,
            fout,
            sort_keys=True,
            indent=2,
            default=numpy_to_json,
        )