__all__ = ["Requirements", "WriteOutput"]

import os

import numpy as np

from cpf.output_formatters.convert_fit_to_unitcell import fits_to_unitcell
from cpf.output_formatters.output_csv import make_header, write_csv
from cpf.settings import get_settings
from cpf.util.io import make_outfile_name
from cpf.util.logging import get_logger

logger = get_logger("cpf.output_formatters.WriteCoefficientTable")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        ##"Output_directory"  # if no direcrtory is specified write to current directory.
        "reflections_to_use": "all",  # -- pick which set of reflections to use for unit cell volume
        "phase": True,  # -- pick whick phase to fit unit cell for.
        "SampleGeometry": "3d",  # -- geometry of the sample for determining the cnetres from. 2D or 3D.
        "SampleDeformation": "compression",  # changes calculation between 'compression' and 'extension'.
        "weighted": True,  # -- weighted fit or not. True/False
        "dp": 5,  # how many decimal points to write out
        "col_width": 15,  # default column width for csv file.
        "ordering_of_output": None,  # Just leave as read -- otherwise list of dataframe headers to order by
    }
    # OptionalParams = [
    #     ##"Output_directory"  # if no direcrtory is specified write to current directory.
    #     "reflections_to_use"  # -- pick which set of reflections to use for unit cell volume
    #     "phase" # -- pick whick phase to fit unit cell for.
    #     "SampleGeometry" # -- geometry of the sample for determining the cnetres from. 2D or 3D.
    #     "weighted" # -- weighted fit or not. True/False
    # ]

    return RequiredParams, OptionalParams


def WriteOutput(
    settings,
    *args,
    **kwargs,
):
    """
    Write unit-cell volumes derived from fitted peak centroids. Writes the values
    to table/csv file.

    Parameters
    ----------
    settings : cpf.settings.Settings() class
        input file or settings class used to make the fits.
    *args : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    """

    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters
    reflections_to_use = settings_class.output_settings.get(
        "reflections_to_use", Requirements()[1]["reflections_to_use"]
    )  # -- pick which set of reflections to use for unit cell volume
    phase = settings_class.output_settings.get("phase", Requirements()[1]["phase"])
    SampleGeometry = settings_class.output_settings.get(
        "SampleGeometry", Requirements()[1]["SampleGeometry"]
    )
    SampleDeformation = settings_class.output_settings.get(
        "SampleDeformation", Requirements()[1]["SampleDeformation"]
    )
    weighted = settings_class.output_settings.get(
        "weighted", Requirements()[1]["weighted"]
    )
    dp = settings_class.output_settings.get("dp", Requirements()[1]["dp"])
    col_width = settings_class.output_settings.get(
        "col_width", Requirements()[1]["col_width"]
    )
    ordering_of_output = settings_class.output_settings.get(
        "ordering_of_output", Requirements()[1]["ordering_of_output"]
    )
    # override with kwargs
    reflections_to_use = kwargs.get("reflections_to_use", reflections_to_use)
    phase = kwargs.get("phase", phase)
    SampleGeometry = kwargs.get("SampleGeometry", SampleGeometry)
    SampleDeformation = kwargs.get("SampleDeformation", SampleDeformation)
    weighted = kwargs.get("weighted", weighted)
    dp = kwargs.get("dp", dp)
    col_width = kwargs.get("col_width", col_width)
    ordering_of_output = kwargs.get("ordering_of_output", ordering_of_output)

    # force all the kwargs that might be needed
    set_params = {  # "temperature": np.nan,
        "reflections_to_use": reflections_to_use,
        "phase": phase,
        "SampleGeometry": SampleGeometry,
        "weighted": weighted,
        "pressure": False,  # to supress pressure, and stresses
    }
    kwargs.update(set_params)

    # get the unit cells from settings.
    df = fits_to_unitcell(settings, **kwargs)
    headers = list(df.columns.values)
    # order the rows
    if ordering_of_output:
        df = df.sort_values(by=ordering_of_output)

    # make filename for output
    base = settings_class.datafile_basename
    if base is None:
        logger.info(
            " ".join(map(str, [("No base filename, using input filename instead.")]))
        )
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    out_file = make_outfile_name(
        base,
        directory=settings_class.output_directory,
        extension=".dat",
        overwrite=True,
        additional_text="unit_cells",
    )

    ## outfile header
    calc_options = {}
    calc_options["Sample Geometry"] = SampleGeometry
    calc_options["Sample Deformation"] = SampleDeformation
    # remove hkls from data frame -- write as a header instead
    cols = [col for col in df.columns if "hkl" in col]
    hkls = {}
    for i in cols:
        hkls[i] = df[i].iloc[0]
        df = df.drop(i, axis=1)
    if hkls:
        key = list(hkls.keys())
        if len(key) != 1:
            # do something with the 'phase' to find the right key
            key = key[0]
        else:
            key = key[0]
        calc_options["peaks used"] = "{" + "} {".join(v for v in hkls[key]) + "}"
    file_header = make_header(
        settings_class, derived="Unit Cells", calc_options=calc_options
    )

    # write file using panda dataframe
    write_csv(out_file, df, headers, file_header, col_width=col_width, dp=dp)
