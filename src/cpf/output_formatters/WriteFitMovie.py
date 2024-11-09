__all__ = ["Requirements", "WriteOutput"]


import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

import cpf.IO_functions as IO
from cpf.BrightSpots import SpotProcess
from cpf.data_preprocess import remove_cosmics as cosmicsimage_preprocess

# from cpf.XRD_FitPattern import logger
from cpf.logger_functions import logger
from cpf.XRD_FitSubpattern import plot_FitAndModel


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = ["fps", "file_types"]

    return RequiredParams, OptionalParams


def WriteOutput(settings_class=None, settings_file=None, debug=False, **kwargs):
    """
    Writes a *.?? file of the fits.

    N.B. this output requires the data files to be present to work.

    Parameters
    ----------
    FitSettings : TYPE
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

    if not "file_types" in kwargs:
        file_types = ".mp4"
    # make sure file_types is a list.
    if isinstance(file_types, str):
        file_types = [file_types]
    if not "fps" in kwargs:
        fps = 10
    elif not isinstance(fps, float):
        raise ValueError("The frames per second needs to be a number.")

    if settings_class is None and settings_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        import cpf.XRD_FitPattern.initiate as initiate

        settings_class = initiate(settings_file)

    # make the base file name
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
    data_class = settings_class.data_class
    data_class.fill_data(
        data_to_fill,
        settings=settings_class,
        debug=debug,
    )

    # get the intensity ranges from the data fits.
    data_range = [[] for i in range(len(settings_class.fit_orders))]
    model_range = [[] for i in range(len(settings_class.fit_orders))]
    resid_range = [[] for i in range(len(settings_class.fit_orders))]
    for z in range(settings_class.image_number):
        settings_class.set_subpattern(z, 0)
        # read fit file
        json_file = IO.make_outfile_name(
            settings_class.subfit_filename,  # diff_files[z],
            directory=settings_class.output_directory,
            extension=".json",
            overwrite=True,
        )
        with open(json_file) as json_data:
            data_fit = json.load(json_data)
        for y in range(len(data_fit)):
            data_range[y].append(data_fit[y]["DataProperties"])
            try:
                # try to see if model range is in the json file. If it is not
                # then just fill with DataProperties.
                model_range[y].append(data_fit[y]["ModelProperties"])
            except:
                model_range[y].append(data_fit[y]["DataProperties"])
            try:
                # try to see if model range is in the json file. If it is not
                # then just fill with DataProperties.
                resid_range[y].append(data_fit[y]["ResidualProperties"])
            except:
                resid_range[y].append({"max": np.nan, "min": np.nan})
    Imax = []
    Imin = []
    Rmax = []
    Rmin = []
    for y in range(len(data_fit)):
        tmp1 = pd.DataFrame(data_range[y], index=list(range(len(data_range[y]))))
        tmp2 = pd.DataFrame(model_range[y], index=list(range(len(model_range[y]))))
        Imax.append(np.max([tmp1["max"].max(), tmp2["max"].max()]))
        Imin.append(np.min([tmp1["min"].min(), tmp2["min"].min()]))
        tmp3 = pd.DataFrame(resid_range[y], index=list(range(len(resid_range[y]))))
        Rmax.append(tmp3["max"].max())
        Rmin.append(tmp3["min"].min())

    duration = (settings_class.image_number) / fps

    for z in range(len(settings_class.fit_orders)):
        y = list(range(settings_class.image_number))

        settings_class.set_subpattern(0, z)

        addd = IO.peak_string(settings_class.subfit_orders, fname=True)
        if settings_class.file_label != None:
            addd = addd + settings_class.file_label
        out_file = IO.make_outfile_name(
            base,
            directory=settings_class.output_directory,
            extension=file_types[0],
            overwrite=True,
            additional_text=addd,
        )
        logger.info(" ".join(map(str, [("Writing %s" % out_file)])))

        # this calls all the iamges and adds them as frames to the video.
        # edited after :https://zulko.github.io/moviepy/getting_started/working_with_matplotlib.html?highlight=matplotlib
        # 4th April 2023.
        def make_frame(t):
            # t scales between 0 and 1.
            # to call each of the images in turn t has to be scaled back
            # into the number of images (here 'y'). And it has to be an integer.
            # logger.info(" ".join(map(str, [(t, int(t*fps), y[int(t*fps)])])))

            # Get diffraction pattern to process.
            data_class.import_image(settings_class.image_list[y[int(t * fps)]])

            if settings_class.datafile_preprocess is not None:
                # needed because image preprocessing adds to the mask and is different for each image.
                data_class.mask_restore()
                if "cosmics" in settings_class.datafile_preprocess:
                    pass  # data_class = cosmicsimage_preprocess(data_class, settings_class)
            else:
                # nothing is done here.
                pass

            # restrict data to the right part.
            sub_data = data_class.duplicate()
            settings_class.set_subpattern(y[int(t * fps)], z)
            sub_data.set_limits(range_bounds=settings_class.subfit_orders["range"])

            # Mask the subpattern by intensity if called for
            if (
                "imax" in settings_class.subfit_orders
                or "imin" in settings_class.subfit_orders
            ):
                sub_data = SpotProcess(sub_data, settings_class)

            # read fit file
            json_file = IO.make_outfile_name(
                settings_class.subfit_filename,  # diff_files[z],
                directory=settings_class.output_directory,
                extension=".json",
                overwrite=True,
            )
            with open(json_file) as json_data:
                data_fit = json.load(json_data)[z]

            # make the plot of the fits.
            fig = plt.figure(1)
            fig = plot_FitAndModel(
                settings_class,
                sub_data,
                # param_lmfit=None,
                params_dict=data_fit,
                figure=fig,
                plot_ColourRange={
                    "max": Imax[z],
                    "min": Imin[z],
                    "rmin": Rmin[z],
                    "rmax": Rmax[z],
                },
            )
            title_str = (
                IO.peak_string(settings_class.subfit_orders)
                + "; "
                + str(y[int(t * fps)])
                + "/"
                + str(settings_class.image_number)
                + "\n"
                + IO.title_file_names(
                    settings_class,
                    num=y[int(t * fps)],
                    image_name=settings_class.subfit_filename,
                )
            )
            if "note" in settings_class.subfit_orders:
                title_str = title_str + " " + settings_class.subfit_orders["note"]
            plt.suptitle(title_str)
            IO.figure_suptitle_space(fig, topmargin=0.4)

            # return the figure
            return mplfig_to_npimage(fig)

        # make the video clip
        animation = VideoClip(make_frame, duration=duration)
        animation.write_videofile(out_file, fps=fps)
        animation.close()
