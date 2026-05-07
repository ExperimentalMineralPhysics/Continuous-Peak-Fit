__all__ = ["Requirements", "WriteOutput"]


import os
from copy import deepcopy
import proglog
import matplotlib.pyplot as plt
import numpy as np
from moviepy.video.VideoClip import VideoClip
# from typing import Literal, Optional
from moviepy import ImageClip
# from moviepy import concatenate
from moviepy import VideoFileClip, concatenate_videoclips
from textwrap import wrap

# import cpf.IO_functions as IO
from  cpf.settings import get_settings
from cpf.IO_functions import make_outfile_name, title_file_names
from cpf.util.logging import get_logger
from cpf.util.output_formatters import mplfig_to_npimage

logger = get_logger("cpf.output_types.WriteCollectionMovie")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = {
        "fps": 10,  # frames per second
        "file_types": ["mp4"],  # movie file type
        "Irange": ["pt1percentile", "99pt9percentile"], # range of colour scale
        "plot data as": "calibrated", 
        "plot type": "default", #"surface",
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

    # FIXME: make so that it can iterate over each range and make a movie of each selected range.

    # make sure settings is a class
    settings_class = get_settings(settings)

    # Parse optional parameters
    fps        = settings_class.output_settings.get("fps", Requirements()[1]["fps"])
    file_types = settings_class.output_settings.get("file_types", Requirements()[1]["file_types"])
    Irange     = settings_class.output_settings.get("Irange", Requirements()[1]["Irange"])
    plot_data_as = settings_class.output_settings.get("plot data as", Requirements()[1]["plot data as"])
    plot_type  = settings_class.output_settings.get("plot type", Requirements()[1]["plot type"])
    #override with kwargs
    fps        = kwargs.get("fps", fps)
    file_types = kwargs.get("file_types", file_types)
    Irange     = kwargs.get("Irange", Irange)
    plot_data_as = kwargs.get("plot_data_as", plot_data_as)
    plot_type  = kwargs.get("plot_type", plot_type)

    # make sure file_types is a list.
    if isinstance(file_types, str):
        file_types = [file_types]
    if not isinstance(fps, float) and not isinstance(fps, int):
        raise ValueError("The frames per second needs to be a number.")
    if plot_data_as != "collected" and plot_data_as != "calibrated":
        raise ValueError("plot_type must be 'collected' or 'calibrated'.")

    # make the base file name
    if settings_class:
        base = settings_class.datafile_basename
    else:
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    if base is None or len(base) == 0:
        logger.info("No base filename, trying ending without extension instead.")
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

    # get the intensity ranges from the data.
    Ipctl = []
    Imin = []

    # to plot maximum inentsity set prctl=100
    prctl = 99.9

    progress = proglog.default_bar_logger("bar")  # shorthand to generate a bar logger
    print("Reading images to get intensity range")
    for z in progress.iter_bar(image=range(settings_class.image_number)):
    # for z in range(settings_class.image_number):
        # read data file
        settings_class.set_subpattern(z, 0)
        data_class.import_image(settings=settings_class)
        Ipctl.append(
            np.nanpercentile(
                np.ma.filled(data_class.intensity, np.nan), prctl
            )
        )
        Imin.append(np.min(data_class.intensity))

    lims = {
        "max": np.max(Ipctl),
        "min": np.min(Imin),
        "rmax": np.max(Ipctl),
        "rmin": np.min(Imin),
    }

    duration = (settings_class.image_number) / fps

    y = list(range(settings_class.image_number))

    fig = plt.figure(figsize=(6, 8))
    ax = fig.add_subplot(1, 1, 1)

    # this calls all the iamges and adds them as frames to the video.
    # edited after :https://zulko.github.io/moviepy/getting_started/working_with_matplotlib.html?highlight=matplotlib
    # 4th April 2023.
    def make_frame(t):
        # t scales between 0 and 1.
        # to call each of the images in turn t has to be scaled back
        # into the number of images (here 'y'). And it has to be an integer.
        # logger.info(" ".join(map(str, [(t, int(t*fps), y[int(t*fps)])])))

        # Get diffraction pattern to process.
        settings_class.set_subpattern(y[int(t * fps)], 0)
        data_class.import_image(settings=settings_class)
        # data_class.import_image(settings_class.image_list[y[int(t * fps)]])

        if settings_class.datafile_preprocess is not None:
            # needed because image preprocessing adds to the mask and is different for each image.
            data_class.mask_restore()
            if "cosmics" in settings_class.datafile_preprocess:
                pass  # data_class = cosmicsimage_preprocess(data_class, settings_class)
        else:
            # nothing is done here.
            pass
        if t==0 and isinstance(t, int):
            # the first time the this function is called by VideoClip t is an integer.
            # everyother time it is a float.
            # use this to determine whether to make the colour bar or not
            cbar = None
        else:
            cbar = False
            
        ax.clear()
        if plot_data_as == "calibrated":
            data_class.plot_calibrated(
                fig_plot=fig, 
                axis_plot=ax, 
                show="intensity", 
                limits=deepcopy(lims),
                plot_type = plot_type,
                cbar_axes=cbar
            )
        else:
            data_class.plot_collected(
                fig_plot=fig, axis_plot=ax, show="intensity", limits=deepcopy(lims),
                cbar_axes=cbar
            )
        ax.set_title("\n".join(wrap(title_file_names(settings_for_fit=settings_class, num=int(t * fps)), 60)))

        # return the figure
        return mplfig_to_npimage(fig)

    # make the video clip
    animation = VideoClip(make_frame, duration=duration)
    for f in range(len(file_types)):
        out_file = make_outfile_name(
            base,
            directory=settings_class.output_directory,
            extension=file_types[f],
            overwrite=True,
        )
        logger.info(" ".join(map(str, [("Writing %s" % out_file)])))
        animation.write_videofile(out_file, fps=fps)
    animation.close()
    
