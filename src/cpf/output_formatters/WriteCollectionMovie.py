__all__ = ["Requirements", "WriteOutput"]


import os
import numpy as np
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
from copy import deepcopy

# import cpf.IO_functions as IO
from cpf.IO_functions import make_outfile_name, title_file_names
from cpf.logging import CPFLogger

logger = CPFLogger("cpf.output_types.WriteCollectionMovie")


def Requirements():
    # List non-universally required parameters for writing this output type.

    RequiredParams = [
        #'apparently none!
    ]
    OptionalParams = [
        "fps",
        "file_types"
    ]

    return RequiredParams, OptionalParams


def WriteOutput(settings_class=None, setting_file=None, debug=False, **kwargs):
    """
    Writes a *.mov file of raw data.
    
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

    # FIXME: make so that it can iterate over each range and make a movie of each selected range.


    if not "file_types" in kwargs:
        file_types = ".mp4"
    # make sure file_types is a list.
    if isinstance(file_types, str):
        file_types = [file_types]
    if not "fps" in kwargs:
        fps = 10
    elif not isinstance(fps, float):
        raise ValueError(
            "The frames per second needs to be a number."
        )
        
    

    if settings_class is None and setting_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif settings_class is None:
        import cpf.XRD_FitPattern.initiate as initiate

        settings_class = initiate(setting_file)

    #make the base file name
    if setting_file:
        base = os.path.splitext(os.path.split(settings_class.settings_file)[1])[0]
    else:
        base = settings_class.datafile_basename
    if base is None or len(base) == 0:
        logger.info(" ".join(map(str, [("No base filename, trying ending without extension instead.")])))
        base = settings_class.datafile_ending
    if base is None:
        logger.info(" ".join(map(str, [("No base filename, using input filename instead.")])))
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
    Ipctl=[]
    Imin=[]
    
    # to plot maximum inentsity set prctl=100
    prctl = 99.9
    
    for z in range(settings_class.image_number):
        # read data file
        data_class.import_image(settings_class.image_list[z])
        Ipctl.append(np.percentile(data_class.intensity[data_class.intensity.mask==False], prctl))
        Imin.append(np.min(data_class.intensity))
        
    lims = {"max": np.max(Ipctl),
            "min": np.min(Imin),
            "rmax": np.max(Ipctl),
            "rmin": np.min(Imin)}
        
    duration = (settings_class.image_number)/fps
    
    y = list(range(settings_class.image_number))
    
    # settings_class.set_subpattern(0, z)


    # addd = IO.peak_string(settings_class.subfit_orders, fname=True)
    # if settings_class.file_label != None:
    #     addd = addd + settings_class.file_label
    out_file = make_outfile_name(
        base, 
        directory=settings_class.output_directory, 
        extension=file_types[0], 
        overwrite=True,
    )
    logger.info(" ".join(map(str, [("Writing %s" % out_file)])))
    
    # this calls all the iamges and adds them as frames to the video.
    # edited after :https://zulko.github.io/moviepy/getting_started/working_with_matplotlib.html?highlight=matplotlib
    # 4th April 2023.
    def make_frame(t):
    
        # t scales between 0 and 1. 
        # to call each of the images in turn t has to be scaled back 
        # into the number of images (here 'y'). And it has to be an integer. 
        #logger.info(" ".join(map(str, [(t, int(t*fps), y[int(t*fps)])])))
                    
        # Get diffraction pattern to process.
        data_class.import_image(settings_class.image_list[y[int(t*fps)]])
                  
        
        if settings_class.datafile_preprocess is not None:
            # needed because image preprocessing adds to the mask and is different for each image.
            data_class.mask_restore()
            if "cosmics" in settings_class.datafile_preprocess:
                pass#data_class = cosmicsimage_preprocess(data_class, settings_class)
        else:
           # nothing is done here.
           pass
            
        fig = plt.figure(figsize=(6, 8))
        ax = fig.add_subplot(1, 1, 1)
        data_class.plot_calibrated(fig_plot=fig, axis_plot=ax, show="intensity", limits=deepcopy(lims))
        plt.title(title_file_names(settings_for_fit=settings_class, num=int(t*fps)))
        
        #return the figure            
        return mplfig_to_npimage(fig)
        
    # make the video clip
    animation = VideoClip(make_frame, duration=duration)
    animation.write_videofile(out_file, fps=fps)
    animation.close()