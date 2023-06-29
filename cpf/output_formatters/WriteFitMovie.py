__all__ = ["Requirements", "WriteOutput"]


import os
import json
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage

import cpf.IO_functions as IO
from cpf.XRD_FitSubpattern import plot_FitAndModel
from cpf.Cosmics import image_preprocess as cosmicsimage_preprocess
from cpf.BrightSpots import SpotProcess


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


def WriteOutput(setting_class=None, setting_file=None, debug=False, **kwargs):
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
        raise ValueError(
            "The frames per second needs to be a number."
        )
        
    

    if setting_class is None and setting_file is None:
        raise ValueError(
            "Either the settings file or the setting class need to be specified."
        )
    elif setting_class is None:
        import cpf.XRD_FitPattern.initiate as initiate

        setting_class = initiate(setting_file)

    #make the base file name
    base = setting_class.datafile_basename
    if base is None or len(base) == 0:
        print("No base filename, trying ending without extension instead.")
        base = setting_class.datafile_ending
    if base is None:
        print("No base filename, using input filename instead.")
        base = os.path.splitext(os.path.split(setting_class.settings_file)[1])[0]
        
    # make the data class.     
    data_to_fill = setting_class.image_list[0]
    data_class = setting_class.data_class
    data_class.fill_data(
        data_to_fill,
        settings=setting_class,
        debug=debug,
    )
    
    
    duration = (setting_class.image_number)/fps
    num_subpatterns = len(setting_class.fit_orders)
    
    for z in range(num_subpatterns):
        
        y = list(range(setting_class.image_number))
        
        setting_class.set_subpattern(0, z)


        addd = IO.peak_string(setting_class.subfit_orders, fname=True)
        if setting_class.file_label != None:
            addd = addd + setting_class.file_label
        out_file = IO.make_outfile_name(
            base, 
            directory=setting_class.output_directory, 
            extension=file_types[0], 
            overwrite=True, 
            additional_text=addd
        )
                    
        # this calls all the iamges and adds them as frames to the video.
        # edited after :https://zulko.github.io/moviepy/getting_started/working_with_matplotlib.html?highlight=matplotlib
        # 4th April 2023.
        def make_frame(t):
        
            # t scales between 0 and 1. 
            # to call each of the images in turn t has to be scaled back 
            # into the number of images (here 'y'). And it has to be an integer. 
            #print(t, int(t*fps), y[int(t*fps)])
                        
            # Get diffraction pattern to process.
            data_class.import_image(setting_class.image_list[y[int(t*fps)]], debug=debug)
                      
            
            if setting_class.datafile_preprocess is not None:
                # needed because image preprocessing adds to the mask and is different for each image.
                data_class.mask_restore()
                if "cosmics" in setting_class.datafile_preprocess:
                    pass#data_class = cosmicsimage_preprocess(data_class, setting_class)
            else:
               # nothing is done here.
               pass
                        
            # restrict data to the right part.           
            sub_data = data_class.duplicate()
            setting_class.set_subpattern(y[int(t*fps)], z)
            sub_data.set_limits(range_bounds=setting_class.subfit_orders["range"])

            # Mask the subpattern by intensity if called for
            if (
                "imax" in setting_class.subfit_orders
                or "imin" in setting_class.subfit_orders
            ):
                sub_data = SpotProcess(sub_data, setting_class)
            
            # read fit file
            json_file = IO.make_outfile_name(
                setting_class.subfit_filename,  # diff_files[z],
                directory=setting_class.output_directory,
                extension=".json",
                overwrite=True,
            )
            with open(json_file) as json_data:
                data_fit = json.load(json_data)[z]
            
            # make the plot of the fits. 
            fig = plt.figure(1)
            fig = plot_FitAndModel(setting_class, 
                                   sub_data, 
                                   #param_lmfit=None, 
                                   params_dict=data_fit, 
                                   figure=fig)
            title_str = (IO.peak_string(setting_class.subfit_orders)+"; " + str(y[int(t*fps)]) +"/" + 
                            str(setting_class.image_number) + "; " +
                            IO.title_file_names(setting_class, num=y[int(t*fps)], image_name=setting_class.subfit_filename)
                            )
            if "note" in setting_class.subfit_orders:
                title_str = title_str + " " + setting_class.subfit_orders["note"]
            plt.suptitle(title_str)   

            #return the figure            
            return mplfig_to_npimage(fig)
        
        # make the video clip
        animation = VideoClip(make_frame, duration=duration)
        animation.write_videofile(out_file, fps=fps)