#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
About
=====
Provides (or will provide) a number of options for smmoothing and noise reduction 
of the images before processing.

This applies filters to the whole image
For example:
    
    data_prepare = {"smooth": {"Gaussian": 2},
                    "cosmics": {},
                    "mask": Calib_mask,
                    "order": {"smooth", "mask"}}


Current useage of intensity limits:    
    "cosmics": {**see cosmics for structure**}
    "smooth": {"smooth": {"filter": "median", "kernel": 10}}
 
Returns: data_class


Implemented:
    - cosmic removal (via another script).
    - median filtering.
    - rolling ball background removal
        
Possible:
    - differernce of gaussians x-ray diffraction


Simon A. Hunt, 2023 - 2024
"""

from skimage import filters, morphology, restoration
from cpf.Cosmics import cosmicsimage
import cpf.settings as settings
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from cpf.XRD_FitPattern import logger

def image_adjust(data, image_process):
    """
    Performs image cleaning processes on diffraction images.
    

    Parameters
    ----------
    data : data class
        Diffraction data to be cleaned.
    image_process : String, list, dictionary or settings class
        String is cleaning process to be performed; list of processes or 
        dictionary containing settings for the cleaning. The setting class must contain 
        the property 'data_prepare'. 
        Examples of licit values are below.

    Raises
    ------
    ValueError
        The type of image processing requested is not know.

    Returns
    -------
    data : data class
        data class containing the cleaned data.
        
    Examples of acceptable inputs:
    image_process = 'cosmics'
    image_process = ['cosmics', 'smooth', 'background']
    image_process = {"cosmics": {"gain": 2.5, "sigclip": 3, "objlim": 3, "sigfrac": 0.3},
                     "smooth": {"filter": "Gaussian", "kernel": 4},
                     "background": {"kernel": 300}
                     "order": ["cosmics", "smooth"]}
                     }
    image_process.data_prepare = ...
    
    If image_process is a dictionary, each sub-dictionary contains the settings for the data cleaning function.
    See each function for full details of the possible variables. 
    
    The image cleaning is performed in the order of the list or dictionary, and in the order of the entires unless
    the dictionary contains the option "order".
    
    The data cleaning only works on image plate data. An error will be created if 
    it is called for non-imageplate X-ray data. 
    
    """

    # issue error if the data type is not correct.
    if data.DispersionType != "AngleDispersive":
        err_str = f"The data cleaning only works on image plate data. This data is {data.DispersionType}"
        raise ValueError(err_str)        

    
    if isinstance(image_process, dict):
        if "order" in image_process:
            processes = image_process["order"]
        else:
            processes = list(image_process.keys())
    elif isinstance(image_process, str):
        processes = [image_process]
    else:
        processes = image_process
        
    for i in range(len(processes)):
        if isinstance(image_process, dict) and processes[i] in image_process:
            options =  image_process[processes[i]]
        else: 
            options=None           

        if processes[i].lower() == "cosmics":
            data = remove_cosmics(data, options=options) 
        elif processes[i].lower() == "smooth":
            data = smooth_image(data, options=options)
        elif processes[i].lower() == "background":
            data = rolling_ball_background(data, options=options)
        elif processes[i].lower() == "scale":
            data = scale_by_background(data, options=options)
        else:
            err_str = "Unknown image pre-process type: %s", processes[i]
            raise ValueError(err_str)
            
    return data
        
        

def remove_cosmics(data, options=None, settings_for_fit=None):
    """
    function to remove the cosmics (very bright spots).

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    settings_for_fit : TYPE
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.

    """
    
    # set defaults
    gain = 2.2
    sigclip = 3.0 # 3 is dioptas default
    objlim = 3.0  # 3 is dioptas default
    sigfrac = 0.3
    params = {"gain": gain, "sigclip": sigclip, "objlim": objlim, "sigfrac": sigfrac}
    if options != None:
        if "gain" in options:
            params["gain"] = options["gain"]
        if "sigclip" in options:
            params["sigclip"] = options["sigclip"]
        if "objlim" in options:
            params["objlim"] = options["objlim"]
        if "sigfrac" in options:
            params["sigfrac"] = options["sigfrac"]
        # FIX ME: use argparse or someway of passing any argument into cosmics.  
    elif settings_for_fit and "cosmics" in settings_for_fit.data_prepare.keys():
        if "gain" in settings_for_fit.data_prepare["cosmics"]:
            params["gain"] = settings_for_fit.data_prepare["cosmics"]["gain"]
        if "sigclip" in settings_for_fit.data_prepare["cosmics"]:
            params["sigclip"] = settings_for_fit.data_prepare["cosmics"]["sigclip"]
        if "objlim" in settings_for_fit.data_prepare["cosmics"]:
            params["objlim"] = settings_for_fit.data_prepare["cosmics"]["objlim"]
        if "sigfrac" in settings_for_fit.data_prepare["cosmics"]:
            params["sigfrac"] = settings_for_fit.data_prepare["cosmics"]["sigfrac"]
        # FIX ME: use argparse or someway of passing any argument into cosmics.
    else: # options == None and settings_for_fit==None:
        pass
    
    test = cosmicsimage(
        data.intensity,
        **params,
    )
    num = 2
    for i in range(num):
        test.lacosmiciteration(True)
        test.clean()
        msk = np.logical_or(
            data.intensity.mask, np.array(test.mask, dtype="bool")
        )

    # apply mask
    data.set_mask(mask=msk)

    return data




def smooth_image(data, options=None, settings_for_fit=None):
    """
    Smooth  the diffraction image with either a median or Gaussian filter. 

    Parameters
    ----------
    data : TYPE
        Diffracion image intensity.
    settings_for_fit : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    data : as imput
        Smoothed diffracion image intensity.

    """
    if isinstance(options, dict):
        if "smooth" in options.keys():
            prep = options
        else:
            prep={}
            prep["smooth"] = options
    elif isinstance(settings_for_fit, settings.settings):
        prep = settings_for_fit.data_prepare
    elif isinstance(settings_for_fit, dict):
        prep = settings_for_fit
    else:
        prep = {"smooth": {"filter": "median", "kernel": 5}}
    
    dt = data.intensity
    
    if prep['smooth']['filter'].lower() == "median":
        med_filt = morphology.disk(prep['smooth']['kernel'])
        #dt = filters.median(dt, med_filt)
        dt = sp.ndimage.median_filter(dt, footprint=med_filt)
    elif prep['smooth']['filter'].lower() == "gaussian":
        sigma = prep['smooth']['kernel']
        dt = filters.gaussian(data.intensity, sigma)
    elif isinstance(prep['smooth']['filter'], str):#.lower() == "nanmedian":
        if "nan" in prep['smooth']['filter']:
            dt.data[dt.mask==True] = np.nan
        if not "kernel" in prep['smooth']:
            prep['smooth']['kernel'] = 5
        try:
            dt = sp.ndimage.generic_filter(dt, getattr(np, prep['smooth']['filter']), [prep['smooth']['kernel'],prep['smooth']['kernel']])
        except:
            err_str = (
                 "'" +prep['smooth']['filter'] + "' is not a numpy function that can be used as a smoothing function."
             )
            raise ValueError(err_str)
        
        # if "nan" in prep['smooth']['filter']:
        #     dt.data[dt.mask==True] = np.nan
        #     print(dt.data[0])       
    else:
        err_str = (
             "Unknown image smoothing type."
         )
        raise ValueError(err_str)
       
    # copy data back into data class
    msk = data.intensity.mask
    dt = np.ma.array(dt, mask=msk)
    data.intensity = dt
            
    return data




def rolling_ball_background(data, options=None, settings_for_fit=None):
    """
    rolling ball background removal is discussed in He's text book on 2D XRD.   

    It is recomended that the image is smoothed before rolling ball is run inorder to 
    negate the effects of any dark pixels in the diffraction image. 

    Parameters
    ----------
    data : data class
        Diffraction data to be cleaned.
    options : dictionary, optional
        Contains options for background removal. 
    settings_for_fit : settings class, optional
        Contains options for background removal. 
        
    Returns
    -------
    data : data_class
        Data class with the updated data.

    """
    
    if isinstance(options, dict):
        if "background" in options.keys():
            prep = options
        else:
            prep={}
            prep["background"] = options
    elif isinstance(settings_for_fit, settings.settings):
        prep = settings_for_fit.data_prepare
    elif isinstance(settings_for_fit, dict):
        prep = settings_for_fit
    else:
        prep = {"background": {"kernel": 50}}    
        
    original_data = data.duplicate()
    
    if "smooth" in prep["background"]:
        # FIXME: we need a median filter that accounts for the mask. 
        # This is either a filter than dignors the mask or we make all the masked values nan and run a nanmedial filter. 
        if prep["background"]["smooth"] == "":
            prep["background"]["smooth"] = {"filter": "median", "kernel": 5}
        data = smooth_image(data, prep["background"]["smooth"])
    elif "smooth" in prep:
        if prep["smooth"] == "":
            prep["smooth"] = {"filter": "median", "kernel": 5}        
        data = smooth_image(data, prep["smooth"])
        
    dt = data.intensity
    dt.data[dt.mask==True] = np.max(dt)
    #FIX ME: the removal here is a metian for no good reason. I am trying to make sure the maksed areas do not bleed back into the iamge proper 
    # I do not know if it is necessary though.
    
    def plot_result(image, background, sigma=None ):
        fig, ax = plt.subplots(nrows=1, ncols=3)
    
        a = ax[0].imshow(image, cmap='jet')
        ax[0].set_title('Original image')
        ax[0].axis('off')
        #plt.colorbar(a,cax=ax[0])
        #ax[1].colorbar()
        fig.colorbar(a, ax=ax[0], orientation='horizontal')
    
        b=ax[1].imshow(background, cmap='jet')
        ax[1].set_title(r'Background, radius=%i' % sigma)
        ax[1].axis('off')
        #plt.colorbar(b,cax=ax[1])
        #ax[1].colorbar()
        fig.colorbar(b, ax=ax[1], orientation='horizontal')
    
        c=ax[2].imshow(image - background, cmap='jet')
        ax[2].set_title('Result')
        ax[2].axis('off')
        #plt.colorbar(c,cax=ax[2])
        #fig.colorbar#add_colorbar(a)
    
        fig.colorbar(c, ax=ax[2], orientation='horizontal')
        fig.tight_layout()
    
    
    background = restoration.rolling_ball(dt, radius=prep['background']['kernel'])
    
    plot_result(original_data.intensity, background, prep['background']['kernel'])
    plt.show()

    no_bg = original_data.intensity - background
    
    # copy data back into data class
    msk = no_bg.mask
    no_bg = np.ma.array(no_bg, mask=msk)
    data.intensity = no_bg
        
    return data

    
    
def scale_by_background(data, options=None, settings_for_fit=None):
    """
    Scales the intsinsities to that of a rolling ball background performed on a smoothed image. 
    
    
    # assume that the backgound is proportional to the transmissivity of the anvils. 
    # Scale the intensities by the background and multiply back up.
    

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    options : TYPE, optional
        DESCRIPTION. The default is None.
    settings_for_fit : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    
    if isinstance(options, dict):
        if "background" in options.keys():
            prep = options
        else:
            prep={}
            prep["background"] = options
    elif isinstance(settings_for_fit, settings.settings):
        prep = settings_for_fit.data_prepare
    elif isinstance(settings_for_fit, dict):
        prep = settings_for_fit
    else:
        prep = {"background": {"kernel": 200, "smooth": {"filter": "median", "kernel": 5}}} 

    original_data = data.duplicate()
    
    background = rolling_ball_background(data, prep)
    
    no_bg = original_data.intensity / (original_data.intensity-background.intensity)
    
    # copy data back into data class
    msk = no_bg.mask
    no_bg = np.ma.array(no_bg, mask=msk)
    data.intensity = no_bg
        
    return data