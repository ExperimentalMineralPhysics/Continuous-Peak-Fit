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


Simon Hunt, March 2023
"""

from skimage import filters, morphology, restoration
#from cpf.Cosmics import image_preprocess as cosmicsimage_preprocess
from cpf.Cosmics import cosmicsimage
import cpf.settings as settings
import numpy as np
import matplotlib.pyplot as plt


def image_adjust(data, settings_for_fit):
                
    
    # image processing routines to prepare the images for fitting. 
    # The order of the processes is delibrate and fixed. 
    # For example: it is sensible to remove the very bright 'cosmics' from the data first before smoothing. 
    
    print(settings_for_fit.data_prepare.items())
    
    for i in range(len(settings_for_fit.data_prepare)):
        
        key = list(settings_for_fit.data_prepare.keys())[i]
        
        if key.lower() == "cosmics":
            data = remove_cosmics(data, options=settings_for_fit.data_prepare["cosmics"])
            
        elif key.lower() == "smooth":
            data = smoothimage(data, settings_for_fit)
            
        elif key.lower() == "background":
            data = rolling_ball_background(data, settings_for_fit)
    
        else:
            err_str = "Unknown image pre-process type: %s", key
            raise ValueError(err_str)
    
        # if "cosmics" in settings_for_fit.data_prepare:
        #     data = cosmicsimage_preprocess(data, settings_for_fit)
            
        # if "smooth" in settings_for_fit.data_prepare:
        #     data = smoothimage(data, settings_for_fit)
            
        # if "background" in settings_for_fit.data_prepare:
        #     data = rolling_ball_background(data, settings_for_fit)
            
    return data
        
        

def remove_cosmics(data_as_class, options=None, settings_for_fit=None):
    """
    function to remove the cosmics (very bright spots).

    Parameters
    ----------
    data_as_class : TYPE
        DESCRIPTION.
    settings_for_fit : TYPE
        DESCRIPTION.

    Returns
    -------
    data_as_class : TYPE
        DESCRIPTION.

    """
    print("Remove Cosmics")
    # set defaults
    gain = 2.2
    sigclip = 1.5  # 3 is dioptas default
    objlim = 3.0  # 3 is dioptas default
    sigfrac = 0.3
    params = {"gain": gain, "sigclip": sigclip, "objlim": objlim, "sigfrac": sigfrac}
    if bool(settings_for_fit.datafile_preprocess["cosmics"]):
        if "gain" in settings_for_fit.datafile_preprocess["cosmics"]:
            params["gain"] = settings_for_fit.datafile_preprocess["cosmics"]["gain"]
        if "sigclip" in settings_for_fit.datafile_preprocess["cosmics"]:
            params["sigclip"] = settings_for_fit.datafile_preprocess["cosmics"]["sigclip"]
        if "objlim" in settings_for_fit.datafile_preprocess["cosmics"]:
            params["objlim"] = settings_for_fit.datafile_preprocess["cosmics"]["objlim"]
        if "sigfrac" in settings_for_fit.datafile_preprocess["cosmics"]:
            params["sigfrac"] = settings_for_fit.datafile_preprocess["cosmics"]["sigfrac"]
        # FIX ME: use argparse or someway of passing any argument into cosmics.
    elif options != None:
        if "gain" in options:
            params["gain"] = options["gain"]
        if "sigclip" in options:
            params["sigclip"] = options["sigclip"]
        if "objlim" in options:
            params["objlim"] = options["objlim"]
        if "sigfrac" in options:
            params["sigfrac"] = options["sigfrac"]
        # FIX ME: use argparse or someway of passing any argument into cosmics.        
    else: # options == None and settings_for_fit==None:
        pass
    
    test = cosmicsimage(
        data_as_class.intensity,
        **params,
    )
    num = 2
    for i in range(num):
        test.lacosmiciteration(True)
        # print(asdf["nnew"], asdf["niter"])
        test.clean()
        msk = np.logical_or(
            data_as_class.intensity.mask, np.array(test.mask, dtype="bool")
        )

    # apply mask
    data_as_class.set_mask(mask=msk)

    return data_as_class




def smoothimage(data, settings_for_fit=None):
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
    
    if isinstance(settings_for_fit, settings.settings):
        prep = settings_for_fit.data_prepare
    elif isinstance(settings_for_fit, dict):
        prep = settings_for_fit
    else:
        prep = {"smooth": {"filter": "median", "kernel": 10}}
        
    if prep['smooth']['filter'].lower() == "median":
       
       print("median")
       
       median_filter = morphology.disk(prep['smooth']['kernel'])
       
       data = filters.median(data, median_filter)
   

    elif prep['smooth']['filter'].lower() == "gaussian":
        print("Gaussian")
        
        # smooth_mean = ndi.correlate(bright_square, mean_kernel)
        
        sigma = prep['smooth']['kernel']
        # smooth = filters.gaussian(bright_square, sigma)
        data = filters.gaussian(data, sigma)
       
    else:
        err_str = (
             "Unknown image smoothing type."
         )
        raise ValueError(err_str)
       
       
    return data




def rolling_ball_background(data, settings_for_fit=None):
    
    # rolling ball background removal is discussed in He's text book on 2D XRD.    
    
    
    # the background removal has to smooth the data before being run inorder to 
    # negate the effects of any dark pixels in the diffraction image. 
    if isinstance(settings_for_fit, settings.settings):
        prep = settings_for_fit.data_prepare
    elif isinstance(settings_for_fit, dict):
        prep = settings_for_fit
    else:
        prep = {"smooth": {"filter": "median", "kernel": 10}}
    
    data = smoothimage(data, prep)
    
    
    # we also need to NaN masked values to repvent bluring of masked parts into the good data. 
    # Not NaN. Make the maseked values big. Then rolling the ball around under the image does not 
    # bleed into the good (unmasked) part of the image. After making the background replace all masked bakgroung values 
    # with either 0, NaN or a mean of the background... 
    # This will bleed the background into the masked areas and not the unmasked parts of the image.
    
    # I don't think we can use NaN. so reaplce all masked values with maximum intensity of image. 
    # see :https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.median
    
    
    def plot_result(image, background, sigma=None ):
        fig, ax = plt.subplots(nrows=1, ncols=3)
    
        a = ax[0].imshow(image, cmap='gray')
        ax[0].set_title('Original image')
        ax[0].axis('off')
        #plt.colorbar(a,cax=ax[0])
        #ax[1].colorbar()
        fig.colorbar(a, ax=ax[0], orientation='horizontal')
    
        b=ax[1].imshow(background, cmap='gray')
        ax[1].set_title(r'Background, radius=%i' % sigma)
        ax[1].axis('off')
        #plt.colorbar(b,cax=ax[1])
        #ax[1].colorbar()
        fig.colorbar(b, ax=ax[1], orientation='horizontal')
    
        c=ax[2].imshow(image - background, cmap='gray')
        ax[2].set_title('Result')
        ax[2].axis('off')
        #plt.colorbar(c,cax=ax[2])
        #fig.colorbar#add_colorbar(a)
    
        fig.colorbar(c, ax=ax[2], orientation='horizontal')
        fig.tight_layout()
    
    if 'kernel' in prep['background']:
        sigma = prep['background']['kernel']
    else:
        sigma = 300
    
    #image = data.coins()
    print("rolling ball")
    background = restoration.rolling_ball(data, radius=sigma, )
    
    plot_result(data, background, sigma)
    plt.show()

    # assume that the backgound is proportional to the transmissivity of the anvils. 
    # Scale the intensities by the background and multiply back up.
    data = data/background
    
    return data