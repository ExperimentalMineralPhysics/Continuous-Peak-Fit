#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file contains the code to make and apply masks to the diffraction data.

The masks are common to all energy dispereive data types. 

The masks can come in a number for formats.
 - an image (or *.mask) file the same size as the diffraction detector
 - a dictionary of polygons, thresholds and ranges (similar to the funcionality of GSAS-II).
    An example of the new dictionary type mask is: 
    Calib_mask     = {"polygon": [[(0,43),(180,43),(205,103),(127,160),(97,247),(80,393),(0,377)],
                                  [(1000,43),(820,43),(795,103),(853,160),(890,247),(910,393),(1000,377)],
                                  [(305,357),(286,375),(286,375),(305,275)],
                                 ],
                      "threshold": [0.25, 8],
                      "image": "file.mask"}



Need functions that can: 
    a. make mask from settings (store as original mask) --- get_mask
    b. make mask of smaller range (azm and tth masks)   --- mask_outside
    c. apply mask to all arrays                         --- mask_apply
    d. restore original mask                            --- mask_restore
    e. remove all masks                                 --- mask_remove
    

Remove mask functions from these class files and put in separate common file. 
"""


import numpy as np
import numpy.ma as ma
# import matplotlib.pyplot as plt
#from matplotlib import gridspec, cm, colors
from PIL import Image, ImageDraw
from cpf.XRD_FitPattern import logger



class _masks():
    """
    Class definiing mask functions. 
    
    These are imported into the detector functions as methods. 
    """
    # def __init__(self, detector_class=None):
    #     """
    #     :param detector_class:
    #     """
    #     self = detector_class
    
    
    def get_mask(self, mask, im_ints=None, debug=False):
        """
        Creates the mask for the diffracion data and returns a boolian image. 
    
    
        Parameters
        ----------
        mask : dictionary or string
            Object that can created into a mask for the data.
        im_ints : array
            Diffraction data array -- so that can build mask around it.
    
        Returns
        -------
        im_mask : boolian array.
            Mask to be applied to the diffraction data.
    
        """
        
        if im_ints is None:
            im_ints = self.intensity
            
        # if the mask is an image string then wrap it in a dictionary. 
        if not isinstance(mask, dict) and mask is not None:
            mask = {"image": mask}
        elif mask is None:
            #make empty dictionary
            mask = {}
            
        # make empty mask
        im_mask = np.zeros(im_ints.shape, dtype='bool')
    
        if "image" in mask:
            #Dioptas mask is compressed Tiff image.
            #Save and load functions within Dioptas are: load_mask and save_mask in dioptas/model/MaskModel.py
            mask_from_image = np.array(Image.open(mask["image"]))
            #im_ints = ma.array(im_ints, mask=im_mask)    
            #im_ints = ma.masked_less(im_ints, 0)
            if mask_from_image.shape == im_ints.shape:
                im_mask = np.asarray(im_mask) | np.asarray(mask_from_image) 
            elif mask_from_image.shape == im_ints.shape[1:2]:
                # the masked image is for a single detector frame and the diffraction data is composed of multiple frames
                # (e.g. ESRFlvp detector).
                #assumes the repitition of the detector frame is in the first dimension of the array
                mask_from_image = np.expand_dims(mask_from_image, axis=0)
                mask_from_image = np.repeat(mask_from_image, mask_from_image.shape[0], axis=0)
            im_mask = np.asarray(im_mask) | ma.asarray(mask_from_image)
    
        if "polygon" in mask:
            polygons = mask['polygon']
            for i in polygons:
                img = Image.new('L', im_ints.shape, 0)
                ImageDraw.Draw(img).polygon(i, outline=1, fill=1)
                im_mask = im_mask +ma.array(img)
    
        if "threshold" in mask:
            threshold = mask['threshold']
            im_mask = np.asarray(im_mask) | ma.masked_outside(im_ints,threshold[0],threshold[1]).mask

        if "detector" in mask:
            # Masks for energy dispersive detector elements for full detector.
            # Works by removing complete detectors.
            for x in range(len(mask["detector"])):
                im_mask[mask["detector"][x] - 1] = True

        if "energy" in mask:
            raise ValueError("'Energy' is not implemented.")
            
        if ("two theta" in mask) or ("twotheta" in mask):
            if ("two theta" in mask):
                lbl_str = 'two theta'
            else:
                lbl_str = 'twotheta'
            limits = mask[lbl_str]
            im_mask = np.asarray(im_mask) | ma.masked_inside(self.tth,limits[0], limits[1]).mask
            
        if ("azm" in mask) or ("azimuth" in mask):
            if ("azm" in mask):
                lbl_str = 'azm'
            else:
                lbl_str = 'azimuth'
            limits = mask[lbl_str]
            im_mask = np.asarray(im_mask) | ma.masked_inside(self.azm,limits[0], limits[1]).mask
    
        # FIX ME: Should also add circles and other polygons as per GSAS-II masks
    
        #mask invalid values
        mask2 = ma.masked_invalid(im_ints).mask
        #mask everything less than 0. 
        mask3 = ma.masked_less(im_ints, 0).mask
    
        #combine masks
        im_mask = np.asarray(im_mask) | np.asarray(mask2) | np.asarray(mask3)
        #im_ints = ma.array(im_ints, mask=im_mask)
    
        """
        if debug:
            # N.B. The plot is a pig with even 1000x1000 pixel images and takes a long time to render.
            if ImTTH.size > 100000:
                logger.info(" ".join(map(str, [("Have patience. The mask plot will appear but it can take its time to render.")])))
            fig_1 = plt.figure()
            ax1 = fig_1.add_subplot(1, 3, 1)
            ax1.scatter(ImTTH, ImAzi, s=1, c=(ImInts.data), edgecolors='none', cmap=plt.cm.jet, vmin=0,
                        vmax=np.percentile(ImInts.flatten(), 98))
            ax1.set_title('All data')
    
            ax2 = fig_1.add_subplot(1, 3, 2)
            ax2.scatter(ImTTH, ImAzi, s=1, c=im_mask, edgecolors='none', cmap='Greys')
            ax2.set_title('Mask')
            ax2.set_xlim(ax1.get_xlim())
    
            ax3 = fig_1.add_subplot(1, 3, 3)
            ax3.scatter(ImTTH, ImAzi, s=1, c=(ImInts), edgecolors='none', cmap=plt.cm.jet, vmin=0,
                        vmax=np.percentile(ImInts.flatten(), 98))
            ax3.set_title('Masked data')
            ax3.set_xlim(ax1.get_xlim())
            # ax2.colorbar()
            plt.show()
    
            plt.close()
        """
        """
        Energy dispersive mask debug
        #     # FIX ME: DMF update to use class plot function!
        #     if debug:
        #         # Plot mask.
        #         # This is left in here for debugging.
        #         fig = plt.figure()
        #         ax = fig.add_subplot(1, 2, 1)
        #         plt.subplot(121)
        #         plt.scatter(
        #             im_two_theta,
        #             im_azimuth,
        #             s=4,
        #             c=im_ints,
        #             edgecolors="none",
        #             cmap=plt.cm.jet,
        #         )
        #         ax.set_ylim([0, 360])
        #         ax = fig.add_subplot(1, 2, 2)
        #         plt.subplot(122)
        #         plt.scatter(
        #             ma.array(im_two_theta, mask=im_mask),
        #             ma.array(im_azimuth, mask=im_mask),
        #             s=4,
        #             c=im_ints,
        #             edgecolors="none",
        #             cmap=plt.cm.jet,
        #         )
        #         ax.set_ylim([0, 360])
        #         plt.colorbar()
        #         plt.show()
        #         plt.close()
        """
        
        self.original_mask = im_mask
        
        return im_mask
    
    
    
    
    def set_mask(self,
        intensity_bounds=[-np.inf, np.inf],
        range_bounds=[-np.inf, np.inf],
        azm_bounds=[-np.inf, np.inf],
        mask=None ):
        """
        Set limits to data
        :param range_bounds:
        :param i_max:
        :param i_min:
        :return:
        """
        
        if mask==None:
            mask = self.intensity.mask
            logger.debug(" ".join(map(str, [("Retreived previous mask")])))
        local_mask = np.where(
            (self.tth >= range_bounds[0]) & (self.tth <= range_bounds[1]) &
            (self.azm >= azm_bounds[0])   & (self.azm <= azm_bounds[1]),
            False, True
        )

        #ma.masked_outside(self.tth,(limits[0]),(limits[1])).mask
        local_mask2 = ma.masked_outside(self.intensity, intensity_bounds[0], intensity_bounds[1]).mask
        
        logger.debug(" ".join(map(str, [("type(mask), type(local_mask)", type(mask), type(local_mask))])))
        logger.debug(" ".join(map(str, [("mask.shape, (local_mask.shape)", mask.shape, (local_mask.shape))])))
        combined_mask = np.ma.mask_or(mask, local_mask)
        combined_mask = np.ma.mask_or(combined_mask, local_mask2)
        NoneType = type(None)
        if not isinstance(mask, NoneType):  # or mask.all() != None:
            combined_mask = np.ma.mask_or(combined_mask, mask)
            
        #apply mask to all arrays
        self.mask_apply(combined_mask)
        
    
    def mask_apply(self, mask, debug=False):
        """
        Applies mask to al the data arrays of the detector class

        Parameters
        ----------
        mask : array
            Mask. The default is None.
        debug : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        self.intensity.mask = mask
        self.tth.mask = mask
        self.azm.mask = mask
        
        if "dspace" in self.__dict__.items():
            self.dspace.mask = mask
        
        if "x" in self.__dict__.items():
            self.x.mask = mask
        if "y" in self.__dict__.items():         
            self.y.mask = mask
        if "z" in self.__dict__.items():
            self.z.mask = mask
        
    


    def mask_restore(self):
        """
        Restores the loaded mask.
        If the image data is still the same size as the original.
        """

        logger.effusive(" ".join(map(str, [("Restore original mask.")])))
        self.mask_apply(self.original_mask)





    def mask_remove(self):
        """
        Revmoes all masks from the data arrays.
        """

        logger.effusive(" ".join(map(str, [("Remove all masks.")])))
        self.mask_apply(None)