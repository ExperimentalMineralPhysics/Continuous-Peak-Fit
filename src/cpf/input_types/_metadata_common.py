#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import numpy as np
import fabio
from datetime import datetime
from dateutil.parser import parse

from cpf import h5_functions
from cpf.util.logging import get_logger

logger = get_logger("cpf.input_types._metadata")



class _metadata_common:
    """
    Class definiing metadata function(s).

    These are imported into the detector functions as methods.
    """
    
    def get_metadata(self, image_name=None, settings=None, metadata_values="default"):
        """
        Gets metadata from the diffraction patterns. 
        The default is to get the timestamps of the images but other data 
        stored in the file can be accessed as well with the correct label.
        
        For h5 files any hdf5 key is allowed but needs '*' wild cards to move through the dataset
        with the images (see h5 documentation).
        For edf and other file formats anything that is a property accessible through fabio image.header is permitted.
        
        Otherfile types are parsed through fabio and behave accordingly.
    
        The implemented options for the metadata parameters are:
        - "time start" -- gives start time of data collection
        - "time mid" -- gives middle time of data collection, accounting for the exposre time
        - "time end" -- gives end time of data collection, accounting for the exposre time
        If no exposure time is determined then only "time_start" is returned.
    
        Parameters
        ----------
        image_name : string, list, optional
            Name of the image set to import. Either this or settings are required.
            The default is None.
        settings : settings class, optional
            cpf settings class. Either this or image_name are required.
            The default is None.
        metadata : string, list, optional
            List the metadata types that are desired from the diffraction data.
            The default is "default".
    
        Returns
        -------
        metadata : dict
            dictionary of the metadata.
    
        """
        
        # check inputs
        if image_name == None and settings.subfit_filename == None:
            raise ValueError("Settings are given but no subpattern is set.")
        if image_name == None:
            # use preset file.
            image_name = settings.subfit_filename
            
        if isinstance(metadata_values, str):
            metadata_values = [metadata_values]
        
        # options for times
        time_opts = ["time_mid" ,"time_start" ,"time_end"]
        
        # make output dictionary
        metadata = {}
        #set default metadata 
        if "default" in metadata_values:
            metadata_values += time_opts
            metadata_values.remove('default')
    
        # get metadata from images
        if (isinstance(image_name, list) or 
            os.path.splitext(os.path.basename(image_name))[1] == ".h5"
            ):
            # then it is a *.h5 file containing data from a the spin of the detector.
            time_location = self.metadata_labels.get('time_label')
            exposure_location = self.metadata_labels.get('exposure_label')
        else: 
            # image(s) are separate tiff, edf, etc. images
            time_location = self.metadata_labels.get('time_label')
            exposure_location = self.metadata_labels.get('exposure_label')
            
            #get list of images (incase there is more than 1)   
            if "*" in image_name:
                imgs_, _ = self._get_sorted_files(image_name, reduce_by=self.reduce_by)
            else:
                imgs_ = [image_name]
        
        #parse metadata_values list 
        discard = []
        if time_location not in metadata_values:
            discard.append(time_location)
        # replace time_opts values that are not in the metadata
        if any(map(lambda v: v in time_opts, metadata_values)):
            replaced = list(set(metadata_values) & set(time_opts))
            metadata_values = list(set(metadata_values) - set(time_opts))
            metadata_values.append(time_location)
            if (("time_end" in replaced or "time_mid" in replaced or "time" in replaced)
                and exposure_location is not None):
                # stop
                if exposure_location not in metadata_values:
                    discard.append(exposure_location)
                metadata_values.append(exposure_location)
        elif metadata_values == ["default"]:
            replaced = {time_location: metadata_values[0]}
            metadata_values = [time_location]
        else:
            replaced = None
            
        
        # get metadata from images
        # look for os level properties first FILE_CREATION and FILE_MODIFIED 
        if "FILE_CREATION" in metadata_values:
            # get file creation time from OS    
            metadata["FILE_CREATION"] = []
            for k,i in enumerate(imgs_):
                metadata["FILE_CREATION"].append(os.path.getctime(image_name))
        if "FILE_MODIFIED" in metadata_values:
            # get file modificaction time from OS   
            metadata["FILE_MODIFIED"] = []
            for k,i in enumerate(imgs_):
                metadata["FILE_MODIFIED"].append(os.path.getmtime(image_name))
        #get information from inside datafiles
        if (isinstance(image_name, list) or
            os.path.splitext(os.path.basename(image_name))[1] == ".h5"):
            # list with h5 image names and keys 
            # or single h5 file.
            
            if not isinstance(image_name, list):
                # single file.
                image_name = [image_name, self.h5_datakey, [0], '0']
            
                err_str = "Metadata extraction not set for single h5 file."
                raise ValueError(err_str)
            
            # check metadata requirements exist 
            # add entries to output dictionary
            for j in metadata_values:
                
                regexp_alphanum = '([-+]?[0-9a-zA-Z-+_]*[.][0-9a-zA-Z-+_]+|[-+]?[0-9a-zA-Z-+_]+)'
                # get iteration number from h5_datakey and image name
                search_term = re.sub('[*]', regexp_alphanum, self.h5_datakey)
                iteration = re.search(search_term, image_name[1])
                
                if len(iteration.groups()) > 1:
                    err_str = f"There is more than 1 wildcard in the h5 key {self.h5_datakey}. This is not implemented here."
                    raise NotImplementedError(err_str)
                else:
                    metadata_key = re.sub('[*]', iteration.groups()[0], j)
                
                #get last index in key as the dictionarry entry label
                ky = j                             
                try:
                    metadata[ky] = h5_functions.get_images([image_name[0], metadata_key, image_name[2], '0'])
                except:                    
                    err_str = f"Metadata type {j} not recognised. Permitted values for this dataset are: any valid h5 key"
                    raise ValueError(err_str)
                    
        else: # image(s) are separate tiff, edf, etc. images
            # check metadata requirements exist and add entries to output dictionary
            im_md = self.get_metadata_dictionary(imgs_[0])
            headers = list(self.get_metadata_dictionary(imgs_[0]))
            for j in metadata_values:
                if j in ["FILE_CREATION", "FILE_MODIFIED"]:
                    pass
                elif j in headers:
                    metadata[j] = []
                else:
                    err_str = f"Metadata type '{j}' not recognised. Permitted values for this dataset are: {headers}."
                    raise ValueError(err_str)
            
            # metadate values from images. 
            for k,i in enumerate(imgs_):
                im_md = self.get_metadata_dictionary(i)
                for j in metadata_values:
                    if j in ["FILE_CREATION", "FILE_MODIFIED"]:
                        pass
                    else:
                        try:
                            metadata[j].append(float(im_md[j]))
                        except:
                            metadata[j].append(im_md[j])

        if replaced != None:
            if "time" in replaced:
                if exposure_location is not None:
                    metadata["time"] = time_combine(metadata[time_location], metadata[exposure_location], second_time_scale=1/2)
                else: 
                    metadata["time"] = "no exposure"
            if "time_mid" in replaced:
                if exposure_location is not None:
                    metadata["time_mid"] = time_combine(metadata[time_location], metadata[exposure_location], second_time_scale=1/2)
                else: 
                    metadata["time"] = "no exposure"
            if "time_start" in replaced:
                metadata["time_start"] = metadata[time_location][0]
            if "time_end" in replaced:
                if exposure_location is not None:
                    metadata["time_end"] = time_combine(metadata[time_location][-1], metadata[exposure_location][-1], second_time_scale=1)
                else: 
                    metadata["time"] = "no exposure"
            for k in discard:
                metadata.pop(k, None)
    
        # collapse everything else down.
        for i in list(metadata):
            if isinstance(metadata[i], list) and len(metadata[i]) > 1:
                try:
                    metadata[i] = np.nanmean(metadata[i])
                except:
                    metadata[i] = metadata[i][np.int_(len(metadata[i])/2)] 
            elif isinstance(metadata[i], list) and len(metadata[i]) == 1:
                metadata[i] = metadata[i][0]
            elif isinstance(metadata[i], np.ndarray) and metadata[i].size > 1:
                metadata[i] = np.nanmean(metadata[i])
            else:
                pass                    
    
        return metadata



def time_combine(first_time, second_time, second_time_scale=None):
    """
    Combine two times and return in the same format (string or time) as first_time

    Parameters
    ----------
    how : TYPE
        DESCRIPTION.
    first_time : TYPE
        DESCRIPTION.
    second_time : TYPE
        DESCRIPTION.
    second_time_scale : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    #format of time string
    as_str = False
    if isinstance(first_time, list):
        as_str = True
        for i in range(len(first_time)):
            first_time[i] = parse(first_time[i]).timestamp()
    first_time = np.nanmean(first_time)
    if isinstance(second_time, list):
        for i in range(len(second_time)):
            if isinstance(second_time[i], str):
                second_time[i] = parse(second_time[i]).timestamp()
    second_time = np.nanmean(second_time)
        
    if second_time_scale is not None:
        out_time = first_time + second_time * second_time_scale
        
    if as_str:
        out_time = f"{datetime.fromtimestamp(out_time):' %Y-%d-%b %H:%M:%S.%f'}"
        
    return out_time
    
    
    
    
    

