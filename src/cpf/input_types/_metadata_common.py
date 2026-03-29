#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import re
import numpy as np
import fabio
from pathlib import Path
from datetime import datetime
from dateutil.parser import parse

# from cpf.settings import is_settings
from cpf import h5_functions
from cpf.util.logging import get_logger

logger = get_logger("cpf.input_types._metadata")



class _metadata_common:
    """
    Class definiing metadata function(s).

    These are imported into the detector functions as methods.
    """
    
    def get_metadata(self, metadata_values="default", report=None):
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
        image : open image file, optional
            Name of the image set to import. Either this or settings are required.
            The default is None.
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
        """
        The meta data is expected to be in the form of a dictionary -- or something that reads to a dictionary.
        
        
        """
        
        no_exposure_message = None#"no exposure"
        
        if isinstance(metadata_values, str):
            metadata_values = [metadata_values]
        
        # options for times
        time_opts = ["time_mid" ,"time_start" ,"time_end", "time_exposure"]
        
        # make output dictionary
        metadata_out = {}
        #set metadata to return 
        if metadata_values == ['all'] or metadata_values == 'all':
            metadata_values = list(self.metadata)
            metadata_values += time_opts
        if "default" in metadata_values:
            metadata_values += time_opts
            metadata_values.remove('default')
        time_location = self.metadata_labels.get('time_label', self._default_metadata_labels['time_label'])
        exposure_location = self.metadata_labels.get('exposure_label', self._default_metadata_labels['exposure_label'])
        
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
            # metadata["FILE_CREATION"] = []
            # for k,i in enumerate(imgs_):
            #     metadata["FILE_CREATION"].append(os.path.getctime(image_name))
            metadata_out["FILE_CREATION"] = self.metadata["FILE_CREATION"]
        if "FILE_MODIFIED" in metadata_values:
            # get file modificaction time from OS   
            # metadata["FILE_MODIFIED"] = []
            # for k,i in enumerate(imgs_):
            #     metadata["FILE_MODIFIED"].append(os.path.getmtime(image_name))
            metadata_out["FILE_MODIFIED"] = self.metadata["FILE_CREATION"]
        #get information from inside datafiles
        if any("*" in x for x in metadata_values) or any("/" in x for x in metadata_values):
            # only hdf5 files should have a "*" as wildcard in the keys.
            # to be sure also check for '/' as a key seperator. 
                        
            #needs --> to be in self.metadata
            # image_name (in the list form) --> get from settings. 
            # self.h5_datakey --> to be used to get wildcard values for metadata keys
            # metadata_values
            
            #get imagename from the meta data
            imagename = self.metadata['image']           
            if not isinstance(imagename, list):
                # single file.
                # imagename = [imagename, self.h5_datakey, [0], '0']
                imagename = [imagename, self.metadata['image'][1], [0], '0']

            # get iteration number from h5_datakey and imagename
            # reverse replacement of the wildcard
            regexp_alphanum = '([-+]?[0-9a-zA-Z-+_]*[.][0-9a-zA-Z-+_]+|[-+]?[0-9a-zA-Z-+_]+)'
            # search_term = re.sub('[*]', regexp_alphanum, self.h5_datakey)
            search_term = re.sub('[*]', regexp_alphanum, self.metadata['h5_datakey'])
            iteration = re.search(search_term, imagename[1])
            
            # check metadata requirements exist 
            # add entries to output dictionary
            for j in metadata_values:
                if len(iteration.groups()) > 1:
                    err_str = f"There is more than 1 wildcard in the h5 key {self.h5_datakey}. This is not implemented here."
                    raise NotImplementedError(err_str)
                else:
                    metadata_key = re.sub('[*]', iteration.groups()[0], j)
                    
                #get last index in key as the dictionarry entry label
                ky = metadata_key.split("/")[-1]
                metadata_out[ky] = h5_functions.get_images([imagename[0], metadata_key, imagename[2], '0'])
                try:
                    metadata_out[ky] = h5_functions.get_images([imagename[0], metadata_key, imagename[2], '0'])
                except:                    
                    err_str = f"Metadata type {metadata_key} not recognised. Permitted values for this dataset are: any valid h5 key"
                    raise ValueError(err_str)
                    
        else: # image(s) are separate tiff, edf, etc. images
            # check metadata requirements exist and add entries to output dictionary
            headers = list(self.metadata)
            for j in metadata_values:
                if j in ["FILE_CREATION", "FILE_MODIFIED"]:
                    pass
                elif j in headers:
                    # metadata_out[j] = []
                    try:
                        metadata_out[j] = float(self.metadata[j])
                    except:
                        metadata_out[j] = self.metadata[j]
                else:
                    err_str = f"Metadata type '{j}' not recognised. Permitted values for this dataset are: {headers}."
                    raise ValueError(err_str)
            
        if replaced != None:
            if "time" in replaced:
                if exposure_location is not None:
                    metadata_out["time"] = time_combine(metadata_out[time_location], metadata_out[exposure_location], second_time_scale=1/2)
                else: 
                    metadata_out["time"] = no_exposure_message
            if "time_mid" in replaced:
                if exposure_location is not None:
                    metadata_out["time_mid"] = time_combine(metadata_out[time_location], metadata_out[exposure_location], second_time_scale=1/2)
                else: 
                    metadata_out["time_mid"] = no_exposure_message
            if "time_start" in replaced:
                metadata_out["time_start"] = time_combine(metadata_out[time_location])
                    # metadata_out["time_start"] = time_combine(self.metadata[time_location])
            if "time_end" in replaced:
                if exposure_location is not None:
                    if isinstance(metadata_out[time_location], list):
                        last = metadata_out[time_location][-1]
                        last_exp = metadata_out[exposure_location][-1]
                    else:
                        last = metadata_out[time_location]
                        last_exp = metadata_out[exposure_location]
                    metadata_out["time_end"] = time_combine(last, last_exp, second_time_scale=1)
                else: 
                    metadata_out["time_end"] = no_exposure_message
            if "time_exposure" in replaced:
                metadata_out["time_exposure"] = metadata_out.get(exposure_location, no_exposure_message)
                # if exposure_location is not None:
                #     metadata_out["time_exposure"] = self.metadata[exposure_location]
                # else: 
                #     metadata_out["time_exposure"] = no_exposure_message
            for k in discard:
                metadata_out.pop(k, None)
            
        # collapse everything else down.
        for i in list(metadata_out):
            if isinstance(metadata_out[i], list) and len(metadata_out[i]) > 1:
                try:
                    metadata_out[i] = np.nanmean(metadata_out[i])
                except:
                    metadata_out[i] = metadata_out[i][np.int_(len(metadata_out[i])/2)] 
            elif isinstance(metadata_out[i], list) and len(metadata_out[i]) == 1:
                metadata_out[i] = metadata_out[i][0]
            elif isinstance(metadata_out[i], np.ndarray) and metadata_out[i].size > 1:
                metadata_out[i] = np.nanmean(metadata_out[i])
            else:
                pass                    
    
        return metadata_out


    def _get_file_created_modified(self, medtadata_dict, image):
        """
        Adds file creation and modification times to the metadata dictionary.

        Parameters
        ----------
        medtadata_dict : dict
            dictionary of matadata.
        image : Pth, str
            location of the image.

        Returns
        -------
        medtadata_dict : dict
            dictionary of matadata.

        """
        # append times to dictionary incase of multiple files. 
        if "FILE_CREATION" not in medtadata_dict:
            medtadata_dict["FILE_CREATION"] = os.path.getctime(image)
        else:
            if not isinstance(medtadata_dict["FILE_CREATION"], list):
                medtadata_dict["FILE_CREATION"] = [medtadata_dict["FILE_CREATION"]]
            medtadata_dict["FILE_CREATION"].append(os.path.getmtime(image))

        if "FILE_MODIFIED" not in medtadata_dict:
            medtadata_dict["FILE_MODIFIED"] = os.path.getmtime(image)
        else:
            if not isinstance(medtadata_dict["FILE_MODIFIED"], list):
                medtadata_dict["FILE_MODIFIED"] = [medtadata_dict["FILE_MODIFIED"]]
            medtadata_dict["FILE_MODIFIED"].append(os.path.getmtime(image))
        
        return medtadata_dict
        
        

def time_combine(first_time, second_time=0, second_time_scale=1):
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
            if isinstance(first_time[i], str):
                first_time[i] = parse(first_time[i]).timestamp()
    elif isinstance(first_time, str):
        as_str = True
        first_time = parse(first_time).timestamp()
    first_time = np.nanmean(first_time)
    if isinstance(second_time, list):
        for i in range(len(second_time)):
            if isinstance(second_time[i], str):
                second_time[i] = parse(second_time[i]).timestamp()
    elif isinstance(second_time, str):
        as_str = True
        second_time = parse(second_time).timestamp()
    second_time = np.nanmean(second_time)
        
    out_time = first_time + second_time * second_time_scale
        
    if as_str:
        out_time = f"{datetime.fromtimestamp(out_time):%Y-%d-%b %H:%M:%S.%f}"
        
    return out_time
    
    
    
    
    

