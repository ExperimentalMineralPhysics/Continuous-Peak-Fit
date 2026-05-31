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


"""
Metadata in continous peak fit.
===============================

Cpf has the options to read metadata from the diffraction files and record it with the fits in the json files and is added to the output files.
Metadata that is required by an output type listed in the Settings() is by default added to the saved metadata.

The 'settings.metadata' is a list of metadata parameter names that will be added to the {fits}.json files as key,value pairs. 
This is formed from:
    (a) the 'metadata' list in the inputfile.
    (b) names required by the outputs, as held in settings.metadata_labels dictionary (see below)
    (c) hdf5 keys in h5_iterate (see below and hdf5 documentation).
    (d) keys returned by metadata_func

For non-hdf5 file types all the metadata is imported from the image file along with the diffraction data by data_class.import_image() 
and is stored in data_class.metadata.
For hdf file image types the number of possible metdata data values is vast, so only the named values are given.

Only the named metadata is added to the saved files. 


Possible metadata names
-----------------------

Possible metadata parameters are anything listed in the data file and accessible by the data_class.  

For hdf5 type files, any valid key is allowed but needs the same number of '*' wild cards to get metadata corresponding to the right images (see h5 documentation).

There is also a set of metadata values designed for use with compound diffraciton images (e.g. ESRF lvp detector, or multiple exposures making a single diffraction collection)
These also work for single frame images and are derived using default "time" and 'exposure' metadata_labels 
The options for this are:
    - "frames start"    -- start time of data collection
    - "frames mid"      -- median time of data collection, accounting for the exposre time
    - "frames end"      -- end time of data collection, accounting for the exposre time
    - "frames exposure" -- total exposure time of data collection, accounting for the frame time
If no exposure time is specified then only "frames start" is returned.

A generic backstop of times is also avaliable in for form of:
    - 'FILE_CREATION'  -- file creation time read from operatiing system
    - 'FILE_MODIFIED'  -- file modification time read from operatiing system


metadata_labels and RequiredParams
----------------------------------

The metadata names or locations are specific to the data type. For example, the time of the data collection trigger or acquition time are called different 
names by different data classes and by different synchroton facilites. 
Therefore parameters called by Write{name of output}.Requirements()["RequiredParams"] are generic pointers, that are speficifed in the data_class.
Each data_class has attribute 'metadata_labels', which is a dictionary that points the RequiredParams to the specific metadata location.
for example:
    data_class.metadata_labels = {"time":        "mean_start_time",
                                  'exposure':    'mean_live_time',
                                  "temperature": "*LVP_tc1_calcs.I"}
These values are set by:
    - initially from data_class._default_metadata_labels*
    - overriding the defaults with values from 'metadata_labels' in input

Dictionary keys from data_class.metadata_labels are added to settings.metadata and saved to json fits files. 


Adding external matadata values -- metadata_read_func
-----------------------------------------------------

metadata_read_func is a method defined by the input that allows interation with external (to the data files) sources of metadata.
It is added to Settings as:
    Settings.metadata_read_func()
from function called `metadata_read_func' in the input.

It returns a dictionary of values to be added to the metadata. The method must:
 - accept 3 optional arguments (settings, image_obj, filename) and kwargs
 - return dictionary of metadata key, value pairs
 - return a dictionary of the same key, value pairs even if the inputs are absent 
        or fail. The empty dict is used to determine the metadata to add to the outputs. 

An example is: 

code block:
        def metadata_read_func(settings=None, image_obj=None, filename=None, **kwargs):
            # example of a metadata_read_func     
            # taken from Example 1
            try:
                hd = image_obj.header
                out = {}
                out["Exposure_time"] = hd["Exposure_time"]
                out["Exposure_period"] = hd["Exposure_period"]
                out['namename'] = image_obj.filename
                out['daft'] = "Needless string"
            except:
                out = {}
                out["Exposure_time"] = None
                out["Exposure_period"] = None
                out['namename'] = None
                out['daft'] = "Needless string"
                
            return out
        
If the outputs of metadata_read_func are required by the outputs, then 'data_calss.metadata_labels' values can be modified to use them.


TO DO:
    1. rename metadata_labels as metadata_names ?? 
    2. write a data_class.metadata_possibilities() method that lists all possible metadata names. 
    3. change definitions of metadata_read_func. strictly it does not need all three of the forced parameters and infact never uses them all 
"""

class _metadata_common:
            
    def get_metadata(self, settings_class= None, metadata_values="default", report=None):
        """
        Gets speficied metadata as from settings_class.metadata and returns as a dictionary.
        
        The default is to get the timestamps of the images but other data 
        stored in the file can be accessed as well with the correct label.
        
        For h5 files any hdf5 key is allowed but needs '*' wild cards to move through the dataset
        with the images (see h5 documentation).
        For edf and other file formats any property accessible through fabio image.header is permitted.
        
        Otherfile types are parsed through fabio and behave accordingly.
    
        The implemented options for the metadata parameters are:
            - "frames start"    -- start time of data collection
                                        min(timestamps)
            - "frames mid"      -- median time of data collection, accounting for the exposre time
                                        np.median(timestamps+exposures)
            - "frames mean"     -- mean time of data collection, accounting for the exposre time  
                                        (timestamps+exposures) / len(timestamps)
            - "frames end"      -- end time of data collection, accounting for the exposre time
                                        max(timestamps+exposures)
            - "frames exposure" -- total exposure time of data collection, accounting for the frame time
                                        max(timestamps+exposures) - min(timestamps)
        If no exposure time is determined then only "time_start" is returned.
       
    
        Parameters
        ----------
        image : open image file, optional
            Name of the image set to import. Either this or settings are required.
            The default is None.
            image_name : string, list, optional
                Name of the image set to import. Either this or settings are required.
                The default is None.
        settings_class : settings class, optional
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
        # parse inputs
        if settings_class and metadata_values=="default":
            metadata_values = settings_class.metadata
        if isinstance(metadata_values, str):
            metadata_values = [metadata_values]
        
        no_exposure_message = None#"no exposure"
        # options for times
        time_opts = ["frames start", "frames mid", "frames mean", "frames end", "frames exposure"]
        # time_opts = ["time_mid" ,"time_start" ,"time_end", "time_exposure"]
        
        # make output dictionary
        metadata_out = {}
        #set metadata to return 
        if metadata_values == ['all'] or metadata_values == 'all':
            metadata_values = list(self.metadata)
            metadata_values += time_opts
        if "default" in metadata_values:
            metadata_values += time_opts
            metadata_values.remove('default')
        time_location = self.metadata_labels.get('time', self._default_metadata_labels.get('time',None))
        exposure_location = self.metadata_labels.get('exposure', self._default_metadata_labels.get('exposure', None))
        
        #parse metadata_values list
        discard = []
        if time_location not in metadata_values:
            discard.append(time_location)
        # replace time_opts values that are not in the metadata
        if any(map(lambda v: v in time_opts, metadata_values)):
            replaced = list(set(metadata_values) & set(time_opts))
            metadata_values = list(set(metadata_values) - set(time_opts))
            metadata_values.append(time_location)
            if (("frames end" in replaced or "frames mid" in replaced or "time" in replaced)
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
        if "FILE_CREATION" in metadata_values and "FILE_CREATION" in self.metadata:
            metadata_out["FILE_CREATION"] = self.metadata["FILE_CREATION"]
        if "FILE_MODIFIED" in metadata_values and "FILE_MODIFIED" in self.metadata:
            metadata_out["FILE_MODIFIED"] = self.metadata["FILE_MODIFIED"]
        #get information from inside datafiles
        if "h5_datakey" in self.metadata:
            # only hdf5 files should have a "h5_datakey" as wildcard in the keys.
            # this is set by cpf.
            
            #if any("*" in x for x in metadata_values) or any("/" in x for x in metadata_values):
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
                if "/" in j:
                    # h5 key must have / in the name. No / neams not h5 key
                    if len(iteration.groups()) > 1:
                        err_str = f"There is more than 1 wildcard in the h5 key {self.h5_datakey}. This is not implemented here."
                        raise NotImplementedError(err_str)
                    elif "*" in j:
                        metadata_key = re.sub('[*]', iteration.groups()[0], j)
                    else:
                        metadata_key = j
                        
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
                if j in ["FILE_CREATION", "FILE_MODIFIED"] or j==None:
                    pass
                elif j in headers:
                    try:
                        metadata_out[j] = float(self.metadata[j])
                    except:
                        metadata_out[j] = self.metadata[j]
                elif "*" in j: # wildcard in header
                    # add all wildcards to header
                    # N.B. this should not be called becuase settings.set_metadata removes all * from metadata.
                    pattern = re.compile(re.sub('[*]', '([0-9a-zA-Z-+_:]*)', j))  
                    matches = [word for word in headers if pattern.match(word)]
                    for k in matches:
                        try:
                            metadata_out[k] = float(self.metadata[k])
                        except:
                            metadata_out[k] = self.metadata[k]
                else:
                    err_str = f"Metadata type '{j}' not recognised. Permitted values for this dataset are: {headers.append(time_opts)}."
                    raise ValueError(err_str)
            
        if replaced != None:
            # if "time" in replaced:
            #     if exposure_location is not None:
            #         metadata_out["time"] = times_combine(metadata_out[time_location], metadata_out[exposure_location])
            #     else:
            #         metadata_out["time"] = times_combine(metadata_out[time_location])
            if "frames start" in replaced:
                metadata_out["frames start"] = metadata_combine(metadata_out, time_location, exposure_location, get="start", as_str=False)
            if "frames mid" in replaced:
                metadata_out["frames mid"] = metadata_combine(metadata_out, time_location, exposure_location, get="mid", as_str=False)
            if "frames mean" in replaced:
                metadata_out["frames mean"] = metadata_combine(metadata_out, time_location, exposure_location, get="mean", as_str=False)
            if "frames end" in replaced:
                metadata_out["frames end"] = metadata_combine(metadata_out, time_location, exposure_location, get="end", as_str=False)
            if "frames exposure" in replaced:
                metadata_out["frames exposure"] = metadata_combine(metadata_out, time_location, exposure_location, get="exposure", as_str=False)
                # metadata_out["frames exposure"] = metadata_out.get(exposure_location, no_exposure_message)
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


    def _get_file_created_modified(self, image):
        """
        Adds file creation and modification times to the metadata dictionary.

        Parameters
        ----------
        image : Pth, str, list
            location of the image or list of images

        Returns
        -------
        medtadata_dict : dict
            dictionary of matadata.

        """
        if not isinstance(image, list):
            image = [image]
            
        metadata_dict={}
        tmp_creations = []
        tmp_modified = []
        for i in image:
            # loop over all images to get data 
            
            # check image input
            if (isinstance(i, str) is False) and (isinstance(i, Path) is False):
                #then image is image object.
                try:
                    i = i.filename
                except:
                    i = i.get_name()
            tmp_creations.append(os.path.getctime(i))
            tmp_modified.append(os.path.getmtime(i))
        
            # # append times to dictionary incase of multiple files. 
            # if "FILE_CREATION" not in metadata_dict:
            #     metadata_dict["FILE_CREATION"] = os.path.getctime(i)
            # else:
            #     if not isinstance(metadata_dict["FILE_CREATION"], list):
            #         metadata_dict["FILE_CREATION"] = [metadata_dict["FILE_CREATION"]]
            #     metadata_dict["FILE_CREATION"].append(os.path.getctime(i))
    
            # if "FILE_MODIFIED" not in metadata_dict:
            #     metadata_dict["FILE_MODIFIED"] = os.path.getmtime(i)
            # else:
            #     if not isinstance(metadata_dict["FILE_MODIFIED"], list):
            #         metadata_dict["FILE_MODIFIED"] = [metadata_dict["FILE_MODIFIED"]]
            #     metadata_dict["FILE_MODIFIED"].append(os.path.getmtime(i))

        metadata_dict["FILE_CREATION"] = np.median(tmp_creations)
        metadata_dict["FILE_MODIFIED"] = np.median(tmp_modified)
        
        return metadata_dict
        

def metadata_combine(metadata, time_location, exposure_location=None, get="mean", as_str=False):
    """
    Combine multiple time stamps and return value in the same format (string or time) as first_time
    
    The values for get are:
        - "start"    -- start time of data collection
                                    min(timestamps)
        - "mid"      -- median time of data collection, accounting for the exposre time
                                    np.median(timestamps+exposures)
        - "mean"     -- mean time of data collection, accounting for the exposre time  
                                    (timestamps+exposures) / len(timestamps)
        - "end"      -- end time of data collection, accounting for the exposre time
                                    max(timestamps+exposures)
        - "exposure" -- total exposure time of data collection, accounting for the frame time
                                    max(timestamps+exposures) - min(timestamps) 
                                    
    Parameters
    ----------
    metadata : TYPE
        DESCRIPTION.
    time_location : str
        key string for time location in metadata dictionary.
    exposure_location : str, optional
        key string for exposure location in metadata dictionary. The default is None.
    get : str, optional
        which value to get. The default is "mean".
    as_str : bool, optional
        Force the timestamps to be strings; if true return in same format as timestamps.
        The default is False.
        
    Raises
    ------
    ValueError
        Unrecognised combination of values.

    Returns
    -------
    out_time : list, float, str
        Require time stamp(s) in the same format as the input.

    """
    
    # parse inputs
    timestamps = metadata[time_location], 
    if exposure_location:
        exposures = metadata[exposure_location]
    else:
        exposures = None
    
    if isinstance(timestamps, str):
        timestamps = [timestamps]
    if exposures == None:
        exposures = 0
    
    #get input format
    if isinstance(timestamps[0], str):
        as_str = True
        
    # force all values to be numbers
    for i in range(len(timestamps)):
        timestamps = [parse(v).timestamp() if isinstance(v, str) else v for v in timestamps]
        
    match get:
        case "start":
             out_time = np.min(np.array(timestamps))
        case "mid":
             out_time = np.median(np.array(timestamps)+np.array(exposures))
        case "mean":
             out_time = np.sum(np.array(timestamps)+np.array(exposures)) / len(timestamps)
        case "end":
             out_time = np.max(np.array(timestamps)+np.array(exposures))
        case "exposure":
             out_time = np.max(np.array(timestamps)+np.array(exposures)) - np.min(np.array(timestamps))
        case _:
            raise ValueError("Unrecognised combination type")
        
    if as_str and get != "exposure":
        out_time = f"{datetime.fromtimestamp(out_time):%Y-%d-%b %H:%M:%S.%f}"
    
    if exposure_location is None:
        return None
    else:
        return out_time
    