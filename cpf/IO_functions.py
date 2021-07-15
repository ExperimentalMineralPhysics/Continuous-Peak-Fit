#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 05:44:13 2021

@author: simon
"""

import numpy as np
import numpy.ma as ma
import os


# Needed for JSON to save fitted parameters.
# Copied from https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array#24375113
# on 13th November 2018
def json_numpy_serialzer(o):
    """
    Serialize numpy types for json
    Parameters:
        o (object): any python object which fails to be serialized by json
    Example:
        >>> import json
        >>> a = np.array([1, 2, 3])
        >>> json.dumps(a, default=json_numpy_serializer)
    """
    numpy_types = (
        np.bool_,
        # np.bytes_, -- python `bytes` class is not json serializable     
        # np.complex64,  -- python `complex` class is not json serializable  
        # np.complex128,  -- python `complex` class is not json serializable
        # np.complex256,  -- special handling below
        # np.datetime64,  -- python `datetime.datetime` class is not json serializable
        np.float16,
        np.float32,
        np.float64,
        # np.float128,  -- special handling below
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        # np.object_  -- should already be evaluated as python native
        np.str_,
        np.timedelta64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.void,
    )

    if isinstance(o, np.ndarray):
        return o.tolist()
    elif isinstance(o, numpy_types):
        return o.item()
    elif isinstance(o, np.float128):
        return o.astype(np.float64).item()
    # elif isinstance(o, np.complex256): -- no python native for np.complex256
    #     return o.astype(np.complex128).item() -- python `complex` class is not json serializable 
    else:
        raise TypeError("{} of type {} is not JSON serializable".format(repr(o), type(o)))


def FileList(FitParameters, FitSettings):
    """
    From the Settings make a list of all the data files.
    This function is called by the output writing scripts to make sure the file names are called consistently.
    :param FitParameters:
    :param FitSettings:
    :return:
    """
    # Define step
    if 'datafile_Step' not in FitParameters:
        Step = 1
    else:
        Step = FitSettings.datafile_Step
    
    # Diffraction patterns -- make list of files
    if 'datafile_Files' not in FitParameters:
        n_diff_files = int(np.floor((FitSettings.datafile_EndNum - FitSettings.datafile_StartNum)/Step + 1))
        diff_files = []
        for j in range(n_diff_files):
            # Make list of diffraction pattern names and no. of pattern
            n = str(FitSettings.datafile_StartNum+j*Step).zfill(FitSettings.datafile_NumDigit)
            # Append diffraction pattern name and directory
            diff_files.append(os.path.abspath(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename
                                              + n + FitSettings.datafile_Ending))
    elif 'datafile_Files' in FitParameters:
        n_diff_files = int(len(FitSettings.datafile_Files)/Step)
        diff_files = []
        for j in range(n_diff_files):
            # Make list of diffraction pattern names and no. of pattern
            n = str(FitSettings.datafile_Files[j*Step]).zfill(FitSettings.datafile_NumDigit)
            # Append diffraction pattern name and directory
            diff_files.append(os.path.abspath(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename
                                              + n + FitSettings.datafile_Ending))
    else:
        n_diff_files = int(len(FitSettings.datafile_Files)/Step + 1)
    return diff_files, n_diff_files




def AnyTermsNull(obj_to_inspect, val_to_find = None, indexpath="", clean=None):
    ''' This function accepts a nested dictionary and list as argument
        and iterates over all values of nested dictionaries and lists.
        If any of the values are "Null" it returns 0
    '''  
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021
    if clean == None:
        clean=1
        
    if isinstance(obj_to_inspect, dict):
        for key, value in obj_to_inspect.items():
            clean = AnyTermsNull(value, val_to_find, indexpath + f"['{key}']", clean=clean)
    
    if isinstance(obj_to_inspect, list):
        for key, value in enumerate(obj_to_inspect):
            clean = AnyTermsNull(value, val_to_find, indexpath + f"[{key}]", clean=clean)
              
    if obj_to_inspect == val_to_find:
        clean = 0
        print(f"Value {val_to_find} found at {indexpath}")

    return clean


def ReplaceNullTerms(obj_to_inspect, val_to_find = None, indexpath="", clean=None):
    ''' This function accepts a nested dictionary and list as argument
        and iterates over all values of nested dictionaries and lists.
        If any of the values are "Null" it returns 0
    '''  
    # copied from https://python-forum.io/thread-24856.html
    # on 26th June 2021
    
    if isinstance(obj_to_inspect, dict):
        for key, value in obj_to_inspect.items():
            obj_to_inspect[key] = ReplaceNullTerms(value, val_to_find, indexpath + f"['{key}']", clean=clean)
    
    if isinstance(obj_to_inspect, list):
        for key, value in enumerate(obj_to_inspect):
            obj_to_inspect[key] = ReplaceNullTerms(value, val_to_find, indexpath + f"[{key}]", clean=clean)
              
    if obj_to_inspect == val_to_find:
        obj_to_inspect = np.nan
        print(f"Value {val_to_find} found at {indexpath}")

    return obj_to_inspect



def peak_string(orders, fname=False):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param fname: str
    :return: string listing peak names
    """
    # FIX ME: this will be a data-type function so should probably move somewhere else in the end.
    p_str = ''
    
    for x in range(len(orders['peak'])):
        if 'phase' in orders['peak'][x]:
            p_str = p_str + orders['peak'][x]['phase']
        else:
            p_str = p_str + "Peak"
        if fname is False:
            p_str = p_str + " ("
        else:
            p_str = p_str + "-"
        if 'hkl' in orders['peak'][x]:
            p_str = p_str + str(orders['peak'][x]['hkl'])
        else:
            p_str = p_str + str(x+1)
        if fname is False:
            p_str = p_str + ")"
                
        if x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
            if fname is False:
                p_str = p_str + " & "
            else:
                p_str = p_str + '_'

    return p_str



def make_outfile_name(basefilename, directory=None, additional_text=None, extension=None, orders=None, overwrite=True):
    ''' Make file names for output files.
    Needs to account for different file name formats 
    '''
        
    # strip directory if it is in the name and there is a new directory
    if directory:
        basefilename = os.path.basename(basefilename)
        
    if basefilename: # catch incase the string does not exist.
        filename,ending = os.path.splitext(basefilename)   
        ending = ending[1:]
        #filename = os.path.splitext(os.path.basename(fnam))[0] + '_'
    else:
        filename = 'Fit2Peak'
        ending = ''
        
    if filename[-1:] == '_': #if the base file name ends in an '_' remove it. 
        filename = filename[0:-1]
    
    # check the iteration is not in the file ending.
    # another way to do this would be to check in the file ending in the input file is empty.
    try:
        int(ending)
        in_ending = True  #the file ending is convertable to a number -- therefore assume the extension contains the iterating parameter.
    except:
        in_ending = False
    if in_ending:
        filename = filename + '_' + ending
    
    if orders: # add phase and hkl to file name
        filename = filename + '__'+ peak_string(orders, fname=True)
    if additional_text: # add additional text -- i.e. notes from orders
        filename = filename + '__' + additional_text
    if directory:
        filename = os.path.join(directory, filename)
    
    if extension and extension[0]==".":
        extension = extension[1:]
    
    if overwrite==False:
        i = 0
        if os.path.exists('{}.{}'.format(filename, extension)):
            i+=1
        while os.path.exists('{}_{:d}.{}'.format(filename, i, extension)):
            i += 1
        if i == 0:
            filename = ('{}.{}'.format(filename, extension))
        else:
            filename = ('{}_{:d}.{}'.format(filename, i, extension))
    else: 
        filename = ('{}.{}'.format(filename, extension))
    
    return filename



def lmfit_fix_int_data_type(fname):
    """
    fixes problem with lmfit save/load model. 
    lmfit load model cannot read int32 data with nulls in it. 
    if replace 'int32' with 'float32' it will read. 
    """
            
    ObjRead = open(fname, "r")   
    txtContent = ObjRead.read();
    ObjRead.close()
    
    txtContent = txtContent.replace('uint', 'float')
    txtContent = txtContent.replace('int', 'float')
    
    print('    Rewriting', fname)
    ObjRead = open(fname, "w")
    ObjRead.write(txtContent)
    ObjRead.close()
    
