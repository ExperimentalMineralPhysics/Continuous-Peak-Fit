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
    if ('datafile_Step' in FitParameters):
        Step = FitSettings.datafile_Step
    elif ('datafile_Files' in FitParameters):
        Step = 1
    else: #if just datafile_StartNum and datafile_EndNum
        if FitSettings.datafile_EndNum > FitSettings.datafile_StartNum:
            Step = 1
        else:
            Step = -1
            
    if not 'datafile_NumDigit' in FitParameters:
        FitSettings.datafile_NumDigit = 1
    
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



def peak_string(orders, fname=False, peak='all'):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param fname: switch indicting if the string is to be a file name or not.
    :param peak: switch to add restricted peaks to the output string
    :return: string listing peak names
    """
    if peak=='all':
        peek = list(range(len(orders['peak'])))
    elif not isinstance(peak, list):
        peek = [int(x) for x in str(peak)]
    else:
        peek=peak
    p_str = ''
    for x in peek:
        if 'phase' in orders['peak'][x]:
            p_str = p_str + orders['peak'][x]['phase']
        else:
            p_str = p_str + "Peak"
        if fname is False:
            p_str = p_str + " ("
        else:
            p_str = p_str + "-"
        if 'hkl' in orders['peak'][x]:
            #p_str = p_str + str(orders['peak'][x]['hkl'])
            p_str = p_str + peak_hkl(orders, peak=x, string=True)[0]
        else:
            p_str = p_str + str(x+1)
        if fname is False:
            p_str = p_str + ")"
                
        if peak=='all' and x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
            if fname is False:
                p_str = p_str + " & "
            else:
                p_str = p_str + '_'
    return p_str


def peak_hkl(orders, peak='all', string=True):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param string: switch indicting if the output is a string or numeric list of hkl.
    :param peak: switch to add restricted peaks to the output string
    :return: string or list of peak hkl.
    """
    #possible hkl input formats:
    #   string -- hkls as a string, generally used for hkls with leading 0 or negative hkl indicy 
    #   int -- hkls all positive and no leading 0
    #   list -- e.g. [h,k,l]. New format. 
    #desired outputs:
    #   string -- hkls as a string using String==True
    #   list -- e.g. [h,k,l]. New format. using String==False
    
    if peak=='all':
        peek = list(range(len(orders['peak'])))
    elif not isinstance(peak, list):
        peek = [int(x) for x in str(peak)]
    else:
        peek=peak
        
    out = []
    for x in peek:
        if 'hkl' in orders['peak'][x]:
            hkl = orders['peak'][x]['hkl']
            if hkl == '0' or hkl == 0:
                hkl = '000'
        else:
            hkl = '000'
        
        if not isinstance(hkl, list) and string==True:
            #string in, string out
            out.append(str(hkl))
        elif isinstance(hkl, list) and string==True:
            #convert numbers to string
            out_tmp = ''
            for y in range(len(hkl)):
                out_tmp=out_tmp+str(hkl[y])
            out.append(out_tmp)
        elif isinstance(hkl, list) and string==False:
            #list in, list out
            out.append(hkl)
        elif not isinstance(hkl, list) and string==False:
            #convert string to array
            pos = 0
            hkl = str(hkl)
            if hkl[0] == '-':
                h = hkl[pos:pos+2]
                pos = pos+2
            else:
                h = hkl[pos:pos+1]
                pos = pos+1
            if hkl[pos] == '-':
                k = hkl[pos:pos+2]
                pos = pos+2
            else:
                k = hkl[pos:pos+1]
                pos = pos+1
            if hkl[pos] == '-':
                l = hkl[pos:pos+2]
                pos = pos+2
            else:
                l = hkl[pos:pos+1]
                pos = pos+1
            out_tmp = [h, k, l]
            if len(hkl)>pos:
                if hkl[pos] == '-':
                    m = hkl[pos:pos+2]
                    pos = pos+2
                else:
                    m = hkl[pos:pos+1]
                    pos = pos+1
                out_tmp.append(m)
            out.append(out_tmp)
                
    return out
    

def peak_phase(orders, peak='all'):
    """
    :param orders: list of peak orders which should include peak names/hkls
    :param string: switch indicting if the output is a string or numeric list of hkl.
    :param peak: switch to add restricted peaks to the output string
    :return: string or list of peak hkl.
    """
    #possible hkl input formats:
    #   string -- hkls as a string, generally used for hkls with leading 0 or negative hkl indicy 
    #   int -- hkls all positive and no leading 0
    #   list -- e.g. [h,k,l]. New format. 
    #desired outputs:
    #   string -- hkls as a string using String==True
    #   list -- e.g. [h,k,l]. New format. using String==False
    
    if peak=='all':
        peek = list(range(len(orders['peak'])))
    elif not isinstance(peak, list):
        peek = [int(x) for x in str(peak)]
    else:
        peek=peak
        
    out = []
    for x in peek:
        if 'phase' in orders['peak'][x]:
            out.append(orders['peak'][x]['phase'])
        else:
            out.append('Unknown')
                
    return out

def make_outfile_name(basefilename, directory=None, additional_text=None, extension=None, orders=None, overwrite=True):
    ''' Make file names for output files.
    Needs to account for different file name formats 
    '''
        
    # strip directory if it is in the name and there is a new directory
    if directory or directory=='':
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
    if additional_text: # add additional text 
        filename = filename + '__' + additional_text
        #FIX ME: Lon lernger used. but need ot check this.
    if orders and 'note' in orders: # add additional text from note in orders.
        filename = filename + '__' + "".join(i for i in orders['note'] if i not in "\/:;*?<>|")          
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
        #if i == 0:
        #    #filename = ('{}.{}'.format(filename, extension))
        #    filename = filename
        #else:
        #    #filename = ('{}_{:d}.{}'.format(filename, i, extension))
        #    filename = ('{}_{:d}'.format(filename, i))
        if i != 0:
            filename = ('{}_{:d}'.format(filename, i))
    #else: 
    #    filename = ('{}'.format(filename, extension))
    
    if extension:
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
    
