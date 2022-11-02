#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import hdf5plugin #have to import it to read the files. not sure why
import h5py
import numpy as np
import cv2   
import matplotlib.pyplot as plt
import os
import cpf.IO_functions as IO
import cpf.XRD_FitPattern as fp

# copied from internet somewhere August 2022
def allkeys(obj):
    "Recursively find all keys in an h5py.Group."
    keys = (obj.name,)
    if isinstance(obj, h5py.Group):
        for key, value in obj.items():
            if isinstance(value, h5py.Group):
                keys = keys + allkeys(value)
            else:
                keys = keys + (value.name,)
    return keys

# copied from internet somewhere August 2022
def getkeys(obj, key_list, level):
    
    keys = list(obj[key_list[0]].keys())
    keys = (obj.name,)
    if isinstance(obj, h5py.Group):
        for key, value in obj.items():
            keys = keys + (value.name,)
            
            # if isinstance(value, h5py.Group):
            #     keys = keys + allkeys(value)
            # else:
                
    return keys



def image_key_validate(h5key_list, key_start, key_end, key_step):
    """
    Validates the lengths and structure of the h5 file key values. 
    Does not validate if all the keys exist. 
    

    Parameters
    ----------
    h5key_list : TYPE
        DESCRIPTION.
    key_start : TYPE, optional
        DESCRIPTION. The default is key_start.
    key_end : TYPE, optional
        DESCRIPTION. The default is key_end.
    key_step : TYPE, optional
        DESCRIPTION. The default is key_step.

    Returns
    -------
    None.

    """
    fail = 0
    
    number_indicies = len(h5key_list)
    
    if len(key_start)!=(1|number_indicies):
        print('There are the wrong number of start positions')
        fail=+1
    
    if len(key_end)!=(1|number_indicies):
        print('There are the wrong number of end positions')
        fail=+1
    
    if len(key_step)!=(1|number_indicies):
        print('There are the wrong number of step lengths')
        fail=+1
    
    if fail != 0:
        stop
    
    return fail
    

def unique_labels(labels, i=3, number_data=None):
    """
    Takes the label array and turns it into strings as sensible lengths.
    The array return will be of unique values if the input array was unique,

    Parameters
    ----------
    labels : TYPE
        the list of labels to be rounded/cleaned up.
    i : number
        The number of decimal places to start the rounding at. 

    Returns
    -------
    Labels rounded

    """
    
    # print(type(labels[0]))
    # print(labels)
    
    number_labels = np.size(labels)
    number_unique = len(np.unique(labels))
    
    #print(number_labels, number_unique, number_data)
    if number_labels == number_unique and (number_data==number_labels or number_labels==1):
        
        done = 0 
        
        while done==0:
            
            #print(i)
                
            try: #catch the unique labels not being numbers
                labels_rounded = np.round(labels, i)
                
                if len(np.unique(labels_rounded)) == number_labels:
                    done=1
                    
                    labels = np.round(labels, i+1)
                    labels = np.array(labels)
                    
                else:
                    i+=1
                #print('pants')
            except:
                done= 1
                #print('ohdear')
                pass
                
    else:
        print('there are insufficient unique labels')
        stop
        
    return labels 
    

    

def get_image_keys(datafile, h5key_list, h5key_names, key_start=0, key_end=-1, key_step=1, key_route='', bottom_level='iterate', index=0, key_str = ""):
    
    #validate the inputs
    #image_key_validate(h5key_list, key_start=key_start, key_end=key_end, key_step=key_step)
    
    
    if isinstance(datafile, str):
        datafile = h5py.File(datafile, 'a')
    
    keys=[]
    sep1 = "_"
    sep2 = "="
    # ...
    #print('called again')
    #print(h5key_list)
    
    if len(h5key_list)==1:
        # we are at the end of the list and need to find out how many values are here, or what the values are
        key = key_route+'/'+h5key_list[0]
        #print('bottom', key)
        
        #print(bottom_level)
        # make sure the key exists
        if key in datafile.keys():
            
            if bottom_level=='sum':
                
                number_data = datafile[key].shape[index]
                #print('number_data', number_data)
                
                if key_end[0] == -1:
                    key_end[0] = number_data+1
                else:
                    key_end[0]+1
                                     
                # if the order is reverse switch the values
                # FIXME: this only works for the input as set needs to be updatedto work for all logics.
                # FIXME: Should use cpf.IO_functions/file_list work for this. 
                if key_step[0] <= -1:
                    key_start[0], key_end[0] = key_end[0], key_start[0]
                    key_end[0]-=1
                    key_start[0]-=1
                
                index_values = list([*range(key_start[0], key_end[0], key_step[0])])
                
                # get the labels.
                # FIXME: this doesnt work if the last enty in the h5key_names list is blank. 
                # FIXME: Should use a function to cover all options and combine with iterate call.           
                labels = key_route+datafile[key_route+'/'+h5key_names[0]][()]
                # FIXME! this is the line to use if the last entry in the list in empty  
                #labels = str(key_str)
                
                return [[key, index_values, labels]]
            
            
            elif bottom_level=='iterate':
                # iterate over the size of the array in the h5 group.
                
                out = []
                
                #print("data set with the key:", datafile[key])
                number_data = datafile[key].shape[index]
                #print('number_data', number_data)
                
                # if we are using all the data make sure we run to the end.
                #print('keyend',key_end)
                if key_end[0] == -1:
                    key_end[0] = number_data
                else:
                    key_end[0]+1
                
                # if we are looking for more images than there are cut the list
                #key_end[0] = np.min(number_data, key_end[0])
                                    
                # if the order is reverse switch the values
                # FIXME: this only works for the input as set needs to be updatedto work for all logics.
                # FIXME: Should use cpf.IO_functions/file_list work for this. 
                if key_step[0] <= -1:
                    key_start[0], key_end[0] = key_end[0], key_start[0]
                    key_end[0]-=1
                    key_start[0]-=1
                                    
                # get the labels
                # FIXME: this doesnt work if the last enty in the h5key_names list is blank. 
                # FIXME: Should use a function to cover all options and combine with iterate call.
                labels = []
                if isinstance(h5key_names[0], list):
                    for i in range(len(h5key_names[0])):
                        #print("key location", str(i), key_route+'/'+h5key_names[0][i])
                        labels_temp = datafile[key_route+'/'+h5key_names[0][i]][()]
                        #print(labels_temp)
                        labels_temp = unique_labels(labels_temp, number_data = number_data)
                        labels.append(labels_temp)
                elif isinstance(h5key_names[0], str) and h5key_names[0] is not "": #it must be a string.
                    #print('string')
                    #print(key_route+'/'+h5key_names[0])
                    #print(datafile[key_route+'/'+h5key_names[0]])
                    labels = [datafile[key_route+'/'+h5key_names[0]][()]]
                    
                    labels = unique_labels(labels, number_data = number_data)
                    
                    h5key_names= [h5key_names]
                        
                else: #the key is empty and so use numbers
                    labels = list(range(number_data))
                    labels = [unique_labels(labels, number_data = number_data)]
                    #sprint('labels', labels)
                    
                    h5key_names= [h5key_names]

                # print(h5key_names)      
                # print(labels)
                #print(type(labels))
                
                #make the list of labels
                for i in [*range(key_start[0], key_end[0], key_step[0])]:
                    lbl_str = ""
                    for j in range(len(h5key_names[0])):
                        # print(i, j, "'",h5key_names[0],"'")
                        if isinstance(h5key_names[0][j], str) and len(h5key_names[0][j])==0:
                            pass
                        else:
                            lbl_str=lbl_str+os.path.basename(os.path.normpath(key_route+'/'+h5key_names[0][j]))
                            
                        if len(h5key_names[0][j]) == 0: #then there is no name to paste into the file name string 
                            if np.size(labels[j])==1:
                                lbl_str = lbl_str+str(labels[j])
                            else:
                                lbl_str = lbl_str+str(labels[j][i])
                        else:
                            if np.size(labels[j])==1:
                                lbl_str = lbl_str+sep2+str(labels[j])
                            else:
                                lbl_str = lbl_str+sep2+str(labels[j][i])
                        if j != len(h5key_names[0])-1:#np.size(labels[j]):
                            lbl_str = lbl_str+sep1
                            
                    #Make label string and replace all illegal filename characters in label_string.
                    # also replace . with pt in the string
                    if len(key_str) != 0:
                        lbl_str = key_str+sep1+lbl_str
                        
                    IO.licit_filename(lbl_str)
                    
                    out.append([key, i, lbl_str])
                return out[:]
                
            else:
                print('oops')
                stop
        else:
            print('oops, key is wrong')
            stop
    else:
        # step inside the key_list
        
        # make current key
        key_str = key_str+'/'+h5key_list[0]
        # print(key_str)
        
        #make list of values
        if isinstance(h5key_names[0], str) and len(h5key_names[0])==0:
            key_lst = list(datafile[h5key_list[0]].keys())
            key_strings = list(range(len(key_lst)))
            # print('1)', len(key_str))
        else:
            key_lst = list(datafile[h5key_list[0]].keys())
            key_strings = key_lst
        #     print('2)')
        # print(key_lst)
        # print(key_strings)
        #cut to what we want ot iterate over
        # if we are using all the data make sure we run to the end.
        if key_end[0] == -1:
            key_end[0] = len(key_lst)
        
        # if we are looking for more images than there are cut the list
        #key_end[0] = np.min(number_data, key_end[0])
                            
        # if the order is reverse switch the values
        if key_step[0] <= -1:
            key_start[0], key_end[0] = key_end[0], key_start[0]
            
        #populate the list
        #print('start loop', [*range(key_start[0], key_end[0], key_step[0])])
        #for i in range(len(key_lst)):
        for i in [*range(key_start[0], key_end[0]+1, key_step[0])]:
        
            #print(key_lst[i])
            #print(h5key_list[1:])
            
            #use the index to pass either the index of the array or the upper level group index
            if bottom_level=='label':
                if index == 0:
                    index = str(key_lst[i])
                else:
                    index = key_lst[i]
            #print('index', index, key_lst[i])
            
            keys.extend(get_image_keys(datafile, h5key_list[1:], h5key_names[1:], key_start=key_start[1:], key_end=key_end[1:], key_step=key_step[1:], key_route=key_lst[i],bottom_level=bottom_level, index = index, key_str=str(key_strings[i])))
            
            #print(keys)
        
    return keys
    


def get_images(image_list=None, settings_file=None, settings_class=None, image_num=None, debug=False):
    """
    Export images from h5 files for use with Dioptas.

    Parameters
    ----------
    image_list : LIST
        This is a list that contains the file names and h5 location of the images. 
    settings_file : TYPE
        DESCRIPTION.
    export_type : TYPE, optional
        DESCRIPTION. The default is "tiff".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    data : numpy array
        The image(s) data array.

    """
    
    if image_list != None:
        if len(image_list) ==2 and len(image_list[1])==3:
            #then the list is for a single file and needs to be wrapped in another list
            image_list = [image_list]
    elif settings_class != None:
        settings_for_fit = settings_class
        image_list = settings_for_fit.image_list
    elif settings_file != None:
        settings_for_fit = fp.initiate(settings_file)
        image_list = settings_for_fit.image_list
    else:
        raise ValueError(
                "There are no settings or image_list. The fitting cannot proceed until a recognised "
                "settings file or class is present."
            )
    
    if image_num== None:
        image_num = list(range(len(image_list)))
    elif isinstance(image_num, int):
        image_num = [image_num]
    #image_num could also be a list -- in which case leave it alone.
        
    #get image
    for i in range(len(image_num)):
        
        n = image_num[i]
        datafile = h5py.File(image_list[n][0], 'a')
        datakey = image_list[n][1][0]
        data_position_in_key = image_list[n][1][1]
        data_tmp = np.array(datafile[datakey][data_position_in_key])
        #print('size of data (',str(i),'): ', data_tmp.shape)
        
        if len(image_num) == 1:
            data = data_tmp
        else:
            axis=2
            if i==0:
                data = data_tmp
                data = np.expand_dims(data, axis=axis)
            else:
                data = np.append(data, np.expand_dims(data_tmp, axis=axis), axis=axis)        
            #print('size of data (',str(i),'): ', data.shape)
        
    return data


def plot_images(settings_file=None, settings_class=None, image_num=None, export_type="tiff", debug=False):
    """
    Export images from h5 files for use with Dioptas.

    Parameters
    ----------
    settings_file : TYPE
        DESCRIPTION.
    export_type : TYPE, optional
        DESCRIPTION. The default is "tiff".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    #FIX ME: there should be an option to override the settings file and just import a subset of the images -- with more precision than just an image_num
    
    #FIX ME: the plotting of the imags should be done through the data class. 
    
    if settings_class != None:
        settings_for_fit = settings_class
    elif settings_file != None:
        settings_for_fit = fp.initiate(settings_file)
    else:
        raise ValueError(
                "There are no settings. The fitting cannot proceed until a recognised "
                "settings file or class is present."
            )
        
    if image_num== None:
        image_num = list(range(settings_for_fit.image_number))
        
    for i in range(len(image_num)):
        
        n = image_num[i]
        im_data = get_images(settings_class=settings_for_fit, image_num=n)
        plt.imshow(np.log10(im_data))
        # plt.imshow(np.log10(im_data[0,:,:]))
        #plt.title(IO.make_outfile_name(os.path.split(settings_for_fit.image_list[n][0])[1],additional_text=settings_for_fit.image_list[n][1][2]))
        plt.title(IO.title_file_names(settings_for_fit=settings_for_fit, num=i))
        plt.show()
        


def save_images(settings_file=None, settings_class=None, image_num=None, export_type="tif", debug=False):
    """
    Export images from h5 files for use with Dioptas.

    Parameters
    ----------
    settings_file : TYPE
        DESCRIPTION.
    export_type : TYPE, optional
        DESCRIPTION. The default is "tiff".
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    #FIX ME: ther should be an option to override the settings file and just import a subset of the images -- with more precision than just an image_num
    
    if settings_class != None:
        settings_for_fit = settings_class
    elif settings_file != None:
        settings_for_fit = fp.initiate(settings_file)
    else:
        raise ValueError(
                "There are no settings. The fitting cannot proceed until a recognised "
                "settings file or class is present."
            )
        
    if image_num== None:
        image_num = list(range(settings_for_fit.image_number))
        
    for i in range(len(image_num)):
        
        n = image_num[i]
        im_data = get_images(settings_class=settings_for_fit, image_num=n)
        
        
        fnam = IO.make_outfile_name(settings_for_fit.image_list[n][0],additional_text=settings_for_fit.image_list[n][1][2], directory=settings_for_fit.output_directory, extension=export_type)

        #print(im_data.shape)
        cv2.imwrite(fnam, im_data[0,:,:])
        print("Writing:",fnam)




if __name__ == "__main__":
    
    '''
    Yunhui says they use two data keys for position.
    Namely: 
    data_key='/27.1/instrument/positioners/dz1' # this entry contains the positions of the detector or sample or somepart to be used as a label.
    data_key='/17.1/instrument/positioners/dy' # Yunsui said this one also but all the values are 0 -- so maybe not.
    
    I am going to use:
    '/' -- the name of it
    'instrument/positioners/dz1' -- the values inside it.
    
    Thw structure of a h5 file is flexible and Yunhui's case has two indicies to iterate over. 
    therefore I am going to make the itertion possible over any number of indicies. By calling a loop rather than assuming how many levels are needed.
    
    The output files will need to be named according to:
        - their positions in the loop
        - an index name or 
        - another index.
    and will need to be named for the number of indicices iterated over. 
    Or a single value that is taken for all indicies. 
    
    We will need to check that all the names are unique.
    
    It might also not be required to process all the data. Therfore we need a:
        - start, end and step for the indicies, or
        - a list of values to process.
    We can copy the normal file pattern.
    But the values have to be a list the same length as the indicies. 
    
    
    
    '''
    x = 1
    y=2
    
    if x==0:
        # IN004
        fname = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/IN004/IN004_0001/IN004_0001.h5'
        # fname = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/IN004/IN004_0001/scan0025/p3_0000.h5'
        # outname = './IN004_0001/scan0025/p3_0000'
        # data_key = "/entry_0000/measurement/data"
    elif x==1:
        # RR006
        fname = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/RR006/RR006_0001/RR006_0001.h5'
        # fname = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/RR006/RR006_0001/scan0007/p3_0000.h5'
        # outname = './RR006_0001/scan0007/p3_0000'
        # data_key = "/entry_0000/measurement/data"
        
        # h5_key_list  = ['/', 'measurement/p3']
        # if 1: #iterate
        #     h5_key_names = ['/', 'instrument/positioners/dy']
        #     h5_key_names = ['/', ['instrument/positioners/dz1','instrument/positioners/dy']]
        #     h5_key_start = [0, 0]
        #     k5_hey_end   = [-1,-1]
        #     h5_key_step  = [1,-1]
        #     h5_data      = 'iterate' # the other options would be "sum". "iterate" only applies to the bottom level.
        # else: #sum
        #     h5_key_names = ['/', 'title']
        #     h5_key_start = [0, 0]
        #     k5_hey_end   = [-1,-1]
        #     h5_key_step  = [1,-1]
        #     h5_data      = 'sum' # the other options would be "sum". "iterate" only applies to the bottom level.
        
    elif x==2:
        # trip_005
        fname = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/trip005/trip005_0001/trip005_0001.h5'
        # fname = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/trip005/trip005_0001/scan0048/p3_0000.h5'
        # outname = './trip005/trip005_0001/scan0048/p3_0000'
        # data_key = "/entry_0000/measurement/data"
    elif x==3:
        datafile_directory = '/Volumes/External/DATA_other_people/UCLengineers_YunhuiChen/data3_ID31/trip005/trip005_0001/'
        datafile_Basename  = 'trip005_0001'
        datafile_Ending    = '.h5'
        #datafile_StartNum  = 1
        #datafile_EndNum    = 10
        datafile_NumDigit  = 0
        datafile_h5=True
        y = 2
    elif x==4:
        datafile_directory = './'
        datafile_Basename  = 'BCC1_2GPa_10s_001_'
        datafile_Ending    = '.tif'
        datafile_StartNum  = 1
        datafile_EndNum    = 10
        datafile_NumDigit  = 5
        datafile_h5=False
        
    else:
        print('no data')
        stop
        
    
    
    h5_key_list  = ['/', 'measurement/p3']
    if y==1: #iterate
        h5_key_names = ['/', 'instrument/positioners/dy']
        h5_key_names = ['/', ['instrument/positioners/dz1','instrument/positioners/dy']]
        h5_key_start = [0, 0]
        k5_hey_end   = [-1,-1]
        h5_key_step  = [1,-1]
        h5_data      = 'iterate' # the other options would be "sum". "iterate" only applies to the bottom level.
    else: #sum
        h5_key_names = ['/', 'title']
        h5_key_start = [0, 0]
        k5_hey_end   = [-1,-1]
        h5_key_step  = [1,-1]
        h5_data      = 'sum' # the other options would be "sum". "iterate" only applies to the bottom level.
    
    
    
        
    
        
    datafile = h5py.File(fname, 'a')
    
    
    # my_list = get_image_keys(datafile, h5_key_list, key_start=h5_key_start, key_end=k5_hey_end, key_step=h5_key_step, bottom_level=h5_data)
    # print(my_list)
    
    my_labels = get_image_keys(datafile, h5_key_list, h5_key_names, key_start=h5_key_start, key_end=k5_hey_end, key_step=h5_key_step, bottom_level=h5_data)
    
    print(my_labels)
