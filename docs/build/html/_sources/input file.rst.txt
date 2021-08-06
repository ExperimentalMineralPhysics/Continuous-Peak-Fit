.. Continuous-Peak-Fit documentation master file


=====================================
Input file structure (long)
=====================================

.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit


The input file for continuous peak fit is a python file that lists and defines most of the parameters needed to run and package on a data set. The settings are as follows:



Data file definitions
=====================================

The file names are assumed to consist of a root or basename (e.g. ```BCC1_...```, an incrementing number that is systematically formatted and a file 'ending' (e.g. .tif). 

The options are listed here: 

===========================  ==============   ================================
Parameter                    Compulsary?      Description
===========================  ==============   ================================
``datafile_directory``       yes              The location of the directory containing the data files.
``datafile_Basename``        yes              The name of the file before the incrementing numbers. Can be an empty string
``datafile_NumDigit``        no               number of digits in the file number. If absent assumed the numbers are not zero padded. This applies to both style a and b below. 
``datafile_Ending``          yes              The ending of the file after the incrementing numbers. Can be an empty string.
``datafile_StartNum``        style a, yes     number of the first diffraction file
``datafile_EndNum``          style a, yes     number of the last diffraction file
``datafile_Step``            style a, no      increment for file numbers. If absent assumed to be 1 (if StartNum > EndNum) or -1 (if StartNum < EndNum) 
``datafile_Files``           style b, yes.    a list of the incrementing part of the file names. It is assumed to be numberic and is zero padded if 'datafile_NumDigit' is present.
===========================  ==============   ================================


Example file name structures:
-----------------------------

- Example 1: style a, from [LINK to Example 1] :

 .. code-block:: python

    datafile_directory  = './Example1-Fe/'     
    datafile_Basename   = 'BCC1_2GPa_10s_001_' 
    datafile_Ending     = '.tif'               
    datafile_StartNum   = 1                    
    datafile_EndNum     = 10                  
    datafile_NumDigit   = 5                   

 - Makes file names::

    ./Example1-Fe/BCC1_2GPa_10s_001_00001.tif,
    ./Example1-Fe/BCC1_2GPa_10s_001_00002.tif,
    ...,
    ./Example1-Fe/BCC1_2GPa_10s_001_00008.tif,
    ./Example1-Fe/BCC1_2GPa_10s_001_00009.tif,
    ./Example1-Fe/BCC1_2GPa_10s_001_00010.tif


- Example 2, style a:

 .. code-block:: python

    datafile_directory  = './'                
    datafile_Basename   = ''        
    datafile_Ending     = 'diffraction.tif' 
    datafile_StartNum   = 15        
    datafile_EndNum     = 3      
    datafile_Step.      = -3    

 - Makes file names::

    ./15_diffraction.tif,
    ./12_diffraction.tif,
    ./9_diffraction.tif,
    ./6_diffraction.tif,
    ./3_diffraction.tif


- Example 3, style a:

 .. code-block:: python

    datafile_directory = 'F:\\remote_server_location\'
    datafile_Basename  = 'M1234.' 
    datafile_Ending    = ''  
    datafile_StartNum  = 1   
    datafile_EndNum    = 30    
    datafile_NumDigit  = 4     

 - Makes file names::

    F:\\remote_server_location\M1234.0001
    F:\\remote_server_location\M1234.0002
    F:\\remote_server_location\M1234.0003
    ...
    F:\\remote_server_location\M1234.0028
    F:\\remote_server_location\M1234.0029
    F:\\remote_server_location\M1234.0030



- Example 4, style b:

 .. code-block:: python

    datafile_directory = '../another/directory/'  
    datafile_Basename  = 'D654_'   
    datafile_Ending    = '.nxs'   
    datafile_Files.    = [6,7, 9,10, 12,13, 115]   
    datafile_NumDigit  = 4     

 - Makes file names::

    ../another/directory/D654_0006.nxs
    ../another/directory/D654_0007.nxs
    ../another/directory/D654_0009.nxs
    ../another/directory/D654_0010.nxs
    ../another/directory/D654_0012.nxs
    ../another/directory/D654_0013.nxs
    ../another/directory/D654_0115.nxs


Calibarion definitions
=====================================

A calibration type and file are required. Other requirements are possible depending on the detector type. ```cpf.XRD_FitPattern.initiate``` (LINK) will inform you if needed definitions are missing. 


==================       =============   ================================
Parameter                Required?            Description
==================       =============   ================================
``Calib_type``           yes             tells the souce for the calibration. Accepted values are "Dioptas" and "MED". 
``Calib_param``          yes             This is the file containing the calibration parameters for the detector. If it is not in the same directory as the imput file then needs to contain a path.
``Calib_mask``           no              Is the mask file for the detector. If it is not in the same directory as the imput file then needs to contain a path.
``Calib_detector``       no              name of the detector. Not used for Dioptas calibration. Might be needed for MED and others though. 
``Calib_data``           no              the name of the data file with the calibration data in. Used in debugging to check the calibration is imported propertly.
``Calib_pixels``         no              needed by GSAS-II calibrations, which currently dont work!!
==================       =============   ================================



Peak definitions and limits
=====================================

Each subpattern is defined by number of parameters contained within a python dictionary. The dictionary for each subpattern is contained with in a list called 'fit_orders'. e.g.:

 .. code-block:: python

  fit_orders = [
         {
           "range": [[3.5,3.8]],
           "background": [1,0],
           "peak": [{
               "phase": "Steel",
               "hkl": '111',
               "d-space": 2,
               "height": 12,
               "profile": 0,
               #"profile_fixed": 1,
               "width": 0,
               #"symmetry": 2
             }],
           },
         {
           "range": [[4.1,4.4]],
           "background": [1,0],
           "peak": [{
               "phase": "Steel",
               "hkl": '200',
               "d-space": 2,
               "height": 12,
               "profile": 0,
               #"profile_fixed": 1,
               "width": 0,
               #"symmetry": 2
             }],
         },
         ]


With in each subpattern dictionary the following parameters are possible or needed: 

==========================  ==========================   ================================
Parameter                   Required?                    Description
==========================  ==========================   ================================
``range``                   yes                          Defines the maximum and minimum of the range to fit with peak in. For angle dispersive diffraction data this is the maximum and minimum two theta angle (in degrees), for energy disperdive diffraction it is the maximum and minimum energy (in keV). The format is ``[[min, max]]``. 
``background``              yes                          This is a list of numbers, see LINK. The length of the list is the order 
``Imax``                    no                           Crops the maximum intensity within the range at this value. Formatted as a single number. 
``Imin``                    no                           Crops the minimum intensity within the range at this value. Formatted as a single number. 
``PeakPositionSelection``   only for multiple peaks
``peak``                    yes                          List of dictionaries. See 
==========================  ==========================   ================================


Peak definition
------------------------------


Output settings
=====================================


blah balh balh


Optional Extras. 
=====================================

Azi_bins


Limits




