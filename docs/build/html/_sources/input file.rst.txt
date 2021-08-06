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

File name construction does not add underscores or fullstops between the components. A fullstop must be contained within the strings used to construct the file name for the file extension to be present. 

Whilst both ``datafile_Basename`` and ``datafile_Ending`` are required, it is permitted for them to be empty strings. See examples 2 and 3 below. 


Example file name structures:
-----------------------------

- Example 1: style a, from [LINK to Example 1]:

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


- Example 3, style b:

 .. code-block:: python

    datafile_directory = 'F:\\remote_server_location\'
    datafile_Basename  = 'M1234.' 
    datafile_Ending    = ''  
    datafile_Files.    = [6,7, 9,10, 12,13, 115]
    datafile_NumDigit  = 4     

 - Makes file names::

    F:\\remote_server_location\M1234.0006
    F:\\remote_server_location\M1234.0007
    F:\\remote_server_location\M1234.0009
    F:\\remote_server_location\M1234.0010
    F:\\remote_server_location\M1234.0012
    F:\\remote_server_location\M1234.0013
    F:\\remote_server_location\M1234.0115



Calibration definitions
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



Subpattern definitions
=====================================

A full description of the structure required for each subpattern is given here: :doc:`subpattern structure`



Output settings
=====================================




blah balh balh


Additional (Optional) Settings 
=====================================

Azimuthal Bins
-------------------------------------
``AziBins`` is the number of bins used in the inital fitting of the data (LINK). For angle dispersive data the default number of bins is 90. This value gives 4 degree azimuthal bins for complete Debye-Scherer rings. For partial diffraction rings a smaller number maybe a better chice. This is set in input file by:

 .. code-block:: python

  AziBins = 45 



Limits
-------------------------------------
The paramters in the fits are bounded (see LINK for details). The default limits are:

 .. code-block:: python

  limits = []


These can be changed by .... 









