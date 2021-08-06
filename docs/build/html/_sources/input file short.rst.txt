.. Continuous-Peak-Fit documentation master file


=====================================
Input file structure (short)
=====================================

.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit


The input file for continuous peak fit is a python file that lists and defines most of the parameters needed to run and package on a data set. 
The input file is a single python file with the following required settings:


.. table:: input table

+-----------------------------------------------------+--------------------------------------------------------------------+
| .. code-block:: python                              | Define the properties of the input files.                          |
|                                                     |                                                                    |
|   # data drive.                                     | Rather than the start and end numbers it is possible to call       |
|   import os                                         | 'datafile_Files' which is a list of the files names without        |
|   drive = "/Directory_string/",                     | endings.                                                           |
|   # properties of the data files.                   |                                                                    |
|   datafile_directory = drive + 'Another_dir_name/', | datafile_Step is optional -- if absent is assumed to be +/-1.      |
|   datafile_Basename = "Pixium_",                    |                                                                    |
|   datafile_Ending = ".tif",                         | For full description see: :ref:`datafile_definitions`              |
|   datafile_StartNum = 1,                            |                                                                    |
|   datafile_EndNum = 1,                              |                                                                    |
|   datafile_NumDigit = 5,                            |                                                                    |
|   datafile_Step = 2,                                |                                                                    |
|   #datafile_Files = [3,5,7,8,9],                    |                                                                    |
+-----------------------------------------------------+--------------------------------------------------------------------+
| .. code-block:: python                              |  Detector Mask and Calibration.                                    |
|                                                     |                                                                    |
|   # Calibration and masking.                        | For full description see: :ref:`calibration_definitions`           |
|   Calib_type = "Dioptas",                           |                                                                    |
|   Calib_detector = 'unknown',                       |                                                                    |
|   Calib_data = drive + '/pixi_00001.tif',           |                                                                    |
|   Calib_param = 'DLS_CeO2_1200mm.poni',             |                                                                    |
|   Calib_mask = 'DLS_CeO2_1200mm.mask',              |                                                                    |
|   Calib_pixels = 296,                               |                                                                    |
+-----------------------------------------------------+--------------------------------------------------------------------+
| .. code-block:: python                              |                                                                    |
|                                                     |                                                                    |
|   # Fitting properties for peaks.,                  | `fit_orders` is a list containing all the parameters needed to fit |
|   fit_orders = [                                    | each subpattern in the data. Each subpattern is described by a set |
|      {"range": [[3.05, 3.25]],                      | of parameters which are contained within a python dictionary.      |
|       "background": [0,0],                          |                                                                    |
|       "peak": [{"phase":"Ti64",                     | Each subpattern covers at least one diffraction peak, the          |
|                 "hkl": '100'                        | numbers for each of the peak properties are the orders of the      |
|                 "d-space": 2,                       | fourier series or splines for that parameter.                      |
|                 "height": 6,                        | The number of coefficients is 2*n+1.                               |
|                 "profile": 0,                       |                                                                    |
|                 #"profile_fixed": 1,                | 'profile_fixed' will fix the pseudo-voigt peak profile to          |
|                 "width": 0,                         | the value in profile.                                              |
|                 "symmetry": 2}]                     |                                                                    |
|      },                                             | Range is in two theta.                                             |
|      {"range": [[3.5, 3.7]],                        |                                                                    |
|       "background": [0,0],                          | For additional peaks just add more dictionaries to the list.       |
|       "peak": [{"phase": "Ti64beta",                |                                                                    |
|                 "hkl": '110',                       | For multiple peaks the peak list has multiple entries and          |
|                 "d-space": 2,                       | ``PeakPositionSelection`` is required.                             |
|                 "height": 6,                        | [[ peak, azimuth, two theta],                                      |
|                 "profile": 0,                       | ...                                                                |
|                 #"profile_fixed": 1,                | ]                                                                  |
|                 "width": 0,                         | A routine exists for making this array, see :ref:`initial_guess`   |
|                 "symmetry": 2,},                    |                                                                    |
|                {"phase": "Ti64alpha",               |                                                                    |
|                 "hkl": '101',                       | A full descrption of the paramters for each subpattern is          |
|                 "d-space": 2,                       | provided here: :doc:`subpattern structure`                         |
|                 "height": 6,                        |                                                                    |
|                 "profile": 0,                       |                                                                    |
|                 #"profile_fixed": 1,                |                                                                    |
|                 "width": 0,                         |                                                                    |
|                 "symmetry": 2,}],                   |                                                                    |
|       "PeakPositionSelection":                      |                                                                    |
|             [[1,-120.5,3.54],                       |                                                                    |
|              [1,-58.997,3.54],                      |                                                                    |
|              [1,59.289,3.54],                       |                                                                    |
|              [1,23.187,3.54],                       |                                                                    |
|              [1,23.212,3.54],                       |                                                                    |
|              [1,23.158,3.54],                       |                                                                    |
|              [1,123.246,3.54],                      |                                                                    |
|              [2,-120.5,3.59],                       |                                                                    |
|              [2,-58.997,3.59],                      |                                                                    |
|              [2,59.289,3.59],                       |                                                                    |
|              [2,23.187,3.59]'                       |                                                                    |
|              [2,23.212,3.59],                       |                                                                    |
|              [2,23.158,3.59],                       |                                                                    |
|              [2,123.246,3.59]],                     |                                                                    |
|      },                                             |                                                                    |
|      ]                                              |                                                                    |
+-----------------------------------------------------+--------------------------------------------------------------------+
| .. code-block:: python                              | Outputs.  The types are 'DifferentialStrain', 'Polydefix',         |
|                                                     | 'multifit' and 'PolyDefixED'.                                      |
|   #Output settings                                  |                                                                    |
|   Output_directory = './'                           | No further options are needed for differential strain.             |
|   Output_type = 'DifferentialStrain'                | But an elastic strain tensor should be provided for the others.    |
|                                                     |                                                                    |
|                                                     | For full description see: :ref:`output_definitions`                |
+-----------------------------------------------------+--------------------------------------------------------------------+

