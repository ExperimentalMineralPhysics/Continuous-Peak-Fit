# Input parameters for BCC1_2GPa_100s_001. 
#
# It is an attempt to collate all the information needed in the fitting routines.
# A 'proper' input file type can be formulated later.

# properties of the data files.
datafile_directory = 'Example1-Fe/'
datafile_Basename  = 'BCC1_2GPa_10s_001_'
datafile_Ending    = '.tif'
datafile_StartNum  = 1
datafile_EndNum    = 10
datafile_NumDigit  = 5
datafile_Step      = 3

'''
# Calibration and masking.
Calib_type   = 'GSAS-II'
Calib_file   = './Calibration2018Oct.imctrl'
Calib_mask   = './Diffraction.immask'
Calib_pixels = 172
'''

# Calibration and masking.
Calib_type   = "Dioptas"
Calib_detector = 'Pilatus1M'
Calib_data     = datafile_directory + 'CeO2_Pil207_E30_2Nov2016_001.tif'
Calib_param    = datafile_directory + 'CeO2_cal_Dioptas.poni'
Calib_mask     = datafile_directory + 'DiffractionMask_Dioptas.mask'
Calib_pixels = 172


# Fitting properties for peaks.
fit_orders = [
      {
        "background": [0],
        "peak": [{
            "d-space": 2, 
            "height": 6, 
            "profile": 0, 
            "width": 1,
"symmetry": 1
          }], 
        "range": [[11.70, 12.03]]
      },
      {
        "background": [0], 
        "peak": [{
            "d-space": 2, 
            "height": 8, 
            "profile": 0, 
            "width": 0,
"symmetry": 1
          }], 
        "range": [[16.6,17.1]]
      },
      {
        "background": [0], 
        "peak": [{
            "d-space": 2, 
            "height": 6, 
            "profile": 0, 
            "width": 0,
"symmetry": 1
          }], 
        "range": [[20.4, 20.9]]
      }, 
    ]

# Number of bins for initial fitting.
AziBins = 90


# Output settings
Output_directory   = 'Example1-Fe/results/'
Output_type        = 'MultiFit'
Output_NumAziWrite = 90
