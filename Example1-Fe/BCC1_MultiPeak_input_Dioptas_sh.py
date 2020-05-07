

# Input parameters for BCC1_2GPa_100s_001.
#
# It is an attempt to collate all the information needed in the fitting routines.
# A 'proper' input file type can be formulated later.

# properties of the data files.
datafile_directory = './Example1-Fe/'
datafile_Basename  = 'BCC1_2GPa_10s_001_'
datafile_Ending    = '.tif'
datafile_StartNum  = 1
datafile_EndNum    = 10
datafile_NumDigit  = 5

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
         "background": [0,0],
         "peak": [{
             "phase": "Other",
             "hkl": '000',
             "d-space": 3,
             "height": 1,
             "profile": 0,
             #"profile_fixed": 1,
             "width": 0,
             "symmetry": 2
           },{
             "phase": "BCC-Fe",
             "hkl": 110,
             "d-space": 3,
             "height": 16,
             "profile": 0,
             #"profile_fixed": 1,
             "width": 0,
             "symmetry": 2
           }],
         "range": [[10.8, 11.67]],
         "PeakPositionSelection": [[1,-170.5,11.19],[1,-150.5,11.19],[1,-120.5,11.19],[1,-58.997,11.19],[1,59.289,11.19],[1,23.187,11.19],[1,23.212,11.19], [1,53.158,11.19], [1,73.246,11.19], [1,93.246,11.19], [1,123.246,11.19], [2,-170, 11.54], [2,-150, 11.54], [2,-110, 11.54], [2,-90, 11.54], [2,-45, 11.54], [2,-0, 11.54], [2,45, 11.54], [2,90, 11.54], [2,120, 11.54], [2,150, 11.54], [2,179, 11.54]]
       },
       {
         "background": [0,0],
         "peak": [{
             "phase": "BCC-Fe",
             "hkl": 200,
             "d-space": 3,
             "height": 8,
             "profile": 0,
             #"profile_fixed": 1,
             "width": 0,
             "symmetry": 2
           }],
         "range": [[16.1,16.6]]
       },
       {
         "background": [0,0],
         "peak": [{
             "phase": "BCC-Fe",
             "hkl": 211,
             "d-space": 3,
             "height": 8,
             "profile": 0,
             #"profile_fixed": 1,
             "width": 0,
             "symmetry": 2
           }],
         "range": [[19.8, 20.27]]
       },
      {
          "background": [1],
          "peak": [{
             "phase": "BCC-Fe",
             "hkl": 220,
              "d-space": 3,
              "height": 8,
              "profile": 0,
              "profile_fixed": 0.5,
              "width": 0,
              "symmetry": 2
            }],
#          "PeakPositionSelection": [[1,-120.5,23.226],[1,-58.997,23.186],[1,59.289,23.208],[1,23.187,120.893],[1,23.212,6.367], [1,23.158,90.025], [1,23.246,25.729]]
          "range": [[23.1, 23.5]]
        },
       {
         "background": [1,1],
         "peak": [{
             "phase": "BCC-Fe",
             "hkl": 310,
             "d-space": 2,
             "height": 1,
             "profile": 0,
             "profile_fixed": 0.5,
             "width": 0,
             "symmetry": 2
           }],
         "range": [[25.7, 26.4]]
       },
    ]

#Number of bins for initial fitting.
AziBins = 90


#Output settings
Output_directory   = './'
Output_type        = 'MultiFit'
Output_NumAziWrite = 90
