# Example 1 -- BCC1_2GPa_100s_001; peaks with no data in 
# 
# This is an empty case; see that everything runs when there are no peaks in the data.
# It only contains required options. 

# properties of the data files.
datafile_directory = "."
datafile_Basename  = 'BCC1_2GPa_10s_001_'
datafile_Ending    = '.tif'
datafile_StartNum  = 1
datafile_EndNum    = 10
datafile_NumDigit  = 5
datafile_Step      = 1

# Calibration and masking.
Calib_type     = "Dioptas"
Calib_detector = 'Pilatus1M'
Calib_data     = "CeO2_Pil207_E30_2Nov2016_001.tif"
Calib_param    = "CeO2_cal_Dioptas.poni"
Calib_mask     = "DiffractionMask_Dioptas.mask"
Calib_pixels = 172

#fit settings
fit_min_peak_intensity = "std"
fit_min_data_intensity = 30

#Output settings
Output_directory   = 'results'
Output_type        = ['CoefficientTable', 'FitMovie']


# define ranges and peaks
fit_orders = [
    {
        "range": [11.95, 12.2],
        "background": [0, 0],
        "peak": [{
            "d-space": 2, 
            "height": 12, 
            "profile": 0, 
            "width": 0,
            "symmetry": 2,
        }], 
    },
    {
        "range": [16.6, 17.1],
        "background": [0, 0],
        "peak": [{
            "d-space": 3,
            "height": 16,
            "profile": 0,
            "width": 0,
            "symmetry": 2,
        }],
    },
    {
        "range": [20.7, 20.95],
        "background": [0, 0],
        "peak": [{
            "d-space": 2,
            "height": 16,
            "profile": 1,
            "width": 0,
            "symmetry": 2,
        }],
    },
    {
        "range": [26.4, 27.0],
        "background": [2, 2],
        "peak": [{
            "d-space": 1,
            "height": 1,
            "profile": 0,
            "width": 0,
            "symmetry": 2,
        }],
    }
]








# # Input parameters for Example 1 -- BCC1_2GPa_100s_001. 

# # properties of the data files.
# datafile_directory = './'
# datafile_Basename  = 'BCC1_2GPa_10s_001_'
# datafile_Ending    = '.tif'
# datafile_StartNum  = 1
# datafile_EndNum    = 10
# datafile_NumDigit  = 5

# # Calibration and masking.
# Calib_type     = "Dioptas"
# Calib_detector = 'Pilatus1M'
# Calib_data     = datafile_directory + 'CeO2_Pil207_E30_2Nov2016_001.tif'
# Calib_param    = datafile_directory + 'CeO2_cal_Dioptas.poni'
# Calib_mask     = datafile_directory + 'DiffractionMask_Dioptas.mask'
# Calib_pixels   = 172

# # Fitting properties for peaks.
# fit_bounds = {
#               "background": ['min', 'max'],
#               "d-space":    ['min', 'max'],
#               "height":     [ 0,    '2*max'],
#               "profile":    [ 0,     1],
#               "width":      [ 'range/(ndata)',  'range/2'],
#               }
# fit_orders = [
#         {
#           "range": [11.95, 12.2],
#           "background": [2, 0],
#           #"background-type": "spline-cubic",
#           "peak": [{
#               "d-space": 3,
#               #"d-space-type": "spline-cubic",
#               "height": 1,
#               #"height-type": "fourier",
#               "profile": 0,
#               "profile_fixed": 0.1,
#               "width": 0,
#               "width_fixed": 0.035,
#               "symmetry": 2
#             }]
#         },
#         {
#           "range": [16.6, 17.1],
#           "background": [0, 0],
#           "peak": [{
#               "d-space": 3,
#               "height": 8,
#               "profile": 0,
#               "profile_fixed": 1,
#               "width": 0,
#               "symmetry": 2
#             }],
#         },
#         {
#           "range": [20.7, 20.95],
#           "background": [0, 0],
#           "peak": [{
#               "d-space": 3,
#               "height": 8,
#               "profile": 1,
#               "profile_fixed": [0.5, 0.25, 0.25],
#               "width": 0,
#               "symmetry": 2
#             }],
#         },
#         {
#            "range": [26.4, 27.0],
#            "background": [2, 2],
#            "peak": [{
#                "d-space": 2,
#                "height": 1,
#                "profile": 0,
#                "profile_fixed": 0.5,
#                "width": 0,
#                "symmetry": 2
#              }],
#            "imax":27,
#         },
#     ]


# #Number of bins for initial fitting.
# AziBins = 90

# #Output settings
# Output_directory   = './results/'
# Output_type        = ['Polydefix', 'DifferentialStrain', 'FitMovie']
# Output_NumAziWrite = 90
# phase              = 'Fe-BCC'
# Output_ElasticProperties = 'Properties_Fe-BCC.txt'