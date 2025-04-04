# Example 1 -- BCC1_2GPa_100s_001; fix series values between peaks. 
# 
# This example demonstrates foring the coeficcients in two peaks in the same range to have the same values.
# It is only applicable to ranges with more than 1 peak.
#
# N.B. This is not the best example date for this feature -- it would be better if there was an example where the peaks overlap and are not distinct.

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

# Fitting properties for peaks.
AziBins = 60 # Number of bins for initial fitting.
fit_bounds = {
              "background": ['min', 'max'],
              "d-space":    ['min', 'max'],
              "height":     [ 0,    '1.2*max'],
              "profile":    [ 0,     1],
              "width":      [ 'range/(ndata)',  'range/2'],
              }

#Output settings
Output_directory   = 'results'
Output_type        = ['Polydefix', 'FitMovie', 'CoefficientTable']
# optional settings for 'Polydefix'.
Output_NumAziWrite = 90 
phase              = 'Fe-BCC' 
Output_ElasticProperties = 'Properties_Fe-BCC.txt'

# define ranges and peaks
fit_orders = [
    {
        "range": [[11.0, 11.7]],
        "background": [2, 0],
        "background-type": "spline-cubic",
        "peak": [{
            "phase": "Other",
            "hkl": '000',
            "d-space": 2,
            "d-space-type": "fourier",
            "height": 1,
            "profile": 0,
            "width": 0,
            "symmetry": 2
        }, {
              "phase": "Fe-BCC",
              "hkl": 110,
              "d-space": [3, "1stOrderEqualPeak0"],
              "height": 16,
              "height-type": "spline-cubic",
              "profile": [0, "0thOrderEqualPeak0"],
              #"profile_fixed": 1,
              "width": [0, "allEqualPeak0"],
              # "width_fixed": 0.035,
              "symmetry": 2
            }],
        "PeakPositionSelection": [[1, -170.5,   11.19],
                                  [1, -150.5,   11.19],
                                  [1, -120.5,   11.19],
                                  [1,  -58.99,  11.19],
                                  [1,   59.289, 11.19],
                                  [1,   23.187, 11.19],
                                  [1,   23.212, 11.19],
                                  [1,   53.158, 11.19],
                                  [1,   73.246, 11.19],
                                  [1,   93.246, 11.19],
                                  [1,  123.246, 11.19],
                                  [2, -170,     11.54],
                                  [2, -150,     11.54],
                                  [2, -110,     11.54],
                                  [2,  -90,     11.54],
                                  [2,  -45,     11.54],
                                  [2,   -0,     11.54],
                                  [2,   45,     11.54],
                                  [2,   90,     11.54],
                                  [2,  120,     11.54],
                                  [2,  150,     11.54],
                                  [2,  179,     11.54]]
    },
    #     {
    #     "range": [16.1, 20.27],
    #     "imax":27,
    #     "background": [0, 0],
    #     "peak": [{
    #         "phase": "Fe-BCC",
    #         "hkl": 200,
    #         "d-space": 3,
    #         "height-type": "spline_linear",
    #         "height": 8,
    #         "profile": 0,
    #         "profile_fixed": 1,
    #         "width": 0,
    #         "symmetry": 2
    #         },{
    #         "phase": "Fe-BCC",
    #         "hkl": 211,
    #         "d-space": 3,
    #         "height": 8,
    #         "height-type": "quadratic",
    #         "profile": 1,
    #         "profile_fixed": [0.5, 0.25, 0.25],
    #         "width": 0,
    #         "symmetry": 2
    #         }],
    # },
]

# fit_orders = [
#         {
#           "range": [[10.8, 11.67]],
#           "background": [2, 0],
#           #"background-type": "spline-cubic",
#           "peak": [{
#               "phase": "Other",
#               "hkl": '000',
#               "d-space": 3,
#               #"d-space-type": "spline-cubic",
#               "height": 1,
#               #"height-type": "fourier",
#               "profile": 0,
#               "profile_fixed": 0.1,
#               "width": 0,
#               # "width_fixed": 0.035,
#               "symmetry": 2
#             }, {
#               "phase": "Fe-BCC",
#               "hkl": 110,
#               "d-space": [3, "1stOrderEqualPeak0"],
#               "height": 16,
#               "height-type": "spline-cubic",
#               "profile": [0, "0thOrderEqualPeak0"],
#               #"profile_fixed": 1,
#               "width": [0, "allEqualPeak0"],
#               # "width_fixed": 0.035,
#               "symmetry": 2
#             }],
#           "PeakPositionSelection": [[1, -170.5,   11.19],
#                                     [1, -150.5,   11.19],
#                                     [1, -120.5,   11.19],
#                                     [1,  -58.99,  11.19],
#                                     [1,   59.289, 11.19],
#                                     [1,   23.187, 11.19],
#                                     [1,   23.212, 11.19],
#                                     [1,   53.158, 11.19],
#                                     [1,   73.246, 11.19],
#                                     [1,   93.246, 11.19],
#                                     [1,  123.246, 11.19],
#                                     [2, -170,     11.54],
#                                     [2, -150,     11.54],
#                                     [2, -110,     11.54],
#                                     [2,  -90,     11.54],
#                                     [2,  -45,     11.54],
#                                     [2,   -0,     11.54],
#                                     [2,   45,     11.54],
#                                     [2,   90,     11.54],
#                                     [2,  120,     11.54],
#                                     [2,  150,     11.54],
#                                     [2,  179,     11.54]]
#         },
#         {
#           "range": [[16.1, 16.6]],
#           "background": [0, 0],
#           "peak": [{
#               "phase": "Fe-BCC",
#               "hkl": 200,
#               "d-space": [0,2],
#               "height": 8,
#               "profile": 0,
#               "profile_fixed": 1,
#               "width": 0,
#               "symmetry": 2
#             }],
#         },
#     ]