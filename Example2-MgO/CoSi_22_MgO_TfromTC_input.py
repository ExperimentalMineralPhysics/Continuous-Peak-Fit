# Input parameters for BCC1_2GPa_10s_001.
#
# It is an attempt to collate all the information needed in the fitting routines.
# A "proper" input file type can be formuated later.

# data drive.
import os 
if os.path.isdir("/Volumes/My Passport"): #portable drive
    drive = "/Volumes/My Passport/"
elif os.path.isdir("/media/simon/My Passport"): #analysis on Linux box
    drive = "/media/simon/My Passport/"
elif os.path.isdir("/volumes/external"): #desk data store
    drive = "/Volumes/External/"
else:
    drive = ""

# properties of the data files.
datafile_directory = "/Users/simon/Documents/WORK/Experiment_Analysis/RelStrength-ol-rw/CoSi_data/COSI_22"
datafile_Basename  = "COSI_022_"
datafile_Ending    = ".med"
datafile_Files     = [6, 9, 11, 13, 15, 17, 
                      19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59,
                      61, 63, 65, 67, #69, 71, 73
                      ]
#datafile_StartNum = 10
#datafile_EndNum =  101
datafile_NumDigit  = 4

# Calibration and masking.
Calib_type     = 'Med'
Calib_detector = '10-element'
Calib_data     = "/Users/simon/Documents/WORK/Experiment_Analysis/RelStrength-ol-rw/CoSi_data/COSI_21/COSI_021_0001.med"
Calib_param    = '/Users/simon/Documents/WORK/Experiment_Analysis/RelStrength-ol-rw/CoSi_data/COSI_21/COSI_021_0001.med'
Calib_mask     = [5,10]
Calib_pixels   = 72

# Fitting properties for peaks.
fit_orders = [
        {
          "range": [[43.5, 46.0]],
          "background": [3],
          "peak": [{
              "phase": "MgO",
              "hkl": 111,
              "d-space": [0,2],
              "height": 3,
              "profile": 0,
              #"profile_fixed": 1,
              "width": 0,
              "symmetry": 2
            }],
        },
        {
          "range": [[51.8, 53.3]],
          "background": [3],
          "peak": [{
              "phase": "MgO",
              "hkl": 200,
              "d-space": [0,2],
              "height": 3,
              "profile": 0,
              #"profile_fixed": 1,
              "width": 0,
              "symmetry": 2
            }],
        },
        {
          "range": [[73.0, 75.5]],
          "background": [3],
          "peak": [{
              "phase": "MgO",
              "hkl": 220,
              "d-space": [0,2],
              "height": 3,
              "profile": 0,
              #"profile_fixed": 1,
              "width": 0,
              "symmetry": 2
            }],
        },
        # {
        #     "range": [[83.0, 95]],
        #     "background": [3],
        #     "peak": [{
        #        "phase": "Other",
        #         "d-space": 2,
        #         "height": 3, 
        #         "profile": 0,
        #         #"profile_fixed": 0.5,
        #         "width": 0,
        #         "symmetry": 2
        #       },{
        #        "phase": "MgO",
        #        "hkl": 311,
        #         "d-space": 2,
        #         "height": 3, 
        #         "profile": 0,
        #         #"profile_fixed": 0.5,
        #         "width": 0,
        #         "symmetry": 2
        #       },{
        #        "phase": "Other",
        #         "d-space": 2,
        #         "height": 3, 
        #         "profile": 0,
        #         #"profile_fixed": 0.5,
        #         "width": 0,
        #         "symmetry": 2
        #       }],
        #     "PeakPositionSelection": [[1,   0.0,  84.7],
        #                               [1,  22.5,  84.7],
        #                               [1,  45.0,  84.7],
        #                               [1,  67.5,  84.7],
        #                               [1,  90.0,  84.7],
        #                               [1, 112.5,  84.7],
        #                               [1, 135.0,  84.7],
        #                               [1, 157.5,  84.7],
        #                               [1, 180.0,  84.7],
        #                               [2,   0.0,  87.4],
        #                               [2,  22.5,  87.4],
        #                               [2,  45.0,  87.4],
        #                               [2,  67.5,  87.4],
        #                               [2,  90.0,  87.4],
        #                               [2, 112.5,  87.4],
        #                               [2, 135.0,  87.4],
        #                               [2, 157.5,  87.4],
        #                               [2, 180.0,  87.4],
        #                               [3,   0.0,  91.8],
        #                               [3,  22.5,  91.8],
        #                               [3,  45.0,  91.8],
        #                               [3,  67.5,  91.8],
        #                               [3,  90.0,  91.8],
        #                               [3, 112.5,  91.8],
        #                               [3, 135.0,  91.8],
        #                               [3, 157.5,  91.8],
        #                               [3, 180.0,  91.8]]
        #   },
        {
          "range": [[89.5, 91.7]],
          "background": [3,3],
          "peak": [{
              "phase": "MgO",
              "hkl": 222,
              "d-space": [0,2],
              "height": 3,
              "profile": 0,
              #"profile_fixed": 0.5,
              "width": 0,
              "symmetry": 2
            }],
        },
    ]

#Output settings
Output_type        = "PolydefixED"
Output_ElasticProperties = "../MgO_elastic_properties.txt"
