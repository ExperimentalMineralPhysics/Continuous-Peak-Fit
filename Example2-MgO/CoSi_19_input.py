# Input parameters for CoSi_19

# data drive.
import os

if os.path.isdir("/Volumes/My Passport"):  # portable drive
    drive = "/Volumes/My Passport/"
elif os.path.isdir("/media/simon/My Passport"):  # analysis on Linux box
    drive = "/media/simon/My Passport/"
elif os.path.isdir("/volumes/external"):  # desk data store
    drive = "/Volumes/External/"
else:
    drive = ""

# properties of the data files.
datafile_directory = "/Users/simon/Desktop/CoSi_data/COSI_19"
datafile_Basename = "COSI_019_"
datafile_Ending = ".med"
datafile_Files = [
    2,
    5,
    7,
    9,
    11,
    13,
    15,
    17,
    19,
    21,
    23,
    25,
    27,
    29,
    31,
    33,
    35,
    37,
    39,
    41,
    43,
    45,
    47,
    49,
    51,
    53,
    55,
    57,
    59,
    61,
    63,
    65,
    67,
    69,
    71,
    73,
    75,
    77,
    79,
]
datafile_NumDigit = 4

# Calibration and masking.
Calib_type = "Med"
Calib_detector = "10-element"
Calib_data = "/Users/simon/Desktop/CoSi_data/COSI_16/COSI_016_0053.med"
Calib_param = "/Users/simon/Desktop/CoSi_data/COSI_16/COSI_016_0053.med"
Calib_mask = [10]  #'../DiffractionMask_49GPa.mask'
Calib_pixels = 72

# Fitting properties for peaks.
fit_orders = [
    {
        "range": [[44.0, 46.5]],
        "background": [3],
        "peak": [
            {
                "phase": "MgO",
                "hkl": 111,
                "d-space": [0, 2],
                "height": 4,
                "profile": 0,
                # "profile_fixed": 1,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
    {
        "range": [[51.0, 53.2]],
        "background": [3],
        "peak": [
            {
                "phase": "MgO",
                "hkl": 200,
                "d-space": [0, 2],
                "height": 4,
                "profile": 0,
                # "profile_fixed": 1,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
    {
        "range": [[72.0, 75.0]],
        "background": [3],
        "peak": [
            {
                "phase": "MgO",
                "hkl": 220,
                "d-space": [0, 2],
                "height": 4,
                "profile": 0,
                # "profile_fixed": 1,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
    {
        "range": [[83.4, 89.0]],
        "background": [3],
        "peak": [
            {
                "phase": "Other",
                "d-space": [0, 2],
                "height": 4,
                "profile": 0,
                "profile_fixed": 0.5,
                "width": 1,
                "symmetry": 2,
            },
            {
                "phase": "MgO",
                "hkl": 311,
                "d-space": [0, 2],
                "height": 4,
                "profile": 0,
                "profile_fixed": 0.9,
                "width": 1,
                "symmetry": 2,
            },
        ],
        "PeakPositionSelection": [
            [1, 0.0, 84.7],
            [1, 22.5, 84.7],
            [1, 45.0, 84.7],
            [1, 67.5, 84.7],
            [1, 90.0, 84.7],
            [1, 112.5, 84.7],
            [1, 135.0, 84.7],
            [1, 157.5, 84.7],
            [1, 180.0, 84.7],
            [2, 0.0, 86.7],
            [2, 22.5, 86.7],
            [2, 45.0, 86.7],
            [2, 67.5, 86.7],
            [2, 90.0, 86.7],
            [2, 112.5, 86.7],
            [2, 135.0, 86.7],
            [2, 157.5, 86.7],
            [2, 180.0, 86.7],
        ],
    },
    {
        "range": [[88.5, 93.0]],
        "background": [3],
        "peak": [
            {
                "phase": "MgO",
                "hkl": 222,
                "d-space": [0, 2],
                "height": 4,
                "profile": 0,
                # "profile_fixed": 0.5,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
]

# Output settings
Output_type = "PolydefixED"
Output_ElasticProperties = "../MgO_elastic_properties.txt"
Output_TemperaturePower = [2.62844e01, 3.06400e00, 4.11279e-03]
