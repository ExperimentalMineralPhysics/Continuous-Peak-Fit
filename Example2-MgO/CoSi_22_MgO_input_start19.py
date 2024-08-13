# Input parameters for CoSi_22. Example for Continuous peak fit.

# properties of the data files.
datafile_directory = "."
datafile_Basename = "COSI_022_"
datafile_Ending = ".med"
datafile_Files = [  # 6, 9, 11, 13, 15, 17,
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
    67,  # 69, #71, 73
]
"""
datafile_StartNum  = 19
datafile_EndNum    = 67
datafile_Step      = 2
"""
datafile_NumDigit = 4

# Calibration and masking.
Calib_type = "Med"
Calib_detector = "10-element"
# Calib_detector = 'W2010_10element_moved.txt'
Calib_data = "COSI_021_0001.med"
Calib_param = "COSI_021_0001.med"
Calib_mask = [5, 10]

# Fitting properties for peaks.
fit_track = True
fit_propagate = True
fit_bounds = {
    "background": ["0.9*min", "max"],
    "d-space": ["min", "max"],
    "height": [0, "1.2*max"],
    "profile": [0, 1],
    "width": ["range/(ndata)", "range/4"],
}
fit_orders = [
    {
        "range": [[43.5, 46.0]],
        "background": [3],
        "background-type": "independent",
        "peak": [
            {
                "phase": "MgO",
                "hkl": 111,
                "d-space": [0, 2],
                "height": 3,
                "height-type": "independent",
                "profile": 0,
                # "profile_fixed": 1,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
    {
        "range": [[51.8, 53.3]],
        "background": [3],
        "background-type": "independent",
        "peak": [
            {
                "phase": "MgO",
                "hkl": 200,
                "d-space": [0, 2],
                "height": 3,
                "height-type": "independent",
                "profile": 0,
                # "profile_fixed": 1,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
    {
        "range": [[73.0, 75.5]],
        "background": [3],
        "background-type": "independent",
        "peak": [
            {
                "phase": "MgO",
                "hkl": 220,
                "d-space": [0, 2],
                "height": 3,
                "height-type": "independent",
                "profile": 0,
                # "profile_fixed": 1,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
    # MgO(311) is overlapped and ignored.
    {
        "range": [[89.5, 91.7]],
        "background": [3, 3],
        "background-type": "independent",
        "peak": [
            {
                "phase": "MgO",
                "hkl": 222,
                "d-space": [0, 2],
                "height": 3,
                "height-type": "independent",
                "profile": 0,
                # "profile_fixed": 0.5,
                "width": 0,
                "symmetry": 2,
            }
        ],
    },
]

# Output settings
Output_directory = "./results/"
Output_type = "PolydefixED"
Output_ElasticProperties = "MgO_elastic_properties.txt"
