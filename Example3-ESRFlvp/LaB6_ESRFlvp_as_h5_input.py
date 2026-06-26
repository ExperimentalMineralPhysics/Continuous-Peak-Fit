# Example 3 -- LaB6_ESRFlvp_as_h5; load h5 data and combine multiple frames to single pattern
# 
# Takes top level h5 file and reads data, combining multiple frames into a single image.

# properties of the data files.
datafile_directory = "./data/"
datafile_Basename  = 'lab6_53keV_3900mm_12steps.h5'

h5_datakey = '/*.1/measurement/p900kw'
h5_iterate = [
        # 1st level -- at the level of the collections in the hdf5 file (the * in h5_datakey)
        {"do":"iterate", 
         "from": 0, 
         "to": -1, 
         "step": 1, 
         "using":"position",
         "label":['pos']},
            # bottom level -- what to do with the data in the collection
            {"do":"combine", 
                     "from": 0, 
                     "to": -1, 
                     "step": 1, 
                     }]

# Calibration and masking.
Calib_type     = "ESRFlvp"
Calib_param    = "./lab6_53keV_3900mm.json"

# Fitting properties for peaks.
reduce_by = 1

#Output settings
Output_directory   = 'results_as_h5'
Output_type        = ['FitMovie', 'CoefficientTable']

# define ranges and peaks
fit_orders = [
    {
        "range": [3.2,3.25],
        "background": [2, 0],
        "peak": [{
            "phase": "LaB6",
            "hkl": '100',
            "d-space": 3,
            "height": 1,
            "profile": 0,
            "width": 0,
            "symmetry": 2
        }, ]
    },
    {
        "range": [4.53,4.6],
        "background": [2, 0],
        "peak": [{
            "phase": "LaB6",
            "hkl": '110',
            "d-space": 3,
            "height": 1,
            "profile": 0,
            "width": 0,
            "symmetry": 2
        }, ]
    },
    {
        "range": [5.55,5.65],
        "background": [2, 0],
        "peak": [{
            "phase": "LaB6",
            "hkl": '111',
            "d-space": 3,
            "height": 1,
            "profile": 0,
            "width": 0,
            "symmetry": 2
        }, ]
    },
    {
        "range": [6.4,6.5],
        "background": [2, 0],
        "peak": [{
            "phase": "LaB6",
            "hkl": '200',
            "d-space": 3,
            "height": 1,
            "profile": 0,
            "width": 0,
            "symmetry": 2
        }, ]
    },
]
