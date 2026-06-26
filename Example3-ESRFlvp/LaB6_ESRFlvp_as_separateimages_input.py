# Example 3 -- LaB6_ESRFlvp_as_separateimages; load separate images and combine multiple frames to single pattern
# 
# Uses wild card to find the images, combining multiple frames into a single image.
# in contrast with 'LaB6_ESRFlvp_as_h5' which uses the top level h5 file to read the images. 

# properties of the data files.
datafile_directory = "./data/scan0001/"
datafile_Basename  = 'p900kw_*.h5'

h5_datakey = '/entry_0000/measurement/data'

# Calibration and masking.
Calib_type     = "ESRFlvp"
Calib_param    = "./lab6_53keV_3900mm.json"

# Fitting properties for peaks.
reduce_by = 1

#Output settings
Output_directory   = 'results_as_separateimages'
Output_type        = ['Polydefix', 'FitMovie', 'CoefficientTable']

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
