# Example 3 -- LaB6_Dioptas_as_sequence; load images as series. 
# 
# Test that Dioptas code can load hdf5 files. 
# It is ingoring the fact that the data files are all from the same collection of the detector. 

# properties of the data files.
datafile_directory = "./data/scan0001/"
datafile_Basename  = 'p900kw_'
datafile_Ending    = '.h5'
datafile_StartNum  = 0
datafile_EndNum    = 330
datafile_NumDigit  = 4
datafile_Step      = 30

# the h5 datakey is standard in Fabio so not required (in example data)
# h5_datakey = '/entry_0000/measurement/data'

# Calibration and masking.
Calib_type     = "Dioptas"
Calib_param    = "./poni_calib_53keV_3900mm_0.poni"

#Output settings
Output_directory   = 'results_as_sequence'
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
