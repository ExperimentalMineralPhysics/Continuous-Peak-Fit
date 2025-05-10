# Example 1 -- BCC1_2GPa_100s_001; basic case 
# 
# This is a baseline example; the processes diffraction patterns for BCC iron.
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

#Output settings
Output_directory   = 'results'
Output_type        = ['Polydefix', 'DifferentialStrain', 'FitMovie', 'CoefficientTable', 'CollectionMovie']


# define ranges and peaks
fit_orders = [
    {
        "range": [11.47, 11.7],
        "background": [0, 0],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 110,
            "d-space": 2, 
            "height": 12, 
            "profile": 0, 
            "width": 0,
        }], 
    },
    {
        "range": [16.1, 16.6],
        "background": [0, 0],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 200,
            "d-space": 3,
            "height": 16,
            "profile": 0,
            "width": 0,
        }],
    },
    {
        "range": [19.8, 20.27],
        "background": [0, 0],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 211,
            "d-space": 2,
            "height": 16,
            "profile": 1,
            "width": 0,
        }],
    },
    {
        "range": [25.7, 26.4],
        "background": [2, 2],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 310,
            "d-space": 1,
            "height": 1,
            "profile": 0,
            "width": 0,
        }],
    }
]
