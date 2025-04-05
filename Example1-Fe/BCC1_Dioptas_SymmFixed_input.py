# Example 1 -- BCC1_2GPa_100s_001; fixed profile, symmetry and intensity limits. 
# 
# This input file has fixed peak profiles and applies a symmetry to the diffraction data.
# The intensity in each range is clipped to being between "imin" and "imax".
# Changes the fit_bounds applied to the fitting
# Change the number of bins used to inital fitting of the Fourier series
# N.B. The restrictions to imin and imax are not (necessarily) useful restrictions to this data. It is for testing purposes

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

# Fitting properties for peaks.
AziBins = 45 # Number of bins for initial fitting.
fit_bounds = {
              "background": ['min', 'max'],
              "d-space":    ['min', 'max'],
              "height":     [ 0,    '1.2*max'],
              "profile":    [ -.1,     1.1],
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
        "range": [11.47, 11.7],
        "imax":500,
        "imin":5,
        "background": [0, 0],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 110,
            "d-space": 2, 
            "height": 12, 
            "profile": 0, 
            "profile_fixed": [0],
            "width": 0,
        }], 
    },
    {
        "range": [16.1, 16.6],
        "imax":27,
        "imin":5,
        "background": [0, 0],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 200,
            "d-space": 3,
            "height": 8,
            "profile": 0,
            "profile_fixed": [1],
            "width": 0,
            "symmetry": 2
        }],
    },
    {
        "range": [19.8, 20.27],
        "imax":27,
        "imin":5,
        "background": [0, 0],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 211,
            "d-space": 3,
            "height": 8,
            "profile": 1,
            "profile_fixed": [0.5, 0.25, 0.25],
            "width": 0,
            "symmetry": 2
        }],
    },
    {
        "range": [25.7, 26.4],
        "imax":27,
        "imin":5,
        "background": [2, 2],
        "peak": [{
            "phase": "Fe-BCC",
            "hkl": 310,
            "d-space": 2,
            "height": 1,
            "profile": 0,
            "profile_fixed": [0.5],
            "width": 0,
            "symmetry": 2
        }],
    }
]