# Input parameters for Example 1 -- BCC1_2GPa_100s_001. 

# properties of the data files.
datafile_directory = './'
datafile_Basename  = 'BCC1_2GPa_10s_001_'
datafile_Ending    = '.tif'
datafile_StartNum  = 1
datafile_EndNum    = 10
datafile_NumDigit  = 5

# Calibration and masking.
Calib_type     = "Dioptas"
Calib_detector = 'Pilatus1M'
Calib_data     = datafile_directory + 'CeO2_Pil207_E30_2Nov2016_001.tif'
Calib_param    = datafile_directory + 'CeO2_cal_Dioptas.poni'
Calib_mask     = datafile_directory + 'DiffractionMask_Dioptas.mask'
Calib_pixels   = 172

# Fitting properties for peaks.
fit_bounds = {
              "background": ['min', 'max'],
              "d-space":    ['min', 'max'],
              "height":     [ 0,    '2*max'],
              "profile":    [ 0,     1],
              "width":      [ 'range/(ndata)',  'range/2'],
              }
fit_orders = [
        {
          "range": [[10.8, 11.67]],
          "background": [2, 0],
          #"background-type": "spline-cubic",
          "peak": [{
              "phase": "Other",
              "hkl": '000',
              "d-space": 3,
              #"d-space-type": "spline-cubic",
              "height": 1,
              #"height-type": "fourier",
              "profile": 0,
              "profile_fixed": 0.1,
              "width": 0,
              # "width_fixed": 0.035,
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
        {
          "range": [[16.1, 16.6]],
          "background": [0, 0],
          "peak": [{
              "phase": "Fe-BCC",
              "hkl": 200,
              "d-space": [0,2],
              "height": 8,
              "profile": 0,
              "profile_fixed": 1,
              "width": 0,
              "symmetry": 2
            }],
        },
    ]

#Number of bins for initial fitting.
AziBins = 90
cascade_number_bins = 100

#Output settings
Output_directory   = './results/'
Output_type        = ['Polydefix', 'DifferentialStrain', 'FitMovie']
Output_NumAziWrite = 90
phase              = 'Fe-BCC'
Output_ElasticProperties = 'Properties_Fe-BCC.txt'