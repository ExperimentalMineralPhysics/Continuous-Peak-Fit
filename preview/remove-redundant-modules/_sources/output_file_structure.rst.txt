.. Continuous-Peak-Fit documentation master file, created by
   sphinx-quickstart on Fri Aug 28 15:27:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



=====================================
Output files:
=====================================

The code produces various json formatted files as default as well as optional image files. The Image files will be not discussed here. The json files are formatted as follows: 

Full fit *.json file
=====================================

Produced by XRD_FitPattern.execute('input_file.py')

Each diffraction pattern when fit produces a *.json file that contains the outputs from the fits. 
These are structured in a similar manner to the 'fit_orders' list of dictionaries in the input file but with additional values. 

There is a dictionary entry for each range in the input file. 
The dictionary contains the following enties. 

* DataProperties
* FitProperties
* PeakLabel
* background
* background_err
* background_type
* correlation_coeffs
* peak
* range

The peak entry is also a list of dictionaries, one for each peak in the input range. 
The parameters in each entry are the fitted values for the named value. 


Each peak parameter has three associated entries in the dictionary:

* parameter -- the fitted coefficients
* parameter_err -- the errors on these coefficients
* parameter_type -- the series type of the coefficients


An example dictionary entry for one range in the json file is: 

 .. code-block:: python
[
  {
    "DataProperties": {
      "max": 4812.0,
      "min": 0.0
    },
    "FitProperties": {
      "ChiSq": 184214466.53235897,
      "RedChiSq": 6144.988542676595,
      "aic": 261587.6361675344,
      "bic": 261654.10405461484,
      "chunks-time": 0.007439851760864258,
      "degree-of-freedom": 29978,
      "function-evaluations": 37,
      "n-data": 29986,
      "n-variables": 8,
      "status": 125,
      "sum-residuals-squared": 184214466.53235897,
      "time-elapsed": 2.0933940410614014
    },
    "PeakLabel": "Ni_FCC_alloy (111)",
    "background": [
      [
        6.853738672982084
      ]
    ],
    "background_err": [
      [
        0.6023096161637345
      ]
    ],
    "background_type": "spline_cubic",
    "correlation_coeffs": "{\"bg_c0_f0\": {\"peak_0_h0\": -0.06549834714563857, \"peak_0_h1\": -0.1666117505895424, \"peak_0_h2\": -0.1852866709498333, \"peak_0_d0\": 0.001256451455142335, \"peak_0_d3\": -0.0005122549843816782, \"peak_0_d4\": -0.0009477254922480384, \"peak_0_w0\": -0.3775816489154083}, \"bg_c0_f_tp\": null, \"peak_0_s0\": null, \"peak_0_h0\": {\"bg_c0_f0\": -0.06549834714563857, \"peak_0_h1\": 0.13653862651020207, \"peak_0_h2\": 0.15482343627632464, \"peak_0_d0\": 0.0014358433026510644, \"peak_0_d3\": 0.00034255031239486183, \"peak_0_d4\": -9.03014679756105e-05, \"peak_0_w0\": -0.4524384746501567}, \"peak_0_h1\": {\"bg_c0_f0\": -0.16661175058954245, \"peak_0_h0\": 0.13653862651020207, \"peak_0_h2\": 0.04582881574689647, \"peak_0_d0\": -0.0044554415717901015, \"peak_0_d3\": 0.00022823349789936358, \"peak_0_d4\": 0.004749905993261898, \"peak_0_w0\": -0.26130121917185667}, \"peak_0_h2\": {\"bg_c0_f0\": -0.1852866709498333, \"peak_0_h0\": 0.15482343627632464, \"peak_0_h1\": 0.045828815746896456, \"peak_0_d0\": 0.008886378058758136, \"peak_0_d3\": -0.0023206655636863567, \"peak_0_d4\": -0.013849007354279543, \"peak_0_w0\": -0.15325046782899288}, \"peak_0_h_tp\": null, \"peak_0_d0\": {\"bg_c0_f0\": 0.0012564514551423366, \"peak_0_h0\": 0.0014358433026510644, \"peak_0_h1\": -0.004455441571790101, \"peak_0_h2\": 0.008886378058758136, \"peak_0_d3\": -0.022822968502537577, \"peak_0_d4\": -0.1159353507330648, \"peak_0_w0\": 0.009292577955583374}, \"peak_0_d1\": null, \"peak_0_d2\": null, \"peak_0_d3\": {\"bg_c0_f0\": -0.0005122549843816782, \"peak_0_h0\": 0.0003425503123948614, \"peak_0_h1\": 0.00022823349789936358, \"peak_0_h2\": -0.0023206655636863567, \"peak_0_d0\": -0.022822968502537577, \"peak_0_d4\": 0.017357152303744103, \"peak_0_w0\": 0.0019994248689470685}, \"peak_0_d4\": {\"bg_c0_f0\": -0.0009477254922480399, \"peak_0_h0\": -9.030146797560995e-05, \"peak_0_h1\": 0.004749905993261897, \"peak_0_h2\": -0.013849007354279538, \"peak_0_d0\": -0.1159353507330648, \"peak_0_d3\": 0.017357152303744103, \"peak_0_w0\": -0.006883271430247167}, \"peak_0_d_tp\": null, \"peak_0_w0\": {\"bg_c0_f0\": -0.37758164891540824, \"peak_0_h0\": -0.4524384746501567, \"peak_0_h1\": -0.26130121917185667, \"peak_0_h2\": -0.15325046782899288, \"peak_0_d0\": 0.009292577955583376, \"peak_0_d3\": 0.0019994248689470685, \"peak_0_d4\": -0.006883271430247167}, \"peak_0_w_tp\": null, \"peak_0_p0\": null, \"peak_0_p_tp\": null}",
    "peak": [
      {
        "d-space": [
          2.0868950217728846,
          0.0,
          0.0,
          -0.0019622103654890477,
          0.0004326382134808848
        ],
        "d-space_err": [
          0.00043161061877024123,
          0,
          0,
          0.0005985061109717787,
          0.0006184947151738631
        ],
        "d-space_type": "fourier",
        "height": [
          49.50312863268371,
          29.817110795113273,
          22.617325696699567
        ],
        "height_err": [
          3.416258208867711,
          3.0443686696480197,
          3.3142337294737043
        ],
        "height_type": "spline_cubic",
        "hkl": 111,
        "phase": "Ni_FCC_alloy",
        "profile": [
          0.7
        ],
        "profile_err": [
          0
        ],
        "profile_type": "fourier",
        "symmetry": 1,
        "width": [
          0.018394009324434288
        ],
        "width_err": [
          0.0014175517606364613
        ],
        "width_type": "fourier"
      }
    ],
    "range": [
      [
        5.000004502946348,
        5.2999923226886
      ],
      [
        2.0201290460258505,
        2.1412477703125647
      ]
    ]
  }, 






Output files: *__chunks.json
=====================================

Produced by cpf.Cascade.execute?("inputfile.py")

The output file contains two top level lists. 
The first list is a list of peak parameters, and the second is the position of the chunks in the first list.

The first list is a list of dictionaries. one dictionary for each range in the input file. 
In has entries:

* peak -- peak name
* chunks -- mean azimuthal position of the chunk
* bg -- background fits
* bg_err -- errors thereon
* d -- dspacing fits
* d_err -- errors thereon
* h -- height fits
* h_err -- errors thereon
* p -- profile fits
* p_err -- errors thereon
* w -- width fits
* w_err -- errors thereon

Each dictionary entry is a list of lists. For the peak parameters (d,h,p,w), the number of the inner lists is the number of peaks in the range. The length of each of these lists is the number of azimuthal chunks that were fit (see ?? for how the azimuthal chunks were defined). 

The background (bg) is also a list. the length of this list of the length of the background list in the input file. It is a polynomial.

chunks lists the mean azimuth of each chunk.

peak labels the peak using the information provided in the input file. 

An example chunk output with a single azimuth (i.e. -180 to 180 degrees) is this:

 .. code-block:: python

[
    {
      "bg": [
        [
          2.691982806955949
        ]
      ],
      "bg_err": [
        [
          0.02297890339346997
        ]
      ],
      "chunks": [
        4.147222281400086
      ],
      "d": [
        [
          3.563848497645126
        ]
      ],
      "d_err": [
        [
          0.00203163076034729
        ]
      ],
      "h": [
        [
          1.2479759854565744
        ]
      ],
      "h_err": [
        [
          0.0630295184727987
        ]
      ],
      "p": [
        [
          0.7
        ]
      ],
      "p_err": [
        [
          0
        ]
      ],
      "peak": "Ni_FCC_alloy_superlattice (100)",
      "w": [
        [
          0.03531636856078305
        ]
      ],
      "w_err": [
        [
          0.0024119947016804956
        ]
      ]
    },


