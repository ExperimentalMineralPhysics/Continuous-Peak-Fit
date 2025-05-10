.. Continuous-Peak-Fit documentation master file, created by
   sphinx-quickstart on Fri Aug 28 15:27:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Subpatterns and peaks
=====================================

Each subpattern is defined by number of parameters contained within a python dictionary. The dictionary for each subpattern is contained with in a list called 'fit_orders'. e.g.:

 .. code-block:: python

  fit_orders = [
         {
           "range": [[3.5,3.8]],
           "background": [1,0],
           "peak": [{
               "phase": "Steel",
               "hkl": '111',
               "d-space": 2,
               "height": 12,
               "profile": 0,
               #"profile_fixed": 1,
               "width": 0,
               #"symmetry": 2
             }],
           },
         {
           "range": [[4.1,4.4]],
           "background": [1,0],
           "peak": [{
               "phase": "Steel",
               "hkl": '200',
               "d-space": 2,
               "height": 12,
               "profile": 0,
               #"profile_fixed": 1,
               "width": 0,
               #"symmetry": 2
             }],
         },
         ]


With in each subpattern dictionary the following parameters are possible or needed: 

==========================  ==========================   ================================
Parameter                   Required?                    Description
==========================  ==========================   ================================
``range``                   yes                          Defines the maximum and minimum of the range to fit with peak in. For angle dispersive diffraction data this is the maximum and minimum two theta angle (in degrees), for energy disperdive diffraction it is the maximum and minimum energy (in keV). The format is ``[[min, max]]``. 
``background``              yes                          This is a list of numbers that defines the background, see LINK.
``peak``                    yes                          List of dictionaries, one for each peak. See LINK
``imax``                    no                           Crops the maximum intensity within the range at this value. Formatted as a single number. 
``imin``                    no                           Crops the minimum intensity within the range at this value. Formatted as a single number. 
``PeakPositionSelection``   only for multiple peaks      a list of [peak_n, azimuth, two theta], that makes a guess for each peak in the range.
==========================  ==========================   ================================



Peak definition
------------------------------

The peaks are a pesugo-Voigt peak with parameters:






Fixing parameters to each other
------------------------------

Parameters can be linked to previously defined paramerters, within the same range. 
This is done by labelling the parameter to a previously labelled parameter.


 .. code-block:: python

         ...
         {
           "range": [[4.1,4.4]],
           "background": [1,0],
           "peak": [{
               "phase": "Steel",
               "hkl": '111',
               "d-space": 2,
               "height": 12,
               "profile": 0,
               #"profile_fixed": 1,
               "width": 0,
               #"symmetry": 2
             },{
               "phase": "Steel",
               "hkl": '200',
               "d-space": [2, **"order2equal0"**],
               "height": 12,
               "profile": **"all"**,
               #"profile_fixed": 1,
               "width": **"less0"**,
               #"symmetry": 2
             }],
          ...



syntax:

The syntax allows parameters to be fixed to the same parameter in a preceeding peak within the same range. 

**is this 0 counted or 1 counted??**
Either way we need to check for min/max out of range


**Syntax definition**

String that goes within list. The string is formed of two parts: 

A. part saying which parts of the series to fix. Commands are:
    
    * all
    
    * order (ideally for Fourier series)

    * n-m, where n and m are integer values (ideally for splines)

B. a part defining its relation to the other parts. It is formed of a string, optionally followed by a number. Possible strings are:

    * equal

Other relations (greater than, less than, times, divide) are not implemented yet.

If either part is missing then it is assumed to be "all" and "equal0". Hence an empty string will make all the series values equal to those of the first peak.



**for fourier series**

"d-space": [2, "2ndOrderEqualPeak0"] --- make 2nd order terms the same betwen this peak and peak 0.

"d-space": [2, "1st2ndOrderEqualPeak0"] --- make 1st and 2nd order terms the same betwen this peak and peak 0.

"d-space": [0,2, "2ndOrderEqualPeak0"] --- no 1st order terms, make 2nd order terms the same betwen this peak and peak 0 (compatible with [0,2]).


"profile": [0, "AllEqual0"] --- make all terms the same as first peak. 

"profile": "AllEqual0" --- make all terms the same as first peak. Same as [*n*, "allwith1"].
"profile": "Fixed0.04" --- fixed as single value of 0.04.
"profile": [1, "Fixed(0.04,0.02,0)", "AllEqual0"]


*strings can be combined. *
"d-space": [2, "1-2OrderEqualPeak0", "0thOrderLessPeak0"] --- make 1st and 2nd order terms the same betwen this peak and peak 0. make 0th order less than peak 0. 

**for splines**

"height": [12, "AllEqualPeak0"] --- make all terms the same as first peak.

"height": [12, "4-13EqualPeak0"] --- make values 4-13 (of 25) the same as first peak.

"height": [12, "4-13GreaterPeak1"] --- make values 4-13 (of 25) greater than those of the second peak.



To fix 3 or more peaks they just have to be assembled in the correct order. 
