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
``background``              yes                          This is a list of numbers, see LINK. The length of the list is the order 
``Imax``                    no                           Crops the maximum intensity within the range at this value. Formatted as a single number. 
``Imin``                    no                           Crops the minimum intensity within the range at this value. Formatted as a single number. 
``PeakPositionSelection``   only for multiple peaks
``peak``                    yes                          List of dictionaries. See 
==========================  ==========================   ================================


Peak definition
------------------------------