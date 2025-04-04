.. Continuous-Peak-Fit documentation master file


=====================================
XRD_FitPattern's examples and the features they test. 
=====================================

.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit


Each input file in the examples is made to demonstrate (and test) a separate feature of the code. 
The test files and their demonstrated feature are as follows:


Example 1: Monochromatic X-ray diffraction
=====================================

The sample in these examples is BCC iron at c. 2 GPa. 

1. **'BCC1_Dioptas_input.py' -- basic test**
    This function tests the basic fitting and outputs with no constraints on the fit. The series are all the default Fourier series. 
    The example creates all the default output types: 'Polydefix', 'DifferentialStrain', 'FitMovie', 'CoefficientTable', 'CollectionMovie'. All with default settings. 

1b. **'BCC1_Dioptas_EmptyRanges_input.py' -- basic test**
    This is an empty case; see that everything runs when there are no peaks in the data.

2. **'BCC1_Dioptas_SymmFixed_input.py' -- symmetry, fixed profile and intensity clipping test**
     This input file has fixed peak profiles and applies a symmetry to the diffraction data. 
     The intensity in each range is clipped to being between "imin" and "imax".
     Changes the fit_bounds applied to the fitting
     Change the number of bins used to inital fitting of the Fourier series
     The example creates output types: 'Polydefix' (with custom settings), 'CoefficientTable', 'CollectionMovie'.
    
     N.B. The restrictions to imin and imax are not (necessarily) useful restrictions to this data. It is for testing purposes
     

3. **'BCC1_Dioptas_MultiPeak_input.py' -- multipeak test**
    This example fits two peaks in the same range. 
    It also applies the PeakPositionSelection to a range with a single peak.

4. **'BCC1_Dioptas_SeriesFunctions_input.py' -- multipeak test**
    Demonstrates all the different series types, except for "independent" (see Example 2).
    Also demonstrates restrictions on the coefficients active in the series. 
    Some series have more than 1 acceptable name, see series documentation for full details.

5. **'BCC1_Dioptas_ConstrainedParams_input.py' -- parameter constraint test**
    This example demonstrates foring the coeficcients in two peaks in the same range to have the same values. It is only applicable to ranges with more than 1 peak.

6. **'BCC1_Dioptas_PeakShape_input.py' -- test of different peak shape test**
    Not implmented

7. **'BCC1_Dioptas_ImPrepare_input.py' -- test of different image pre-filtering functions**
    Not completed. 


Example 2: Angle dispersive X-ray diffraction
=====================================

The sample in these examples is MgO at c. 4 GPa. 
The data was collected at the X17B2 beamline at the NSLS, using the 10-element detector of Weidner et al (2010). 

1. **'BCC1_MED_input.py' -- basic test**
    This function tests the basic fitting and outputs with no constraints on the fit. The series are all the default Fourier series. 




CoSi_22_MgO_TfromTC_input.py
CoSi_22_MgO_input_start17down.py
CoSi_22_MgO_input_start19.py
CoSi_22_MgO_input_start19_fourier.py