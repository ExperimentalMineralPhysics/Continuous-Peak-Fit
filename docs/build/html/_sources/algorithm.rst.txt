.. Continuous-Peak-Fit documentation master file, created by
   sphinx-quickstart on Fri Aug 28 15:27:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


===============================================
Algorithm
===============================================


Many data sets contain 2D peaks that are generally treated in data processing as a series of indepentent, unrelated peaks. There are good reasons for this e.g. keeping errors independent. However, in many cases the data are explicitally related in which case the data can and should be treated as a single 2-dimensional data set. Treating thedata as a two D data set requires that each paramter is a function of the dispersion -- for X-ray diffraction, azimuth, time series dsata time. 

The continoous Peak Fit algorithm is contained within the script `XRD_FitPattern.py`. The potential complexity of 2D data requires a number of sequential steps each of which produces a more refined model of the data. The steps are:

* 1D peak fitting and contonuous function fitting 
* Refinement of dispersion parameters
* Final fit of all parameters. 


Without the incremental steps the minimiser is very unlikely to find the optimal solution. However, a solution for one data can be propagated onto the next if the two data are similar enough (e.g. XRD data).



1D peak fitting and contonuous function fitting
-----------------------------------------------


In a 1D diffraction pattern the intensity of a single peak (Ihkl) as a function of :math:`2\theta` can be approximated by a Gaussian, Lorentzian or a convolution of the two (usually a pesudoVoigt function). The peak has height H and width W. The intensity with :math:`2\theta` for a Gaussian peak is:

􏰀    `(2\theta − 2\theta_0)/2 − 2\sigma^2`

Ihkl(2θ) = H exp
where: H is the height of the peak, 2θ0 is the centroid and σ is the width of the distribution. For a Lorentzian peak:
where: γ is the HWHM of the distribution.

􏰂   `γ2 􏰃Ihkl(2θ) = H (2θ − 2θ0)2 + γ2 (β.2)`

the PseudoVoigt function is a linear combination of Gaussian and Lorentz peak shapes.
Ihkl(2θ) = H[ηL(2θ − 2θ0) + (1 − η)G(2θ − 2θ0)] (β.3)
where η is the fraction of Lorentz (L) and Gaussian (G) functions used.



Refinement of dispersion parameters
-----------------------------------

text goes here

Final fit of all parameters
---------------------------

text goes here.

