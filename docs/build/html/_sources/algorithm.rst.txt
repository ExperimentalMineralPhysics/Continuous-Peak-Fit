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





Meaning of 'Fit Status' 
---------------------------
The status in the output json file is a record of all the steps that the fitting passed through. 
The assignent of each value and that path to get the value is as follows:


*Initial 'step' value*:

   step = 0 [new fit, do chunks]

   step = 5 [fit using previous parameters]


*Initiate parameters and perform inital fit to chunks*. At the end of which:

   step = step + 10

Giving:

   step = 10 [new series fits, using chunked fits]

   step = 15 [series fits from previous parameters]

*Refine the peak parameters one at a time.*

if refinementments are used:

   step = step + 10

otherwise:

   step = step + 11

Step values at the end of the first attempt at refinements:

   step = 20 [new series fits, refined]

   step = 21 [new series fits, not refined]

   step = 25 [series fits from previous parameters, refined]

   step = 26 [series fits from previous parameters, not refined]

Other values (22 or 27) are possible if the peak parameters have been refined more than once. 


*Perform a global fit to all parameters*.


if the step value is 20, 21, 22, 25, 26, or 27

   use default maximum function evaluaitons (max_f_eval)

otherwise if the step value is 23 or 28

   use 2 time max_f_eval

otherwise for step values of 24 or 29

   use inf for  max_f_eval


*Determine the state of the fit and act accordingly*

if the fit is sucsessful (i.e. lmfit sucsess==1) and there are no problems with the fit

   step = step+100

alternatively if the fit is sucsessful but there are problmes with the fit 

   either 

      try changing some of the setttings and trying again

   or 

      Dicscard everything and start again


The possible final status values are:

* **-176** -- the fitting ended up in a state that is not going to work. Exit neatly.
* **-11** -- no valid data, failed, exit neatly.
* **120** -- no previous fits, refined, finished
* **121** -- no previous fits, not refined, finished
* **122** -- no previous fits, not refined, failed, refined, finished
* **123** -- (no previous fits, refined, failed, refined, finished) or (no previous fits, not refined, failed, refined, failed, refined, finished)
* **124** -- not possible, raises an error
* **125** -- previous fits, refiend, finished
* **126** -- previous fits, not refined, finished
* **127** -- previous fits, not refined, failed, refined, finished
* **128** -- (previous fits, refined, failed, refined, finished) or (previous fits, not refined, failed, refined, failed, refined, finished)
* **129** -- not possible as final status. But it passes through everything and starts again: process is previous fits, (not refined), failed, refined, failed, refined, failed, discard all and start again with no previous fits. 
