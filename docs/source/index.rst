.. Continuous-Peak-Fit documentation master file

2-dimensional peak fitting 
=====================================


.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit
.. _Fit2D:   https://www.esrf.fr/computing/scientific/FIT2D/
.. _pyFAI: https://pyfai.readthedocs.io/en/master/index.html

`Continuous-peak-fit` is a python package for fitting peaks in 2-dimensional (X, Y + intensity) data. It was originally developed to improve the precision of azimuthally dependent properties determined from X-ray diffraction data collected on area detectors. With suitable input and output files, it is applicable to time-series and any other 2D data (e.g. radio astronimy data).

MAKE FIGURE TO GO HERE

The core idea of `continuous peak fit` is that peaks in two dimensional data are often inherently continuous -- the parameters describing a peak (in X: position, height, width, shape) are related to those describing the peak in the immedately adjacent position (Y+dY, where Y is azimuth, time, etc.).
For the most part previous data fitting routines have ignored this relationship, and either collapsed all dispersed data in to a single peak loosing the dispersion information, or treated the data as a series of independent peaks.

`Continuous peak fit` package was develped to overcome these limitations by treating peaks as a single two-dimensional objects. The properties of the peak are numerated as fourier series or spline functions in the Y. Treating the data as continuous improves the resolution of the fitted parameters as well as enabling fitting of previously un-fittable data.


==========
X-ray data
==========

Historically, 2D X-ray diffraction data has been investigated for azimuthal properties (e.g. differential stress and fabric) by transforming the data into a number of 1D diffraction patterns using software such as `Fit2D`_ and `pyFAI`_. Each of these 1D diffraction patterns represents a small range of azimuth (χ). The peaks in each slice are fit independantly and from these fits azimuthal properties (e.g. differential strain, stress and fabric) are calcualted. 
This approach has been used in the analysis of the azimuthally changing properties in diffraction patterns but it has its limitations. The drawback is that resampling of the 2D diffraction pattern reduces the number of data and necessarily applies a degree of averaging. The resampaling also makes the data points equally dense over the entire 2θ range, something that is not true of the original area detector data.

Rather than resampling the data `Continuous peak fit` uses the original detector observations and describes the diffraction peaks as continuous finctions in azimuth (χ). This description enables the fitting of weak, spotty and incomplete peaks or samples with extremally strong LPO that was not possible previously. Moreover, this method better accounts for the exclusions required by certain detectors and diffraction geometries.


============
Radio data
============

Noise in radio data looks a bit like the signal in X-ray diffraction patterns.

More Blurb. 


=========================
Software
=========================


The `Continuous-peak-fit` package is Free software, using the MIT open source license. Both the software and this documentation are works in progress. Please provide any feedback (good or bad) to the developers via GitHub (`cpf GitHub repository`_). 



===========
Citation
===========
If you use continuous peak fit for your research, please reference the manuscript below:

.. please reference one (or more) of the following papers that best fits your application:

.. literalinclude:: cpf.bib
   :language: latex



===============
Acknowledgments
===============
This package was developed while the authors were funded by NERC and SKA. 


===================
Documentation Parts
===================

.. toctree::
   :maxdepth: 2

   algorithm
   installation
   running cpf
   input file
   input file short
   whatsnew



==================
Indices and tables
==================
* :ref:‘index‘
* :ref:‘modindex‘
* :ref:‘search‘




