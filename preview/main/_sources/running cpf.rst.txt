.. Continuous-Peak-Fit documentation master file

*Continuous peak fit*, via the command line
==================================================

.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit
.. _Example1: file:///../Fourier-Peak-Fit/Example1-Fe/
.. _Dioptas: http://www.clemensprescher.com/programs/dioptas
.. _Dioptas documentation : https://dioptas.readthedocs.io/en/stable/
.. _Dioptas integration : https://dioptas.readthedocs.io/en/stable/integration.html

.. _installation: installation.html

.. role:: python(code)
   :language: python


`Continuous-Peak-Fit` is currently run from the python command line. We will make a GUI and OS-specific executables but for now this does not exist. 

This tutorial uses the angle dispersive X-ray diffraction data contained in cpf's Example 1. The data is for BCC Iron from a DAC experiment and is the first 10 images of a data set of 600 consecutive images. 

`Continuous-Peak-Fit` is excuted by calling an input file that contains all the information about the data, calibration and what fitting is desired. Examples of these input files are contained in the `Example 1 directory <Example1_>`_. 

The full package consists of a nuber of subfunctions :python:`cpf.XRD_Fitpattern`. Although each stage is not required it is recommended to ease use of the package. 


=====================================
Calibration
=====================================

The detector calibrtation is undertaken in external software. The calibration and mask are saved for import to `cpf`. `Dioptas`_ is used for the calibrations and masks of angle dispersive diffraction data. Other software (e.g. GSAS-II or Fit2D) are implementable if desired. 

From Dioptas `Continuous peak fit` reuqires a calibration file and an optional image mask. See the `Dioptas documentation`_ for how to mask and calibrate and the detector. 



=====================================
Loading `continuous peak fit`.
=====================================
See :doc:`installation`. Once installed, run:

 .. code-block:: python    

   import cpf



=====================================
Input file creation
=====================================

The input file is structured as outlined here: :doc:`input file` and here: :doc:`input file short`. 
If is a python file containing a number of dictionaries and lists.

In the following documentation this file is refered to as `input_file`. The file name may have a py ending.



=====================================
Initiate
=====================================
Initiating the input checks that the input file is correctly formatted and if something is missing returns an error message. Run:

 .. code-block:: python    

   cpf.XRD_Fitpattern.initiate('input_file')

If this completes without returning any errors or other messages the input file is *probably* formatted correctly. It does not though check that the selections are sensible. 


=====================================
Set range
=====================================

The two-theta vs. azimuth data plotted in the `Dioptas integration`_ window is resampled onto an orthonormal two-theta--azimuth grid and consequently somewhat smoothed. `Continuous peak fit` uses the raw data from the detector, therefore is it worth checking that the data range is acceptable and the data does not contain too many dark or bright comsic spots that would distort the fitting. 

To check the data ranges, run:

 .. code-block:: python    

   cpf.XRD_Fitpattern.setrange('input_file')

This will make range plots for all the selected peaks. A typical figure will look something like this:

RANGE FIGURE. 

In the left hand column is the raw, unmasked data, plotted as two theta vs. azimuth and two theta vs. intensity. The middle column plots the mask (top) and the Intineisties cumulative distribution function (bottom). The right hand column contains the masked data again plotted as two theta vs. azimuth and two theta vs. intensity. 

The extent of the selection is adjusted by changing the max and min ``range`` values in the input file (LINK). It is sometimes necessary to limit the maximum intenstiy in each peak selection, which is done via the ``Imax`` peak definition (LINK).

Cheching the range for a single peak is done by:

 .. code-block:: python    

   cpf.XRD_Fitpattern.setrange('input_file', peak_selection = .....)

The peak_selection is zero counted, so the first peak is 0. 

`input_file`\_mask.tif files are saved in the output directory.


.. _initial_guess:

==========================================================================
Initial guess selection for multiple, weak or obscured peaks
==========================================================================

The strong peaks coninuous peak fit is able to guess the position of a single strong peak (LINK). But for multiple, weak or spotty peaks some inital guidance is needed.

Initial guesses can be derived from Dioptas manually by using the two theta and azimuth position of the curser on the `Dioptas integration`_. Alternatively, the peak guesses can be selected using the mouse +/- keyboard by running: 

 .. code-block:: python    

   cpf.XRD_Fitpattern.initialpeakposition('input_file', [peak_selection = n])

where ``n`` is the subpattern range to be selected and is an optional call. As with set range [LINK] the subpatterns are 0 counted. Leaving the ``peak selection`` out works through the initial guesses for all the peak ranges. 

This presents a file of the masked data like so:

IMAGE

The position of the peaks are selected using the mouse or keyboard. peaks should be selected 1 - n. The left mouse botton returns '1' and the right mouse button is mapped to '2'. More peaks have to be selected using the keyboard. When the figure is closed the peak selection is printed in the command window.

It is formatted as:


 .. code-block:: python    

   [
   [ peak num, azimuth, two theta],
   ...
   ]


The peaks are 1 counted. This array needs to be copied into the input file manually. The array is the content of the "InitialPeakPosition" ??? of the peak list (LINK). 

An example of a properly structured input is: 

COPY SOMETHING FROM AN INPUT FILE



=====================================
Order search
=====================================


The optimal number of coefficients (the order) for each parameter is not necessarily obvious. For smooth but structured diffraction data a good guess for the order of the peak height is the number of peaks in the diffraction ring. But this is not necessarily optimal, and there is often a ply off between quality of fit and execution time.

The parameter that is most strongly affected by the order is the peak height because for X-ray data the other coefficients have orders generally <= 2. To run a sweep over the possible height orders run:

 .. code-block:: python    

   cpf.XRD_Fitpattern.ordersearch('input_file')

This runs a seach of heights over all the subpatterns and peaks. 

To change the parameter or peak searched over add the optional calls, formatted as option=set:

===================  ==========================   ================================
option               Values                       default, notes
===================  ==========================   ================================
search_parameter     string of peak parameter     default 'height' (e.g. 'background', 'height', 'd-space'... )
search_over          [min,max]                    [0,20]
subpattern           'all' or list of numbers     'all'     no        
search_peak          number                       0, zero counted, can't be greater than the number of peaks, -1 for last peak is permitted. 
search_series        series string                ['fourier', 'spline'] 
===================  ==========================   ================================




=====================================
Execution
=====================================

Running the fitting is done by calling: 

 .. code-block:: python    

   cpf.XRD_Fitpattern.execute('input_file')

This fits all the files in turn, propagating the fit from one file to the next (LINK). All possible swithes are listed here: LINK.
The most likely ones (and their defaulr settings) are:

 .. code-block:: python    

   cpf.XRD_Fitpattern.execute('input_file', 
                                  debug = False, 
                                  propagate = True, 
                                  parallel = True, 
                                  save_all=False)

debug makes some more figures.

save_all makes even more files than debug, saving an image for each subpattern fitted.

propagate uses the fit from the previous data file as a initial guess for the current one. 

parallel uses the parallel options (if installed), speeding up the code exection.


A json file is created for each diffraction pattern which contains the fit parameters.
For the first file in the sequence a figure of the fit for each region is also saved, e.g.:

FIT FIGURE. 

We also save a \*.sav file which contains all the lmfit fit object (LINK). 


=====================================
Creating output files. 
=====================================
The execute call will generate the output files listed in  ``input_file``. But the outputs can be generated directly, from the json files, by calling:

 .. code-block:: python    

   cpf.XRD_Fitpattern.write_output('input_file')

The optional call ``write_output = ...`` allows for the creation of different output types without editing the input file. 

The currently recognised output types are listed here: LINK TO OUTPUT TYPES. 

Additional types can be generated by editing the cpf.Write... file to make a file formated as you wish. 






