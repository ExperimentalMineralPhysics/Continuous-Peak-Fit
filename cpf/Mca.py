"""
This module defines a device-independent MultiChannel Analyzer (MCA) class,
and a number of support classes.

Author:
   Mark Rivers
   
Created:
   Sept. 16, 2002.  Based on my earlier IDL code.

Modifications:
   Sept. 24, 2002 MLR
      - Fixed bug in saving ROIs in Meds
      - Fixed bug reading environment variables
   Sept. 25, 2002
      - Changed McaPeak.get_calibration to McaPeak.update, make more fields
        consistent.
      - Added McaPeak.ignore field to fix problem with other fields getting
        clobbered in fitPeaks.
      - Fixed serious bug in .d_to_channel()
   SAH - June 2019.
     - removed background and peak fitting functions because not required in integrated code.
     - removed import of CARSmath Xrf and fitPeaks
"""

import numpy as Numeric
#import string
import copy
import math
import sys
import time
#import CARSMath
#import Xrf
#import fitPeaks

########################################################################
class McaBackground:
   """
   Defines the parameters for fitting backgrounds in Mca objects.
   Fields and default values:
      .exponent     = 2
      .top_width    = 0.
      .bottom_width = 4.
      .tangent      = 0
      .compress     = 4
      
   See the documentation on fit_background() for information on the meaning
   of these fields.
   """
   def __init__(self):
      self.exponent     = 2
      self.top_width    = 0.
      self.bottom_width = 4.
      self.tangent      = 0
      self.compress     = 4

########################################################################
class McaFit:
   """
   Class used for the global input and output parameters for fit_peaks.
   """
   def __init__(self, mca=None):
      if (mca == None): mca=Mca()
      calibration = mca.get_calibration()
      self.npeaks =       0      # Number of peaks to fit
      self.first_chan =   0      # First channel to fit
      self.nchans =       len(mca.data)   # Number of channels 
      self.last_chan =    self.nchans-1   # Last channel to fit
      self.nparams =      0      # Number of fit parameters
      self.initial_energy_offset = calibration.offset
      self.initial_energy_slope =  calibration.slope
      self.energy_offset =  0.   # Energy calibration offset
      self.energy_slope =   0.   # Energy calibration slope
      self.energy_flag =    1    # Energy flag
                            #   0 = Fix energy calibration coefficients
                            #   1 = Optimize energy calibration coefficients
      self.initial_fwhm_offset = .15 # FWHM offset = 150 eV
      self.initial_fwhm_slope =  0.
      self.fwhm_offset =    0.   # FWHM offset
      self.fwhm_slope =     0.   # FWHM slope
      self.fwhm_flag =      1    # Fwhm flag
                            #   0 = Fix FWHM coefficients
                            #   1 = Optimize FWHM coefficients
      self.chi_exp =      0.     # Exponent of chi
      self.max_eval =     0      # Maximum # function evaluations
      self.n_eval =       0      # Actual number of function evalutions
      self.max_iter =     20     # Maximum number of iterations
      self.n_iter =       0      # Actual number of iterations
      self.tolerance =    1.e-4  # Convergence tolerance
      self.chisqr =       0.     # Chi-squared on output
      self.status =       0      # Output status code
      self.err_string =   ''     # Output error string
      self.debug =        0      # 0 (no debug) or 1
      self.results_file = 'fit_results.txt'  # Name of file for results
      self.spreadsheet_file = 'fit_spreadsheet.txt'  # Name of spreadsheet
      self.background =   McaBackground()  # Background object
      self.peaks =        []     # List of McaPeak() objects

   def update(self, mca):
      """
      Updates the McaFit object to be consistent with the Mca object.
      The calibration offset and slope are copied to the initial values,
      and the number of channels is set.

      Inputs:
         mca:
            An Mca object to copy the calibration from
      """
      cal = mca.get_calibration()
      self.initial_energy_offset = cal.offset
      self.initial_energy_slope = cal.slope
      data = mca.get_data()
      self.nchans = len(data)
      self.last_chan = self.nchans-1

########################################################################
class McaPeak:
   """
   Class for definiing the input and output parameters for each peak in
   fit_peaks().
   Input fields set bfore calling fit_peaks(), defaults and descriptions
      .label =       ""      # Peak label
      .self.energy_flag = 0  # Flag for fitting energy
                             #   0 = Fix energy 
                             #   1 = Optimize energy
      .fwhm_flag =   0       # Flag for fitting FWHM
                             #   0 = Fix FWHM to global curve
                             #   1 = Optimize FWHM
                             #   2 = Fix FWHM to input value
      .ampl_factor = 0.      # Fixed amplitude ratio to previous peak
                             #   0.  = Optimize amplitude of this peak
                             #   >0. = Fix amplitude to this value relative
                             #         to amplitude of previous free peak
                             #  -1.0 = Fix amplitude at 0.0
      .initial_energy = 0.   # Peak energy
      .initial_fwhm =   0.   # Peak FWHM
      .initial_ampl =   0.   # Peak amplitude
      
   Output fields returned by fit_peaks(), defaults and descriptions
      .energy =         0.   # Peak energy
      .fwhm =           0.   # Peak FWHM
      .ampl =           0.   # Peak amplitude
      .area =           0.   # Area of peak
      .bgd =            0.   # Background under peak
   """
   def __init__(self):
      self.label =       ""  # Peak label
      self.energy_flag = 0   # Flag for fitting energy
                             #   0 = Fix energy 
                             #   1 = Optimize energy
      self.fwhm_flag =   0   # Flag for fitting FWHM
                             #   0 = Fix FWHM to global curve
                             #   1 = Optimize FWHM
                             #   2 = Fix FWHM to input value
      self.ampl_factor = 0.  # Fixed amplitude ratio to previous peak
                             #   0.  = Optimize amplitude of this peak
                             #   >0. = Fix amplitude to this value relative
                             #         to amplitude of previous free peak
                             #  -1.0 = Fix amplitude at 0.0
      self.initial_energy = 0.  # Peak energy
      self.energy =         0.  # Peak energy
      self.initial_fwhm =   0.  # Peak FWHM
      self.fwhm =           0.  # Peak FWHM
      self.initial_ampl =   0.  # Peak amplitude
      self.ampl =           0.  # Peak amplitude
      self.area =           0.  # Area of peak
      self.bgd =            0.  # Background under peak
      self.ignore =         0   # Don't fit peak

########################################################################
class McaROI:
   """
   Class that defines a Region-Of-Interest (ROI)
   Fields
      .left      # Left channel or energy
      .right     # Right channel or energy
      .centroid  # Centroid channel or energy
      .fwhm      # Width
      .bgd_width # Number of channels to use for background subtraction
      .use       # Flag: should the ROI should be used for energy calibration
      .preset    # Is this ROI controlling preset acquisition
      .label     # Name of the ROI
      .d_spacing # Lattice spacing if a diffraction peak
      .energy    # Energy of the centroid for energy calibration
   """
   def __init__(self, left=0., right=0., centroid=0., fwhm=0., bgd_width=0,
                      use=1, preset=0, label='', d_spacing=0., energy=0.):
      """
      Keywords:
         There is a keyword with the same name as each attribute that can be
         used to initialize the ROI when it is created.
      """
      self.left = left
      self.right = right
      self.centroid = centroid
      self.fwhm = fwhm
      self.bgd_width = bgd_width
      self.use = use
      self.preset = preset
      self.label = label
      self.d_spacing = d_spacing
      self.energy = energy
   def __cmp__(self, other):
      """
      Comparison operator.  The .left field is used to define ROI ordering
      """ 
      return (self.left - other.left)

########################################################################
class McaCalibration:
   """
   Class defining an Mca calibration.  The calibration equation is
      energy = .offset + .slope*channel + .quad*channel**2
   where the first channel is channel 0, and thus the energy of the first
   channel is .offset.
   
   Fields:
      .offset    # Offset
      .slope     # Slope
      .quad      # Quadratic
      .units     # Calibration units, a string
      .two_theta # 2-theta of this Mca for energy-dispersive diffraction
   """
   def __init__(self, offset=0., slope=1.0, quad=0., units='keV', 
                      two_theta=10.):
      """
      There is a keyword with the same name as each field, so the object can
      be initialized when it is created.
      """
      self.offset = offset
      self.slope = slope
      self.quad = quad
      self.units = units
      self.two_theta = two_theta

########################################################################
class McaElapsed:
   """
   The elapsed time and counts for an Mca.
   
   Fields:
      .start_time   # Start time and date, a string
      .live_time    # Elapsed live time in seconds
      .real_time    # Elapsed real time in seconds
      .read_time    # Time that the Mca was last read in seconds
      .total_counts # Total counts between the preset start and stop channels
   """
   def __init__(self, start_time='', live_time=0., real_time=0., 
                      read_time=0., total_counts=0.):
      self.start_time = start_time
      self.live_time = live_time
      self.real_time = real_time
      self.read_time = read_time
      self.total_counts = total_counts

########################################################################
class McaPresets:
   """
   The preset time and counts for an Mca.
   
   Fields:
      .live_time       # Preset live time in seconds
      .real_time       # Preset real time in seconds
      .read_time       # Time that the Mca was last read in seconds
      .total_counts    # Preset total counts between the preset
                       #    start and stop channels
      .start_channel   # Start channel for preset counts
      .end_channel     # End channel for preset counts
      .dwell           # Dwell time per channel for MCS devices
      .channel_advance # Channel advance source for MCS hardware:
                       #    0=internal, 1=external
      .prescale        # Prescaling setting for MCS hardware
   """
   def __init__(self):
      self.live_time = 0.
      self.real_time = 0.
      self.total_counts = 0.
      self.start_channel = 0
      self.end_channel = 0
      self.dwell = 0.
      self.channel_advance = 0
      self.prescale = 0

########################################################################
class McaEnvironment:
   """
   The "environment" or related parameters for an Mca.  These might include
   things like motor positions, temperature, anything that describes the
   experiment.

   An Mca object has an associated list of McaEnvironment objects, since there
   are typically many such parameters required to describe an experiment.

   Fields:
      .name         # A string name of this parameter, e.g. "13IDD:m1"
      .value        # A string value of this parameter,  e.g. "14.223"
      .description  # A string description of this parameter, e.g. "X stage"
   """
   def __init__(self, name='', value='', description=''):
      self.name = name
      self.value = value
      self.description = description

########################################################################
class Mca(object):
   """ Device-independent MultiChannel Analyzer (MCA) class """
   def __init__(self, file=None, **filekw):
      """
      Creates new Mca object.  The data are initially all zeros, and the number
      of channels is 2048.
      
      Keywords:
         file:
            Name of a file to read into the new Mca object with read_file()
            
      Example:
         m = Mca('my_spectrum.dat')
      """
      self.name = ''
      self.n_detectors = 1
      nchans = 2048
      self.data = Numeric.zeros(nchans)
      self.rois = []
      self.calibration = McaCalibration()
      self.elapsed = McaElapsed()
      self.presets = McaPresets()
      self.environment = []
      if (file != None):
         self.read_file(file, **filekw)

   ########################################################################
   def __copy__(self):
      """
      Makes a "shallow" copy of an Mca instance, using copy.copy() on all of
      the attributes of the Mca instance.  The .rois and .environment attributes
      will still point to the same values, because they are lists.
      """
      new = Mca()
      new.data = copy.copy(self.data)
      new.rois = copy.copy(self.rois)
      new.elapsed = copy.copy(self.elapsed)
      new.calibration = copy.copy(self.calibration)
      new.presets = copy.copy(self.presets)
      new.environment = copy.copy(self.environment)
      return(new)


   ########################################################################
   def __deepcopy__(self, visit):
      """
      Makes a "deep" copy of an Mca instance, using copy.deepcopy() on all of
      the attributes of the Mca instance. All of the attribute will point to
      new objects.
      """
      new = Mca()
      new.data = copy.copy(self.data)
      new.rois = copy.deepcopy(self.rois, visit)
      new.elapsed = copy.copy(self.elapsed)
      new.calibration = copy.copy(self.calibration)
      new.presets = copy.copy(self.presets)
      new.environment = copy.deepcopy(self.environment, visit)
      return(new)

   ########################################################################
   def get_calibration(self):
      """ Returns the Mca calibration, as an McaCalibration object """
      return self.calibration

   ########################################################################
   def set_calibration(self, calibration):
      """
      Sets the Mca calibration.
      
      Inputs:
         calibration:
            An McaCalibration object
      """
      self.calibration = calibration

   ########################################################################
   def get_presets(self):
      """ Returns the Mca presets, as an McaCPresets object """
      return self.presets

   ########################################################################
   def set_presets(self, presets):
      """
      Sets the Mca presets.
      
      Inputs:
         presets:
            An McaPresets object
      """
      self.presets = presets

   ########################################################################
   def get_elapsed(self):
      """ Returns the Mca elapsed parameters, as an McaElapsed object """
      return self.elapsed

   ########################################################################
   def set_elapsed(self, elapsed):
      """
      Sets the Mca elapsed parameters.
      
      Inputs:
         elapsed:
            An McaElapsed object
      """
      self.elapsed = elapsed

   ########################################################################
   def get_name(self):
      """ Returns the Mca name as a string """
      return self.name

   ########################################################################
   def set_name(self, name):
      """
      Sets the Mca name.
      
      Inputs:
         name:
            A string
      """
      self.name = name

   ########################################################################
   def get_rois(self, energy=0):
      """ Returns the Mca ROIS, as a list of McaROI objects """
      rois = copy.copy(self.rois)
      if (energy != 0):
         for roi in rois:
            roi.left = self.channel_to_energy(roi.left)
            roi.right = self.channel_to_energy(roi.right)
      return rois

   ########################################################################
   def set_rois(self, rois, energy=0):
      """
      Sets the region-of-interest parameters for the MCA.
      The rois information is contained in an object of class McaRoi.
      This routine is not needed if the information in the McaRoi instance
      is already in channel units.  It is needed if the information in the
      .left and .right fields is in terms of energy.

      Inputs:
         rois:
            A list of objects of type McaROI
            
      Keywords:
         energy:
            Set this flag to indicate that the .left and .right fields
            of rois are in units of energy rather than channel number.
            
      Example:
        mca = Mca('mca.001')
        r1 = McaROI()
        r1.left = 5.4
        r1.right = 5.6
        r2 = McaROI()
        r2.left = 6.1
        r2.right = 6.2
        mca.set_rois([r1,r2], energy=1)
      """
      self.rois = []
      for roi in rois:
         r = roi
         if (energy == 1):
            r.left =  self.energy_to_channel(r.left, clip=1)
            r.right = self.energy_to_channel(r.right, clip=1)
         self.rois.append(r)

   ########################################################################
   def get_roi_counts(self, background_width=1):
      """
      Returns a tuple (total, net) containing the total and net counts of
      each region-of-interest in the MCA.

      Kyetwords:
         background_width:
            Set this keyword to set the width of the background region on either
            side of the peaks when computing net counts.  The default is 1.
            
      Outputs:
          total:  The total counts in each ROI.
          net:    The net counts in each ROI.

          The dimension of each list is NROIS, where NROIS
          is the number of currently defined ROIs for this MCA.  It returns
          and empty list for both if NROIS is zero.
          
      Example:
         mca = Mca('mca.001')
         total, net = mca.get_roi_counts(background_width=3)
         print 'Net counts = ', net
      """
      total = []
      net = []
      nchans = len(self.data)
      for roi in self.rois:
         left = roi.left
         ll = max((left-background_width+1), 0)
         if (background_width > 0):
             bgd_left = sum(self.data[ll:(left+1)]) / (left-ll+1)
         else: bgd_left = 0.
         right = roi.right
         rr = min((right+background_width-1), nchans-1)
         if (background_width > 0):
             bgd_right = sum(self.data[right:rr+1]) / (rr-right+1)
         else: bgd_right = 0.
         total_counts = self.data[left:right+1]
         total.append(sum(total_counts))
         n_sel        = right - left + 1
         bgd_counts   = bgd_left + Numeric.arange(n_sel,Float)/(n_sel-1) * \
                                  (bgd_right - bgd_left)
         net_counts   = total_counts - bgd_counts
         net.append(sum(net_counts))
      return (total, net)

   ########################################################################
   def get_environment(self):
      """
      Returns a list of McaEnvironment objects that contain the environment
      parameters of the Mca.
      """
      return self.environment

   ########################################################################
   def set_environment(self, environment):
      """
      Copies a list of McaEnvironment objects to the Mca object.

      Inputs:
         environment:
            A list of McaEnvironment objects.
      """
      self.environment = environment

   ########################################################################
   def get_data(self):
      """ Returns the data (counts) from the Mca """
      return self.data

   ########################################################################
   def set_data(self, data):
      """
      Copies an array of data (counts) to the Mca.

      Inputs:
         data:
            A Numeric array of data (counts).
      """
      self.data = data

   ########################################################################
   def get_energy(self):
      """
      Returns a list containing the energy of each channel in the MCA spectrum.

      Procedure:
         Simply returns mca.channel_to_energy() for each channel
         
      Example:
          from Mca import *
          mca = Mca('mca.001')
          energy = mca.get_energy()
      """
      channels = Numeric.arange(len(self.data))
      return self.channel_to_energy(channels)

   ########################################################################
   def initial_calibration(self, energy):
      """
      Performs an initial coarse energy calibration of the Mca, setting only
      the slope, and setting the offset parameter to 0.

      Inputs:
         energy: The energy of the biggest peak in the MCA spectrum.
         
      Procedure:
         This routine does the following:
            1) Sets the offset coefficient to 0.0
            2) Sets the quadratic coefficient to 0.0
            3) Determines which channel contains the most counts, PEAK_CHAN
            4) Sets the slope equal to the input energy divided by PEAK_CHAN
            
      Example:
         from Mca import *
         mca = Mca('mca.001')
         mca.initial_calibration(20.1)
      """
      peak_chan = Numeric.argmax(self.data)
      peak_chan = max(peak_chan,1)
      self.calibration.offset = 0.
      self.calibration.slope = float(energy)/peak_chan
      self.calibration.quad = 0.

   ########################################################################
   def add_roi(self, roi, energy=0):
      """
      This procedure adds a new region-of-interest to the MCA.

      Inputs:
         roi: An object of type mcaROI.
         
      Kyetwords:
         energy:
            Set this flag to 1 to indicate that the .left and .right 
            fields of roi are in units of energy rather than channel 
            number.
            
      Example:
         mca = Mca('mca.001')
         roi = McaROI()
         roi.left = 500
         roi.right = 600
         roi.label = 'Fe Ka'
         mca,add_roi(roi)
      """
      r = copy.copy(roi)
      if (energy == 1):
         r.left = self.energy_to_channel(r.left, clip=1)
         r.right = self.energy_to_channel(r.right, clip=1)
      self.rois.append(r)

      # Sort ROIs.  This sorts by left channel.
      self.rois.sort()

   ########################################################################
   def find_roi(self, left, right, energy=0):
      """
      This procedure finds the index number of the ROI with a specified
      left and right channel number.

      Inputs:
         left:
            Left channel number (or energy) of this ROI
            
         right:
            Right channel number (or energy) of this ROI
         
      Keywords:
         energy:
            Set this flag to 1 to indicate that Left and Right are in units
            of energy rather than channel number.
            
      Output:
         Returns the index of the specified ROI, -1 if the ROI was not found.
         
      Example:
         mca = Mca('mca.001')
         index = mca.find_roi(100, 200)
      """
      l = left
      r = right
      if (energy == 1):
         l = self.energy_to_channel(l, clip=1)
         r = self.energy_to_channel(r, clip=1)
      index = 0
      for roi in self.rois:
         if (l == roi.left) and (r == roi.right): return index
         index = index + 1
      return -1

   ########################################################################
   def delete_roi(self, index):
      """
      This procedure deletes the specified region-of-interest from the MCA.

      Inputs:
         index:  The index of the ROI to be deleted, range 0 to len(mca.rois)
         
      Example:
        mca = Mca('mca.001')
        mca.delete_roi(2)
      """
      del self.rois[index]

   ########################################################################
   def channel_to_energy(self, channels):
      """
      Converts channels to energy using the current calibration values for the
      Mca.  This routine can convert a single channel number or an array of
      channel numbers.  Users are strongly encouraged to use this function
      rather than implement the conversion calculation themselves, since it
      will be updated if additional calibration parameters (cubic, etc.) are
      added.

      Inputs:
         channels:
            The channel numbers to be converted to energy.  This can be
            a single number or a sequence of channel numbers.
            
      Outputs:
         This function returns the equivalent energy for the input channels.
         
      Example:
         mca = Mca('mca.001')
         channels = [100, 200, 300]
         energy = mca.channel_to_energy(channels) # Get the energy of these
      """

      c = Numeric.asarray(channels)
      return self.calibration.offset + \
             self.calibration.slope * c + \
             self.calibration.quad * Numeric.power(c, 2)

   ########################################################################
   def channel_to_d(self, channels):
      """
      Converts channels to "d-spacing" using the current calibration values for
      the Mca.  This routine can convert a single channel number or an array of
      channel numbers.  Users are strongly encouraged to use this function
      rather than implement the conversion calculation themselves, since it
      will be updated if additional calibration parameters are added.  This
      routine is useful for energy dispersive diffraction experiments.  It uses
      both the energy calibration parameters and the "two-theta" calibration
      parameter.

      Inputs:
         channels:
            The channel numbers to be converted to "d-spacing".
            This can be a single number or a list of channel numbers.
            
      Outputs:
         This function returns the equivalent "d-spacing" for the input channels.
         The output units are in Angstroms.
         
      Restrictions:
         This function assumes that the units of the energy calibration are keV
         and that the units of "two-theta" are degrees.
         
      Example:
         mca = Mca('mca.001')
         channels = [100,200,300]
         d = mca.channel_to_d(channels)       # Get the "d-spacing" of these
      """
      e = self.channel_to_energy(channels)
      return 12.398 / (2. * e * math.sin(self.calibration.two_theta/2.*math.pi/180.))
  
    #FIXED : added divide by 2 to conversion

   ########################################################################
   def energy_to_channel(self, energy, clip=0):
      """
      Converts energy to channels using the current calibration values for the
      Mca.  This routine can convert a single energy or an array of energy
      values.  Users are strongly encouraged to use this function rather than
      implement the conversion calculation themselves, since it will be updated
      if additional calibration parameters are added.

      Inputs:
         energy:
            The energy values to be converted to channels. This can be a
            single number or a sequence energy values.
            
      Keywords:
         clip:
            Set this flag to 1 to clip the returned values to be between
            0 and nchans-1.  The default is not to clip.
            
      Outputs:
         This function returns the closest equivalent channel for the input
         energy.  Note that it does not generate an error if the channel number
         is outside the range 0 to (nchans-1), which will happen if the energy
         is outside the range for the calibration values of the Mca.
         
      Example:
         mca = Mca('mca.001')
         channel = mca.energy_to_channel(5.985)
      """
      if (self.calibration.quad == 0.0):
         channel = ((energy-self.calibration.offset) /
                    self.calibration.slope)
      else:
         # Use the quadratic formula, use some shorthand
         a = self.calibration.quad
         b = self.calibration.slope
         c = self.calibration.offset - energy
         # There are 2 roots.  I think we always want the "+" root?
         channel = (-b + Numeric.sqrt(b**2 - 4.*a*c))/(2.*a)
      channel = Numeric.around(channel)
      if (clip != 0): 
         nchans = len(self.data)
         channel = Numeric.clip(channel, 0, nchans-1)
      if (type(channel) == Numeric.ArrayType): 
         return channel.astype(Numeric.Int)
      else:
         return int(channel)

   ########################################################################
   def d_to_channel(self, d, clip=0):
      """
      Converts "d-spacing" to channels using the current calibration values
      for the Mca.  This routine can convert a single "d-spacing" or an array
      of "d-spacings".  Users are strongly encouraged to use this function
      rather than implement the conversion calculation themselves, since it
      will be updated if additional calibration parameters are added.
      This routine is useful for energy dispersive diffraction experiments.
      It uses both the energy calibration parameters and the "two-theta"
      calibration parameter.

      Inputs:
         d:
            The "d-spacing" values to be converted to channels.
            This can be a single number or an array of values.
            
      Keywords:
         clip:
            Set this flag to 1 to clip the returned values to be between
            0 and nchans-1.  The default is not to clip.
            
      Outputs:
         This function returns the closest equivalent channel for the input
         "d-spacing". Note that it does not generate an error if the channel
         number is outside the range 0 to (nchans-1), which will happen if the
         "d-spacing" is outside the range for the calibration values of the Mca.
         
      Example:
         mca = Mca('mca.001')
         channel = mca.d_to_chan(1.598)
      """
      e = 12.398 / (2. * d * math.sin(self.calibration.two_theta*math.pi/180./2.))
      return self.energy_to_channel(e, clip=clip)

   ########################################################################
   def write_file(self, file, netcdf=0):
      """
      Writes Mca or Med objects to a disk file.
      
      It calls Mca.write_netcdf_file if the netcdf keyword flg is set,

      Note that users who want to read such files with Python are strongly
      encouraged to use Mca.read_file()

      Inputs:
         file:
            The name of the disk file to write.
            
      Keywords:
         netcdf:
            Set this flag to write the file in netCDF format, otherwise
            the file is written in ASCII format.  See the documentation
            for Mca.write_ascii_file and Mca.write_netcdf_file for 
            information on the formats.
 
      Example:
         mca = Mca()
         mca.write_file('mca.001')
      """
      # Call the get_xxx() methods to make sure things are up to date
      data = self.get_data()
      calibration = self.get_calibration()
      elapsed = self.get_elapsed()
      presets = self.get_presets()
      rois = self.get_rois()
      environment = self.get_environment()

      if (netcdf != 0):
        write_netcdf_file(file, data, calibration, elapsed, presets, rois,
                          environment)
      else:
        write_ascii_file(file, data, calibration, elapsed, presets, rois,
                         environment)



   ########################################################################
   def read_file(self, file, netcdf=0, detector=0):
      """
      Reads a disk file into an MCA object.  If the netcdf=1 flag is set it
      reads a netcdf file, else it assumes the file is ASCII.
      If the data file has multiple detectors then the detector keyword can be
      used to specify which detector data to return.

      Inputs:
         file:
            The name of the disk file to read.
            
      Keywords:
         netcdf:
            Set this flag to read files written in netCDF format, otherwise
            the routine assumes that the file is in ASCII format.
            See the documentation for Mca.write_ascii_file and
            Mca.write_netcdf_file for information on the formats.

         detector:
            Specifies which detector to read if the file has multiple detectors.
            
      Example:
         mca = Mca()
         mca.read_file('mca.001')
      """
      if (netcdf != 0):
         r = read_netcdf_file(file)
      else:
         r = read_ascii_file(file)
      self.name = file
      self.calibration = r['calibration'][detector]
      self.data = r['data'][detector]
      self.elapsed = r['elapsed'][detector]
      self.rois = r['rois'][detector]
      self.environment = r['environment']



   ########################################################################
#   def fit_background(self, bottom_width=4., top_width=0., exponent=2, 
#                      tangent=0, compress=4):
#      """
#      This function fits a background to an MCA spectrum. The background is
#      fitted using an enhanced version of the algorithm published by
#      Kajfosz, J. and Kwiatek, W .M. (1987)  "Non-polynomial approximation of
#      background in x-ray spectra." Nucl. Instrum. Methods B22, 78-81.
#      
#      Keywords:
#         top_width:
#            Specifies the width of the polynomials which are concave upward.
#            The top_width is the full width in energy units at which the
#            magnitude of the polynomial is 100 counts. The default is 0, which
#            means that concave upward polynomials are not used.
#
#         bottom_width:
#            Specifies the width of the polynomials which are concave downward.
#            The bottom_width is the full width in energy units at which the
#            magnitude of the polynomial is 100 counts. The default is 4.
#
#         exponent:
#            Specifies the power of polynomial which is used. The power must be
#            an integer. The default is 2, i.e. parabolas. Higher exponents,
#            for example EXPONENT=4, results in polynomials with flatter tops
#            and steeper sides, which can better fit spectra with steeply
#            sloping backgrounds.
#
#         tangent:
#            Specifies that the polynomials are to be tangent to the slope of the
#            spectrum. The default is vertical polynomials. This option works
#            best on steeply sloping spectra. It has trouble in spectra with
#            big peaks because the polynomials are very tilted up inside the
#            peaks.
#
#         compress:
#            Compression factor to apply before fitting the background.
#            Default=4, which means, for example, that a 2048 channel spectrum
#            will be rebinned to 512 channels before fitting.
#            The compression is done on a temporary copy of the input spectrum,
#            so the input spectrum itself is unchanged.
#            The algorithm works best if the spectrum is compressed before it
#            is fitted. There are two reasons for this. First, the background
#            is constrained to never be larger than the data itself. If the
#            spectrum has negative noise spikes they will cause the fit to be
#            too low. Compression will smooth out such noise spikes.
#            Second, the algorithm requires about 3*N^2 operations, so the time
#            required grows rapidly with the size of the input spectrum. On a
#            200 MHz Pentium it takes about 3 seconds to fit a 2048 channel
#            spectrum with COMPRESS=1 (no compression), but only 0.2 seconds
#            with COMPRESS=4 (the default).
#            
#     Procedure:
#         1) At each channel "i" an n'th degree polynomial which is concave up
#         is fitted. Its equation is
#
#                                       n
#                          (e(i) - e(j))
#         f(j,i) = y(i) + --------------
#                                    n
#                           top_width
#
#         where f(j,i) is the fitted counts in channel j for the polynomial
#         centered in channel i. y(i) is the input counts in channel "i", e(i) is
#         the energy of channel i, e(j) is the energy of channel j, and
#         "top_width" and "n" are user-specified parameters. The background count
#         in channel "j", b(j) is defined as
#
#         b(j) = min ((f(j,i), y(j))
#                 i
#
#         b(j) is thus the smallest fitted polynomial in channel j, or the raw
#         data, whichever is smaller.
#
#         2) After the concave up polynomials have been fitted, a series of
#         concave down polynomials are constructed. At each channel "i" an n'th
#         degree polynomial which is concave up is fitted. The polynomial is slid
#         up from below until it "just touches" some channel of the spectrum. Call
#         this channel "i". The maximum height of the polynomial is thus
#
#                                                  n
#                                     (e(i) - e(j))
#         height(j) = max ( b(j) +  --------------  )
#                      i                          n
#                                     bottom_width
#
#         where bottom_width is a user_specified parameter.
#
#         3) Once the value of height(i) is known the polynomial is fitted. The
#         background counts in each channel are then determined from:
#
#                                                  n
#                                     (e(i) - e(j))
#         bgd(j) = max ( height(i) + --------------
#                   i                             n
#                                     bottom_width
#
#         bgd(j) is thus the maximum counts for any of the concave down
#         polynomials passing though channel j.
#
#         Before the concave-down polynomials are fitted the spectrum at each
#         channel it is possible to subtract out a straight line which is
#         tangent to the spectrum at that channel. Use the /TANGENT qualifier to
#         do this. This is equivalent to fitting a "tilted" polynomial whose
#         apex is tangent to the spectrum at that channel. By fitting
#         polynomials which are tangent rather than vertical the background fit
#         is much improved on spectra with steep slopes.
#
#      Outputs:
#         This function returns an MCA object which is identical to the calling
#         object, except that the data have been replaced by the background fit.
#
#      Example:
#        mca = Mca()
#        mca.read_file('mca.001')
#        bgd = mca.fit_background(mca, bottom=6, exponent=4)
#      """
#      REFERENCE_AMPL=100.
#      TINY = 1.E-20
#      HUGE = 1.E20
#      MAX_TANGENT=2
#
#      bgd = copy.copy(self)
#      nchans = len(bgd.data)
#      calibration = bgd.get_calibration()
#      scratch = copy.copy(bgd.get_data())
#      slope = calibration.slope
#
#      # Compress scratch spectrum
#      if (compress > 1):
#         scratch = CARSMath.compress_array(scratch, compress)
#         slope = slope * compress
#         nchans = nchans / compress
#
#      # Copy scratch spectrum to background spectrum
#      bckgnd = copy.copy(scratch)
#
#      # Find maximum counts in input spectrum. This information is used to
#      # limit the size of the function lookup table
#      max_counts = max(scratch)
#
#      #  Fit functions which come down from top
#      if (top_width > 0.):
#         #   First make a lookup table of this function
#         chan_width = top_width / (2. * slope)
#         denom = chan_width**exponent
#         indices = Numeric.arange(float(nchans*2+1)) - nchans
#         power_funct = indices**exponent * (REFERENCE_AMPL / denom)
#         power_funct = Numeric.compress((power_funct <= max_counts), power_funct)
#         max_index = len(power_funct)/2 - 1
#
#         for center_chan in range(nchans):
#            first_chan = max((center_chan - max_index), 0)
#            last_chan = min((center_chan + max_index), (nchans-1))
#            f = first_chan - center_chan + max_index
#            l = last_chan - center_chan + max_index
#            test = scratch[center_chan] + power_funct[f:l+1]
#            sub = bckgnd[first_chan:last_chan+1] 
#            bckgnd[first_chan:last_chan+1] = Numeric.maximum(sub, test)
#
#      # Copy this approximation of background to scratch
#      scratch = copy.copy(bckgnd)
#
#      # Find maximum counts in scratch spectrum. This information is used to
#      #   limit the size of the function lookup table
#      max_counts = max(scratch)
#
#      # Fit functions which come up from below
#      bckgnd = Numeric.arange(float(nchans)) - HUGE
#
#      # First make a lookup table of this function
#      chan_width = bottom_width / (2. * slope)
#      if (chan_width == 0.): denom = TINY
#      else: denom = chan_width**exponent
#      indices = Numeric.arange(float(nchans*2+1)) - nchans
#      power_funct = indices**exponent  * (REFERENCE_AMPL / denom)
#      power_funct = Numeric.compress((power_funct <= max_counts), power_funct)
#      max_index = len(power_funct)/2 - 1
#
#      for center_chan in range(nchans-1):
#         tangent_slope = 0.
#         if (tangent):
#            # Find slope of tangent to spectrum at this channel
#            first_chan = max((center_chan - MAX_TANGENT), 0)
#            last_chan = min((center_chan + MAX_TANGENT), (nchans-1))
#            denom = center_chan - Numeric.arange(float(last_chan - first_chan + 1))
#            tangent_slope = (scratch[center_chan] - 
#                             scratch[first_chan:last_chan+1]) / max(denom, 1)
#            tangent_slope = Numeric.sum(tangent_slope) / (last_chan - first_chan)
#
#         first_chan = max((center_chan - max_index), 0)
#         last_chan = min((center_chan + max_index), (nchans-1))
#         last_chan = max(last_chan, first_chan)
#         nc = last_chan - first_chan + 1
#         lin_offset = scratch[center_chan] + \
#                     (Numeric.arange(float(nc)) - nc/2) * tangent_slope
#
#         # Find the maximum height of a function centered on this channel
#         # such that it is never higher than the counts in any channel
#
#         f = first_chan - center_chan + max_index
#         l = last_chan - center_chan + max_index
#         test = scratch[first_chan:last_chan+1] - lin_offset + \
#                                                            power_funct[f:l+1]
#         height = min(test)
#
#         # We now have the function height. Set the background to the
#         # height of the maximum function amplitude at each channel
#
#         test = height + lin_offset - power_funct[f:l+1]
#         sub = bckgnd[first_chan:last_chan+1]
#         bckgnd[first_chan:last_chan+1] = Numeric.maximum(sub, test)
#
#      # Expand spectrum
#      if (compress > 1): bckgnd = CARSMath.expand_array(bckgnd, compress)
#      bgd.set_data(bckgnd.astype(Numeric.Int))
#      return bgd
#

   ########################################################################
#   def fit_peaks(self, peaks, fit=None, background=None,
#                 output='', spreadsheet=None, 
#                 append=1, **background_kw):
#      """
#      Fits the peaks in the MCA spectrum. It provides a convenient interface to
#      the fitPeaks() function.
#   
#      Inputs:
#         peaks:  A list of McaPeak objects.  See fitPeaks() and read_peaks() for
#         more information.
#
#      Keywords:
#         fit:
#            An object of type McaFit which can be used to control the
#            peak fitting.  If this keyword is omitted then the fit structure
#            is created with McaFit()
#            
#         background:
#            An Mca object containing the fitted background.  If this keyword
#            is omitted then this function will call Mca.fit_background() before
#            calling fitPeaks().
#            
#         output:
#            The name of an output file to receive the ASCII printout of the
#            fit results.  This keyword is simply passed to fit_peaks_report().
#            
#         spreadsheet:
#            The name of an output file to receive the ASCII output of the
#            fit results in spreadsheet format.  This keyword is simply passed to
#            fit_peaks_report().
#            
#         append:
#            Flag indicating whether the output and spreadsheet files should be
#            appended to or overwritten.  This keyword is simply passed to
#            fit_peaks_report9).
#
#         In addition to these keywords, all keywords accepted by the
#         fit_background() function are accepted if the background keyword is
#         not present, i.e. if this function will be calling fit_background().
#
#      Outputs:
#         This function returns an Mca object which is identical to the calling
#         object, except that the data have been replaced by the peak fit.
#
#      Procedure:
#         The function does the following:
#            - Creates the Fit structure with mca->FIT_INITIALIZE() if Fit
#              was not passed as a keyword parameter.
#            - Fits the background using fit_background() if background
#              was not passed as a keyword parameter.
#            - Extracts the data from the input spectrum and the background
#              spectrum.
#            - Calls fitPeaks() with the background subtracted data.
#            - Calls fit_peaks_report()
#            - Creates a new Mca object using Mca.deepcopy and stores the output
#              of fitPeaks() in this new object with set_data().  It then
#              returns this new MCA object as the function return value.
#
#       Example:
#             mca = Mca()
#             mca.read_file('mca.001')
#             peaks = read_peaks('mypeaks.pks')
#             fit = mca.fit_peaks(peaks, bottom=6, exponent=4)
#      """
#      # If fit is not an object of type McaFit then initialize it
#      if (not isinstance(fit, McaFit)):
#         fit = McaFit(self)
#
#      # Copy the calibration parameters to fit
#      fit.update(self)
#
#      # If background is not an object of type Mca then initialize it
#      if (not isinstance(background, Mca)):
#          background = self.fit_background(**background_kw)
#
#      fit.npeaks = len(peaks)
#      background_counts = background.get_data()
#      observed_counts = self.get_data()
#      t0 = time.time()
#      [fit, peaks, fit_counts] = fitPeaks.fitPeaks(fit, peaks,
#                                      observed_counts - background_counts)
#      t1 = time.time()
#      fit_counts = fit_counts + background_counts
#      self.fit_peaks_report(fit, peaks, background, output=output, 
#                       spreadsheet=spreadsheet, append=append, time=t1-t0)
#      fit_mca = copy.copy(self)
#      fit_mca.set_data(fit_counts)
#      return([fit, peaks, fit_mca])
#
#

#   ########################################################################
#   def fit_peaks_report(self, fit, peaks, background, output='',
#                        spreadsheet=None, append=1, time=None):
#      """
#      Prints out the results from <A HREF="mca_utility_routines.html#FIT_PEAKS">FIT_PEAKS</A>
#
#      Inputs:
#         fit:
#            An McaFit object with the global fitting parameters.
#
#         peaks:
#            A list of McaPeak objects with the fit results for each peak.
#
#         See fit_peaks for more information on fit and peaks.
#
#      background:
#         An Mca object containing the fitted background spectrum.
#
#      Keywords:
#         output:
#            The name of an output file to receive the ASCII printout of the
#            fit results.  If this keyword is omitted then the output will be
#            written to stdout, i.e. the IDL output window.  If the Output file
#            already exists then the new information will (by default) be appended
#            to the file.
#
#         spreadsheet:
#            The name of an output file to receive the ASCII output of the
#            fit results in a format easily imported into a spreadsheet.  If this
#            keyword is omitted then no spreadsheet output will be generated.
#            written to stdout, i.e. the IDL output window.
#            If the spreadhseet file already exists then the new information will
#            (by default) be appended to the file.
#
#         append:
#            Set this keyword to 0 to overwrite the output and spreadsheet files
#            rather than to append to them, which is the default behavior.
#
#      Example:
#         mca = Mca(file='mca.001')
#         peaks = read_peaks('mypeaks.pks')
#         [fit, peaks, predicted] = mca.fit_peaks(peaks, fit,
#                                                 bottom=6, exponent=4)
#         mca.fit_peaks_report(fit, peaks, background, output='mca.001_out')
#      """
#
#      if (append): mode = 'a'
#      else: mode = 'w'
#      if (output == ''): out_fp = sys.stdout
#      else: out_fp = open(output, mode)
#      if (spreadsheet != None): spread_fp = open(spreadsheet, mode)
#      else: spread_fp = None
#
#      SIGMA_TO_FWHM = 2.*Numeric.sqrt(2.*Numeric.log(2.))
#
#      # Compute backgrounds
#      background_counts = background.get_data()
#      nchans = len(self.get_data())
#      for peak in peaks:
#         low = background.energy_to_channel(peak.energy - 
#                 2.*peak.fwhm / SIGMA_TO_FWHM)
#         low = min(max(low, 0), nchans-3)
#         hi  = background.energy_to_channel(peak.energy + 
#                 2.*peak.fwhm / SIGMA_TO_FWHM)
#         hi = min(max(hi, low+1), nchans-1)
#         peak.bgd = Numeric.sum(background_counts[low:hi+1])
#
#      if (out_fp != None):
#         out_fp.write('\n')
#         out_fp.write('*******************************************************\n')
#         out_fp.write( '    Fit of ' + self.name + '\n')
#         out_fp.write('\n')
#         elapsed = self.get_elapsed()
#         out_fp.write('Real time (seconds):            ' +
#                        ('%10.2f' % elapsed.real_time) + '\n')
#         out_fp.write('Live time (seconds):            ' +
#                        ('%10.2f' % elapsed.live_time) + '\n')
#         out_fp.write('Initial FWHM offset, slope:     ' +
#                        ('%10.6f' % fit.initial_fwhm_offset) + 
#                        ('%10.6f' % fit.initial_fwhm_slope) + '\n')
#         out_fp.write('Optimized FWHM offset, slope:   ' +
#                        ('%10.6f' % fit.fwhm_offset) +
#                        ('%10.6f' % fit.fwhm_slope) +  '\n')
#         out_fp.write('Initial energy offset, slope:   ' +
#                        ('%10.6f' % fit.initial_energy_offset) +
#                        ('%10.6f' % fit.initial_energy_slope) + '\n')
#         out_fp.write('Optimized energy offset, slope: ' +
#                        ('%10.6f' % fit.energy_offset) +
#                        ('%10.6f' % fit.energy_slope) + '\n')
#         out_fp.write('# Iterations, function evals:   ' +
#                        ('%10d' % fit.n_iter) +
#                        ('%10d' % fit.n_eval) + '\n')
#         out_fp.write('Chi squared:                    ' +
#                        ('%.6g' % fit.chisqr) + '\n')
#         out_fp.write('Status code:                    ' +
#                        ('%d' % fit.status) + '\n')
#         if (fit.status <= 0):
#            out_fp.write('Error message:               ' +
#                           fit.err_string)
#         if (time != None):
#            out_fp.write('Time to fit:                 ' +
#                        ('%.3f' % time) + '\n')
#             
#         out_fp.write('\n')
#         out_fp.write('        Peak       Energy    FWHM      Area     ' + 
#                'Background   Area/MDL   Area/Bkg\n')
#         out_fp.write('\n')
#         for peak in peaks:
#            if (peak.energy_flag == 0):  esym=' '
#            else: esym='*'
#            if (peak.fwhm_flag == 0):    fsym=' '
#            else: fsym='*'
#            if (peak.ampl_factor == 0.): asym=' '
#            else: asym='*'
#            out_fp.write(
#                  ('%15s' % peak.label) +
#                  ('%10.3f%1s' % (peak.energy, esym)) +
#                  ('%10.4f%1s' % (peak.fwhm, fsym)) +
#                  ('%10.1f%1s' % (peak.area, asym)) +
#                  ('%10.1f' % peak.bgd) +
#                  ('%10.1f' % (peak.area/
#                               max((3.*Numeric.sqrt(peak.bgd)), 1.0))) +
#                  ('%10.1f' % (peak.area/max(peak.bgd, 1.0)))+'\n')
#
#      if (spread_fp != None):
#         el = ('%10.2f' % elapsed.live_time)
#         er = ('%10.2f' % elapsed.real_time)
#         spread_fp.write(self.name+'#Labels#Live time#Real time#')
#         for peak in peaks:
#            spread_fp.write(peak.label+'#')
#         spread_fp.write('\n')
#         spread_fp.write(self.name+'#Energy#'+el+'#'+er+'#')
#         for peak in peaks:
#            spread_fp.write(('%10.3f' % peak.energy)+'#')
#         spread_fp.write('\n')
#         spread_fp.write(self.name+'#FWHM#'+el+'#'+er+'#')
#         for peak in peaks:
#            spread_fp.write(('%10.4f' % peak.fwhm)+'#')
#         spread_fp.write('\n')
#         spread_fp.write(self.name+'#Area#'+el+'#'+er+'#')
#         for peak in peaks:
#            spread_fp.write(('%10.1f' % peak.area)+'#')
#         spread_fp.write('\n')
#         spread_fp.write(self.name+'#Background#'+el+'#'+er+'#')
#         for peak in peaks:
#            spread_fp.write(('%10.1f' % peak.bgd)+'#')
#         spread_fp.write('\n')
#         spread_fp.close()

#######################################################################
def write_ascii_file(file, data, calibration, elapsed, presets, rois,
                     environment):
   """
   Writes Mca or Med data to a disk file.  The file 
   format is a tagged ASCII format.  The file contains the information 
   from the Mca object which it makes sense to store permanently, but 
   does not contain all of the internal state information for the Mca.  
   Files written with this routine can be read with read_ascii_file(), which
   is called by Mca.read_file() if the netcdf flag is 0.

   This procedure is typically not called directly, but is called
   by Mca.write_file if the netcdf=1 keyword is not used.

   This function can be used for writing for Mca objects, in which case
   each input parameter is an object, such as McaElapsed, etc.
   It can also be used for writing Med objects, in which case each input
   parameter is a list.
   
   If the rank of data is 2 then this is an Med, and the number of detectors
   is the first dimension of data

   Inputs:
      file:
         The name of the disk file to write.
         
      data:
         The data to write.  Either 1-D array or list of 1-D arrays.

      calibration:
         An object of type McaCalibration, or a list of such objects.

      elapsed:
         An object of type McaElapsed, or a list of such objects.

      presets:
         An object of type McaPresets, or a list of such objects.

      rois:
         A list of McaROI objects, or a list of lists of such objects.
      
      environment:
         A list of McaEnvironment objects, or a list of lists of such objects.
   """
   if (Numeric.rank(data) == 2):
      n_det = len(data)
   else:
      n_det = 1
   fformat = '%f ' * n_det
   eformat = '%e ' * n_det
   iformat = '%d ' * n_det
   sformat = '%s ' * n_det
   if (n_det == 1):
      # For convenience we convert all attributes to lists
      data = [data]
      rois = [rois]
      calibration = [calibration]
      presets = [presets]
      elapsed = [elapsed]
   nchans = len(data[0])
   start_time = elapsed[0].start_time

   fp = open(file, 'w')
   fp.write('VERSION:    '+'3.1'+'\n')
   fp.write('ELEMENTS:   '+str(n_det)+'\n')
   fp.write('DATE:       '+str(start_time)+'\n')
   fp.write('CHANNELS:   '+str(nchans)+'\n')

   nrois = []
   for roi in rois:
      nrois.append(len(roi))
   fp.write('ROIS:       '+(iformat % tuple(nrois))+'\n')
   real_time=[]; live_time=[]
   for e in elapsed:
      real_time.append(e.real_time)
      live_time.append(e.live_time)
   fp.write('REAL_TIME:  '+(fformat % tuple(real_time))+'\n')
   fp.write('LIVE_TIME:  '+(fformat % tuple(live_time))+'\n')
   offset=[]; slope=[]; quad=[]; two_theta=[]
   for c in calibration:
      offset.append(c.offset)
      slope.append(c.slope)
      quad.append(c.quad)
      two_theta.append(c.two_theta)
   fp.write('CAL_OFFSET: '+(eformat % tuple(offset))+'\n')
   fp.write('CAL_SLOPE: '+(eformat % tuple(slope))+'\n')
   fp.write('CAL_QUAD: '+(eformat % tuple(quad))+'\n')
   fp.write('TWO_THETA: '+(fformat % tuple(two_theta))+'\n')

   for i in range(max(nrois)):
      num = str(i)
      left=[]; right=[]; label=[]
      for d in range(n_det):
         if (i < nrois[d]):
            left.append(rois[d][i].left)
            right.append(rois[d][i].right)
            label.append(rois[d][i].label + '&')
         else:
            left.append(0)
            right.append(0)
            label.append(' &')
      fp.write('ROI_'+num+'_LEFT:  '+(iformat % tuple(left))+'\n')
      fp.write('ROI_'+num+'_RIGHT:  '+(iformat % tuple(right))+'\n')
      fp.write('ROI_'+num+'_LABEL:  '+(sformat % tuple(label))+'\n')
   for e in environment:
     fp.write('ENVIRONMENT: '       + str(e.name) +
                              '="'  + str(e.value) +
                              '" (' + str(e.description) + ')\n')
   fp.write('DATA: \n')
   counts = Numeric.zeros(n_det)
   for i in range(nchans):
      for d in range(n_det):
         counts[d]=data[d][i]
      fp.write((iformat % tuple(counts))+'\n')
   fp.close()

########################################################################
def read_ascii_file(file):
   """
   Reads a disk file.  The file format is a tagged ASCII format.
   The file contains the information from the Mca object which it makes sense
   to store permanently, but does not contain all of the internal state
   information for the Mca.  This procedure reads files written with
   write_ascii_file().

   Inputs:
      file:
         The name of the disk file to read.
         
   Outputs:
      Returns a dictionary of the following type:
      'n_detectors': int,
      'calibration': [McaCalibration()],
      'elapsed':     [McaElapsed()],
      'rois':        [[McaROI()]]
      'data':        [Numeric.array]
      'environment': [[McaEnvironment()]]
      
   Example:
      m = read_ascii_file('mca.001')
      m['elapsed'][0].real_time
   """
   fp = open(file, 'r')
   line = ''
   start_time = ''
   data = None
 
   environment = []
   n_detectors = 1  # Assume single element data
   elapsed = [McaElapsed()]
   calibration = [McaCalibration()]
   rois = [[]]
   while(1):
      line = fp.readline()
      if (line == ''): break
      pos = line.find(' ') #pos = string.find(line, ' ')
      if (pos == -1): pos = len(line)
      tag = line[0:pos]
      value = str.strip(line[pos:])#line[pos:].strip #value = string.strip(line[pos:])
      values = str.split(value) #values = string.split(value)
      if (tag == 'VERSION:'):
          pass
      elif (tag == 'DATE:'):  
         start_time = value
      elif (tag == 'ELEMENTS:'):
         n_detectors  = int(value)
         for det in range(1, n_detectors):
             elapsed.append(McaElapsed())
             calibration.append(McaCalibration())
             rois.append([])
      elif (tag == 'CHANNELS:'):
         nchans = int(value)
      elif (tag == 'ROIS:'):
         nrois = []
         for d in range(n_detectors):
            nrois.append(int(values[d]))
         max_rois = max(nrois)
         for d in range(n_detectors):
            for r in range(nrois[d]):
               rois[d].append(McaROI())
      elif (tag == 'REAL_TIME:'):
         for d in range(n_detectors):
            elapsed[d].start_time = start_time
            elapsed[d].real_time = float(values[d])
      elif (tag == 'LIVE_TIME:'):  
         for d in range(n_detectors):
            elapsed[d].live_time = float(values[d])
      elif (tag == 'CAL_OFFSET:'):
         for d in range(n_detectors):
            calibration[d].offset = float(values[d])
      elif (tag == 'CAL_SLOPE:'):
         for d in range(n_detectors):
            calibration[d].slope = float(values[d])
      elif (tag == 'CAL_QUAD:'):  
         for d in range(n_detectors):
            calibration[d].quad = float(values[d])
      elif (tag == 'TWO_THETA:'):
         for d in range(n_detectors):
            calibration[d].two_theta = float(values[d])
      elif (tag == 'ENVIRONMENT:'):
         env = McaEnvironment()
         p1 = str.find(value, '=')
         env.name = value[0:p1]
         p2 = str.find(value[p1+2:], '"')
         env.value = value[p1+2: p1+2+p2]
         env.description = value[p1+2+p2+3:-1]
         environment.append(env)
      elif (tag == 'DATA:'):
         data = []
         for d in range(n_detectors):
            data.append(Numeric.zeros(nchans, 'i'))
         for chan in range(nchans):
            line = fp.readline()
            counts = str.split(line)
            for d in range(n_detectors):
               data[d][chan]=int(counts[d])
      else:
         for i in range(max_rois):
             roi = 'ROI_'+str(i)+'_'
             if (tag == roi+'LEFT:'):
                for d in range(n_detectors):
                   if (i < nrois[d]):
                       rois[d][i].left = int(values[d])
                break
             elif (tag == roi+'RIGHT:'):
                for d in range(n_detectors):
                   if (i < nrois[d]):
                      rois[d][i].right = int(values[d])
                break
             elif (tag == roi+'LABEL:'):
                labels = str.split(value, '&')
                for d in range(n_detectors):
                   if (i < nrois[d]):
                      rois[d][i].label = str.strip(labels[d])
                break
         else:
            print('Unknown tag = '+tag+' in file: ' + file + '.')

   # Make sure DATA array is defined, else this was not a valid data file
   if (data == None): print('Not a valid data file: ' + file + '.')
   fp.close()
   # Built dictionary to return
   r = {}
   r['n_detectors'] = n_detectors
   r['calibration'] = calibration
   r['elapsed'] = elapsed
   r['rois'] = rois
   r['data'] = data
   r['environment'] = environment
   return r

########################################################################
#def read_peaks(file):
#   """
#   Reads a disk file into an array of structures of type 
#   McaPeak.  This routine is typically called before calling fit_peaks.
#   The routine also returns a structure of type McaBackground.  This structure
#   may or may not actually be defined in the file (older peaks files lacked it),
#   but a reasonable default value will always be returned.
#
#   Inputs:
#      file:  
#         The name of a disk file containing the peak definitions.
#         
#   Outputs:
#        This function returns a dictionary:
#        'peaks': [McaPeak]          # A list of McaPeak objects
#        'background': McaBackground # An McaBackground object
#        
#   The format of the disk file is as follows:
#      - Lines containing the parameters for the background (fields in
#        McaBackground class) have the following format:
#           Background_exponent, 4
#           Background_top_width, 0
#           Background_bottom_width, 4
#           Background_tangent, 0
#           Background_compress, 8
#      - There is one line in the file for each peak
#           - Each line consists of the following fields, separated by commas:
#             energy, energy_flag, fwhm, fwhm_flag, ampl_factor, label
#           - All fields except the first, "energy", are optional and default
#             values of 0 or blank.
#           - Field descriptions
#               - energy:       This field can either be an energy in keV or a
#                               string which can be parsed by Xrf.lookup_xrf_line()
#                               or Xrf.lookup_gamma_line()
#                               The energy in keV of the peak is put in the 
#                               .energy value for the peak.
#               - energy_flag:  The .energy_flag value for the peak. 0 or 1.
#               - fwhm:         The .initial_fwhm value for the peak.  This 
#                               can be 0 if the .fwhm_flag is 0.
#               - fwhm_flag:    The .fwhm_flag value for the peak. 0, 1 or 2.
#               - ampl_factor:  The .ampl_factor for the peak.
#               - label:        The .label string for the peak.  If the energy
#                               field is a string, and the label field is blank
#                               then the energy string will be put in .label.
#
#   The following is an example peak file. Most peaks in this example use 
#   the default (0) values for the energy_flag, fwhm_flag and
#   amplitude_factor fields.
#      4.660,,,,,Fe escape  ! This peak uses a numeric energy value
#      Fe Ka                ! This peak uses a string energy value
#      Fe Kb,,,,.5          ! Fe Kb is constrained to be 0.5 of Fe Ka
#      Ni Ka                ! These peaks energies are found with
#      Ni Kb                !  LOOKUP_XRF_LINE
#      Co57 G2              ! These peak energies are found with
#      Cd109 G1             !  LOOKUP_GAMMA_LINE
#      15.9,1,.3,1,,Diffraction ! Fit both the energy and fwhm of this peak
#      17.443,1,,,,Unknown  ! Fit the energy, but not fwhm, of this peak
#   """
#
#   fp = open(file, 'r')
#   lines = fp.readlines()
#   fp.close()
#   peaks = []
#   background = McaBackground()
#   for s in lines:
#      s = s.strip()
#      q = s.find(',')
#      if (q >= 0): 
#         label = s[0:q].upper()
#         if (label == "BACKGROUND_EXPONENT"):
#            background.exponent = int(float(s[q+1:]))
#         elif (label == "BACKGROUND_TOP_WIDTH"):
#            background.top_width = float(s[q+1:])
#         elif (label == "BACKGROUND_BOTTOM_WIDTH"):
#            background.bottom_width = float(s[q+1:])
#         elif (label == "BACKGROUND_TANGENT"):
#            background.tangent = int(s[q+1:])
#         elif (label == "BACKGROUND_COMPRESS"):
#            background.compress = int(s[q+1:])
#         else:
#            peaks.append(parse_peak(s))
#      else:
#         peaks.append(parse_peak(s))
#   return {'peaks': peaks, 'background': background}
#
#########################################################################
#def parse_peak(text):
#   """ Private function """
#   peak = McaPeak()
#   params = text.split(',')
#   n_params = len(params)
#   for i in range(n_params):
#      params[i] = params[i].strip()
#   if (n_params == 0): return peak
#
#   try:
#      peak.initial_energy = float(params[0])
#   except:
#      peak.label = params[0]
#      e = Xrf.lookup_xrf_line(params[0])
#      if (e == None):
#         e = Xrf.lookup_gamma_line(params[0])
#      if (e == None): e = 0.
#      peak.initial_energy = e
#   if (n_params == 1): return peak
#
#   if (len(params[1]) > 0):
#      peak.energy_flag = int(float(params[1]))
#   if (n_params == 2): return peak
#
#   if (len(params[2]) > 0):
#      peak.intial_fwhm = float(params[2])
#   if (n_params == 3): return peak
#
#   if (len(params[3]) > 0):
#      peak.fwhm_flag = int(float(params[3]))
#   if (n_params == 4): return peak
#
#   if (len(params[4]) > 0):
#      peak.ampl_factor = float(params[4])
#   if (n_params == 5): return peak
#
#   if (len(params[5]) > 0):
#      peak.label = params[5]
#   return peak

########################################################################
def write_peaks(file, peaks, background=None):
   """
   Writes a list of obejcts of type McaPeak to a disk file.
   If the background parameter is present it also writes the background
   structure to the file.
   
   Inputs:
      file:
         The name of a disk file to be written ;
      peaks:
         A list of McaPeak objects
         
   Keywords:
      background:
         An object of type McaBackground
         
   Example:
      r = read_peaks('my_peaks.pks')
      peaks = r['peaks']
      peaks[1].initial_energy = 6.4
      write_peaks('mypeaks.pks', peaks)
   """
   lines = []
   if (background != None):
      lines.append('Background_exponent,'     + str(background.exponent)+'\n')
      lines.append('Background_top_width,'    + str(background.top_width)+'\n')
      lines.append('Background_bottom_width,' + str(background.bottom_width)+'\n')
      lines.append('Background_tangent,'      + str(background.tangent)+'\n')
      lines.append('Background_compress,'     + str(background.compress)+'\n')
   for peak in peaks:
      lines.append(str(peak.initial_energy) + ', ' + 
                   str(peak.energy_flag)    + ', ' + 
                   str(peak.initial_fwhm)   + ', ' + 
                   str(peak.fwhm_flag)      + ', ' + 
                   str(peak.ampl_factor)    + ', ' + 
                   str(peak.label) + '\n')
   fp = open(file, 'w')
   fp.writelines(lines)
   fp.close()

