"""
Support for Multi-Element Detectors (Med).

Author:        Mark Rivers
Created:       Sept. 18, 2002.  Based on earlier IDL code.
Modifications:
    Edited by SAH (June 2019) to suit analysis code. 
    Stripped out parts that are not required:
        - spline function and option to align the called data. 
"""
from cpf import Mca
import numpy as Numeric
#import spline

#########################################################################
class Med(Mca.Mca):
   """
   The MED class is basically a collection of Mca objects.
   
   This class is device-independent.
   
   Its methods generally simply apply the Mca class methods to each Mca object
   in the collection. The Med class itself is most commonly used for reading
   data from disk files. More importantly, this class is the superclass of the
   epicsMed class.
   """

   def __init__(self, n_detectors=16, file=None):
      """
      Initialization code for creating a new Med object.

      Keywords:
         n_detectors:
            The number of detectors (Mca objects) in the Med.

         file:
            The name of a disk file to read into the Med after it is created.
            The number of detectors in the Med will be changed if the number of
            Mca objects in the disk file is different from the Med.
      """
      Mca.Mca.__init__(self)  # Invoke base class initialization
      self.n_detectors = n_detectors
      self.mcas = []
      for i in range(n_detectors):
         self.mcas.append(Mca.Mca())
      if (file != None): self.read_file(file)


   #########################################################################
   def initial_calibration(self, energy):
      """
      Performs an initial energy calibration for each Mca in the Med.

      Inputs:
         energy:
            The energy of the largest peak in the spectrum.
            
      See the documentation for Mca.initial_calibration() for more information.
      """
      # Read the data first in case this is a hardware Mca
      junk = self.get_data()
      for mca in self.mcas:
         mca.initial_calibration(energy)

   #########################################################################
   def final_calibration(self, peaks):
      """
      Performs a final energy calibration for each Mca in the Med.

      Inputs:
         peaks:
            A list of McaPeak objects. This list is typically read from a
            disk file with function Mca.read_peaks().
            
      See the documentation for Mca.final_calibration() for more information.
      """
      # Read the data first in case this is a hardware Mca
      junk = self.get_data()
      for mca in self.mcas:
         mca.final_calibration(peaks)

   #########################################################################
   def get_energy(self):
      """
      Returns a list of energy arrays, one array for each Mca in the Med.
      See the documentation for Mca.get_energy() for more information.
      """
      energy = []
      for mca in self.mcas:
         energy.append(mca.get_energy())
      return energy

   #########################################################################
   def get_mcas(self):
      """
      Returns a list of Mca objects from the Med.
      """
      return self.mcas

   #########################################################################
   def get_calibration(self):
      """
      Returns a list of McaCalibration objects, one for each Mca in the Med.
      """
      calibration = []
      for mca in self.mcas:
         calibration.append(mca.get_calibration())
      self.calibration = calibration
      return calibration

   #########################################################################
   def set_calibration(self, calibration):
      """
      This procedure sets the calibration parameters for the Med.
      The calibration information is contained in an object or list of 
      objects of type McaCalibration.
      
      Inputs:
         calibration:
            A single object or a list of objects of type McaCalibration
            containing the calibration parameters for each Mca.
            If a single object is passed then this is written to each Mca.
            If a list of objects is passed then calibration[i] is written to
            Mca[i].
      """
      if (isinstance(calibration, Mca.McaCalibration)):
         for mca in self.mcas:
            mca.set_calibration(calibration)
      else:  # Assume it is a list or tuple
         for i in range(self.n_detectors):
            self.mcas[i].set_calibration(calibration[i])


   #########################################################################
   def get_elapsed(self):
      """
      Returns the elapsed parameters for the Med.
      The elapsed information is contained in a list of structures of type
      McaElapsed.
      
      Outputs:
         Returns a list of structures of type McaElapsed.
         
      Procedure:
         This function simply invokes Mca.get_elapsed for each Mca in the Med
         and stores the results in the returned list.
      """
      elapsed = []
      for mca in self.mcas:
         elapsed.append(mca.get_elapsed())
      return elapsed

   #########################################################################
   def set_elapsed(self, elapsed):
      """
      Sets the elapsed parameters for the Med.
      The elapsed information is contained in an object or list of 
      objects of type McaElapsed.

      Inputs:
         elapsed:
            A single structure or a list of structures of type McaElapsed
            containing the elapsed parameters for each Mca.
            If a single object is passed then this is written to each Mca.
            If a list of objects is passed then elapsed[i] is written to Mca[i].
      """
      if (isinstance(elapsed, Mca.McaElapsed)):
         for mca in self.mcas:
            mca.set_elapsed(elapsed)
      else:  # Assume it is a list or tuple
         for i in range(self.n_detectors):
            self.mcas[i].set_elapsed(elapsed[i])


   #########################################################################
   def get_presets(self):
      """
      Returns the preset parameters for the Med.
      The preset information is contained in a list of objects of type
      McaPresets.
      
      Outputs:
         Returns a list of structures of type McaPresets.
         
      Procedure:
         This function simply invokes Mca.get_presets() for each Mca in the Med
         and stores the results in the returned list.
      """
      presets = []
      for mca in self.mcas:
         presets.append(mca.get_presets())
      return presets

   #########################################################################
   def set_presets(self, presets):
      """
      This procedure set the preset parameters for the Med.
      The elapsed information is contained in an object or list of 
      objects of type McaPresets.

      Inputs:
         presets:
            A single object or a list of objects of type McaPresets containing
            the preset parameters for each Mca.
            If a single object is passed then this is written to each Mca.
            If a list of objects is passed then presets[i] is written to Mca[i].
      """
      if (isinstance(presets, Mca.McaPresets)):
         for mca in self.mcas:
            mca.set_presets(presets)
      else:  # Assume it is a list or tuple
         for i in range(self.n_detectors):
            self.mcas[i].set_presets(presets[i])

   #########################################################################
   def get_rois(self):
      """
      Returns the region-of-interest information for each Mca in the Med.

      Outputs:
         Returns a list of list of lists of McaRoi objects.
         The length of the outer list is self.n_detectors, the length of the
         list for each Mca is the number of ROIs defined for that Mca.
      """
      rois = []
      for mca in self.mcas:
         rois.append(mca.get_rois())
      return rois

   #########################################################################
   def get_roi_counts(self, background_width=1):
      """
      Returns the net and total counts for each Roi in each Mca in the Med.

      Outputs:
         Returns a tuple (total, net).  total and net are lists of lists
         containing the total and net counts in each ROI.  The length of the
         outer list is self.n_detectors, the length of the total and net lists
         list for each Mca is the number of ROIs defined for that Mca.
      """
      total = []
      net = []
      for mca in self.mcas:
         t, n = mca.get_roi_counts(background_width)
         total.append(t)
         net.append(n)
      return (total, net)

   #########################################################################
   def set_rois(self, rois, energy=0):
      """
      This procedure sets the ROIs for the Med.
      The elapsed information is contained in a list of McaRoi objects,
      or list of such lists.

      Inputs:
         rois:
            A single list or a nested list of objects McaROI objects.
            If a single list is passed then this is written to each Mca.
            If a list of lists is passed then rois[i][*] is written to Mca[i].
      """
      if (len(rois) <= 1):
         for mca in self.mcas:
            mca.set_rois(rois, energy=energy)
      else:
         for i in range(self.n_detectors):
            self.mcas[i].set_rois(rois[i], energy=energy)

   #########################################################################
   def add_roi(self, roi, energy=0):
      """
      This procedure adds an ROI to each Mca in the Med.

      Inputs:
         roi:
            A single McaROI to be added.
      """
      for mca in self.mcas:
         mca.add_roi(roi, energy=energy)
         
   #########################################################################
   def delete_roi(self, index):
      """
      This procedure deletes the ROI at position "index" from each Mca in the
      Med.

      Inputs:
         index:  The index number of the ROI to be deleted.
      """
      for mca in self.mcas:
         mca.delete_roi(index)

   #########################################################################
   def copy_rois(self, source_mca=0, energy=0):
      """
      This procedure copies the ROIs defined for one Mca in the Med to all of
      the other Mcas.

      Inputs:
         source_mca:
            The index number of the Mca from which the ROIs are to
            be copied.  This number ranges from 0 to self.n_detectors-1.
            The default is the first Mca (index=0).
            
      Keywords:
         energy:
            Set this keyword if the ROIs should be copied by their position
            in energy rather than in channels. This is very useful when 
            copying ROIs when the calibration parameters for each Mca in 
            the Med are not identical.
      """
      rois = self.mcas[source_mca].get_rois(energy=energy)
      self.set_rois(rois, energy=energy)

   #########################################################################
   def get_data(self, total=0, align=0):
      """
      Returns the data from each Mca in the Med as a 2-D Numeric array
      
      Keywords:
         total:
            Set this keyword to return the sum of the spectra from all
            of the Mcas as a 1-D Numeric array.
            
         align:
            Set this keyword to return spectra which have been shifted and
            and stretched to match the energy calibration parameters of the
            first detector.  This permits doing arithmetic on a
            "channel-by-channel" basis. This keyword can be used alone
            or together with the TOTAL keyword, in which case the data
            are aligned before summing.
            
      Outputs:
         By default this function returns a long 2-D array of counts dimensioned
         [nchans, self.n_detectors]
         If the "total" keyword is set then the function returns a long 1-D
         array dimensioned [nchans].
      """
      temp = self.mcas[0].get_data()
      nchans = len(temp)
      data = Numeric.zeros((self.n_detectors, nchans))
      for i in range(self.n_detectors):
         data[i,:] = self.mcas[i].get_data()
      if (align != 0):
          print('This option has been removed')
          stop
#         ref_energy = self.mcas[0].get_energy()
#         for i in range(self.n_detectors):
#            energy = self.mcas[i].get_energy()
#            temp = spline.spline_interpolate(energy, data[i,:], ref_energy)
#            data[i,:] = (temp+.5).astype(Numeric.Int)
      if (total != 0):
         d = Numeric.sum(data)
         self.data = d
         return d
      else:
         self.data = data
         return data

   #########################################################################
   def read_file(self, file, netcdf=0):
      """
      Reads a disk file into an Med object. The file contains the information
      from the Med object which it makes sense to store permanently, but does
      not contain all of the internal state information for the Med.

      Inputs:
         file:
            The name of the disk file to read.
      """
      if (netcdf != 0):
         r = Mca.read_netcdf_file(file)
      else:
         r = Mca.read_ascii_file(file)
      self.name = file
      self.n_detectors = r['n_detectors']
      self.mcas = []
      for i in range(self.n_detectors):
         self.mcas.append(Mca.Mca())
         self.mcas[i].set_rois(r['rois'][i])
         self.mcas[i].set_data(r['data'][i])
         self.mcas[i].set_name(self.name + ':' + str(i+1))
      self.set_elapsed(r['elapsed'])
      self.set_calibration(r['calibration'])
      self.set_environment(r['environment'])
