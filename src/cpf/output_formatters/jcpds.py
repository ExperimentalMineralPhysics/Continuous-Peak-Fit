# -*- coding: utf-8 -*-
"""
Support for calculating D spacing for powder diffraction lines as
as function of pressure and temperature, given symmetry, zero-pressure lattice
constants and equation of state parameters.

Author:
  Mark Rivers

Created:
   Sept. 10, 2002 from older IDL version

Modifications:
    Sept. 26, 2002 MLR
        - Implemented Birch-Murnaghan solver using CARSnp.newton root finder
    Mai 27, 2014 Clemens Prescher
        - changed np function to numpy versions,
        - using scipy optimize for solving the inverse Birch-Murnaghan problem
        - fixed a bug which was causing a gamma0 to be 0 for cubic unit cell
    August 22, 2014 Clemens Prescher
        - calculation of d spacings is now done by using arrays
        - added several new utility function -- calculate_d0, add_reflection
        - updated the write_file function to be able to use new standard
    August 26, 2014 Clemens Prescher
        - added sorting functions
        - fixed the d spacing calculation for triclinic structure - equation used was wrong...
    August 27, 2014 Clemens Prescher
        - added modified flag and the surrounding functions. When an attribute is changed, it will set it to true and the
          filename and name will have an asterisk appended to indicate that this is not the original jcpds loaded
        - added a reload function
        - renamed read and write to load and save
        - the load function will now reset all parameters (previously parameters not set in the newly loaded file, were
          taken over from the previous state of the object)
    September 2025 Simon Hunt
        - added "all" option to delete_reflection
        - added functions to fit unitcell to observations

"""
import logging
logger = logging.getLogger(__name__)

import os
import string
import numpy as np
# from scipy.optimize import minimize
import lmfit
from uncertainties import ufloat, ufloat_fromstr
from uncertainties.umath import sin, cos, tan, sqrt

# import time

class jcpds_reflection:
    """
    Class that defines a reflection.
    Attributes:
       d0:     Zero-pressure lattice spacing
       d:      Lattice spacing at P and T
       dobs:   Observed lattice spacing (possiblt at P and T)
       inten:  Relative intensity to most intense reflection for this material
       h:      H index for this reflection
       k:      K index for this reflection
       l:      L index for this reflection

    """

    def __init__(self, h=0., k=0., l=0., intensity=0., d=0., dobs=None):
        self.d0 = d
        self.d = d
        self.dobs = dobs
        self.intensity = intensity
        self.h = h
        self.k = k
        self.l = l


class jcpds(object):
    def __init__(self):
        self._filename = ' '
        self._name = ' '
        self.version = 0
        self.comments = []
        self.symmetry = ''
        self.k0 = ufloat(0.,1E-10)
        self.k0p0 = ufloat(0.,1E-10)  # k0p at 298K
        self.k0p = ufloat(0.,1E-10)  # k0p at high T
        self.dk0dt = ufloat(0.,1E-10)
        self.dk0pdt = ufloat(0.,1E-10)
        self.alpha_t0 = ufloat(0.,1E-10)  # alphat at 298K
        self.alpha_t = ufloat(0.,1E-10)  # alphat at high temp.
        self.d_alpha_dt = ufloat(0.,1E-10)
        self.a0 = 0.
        self.b0 = 0.
        self.c0 = 0.
        self.alpha0 = 0.
        self.beta0 = 0.
        self.gamma0 = 0.
        self.v0 = 0.
        self.a = ufloat(0.,1E-10)
        self.b = ufloat(0.,1E-10)
        self.c = ufloat(0.,1E-10)
        self.alpha = ufloat(0.,1E-10)
        self.beta = ufloat(0.,1E-10)
        self.gamma = ufloat(0.,1E-10)
        self.v = ufloat(0.,1E-10)
        self.pressure = ufloat(0.,1E-10)
        self.temperature = 298.
        self.reflections = []
        self.modified = False

    def load_file(self, filename):
        """
        Reads a JCPDS file into the JCPDS object.

        Inputs:
           file:  The name of the file to read.

        Procedure:
           This procedure read the JCPDS file.  There are several versions of the
           formats used for JCPDS files.  Versions 1, 2 and 3 used a fixed
           format, where a particular entry had to be in a specific location on
           a specific line.  Versions 2 and 3 were used only by Dan Shim.
           This routine can read these old files, but no new files should be
           created in this format, they should be converted to Version 4.
           Version 4 is a "keyword" driven format.  Each line in the file is of
           the form:
           KEYWORD: value
           The order of the lines is not important, except that the first line of
           the file must be "VERSION: 4".
           The following keywords are currently supported:
               COMMENT:    Any information describing the material, literature
                           references, etc.  There can be multiple comment lines
                           per file.
               K0:         The bulk modulus in GPa.
               K0P:        The change in K0 with pressure, for Birch-Murnaghan
                           equation of state.  Dimensionless.
               DK0DT:      The temperature derivative of K0, GPa/K.
               DK0PDT:     The temperature derivative of K0P, 1/K.
               SYMMETRY:   One of CUBIC, TETRAGONAL, HEXAGONAL, RHOMBOHEDRAL,
                           ORTHORHOMBIC, MONOCLINIC or TRICLINIC
               A:          The unit cell dimension A
               B:          The unit cell dimension B
               C:          The unit cell dimension C
               ALPHA:      The unit cell angle ALPHA
               BETA:       The unit cell angle BETA
               GAMMA:      The unit cell angle GAMMA
               VOLUME:     The unit cell volume
               ALPHAT:     The thermal expansion coefficient, 1/K
               DALPHADT:   The temperature derivative of the thermal expansion
                           coefficient, 1/K^2
               DIHKL:      For each reflection, the D spacing in Angstrom, the
                           relative intensity (0-100), and the H, K, L indices.

           This procedure calculates the D spacing of each relfection, using the
           symmetry and unit cell parameters from the file.  It compares the
           calculated D spacing with the input D spacing for each line.  If they
           disagree by more than 0.1% then a warning message is printed.
           The following is an example JCPDS file in the Version 4 format:
               VERSION:  4
               COMMENT: Alumina (JCPDS 0-173, EOS n/a)
               K0:          194.000
               K0P:           5.000
               SYMMETRY: HEXAGONAL
               A:            4.758
               C:            12.99
               VOLUME:        22.0640
               ALPHAT:    2.000e-6
               DIHKL:        3.4790      75.0   0   1   2
               DIHKL:        2.5520      90.0   1   0   4
               DIHKL:        2.3790      40.0   1   1   0
               DIHKL:        2.0850     100.0   1   1   3
               DIHKL:        1.7400      45.0   0   2   4
               DIHKL:        1.6010      80.0   1   1   6
               DIHKL:        1.4040      30.0   2   1   4
               DIHKL:        1.3740      50.0   3   0   0
               DIHKL:        1.2390      16.0   1   0  10

           Note that B and ALPHA, BETA and GAMMA are not present, since they are
           not needed for a hexagonal material, and will be simple ignorred if
           they are present.
        """
        self.__init__()
        # Initialize variables
        self._filename = filename
        # Construct base name = file without path and without extension
        name = os.path.basename(filename)
        pos = name.find('.')
        if (pos >= 0): name = name[0:pos]
        self._name = name
        line = ''
        version = 0.
        self.comments = []
        nd = 0
        self.reflections = []

        # Determine what version JCPDS file this is
        # In current files have the first line starts with the string VERSION:
        fp = open(filename, 'r')
        line = fp.readline()
        pos = line.index(' ')
        tag = line[0:pos].upper()
        value = line[pos:].strip()
        if (tag == 'VERSION:'):
            self.version = value
            # This is the current, keyword based version of JCPDS file
            while (1):
                line = fp.readline()
                if (line == ''): break
                pos = line.index(' ')
                tag = line[0:pos].upper()
                value = line[pos:].strip()
                if (tag == 'COMMENT:'):
                    self.comments.append(value)
                elif (tag == 'K0:'):
                    self.k0 = ufloat_fromstr(value)
                elif (tag == 'K0P:'):
                    self.k0p0 = ufloat_fromstr(value)
                elif (tag == 'DK0DT:'):
                    self.dk0dt = ufloat_fromstr(value)
                elif (tag == 'DK0PDT:'):
                    self.dk0pdt = ufloat_fromstr(value)
                elif (tag == 'SYMMETRY:'):
                    self.symmetry = value.upper()
                elif (tag == 'A:'):
                    self.a0 = float(value)
                elif (tag == 'B:'):
                    self.b0 = float(value)
                elif (tag == 'C:'):
                    self.c0 = float(value)
                elif (tag == 'ALPHA:'):
                    self.alpha0 = float(value)
                elif (tag == 'BETA:'):
                    self.beta0 = float(value)
                elif (tag == 'GAMMA:'):
                    self.gamma0 = float(value)
                elif (tag == 'VOLUME:'):
                    self.v0 = float(value)
                elif (tag == 'ALPHAT:'):
                    self.alpha_t0 = ufloat_fromstr(value)
                elif (tag == 'DALPHADT:'):
                    self.d_alpha_dt = ufloat_fromstr(value)
                elif (tag == 'DIHKL:'):
                    dtemp = value.split()
                    dtemp = list(map(float, dtemp))
                    reflection = jcpds_reflection()
                    reflection.d0 = dtemp[0]
                    reflection.intensity = dtemp[1]
                    reflection.h = int(dtemp[2])
                    reflection.k = int(dtemp[3])
                    reflection.l = int(dtemp[4])
                    self.reflections.append(reflection)
        else:
            # This is an old format JCPDS file
            self.version = 1.
            header = ''
            self.comments.append(line)  # Read above
            line = fp.readline()
            # Replace any commas with blanks, split at blanks
            temp = string.split(line.replace(',', ' '))
            temp = list(map(float, temp[0:5]))
            # The symmetry codes are as follows:
            #   1 -- cubic
            #   2 -- hexagonal
            if (temp[0] == 1):
                self.symmetry = 'CUBIC'
            elif (temp[0] == 2):
                self.symmetry = 'HEXAGONAL'
            self.a0 = temp[1]
            self.k0 = temp[2]
            self.k0p0 = temp[3]
            c0a0 = temp[4]
            self.c0 = self.a0 * c0a0
            line = fp.readline()  # Ignore, just column labels

            while (1):
                line = fp.readline()
                if (line == ''): break
                dtemp = line.split()
                dtemp = list(map(float, dtemp))
                reflection = jcpds_reflection()
                reflection.d0 = dtemp[0]
                reflection.intensity = dtemp[1]
                reflection.h = int(dtemp[2])
                reflection.k = int(dtemp[3])
                reflection.l = int(dtemp[4])
                self.reflections.append(reflection)

        fp.close()
        self.compute_v0()
        self.a = self.a0
        self.b = self.b0
        self.c = self.c0
        self.alpha = self.alpha0
        self.beta = self.beta0
        self.gamma = self.gamma0
        self.v = self.v0
        # Compute D spacings, make sure they are consistent with the input values

        self.compute_d()
        for reflection in self.reflections:
            reflection.d0 = reflection.d

        self.modified = False

        ## we just removed this check because it should be better to care more about the actual a,b,c values than
        # individual d spacings
        # reflections = self.get_reflections()
        # for r in reflections:
        #     diff = abs(r.d0 - r.d) / r.d0
        #     if (diff > .001):
        #         logger.info(('Reflection ', r.h, r.k, r.l, \
        #             ': calculated D ', r.d, \
        #             ') differs by more than 0.1% from input D (', r.d0, ')'))


    def save_file(self, filename):
        """
        Writes a JCPDS object to a file.

        Inputs:
           filename:  The name of the file to written.

        Procedure:
           This procedure writes a JCPDS file.  It always writes files in the
           current, keyword-driven format (Version 4).  See the documentation for
           read_file() for information on the file format.

        Example:
           This reads an old format file, writes a new format file.
           j = jcpds.jcpds()
           j.read_file('alumina_old.jcpds')
           j.write_file('alumina_new.jcpds')
        """
        fp = open(filename, 'w')
        fp.write('VERSION:   4\n')
        for comment in self.comments:
            fp.write('COMMENT: ' + comment+'\n')
        fp.write('K0:       ' + str(self.k0)+'\n')
        fp.write('K0P:      ' + str(self.k0p0)+'\n')
        fp.write('DK0DT:    ' + str(self.dk0dt)+'\n')
        fp.write('DK0PDT:   ' + str(self.dk0pdt)+'\n')
        fp.write('SYMMETRY: ' + self.symmetry+'\n')
        fp.write('A:        ' + str(self.a0)+'\n')
        fp.write('B:        ' + str(self.b0)+'\n')
        fp.write('C:        ' + str(self.c0)+'\n')
        fp.write('ALPHA:    ' + str(self.alpha0)+'\n')
        fp.write('BETA:     ' + str(self.beta0)+'\n')
        fp.write('GAMMA:    ' + str(self.gamma0)+'\n')
        fp.write('VOLUME:   ' + str(self.v0)+'\n')
        fp.write('ALPHAT:   ' + str(self.alpha_t0)+'\n')
        fp.write('DALPHADT: ' + str(self.d_alpha_dt)+'\n')
        reflections = self.get_reflections()
        for r in reflections:
            fp.write('DIHKL:    {:g}\t{:g}\t{:g}\t{:g}\t{:g}\n'.format(r.d0, r.intensity, r.h, r.k, r.l))
        fp.close()

        self._filename = filename
        name = os.path.basename(filename)
        pos = name.find('.')
        if (pos >= 0): name = name[0:pos]
        self._name = name

        self.modified = False

    def reload_file(self):
        pressure = self.pressure
        temperature = self.temperature
        self.load_file(self._filename)
        self.pressure = pressure
        self.temperature = temperature
        self.compute_d()

    def __setattr__(self, key, value):
        if key in ['comments', 'a0', 'b0', 'c0', 'alpha0', 'beta0', 'gamma0',
                   'symmetry', 'k0', 'k0p0', 'dk0dt', 'dk0pdt',
                    'alpha_t0', 'd_alpha_dt', 'reflections']:
            self.modified = True
        super(jcpds, self).__setattr__(key, value)

    @property
    def filename(self):
        if self.modified:
            return self._filename + '*'
        else:
            return self._filename

    @filename.setter
    def filename(self, value):
        self._filename = value

    @property
    def name(self):
        if self.modified:
            return self._name + '*'
        else:
            return self._name

    @name.setter
    def name(self, value):
        self._name = value


    def compute_v0(self):
        """
        Computes the unit cell volume of the material at zero pressure and
        temperature from the unit cell parameters.

        Procedure:
           This procedure computes the unit cell volume from the unit cell
           parameters.

        Example:
           Compute the zero pressure and temperature unit cell volume of alumina
           j = jcpds()
           j.read_file('alumina.jcpds')
           j.compute_v0()
        """
        if (self.symmetry == 'CUBIC'):
            self.b0 = self.a0
            self.c0 = self.a0
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 90.

        elif (self.symmetry == 'TETRAGONAL'):
            self.b0 = self.a0
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 90.

        elif (self.symmetry == 'ORTHORHOMBIC'):
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 90.

        elif (self.symmetry == 'HEXAGONAL'):
            self.b0 = self.a0
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 120.

        elif (self.symmetry == 'RHOMBOHEDRAL'):
            self.b0 = self.a0
            self.c0 = self.a0
            self.beta0 = self.alpha0
            self.gamma0 = self.alpha0

        elif (self.symmetry == 'MONOCLINIC'):
            self.alpha0 = 90.
            self.gamma0 = 90.

        elif (self.symmetry == 'TRICLINIC'):
            pass

        dtor = np.pi / 180.
        self.v0 = (self.a0 * self.b0 * self.c0 *
                   np.sqrt(1. -
                           np.cos(self.alpha0 * dtor) ** 2 -
                           np.cos(self.beta0 * dtor) ** 2 -
                           np.cos(self.gamma0 * dtor) ** 2 +
                           2. * (np.cos(self.alpha0 * dtor) *
                                 np.cos(self.beta0 * dtor) *
                                 np.cos(self.gamma0 * dtor))))

    def compute_volume(self, pressure=None, temperature=None):
        """
        Computes the unit cell volume of the material.
        It can compute volumes at different pressures and temperatures.

        Keywords:
           pressure:
              The pressure in GPa.  If not present then the pressure is
              assumed to be 0.

           temperature:
              The temperature in K.  If not present or zero, then the
              temperature is assumed to be 298K, i.e. room temperature.

        Procedure:
           This procedure computes the unit cell volume.  It starts with the
           volume read from the JCPDS file or computed from the zero-pressure,
           room temperature lattice constants.  It does the following:
              1) Corrects K0 for temperature if DK0DT is non-zero.
              2) Computes volume at zero-pressure and the specified temperature
                 if ALPHAT0 is non-zero.
              3) Computes the volume at the specified pressure if K0 is non-zero.
                 The routine uses the IDL function FX_ROOT to solve the third
                 order Birch-Murnaghan equation of state.

        Example:
           Compute the unit cell volume of alumina at 100 GPa and 2500 K.
           j = jcpds()
           j.read_file('alumina.jcpds')
           j.compute_volume(100, 2500)

        """
        if pressure == None:
            pressure = self.pressure
        else:
            self.pressure = pressure

        if temperature == None:
            temperature = self.temperature
        else:
            self.temperature = temperature

        # Assume 0 K really means room T
        if (temperature == 0): temperature = 298.
        # Compute values of K0, K0P and alphat at this temperature
        self.alpha_t = self.alpha_t0 + self.d_alpha_dt * (temperature - 298.)
        self.k0p = self.k0p0 + self.dk0pdt * (temperature - 298.)

        if (pressure == 0.):
            self.v = self.v0 * (1 + self.alpha_t * (temperature - 298.))
        else:
            if (self.k0 <= 0.):
                logger.info('K0 is zero, computing zero pressure volume')
                self.v = self.v0
            else:
                self.mod_pressure = pressure - \
                                    self.alpha_t * self.k0 * (temperature - 298.)
                
                
                
                def tempfunc(self, var):
                    return self.bm3_inverse(var)
                
                fit_params = lmfit.create_params(P=1)
                out = lmfit.minimize(self.bm3_inverse, fit_params)
                res = ufloat(out.params["P"].value, out.params["P"].stderr)
                # print
                # res = minimize(self.bm3_inverse, 1.)
                self.v = self.v0 / res
                
                # res = minimize(self.bm3_inverse, 1.)
                # self.v = self.v0 / float(res.x)


    def compute_unitcell_volume(self):
        """
        Computes the unit cell volume of the material from given unit cell parameters.
    
        Procedure:
           This procedure computes the unit cell volume from the unit cell
           parameters.
    
        Example:
           Compute the zero pressure and temperature unit cell volume of alumina
           j = jcpds()
           j.read_file('alumina.jcpds')
           j.compute_unitcell_volume()
        """
        self.apply_symmetry()
        dtor = np.pi / 180.
        self.v = (self.a * self.b * self.c *
                   np.sqrt(1. -
                           np.cos(self.alpha * dtor) ** 2 -
                           np.cos(self.beta * dtor) ** 2 -
                           np.cos(self.gamma * dtor) ** 2 +
                           2. * (np.cos(self.alpha * dtor) *
                                 np.cos(self.beta * dtor) *
                                 np.cos(self.gamma * dtor))))


    def apply_symmetry(self):
        """
        Updates the unit cell parameters to correspond to the unit cell symmetry
        """
        
        if (self.symmetry == 'CUBIC'):
            self.b0 = self.a0
            self.c0 = self.a0
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 90.
            self.b = self.a
            self.c = self.a
            self.alpha = 90.
            self.beta = 90.
            self.gamma = 90.

        elif (self.symmetry == 'TETRAGONAL'):
            self.b0 = self.a0
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 90.
            self.b = self.a
            self.alpha = 90.
            self.beta = 90.
            self.gamma = 90.

        elif (self.symmetry == 'ORTHORHOMBIC'):
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 90.
            self.alpha = 90.
            self.beta = 90.
            self.gamma = 90.

        elif (self.symmetry == 'HEXAGONAL'):
            self.b0 = self.a0
            self.alpha0 = 90.
            self.beta0 = 90.
            self.gamma0 = 120.
            self.b = self.a
            self.alpha = 90.
            self.beta = 90.
            self.gamma = 120.

        elif (self.symmetry == 'RHOMBOHEDRAL'):
            self.b0 = self.a0
            self.c0 = self.a0
            self.beta0 = self.alpha0
            self.gamma0 = self.alpha0
            self.b = self.a
            self.c = self.a
            self.beta = self.alpha
            self.gamma = self.alpha

        elif (self.symmetry == 'MONOCLINIC'):
            self.alpha0 = 90.
            self.gamma0 = 90.
            self.alpha = 90.
            self.gamma = 90.

        elif (self.symmetry == 'TRICLINIC'):
            pass
        
    
    def get_unique_unitcell_params(self):
        """
        Get indiepende values for unit cell given symmetry

        """
        if (self.symmetry == 'CUBIC'):
            return ["a"]

        elif (self.symmetry == 'TETRAGONAL'):
            return ["a", "c"]

        elif (self.symmetry == 'ORTHORHOMBIC'):
            return ["a", "b", "c"]

        elif (self.symmetry == 'HEXAGONAL'):
            return ["a", "c"]

        elif (self.symmetry == 'RHOMBOHEDRAL'):
            return ["a", "alpha"]
        
        elif (self.symmetry == 'MONOCLINIC'):
            return ["a", "b", "c", "beta"]

        elif (self.symmetry == 'TRICLINIC'):
            return ["a", "b", "c", "alpha", "beta", "gamma"]
            
        

    def fit_lattice_parameters(self, weighted=True):
        """
        Fits the unit cell volume of the material to the reflection d values.

        Procedure:
           This procedure computes the unit cell volume from the observed d spacing 
           in the reflections. 
           
        Example:
           Compute the unit cell volume of from observations.
           
           jcpds_obj = jcpds()
           jcpds_obj.load_file(jcpds_file)
           jcpds_obj.remove_reflection("all")
           
           jcpds_obj.add_reflection(h,k,l,dobs=obshkl)
           
           jcpds_obj.compute_d0() # compute lattice parameters for unit cell from jcpds
           jcpds_obj.fit_lattice_parameters
          
            
        Returns
        -------
        out : lmfit object
            fit object from lmfit.
 
          
        """
        use = self.get_unique_unitcell_params()
        params = lmfit.Parameters()
        for i in range(len(use)):
            val = getattr(self,use[i])
            try:
                val=getattr(self,use[i]).nominal_value
            except:
                val=val
            params.add(use[i], val, vary=True, min=0)
            
        #get weights from observation errors
        weights = []
        obs = np.array([])
        for i in self.get_reflections():
            try:
                weights.append(1/i.dobs.std_dev**2)
            except:
                weights.append(np.nan)
            try:
                obs = np.append(obs, i.dobs.nominal_value)
            except:
                obs = np.append(obs, i.dobs)
            
        # out = lmfit.minimize(self._resid, params)
        cmodel = lmfit.Model(self._lattice_params_model, independent_vars=["jcpds"])
        
        if weighted is True:
            out = cmodel.fit(obs, params, jcpds=None, weights=weights)   
        else:
            out = cmodel.fit(obs, params, jcpds=None)   
            
        uc_parts = self.get_unique_unitcell_params()
        uc_parms = {}
        for ind in range(len(uc_parts)):
            setattr(self, uc_parts[ind], ufloat(out.params[uc_parts[ind]].value, out.params[uc_parts[ind]].stderr) )
            # self.[uc_parts[ind]] = getattr(jcpds_obj, uc_parts[ind])    
        # uc_parms["volume"] = jcpds_obj.v        
        
        #update unitcell volume
        self.compute_unitcell_volume()
        
        return out
        
    
    def _lattice_params_model(self, a=None, b=None, c=None, alpha=None, beta=None, gamma=None, jcpds=None):
         """
         Rerust differences between calculated and observed reflcations
    
         Parameters
         ----------
         params : lmfit.Parameters()
             unitcell parameter set for calculating d-spacings.
    
         Returns
         -------
         list
             Difference in d-spacing between calculated and observed reflections.
    
         """
         use = self.get_unique_unitcell_params()
         lattice_parameters ={}
         for i in use:
             if i in use:
                 lattice_parameters[i] = locals()[i]
             else: 
                 try:
                     lattice_parameters[i] = getattr(self, i).nominal_value
                 except:
                     lattice_parameters[i] = getattr(self, i)
         self.compute_d(lattice_parameters=lattice_parameters)  
         r = self.get_reflections()
         calc = []
         obs = []
         for i in range(len(r)):
             if r[i].dobs != None:
                 obs.append(r[i].dobs)
                 calc.append(r[i].d)
         return np.array(calc)


    # def _resid(self, params, weights=None):
    #     """
    #     Rerust differences between calculated and observed reflcations

    #     Parameters
    #     ----------
    #     params : lmfit.Parameters()
    #         unitcell parameter set for calculating d-spacings.

    #     Returns
    #     -------
    #     list
    #         Difference in d-spacing between calculated and observed reflections.

    #     """
    #     params = params.valuesdict()
    #     self.compute_d(lattice_parameters=params)  
    #     r = self.get_reflections()
    #     calc = []
    #     obs = []
    #     for i in range(len(r)):
    #         if r[i].dobs != None:
    #             obs.append(r[i].dobs)
    #             calc.append(r[i].d)
    #     out = np.array(calc) - np.array(obs)
    #     try:
    #         out = [i.nominal_value for i in out]
    #     except:
    #         out = out
            
    #     if weights is not None:
    #         return out/weights
    #     return out

    def bm3_inverse(self, v0_v_in):
        """
        Returns the value of the third order Birch-Murnaghan equation minus
        pressure.  It is used to solve for V0/V for a given
           P, K0 and K0'.

        Inputs:
           v0_v:   The ratio of the zero pressure volume to the high pressure
                   volume
        Outputs:
           This function returns the value of the third order Birch-Murnaghan
           equation minus pressure.  \

        Procedure:
           This procedure simply computes the pressure using V0/V, K0 and K0',
           and then subtracts the input pressure.

        Example:
           Compute the difference of the calculated pressure and 100 GPa for
           V0/V=1.3 for alumina
           jcpds = obj_new('JCPDS')
           jcpds->read_file,  'alumina.jcpds'
           common bm3_common mod_pressure, k0, k0p
           mod_pressure=100
           k0 = 100
           k0p = 4.
           diff = jcpds_bm3_inverse(1.3)
        """
        if isinstance(v0_v_in, type(lmfit.create_params(P=1))):
            v0_v = v0_v_in["P"].value
        else:
            v0_v = v0_v_in

        out = (1.5 * self.k0 * (v0_v ** (7. / 3.) - v0_v ** (5. / 3.)) *
                (1 + 0.75 * (self.k0p - 4.) * (v0_v ** (2. / 3.) - 1.0)) -
                self.mod_pressure) ** 2
        
        # the output has to be a single number for the minimisers. 
        # so try to strip it down. 
        try:
            out = out.nominal_value
        except: 
            out = out[0].nominal_value
        finally:
            out=out
        return out

    def compute_d0(self):
        """
        computes d0 values for the based on the the current lattice parameters
        """
        a = self.a0
        b = self.b0
        c = self.c0
        degree_to_radians = np.pi / 180.
        alpha = self.alpha0 * degree_to_radians
        beta = self.beta0 * degree_to_radians
        gamma = self.gamma0 * degree_to_radians

        h = np.zeros(len(self.reflections))
        k = np.zeros(len(self.reflections))
        l = np.zeros(len(self.reflections))

        for ind, reflection in enumerate(self.reflections):
            h[ind] = reflection.h
            k[ind] = reflection.k
            l[ind] = reflection.l

        if (self.symmetry == 'CUBIC'):
            d2inv = (h ** 2 + k ** 2 + l ** 2) / a ** 2
        elif (self.symmetry == 'TETRAGONAL'):
            d2inv = (h ** 2 + k ** 2) / a ** 2 + l ** 2 / c ** 2
        elif (self.symmetry == 'ORTHORHOMBIC'):
            d2inv = h ** 2 / a ** 2 + k ** 2 / b ** 2 + l ** 2 / c ** 2
        elif (self.symmetry == 'HEXAGONAL'):
            d2inv = (h ** 2 + h * k + k ** 2) * 4. / 3. / a ** 2 + l ** 2 / c ** 2
        elif (self.symmetry == 'RHOMBOHEDRAL'):
            d2inv = (((1. + np.cos(alpha)) * ((h ** 2 + k ** 2 + l ** 2) -
                                              (1 - np.tan(0.5 * alpha) ** 2) * (h * k + k * l + l * h))) /
                     (a ** 2 * (1 + np.cos(alpha) - 2 * np.cos(alpha) ** 2)))
        elif (self.symmetry == 'MONOCLINIC'):
            d2inv = (h ** 2 / np.sin(beta) ** 2 / a ** 2 +
                     k ** 2 / b ** 2 +
                     l ** 2 / np.sin(beta) ** 2 / c ** 2 +
                     2 * h * l * np.cos(beta) / (a * c * np.sin(beta) ** 2))
        elif (self.symmetry == 'TRICLINIC'):
            V = (a  * b  * c  *
                 np.sqrt(1. - np.cos(alpha) ** 2 - np.cos(beta) ** 2 -
                  np.cos(gamma) ** 2 +
                  2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma)))
            s11 = b ** 2 * c ** 2 * np.sin(alpha) ** 2
            s22 = a ** 2 * c ** 2 * np.sin(beta) ** 2
            s33 = a ** 2 * b ** 2 * np.sin(gamma) ** 2
            s12 = a * b * c ** 2 * (np.cos(alpha) * np.cos(beta) -
                                    np.cos(gamma))
            s23 = a ** 2 * b * c * (np.cos(beta) * np.cos(gamma) -
                                    np.cos(alpha))
            s31 = a * b ** 2 * c * (np.cos(gamma) * np.cos(alpha) -
                                    np.cos(beta))
            d2inv = (s11 * h ** 2 + s22 * k ** 2 + s33 * l ** 2 +
                     2. * s12 * h * k + 2. * s23 * k * l + 2. * s31 * l * h) / V ** 2
        else:
            logger.error(('Unknown crystal symmetry = ' + self.symmetry))
        d_spacings = np.sqrt(1. / d2inv)

        for ind in range(len(self.reflections)):
            self.reflections[ind].d0 = d_spacings[ind]



    def compute_d(self, lattice_parameters=None, pressure=None, temperature=None):
        """
        Computes the D spacings of the material.
        It can compute D spacings at different pressures and temperatures.

        Keywords:
           lattice_parameters:
               dictionary of lattice parameters for the unit cell symmetry.
                              
           pressure:
              The pressure in GPa.  If not present then the pressure is
              assumed to be 0.

           temperature:
              The temperature in K.  If not present or zero, then the
              temperature is assumed to be 298K, i.e. room temperature.

        Outputs:
           None.  The D spacing information in the JCPDS object is calculated.

        Procedure:
            This procedure first calls jcpds.compute_volume().
            It then assumes that each lattice dimension fractionally changes by
            the cube root of the fractional change in the volume.
            Using the equations for the each symmetry class it then computes the
            change in D spacing of each reflection.

        Example:
           Compute the D spacings of alumina at 100 GPa and 2500 K.
           j=jcpds()
           j.read_file('alumina.jcpds')
           j.compute_d(100, 2500)
           refl = j.get_reflections()
           for r in refl:
              # Print out the D spacings at ambient conditions
              print, r.d0
              # Print out the D spacings at high pressure and temperature
              print, r.d
        """
        
        if lattice_parameters == None:
            self.compute_volume(pressure, temperature)
    
            # Assume each cell dimension changes by the same fractional amount = cube
            # root of volume change ratio
            ratio = ((self.v / self.v0) ** (1.0 / 3.0))
            self.a = self.a0 * ratio
            self.b = self.b0 * ratio
            self.c = self.c0 * ratio
    
        else:
            # update unit cell
            if "a" in list(lattice_parameters):
                self.a = lattice_parameters["a"]
            if "b" in list(lattice_parameters):
                self.b = lattice_parameters["b"]
            if "c" in list(lattice_parameters):
                self.c = lattice_parameters["c"]
            if "alpha" in list(lattice_parameters):
                self.alpha = lattice_parameters["alpha"]
            if "beta" in list(lattice_parameters):
                self.beta = lattice_parameters["beta"]
            if "gamma" in list(lattice_parameters):
                self.gamma = lattice_parameters["gamma"]
            self.apply_symmetry()
            
        a = self.a
        b = self.b
        c = self.c
        dtor = np.pi / 180.
        alpha = self.alpha0 * dtor # in radians
        beta = self.beta0 * dtor   # in radians
        gamma = self.gamma0 * dtor # in radians

        h = np.zeros(len(self.reflections))
        k = np.zeros(len(self.reflections))
        l = np.zeros(len(self.reflections))

        for ind, reflection in enumerate(self.reflections):
            h[ind] = reflection.h
            k[ind] = reflection.k
            l[ind] = reflection.l

        if (self.symmetry == 'CUBIC'):
            d2inv = (h ** 2 + k ** 2 + l ** 2) / a ** 2
        elif (self.symmetry == 'TETRAGONAL'):
            d2inv = (h ** 2 + k ** 2) / a ** 2 + l ** 2 / c ** 2
        elif (self.symmetry == 'ORTHORHOMBIC'):
            d2inv = h ** 2 / a ** 2 + k ** 2 / b ** 2 + l ** 2 / c ** 2
        elif (self.symmetry == 'HEXAGONAL'):
            d2inv = (h ** 2 + h * k + k ** 2) * 4. / 3. / a ** 2 + l ** 2 / c ** 2
        elif (self.symmetry == 'RHOMBOHEDRAL'):
            d2inv = (((1. + cos(alpha)) * ((h ** 2 + k ** 2 + l ** 2) -
                                              (1 - tan(0.5 * alpha) ** 2) * (h * k + k * l + l * h))) /
                     (a ** 2 * (1 + cos(alpha) - 2 * cos(alpha) ** 2)))
        elif (self.symmetry == 'MONOCLINIC'):
            d2inv = (h ** 2 / sin(beta) ** 2 / a ** 2 +
                     k ** 2 / b ** 2 +
                     l ** 2 / sin(beta) ** 2 / c ** 2 +
                     2 * h * l * cos(beta) / (a * c * sin(beta) ** 2))
        elif (self.symmetry == 'TRICLINIC'):
            V = (  a * b  * c  *
                 sqrt(1. - cos(alpha) ** 2 - cos(beta) ** 2 -
                  cos(gamma) ** 2 +
                  2 * cos(alpha) * cos(beta) * cos(gamma)))
            s11 = b ** 2 * c ** 2 * sin(alpha) ** 2
            s22 = a ** 2 * c ** 2 * sin(beta) ** 2
            s33 = a ** 2 * b ** 2 * sin(gamma) ** 2
            s12 = a * b * c ** 2 * (cos(alpha) * cos(beta) -
                                    cos(gamma))
            s23 = a ** 2 * b * c * (cos(beta) * cos(gamma) -
                                    cos(alpha))
            s31 = a * b ** 2 * c * (cos(gamma) * cos(alpha) -
                                    cos(beta))
            d2inv = (s11 * h ** 2 + s22 * k ** 2 + s33 * l ** 2 +
                     2. * s12 * h * k + 2. * s23 * k * l + 2. * s31 * l * h) / V ** 2
        else:
            logger.error(('Unknown crystal symmetry = ' + self.symmetry))
        d_spacings = (1. / d2inv)**(1/2)
        for ind in range(len(self.reflections)):
            self.reflections[ind].d = d_spacings[ind]

    def add_reflection(self, h=0., k=0., l=0., intensity=0., d=0., dobs=None, dobs_err=None):
        if dobs_err is not None:
            dobs = ufloat(dobs, dobs_err)
        new_reflection = jcpds_reflection(h,k,l, intensity, d, dobs=dobs)
        self.reflections.append(new_reflection)
        self.modified = True

    def remove_reflection(self, ind):
        if ind == "all":
            self.reflections = []
        else:
            del self.reflections[ind]
        self.modified = True

    def get_reflections(self):
        """
        Returns the information for each reflection for the material.
        This information is an array of elements of class jcpds_reflection
        """
        return self.reflections

    def reorder_reflections_by_index(self, ind_list, reversed=False):
        if reversed:
            ind_list= ind_list[::-1]
        new_reflections = []
        for ind in ind_list:
            new_reflections.append(self.reflections[ind])

        modified_flag = self.modified
        self.reflections = new_reflections
        self.modified = modified_flag

    def sort_reflections_by_h(self, reversed = False):
        h_list = []
        for reflection in self.reflections:
            h_list.append(reflection.h)
        sorted_ind = np.argsort(h_list)
        self.reorder_reflections_by_index(sorted_ind, reversed)

    def sort_reflections_by_k(self, reversed = False):
        k_list = []
        for reflection in self.reflections:
            k_list.append(reflection.k)
        sorted_ind = np.argsort(k_list)
        self.reorder_reflections_by_index(sorted_ind, reversed)

    def sort_reflections_by_l(self, reversed = False):
        l_list = []
        for reflection in self.reflections:
            l_list.append(reflection.l)
        sorted_ind = np.argsort(l_list)
        self.reorder_reflections_by_index(sorted_ind, reversed)

    def sort_reflections_by_intensity(self, reversed = False):
        intensity_list = []
        for reflection in self.reflections:
            intensity_list.append(reflection.intensity)
        sorted_ind = np.argsort(intensity_list)
        self.reorder_reflections_by_index(sorted_ind, reversed)

    def sort_reflections_by_d(self, reversed = False):
        d_list = []
        for reflection in self.reflections:
            d_list.append(reflection.d0)
        sorted_ind = np.argsort(d_list)
        self.reorder_reflections_by_index(sorted_ind, reversed)

    def has_thermal_expansion(self):
        return (self.alpha_t0!=0) or (self.d_alpha_dt!=0)


def lookup_jcpds_line(in_string,
                      pressure=0.,
                      temperature=0.,
                      path=os.getenv('JCPDS_PATH')):
    """
    Returns the d-spacing in Angstroms for a particular lattice plane.

    Inputs:
       Diffaction_plane: A string of the form 'Compound HKL', where Compound
       is the name of a material (e.g. 'gold', and HKL is the diffraction
       plane (e.g. 220).
       There must be a space between Compound and HKL.
         Examples of Diffraction_plane:
             'gold 111' - Gold 111 plane
             'si 220'   - Silicon 220 plane

    Keywords:
       path:
          The path in which to look for the file 'Compound.jcpds'.  The
          default is to search in the directory pointed to by the
          environment variable JCPDS_PATH.

       pressure:
          The pressure at which to compute the d-spacing.  Not yet
          implemented, zero pressure d-spacing is always returned.

       temperature:
           The temperature at which to compute the d-spacing.  Not yet
           implemented.  Room-temperature d-spacing is always returned.

    Outputs:
       This function returns the d-spacing of the specified lattice plane.
       If the input is invalid, e.g. non-existent compound or plane, then the
       function returns None.

    Restrictions:
       This function attempts to locate the file 'Compound.jcpds', where
       'Compound' is the name of the material specified in the input parameter
       'Diffraction_plane'.  For example:
           d = lookup_jcpds_line('gold 220')
       will look for the file gold.jcpds.  It will either look in the file
       specified in the PATH keyword parameter to this function, or in the
       the directory pointed to by the environtment variable JCPDS_PATH
       if the PATH keyword is not specified.  Note that the filename will be
       case sensitive on Unix systems, but not on Windows.

       This function is currently only able to handle HKL values from 0-9.
       The parser will need to be improved to handle 2-digit values of H,
       K or L.

    Procedure:
       This function calls jcpds.read_file() and searches for the specified HKL plane
       and returns its d-spacing;

    Example:
       d = lookup_jcpds_line('gold 111')   # Look up gold 111 line
       d = lookup_jcpds_line('quartz 220') # Look up the quartz 220 line
    """

    temp = in_string.split()
    if (len(temp) < 2): return None
    file = temp[0]
    nums = temp[1].split()
    n = len(nums)
    if (n == 1):
        if (len(nums[0]) == 3):
            try:
                hkl = (int(nums[0][0]), int(nums[0][1]), int(nums[0][2]))
            except:
                return None
        else:
            return None
    elif (n == 3):
        hkl = list(map(int, nums))
    else:
        return None

    full_file = path + file + '.jcpds'
    try:
        j = jcpds()
        j.load_file(full_file)
        refl = j.get_reflections()
        for r in refl:
            if (r.h == hkl[0] and
                        r.k == hkl[1] and
                        r.l == hkl[2]): return r.d0
        return None
    except:
        return None