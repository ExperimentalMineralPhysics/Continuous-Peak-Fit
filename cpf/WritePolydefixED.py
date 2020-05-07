
__all__ = ['Requirements', 'WriteOutput']


import os

import numpy as np
import cpf.PeakFunctions as ff
import json


def Requirements():
    #List non-universally required parameters for writing this output type.
    
    RequiredParams = [
            #'apparently none!
            ]
    OptionalParams = [
            'ElasticProperties', #FIX ME: this needs to be included
            'Output_directory' # if no direcrtory is specified write to current directory.
            ]
    
    return RequiredParams
    
    
def WriteOutput(FitSettings, parms_dict):
# writes *.exp files required by polydefixED.
# N.B. this is a different file than that required by polydefix for monochromatic diffraction.
# This file contains a list of all the diffraction information. Hence it has to be written after the fitting as a single operation.


    FitParameters = dir(FitSettings)
    
    #Parse the required inputs. 
    base_file_name = FitSettings.datafile_Basename
    
    # diffraction patterns 
    # Identical to code in XRD_FitPatterns
    if not 'datafile_Files' in FitParameters:
        n_diff_files = FitSettings.datafile_EndNum - FitSettings.datafile_StartNum + 1
        diff_files = []
        for j in range(n_diff_files):
            # make list of diffraction pattern names
    
            #make number of pattern
            n = str(j+FitSettings.datafile_StartNum).zfill(FitSettings.datafile_NumDigit)
            
            #append diffraction pattern name and directory
            diff_files.append(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n + FitSettings.datafile_Ending)
    elif 'datafile_Files' in FitParameters:
        n_diff_files = len(FitSettings.datafile_Files)
        diff_files = []
        for j in range(n_diff_files):
            # make list of diffraction pattern names
    
            #make number of pattern
            n = str(FitSettings.datafile_Files[j]).zfill(FitSettings.datafile_NumDigit)
            
            #append diffraction pattern name and directory
            diff_files.append(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n + FitSettings.datafile_Ending)
    else:
        n_diff_files = len(FitSettings.datafile_Files)

    if 'Output_directory' in FitParameters:
        out_dir = FitSettings.Output_directory
    else:
        out_dir = '../'
            

    # create output file name from passed name
    path, filename = os.path.split(base_file_name)
    base, ext = os.path.splitext(filename)
    if base[-1:] == '_': #if the base file name ends in an '_' remove it. 
        base = base[0:-1]
    out_file = out_dir + base + '.exp'

    text_file = open(out_file, "w")
    print('Writing', out_file)

    #headers. set file version to be 1.
    text_file.write("# Experiment analysis file. to be used with PolydefixED\n")
    text_file.write("# For more information: http://merkel.zoneo.net/Polydefix/\n" )
    text_file.write("# File Created by FPF_WritePolydefixED function in ???\n")
    text_file.write("# For more information: http://www.github.com/me/something\n" )
    text_file.write("# File version\n" )
    text_file.write("     1\n" )


    # write two theta angle.
    # Writes two theta from the first detector.
    # FIX ME: should this be the average of all the two theat angles?
    text_file.write("# 2 theta angle (in degrees)\n")
    text_file.write("%12.5f\n" % parms_dict['calibs'].mcas[0].calibration.two_theta)                
       
    # Write detector properties.             
    text_file.write("# Number of detector positions\n")
    text_file.write("%6i\n" % len(parms_dict['azimuths']))
    text_file.write("# Angles for detector positions: Number. Use (1/0). Angle. \n")
    for x in range(len(parms_dict['azimuths'])):
        
        # determine if detector masked
        if FitSettings.Calib_mask:
            if x+1 in FitSettings.Calib_mask:
                use = 0
            else:
                use = 1
                
        # print detector number, if it is being used, then angle, lastly just a number.
        text_file.write("%6i  %i  %7.3f  %i \n" % (x+1, use, parms_dict['azimuths'][x], 1))
        
        # FIX ME: if 10 element detector then turn detector 10 off.
        
        
    #Write number of peaks
    num_subpatterns = len(FitSettings.fit_orders)
    num_peaks = 0
    for x in range(num_subpatterns):
        num_peaks = num_peaks + len(FitSettings.fit_orders[x]['peak'])
    text_file.write("# Number of peaks\n")
    text_file.write("%6i \n" % num_peaks)
    
    #Write peak properties
    text_file.write("# Peak information \n")
    text_file.write("# Number. use (1/0). h. k. l\n")
    peak = 1
    for x in range(num_subpatterns):
        for y in range(len(FitSettings.fit_orders[x]['peak'])):
            if 'hkl' in FitSettings.fit_orders[x]['peak'][y]:
                hkl = str(FitSettings.fit_orders[x]['peak'][y]['hkl'])
                
                #check if the d-spacing fits are NaN or not. if NaN switch off (0) otherwise use (1).
                # Got first file name        
                filename = os.path.splitext(os.path.basename(diff_files[0]))[0]
                filename = filename+'.json'
                # Read JSON data from file
                with open(filename) as json_data:
                    fit = json.load(json_data)
                #check if the d-spacing fits are NaN or not. if NaN switch off.
                if np.isnan(fit[x]['peak'][y]['d-space'][0]): 
                    use = 0
                else:
                    use = 1
            else:
                hkl = '000'
                use = 0
            pos = 0
            if hkl[0] == '-':
                h = hkl[pos:pos+2]
                pos = pos+2
            else:
                h = hkl[pos:pos+1]
                pos = pos+1
            if hkl[pos] == '-':
                k = hkl[pos:pos+2]
                pos = pos+2
            else:
                k = hkl[pos:pos+1]
                pos = pos+1
            if hkl[pos] == '-':
                l = hkl[pos:pos+2]
                pos = pos+2
            else:
                l = hkl[pos:pos+1]
                pos = pos+1
            text_file.write(" %5i %5i    %s    %s    %s\n" % (peak, use, h, k, l))
            peak = peak+1
    
    # Write fetting behaviour
    text_file.write("# Fit offset for maximum stress 1 for yes. 0 for no\n")
    text_file.write("     1\n")
    text_file.write("# Starting offset value. in degrees\n")
    text_file.write("%6i \n" % parms_dict['azimuths'][0])
    
    #Write material properties
    # FIX ME: need to be able to import material properties and write them to the file without messing about.
    text_file.write("# Material properties set (1/0)\n")
    if 'Material' in parms_dict:
        text_file.write("     1\n")
    else:
        text_file.write("     1\n")
        
        # FIX ME: here I am just writing the elastic properties of olivine with no regard for the structure of the file or the avaliable data. 
        #         This should really be fixed. 
        text_file.write("# Material properties\n")
        text_file.write("# Version\n")
        text_file.write("2       \n")
        text_file.write("# Name\n")
        text_file.write("Olivine\n")
        text_file.write("# Symmetry\n")
        text_file.write("ortho  \n")
        text_file.write("# EOS stuff (v0. k0. dk0/P. dK0/dT) (Li Li2006)\n")
        text_file.write("# Thermal expansion coefficients (a. b. c)\n")
        text_file.write("# 2.69000e-05 2.12000e-08 0.00000\n")
        text_file.write("# EOS stuff (v0. k0. dk0/P. dK0/dT)       \n")
        text_file.write(" 292.58 128.800 4.2000 0.00000         \n")
        text_file.write("# Thermal expansion coefficients (a. b. c) (ISAAK ET AL.. 1989)   \n")
        text_file.write(" 3.40e-05 0.7100e-08 0.00000          \n")
        text_file.write("# Elastic model           \n")
        text_file.write("1             \n")
        text_file.write("# Parameters for isotropic elastic model (g0. g1. g2. g3. g4)   \n")
        text_file.write("# G = G0 + G1*P + G2*P*P + G3*(T-300) + G4*(T-300)*(T-300)  \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000        \n")
        text_file.write("# Parameters for anisotropic elastic model (Cij) (ISAAK ET AL.. 1989)  \n")
        text_file.write(" 323.700 66.4000 71.6000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 66.4000 197.600 67.2300 0.00000 0.00000 0.00000       \n")
        text_file.write(" 71.6000 67.2300 235.100 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 64.6200 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 78.0500 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 79.0400       \n")
        text_file.write("# Parameters for anisotropic elastic model (dCij/dp) (ZHA. 1996)     \n")
        text_file.write(" 7.98000 4.74000 4.48000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 4.74000 6.37000 3.76000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 4.48000 3.76000 6.38000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 2.17000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 1.64000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 2.31000       \n")
        text_file.write("# Parameters for anisotropic elastic model (d2Cij/dp2)       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write("# Parameters for anisotropic elastic model (dCij/dT) (ISAAK ET AL.. 1989)   \n")
        text_file.write(" -0.0340000 -0.0105000 -0.00940000 0.00000 0.00000 0.00000       \n")
        text_file.write(" -0.0105000 -0.0285000 -0.00510000 0.00000 0.00000 0.00000       \n")
        text_file.write(" -0.00940000 -0.00510000 -0.0286000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 -0.0128000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 -0.0130000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 -0.0157000       \n")
        text_file.write("# Parameters for anisotropic elastic model (d2Cij/dT2)       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000       \n")
        text_file.write(" 0.00000 0.00000 0.00000 0.00000 0.00000 0.00000      \n")
    
    
    # Write details of experiment.
    
    # write number of steps
    text_file.write("# Number of time steps in the experiment\n")
    text_file.write("%i \n" %n_diff_files)
    text_file.write("# Information on time step       \n")
    text_file.write("# Step number.  Step name.         Step time.    Step temperature (K).    taux de def E\n")
    for x in range(n_diff_files):
        #FIX ME: if exists write the input for strain or temperature.
        
        #strip directory and ending from filename
        diff_name = os.path.splitext(os.path.basename(diff_files[x]))[0]
        
        text_file.write("%8i        %s  %3.4g                  %3.4g\n" %  (x+1, "{:<18}".format(diff_name), x, 25+273))
        

    #Write the experimental data fits.
    text_file.write("# Experimental data \n")
    text_file.write("# Peak number.  h.  k.  l.  d-spacing.  intensity.  detector number. step number E t \n")
    for z in range(n_diff_files):
        
        filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        filename = filename+'.json'
        
        # Read JSON data from file
        with open(filename) as json_data:
            fit = json.load(json_data)
        
        peak = 1
        for x in range(num_subpatterns):
            for y in range(len(FitSettings.fit_orders[x]['peak'])):
                
                # Write the following structure to the data file.    
                #[peak number     H     K     L     d-spacing     Intensity     detector number     step number   ]  repeat for the number of 
                #...                                                                                ]  peaks*number detectors*number patterns
                #[peak number     H     K     L     d-spacing     Intensity     detector number     step number   ]  in data set.

                if 'hkl' in FitSettings.fit_orders[x]['peak'][y]:
                    hkl = str(FitSettings.fit_orders[x]['peak'][y]['hkl'])
                    use = 1
                else:
                    hkl = '000'
                    use = 0
                pos = 0
                if hkl[0] == '-':
                    h = hkl[pos:pos+2]
                    pos = pos+2
                else:
                    h = hkl[pos:pos+1]
                    pos = pos+1
                if hkl[pos] == '-':
                    k = hkl[pos:pos+2]
                    pos = pos+2
                else:
                    k = hkl[pos:pos+1]
                    pos = pos+1
                if hkl[pos] == '-':
                    l = hkl[pos:pos+2]
                    pos = pos+2
                else:
                    l = hkl[pos:pos+1]
                    pos = pos+1
                                
                if 'symmetry' in FitSettings.fit_orders[x]['peak'][y]:
                    sym = FitSettings.fit_orders[x]['peak'][y]['symmetry']
                else:
                    sym = 1
                    
                for w in range(len(parms_dict['calibs'].mcas)):
                    
                    # determine if detector masked
                    if FitSettings.Calib_mask:
                        if w+1 in FitSettings.Calib_mask:
                            use = 0
                        else:
                            use = 1
                    else:
                            use = 1
                            
                    #only write output for unmasked detectors
                    if use == 1:
                        az = parms_dict['azimuths'][w]
                        
                        peak_d = ff.Fourier_expand((az), fit[x]['peak'][y]['d-space'])
                        #peak_tth = 2.*np.degrees(np.arcsin(wavelength/2/peak_d))
                        peak_i = ff.Fourier_expand((az)*sym, fit[x]['peak'][y]['height'])  #FIX ME - Is this the height of the peak or the integral under it?
                        #peak_w = ff.Fourier_expand((az)*sym, fit['peak'][y]['width'])   #FIX ME - is this the correct half width?
                        #peak_p = ff.Fourier_expand((az)*sym, fit['peak'][y]['profile']) #FIX ME - is this 1 or 0 for Gaussian?
                        text_file.write("%8i        %s   %s   %s   %8.4f  %10.3f       %3i           %3i\n" % (peak, h, k, l, peak_d, peak_i, w+1,  z+1))
                
                peak = peak+1
        
    text_file.close()