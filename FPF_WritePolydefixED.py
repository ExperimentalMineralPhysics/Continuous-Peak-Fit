
import os

import numpy as np
import FPF_PeakFourierFunctions as ff
import json


def WritePolydefixED(FitSettings, parms_dict):
# writes *.exp files required by polydefixED.
# N.B. this is a different file than that required by polydefix for monochromatic diffraction.
# This file contains a list of all the diffraction information. Hence it has to be written after the fitting as a single operation.

    
    #Parse the required inputs. 
    base_file_name = FitSettings.datafile_Basename
    
    # diffraction patterns 
    # Identical to code in XRD_FitPatterns
    if not FitSettings.datafile_Files:
        n_diff_files = FitSettings.datafile_EndNum - FitSettings.datafile_StartNum + 1
        diff_files = []
        for j in range(n_diff_files):
            # make list of diffraction pattern names
    
            #make number of pattern
            n = str(j+FitSettings.datafile_StartNum).zfill(FitSettings.datafile_NumDigit)
            
            #append diffraction pattern name and directory
            diff_files.append(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n + FitSettings.datafile_Ending)
    elif FitSettings.datafile_Files:
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



    # create output file name from passed name
    path, filename = os.path.split(base_file_name)
    base, ext = os.path.splitext(filename)
    if base[-1:] == '_': #if the base file name ends in an '_' remove it. 
        base = base[0:-1]
    out_file = base + '.exp'

    text_file = open(out_file, "w")


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
        # print detector number, if it is being used, then angle, lastly just a number.
        text_file.write("%6i  " % (x+1))
        text_file.write("%i  " % (1))
        text_file.write("%7.3f " % (parms_dict['azimuths'][x]))
        if FitSettings.Calib_mask:
            if x+1 == 10:
                use = 0
            else:
                use = 1
        text_file.write("%i\n" % use )
        #text_file.write('%i  %i  %7.3f i% \n' % (x+1, 1, parms_dict['azimuths'][x], 1))

        # FIX ME: make into a single line for writing. I dont know why it doesnt work
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
                    
                    az = parms_dict['azimuths'][w]
                    
                    peak_d = ff.Fourier_expand((az), fit[x]['peak'][y]['d-space'])
                    #peak_tth = 2.*np.degrees(np.arcsin(wavelength/2/peak_d))
                    peak_i = ff.Fourier_expand((az)*sym, fit[x]['peak'][y]['height'])  #FIX ME - Is this the height of the peak or the integral under it?
                    #peak_w = ff.Fourier_expand((az)*sym, fit['peak'][y]['width'])   #FIX ME - is this the correct half width?
                    #peak_p = ff.Fourier_expand((az)*sym, fit['peak'][y]['profile']) #FIX ME - is this 1 or 0 for Gaussian?
    
        
                        
                    text_file.write("%8i        %s   %s   %s   %8.4f  %10.3f       %3i           %3i\n" % (peak, h, k, l, peak_d, peak_i, w+1,  z+1))
                
                peak = peak+1
        
    

    text_file.close()



# =============================================================================
# This is the test case for the monochromatic output file.
#
# #test cases for wirting files
# test = 0
# 
# if test==1:
#     Fouriers_to_write = [{'peak': [{
#     'profile': 0, 
#     'width': [ 0.03296637], 
#     'd-space': [  1.45521770e+02,   3.15812130e-03,   3.63949255e-03,   -1.43496474e+02,  -1.43500658e+02], 
#     'height': [ -1.61854953e+08,  -7.90551608e+01,  -1.86863466e+02,
#      1.84245125e+08,  -5.01600833e+05,   8.75674029e+02,
#      3.71825217e+03,  -5.17425137e+06,   3.65458422e+08,
#     -3.90556882e+03,  -2.25773282e+04,   1.86844720e+08,
#     -2.28359051e+08,   8.22152351e+03,   5.67280685e+04,
#     -3.07040812e+08,   9.14094789e+06,  -8.05594617e+03,
#     -6.26194032e+04,   1.24648536e+08,  -5.65124620e+06,
#      2.94882628e+03,   2.51878566e+04,  -2.16682847e+07,
#      2.17678521e+07]}], 
#     'background': [4.0],
#     'range': [[11.700006, 12.079979], [1.9937644, 2.0582786]]}, 
#     {'peak': [{
#     'profile': 0, 
#     'width': [ 0.06473334], 
#     'd-space': [ -1.48371656e+02,   1.32237912e-03,   1.75463473e-03,   1.49805560e+02,   1.49801480e+02], 
#     'height': [  2.05724056e+07,  -9.25144706e+00,  -5.38376022e+00,
#     -2.10599065e+07,  -6.44172806e+06,   1.50563973e+02,
#      9.46366683e+01,  -6.81122757e+05,  -3.99125724e+07,
#     -7.06638114e+02,  -4.55423181e+02,  -6.71592503e+06,
#      3.26600155e+07,   1.31413614e+03,   9.36273739e+02,
#     -1.63116330e+06,   3.30713638e+06,  -1.04632954e+03,
#     -9.05402148e+02,   1.44412561e+07,  -1.51098537e+07,
#      2.97575768e+02,   3.39659190e+02,  -4.92552756e+06,
#      4.92459547e+06]}], 
#     'background': [4.0],
#     'range': [[16.600002, 17.099998], [1.4110823, 1.4532731]]}, 
#      {'peak': [{
#      'profile': 0, 
#      'width': [ 0.05999371], 
#      'd-space': [ -1.38852649e+02,   8.75848195e-04,   1.30742495e-03,  1.40024571e+02,   1.40021371e+02], 
#      'height': [ -2.94680617e+05,  -7.51087348e+00,  -5.01955745e-01,
#      1.71972262e+05,   3.01886458e+06,   1.47013354e+02,
#     -3.97875774e+01,   1.00812116e+06,   1.31668009e+06,
#     -8.64958125e+02,   4.82478179e+02,  -3.80495454e+07,
#      2.73111970e+06,   2.04088806e+03,  -1.76521792e+03,
#      9.43916850e+07,  -4.37268233e+07,  -2.05402403e+03,
#      2.45598540e+03,  -8.07779850e+07,   6.04951412e+07,
#      7.38560951e+02,  -1.13116306e+03,   2.35504567e+07,
#     -2.35403011e+07]}], 
#      'background': [4.0],
#      'range': [[20.400015, 20.89999], [1.1566441, 1.1846806]]}]
# 
#     file_name = "Output1.txt"
# 
#     Azimuths = 360.
# 
#     WriteMultiFit(file_name, Fouriers_to_write, Azimuths)
# 
# 
#     #Write the fits to a temporary file
#     TempFilename = open("Fit_previous_JSON.dat", "w")
# 
#     # Write a JSON string into the file.
#     json_string = json.dumps(Fouriers_to_write, TempFilename, sort_keys=True, indent=4, separators=(',', ': '))
#     TempFilename.write(json_string)
# 
#     
# =============================================================================
