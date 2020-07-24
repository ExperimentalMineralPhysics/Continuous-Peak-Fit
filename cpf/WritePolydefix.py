
__all__ = ['Requirements', 'WriteOutput']


import os

import numpy as np
#import cpf.PeakFunctions as ff
import cpf.WriteMultiFit as WriteMultiFit
import json
#import re
#import datetime

def Requirements():
    #List non-universally required parameters for writing this output type.
    
    RequiredParams = [
            #'apparently none!
            ]
    OptionalParams = [
            'ElasticProperties', #FIX ME: this needs to be included
            'tc', # which thermocoule to include from 6BMB/X17B2 collection system. default to 1. #FIX ME: this needs to be included
            'Output_directory' # if no direcrtory is specified write to current directory.
            ]
    
    return RequiredParams
    
    
def WriteOutput(FitSettings, parms_dict):
# writes *.fit files produced by multifit (as input for polydefix)
# writes *.exp files required by polydefix.
# N.B. this is a different file than that required by polydefix for energy dispersive diffraction.

    
    #Write fit files
    WriteMultiFit.WriteOutput(FitSettings, parms_dict)



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
        out_dir = './'
            

    # create output file name from passed name
    path, filename = os.path.split(base_file_name)
    base, ext = os.path.splitext(filename)
    if base[-1:] == '_': #if the base file name ends in an '_' remove it. 
        base = base[0:-1]
    out_file = out_dir + base + '.exp'

    text_file = open(out_file, "w")
    print('Writing', out_file)

    #headers. set file version to be 1.
    text_file.write("# Experiment analysis file. to be used with Polydefix\n")
    text_file.write("# For more information: http://merkel.zoneo.net/Polydefix/\n" )
    text_file.write("# File Created by FPF_WritePolydefix function in Continuous-Peak-Fit\n")
    text_file.write("# For more information: http://www.github.com/me/something\n" )
    
    
    text_file.write("# File version\n" )
    text_file.write("     2\n" )

    #Write data properties
    text_file.write("# Directory with FIT files\n" )
    # needs absolute path for the data files and the string has to end with a '/' e.g. Users/me/data/BCC1_2GPa_10s_e/
    text_file.write("     %s/\n" %  os.path.abspath(os.getcwd() + '/' + FitSettings.Output_directory) ) 
    text_file.write("# Basename for FIT files\n" )
    #if last symbol in the file name is '_' then we need to strip it from the name
    text_file.write("     %s\n" %  FitSettings.datafile_Basename.strip('_') )
    text_file.write("# First index for FIT files\n" )
    text_file.write("     %i\n" %  FitSettings.datafile_StartNum )
    text_file.write("# Last index for FIT files\n" )
    text_file.write("     %i\n" %  FitSettings.datafile_EndNum ) 
    text_file.write("# Number of digits for FIT files\n" )
    text_file.write("     %i\n" %  FitSettings.datafile_NumDigit )
    text_file.write("# Wavelength\n" )
    text_file.write("     %8.7g\n" %  parms_dict['conversion_constant'] )
    text_file.write("# Fit offset for maximum stress 1 for yes, 0 for no\n" )
    text_file.write("     %i\n" %  1 )
    text_file.write("# Starting offset value, in degrees\n" )
    text_file.write("     %8.4g\n" %  10.0000 )
    text_file.write("# Fit beam center; 1 for yes, 0 for no\n" )
    text_file.write("     %i\n" %  0 )
    text_file.write("# Material properties set (1/0)\n" )
    text_file.write("     %i\n" %  1 )
    text_file.write("# Peaks properties set (1/0)\n" )
    text_file.write("     %i\n" %  1 )
    
    #write number of peaks and hkls.
    numpeaks = 0
    for x in range(len(FitSettings.fit_orders)):
        numpeaks = numpeaks + len(FitSettings.fit_orders[x]['peak'])
    text_file.write("# Number of peaks\n" )
    text_file.write("     %i\n" %  numpeaks )
    text_file.write("# Peaks info (use, h, k, l)\n" )
    for x in range(len(FitSettings.fit_orders)):
        for y in range(len(FitSettings.fit_orders[x]['peak'])):
            if 'hkl' in FitSettings.fit_orders[x]['peak'][y]:
                hkl = str(FitSettings.fit_orders[x]['peak'][y]['hkl'])
                
                if hkl == '0' or hkl == 0:
                    hkl = '000'
                    
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
                    
                #check if the  phase name matches the material (if exists)
                if 'Phase' or 'Output_ElasticProperties' in FitParameters: # FIX ME: should catch differences between output elastic properties and phases. 
                    if FitSettings.phase == FitSettings.fit_orders[x]['peak'][y]['phase']:
                        use = use
                    else:
                        use = 0
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
            text_file.write(" %5i    %s    %s    %s\n" % (use, h, k, l))
            
    #material properties
    text_file.write("# Material properties\n" )
    if 'Material' in FitParameters:
        
        raise ValueError("''Material'' is depreciated as an option. Change your input file to have a ''Phase'' instead.")
        
    elif 'Output_ElasticProperties' in FitParameters:
        
        fid = open(FitSettings.Output_ElasticProperties, 'r')
        
        #pipe ealstic properties to the output file.
        text_file.write(fid.read())
    
    elif 'phase' in FitParameters:
        
        try: 
            fid = open((FitSettings.phase + '_elastic_properties.txt'), 'r')
            #pipe ealstic properties to the output file.
            text_file.write(fid.read())
        except:
            raise ValueError("The file " + FitSettings.phase + "_elastic_properties.txt is not found. Edit input file and try again.")
            
    else:
        
        # FIX ME: here I am just writing the elastic properties of BCC iron with no regard for the structure of the file or the avaliable data. 
        #         This should really be fixed. 
        text_file.write("# Name:  (N.B. this is the default phase)\n")
        text_file.write("BCC-Fe\n")
        text_file.write("# Symmetry\n")
        text_file.write("cubic  \n")
        text_file.write("# EOS stuff (v0, k0, k'0)\n")
        text_file.write("      23.5530      159.900      6.52000\n")
        text_file.write("# Elastic model\n")
        text_file.write("       1\n")
        text_file.write("# Parameters for isotropic elastic model (k0 k1 k2, g0, g1, g2)\n")
        text_file.write("      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000\n")
        text_file.write("# Parameters for anisotropic elastic model (Cij)\n")
        text_file.write("      223.000      127.000      127.000      0.00000      0.00000      0.00000\n")
        text_file.write("      127.000      223.000      127.000      0.00000      0.00000      0.00000\n")
        text_file.write("      127.000      127.000      223.000      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      122.000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      122.000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      122.000\n")
        text_file.write("# Parameters for anisotropic elastic model (dCij/dp)\n")
        text_file.write("      6.12245      4.08163      4.08163      0.00000      0.00000      0.00000\n")
        text_file.write("      4.08163      6.12245      4.08163      0.00000      0.00000      0.00000\n")
        text_file.write("      4.08163      4.08163      6.12245      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      2.44898      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      2.44898      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      2.44898\n")
        text_file.write("# Parameters for anisotropic elastic model (d2Cij/dp2)\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n")
        text_file.write("      0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n")
        text_file.write("     0.00000      0.00000      0.00000      0.00000      0.00000      0.00000\n")

        
    text_file.close()
   
   