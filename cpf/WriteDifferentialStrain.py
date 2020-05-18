
__all__ = ['Requirements', 'WriteOutput']


import os

import numpy as np
import cpf.PeakFunctions as ff
import json
import lmfit
from lmfit.model import load_modelresult

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
            

    for z in range(n_diff_files):
        num_subpatterns = len(FitSettings.fit_orders)
        for y in range(num_subpatterns):
            
            orders = FitSettings.fit_orders[y]
            
            filename = os.path.splitext(os.path.basename(diff_files[z]))[0] + '_'
            
            for x in range(len(orders['peak'])):
                if 'phase' in orders['peak'][x]:
                    filename = filename + orders['peak'][x]['phase']
                else:
                    filename = filename + "Peak"
                if 'hkl' in orders['peak'][x]:
                    filename = filename + str(orders['peak'][x]['hkl'])
                else:
                    filename = filename + x
                        
                if x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
                    filename = filename + '_'
                        
            gmodel = load_modelresult(filename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
        
            corr = gmodel.params['peak_0_d3'].correl['peak_0_d4']
            print(filename, ': correlation of differential stress coefficents = ', corr)
            #ci, trace = lmfit.conf_interval(mini, gmodel, sigmas=[1, 2], trace=True)
            #lmfit.printfuncs.report_ci(ci)

    

