
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
    

def lmfit_fix_int_data_type(fname):
    """
    fixes problem with lmfit save/load model. 
    lmfit load model cannot read int32 data with nulls in it. 
    if replace 'int32' with 'float32' it will read. 
    """
            
    ObjRead = open(fname, "r")   
    txtContent = ObjRead.read();
    ObjRead.close()
    
    txtContent = txtContent.replace('int', 'float')
    
    print('    Rewriting', fname)
    ObjRead = open(fname, "w")
    ObjRead.write(txtContent)
    ObjRead.close()
    

    
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
        out_dir = './'
            

    # create output file name from passed name
    path, filename = os.path.split(base_file_name)
    base, ext = os.path.splitext(filename)
    if not base:
        print("No base filename, using input filename instead.")
        base =  os.path.splitext(FitSettings.inputfile)[0]
        
    if base[-1:] == '_': #if the base file name ends in an '_' remove it. 
        base = base[0:-1]
        
    out_file = out_dir + base + '.dat'

    text_file = open(out_file, "w")
    print('Writing', out_file)
    
    text_file.write("# Summary of fits produced by continuous_peak_fit for input file: %s.\n" % FitSettings.inputfile)
    text_file.write("# For more information: http://www.github.com/me/something\n" )
    text_file.write("# File version: %i \n" % 1 )
    text_file.write("# \n")
    
    #write header
    col_width = 12
    text_file.write(("# {0:<"+str(col_width+5)+"}").format("Data File"+","))
    text_file.write(("{0:<"+str(col_width)+"}").format("Peak"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d0"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d0_err"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d2cos"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d2cos_err"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d2sin"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d2sin_err"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("corr coef"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("diff strain"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("diff s err"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("orientation"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("orient err"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d_max"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("d_min"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("h0"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("h0_err"+",")) 
    text_file.write(("{0:>"+str(col_width)+"}").format("w0"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("w0_err"+",")) 
    text_file.write(("{0:>"+str(col_width)+"}").format("p0"+","))
    text_file.write(("{0:>"+str(col_width)+"}").format("p0_err"+","))
    text_file.write("\n")
    
    for z in range(n_diff_files):
        
        filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        filename = filename+'.json'
        
        # Read JSON data from file
        with open(filename) as json_data:
            fit = json.load(json_data)
        
        num_subpatterns = len(FitSettings.fit_orders)
        for y in range(num_subpatterns):
            
            orders = FitSettings.fit_orders[y]
            
            subfilename = os.path.splitext(os.path.basename(diff_files[z]))[0] + '_'
            
            for x in range(len(orders['peak'])):
                if 'phase' in orders['peak'][x]:
                    subfilename = subfilename + orders['peak'][x]['phase']
                else:
                    subfilename = subfilename + "Peak"
                if 'hkl' in orders['peak'][x]:
                    subfilename = subfilename + str(orders['peak'][x]['hkl'])
                else:
                    subfilename = subfilename + x
                        
                if x < len(orders['peak']) - 1 and len(orders['peak']) > 1:
                    subfilename = subfilename + '_'
                        
            print('  Incorporating ' + subfilename)
            try:
                gmodel = load_modelresult(subfilename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
            except:
                lmfit_fix_int_data_type(subfilename+'.sav')
                try:
                    gmodel = load_modelresult(subfilename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
                except:
                    raise FileNotFoundError     
        
            # FIX ME: this will only work for one peak. Needs fixing if more then one peak in subpattern
            corr = gmodel.params['peak_0_d3'].correl['peak_0_d4']
            #ci, trace = lmfit.conf_interval(mini, gmodel, sigmas=[1, 2], trace=True)
            #lmfit.printfuncs.report_ci(ci)

            #get parameters to write to output file.
            #file name
            out_name = os.path.splitext(os.path.basename(diff_files[z]))[0]
            
            for x in range(len(orders['peak'])):
                
                out_peak = []
                
                #peak
                if 'phase' in orders['peak'][x]:
                    out_peak = orders['peak'][x]['phase']
                else:
                    out_peak = out_peak + "Peak"
                if 'hkl' in orders['peak'][x]:
                    out_peak = out_peak + ' (' + str(orders['peak'][x]['hkl']) + ')'
                else:
                    out_peak = out_peak + ' ' + str(x)
                
                # d0
                out_d0       = fit[y]['peak'][x]['d-space'][0]
                out_d0err    = fit[y]['peak'][x]['d-space_err'][0]
                #differential coefficients, errors and covarience
                out_dcos2    = fit[y]['peak'][x]['d-space'][3]
                out_dsin2    = fit[y]['peak'][x]['d-space'][4]
                out_dcos2err = fit[y]['peak'][x]['d-space_err'][3]
                out_dsin2err = fit[y]['peak'][x]['d-space_err'][4]
                out_dcorr    = gmodel.params['peak_0_d3'].correl['peak_0_d4']
                
                #differential strain
                # a.cos(2t) + b.sin(2t) = (a^2+b^2)^(1/2) . cos(2t - atan(b/a))
                # out_dd  = (a^2+b^2)^(1/2)
                out_dd    = (fit[y]['peak'][x]['d-space'][3]**2 + fit[y]['peak'][x]['d-space'][4]**2)**(1/2)
                #out_dderr= [(2.a.da.)^2 + (2.b.db)^2]^(1/2)]^(1/2)
                out_dderr = ((2 * fit[y]['peak'][x]['d-space'][3] * fit[y]['peak'][x]['d-space_err'][3])**2 + 
                             (2 * fit[y]['peak'][x]['d-space'][4] * fit[y]['peak'][x]['d-space_err'][4])**2 ) ** (1/4)
                
                #angle
                # out_ang  = atan(b/a)
                out_ang    = (np.arctan(fit[y]['peak'][x]['d-space'][4] / fit[y]['peak'][x]['d-space'][3]))
                # d (atan(c))/dc = 1/(c^2+1). c = b/a. dc = c.((da/a)^2 + (db/b)^2)^(1/2)
                # out_angerr = dc. 
                out_angerr = 1 / ((fit[y]['peak'][x]['d-space'][4] / fit[y]['peak'][x]['d-space'][3])**2 + 1) * (np.abs(out_ang) * ((fit[y]['peak'][x]['d-space_err'][3]/fit[y]['peak'][x]['d-space'][3])**2 + (fit[y]['peak'][x]['d-space_err'][4]/fit[y]['peak'][x]['d-space'][4])**2)**(1/2))
                # correction to make angle correct (otherwise out by pi)
                if fit[y]['peak'][x]['d-space'][3] < 0:
                    out_ang = out_ang + np.pi
                # force angles to be between -pi and pi
                if out_ang > np.pi/2:
                    out_ang = out_ang - np.pi
                elif out_ang < -np.pi/2:
                    out_ang = out_ang + np.pi
                # convert angles to degrees.
                out_ang    = np.rad2deg(out_ang)
                out_angerr = np.rad2deg(out_angerr)
                
                
                #differential max
                out_dmax    = out_dd * np.cos(2*np.pi/4 - np.deg2rad(out_ang)) + out_d0
                #differential min 
                out_dmin    = out_dd * np.cos(2*3*np.pi/4 - np.deg2rad(out_ang)) + out_d0
                
                #height mean
                out_h0    = fit[y]['peak'][x]['height'][0]
                out_h0err = fit[y]['peak'][x]['height_err'][0]
                #width mean
                out_w0    = fit[y]['peak'][x]['width'][0]
                out_w0err = fit[y]['peak'][x]['width_err'][0]
                #profile mean
                out_p0    = fit[y]['peak'][x]['profile'][0]
                out_p0err = fit[y]['peak'][x]['profile_err'][0]
            
            
            
                #write numbers to file
                dp = 5
                text_file.write(("{0:<"+str(col_width+7)+"}").format(out_name+","))
                text_file.write(("{0:<"+str(col_width)+"}").format(out_peak+","))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_d0))
                #text_file.write(" %10.5f," % out_d0)
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_d0err))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dcos2))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dcos2err))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dsin2))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dsin2err))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dcorr))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dd))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dderr))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_ang))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_angerr))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dmax))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_dmin))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_h0))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_h0err))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_w0))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_w0err))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_p0))
                text_file.write(("{0:"+str(col_width-1)+"."+str(dp)+"f},").format(out_p0err))
                text_file.write("\n")
            
            
            
    # col_width = 10
    # text_file.write(("# {0:<"+str(col_width+5)+"}").format("Data File"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("Peak"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("d0"+","))
    # text_file.write(("{0:<"+str(col_width)+"}").format("d0_err"+","))
            
    text_file.close()
            
            


            
            
            

