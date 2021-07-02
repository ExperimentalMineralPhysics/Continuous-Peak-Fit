
__all__ = ['Requirements', 'WriteOutput']


import os

import numpy as np
import cpf.PeakFunctions as ff
import json
import lmfit
from lmfit.model import load_modelresult
import cpf.XRD_FitPattern as XRD_FP

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
    
    txtContent = txtContent.replace('uint', 'float')
    txtContent = txtContent.replace('int', 'float')
    
    print('    Rewriting', fname)
    ObjRead = open(fname, "w")
    ObjRead.write(txtContent)
    ObjRead.close()
    

    
def WriteOutput(FitSettings, parms_dict, debug=True, **kwargs):
# writes some of the fitted coeficients to a table. With a focus on the differential strain coefficents


    FitParameters = dir(FitSettings)
    
    #Parse the required inputs. 
    base_file_name = FitSettings.datafile_Basename
    
    # diffraction patterns 
    diff_files, n_diff_files = XRD_FP.FileList(FitParameters, FitSettings)
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
    text_file.write("# File version: %i \n" % 1.1 )
    # version 1.1 has elasped time and reduced chi squared added to the output table.
    text_file.write("# \n")
    
    #write header
    width_col = 12
    dp = 5 #used later or the number of decial places in the numbers -- but set here so all settings are in the same place.
    width_fnam = 25
    width_hkl = 15
    text_file.write(("# {0:<"+str(width_fnam-2)+"}").format("Data File"+","))
    text_file.write(("{0:<"+str(width_hkl)+"}").format("Peak"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d0"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d0_err"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d2cos"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d2cos_err"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d2sin"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d2sin_err"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("corr coef"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("diff strain"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("diff s err"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("orientation"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("orient err"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d_max"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("d_min"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("mean h"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("h_err"+",")) 
    text_file.write(("{0:>"+str(width_col)+"}").format("mean w"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("w_err"+",")) 
    text_file.write(("{0:>"+str(width_col)+"}").format("mean p"+","))
    text_file.write(("{0:>"+str(width_col)+"}").format("p0_err"+","))
    if debug:
        text_file.write(("{0:>"+str(width_col)+"}").format("Time taken"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Chunk time"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Sum Resid^2"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Status"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Func eval"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Num vars"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Num data"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Deg Freedom"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("ChiSq"+","))
        text_file.write(("{0:>"+str(width_col)+"}").format("Red. ChiSq"+","))
        text_file.write(("{0:<"+str(width_col)+"}").format("Akaike Information Criterion"+","))
        text_file.write(("{0:<"+str(width_col)+"}").format("Bayesian Information Criterion"+","))
    text_file.write("\n")
    
    for z in range(n_diff_files):
        
        filename = os.path.splitext(os.path.basename(diff_files[z]))[0]
        filename = filename+'.json'
        
        
        if os.path.isfile(filename):
            # Read JSON data from file
            with open(filename) as json_data:
                fit = json.load(json_data)
            
            num_subpatterns = len(fit)
            for y in range(num_subpatterns):
                
                try:
                    orders = FitSettings.fit_orders[y]
                except:
                    orders = []    
                
                subfilename = os.path.splitext(os.path.basename(diff_files[z]))[0] + '_'
                
                for x in range(len(fit[y]['peak'])):
                    if 'phase' in fit[y]['peak'][x]:
                        subfilename = subfilename + fit[y]['peak'][x]['phase']
                    elif 'phase' in orders['peak'][x]:
                        subfilename = subfilename + orders['peak'][x]['phase']
                    else:
                        subfilename = subfilename + "Peak"
                    if 'hkl' in fit[y]['peak'][x]:
                        subfilename = subfilename + str(fit[y]['peak'][x]['hkl'])
                    elif 'hkl' in orders['peak'][x]:
                        subfilename = subfilename + orders['peak'][x]['hkl']
                    else:
                        subfilename = subfilename + x
                            
                    if x < len(fit[y]['peak']) - 1 and len(fit[y]['peak']) > 1:
                        subfilename = subfilename + '_'
                            
                print('  Incorporating ' + subfilename)
                if os.path.isfile(subfilename+'.sav'):
                    try:
                        gmodel = load_modelresult(subfilename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
                    except:
                        lmfit_fix_int_data_type(subfilename+'.sav')
                        try:
                            gmodel = load_modelresult(subfilename+'.sav', funcdefs={'PeaksModel': ff.PeaksModel})
                        except:
                            raise FileNotFoundError  
                    # FIX ME: this will only work for one peak. Needs fixing if more then one peak in subpattern
                    try:
                        corr = gmodel.params['peak_0_d3'].correl['peak_0_d4']
                    except:
                        corr = np.nan
                    #ci, trace = lmfit.conf_interval(mini, gmodel, sigmas=[1, 2], trace=True)
                    #lmfit.printfuncs.report_ci(ci)
    
                else:
                    #print('  No file ' + subfilename)
                    corr = np.nan
            
                #get parameters to write to output file.
                #file name
                out_name = os.path.splitext(os.path.basename(diff_files[z]))[0]

                width_fnam = np.max((width_fnam, len(out_name)))
            
                for x in range(len(fit[y]['peak'])):
                    
                    out_peak = []
                    
                    if fit[y]['peak'][x]['d-space-type'] != 'fourier':
                        raise ValueError('This output type is not setup to process non-Fourier peak centroids.')
                    
                    #peak
                    if 'phase' in fit[y]['peak'][x]:
                        out_peak = fit[y]['peak'][x]['phase']
                    elif 'phase' in orders['peak'][x]:
                        out_peak = orders['peak'][x]['phase']
                    else:
                        out_peak = out_peak + "Peak"
                    if 'hkl' in fit[y]['peak'][x]:
                        out_peak = out_peak + ' (' + str(fit[y]['peak'][x]['hkl']) + ')'
                    elif 'hkl' in orders['peak'][x]:
                        out_peak = out_peak + orders['peak'][x]['hkl']
                    else:
                        out_peak = out_peak + ' ' + str(x)
                    
                    width_hkl = np.max((width_hkl, len(out_peak)))
                    
                    # d0
                    out_d0       = fit[y]['peak'][x]['d-space'][0]
                    out_d0err    = fit[y]['peak'][x]['d-space_err'][0]
                    if out_d0err is None:  #catch  'null' as an error 
                        out_d0err = np.nan
                        
                    #differential coefficients, errors and covarience
                    try:
                        out_dcos2    = fit[y]['peak'][x]['d-space'][3]
                        out_dsin2    = fit[y]['peak'][x]['d-space'][4]
                        out_dcos2err = fit[y]['peak'][x]['d-space_err'][3]
                        out_dsin2err = fit[y]['peak'][x]['d-space_err'][4]
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
                    except:
                        out_dcos2    = np.nan
                        out_dsin2    = np.nan
                        out_dcos2err = np.nan
                        out_dsin2err = np.nan
                        out_dd       = np.nan
                        out_dderr    = np.nan
                        out_ang      = np.nan
                        out_angerr   = np.nan
                        
                    out_dcorr    = corr#gmodel.params['peak_0_d3'].correl['peak_0_d4']
                    
                    
                    
                    #differential max
                    out_dmax    = out_dd * np.cos(2*np.pi/4 - np.deg2rad(out_ang)) + out_d0
                    #differential min 
                    out_dmin    = out_dd * np.cos(2*3*np.pi/4 - np.deg2rad(out_ang)) + out_d0
                    
                    #height mean
                    if fit[y]['peak'][x]['height-type'] == 'fourier':
                        out_h0    = fit[y]['peak'][x]['height'][0]
                        out_h0err = fit[y]['peak'][x]['height_err'][0]
                    else:
                        tot    = np.sum(fit[y]['peak'][x]['height'])
                        errsum = np.sqrt(np.sum(np.array(fit[y]['peak'][x]['height_err'])**2))
                        num    = np.shape(fit[y]['peak'][x]['height'])
                        out_h0 = float(tot/num)
                        out_h0err = float(errsum/num)
                    if out_h0 is None:  #catch  'null' as an error 
                        out_h0 = np.nan
                    if out_h0err is None:  #catch  'null' as an error 
                        out_h0err = np.nan
                        
                    #width mean
                    if fit[y]['peak'][x]['width-type'] == 'fourier':
                        out_w0    = fit[y]['peak'][x]['width'][0]
                        out_w0err = fit[y]['peak'][x]['width_err'][0]
                    else:
                        tot    = np.sum(fit[y]['peak'][x]['width'])
                        errsum = np.sqrt(np.sum(np.array(fit[y]['peak'][x]['width_err'])**2))
                        num    = np.shape(fit[y]['peak'][x]['width'])
                        out_w0 = float(tot/num)
                        out_w0err = float(errsum/num)
                    if out_w0 is None:  #catch  'null' as an error 
                        out_w0 = np.nan
                    if out_w0err is None:  #catch  'null' as an error 
                        out_w0err = np.nan
               
                    #profile mean
                    if fit[y]['peak'][x]['profile-type'] == 'fourier':
                        out_p0    = fit[y]['peak'][x]['profile'][0]
                        out_p0err = fit[y]['peak'][x]['profile_err'][0]
                    else:
                        tot    = np.sum(fit[y]['peak'][x]['profile'])
                        errsum = np.sqrt(np.sum(np.array(fit[y]['peak'][x]['profile_err'])**2))
                        num    = np.shape(fit[y]['peak'][x]['profile'])
                        out_p0 = float(tot/num)
                        out_p0err = float(errsum/num)
                    if out_p0 is None:  #catch  'null' as an error 
                        out_p0 = np.nan
                    if out_p0err is None:  #catch  'null' as an error 
                        out_p0err = np.nan
                
                    #write numbers to file
                    text_file.write(("{0:<"+str(width_fnam)+"}").format(out_name+","))
                    text_file.write(("{0:<"+str(width_hkl)+"}").format(out_peak+","))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_d0))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_d0err))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dcos2))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dcos2err))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dsin2))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dsin2err))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dcorr))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dd))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dderr))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_ang))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_angerr))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dmax))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_dmin))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_h0))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_h0err))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_w0))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_w0err))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_p0))
                    text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(out_p0err))
                    if debug:
                        # include properties from the lmfit output that were passed with the fits.
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(fit[y]['FitProperties']['time-elapsed']))
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(fit[y]['FitProperties']['chunks-time']))
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp-1)+"e},").format(fit[y]['FitProperties']["sum-residuals-squared"]))
                        text_file.write(("{0:"+str(width_col-1)+"d},").format(fit[y]['FitProperties']['status']))
                        text_file.write(("{0:"+str(width_col-1)+"d},").format(fit[y]['FitProperties']['function-evaluations']))
                        text_file.write(("{0:"+str(width_col-1)+"d},").format(fit[y]['FitProperties']["n-variables"]))
                        text_file.write(("{0:"+str(width_col-1)+"d},").format(fit[y]['FitProperties']["n-data"]))
                        text_file.write(("{0:"+str(width_col-1)+"d},").format(fit[y]['FitProperties']["degree-of-freedom"]))
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp-1)+"e},").format(fit[y]['FitProperties']["ChiSq"]))
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(fit[y]['FitProperties']['RedChiSq']))
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(fit[y]['FitProperties']["aic"]))
                        text_file.write(("{0:"+str(width_col-1)+"."+str(dp)+"f},").format(fit[y]['FitProperties']["bic"]))
                    text_file.write("\n")
                
            
    text_file.close()