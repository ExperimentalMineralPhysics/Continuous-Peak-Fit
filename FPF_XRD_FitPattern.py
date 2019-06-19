# /usr/bin/python

import numpy as np
import numpy.ma as ma
import copy, os, sys
import json
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.path as mlp
from scipy.optimize import curve_fit
np.set_printoptions(threshold='nan')

from FPF_XRD_FitSubpattern import FitSubpattern
import FPF_PeakFourierFunctions as ff
import FPF_WriteMultiFit as wr

# Needed for JSON to save fitted parameters. 
# copied from https://stackoverflow.com/questions/3488934/simplejson-and-numpy-array#24375113
# on 13th November 2018
def json_numpy_serialzer(o):
    """ Serialize numpy types for json

    Parameters:
        o (object): any python object which fails to be serialized by json

    Example:

        >>> import json
        >>> a = np.array([1, 2, 3])
        >>> json.dumps(a, default=json_numpy_serializer)

    """
    numpy_types = (
        np.bool_,
        # np.bytes_, -- python `bytes` class is not json serializable     
        # np.complex64,  -- python `complex` class is not json serializable  
        # np.complex128,  -- python `complex` class is not json serializable
        # np.complex256,  -- special handling below
        # np.datetime64,  -- python `datetime.datetime` class is not json serializable
        np.float16,
        np.float32,
        np.float64,
        # np.float128,  -- special handling below
        np.int8,
        np.int16,
        np.int32,
        np.int64,
        # np.object_  -- should already be evaluated as python native
        np.str_,
        np.timedelta64,
        np.uint8,
        np.uint16,
        np.uint32,
        np.uint64,
        np.void,
    )

    if isinstance(o, np.ndarray):
        return o.tolist()
    elif isinstance(o, numpy_types):        
        return o.item()
    elif isinstance(o, np.float128):
        return o.astype(np.float64).item()
    # elif isinstance(o, np.complex256): -- no python native for np.complex256
    #     return o.astype(np.complex128).item() -- python `complex` class is not json serializable 
    else:
        raise TypeError("{} of type {} is not JSON serializable".format(repr(o), type(o)))




### Inputs ###
# and setup
        
# Load settings fit settings file.
sys.path.append(os.getcwd())
print '  '
print "The name of the settings file is: ", sys.argv[1]
FitSettings = __import__(sys.argv[1])
FitParameters = dir(FitSettings)

# List all the required parameters.
# NOTE: No doubt this will change hence denoted as a list that is worked over. 
# 1. Have option to have list of files rather than sequence. This list needs to be called diff_files.
RequiredList = ['AziBins',
                'Calib_type', 
                'Calib_data', 
                'Calib_param', 
                'Calib_pixels',   #this is only needed for GSAS-II and should be read from image file.
                #'Calib_mask',    #a mask file is not strictly required.
                'Output_NumAziWrite', 
                'Output_directory', 
                'Output_type', 
                'datafile_Basename', 
                'datafile_Ending', 
                'datafile_EndNum', 
                'datafile_NumDigit', 
                'datafile_StartNum', 
                'datafile_directory', 
                'fit_orders']

#check all the required values are present.
all_present = 1
# properties of the data files.
for x in range(len(RequiredList)):
    if RequiredList[x] in FitParameters:
        print 'Got: ', RequiredList[x]
    else:
        print 'The settings file requires a parameter called  \'',RequiredList[x],'\''
        all_present = 0
# exit if all parameters are not present
if all_present == 0:
    sys.exit('The highlighed Settings are missing from the input file. Fitting cannot proceed until they are all present.')



# load the detector functions here.
calib_mod = 'FPF_' + FitSettings.Calib_type + 'Functions'
det = __import__(calib_mod)


# Detector check.
# Check if the detector is required and if so if it is present and recognised.
if 'Calib_detector' in FitParameters:
    det_name = FitSettings.Calib_detector
else:
    det_name = None
FitSettings.Calib_detector = det.DetectorCheck(FitSettings.Calib_data, det_name)

    

# define locally required names.
temporary_data_file = 'PreviousFit_JSON.dat'

backg_type = 'flat'
backg = [[4]]





#### Main code ####


## Load files ##

# calibration parameter file #
parms_dict = det.GetCalibration(FitSettings.Calib_param)


### Sort background out ####
if backg_type == 'coeffs':
    ## coefficients provided, figure out height using values
    tempbackg = []
    bg_order = []
    if isinstance(backg, list):
        for val in backg:
            if not isinstance(val, list):
                tempbackg.append([val])
                bg_order.append(0)
            else:
                tempbackg.append(val)
                bg_order.append((len(val)-1)/2)
    else: 
        tempbackg = [[backg]]
    backg = tempbackg
            
elif backg_type == 'order':
    ## orderes requested, figure out height using best guess
    print 'TBD'
    tempbackg = []
    bg_order = []
    if isinstance(backg, list):
        if len(backg) == 2:
            nfourier = backg[0]
            npoly = backg[1]
        else:
            nfourier = backg[0]
            npoly = 0
    else:
        nfourier = backg
        npoly = 0
    for i in xrange(npoly+1):
        tempbackg_f =  [1 for j in xrange(nfourier+1)]
        tempbackg.append(tempbackgf)
        bg_order.append(nfourier)
    backg = tempbackg

elif backg_type == 'flat':
    backg = [[backg]]



# diffraction patterns #
if not 'diff_files' in locals():
    n_diff_files = FitSettings.datafile_EndNum - FitSettings.datafile_StartNum + 1
    diff_files = []
    for j in range(n_diff_files):
        # make list of diffraction pattern names

        #make number of pattern
        n = str(j+FitSettings.datafile_StartNum).zfill(FitSettings.datafile_NumDigit)
        
        #append diffraction pattern name and directory
        diff_files.append(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n + FitSettings.datafile_Ending)
else:
    n_diff_files = len(diff_files)


  

#get calibration for image set.
#get array for azimuth, two theta and d-space the same size as the image.
azimu    = det.GetAzm(FitSettings.Calib_data, parms_dict, FitSettings.Calib_pixels, FitSettings.Calib_detector)
twotheta = det.GetTth(FitSettings.Calib_data, parms_dict, FitSettings.Calib_pixels, FitSettings.Calib_detector)
dspace   = det.GetDsp(FitSettings.Calib_data, parms_dict, FitSettings.Calib_pixels, FitSettings.Calib_detector)
#FIX ME: Should be able to send diffraction patter, mask and calibration parameters to a function and get back three masked arrays.

# load calibration data/image -- to get intensities for mask.
im     = det.ImportImage(FitSettings.Calib_data)
intens = ma.array(im)

    
# create mask from mask file if present.
# if not make all values valid
if 'Calib_mask' in FitParameters:
    mask_file = FitSettings.Calib_mask
    intens = det.GetMask(mask_file, intens, twotheta, azimu, FitSettings.Calib_data, FitSettings.Calib_pixels)
else:
    #print np.shape(intens), type(intens)
    #print np.shape(intens.mask), type(intens)
    ImMask = np.zeros_like(intens)
    intens = ma.array(intens, mask=ImMask)
    intens = ma.masked_outside(intens,0,np.inf)
    #print np.shape(intens), type(intens)
    #print np.shape(intens.mask), type(intens)
    #stop
    
#apply mask to arrays
azimu    = ma.array(azimu, mask=intens.mask)
twotheta = ma.array(twotheta,mask=intens.mask)
dspace   = ma.array(dspace,mask=intens.mask)



## plot calibration file
if 1:
    fig = plt.figure()
    #plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 1000)  #FIX ME: need variable for vmax.
    plt.scatter(twotheta, azimu, s=intens/10, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 1000)  #FIX ME: need variable for vmax.
    plt.colorbar()
    plt.show()
    plt.close()
    
    quit()
    

# Process the diffraction patterns #

for j in range(n_diff_files):

    print 'Process ',diff_files[j]
    
    # Get diffraction pattern to process.
    im = det.ImportImage(diff_files[j])
    
    #get intensity array and mask it
    intens    = ma.array(im, mask=azimu.mask)
    
    ## plot input file
    if 0:
        fig = plt.figure()
        plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 4000)  #FIX ME: need variable for vmax.
        plt.colorbar()
        plt.show()
        plt.close()
        
        quit()
        

    #get previous fit (if exists)
    # load temporary fit file - if it exists. 
    if os.path.isfile(temporary_data_file):

        # Read JSON data from file
        with open(temporary_data_file) as json_data:
            previous_fit = json.load(json_data)

    
    ## Pass each subpatern to FPF_Fit_Subpattern for fitting in turn.
    #get number of subpatterns to be fitted
    n_subpats = len(FitSettings.fit_orders)
    Fitted_param = []
    
    for i in range(n_subpats):

        tthRange = FitSettings.fit_orders[i]['range'][0]

        # get subsections of data to pass
        subpat       = np.where((twotheta>=tthRange[0])&(twotheta<=tthRange[1]))
        twotheta_sub = twotheta[subpat]
        dspacing_sub = dspace[subpat]
        azimu_sub    = azimu[subpat]
        intens_sub   = intens[subpat]

        if 'previous_fit' in locals():
            params = previous_fit[i]
        else:
            params = []

        if FitSettings.fit_orders:
            orders = FitSettings.fit_orders[i]
            orders['AziBins'] = FitSettings.AziBins
        #print orders

        #fit the subpattern
        Fitted_param.append(FitSubpattern([twotheta_sub, dspacing_sub, parms_dict], azimu_sub, intens_sub, orders, params, DetFuncs=calib_mod))
        #Fitted_param.append(FitSubpattern([twotheta_sub, dspacing_sub, parms_dict['conversion_constant']], azimu_sub, intens_sub, orders, params))

    #stop

    #write output files
    
    # store all the fit information as a JSON file. 
    #filename, file_extension = os.path.splitext(diff_files[j])
    filename = os.path.splitext(os.path.basename(diff_files[j]))[0]
    filename = filename+'.json'
    print filename
    with open(filename, 'w') as TempFile:
        # Write a JSON string into the file.
        json_string = json.dump(Fitted_param, TempFile, sort_keys=True, indent=2, default=json_numpy_serialzer)
    
    # FIX ME: Need to add optins for choice of output.
    wr.WriteMultiFit(diff_files[j], Fitted_param, FitSettings.Output_NumAziWrite, parms_dict['conversion_constant'])

    # #Write the fits to a temporary file
    # TempFilename = open(temporary_data_file, 'w')
    # # Write a JSON string into the file.
    # json_string = json.dumps(Fitted_param, TempFilename, sort_keys=True, indent=4, default=json_numpy_serialzer)
    # #json_string = json.dumps(Fitted_param, TempFilename, sort_keys=True, indent=4, separators=(',', ': '))
    # TempFilename.write(json_string)

    #print json_string
    with open(temporary_data_file, 'w') as TempFile:

        # Write a JSON string into the file.
        json_string = json.dump(Fitted_param, TempFile, sort_keys=True, indent=2, default=json_numpy_serialzer)
