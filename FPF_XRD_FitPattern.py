# /usr/bin/python

import numpy as np
import numpy.ma as ma
import copy, os, sys
import json
#from PIL import Image
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.path as mlp
from scipy.optimize import curve_fit
np.set_printoptions(threshold='nan')

from FPF_XRD_FitSubpattern import FitSubpattern
import FPF_PeakFourierFunctions as ff
import FPF_WriteMultiFit as wr

# load as det (detector) so that can readilt replace the GSASII functions with e.g. diaoptas without confusing names in this script
#calib_mod = 'FPF_DioptasFunctions'
#det = __import__(calib_mod)
#import calib_mod as det



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
sys.path.append(os.getcwd())
print '  '
print "The name of the settings file is: ", sys.argv[1]
FitSettings = __import__(sys.argv[1])
#from BCC1_input import *
FitParameters = dir(FitSettings)

#list all the required parameters..
# NOTE: No doubt this will change hence denoted as a list that is worked over. 
# 1. Have option to have list of files rather than sequence. This list needs to be called diff_files.
RequiredList = ['AziBins',
                'Calib_type', 
                'Calib_data', 
                'Calib_param', 
                'Calib_pixels',   #this is only needed for GSAS-II and should be read from image file.
                #'Calib_mask',    #a mask file is not stricktly required.
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




# load the detetor functions here. 
calib_mod = 'FPF_' + FitSettings.Calib_type + 'Functions'
det = __import__(calib_mod)
#import calib_mod as det




# define locally required names.
temporary_data_file = 'PreviousFit_JSON.dat'

backg_type = 'flat'
backg = [[4]]



# tthRange -- is now defined in 'fit_orders' 



# calib_file = './Example1-Fe/CeO2_Pil207_E30_2Nov2016_001.imctrl'
# calib_file = './Example1-Fe/Calibration2018Oct.imctrl'
# pix = 172 #microns ##SAH: edit 172 microns not 17.2

# opt = 2

# if opt==1:
#     diff_files = './Example1-Fe/CeO2_Pil207_E30_2Nov2016_001.tif'
#     mask_file = './Example1-Fe/Diffraction.immask'
#     tthRange = [7.6,7.8]

#     total = 25
    
#     backg_type = 'flat'
#     backg = 10.

#     h_order = 0  ##  no formal limit
#     w_order = 5  ## no formal limit
#     d0_order = 0 ##  probably limited to 2 -- generalisation might work better for radio data.

# elif opt==2:
#     diff_files = ['./Example1-Fe/BCC1_2GPa_10s_001_00001.tif','./Example1-Fe/BCC1_2GPa_10s_001_00002.tif']
#     mask_file = './Example1-Fe/Diffraction.immask'
#     tthRange = [[11.7, 12.03], [16.6,17.1], [20.4, 20.9]]
#     #tthRange = [[16.6,17.1] ]

#     AziBins = 45
#     backg_type = 'coeffs'
#     backg = [[4.,1.,1.]] #have to make sure this number is the right size for the bg order.

#     h_order  = 12  ##  no formal limit
#     w_order  = 1   ## no formal limit
#     d0_order = 2 ##  probably limited to 2 -- generalisation might work better for radio data.
#     bg_order = 0 ## no limit.

#     NumAziWrite = 180


#     fit_orders = [
#           {
#             "background": [1], 
#             "peak": [{
#                 "d-space": 6, 
#                 "height": 20, 
#                 "profile": 0, 
#                 "width": 8
#               }], 
#             "range": [[11.70, 12.03]]
#           }, 
#           {
#             "background": [1], 
#             "peak": [{
#                 "d-space": 2, 
#                 "height": 12, 
#                 "profile": 0, 
#                 "width": 1
#               }], 
#             "range": [[16.6,17.1]]
#           },  
#           {
#             "background": [1], 
#             "peak": [{
#                 "d-space": 2, 
#                 "height": 12, 
#                 "profile": 0, 
#                 "width": 1
#               }], 
#             "range": [[20.4, 20.9]]
#           }, 
#         ]




# GUI needed to set inputs? ###
# SAH: A GUI would be good for the initial inputs -- most geoscientist types find it is easier to select areas with a mouse than the keyboard.
# SAH: The GUI should save an output file that is then used as the input file for batch processing lots of files.



#### Main code ####


## Load files ##

# calibration parameter file #
parms_dict = det.GetCalibration(FitSettings.Calib_param)

#parms_file = open(FitSettings.Calib_param,'rb')
#parms_dict = det.GetCalibration(parms_file)

## may be a mask file?
## SAH: yes mask files -- they will either be a list of regions (x,y) or a binary image.
## SAH: I will generate you a mask file and send it to you.




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





# diffraction pattern #

if not 'diff_files' in locals():

    n_diff_files = FitSettings.datafile_EndNum - FitSettings.datafile_StartNum + 1
    diff_files = []
    for j in range(n_diff_files):
        # make list of diffraction pattern names

        #make number of pattern
        n = str(j+1).zfill(FitSettings.datafile_NumDigit)

        #append diffraction pattern name and directory
        diff_files.append(FitSettings.datafile_directory + os.sep + FitSettings.datafile_Basename + n + FitSettings.datafile_Ending)

else:
    n_diff_files = len(diff_files)


#for now redfine previously used parameters from input file.
pix = FitSettings.Calib_pixels


#get calibration for image set.
#get array for two theta and d-space the same size as the image.
azimu    = det.GetAzm(FitSettings.Calib_data, parms_dict, pix)
twotheta = det.GetTth(FitSettings.Calib_data, parms_dict, pix)
dspace   = det.GetDsp(FitSettings.Calib_data, parms_dict, pix)

# load calibration data/image -- to get intensities for mask.
im     = det.ImportImage(FitSettings.Calib_data)
intens = ma.array(im)

# create mask from mask file if present.
# if not make all values valid
if 'Calib_mask' in FitParameters:
    
    #for now redfine previously used parameters from input file.
    mask_file = FitSettings.Calib_mask

    intens = det.GetMask(mask_file, intens, twotheta, azimu, FitSettings.Calib_data, pix)
else:
    
#    #intens = ma.masked_array(intens, mask=np.zeros(intens.shape))
    intens = ma.masked_outside(intens,0,np.inf)
    
    
#apply mask to arrays
azimu    = ma.array(azimu, mask=intens.mask)
twotheta = ma.array(twotheta,mask=azimu.mask)
dspace   = ma.array(dspace,mask=azimu.mask)




# Process the diffraction patterns #

for j in range(n_diff_files):

    print 'Process ',diff_files[j]
    
    im = det.ImportImage(diff_files[j])
    
    #get intensity array and mask it
    intens    = ma.array(im, mask=azimu.mask)



    #mask negative intensities
    #intens = ma.masked_outside(intens,int(0),np.inf)
    #intens   = ma.array(intens, mask=(intens.mask*azimu.mask))
    #azimu    = ma.array(azimu, mask=intens.mask)
    #twotheta = ma.array(twotheta,mask=intens.mask)
    #dspace   = ma.array(dspace,mask=intens.mask)
    
    
    #im = Image.open(diff_file) ##always tiff?
    ## SAH: No. could be propreity formats. Possible formats for GSAS-II diffracion import are: tif, ADSC, cbf, png, edf, ge, hdf5, mar, rigaku, sfrm, 'gsas-ii image file' and guess.
    ## SAH: GSAS-II import files are located at https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/ with the file names 'G2img_*.py'.
    # im.show()
    #imarray = np.array(im)
    #print imarray.shape,imarray.max()

    ## Convert data ##

    #The bits of the calibration that we need to convert X,Y to 2theta, azimuth or d,azimuth  are: Wavelength (Angstroms), distance (mm), center (X,Y, mm), tilt (degrees), rotation (degrees) and DetDepth (called penetration in screen shot).
    #x = np.arange(imarray.shape[1])
    #y = np.arange(imarray.shape[0])[0:imarray.shape[1]]

    #gd = np.mgrid[0:imarray.shape[0],0:imarray.shape[1]]       ##SAH: edit
    #print gd.shape
    #y = gd[0,:,:]+1       ##SAH: edit
    #x = gd[1,:,:]+1       ##SAH: edit
    #print y.shape, x.shape
    #print x[2]
    #y = (y) * pix / 1e3   ##SAH: edit
    #x = (x) * pix / 1e3   ##SAH: edit
    #print x             ##SAH: edit
    #print y             ##SAH: edit
    
    
    '''
    ## SAH: would linspace+meshgrid be a better way of making the arrays?
    ## DMF: Could do it this way...
    xn = np.linspace(0, imarray.shape[0],num=imarray.shape[0],endpoint=False,dtype=int)#nx)
    yn = np.linspace(0, imarray.shape[1],num=imarray.shape[1],endpoint=False,dtype=int)#ny)
    xv, yv = np.meshgrid(xn, yn)

    #print xv.shape,yv.shape,xv[1],yv[1]

    '''

    
    ## Convert to twotheta and azimuth arrays ##
    #twotheta = det.GetTth(x,y,parms_dict)
    #azimu    = det.GetAzm(x,y,parms_dict)
    #dspace   = det.GetDsp(x,y,parms_dict)
    #intens   = imarray
    
    #azimu = (azimu + 180) % 360


    # checks via prints and plots #
    #print twotheta.shape,azimu.shape,dspace.shape,twotheta[0]
    #colors = ("red", "green", "blue")
    #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
    # plt.scatter(twotheta, azimu, s=1, c=intens, edgecolors='none', cmap=plt.cm.jet)
    # plt.colorbar()
    # plt.show()
    # plt.close()


    ## selected limits via gui - 2theta_min, 2theta_max, n_peaks (each peak 2theta min+max)
    ## for each peak....
    

    # create mask from mask file if present.
    # if not make all values valid
    
    #apply mask here incase there is an intensity threshold in the mask.
    #if 'Calib_mask' in FitParameters:
    #            
    #    #for now redfine previously used parameters from input file.
    #    #pix = FitSettings.Calib_pixels
    #    mask_file = FitSettings.Calib_mask
    #
    #    intens = det.GetMask(mask_file, intens, twotheta, azimu, FitSettings.Calib_file, pix)
    #else:
    #intens = ma.masked_array(intens, mask=np.zeros(intens.shape))
    #intens = ma.masked_outside(intens,int(0),np.inf)
    
    #make sure there are no 0 or negative pixels.   
    
    
        
    #apply mask to other arrays
    
    intens   = ma.array(intens, mask=(azimu.mask))
    #azimu    = ma.array(azimu, mask=intens.mask)
    #twotheta = ma.array(twotheta,mask=intens.mask)
    #dspace   = ma.array(dspace,mask=intens.mask)
    
    


    ## plot input file
    if 0:
        fig = plt.figure()
        plt.scatter(twotheta, azimu, s=4, c=(intens), edgecolors='none', cmap=plt.cm.jet, vmin = 0, vmax = 300)
        #plt.scatter(dspace.flatten()[tthchunk],azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
        plt.colorbar()
        plt.show()
        plt.close()

        quit()

    ## Pass each subpatern to FPF_Fit_Subpattern for fitting in turn.

    ##set up things to be fed to subpattern fit.
    #if 0:
        # pass previous coefficients to the fit.
    #elif:
        #something here

    '''
    print np.size(tthRange,0)
    print np.size(tthRange,1)
    print tthRange
    print type(tthRange)
    '''

    #get previous fit (if exists)
    # load temporary fit file - if it exists. 
    if os.path.isfile(temporary_data_file):

        # Read JSON data from file
        with open(temporary_data_file) as json_data:
            previous_fit = json.load(json_data)
        #print previous_fits
        #print previous_fits[0]['peak'][0]['profile']


    #n_subpats = np.size(tthRange,0)
    n_subpats = len(FitSettings.fit_orders)
    Fitted_param = []

    for i in range(n_subpats):

        tthRange = FitSettings.fit_orders[i]['range'][0]
        #print tthRange

        # get subsections of data to pass
        temptth = twotheta.flatten()
        #subpat       = np.where((twotheta>=tthRange[i][0])&(twotheta<=tthRange[i][1]))
        subpat       = np.where((twotheta>=tthRange[0])&(twotheta<=tthRange[1]))
        twotheta_sub =  twotheta[subpat]
        dspacing_sub =  dspace[subpat]
        azimu_sub    =  azimu[subpat]
        intens_sub   =  intens[subpat]

        if 'previous_fit' in locals():
            params = previous_fit[i]
        else:
            params = []


        if FitSettings.fit_orders:
            orders = FitSettings.fit_orders[i]
            orders['AziBins'] = FitSettings.AziBins

        #print orders
        #params = []
        #if backg:
        #    params = {'background': backg, 'backgrnd_type': backg_type}

        #get fourier order (if required). The last parameter is the number of bins to use in initial fitting 
        #orders = [d0_order, h_order, w_order, bg_order, total]

        #fit the subpattern

        #print parms_dict
        Fitted_param.append(FitSubpattern([twotheta_sub, dspacing_sub, parms_dict['wavelength']], azimu_sub, intens_sub, orders, params))

        #print Fitted_param

    stop

    #write output files


    wr.WriteMultiFit(diff_files[j], Fitted_param, FitSettings.Output_NumAziWrite, parms_dict['wavelength'])

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
