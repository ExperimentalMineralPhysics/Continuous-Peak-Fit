# /usr/bin/python

import numpy as np
import numpy.ma as ma
import copy, os, sys
from PIL import Image
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
import FPF_GSASIIFunctions as det


### Inputs ###
# and setup

#calib_file = '/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/python/CeO2_Pil207_E30_2Nov2016_001.imctrl'
#calib_file = '/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/Calibration2018Oct.imctrl'
calib_file = '/local/shunt/CeO2_Pil207_E30_2Nov2016_001.imctrl'
calib_file = '/local/shunt/test.imctrl'
pix = 172 #microns ##SAH: edit 172 microns not 17.2

opt = 2

if opt==1:
    diff_files = '/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/python/CeO2_Pil207_E30_2Nov2016_001.tif'
    mask_file = '/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/python/Diffraction.immask'
    tthRange = [7.6,7.8]

    total = 25

    backg = 10.

    h_order = 0  ##  no formal limit
    w_order = 5  ## no formal limit
    d0_order = 0 ##  probably limited to 2 -- generalisation might work better for radio data.

elif opt==2:
    # diff_files = ['/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/python/BCC1_2GPa_10s_001_00001.tif','/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/python/BCC1_2GPa_10s_001_00001.tif']
    # mask_file = '/Users/simon/Documents/WORK/Experiment_Analysis/DebyeFitting/python/Diffraction.immask'
    diff_files = ['/local/shunt/BCC1_2GPa_10s_001_00001.tif','/local/shunt/BCC1_2GPa_10s_001_00001.tif']
    mask_file = '/local/shunt/Diffraction.immask'
    tthRange = [[11.7, 12.08], [16.6,17.1], [20.4, 20.9]]

    total = 72
    backg = [4.] #have to make sure this number is the right size for the bg order.

    h_order  = 12  ##  no formal limit
    w_order  = 0   ## no formal limit
    d0_order = 2 ##  probably limited to 2 -- generalisation might work better for radio data.
    bg_order = 0 ## no limit.

# GUI needed to set inputs? ###
# SAH: A GUI would be good for the initial inputs -- most geoscientist types find it is easier to select areas with a mouse than the keyboard.
# SAH: The GUI should save an output file that is then used as the input file for batch processing lots of files.




#### Main code ####



## Load files ##

# calibration parameter file #
parms_file = open(calib_file,'rb')
parms_dict = det.load_inputs_file(parms_file)
#print parms_dict
## may be a mask file?
## SAH: yes mask files -- they will either be a list of regions (x,y) or a binary image.
## SAH: I will generate you a mask file and send it to you.


# diffraction pattern #


n_diff_files = len(diff_files)

for j in range(n_diff_files):


    im = det.ImportImage(diff_files[j])


    #im = Image.open(diff_file) ##always tiff?
    ## SAH: No. could be propreity formats. Possible formats for GSAS-II diffracion import are: tif, ADSC, cbf, png, edf, ge, hdf5, mar, rigaku, sfrm, 'gsas-ii image file' and guess.
    ## SAH: GSAS-II import files are located at https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/ with the file names 'G2img_*.py'.
    # im.show()
    imarray = np.array(im)
    #print imarray.shape,imarray.max()

    ## Convert data ##

    #The bits of the calibration that we need to convert X,Y to 2theta, azimuth or d,azimuth  are: Wavelength (Angstroms), distance (mm), center (X,Y, mm), tilt (degrees), rotation (degrees) and DetDepth (called penetration in screen shot).
    #x = np.arange(imarray.shape[1])
    #y = np.arange(imarray.shape[0])[0:imarray.shape[1]]

    gd = np.mgrid[0:imarray.shape[0],0:imarray.shape[1]]       ##SAH: edit
    #print gd.shape
    y = gd[0,:,:]       ##SAH: edit
    x = gd[1,:,:]       ##SAH: edit
    #print y.shape, x.shape
    #print x[2]
    y = y * pix / 1e3   ##SAH: edit
    x = x * pix / 1e3   ##SAH: edit
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

    twotheta = det.GetTth(x,y,parms_dict)
    azimu    = det.GetAzm(x,y,parms_dict)
    dspace   = det.GetDsp(x,y,parms_dict)
    intens   = imarray




    # print twotheta.shape,azimu.shape,dspace.shape,intens.shape

    #checks#
    # print twotheta.shape,azimu.shape,dspace.shape,twotheta[0]
    colors = ("red", "green", "blue")
    #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
    # plt.scatter(twotheta, azimu, s=1, c=intens, edgecolors='none', cmap=plt.cm.jet)
    # plt.colorbar()
    # plt.show()


    # plt.close()


    ## selected limits via gui - 2theta_min, 2theta_max, n_peaks (each peak 2theta min+max)
    ## for each peak....


    ## Example introduction of mask via intensity cut. Mask created
    ## for intens array then applied to all arrays.
    # intens = ma.asarray(intens) 
    # intens = ma.masked_less(intens,0)
    # #print azimu[~intens.mask]
    # #print intens.shape
    # #newazimu = ma.asarray(azimu)
    # azimu = ma.array(azimu, mask=intens.mask)
    # #print azimu.shape
    # twotheta = ma.array(twotheta,mask=intens.mask)
    # dspace = ma.array(dspace,mask=intens.mask)



    # create mask from GSAS-II mask file.
    intens = det.CreateGSASIIMask(mask_file, intens, gd.shape, twotheta, azimu, y, x)

    #apply mask to other arrays
    azimu    = ma.array(azimu, mask=intens.mask)
    twotheta = ma.array(twotheta,mask=intens.mask)
    dspace   = ma.array(dspace,mask=intens.mask)




    # plot input file


    # fig = plt.figure()
    # plt.scatter(twotheta, azimu, s=4, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
    # #plt.scatter(dspace.flatten()[tthchunk],azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)

    # plt.colorbar()
    # plt.show()

    # plt.close()



    # Pass each subpatern to FPF_Fit_Subpattern for fitting in turn.

    ##set up things to be fed to subpattern fit.
    #if 0:
        # pass previous coefficients to the fit.
    #elif:
        #something here

    print np.size(tthRange,0)
    print np.size(tthRange,1)
    print tthRange
    print type(tthRange)

    n_subpats = np.size(tthRange,0)
    Fitted_param = []

    for i in range(n_subpats):

        # get subsections of data to pass
        temptth = twotheta.flatten()
        subpat       = np.where((twotheta>=tthRange[i][0])&(twotheta<=tthRange[i][1]))
        twotheta_sub =  twotheta[subpat]
        dspacing_sub =  dspace[subpat]
        azimu_sub    =  azimu[subpat]
        intens_sub   =  intens[subpat]

        #get previous fit (if exists)
        params = []
        if backg:
            params = {'background': backg}

        #get fourier order (if required). The last parameter is the number of bins to use in initial fitting 
        orders = [d0_order, h_order, w_order, bg_order, total]

        #fit the subpattern

        Fitted_param.append(FitSubpattern([twotheta_sub, dspacing_sub, parms_dict['wavelength']], azimu_sub, intens_sub, orders, params))

        print Fitted_param



    #write output files


    wr.WriteMultiFit(diff_files[j], Fitted_param)





