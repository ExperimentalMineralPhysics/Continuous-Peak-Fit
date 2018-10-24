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

import FPF_PeakFourierFunctions as ff

# load as det (detector) so that can readilt replace the GSASII functions with e.g. diaoptas without confusing names in this script
import FPF_GSASIIFunctions as det


### Inputs ###
# and setup

calib_file = 'CeO2_Pil207_E30_2Nov2016_001.imctrl'
pix = 172 #microns ##SAH: edit 172 microns not 17.2

opt = 2

if opt==1:
    diff_file = 'CeO2_Pil207_E30_2Nov2016_001.tif'
    mask_file = 'Calibration.immask'
    tthRange = [7.6,7.8]

    total = 25

    backg = 10.

    h_order = 0  ##  no formal limit
    w_order = 5  ## no formal limit
    d0_order = 0 ##  probably limited to 2 -- generalisation might work better for radio data.

elif opt==2:
    diff_file = 'BCC1_2GPa_10s_001_00001.tif'
    mask_file = 'Diffraction.immask'
    tthRange = [16.6,17.1]

    total = 72
    not_use = 5
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

im = Image.open(diff_file) ##always tiff?
## SAH: No. could be propreity formats. Possible formats for GSAS-II diffracion import are: tif, ADSC, cbf, png, edf, ge, hdf5, mar, rigaku, sfrm, 'gsas-ii image file' and guess.
## SAH: GSAS-II import files are located at https://subversion.xray.aps.anl.gov/pyGSAS/trunk/imports/ with the file names 'G2img_*.py'.
#im.show()
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
#plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
#plt.colorbar()
#plt.show()




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


#fig = plt.figure()
#plt.scatter(twotheta, azimu, s=4, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
##plt.scatter(dspace.flatten()[tthchunk],azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
#
#plt.colorbar()
#plt.show()
#
#plt.close()


## Split data into chunks in azimuth, replace azimuths with avg.
## doesn't have to be 100

print '\nGroup fitting in azimuth...\n'


azmax = azimu.max()
azmin = azimu.min()
binsize = (azmax-azmin)/total
chunks = []
azichunks = []
for i in xrange(total):
    end = (i+1)*binsize
    start = i*binsize
    temptth = twotheta.flatten()
    tthchunk = np.where((temptth>=tthRange[0])&(temptth<=tthRange[1]))
    tempazi = azimu.flatten()
    azichunk = np.where((tempazi>start)&(tempazi<=end))
    chunkels = np.intersect1d(tthchunk,azichunk)

    azichunks.append(((end-start)/2)+start)
    chunks.append(chunkels)


##### Final output list of azimuths with corresponding twotheta_0,h,w

wavelength=parms_dict['wavelength']
newd0=[]
newd0Err=[]
newHall=[]
newHallErr=[]
newWall=[]
newWallErr=[]
newBGall=[]
newBGallErr=[]

for j in range(len(chunks)):


    # print '\nFitting azimuth chunk....\n'

    # print intens.flatten()[chunks[j]].dtype
    hguess_I=intens.flatten()[chunks[j]].argmax()
    #print hguess_I
    hguess=intens.flatten()[chunks[j]][hguess_I]-backg[0]
    #print hguess
    dguess=dspace.flatten()[chunks[j]][hguess_I]
    bgguess = backg
    wguess = 0.005
    print hguess_I,hguess,dguess,backg,'hguess'
    po,pc= ff.singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[j]],d0s=dguess,heights=hguess,widths=wguess,wavelength=np.float(parms_dict['wavelength']),bg=bgguess)
    #print po
    AllErr = np.sqrt(np.abs(np.diag(pc)))
    #print AllErr

    #print newd0Err,'newd0Err'
    newd0.append(po[0])
    newd0Err.append(AllErr[0])
    newHall.append(po[1])
    newHallErr.append(AllErr[1])
    newWall.append(po[2])
    newWallErr.append(AllErr[2])
    newBGall.append(po[3])
    newBGallErr.append(AllErr[3])

    print po, 'p0'

    '''
    newd0.append(po[0])
    newd0Err.append(np.sqrt(pc[0,0]))
    newHall.append(po[1])
    newHallErr.append(np.sqrt(pc[1,1]))
    newWall.append(po[2])
    newWallErr.append(np.sqrt(pc[2,2]))
    newBGall.append(po[3])
    newBGallErr.append(np.sqrt(pc[3,3]))
    '''

    # #plot the fits.
    # asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
    # asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
    # asdf2 = intens.flatten()[chunks[j]][asdf_ord]
    # asdf3 = singleInt(asdf1,tuple(po))
    # #asdf3 = singleInt(asdf1,po[0],po[1],po[2],po[3])


    # plt.plot(asdf1, asdf2,'.')
    # plt.plot(asdf1, asdf3, '-')
    # plt.show()




### Feed each d_0,h,w into fourier expansion function to get fit for fourier component
### parameters as output.

# print '\n Doing d0 fourier fit...', d0_order

#print 'newd0Err', newd0Err

dfour = ff.Fourier_fit(np.array(azichunks),np.array(newd0),terms=d0_order, errs=np.array(newd0Err))
# print '\n Doing h fourier fit...'
hfour = ff.Fourier_fit(np.array(azichunks),np.array(newHall),terms=h_order, errs=np.array(newHallErr))
wfour = ff.Fourier_fit(np.array(azichunks),np.array(newWall),terms=w_order, errs=np.array(newWallErr))

'''
print 'shapes'
print dfour[0].shape
print hfour[0].shape
print wfour[0].shape
#print 'hfour'
#print hfour
#print 'wfour'
#print wfour
'''

#plot output of fourier fits....

# fig = plt.figure()
# ax = fig.add_subplot(3,1,1)
# plt.subplot(311)
# plt.plot(azichunks,newd0, 'bo')
# plt.title('D-spacing')
# # print len(azichunks), 'len azichunks'
# print dfour[0]
# testout = Fourier_expand(np.array(azichunks), dfour[0])
# print testout, 'testout'
# #print stop
# plt.plot(azichunks,Fourier_expand(np.array(azichunks), dfour[0]), 'r-')
# plt.subplot(312)
# plt.plot(azichunks,newHall, 'bo')
# plt.title('height')
# plt.plot(azichunks,Fourier_expand(np.array(azichunks), hfour[0]), 'r-')
# plt.subplot(313)
# plt.plot(azichunks,newWall, 'bo')
# plt.title('width')
# plt.plot(azichunks,Fourier_expand(np.array(azichunks), np.array(wfour[0])), 'r-')
# plt.show()

# plt.close()

# print stop


### Feed each d,h,w into versions of main equation to fit only or one at a time
### including the new paramters from the Fourier fits above as initial inputs
## Slice arrays for twothetarange

print '\nRe-fitting for d, h, w separately...\n'


temptth = twotheta.flatten()
tthchunk = np.where((temptth>=tthRange[0])&(temptth<=tthRange[1]))
#print 'check point'

## Fit for d0 with existing h and w

heights  = ff.Fourier_expand(azimu.flatten()[tthchunk],hfour[0])
widths   = ff.Fourier_expand(azimu.flatten()[tthchunk],wfour[0])
#print dfour[0]
newdfour = ff.dfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],dspace.flatten()[tthchunk],dfour[0],heights,widths,wavelength,bg=backg)
print 'Old d coefficients: ', dfour[0]
print 'New d coefficients: ', newdfour[0]

## Fit for h with new d and existing w

newd0s = ff.Fourier_expand(azimu.flatten()[tthchunk],newdfour[0])
widths = ff.Fourier_expand(azimu.flatten()[tthchunk],wfour[0])
#print hfour[0]
newhfour = ff.hfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],newd0s,widths,hfour[0],wavelength,bg=backg)
print 'Old h coefficients: ', hfour[0]
print 'New h coefficients: ', newhfour[0]

## Fit for w with new d and h

newd0s     = ff.Fourier_expand(azimu.flatten()[tthchunk],newdfour[0])
newheights = ff.Fourier_expand(azimu.flatten()[tthchunk],newhfour[0])
#print wfour[0]
newwfour   = ff.wfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],newd0s,newheights,wfour[0],wavelength,bg=backg)
print 'Old w coefficients: ', wfour[0]
print 'New w coefficients: ', newwfour[0]


### Refit through full equation with all data for d,h,w independently

print '\nFinal fit solving for all parms...\n'

newwidths = ff.Fourier_expand(azimu.flatten()[tthchunk],newwfour[0])

Finalparms = ff.Allfit(intens.flatten()[tthchunk],twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],newdfour[0],newhfour[0],newwfour[0],wavelength,bg=backg)

starth = len(newdfour[0])
startw = len(newdfour[0])+len(newhfour[0])
end = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])

fin_d = Finalparms[0][0:starth]
fin_h = Finalparms[0][starth:startw]
fin_w = Finalparms[0][startw:end]


#print len(fin_d)

print 'Final d0 coefficients: ', Finalparms[0][0:starth]
print 'Final h coefficients: ', Finalparms[0][starth:startw]
print 'Final w coefficients: ', Finalparms[0][startw:end]





### Plot results to check

print '\nPlotting results for fit...\n'

tthchunk = np.where((twotheta.flatten()>=tthRange[0])&(twotheta.flatten()<=tthRange[1]))
# print Finalparms[0]
fit_intens = ff.Allchange((twotheta.flatten()[tthchunk],azimu.flatten()[tthchunk],np.array([len(fin_d),len(fin_h),len(fin_w)]), wavelength, backg), *Finalparms[0].squeeze())
fullfit_intens = np.ones((intens.shape[0],intens.shape[1]),float)
oldshape = fullfit_intens.shape
fullfit_intens = fullfit_intens.flatten()
fullfit_intens[tthchunk] = fit_intens
#print fullfit_intens.max()
#print fullfit_intens.shape
#print oldshape
fullfit_intens = np.reshape(fullfit_intens,(oldshape[0],oldshape[1]))

ValMax1 = np.max(intens.flatten()[tthchunk])
ValMax2 = np.max(fullfit_intens.flatten()[tthchunk])
ValMax3 = np.max((intens.flatten()[tthchunk]-fullfit_intens.flatten()[tthchunk]))
ValMax  = 40#np.max([ValMax1, ValMax2, ValMax3])
ValMin1 = np.min(intens.flatten()[tthchunk])
ValMin2 = np.min(fullfit_intens.flatten()[tthchunk])
ValMin3 = np.min((intens.flatten()[tthchunk]-fullfit_intens.flatten()[tthchunk]))
ValMin  = -3#np.min([ValMin1, ValMin2, ValMin3])

#print fullfit_intens.shape
#print twotheta.shape
#print azimu.shape
#print intens.shape
#print tthchunk
#print twotheta[0,]
#print np.where((twotheta[0]>=tthRange[0])&(twotheta[0]<=tthRange[1]))

if 1:
    fig = plt.figure()
    ax = fig.add_subplot(1,3,1)
    plt.subplot(131)
    #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
    plt.scatter(dspace.flatten()[tthchunk],azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
    plt.subplot(132)
    #plt.scatter(twotheta, azimu, s=1, c=np.log(fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
    plt.scatter(dspace.flatten()[tthchunk], azimu.flatten()[tthchunk], s=4, c=(fullfit_intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
    plt.colorbar()
    plt.subplot(133)
    #plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
    plt.scatter(dspace.flatten()[tthchunk], azimu.flatten()[tthchunk], s=4, c=(intens.flatten()[tthchunk]-fullfit_intens.flatten()[tthchunk]), edgecolors='none', cmap=plt.cm.jet)
    plt.colorbar()
    plt.show()

    plt.close()



    # print '\nFitting azimuth chunk....\n'

  #  # print intens.flatten()[chunks[j]].dtype
  #  hguess_I=intens.flatten()[chunks[j]].argmax()
  #  hguess=intens.flatten()[chunks[j]][hguess_I]
  #  dguess=dspace.flatten()[chunks[j]][hguess_I]
  #  #print hguess_I,hguess,dguess
  #  po,pc=singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[#j]],d0s=dguess,heights=hguess,widths=0.005,wavelength=np.float(parms_dict['wavelength']),bg=10.)
  #  print po,p#c
  #  newd0.append(po[0])
  #  newHall.append(po[1])
  #  newWall.append(po[2])

  #  #plot the fits.
  #  asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
  #  asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
  #  asdf2 = intens.flatten()[chunks[j]][asdf_ord]
  #  asdf3 = singleInt(asdf1,po[0],po[1],po[2])

  #  plt.plot(asdf1, asdf2,'x--')
   # plt.plot(asdf1, asdf3, '-')
  #  plt.show()












## Specify ranges to search/starting parameters ##



#Input from GUI/file?
## SAH: both. I think I would have two functions.
## SAH: A GUI setup function and a wrapper that calls this function with a parameter file or if it does not exist then calls the GUI function.

#Load plot to screen and choose regions to fit? ## SAH: yes.
#can then slice arrays and calcualte intial guesses from peaks and avgs.



## Setup/Choose function to fit ##


## make arrays. -- not line arrays n by 2. possibly better as a line array becuase then we are not carrying around an unused number.
#heights = ones(2,h_order) ## better to be intellegent guesses
#widths = ones(2,w_order)  ## better to be intellegent guesses
#d0s = ones(2,d0_order)    ## better to be intellegent guesses


## Now have starting model


## refit starting model via full intensity equation but hold other parameters fixed e.g. fit for d, not h,w,
## do for each in turn i.e. if strat with d, then h, then w


## then continue fitting with full equation, varying all parms.





'''
#start_model
Segregate into chunks in azimuth
get avg. azimuth
flatten intens and twotheta arrays to pass to fit
pass to gauss fit to get first approx of 2th_0,h,w with nterms =0 for w,h,2th for each block

#start_fourier
take array of 2th_0,h,w over azimu blocks
choose nterms for each for now
fit each with Fourier expansion of nterms

#main_fit
use all data
take fourier expansion with same nterms and input parms from start_fourier
send to allfit to allow all parms to fit

#inside main fit
redo 2th_0/d_0,h,w fit iteratively
then input to allfit to fit parms in gauss function
increase nterms?

generally think will have to pass flattened azimu,2th and intens arrays.

'''
