# /usr/bin/python

# FPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier order

# The script should not need to know what the limits are - it should only see the data that it needs to fit. 
# The two theta 

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
#import FPF_GSASIIFunctions as det



def FitSubpattern(TwoThetaAndDspacings, azimu, intens, orders=None, params=None):


    #print np.shape(TwoThetaAndDspacings)[0]


    # sort inputs

    # two theta and dspacing are both 'x' dimensions. therefore import them as single array.
    # split them up here if needed. 
    if np.shape(TwoThetaAndDspacings)[0] == 3:
        twotheta = TwoThetaAndDspacings[0]
        dspace = TwoThetaAndDspacings[1]
        wavelength = TwoThetaAndDspacings[2]
    else:
        twotheta = TwoThetaAndDspacings
        #dspacing = twotheta


    # if no parameters
    #     get orders of arrays
    # else if no orders
    #     use parameters as given
    # else
    #     change the size of the parameter arrays to match orders

    if params:
        backg = params['background']
        backg_type =params['backgrnd_type']
        lenbg=[]
        singlelenbg=[]
        for val in backg:
            lenbg.append(len(val))
            singlelenbg.append(1)
        lenbg=np.array(lenbg)
        singlelenbg=np.array(singlelenbg)


    if orders:
        d0_order, h_order, w_order, bg_order, total = orders


    



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
        #temptth = twotheta.flatten()
        #tthchunk = np.where((temptth>=tthRange[0])&(temptth<=tthRange[1]))
        tempazi = azimu.flatten()
        azichunk = np.where((tempazi>start)&(tempazi<=end))
        #chunkels = np.intersect1d(tthchunk,azichunk)

        azichunks.append(((end-start)/2)+start)
        #chunks.append(chunkels)
        chunks.append(azichunk)


    ##### Final output list of azimuths with corresponding twotheta_0,h,w

    #wavelength=parms_dict['wavelength']
    newd0=[]
    newd0Err=[]
    newHall=[]
    newHallErr=[]
    newWall=[]
    newWallErr=[]
    newBGall=[]
    newBGallErr=[]
    newAziChunks=[]

    for j in range(len(chunks)):


        # print '\nFitting azimuth chunk....\n'

        #only compute the fit and save it if there are 'enough' data
        if np.sum(~intens.flatten()[chunks[j]].mask) >= 20: #FIX ME: this switch is.a crude way of finding if there are less than 20 data in the bin. It should also be a variable. 

            ## background estimates ##
            if backg_type == 'coeffs':
                backg_base = backg[0][0]
            elif backg_type == 'order':
                print 'uh oh!'
            elif backg_type == 'flat':
                backg_base = backg[0][0]
            hguess_I=intens.flatten()[chunks[j]].argmax()
            hguess=intens.flatten()[chunks[j]][hguess_I]-backg_base
            # hguess=(intens.flatten()[chunks[j]][hguess_I]-backg[0])*0.8
            dguess=dspace.flatten()[chunks[j]][hguess_I]
            dguess= (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/2
            wguess = (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/100
            #print hguess_I,hguess,dguess,bgguess,'hguess'

            #bgguess = backg
            singleBackg = [[] for i in range(len(backg))]
            for i in xrange(len(backg)):
                singleBackg[i] = [backg[i][0]]
            bgguess = singleBackg

            po,pc= ff.singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[j]],d0s=dguess,heights=hguess,widths=wguess,wavelength=wavelength,bg=bgguess)
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
            newBGall.append(list(po[3:]))
            newBGallErr.append(list(AllErr[3:]))

            newAziChunks.append(azichunks[j])


        #plot the fits.
        if 0:
            asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
            asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
            asdf2 = intens.flatten()[chunks[j]][asdf_ord]
            asdf4 = azimu.flatten()[chunks[j]][asdf_ord]
            asdf3 = ff.singleInt((asdf1,asdf4,singlelenbg,wavelength),tuple(po))
            #asdf3 = ff.singleInt((asdf1,asdf4,singlelenbg,wavelength),po[0],po[1],po[2],po[3])


            plt.plot(asdf1, asdf2,'.')
            plt.plot(asdf1, asdf3, '-')
            plt.show()



    ### Feed each d_0,h,w into fourier expansion function to get fit for fourier component
    ### parameters as output.

    # print '\n Doing d0 fourier fit...', d0_order

    dfour = ff.Fourier_fit(np.array(newAziChunks),np.array(newd0),terms=d0_order, errs=np.array(newd0Err))
    hfour = ff.Fourier_fit(np.array(newAziChunks),np.array(newHall),terms=h_order, errs=np.array(newHallErr))
    wfour = ff.Fourier_fit(np.array(newAziChunks),np.array(newWall),terms=w_order, errs=np.array(newWallErr))

    newBGall = np.array(newBGall)
    newBGallErr = np.array(newBGallErr)
    bgfour=[]
    for i in xrange(len(lenbg)):
        tempbg,tempbgc = ff.Fourier_fit(np.array(newAziChunks),np.array(newBGall[:,i]),terms=bg_order[i], errs=np.array(newBGallErr[:,i]))
        bgfour.append(tempbg)



    #plot output of fourier fits....
    if 1:
        fig = plt.figure()
        ax = fig.add_subplot(3,1,1)
        plt.subplot(411)
        plt.plot(newAziChunks,newd0, 'bo')
        plt.title('D-spacing')
        plt.plot(newAziChunks,ff.Fourier_expand(np.array(newAziChunks), dfour[0]), 'r-')
        plt.subplot(412)
        plt.plot(newAziChunks,newHall, 'bo')
        plt.title('height')
        plt.plot(newAziChunks,ff.Fourier_expand(np.array(newAziChunks), hfour[0]), 'r-')
        plt.subplot(413)
        plt.plot(newAziChunks,newWall, 'bo')
        plt.title('width')
        plt.plot(newAziChunks,ff.Fourier_expand(np.array(newAziChunks), np.array(wfour[0])), 'r-')
        plt.subplot(414)
        plt.plot(newAziChunks,newBGall, 'bo')
        plt.title('Background')
        plt.plot(newAziChunks,ff.Fourier_expand(np.array(newAziChunks), np.array(bgfour[0])), 'r-')
        plt.show()
        plt.close()



    ### Feed each d,h,w into versions of main equation to fit only or one at a time
    ### including the new paramters from the Fourier fits above as initial inputs
    ## Slice arrays for twothetarange

    print '\nRe-fitting for d, h, w separately...\n'



    ## Fit for background with existing d0, h ,w

    d0s = ff.Fourier_expand(azimu.flatten(),dfour[0])
    heights  = ff.Fourier_expand(azimu.flatten(),hfour[0])
    widths   = ff.Fourier_expand(azimu.flatten(),wfour[0])
    newbgfit = ff.bgfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),d0s,heights,widths,wavelength,bg=bgfour)
    newbg = ff.update_backgrnd(newbgfit[0],lenbg,0)


    ## Fit for d0 with existing h and w

    newdfour = ff.dfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),dfour[0],heights,widths,wavelength,bg=backg)
    print 'Old d coefficients: ', dfour[0]
    print 'New d coefficients: ', newdfour[0]

    ## Fit for h with new d and existing w

    newd0s = ff.Fourier_expand(azimu.flatten(),newdfour[0])
    # widths = ff.Fourier_expand(azimu.flatten(),wfour[0])
    #print hfour[0]
    newhfour = ff.hfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newd0s,widths,hfour[0],wavelength,bg=backg)
    print 'Old h coefficients: ', hfour[0]
    print 'New h coefficients: ', newhfour[0]

    ## Fit for w with new d and h

    # newd0s     = ff.Fourier_expand(azimu.flatten(),newdfour[0])
    newheights = ff.Fourier_expand(azimu.flatten(),newhfour[0])
    #print wfour[0]
    newwfour   = ff.wfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newd0s,newheights,wfour[0],wavelength,bg=backg)
    print 'Old w coefficients: ', wfour[0]
    print 'New w coefficients: ', newwfour[0]


    ### Refit through full equation with all data for d,h,w,bg independently

    print '\nFinal fit solving for all parms...\n'

    newwidths = ff.Fourier_expand(azimu.flatten(),newwfour[0])

    Finalparms = ff.Allfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newdfour[0],newhfour[0],newwfour[0],wavelength,bg=backg)

    starth = len(newdfour[0])
    startw = len(newdfour[0])+len(newhfour[0])
    end = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])

    fin_d = Finalparms[0][0:starth]
    fin_h = Finalparms[0][starth:startw]
    fin_w = Finalparms[0][startw:end]
    # reconstruct background coeffs list
    backg = ff.update_backgrnd(Finalparms[0],lenbg,end)

    #print len(fin_d)

    print 'Final d0 coefficients: ', Finalparms[0][0:starth]
    print 'Final h coefficients: ', Finalparms[0][starth:startw]
    print 'Final w coefficients: ', Finalparms[0][startw:end]
    print 'Final background coefficients: ', backg




    ### Plot results to check
    if 1:
        print '\nPlotting results for fit...\n'

        #tthchunk = np.where((twotheta.flatten()>=tthRange[0])&(twotheta.flatten()<=tthRange[1]))
        # print Finalparms[0]
        fullfit_intens = ff.Allchange((twotheta.flatten(),azimu.flatten(),lenbg,np.array([len(fin_d),len(fin_h),len(fin_w)]), wavelength), *Finalparms[0].squeeze())

        # Set the maximum to the n+1th value in the array
        # FIX ME: the max find is not necessarily efficient. A quicker way should be found if it exists.
        max_pos = 5
        
        ValMax1 = -np.sort(-intens.flatten())[max_pos]#np.max(intens.flatten())
        ValMax2 = -np.sort(-fullfit_intens.flatten())[max_pos]#np.max(fullfit_intens.flatten())
        ValMax3 = -np.sort(-(intens.flatten()-fullfit_intens.flatten()))[max_pos]#np.max((intens.flatten()-fullfit_intens.flatten()))
        ValMax  = np.max([ValMax1, ValMax2])
        ValMin1 = np.min(intens.flatten())
        ValMin2 = np.min(fullfit_intens.flatten())
        ValMin3 = np.min((intens.flatten()-fullfit_intens.flatten()))
        ValMin  = np.min([ValMin1, ValMin2])


        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        plt.subplot(131)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(dspace.flatten(),azimu.flatten(), s=4, c=(intens.flatten()), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
        plt.subplot(132)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(dspace.flatten(), azimu.flatten(), s=4, c=(fullfit_intens.flatten()), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin, vmax=ValMax)
        plt.colorbar()
        plt.subplot(133)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(dspace.flatten(), azimu.flatten(), s=4, c=(intens.flatten()-fullfit_intens.flatten()), edgecolors='none', cmap=plt.cm.jet, vmin=ValMin3, vmax=ValMax3)
        plt.colorbar()
        plt.show()

        plt.close()


    '''
    collated = {"d-space":    fin_d,
                "height":     fin_h,
                "width":      fin_w,
                "background": backg
    }
    '''
    
    tth_min = twotheta.min()
    tth_max = twotheta.max()
    d_min = dspace.min()
    d_max = dspace.max()
    extent = [[tth_min, tth_max],[d_min, d_max]]


    peaks = []
    peaks.append({"d-space":    fin_d,
                "height":     fin_h,
                "width":      fin_w,
                "profile":      0
    })
    collated = {"background": backg, "range": extent, "peak":peaks}


    return collated
    


