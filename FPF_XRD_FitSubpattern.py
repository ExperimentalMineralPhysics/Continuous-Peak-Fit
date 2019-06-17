# /usr/bin/python

# FPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier order

# The script should not need to know what the limits are - it should only see the data that it needs to fit. 
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



def FitSubpattern(TwoThetaAndDspacings, azimu, intens, orders=None, PreviousParams=None, DetFuncs='FPF_DioptasFunctions'):
    
    #load detector functions to use where needed. Mostly for converstion of tth to dspacing.
    #This is loaded from detector functions file to avoid making any assumptions about the conversion between values (if indeed any is required). 
    det = __import__(DetFuncs)

    # sort inputs

    # two theta and dspacing are both 'x' dimensions. therefore import them as single array.
    # split them up here if needed. 
    if np.shape(TwoThetaAndDspacings)[0] == 3:
        twotheta = TwoThetaAndDspacings[0]
        dspace = TwoThetaAndDspacings[1]
        conversion_factor = TwoThetaAndDspacings[2]
    else:
        twotheta = TwoThetaAndDspacings
        #FIX ME: need to seta value for 'conversion_factor' in here and propagate its none-use through the code.

    # if no parameters
    #     get orders of arrays
    # else if no orders
    #     use parameters as given
    # else
    #     change the size of the parameter arrays to match orders

    #define the number of peaks to be fitted
    peeks = len(orders['peak'])
    
    if orders:
        #print orders
        total = orders['AziBins']
        backg_type = 'order'
        bg_order = orders['background']

        #FIX ME: assumes a fixed number of peaks.
        if 'symmetry' in orders['peak'][0]:
            symmetry = orders['peak'][0]['symmetry']
        else:
            symmetry = 1

        #print bg_order
        backg=[]
        for j in xrange(len(bg_order)):
            backg.append(np.zeros(len(bg_order)*2+1))

    if PreviousParams and orders:
        #need to check sizes of both arrays here.
        #orders takes presidence over PreviousParams so we add or remove coefficients as necessary
        
        for y in range(peeks):

            #loop over parameters
            for x in range(4):

                #only 4 parameters describe each peak
                if x == 0:
                    Parm = 'd-space'
                elif x ==1:
                    Parm = 'height'
                elif x ==2:
                    Parm = 'width'
                elif x ==3:
                    Parm = 'profile'
                elif x ==4:
                    Parm = 'profile_fixed'

                if orders['peak'][y][Parm] > ff.Fourier_order(PreviousParams['peak'][y][Parm]):
                    #print 'longer'
                    changeby = orders['peak'][y][Parm]*2 + 1 - np.size(PreviousParams['peak'][y][Parm])
                    #PreviousParams['peak'][y][Parm] = np.pad(PreviousParams['peak'][y][Parm], (0,changeby), mode='constant', constant_values=0)
                    PreviousParams['peak'][y][Parm] = (PreviousParams['background'][y] + [0]*changeby)

                elif orders['peak'][y][Parm] < ff.Fourier_order(PreviousParams['peak'][y][Parm]):
                    #print 'smaller'
                    changeby = np.size(PreviousParams['peak'][y][Parm]) - ((orders['peak'][y][Parm])*2 + 1)
                    PreviousParams['peak'][y][Parm] = PreviousParams['peak'][y][Parm][1:-(changeby-1)]

                #else:
                    #print 'same'

        # loop for background orders/size
        for y in range( np.max([len(orders['background']), len(PreviousParams['background'])]) ):
            #loop over the degree of the background

            if (len(PreviousParams['background'])-1 >= y and len(orders['background'])-1 >= y):
                #if present in both arrays make sure it is the right size.

                if orders['background'][y] > ff.Fourier_order(PreviousParams['background'][y]):
                    #print 'longer'
                    changeby = orders['background'][y]*2 + 1 - np.size(PreviousParams['background'][y])
                    #print changeby
                    #PreviousParams['background'][y] = np.pad(PreviousParams['background'][y], (0,changeby), mode='constant', constant_values=0)
                    PreviousParams['background'][y] = (PreviousParams['background'][y] + [0]*changeby)

                elif orders['background'][y] < ff.Fourier_order(PreviousParams['background'][y]):
                    #print 'smaller'
                    changeby = np.size(PreviousParams['background'][y]) - ((orders['background'][y])*2 + 1)
                    PreviousParams['background'][y] = PreviousParams['background'][y][:-(changeby)]

            elif len(PreviousParams['background'])-1 >= y and len(orders['background'])-1 < y:
                #if previous params has too many degrees remove the higher ones
                PreviousParams['background'][y] = []

            elif len(PreviousParams['background'])-1 < y and len(orders['background'])-1 >= y:
                #if previous parameters does not have enough degrees add more.
                PreviousParams['background'].append([0]*(orders['background'][y]*2+1))

        #tidy up from reductions in length of PreviousParams
        PreviousParams['background'] = [x for x in PreviousParams['background'] if x != []]
        

    if PreviousParams:
        backg = PreviousParams['background']
        backg_type = 'coeffs'


    lenbg=[]
    singlelenbg=[]
    for val in backg:
        lenbg.append(len(val))
        singlelenbg.append(1)
    lenbg=np.array(lenbg)
    singlelenbg=np.array(singlelenbg)



    ## If there is not previous fit -- Fit data in azimuthal chunks.
    if not PreviousParams:
        
        if not 'PeakPositionSelection' in orders:
            print '\nNo previous fits or selections. Assuming one peak and group fitting in azimuth...\n'
            
        else:
            print '\nUsing manual selections for initial guesses...\n'
            
            #FIX ME: Need to cheak that the number of peaks for which we have parameters is the same as the number of peaks guessed at. 
            
            tthguesses = np.array(orders['PeakPositionSelection'])
            
#            peeks = np.unique(tthguesses[:,0])
#            print peeks 
#            peeks = len(orders['peak'])
#            print peeks
#            print range(peeks)
#            stop
            
            #fit fourier series to two theta/d-spacing for each peak
            # FIX ME: this fitting should be done to the dspacing not the tth/energy. Probably need to import the detector functions to make this happen.
            dfour = []
            for j in range(peeks):
                
                peek = j+1
                tthguess = tthguesses[tthguesses[:,0]==peek,1:]
                
                dfour.append(ff.Fourier_fit(tthguess[:,0],det.Conversion(tthguess[:,1], conversion_factor),terms=orders['peak'][j]['d-space']))
            
                #print 'tthguess:', tthguess
                #print 'dfour:',dfour
        
        azmax = azimu.max()
        azmin = azimu.min()
        binsize = (azmax-azmin)/total
        chunks = []
        azichunks = []


        for i in xrange(total):
            
            end = azmin + (i+1)*binsize
            start = azmin + i*binsize
            
            tempazi = azimu.flatten()
            azichunk = np.where((tempazi>start)&(tempazi<=end))
            
            if np.unique(tempazi[azichunk]).size == 1: #FIXME: this is a quick way of distinguishing between angle and energy dispersive data without having to know anything about hte number of detetors.... should be made more robust because the azimuth is needed to convert energy to dspacing and wont work if it is wrong.
                azichunks.append(tempazi[azichunk[0][0]])
            else:
                azichunks.append(((end-start)/2)+start)
            chunks.append(azichunk)
            
        ##### Final output list of azimuths with corresponding twotheta_0,h,w

        #wavelength=parms_dict['wavelength']
        newd0      = [[] for _ in range(peeks)] #[[]]*peeks
        newd0Err   = [[] for _ in range(peeks)] #[[]]*peeks
        newHall    = [[] for _ in range(peeks)] #[[]]*peeks
        newHallErr = [[] for _ in range(peeks)] #[[]]*peeks
        newWall    = [[] for _ in range(peeks)] #[[]]*peeks
        newWallErr = [[] for _ in range(peeks)] #[[]]*peeks
        newPall    = [[] for _ in range(peeks)] #[[]]*peeks
        newPallErr = [[] for _ in range(peeks)] #[[]]*peeks
        newBGall    = [[] for _ in range(len(bg_order))] #[[]]*len(bg_order)
        newBGallErr = [[] for _ in range(len(bg_order))] #[[]]*len(bg_order)
        newAziChunks=[]
        
        for j in range(len(chunks)):

            #print '\nFitting azimuth chunk....\n'
            
            #only compute the fit and save it if there are 'enough' data
            if np.sum(~intens.flatten()[chunks[j]].mask) >= 20: #FIX ME: this switch is.a crude way of finding if there are less than 20 data in the bin. It should be a variable. 

                ## background estimates ##
                if backg_type == 'coeffs':
                    backg_base = backg[0][0]
                elif backg_type == 'order':
                    #print 'uh oh!'
                    backg_base = intens.flatten()[chunks[j]].argmin()
                    backg_base = intens[backg_base]
                    if not isinstance(backg_base, int):
                        backg_base = 0
                    backg[0][0] = backg_base

                elif backg_type == 'flat':
                    backg_base = backg[0][0]
                

                #Organise guesses to be refined.
                dguess = []
                hguess = []
                wguess = []
                pguess = []
                pfixed = 0
                if 'dfour' in locals():
                    #if the positions have been pre guessed extract values from dfour.
                                        
                    for k in range(peeks):
                        
                        #guess returns d-spacing then converted to two theta. 
                        #guess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0])
                        #tthguess = (det.Conversion(guess, wavelength, reverse=1))
                        
                        dguess.append(ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0]))
                        tthguess = (det.Conversion(dguess[k], conversion_factor, reverse=1))
                        print 'tthguess', tthguess
                    
                        tthind = ((np.abs(twotheta.flatten()[chunks[j]] - tthguess)).argmin() )
                        #FIX ME: this takes the closest tth value to get the inensity. The mean of a number of the smallest values would be more stable.
                       
                        hguess.append(intens.flatten()[chunks[j]][tthind]-backg_base)  #FIX ME: This is very crude. Need a better way to do it.                 
                        wguess.append( (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/40)  #FIX ME: This is very crude. Need a better way to do it. \
                        
                        if 'profile_fixed' in orders['peak'][k]:
                            pguess.append( orders['peak'][k]['profile_fixed'] ) #FIX ME: this should be a Fourier expansion because it is fixed....
                            pfixed = 1
                        else:
                            pguess.append(0.5)
                        
  
                else:
                    # guess guesses from data slice - assuming a single peak
                    hguess_I=intens.flatten()[chunks[j]].argmax()
                    hguess=intens.flatten()[chunks[j]][hguess_I]*.9-backg_base  #FIX ME: This is very crude. Need a better way to do it. 
                    
                    #hguess=(intens.flatten()[chunks[j]][hguess_I]-backg[0])*0.8
                    n = 5
                    idx = (-intens.flatten()[chunks[j]]).argsort()[:n]
                    #print (intens.flatten()[chunks[j]])[idx]
                    #print idx

                    hguess = (intens.flatten()[chunks[j]])[idx[n-1]]-5
                    
                    dguess=dspace.flatten()[chunks[j]][hguess_I]
                    
                    dguess = dspace.flatten()[chunks[j]][idx[n-1]]
                    
                    #dguess= (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/2
                
                    wguess = 0.02#(np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/40  #FIX ME: This is very crude. Need a better way to do it. \
                    
                    #if profile_fixed exists then set pguess to fixed value
                    #then set pfixed, 1 if fixed else 0
                    pguess = 1 #guess half if solving for profile
                    pfixed = 0   #set to 0 so solvering unless profile_fixed exists.
                    if 'profile_fixed' in orders['peak'][0]:
                        pguess = orders['peak'][0]['profile_fixed']   #FIX ME: this should be a Fourier expansion because it is fixed....
                        pfixed = 1
                        
                #print 'dguess',dguess
                #print 'hguess',hguess
                #print 'wguess',wguess
                #print 'pguess',pguess
                     
                singleBackg = [[] for i in range(len(backg))]
                for i in xrange(len(backg)):
                    singleBackg[i] = [backg[i][0]]
                bgguess = singleBackg
                bgguess[0][0] = 5

                if pfixed == 1:
                    Chng = ff.ParamFitArr(peeks, len(backg), 1, p=[1,1,1,0])
                else:
                    Chng = ff.ParamFitArr(peeks, len(backg), 1)

                #Chng[0][0][3]=0
                
                d_,h_,w_,p_,bg__,d0s_err, Hs_err, Ws_err, Ps_err, bg__err, po,pc=ff.FitModel(intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]], azichunks[j], Chng, [[dguess]], [[hguess]], [[wguess]], [[pguess]], bgguess, fixed=pfixed, Conv=conversion_factor)
                #AllErr = np.sqrt(np.abs(np.diag(pc)))
                
                for i in range(peeks):
                    newd0[i].append(d_[i][0])
                    newd0Err[i].append(d0s_err[i][0])
                    newHall[i].append(h_[i][0])
                    newHallErr[i].append(Hs_err[i][0])
                    newWall[i].append(w_[i][0])
                    newWallErr[i].append(Ws_err[i][0])
                    if pfixed==1:
                        newPall[i].append(pguess)
                        newPallErr[i].append(0.0001)
                    else:
                        newPall[i].append(p_[i][0])
                        newPallErr[i].append(Ps_err[i][0])
                #print bg__
                #print newBGall
                #print len(lenbg)
                for i in range(len(lenbg)):
                    #print bg__[i][0]
                    #print newBGall[i]
                    app = bg__[i][0]
                    #print 'app', app
                    newBGall[i].append(app)
                    #print 'newBGall',newBGall
                    newBGallErr[i].append(bg__err[i][0])
                   
                newAziChunks.append(azichunks[j])
                
#                po,pc= ff.singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[j]], d0s=dguess,heights=hguess,widths=wguess, profile=pguess,fixed=pfixed, wavelength=wavelength, bg=bgguess)
#                AllErr = np.sqrt(np.abs(np.diag(pc)))
#
#                #print po
#            
#
#
#                newd0.append(po[0:peeks])
#                newd0Err.append(AllErr[0:peeks])
#                newHall.append(po[peeks:2*peeks])
#                newHallErr.append(AllErr[peeks:2*peeks])
#                newWall.append(po[2*peeks:3*peeks])
#                newWallErr.append(AllErr[2*peeks:3*peeks])
#                if pfixed==1:
#                    newBGall.append(list(po[3*peeks:]))
#                    newBGallErr.append(list(AllErr[3*peeks:]))
#                    newProfileAll.append(pguess)
#                    newProfileAllErr.append([0] * peeks)
#                else:
#                    newBGall.append(list(po[3*peeks:np.size(po)-peeks])) #FIX ME: need right length for profile.
#                    newBGallErr.append(list(AllErr[3*peeks:np.size(po)-peeks]))
#                    newProfileAll.append(po[-peeks])
#                    newProfileAllErr.append(AllErr[-peeks])
#                    
#                newAziChunks.append(azichunks[j])

                #plot the fits.
                if 1:
                    asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
                    asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
                    asdf2 = intens.flatten()[chunks[j]][asdf_ord]
                    asdf4 = azimu.flatten()[chunks[j]][asdf_ord]
                    if 'profile_fixed' in orders['peak'][0]:
                        print 'pguess'
                        ins = (asdf1,asdf4,singlelenbg,conversion_factor,pguess)
                    else:
                        print 'no pguess'
                        #profile in the fitted parameter list
                        ins = (asdf1,asdf4,singlelenbg,conversion_factor)
                    #print 'ins and po', ins, po
                    #asdf3 = ff.singleInt(ins,tuple(po))
                    #asdf3 = ff.singleInt((asdf1,asdf4,singlelenbg,conversion_factor,pguess),tuple(po))
                    asdf5 = ff.PeaksModel(asdf1,asdf4,d_,h_,w_,p_,bg__,Conv=conversion_factor)
                    #asdf3 = ff.singleInt((asdf1,asdf4,singlelenbg,conversion_factor),po[0],po[1],po[2],po[3])
    
                    plt.plot(asdf1, asdf2,'.')
                    #plt.plot(asdf1, asdf3, 'g-')
                    plt.plot(asdf1, asdf5[0], 'r.-')
                    plt.show()

               
        ### Feed each d_0,h,w into fourier expansion function to get fit for fourier component
        ### parameters as output.

        # print '\n Doing d0 fourier fit...', d0_order
        
        dfour = []
        hfour = []
        wfour = []
        pfour = []
        for j in range(peeks):
            
            d0_order = orders['peak'][j]['d-space']
            h_order  = orders['peak'][j]['height']
            w_order  = orders['peak'][j]['width']
            p_order  = orders['peak'][j]['profile']
            dfour.append( ff.Fourier_fit(np.array(newAziChunks),         np.array(newd0[j]),  terms=d0_order, errs=np.array(newd0Err[j])))            
            hfour.append( ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newHall[j]),terms=h_order,  errs=np.array(newHallErr[j])))
            wfour.append( ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newWall[j]),terms=w_order,  errs=np.array(newWallErr[j])*np.array(newWallErr[j])))
            pfour.append( ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newPall[j]),terms=p_order,  errs=np.array(newPallErr[j])*np.array(newPallErr[j])))
        
#        
#        dfour = ff.Fourier_fit(np.array(newAziChunks),np.array(newd0),terms=d0_order, errs=np.array(newd0Err))
#        hfour = ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newHall),terms=h_order, errs=np.array(newHallErr))
#        wfour = ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newWall),terms=w_order, errs=np.array(newWallErr)*np.array(newWallErr))
#        if 'profile_fixed' in orders['peak'][0]:
#            pfour = orders['peak'][0]['profile_fixed']
#        else:
#            pfour = ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newPall),terms=p_order, errs=np.array(newProfileAllErr)*np.array(newProfileAllErr))
#        
        newBGall = np.array(newBGall)
        newBGallErr = np.array(newBGallErr)
        
        bgfour=[]
        
        for i in xrange(len(lenbg)):
            tempbg,tempbgc = ff.Fourier_fit(np.array(newAziChunks),np.array(newBGall[i]),terms=bg_order[i], errs=np.array(newBGallErr[i]))
            bgfour.append(tempbg)

        #plot output of fourier fits....
        if 1:

            y_lims  = np.array([np.min(newAziChunks), np.max(newAziChunks)])
            y_lims  = np.around(y_lims/180)*180
            AziPlot = range(np.int(y_lims[0]),np.int(y_lims[1]),2)
            
            
            fig = plt.figure()
            ax = fig.add_subplot(5,1,1)
            
            plt.subplot(511)
            plt.title('D-spacing')
            for j in range(peeks):
                #plt.scatter(newAziChunks,newd0[j], s=20, facecolors='none', edgecolors='b') #FIX ME: should have a different colour for each peak in the subpattern
                plt.errorbar(newAziChunks,newd0[j],newd0Err[j], linestyle="None") #FIX ME: should have a different colour for each peak in the subpattern
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot), dfour[j][0]), 'r-')
            
            plt.subplot(512)
            plt.title('height')
            for j in range(peeks):
                plt.scatter(newAziChunks,newHall[j], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symmetry, hfour[j][0]), 'r-')
            
            plt.subplot(513)
            plt.title('width')
            for j in range(peeks):
                plt.scatter(newAziChunks,newWall[j], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symmetry, np.array(wfour[j][0])), 'r-')
            
            plt.subplot(514)
            plt.title('Profile')
            for j in range(peeks):
                plt.scatter(newAziChunks,newPall[j], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symmetry, np.array(pfour[j][0])), 'r-')
            
            for k in range(len(lenbg)):
                x_plt = len(lenbg)
                y_plt = len(lenbg)*4 + k +1
                plt.subplot(5,x_plt,y_plt)
                plt.title('Background')
                plt.scatter(newAziChunks,newBGall[k], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot), np.array(bgfour[k][0])), 'r-')
                
            plt.show()
            plt.close()


        #put fits into data structure.
        #FIX ME: N.B. As written now this only works for a single peak.
        #FIX ME: this has been changed but not verified.
        for j in range(peeks):
            peaks = []
            peaks.append({"d-space":    dfour[j][0],
                            "height":   hfour[j][0],
                            "width":    wfour[j][0],
                            "profile":  pfour[j][0]})
        NewParams = []
        NewParams = {"background": bgfour, "peak":peaks}

        #print NewParams
        ### Feed each d,h,w into versions of main equation to fit only or one at a time
        ### including the new paramters from the Fourier fits above as initial inputs
        ## Slice arrays for twothetarange

        refine = 1
        iterations = 1
        if refine:
            print '\nRe-fitting for d, h, w, bg separately...\n'

            #iterate through the variables.
            for j in range(iterations):
                
                #refine parameters for each peak.
                for k in range(peeks):
                
                    if 'profile_fixed' in orders['peak'][0]:
                        refine_parms = 3
                    else:
                        refine_parms = 4
        
                    #loop over parameters
                    for l in range(refine_parms):
                    
                        
                        Chng = ff.ParamFitArr(peeks, len(backg),0)
                        Chng[0][k][l] = 1
                            
                        print Chng
                        
                        
                        #only 4 parameters describe each peak
#                        # set paramter to change to 1
#                        if l == 0:
#                            Parm = 'd-space'
#                        elif l == 1:
#                            Parm = 'height'
#                        elif l == 2:
#                            Parm = 'width'
#                        elif l == 3:
#                            Parm = 'profile'
#                                
#                    
#                       
#                   
                        print peaks
                        print peaks[k]['d-space']
                        print dfour
                        stop
                        dfour, hfour, wfour, pfour, bgfour, d0s_err, Hs_err, Ws_err, Ps_err, bg__err, po,pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng, dfour, hfour, wfour, pfour, bgfour, fixed=pfixed, Conv=conversion_factor, symm=symmetry )
                        
                
            
            stop


        refine = 1
        if refine:
            print '\nRe-fitting for d, h, w, bg separately...\n'
            
            #set peek for a loop over peaks that doesn't yet exist. 
            peek = 0

            d0s        = ff.Fourier_expand(azimu.flatten(),NewParams['peak'][peek]['d-space'])
            heights    = ff.Fourier_expand(azimu.flatten()*symmetry,NewParams['peak'][peek]['height'])
            widths     = ff.Fourier_expand(azimu.flatten()*symmetry,NewParams['peak'][peek]['width'])
            profiles   = ff.Fourier_expand(azimu.flatten()*symmetry,NewParams['peak'][peek]['profile'])
            background = ff.Fourier_expand(azimu.flatten(),NewParams['background'])

            print  NewParams['peak'][peek]['profile']
            ## Fit for background with existing d0, h ,w

            #newbgfit = ff.bgfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),d0s,heights,widths,wavelength,bg=bgfour)
            bgexp = ff.update_backgrnd(bgfour[0],lenbg,0)



            if 'profile_fixed' in orders['peak'][0]:
                refine_parms = 3
            else:
                refine_parms = 4

            #loop over parameters
            for x in range(refine_parms):

                #only 4 parameters describe each peak
                if x == 0:
                    Parm = 'd-space'
                    sym  = 1
                elif x ==1:
                    Parm = 'height'
                    sym  = symmetry
                elif x ==2:
                    Parm = 'width'
                    sym  = symmetry
                elif x ==3:
                    Parm = 'profile'
                    sym  = symmetry

                ## Fit for each parameter keeping others constant
                print 'Old ', Parm, ' coefficients: ', NewParams['peak'][peek][Parm]

                print conversion_factor

                newfour = ff.ParamFit(Parm, NewParams['peak'][peek][Parm], intens.flatten(),twotheta.flatten(),azimu.flatten()*sym,dspace.flatten(),d0s,heights,widths,profiles,conversion_factor['conversion_constant'],bg=bgexp)

                NewParams['peak'][peek][Parm] = newfour[0]
                print 'New ', Parm, ' coefficients: ', NewParams['peak'][peek][Parm]

                #copy new fit back into correct array
                NewParmsExp = ff.Fourier_expand(azimu.flatten()*sym,NewParams['peak'][peek][Parm])
                if x == 0:
                    d0s      = NewParmsExp
                elif x ==1:
                    heights  = NewParmsExp
                elif x ==2:
                    widths   = NewParmsExp
                elif x ==3:
                    profiles = NewParmsExp

                if 0:
                    AziPlot = range(0,360,2)
                    fig = plt.figure()
                    ax = fig.add_subplot(1,1,1)
                    plt.subplot(111)
                    plt.plot(azimu.flatten(),NewParmsExp, 'r.')
                    plt.plot(np.array(range(360)),ff.Fourier_expand(np.array(range(360))*sym,NewParams['peak'][peek][Parm]),'r-')
                    plt.title(Parm)
                    plt.show()
                    plt.close()

            ## Fit for background keeping rest constant
            Parm = 'background'
            print 'Old ', Parm, ' coefficients: ', NewParams[Parm]

            newfour = ff.ParamFit(Parm, NewParams[Parm], intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),d0s,heights,widths,profiles,conversion_factor['conversion_constant'],bg=bgexp)

            NewParams[Parm] = newfour[0]
            print 'New ', Parm, ' coefficients: ', NewParams[Parm]
            bgexp = ff.Fourier_expand(azimu.flatten(),NewParams[Parm])


            print  NewParams['peak'][peek]['profile']
 

    elif orders['PeakPositionSelection']:

        print '\nUse position selections for initial guess...\n'
            
        stop
        

    else:
        NewParams = PreviousParams
        

    ### Refit through full equation with all data for d,h,w,bg independently

    print '\nFinal fit solving for all parms...\n'


    #setup the parameters
    newdfour = []
    newdfour.append(NewParams['peak'][0]['d-space'])
    newhfour = []
    newhfour.append(NewParams['peak'][0]['height'])
    newwfour = []
    newwfour.append(NewParams['peak'][0]['width'])
    newpfour = []
    newpfour.append(NewParams['peak'][0]['profile'])

    print newpfour[0]
    print newwfour[0]


    backg = ff.update_backgrnd(NewParams['background'],lenbg,0)

    #Finalparms = ff.AllfitNew(intens.flatten(),twotheta.flatten(),azimu.flatten(),NewParams,wavelength)

    if 'profile_fixed' in orders['peak'][0]:
        fixed = orders['peak'][0]['profile_fixed']
        starth = len(newdfour[0])
        startw = len(newdfour[0])+len(newhfour[0])
        end    = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])

    else:
        fixed = None
        starth = len(newdfour[0])
        startw = len(newdfour[0])+len(newhfour[0])
        startp = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])
        end    = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])+len(newpfour[0])

    Finalparms = ff.Allfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newdfour[0],newhfour[0],newwfour[0],newpfour[0],conversion_factor['conversion_constant'],bg=backg, symm=symmetry, fix=fixed)

    #starth = len(newdfour[0])
    #startw = len(newdfour[0])+len(newhfour[0])
    #startp = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])
    #end    = len(newdfour[0])+len(newhfour[0])+len(newwfour[0])+len(newpfour[0])

    fin_d = Finalparms[0][0:starth]
    fin_h = Finalparms[0][starth:startw]
    if 'profile_fixed' in orders['peak'][0]:
        fin_w = Finalparms[0][startw:end]
        fin_p = orders['peak'][0]['profile_fixed']
    else:
        fin_w = Finalparms[0][startw:startp]
        fin_p = Finalparms[0][startp:end]
    # reconstruct background coeffs list
    backg = ff.update_backgrnd(Finalparms[0],lenbg,end)

    #print len(fin_d)

    print 'Final d0 coefficients: ', fin_d
    print 'Final h coefficients: ',  fin_h
    print 'Final w coefficients: ',  fin_w
    print 'Final p coefficients: ',  fin_p
    print 'Final background coefficients: ', backg


    #Calculate SSD for fit. 

    if 'profile_fixed' in orders['peak'][0]:
        inp = (twotheta.flatten(),azimu.flatten(),lenbg,np.array([len(fin_d),len(fin_h),len(fin_w)]), conversion_factor['conversion_constant'], symmetry, newpfour[0])
    else:
        inp = (twotheta.flatten(),azimu.flatten(),lenbg,np.array([len(fin_d),len(fin_h),len(fin_w),len(fin_p)]), conversion_factor['conversion_constant'], symmetry)
    fullfit_intens = ff.Allchange(inp, *Finalparms[0].squeeze())
    SSD = np.sum( (intens.flatten() - fullfit_intens.flatten())**2 ) 
    print SSD

    #We choose the degree of polynomial for which the variance as computed by
    #Sr(m)/(n-m-1)
    #is a minimum or when there is no significant decrease in its value as the degree of polynomial is increased. In the above formula,
    #Sr(m) = sum of the square of the residuals for the mth order polynomial
    #n= number of data points
    #m=order of polynomial (so m+1 is the number of constants of the model)

    # nunmber elements
    n_points = intens.flatten().size
    print 'number data', n_points

    #degrees of freedom
    # number of parameters solved for.
    deg_freedom = len(Finalparms[0])
    print 'degrees of freedom', deg_freedom


    SSD_var = SSD/(n_points + deg_freedom -1)
    print 'variance', SSD_var


    ### Plot results to check
    if 1:
        print '\nPlotting results for fit...\n'

        #tthchunk = np.where((twotheta.flatten()>=tthRange[0])&(twotheta.flatten()<=tthRange[1]))
        # print Finalparms[0]
        #fullfit_intens = ff.Allchange((twotheta.flatten(),azimu.flatten(),lenbg,np.array([len(fin_d),len(fin_h),len(fin_w)]), wavelength), *Finalparms[0].squeeze())

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


        y_lims = np.array([np.min(azimu.flatten()), np.max(azimu.flatten())])
        y_lims = np.around(y_lims/180)*180
        
        dot_size = 16

        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        ax1 = plt.subplot(131)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(),azimu.flatten(), s=dot_size, c=(intens.flatten()), edgecolors='none', cmap=plt.cm.magma_r, vmin=ValMin, vmax=ValMax)
        ax1.set_title('Data')
        ax1.set_xlabel('2 Theta (deg)')
        ax1.set_ylabel('Azimuth (deg)')
        ax1.set_xlim([np.min(twotheta.flatten()),np.max(twotheta.flatten())])
        ax1.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        plt.colorbar()
        ax2 = plt.subplot(132)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(fullfit_intens.flatten()), edgecolors='none', cmap=plt.cm.magma_r, vmin=ValMin, vmax=ValMax)
        plt.colorbar()
        ax2.set_title('Model')
        ax2.set_xlabel('2 Theta (deg)')
        ax2.set_ylabel('Azimuth (deg)')
        ax2.set_xlim([np.min(twotheta.flatten()),np.max(twotheta.flatten())])
        ax2.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        ax3 = plt.subplot(133)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(intens.flatten()-fullfit_intens.flatten()), edgecolors='none', cmap=plt.cm.coolwarm, vmin=ValMin3, vmax=ValMax3)
        ax3.set_title('Residuals (data - model)')
        ax3.set_xlabel('2 Theta (deg)')
        ax3.set_ylabel('Azimuth (deg)')
        ax3.set_xlim([np.min(twotheta.flatten()),np.max(twotheta.flatten())])
        ax3.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        plt.colorbar()

        plt.tight_layout()

        # save figures without overwriting old names
        filename = 'Fit2Peak'
        i = 1
        while os.path.exists('{}{:d}.png'.format(filename, i)):
            i += 1
        plt.savefig('{}{:d}.png'.format(filename, i))

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
    extent = [[tth_min, tth_max],[d_min, d_max]] #FIX ME: Tis is taking the maximum and hte minimum of the data not the 'range' itself.


    peaks = []
    peaks.append({"d-space":     fin_d,
                "height":        fin_h,
                "width":         fin_w,
                "profile":       fin_p, 
                "profile_fixed": fin_p, #FIX ME: this is wrong here! it might not exist.
                "symmetry":      symmetry
    })
    collated = {"background": backg, "range": extent, "peak":peaks}

    return collated
    


