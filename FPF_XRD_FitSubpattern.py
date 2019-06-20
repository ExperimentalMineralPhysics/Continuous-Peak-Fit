# /usr/bin/python

# FPF_XRD_FitSubpattern
# Script fits subset of the data with peaks of pre-defined Fourier order

import numpy as np
#import numpy.ma as ma
import os
import matplotlib.pyplot as plt
np.set_printoptions(threshold='nan')
import FPF_PeakFourierFunctions as ff


# FIX ME;
# The script should not need to know what the limits are - it should only see the data that it needs to fit. 


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

    # FIX ME: need to match orders of arrays to previous numbers.
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

#        #FIX ME: assumes a fixed number of peaks.
#        if 'symmetry' in orders['peak'][0]:
#            symmetry = orders['peak'][0]['symmetry']
#        else:
#            symmetry = 1

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
            
            #FIX ME: Need to check that the number of peaks for which we have parameters is the same as the number of peaks guessed at. 
            
            tthguesses = np.array(orders['PeakPositionSelection'])
            
            #fit fourier series to two theta/d-spacing for each peak
            # FIX ME: this fitting should be done to the dspacing not the tth/energy. Probably need to import the detector functions to make this happen.
            dfour = []
            for j in range(peeks):
                
                peek = j+1
                tthguess = tthguesses[tthguesses[:,0]==peek,1:]
                
                dfour.append(ff.Fourier_fit(tthguess[:,0],det.Conversion(tthguess[:,1], conversion_factor),terms=orders['peak'][j]['d-space']))
            
                print 'tthguess:', tthguess
                print 'dfour:',dfour
        
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
        newd0      = [[] for _ in range(peeks)] 
        newd0Err   = [[] for _ in range(peeks)] 
        newHall    = [[] for _ in range(peeks)] 
        newHallErr = [[] for _ in range(peeks)] 
        newWall    = [[] for _ in range(peeks)] 
        newWallErr = [[] for _ in range(peeks)] 
        newPall    = [[] for _ in range(peeks)] 
        newPallErr = [[] for _ in range(peeks)] 
        newBGall    = [[] for _ in range(len(bg_order))] 
        newBGallErr = [[] for _ in range(len(bg_order))] 
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
                                
                peaks = []
                    
                if 'dfour' in locals():
                    #if the positions have been pre guessed extract values from dfour.
                                        
                    for k in range(peeks):
                        
                        #guess returns d-spacing then converted to two theta.                         
                        dguess = ff.Fourier_expand(np.mean(azimu.flatten()[chunks[j]]), dfour[k][0])
                        tthguess = (det.Conversion(dguess, conversion_factor, reverse=1))
                    
                        tthind = ((np.abs(twotheta.flatten()[chunks[j]] - tthguess)).argmin() )
                        #FIX ME: this takes the closest tth value to get the inensity. The mean of a number of the smallest values would be more stable.
                       
                        hguess = (intens.flatten()[chunks[j]][tthind]-backg_base)  #FIX ME: This is very crude. Need a better way to do it.                 
                        wguess = ( (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/40)  #FIX ME: This is very crude. Need a better way to do it. \
                        
                        pguess = 0.5 #guess half if solving for profile
                        pfixed = 0   #set to 0 so solvering unless profile_fixed exists.
                        if 'profile_fixed' in orders['peak'][k]:
                            pguess = orders['peak'][k]['profile_fixed']  #FIX ME: this should be a Fourier expansion because it is fixed....
                            pfixed = 1

                        peaks.append({"d-space": [dguess], 
                                      "height": [hguess], 
                                      "profile": [pguess], 
                                      "width": [wguess]
                                      })  
            
                else:
                    # guess guesses from data slice - assuming a single peak
                    hguess_I=intens.flatten()[chunks[j]].argmax()
                    hguess=intens.flatten()[chunks[j]][hguess_I]*.9-backg_base  #FIX ME: This is very crude. Need a better way to do it. 
                    
                    n = 5
                    idx = (-intens.flatten()[chunks[j]]).argsort()[:n]
                    hguess = (intens.flatten()[chunks[j]])[idx[n-1]]-5
                    
                    dguess=dspace.flatten()[chunks[j]][hguess_I]
                    
                    dguess = dspace.flatten()[chunks[j]][idx[n-1]]
                    
                    wguess = 0.02#(np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/40  #FIX ME: This is very crude. Need a better way to do it. \
                    
                    #if profile_fixed exists then set pguess to fixed value
                    #then set pfixed, 1 if fixed else 0
                    pguess = 0.5 #guess half if solving for profile
                    pfixed = 0   #set to 0 so solvering unless profile_fixed exists.
                    if 'profile_fixed' in orders['peak'][0]:
                        pguess = orders['peak'][0]['profile_fixed']   #FIX ME: this should be a Fourier expansion because it is fixed....
                        pfixed = 1
                        
                    peaks.append({"d-space": [dguess], 
                                  "height": [hguess], 
                                  "profile": [pguess], 
                                  "width": [wguess]
                                  })
                                             
                singleBackg = [[] for i in range(len(backg))]
                for i in xrange(len(backg)):
                    singleBackg[i] = [backg[i][0]]
                bgguess = singleBackg
                bgguess[0][0] = 5

                if pfixed == 1:
                    Chng = ff.ParamFitArr(peeks, len(backg), 1, p=[1,1,1,0])
                else:
                    Chng = ff.ParamFitArr(peeks, len(backg), 1)

                Guesses={  "background": bgguess, 
                    "peak": peaks }
                Guesses, po,pc = ff.FitModel(intens.flatten()[chunks[j]], twotheta.flatten()[chunks[j]], azichunks[j], Chng, Guesses, fixed=pfixed, Conv=conversion_factor)
                for i in range(peeks):
                    newd0[i].append(Guesses['peak'][i]['d-space'][0])
                    newd0Err[i].append(Guesses['peak'][i]['d-space_err'][0])
                    newHall[i].append(Guesses['peak'][i]['height'][0])
                    newHallErr[i].append(Guesses['peak'][i]['height_err'][0])
                    newWall[i].append(Guesses['peak'][i]['width'][0])
                    newWallErr[i].append(Guesses['peak'][i]['width_err'][0])
                    newPall[i].append(Guesses['peak'][i]['profile'][0])
                    if pfixed==1:
                        newPallErr[i].append(0.0001)
                    else:
                        newPallErr[i].append(Guesses['peak'][i]['profile_err'][0])
                for i in range(len(Guesses['background'])):
                    newBGall[i].append(Guesses['background'][i][0])
                    newBGallErr[i].append(Guesses['background_err'][i][0])
                   
                newAziChunks.append(azichunks[j])

                #plot the fits.
                if 0:
                    asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
                    tth_plot = twotheta.flatten()[chunks[j]][asdf_ord]
                    int_plot = intens.flatten()[chunks[j]][asdf_ord]
                    azm_plot = azimu.flatten()[chunks[j]][asdf_ord]
                    mod_plot = ff.PeaksModel(tth_plot,azm_plot,Guesses,Conv=conversion_factor)
    
                    plt.plot(tth_plot, int_plot,'.')
                    plt.plot(tth_plot, mod_plot[0], 'r.-')
                    plt.show()

               
        ### Feed each d_0,h,w into fourier expansion function to get fit for fourier component
        ### parameters as output.

        # print '\n Doing d0 fourier fit...', d0_order
        
        dfour = []
        hfour = []
        wfour = []
        pfour = []
        for j in range(peeks):
            
            #symmety 
            if 'symmetry' in orders['peak'][j].keys():
                symm = orders['peak'][j]['symmetry']
            else:
                symm = 1
            
            d0_order = orders['peak'][j]['d-space']
            h_order  = orders['peak'][j]['height']
            w_order  = orders['peak'][j]['width']
            p_order  = orders['peak'][j]['profile']
            dfour.append( ff.Fourier_fit(np.array(newAziChunks),      np.array(newd0[j]),  terms=d0_order, errs=np.array(newd0Err[j])))            
            hfour.append( ff.Fourier_fit(np.array(newAziChunks)*symm, np.array(newHall[j]),terms=h_order,  errs=np.array(newHallErr[j])))
            wfour.append( ff.Fourier_fit(np.array(newAziChunks)*symm, np.array(newWall[j]),terms=w_order,  errs=np.array(newWallErr[j])*np.array(newWallErr[j])))
            pfour.append( ff.Fourier_fit(np.array(newAziChunks)*symm, np.array(newPall[j]),terms=p_order,  errs=np.array(newPallErr[j])*np.array(newPallErr[j])))
        
        newBGall = np.array(newBGall)
        newBGallErr = np.array(newBGallErr)
        
        bgfour=[]
        
        for i in xrange(len(lenbg)):
            tempbg,tempbgc = ff.Fourier_fit(np.array(newAziChunks),np.array(newBGall[i]),terms=bg_order[i], errs=np.array(newBGallErr[i]))
            bgfour.append(tempbg.tolist())

        #plot output of fourier fits....
        if 0:

            y_lims  = np.array([np.min(newAziChunks), np.max(newAziChunks)])
            y_lims  = np.around(y_lims/180)*180
            AziPlot = range(np.int(y_lims[0]),np.int(y_lims[1]),2)
            
            
            fig = plt.figure()
            ax = fig.add_subplot(5,1,1)
            
            plt.subplot(511)
            plt.title('D-spacing')
            for j in range(peeks):
                plt.scatter(newAziChunks,newd0[j], s=20, facecolors='none', edgecolors='b') #FIX ME: should have a different colour for each peak in the subpattern
                #plt.errorbar(newAziChunks,newd0[j],newd0Err[j], linestyle="None") #FIX ME: should have a different colour for each peak in the subpattern
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot), dfour[j][0]), 'r-')
            
            plt.subplot(512)
            plt.title('height')
            for j in range(peeks):
                if 'symmetry' in orders['peak'][j].keys():
                    symm = orders['peak'][j]['symmetry']
                    print "blah"
                else:
                    symm = 1
                    print 'bummer'
                plt.scatter(newAziChunks,newHall[j], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symm, hfour[j][0]), 'r-')
            
            plt.subplot(513)
            plt.title('width')
            for j in range(peeks):
                if 'symmetry' in orders['peak'][j].keys():
                    symm = orders['peak'][j]['symmetry']
                    print "blah"
                else:
                    symm = 1
                    print 'bummer'
                plt.scatter(newAziChunks,newWall[j], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symm, np.array(wfour[j][0])), 'r-')
            
            plt.subplot(514)
            plt.title('Profile')
            for j in range(peeks):
                if 'symmetry' in orders['peak'][j].keys():
                    symm = orders['peak'][j]['symmetry']
                    print "blah"
                else:
                    symm = 1
                    print 'bummer'
                plt.scatter(newAziChunks,newPall[j], s=20, facecolors='none', edgecolors='b')
                plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symm, np.array(pfour[j][0])), 'r-')
            
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
        peaks = []
        for j in range(peeks):
            if 'symmetry' in orders['peak'][j]: #FIX ME: the symmetry part of this is a horrible way of doing it.
                peaks.append({"d-space":    dfour[j][0].tolist(),
                                "height":   hfour[j][0].tolist(),
                                "width":    wfour[j][0].tolist(),
                                "profile":  pfour[j][0].tolist(),
                                "symmetry": orders['peak'][j]['symmetry']})
            else:
                peaks.append({"d-space":    dfour[j][0].tolist(),
                                "height":   hfour[j][0].tolist(),
                                "width":    wfour[j][0].tolist(),
                                "profile":  pfour[j][0].tolist()})
        NewParams = []
        NewParams = {"background": bgfour, "peak":peaks}
        
        
        ### Interate over each fourier series in turn.
        
        refine = 1         #FIX ME: should be a variable
        iterations = 1     #FIX ME: should be a variable
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
                        Chng[0][k][l] = 1 #necessary?
                        # only 4 parameters describe each peak
                        # set paramter to change to length of array
                        # FIX ME: should this just be 1/0 as a switch or the order of the array - then we pass the order array around and match things up later?
                        if l == 0:
                            Chng[0][k][l] = orders['peak'][k]['d-space']*2 +1
                            Parm = 'd-space'
                        elif l == 1:
                            Chng[0][k][l] = orders['peak'][k]['height']*2 +1
                            Parm = 'height'
                        elif l == 2:
                            Chng[0][k][l] = orders['peak'][k]['width']*2 +1
                            Parm = 'width'
                        elif l == 3:
                            Chng[0][k][l] = orders['peak'][k]['profile']*2 +1
                            Parm = 'profile'
                            
                        print Chng                        
                        NewParams, po,pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng, NewParams, fixed=pfixed, Conv=conversion_factor)
                        #NewParams = ModParams
                        #print 'Parms out', NewParams, '\n'
                    
                #refine background
                Chng = ff.ParamFitArr(peeks, len(backg),0)
                for k in range(len(orders['background'])):
                    Chng[1][k] = orders['background'][k]*2 +1
                print Chng
                NewParams, po,pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng, NewParams, fixed=pfixed, Conv=conversion_factor)    
                    

    else:
        NewParams = PreviousParams
        pfixed = 'something'



    ### Refit through full equation with all data for d,h,w,bg independently
    print '\nFinal fit solving for all parms...\n'
    Chng = [[],[]]
    for x in range(len(orders['peak'])):
        if 'profile_fixed' in orders['peak'][x]:
            Chng[0].append([orders['peak'][x]['d-space']*2+1,orders['peak'][x]['height']*2+1,orders['peak'][x]['width']*2+1,0])
        else:
            Chng[0].append([orders['peak'][x]['d-space']*2+1,orders['peak'][x]['height']*2+1,orders['peak'][x]['width']*2+1,orders['peak'][x]['profile']*2+1])
    for x in range(len(orders['background'])):
        Chng[1].append(orders['background'][x]*2+1)
    print Chng
    NewParams, po,pc = ff.FitModel(intens.flatten(), twotheta.flatten(), azimu.flatten(), Chng, NewParams, fixed=pfixed, Conv=conversion_factor)


    print 'Final Coefficients\n', NewParams


    #Calculate SSD for fit. 
    inp,pen = ff.PeaksModel(twotheta.flatten(), azimu.flatten(), NewParams, Conv=conversion_factor)
    SSD = np.sum( (intens.flatten() - inp.flatten())**2 ) 
    print 'SSD', SSD    
    
    #We choose the degree of polynomial for which the variance as computed by
    #Sr(m)/(n-m-1)
    #is a minimum or when there is no significant decrease in its value as the degree of polynomial is increased. In the above formula,
    #Sr(m) = sum of the square of the residuals for the mth order polynomial
    #n= number of data points
    #m=order of polynomial (so m+1 is the number of constants of the model)

    # nunmber elements
    n_points = intens.flatten().size
    print 'number data', n_points

    #degrees of freedom (number of parameters solved for).
    deg_freedom = sum(ff.flatten(Chng))
    print 'degrees of freedom', deg_freedom


    SSD_var = SSD/(n_points + deg_freedom -1)
    print 'variance', SSD_var


    ### Plot results to check
    if 1:
        print '\nPlotting results for fit...\n'

        fullfit_intens = inp
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
        
        if conversion_factor['DispersionType']  == 'EnergyDispersive':
            xlabel_str = 'Energy (keV)'
        elif conversion_factor['DispersionType']  == 'AngleDispersive':
            xlabel_str = 'Two theta (deg)'
        else:
            xlabel_str = 'Dispersive units'
            
        dot_size = 16

        fig = plt.figure()
        ax = fig.add_subplot(1,3,1)
        ax1 = plt.subplot(131)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(),azimu.flatten(), s=dot_size, c=(intens.flatten()), edgecolors='none', cmap=plt.cm.magma_r, vmin=ValMin, vmax=ValMax)
        ax1.set_title('Data')
        ax1.set_xlabel(xlabel_str)
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
        ax2.set_xlabel(xlabel_str)
        ax2.set_ylabel('Azimuth (deg)')
        ax2.set_xlim([np.min(twotheta.flatten()),np.max(twotheta.flatten())])
        ax2.set_ylim(y_lims)
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        ax3 = plt.subplot(133)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(intens.flatten()-fullfit_intens.flatten()), edgecolors='none', cmap=plt.cm.coolwarm, vmin=ValMin3, vmax=ValMax3)
        ax3.set_title('Residuals (data - model)')
        ax3.set_xlabel(xlabel_str)
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
    extent = [[tth_min, tth_max],[d_min, d_max]] #FIX ME: This is taking the maximum and hte minimum of the data not the 'range' itself.

    NewParams.update({'range': extent})

    return NewParams