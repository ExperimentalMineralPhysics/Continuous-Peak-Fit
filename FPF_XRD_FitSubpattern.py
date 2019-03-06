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



def FitSubpattern(TwoThetaAndDspacings, azimu, intens, orders=None, PreviousParams=None):

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

    if orders:
        print orders
        total = orders['AziBins']
        backg_type = 'order'
        bg_order = orders['background']

        if 'symmetry' in orders['peak'][0]:
            symmetry = orders['peak'][0]['symmetry']
        else:
            symmetry = 1

        print bg_order
        backg=[]
        for j in xrange(len(bg_order)):
            backg.append(np.zeros(len(bg_order)*2+1))

    if PreviousParams and orders:
        #need to check sizes of both arrays here.
        #orders takes presidence over PreviousParams so we add or remove coefficients as necessary

        # print 'orders'
        # print orders
        # print 'Previous Fit'
        # print PreviousParams
        # print ff.Fourier_order(PreviousParams['peak'][0]['d-space'])
        # print ff.Fourier_order(PreviousParams['peak'][0]['height'])
        # print ff.Fourier_order(PreviousParams['peak'][0]['profile'])


        peeks = 1 #FIX ME: replace with loop over length of peaks.
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


    # if orders:
    #     d0_order, h_order, w_order, bg_order, total = orders


    ## Split data into chunks in azimuth, replace azimuths with avg.
    ## doesn't have to be 100

    if not PreviousParams:

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
        newProfileAll=[]
        newProfileAllErr=[]
        newAziChunks=[]

        for j in range(len(chunks)):


            # print '\nFitting azimuth chunk....\n'

            #only compute the fit and save it if there are 'enough' data
            if np.sum(~intens.flatten()[chunks[j]].mask) >= 20: #FIX ME: this switch is.a crude way of finding if there are less than 20 data in the bin. It should also be a variable. 

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
                hguess_I=intens.flatten()[chunks[j]].argmax()
                hguess=intens.flatten()[chunks[j]][hguess_I]-backg_base
                # hguess=(intens.flatten()[chunks[j]][hguess_I]-backg[0])*0.8
                dguess=dspace.flatten()[chunks[j]][hguess_I]
                #dguess= (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/2
                wguess = (np.max(dspace.flatten()[chunks[j]]) + np.min(dspace.flatten()[chunks[j]]))/40

                #if profile_fixed exists then set pguess to fixed value
                #then set pfixed, 1 if fixed else 0
                pguess = 0.5 #guess half if solving for profile
                pfixed = 0   #set to 0 so solvering unless profile_fixed exists.
                if 'profile_fixed' in orders['peak'][0]:
                    pguess = orders['peak'][0]['profile_fixed'][0]   #FIX ME: this should be a Fourier expansion because it is fixed....
                    pfixed = 1

                #print hguess_I,hguess,dguess,bgguess,'hguess'
                #print wguess
                #bgguess = backg
                singleBackg = [[] for i in range(len(backg))]
                for i in xrange(len(backg)):
                    singleBackg[i] = [backg[i][0]]
                bgguess = singleBackg
                #print 'bgguess ', bgguess
                po,pc= ff.singleFit(intens.flatten()[chunks[j]],twotheta.flatten()[chunks[j]],azichunks[j],dspace.flatten()[chunks[j]],d0s=dguess,heights=hguess,widths=wguess, profile=pguess,fixed=pfixed,wavelength=wavelength,bg=bgguess)
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
                if pfixed==1:
                    newBGall.append(list(po[3:]))
                    newBGallErr.append(list(AllErr[3:]))
                    newProfileAll.append(pguess)
                    newProfileAllErr.append(0)
                else:
                    newBGall.append(list(po[3:np.size(po)-1])) #FIX ME: need right length for profile.
                    newBGallErr.append(list(AllErr[3:np.size(po)-1]))
                    newProfileAll.append(po[-1])
                    newProfileAllErr.append(AllErr[-1])

                #print 'newWall ',newWall
                newAziChunks.append(azichunks[j])
                #print 'newBGall ',newBGall

            #plot the fits.
            if 0:
                asdf_ord = np.argsort(twotheta.flatten()[chunks[j]])
                asdf1 = twotheta.flatten()[chunks[j]][asdf_ord]
                asdf2 = intens.flatten()[chunks[j]][asdf_ord]
                asdf4 = azimu.flatten()[chunks[j]][asdf_ord]
                if 'profile_fixed' in orders['peak'][0]:
                    ins = (asdf1,asdf4,singlelenbg,wavelength,pguess)
                else:
                    ins = (asdf1,asdf4,singlelenbg,wavelength)
                asdf3 = ff.singleInt(ins,tuple(po))
                asdf3 = ff.singleInt((asdf1,asdf4,singlelenbg,wavelength,pguess),tuple(po))
                #asdf3 = ff.singleInt((asdf1,asdf4,singlelenbg,wavelength),po[0],po[1],po[2],po[3])


                plt.plot(asdf1, asdf2,'.')
                plt.plot(asdf1, asdf3, '-')
                plt.show()



        ### Feed each d_0,h,w into fourier expansion function to get fit for fourier component
        ### parameters as output.

        # print '\n Doing d0 fourier fit...', d0_order

        d0_order = orders['peak'][0]['d-space']
        h_order = orders['peak'][0]['height']
        w_order = orders['peak'][0]['width']
        p_order = orders['peak'][0]['profile']

        dfour = ff.Fourier_fit(np.array(newAziChunks),np.array(newd0),terms=d0_order, errs=np.array(newd0Err))
        hfour = ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newHall),terms=h_order, errs=np.array(newHallErr))
        wfour = ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newWall),terms=w_order, errs=np.array(newWallErr)*np.array(newWallErr))
        if 'profile_fixed' in orders['peak'][0]:
            pfour = orders['peak'][0]['profile_fixed']
        else:
            pfour = ff.Fourier_fit(np.array(newAziChunks)*symmetry,np.array(newProfileAll),terms=p_order, errs=np.array(newProfileAllErr)*np.array(newProfileAllErr))
        
        newBGall = np.array(newBGall)
        newBGallErr = np.array(newBGallErr)
        bgfour=[]
        for i in xrange(len(lenbg)):
            #print bg_order[i]
            #print np.array(newBGall[:,i])
            #print np.array(newBGallErr[:,i])
            tempbg,tempbgc = ff.Fourier_fit(np.array(newAziChunks),np.array(newBGall[:,i]),terms=bg_order[i], errs=np.array(newBGallErr[:,i]))
            bgfour.append(tempbg)

        #plot output of fourier fits....
        if 1:

            AziPlot = range(0,360,2)
            fig = plt.figure()
            ax = fig.add_subplot(4,1,1)
            plt.subplot(511)
            plt.plot(newAziChunks,newd0, 'bo')
            plt.title('D-spacing')
            plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot), dfour[0]), 'r-')
            plt.subplot(512)
            plt.plot(newAziChunks,newHall, 'bo')
            plt.title('height')
            plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symmetry, hfour[0]), 'r-')
            plt.subplot(513)
            plt.plot(newAziChunks,newWall, 'bo')
            plt.title('width')
            plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symmetry, np.array(wfour[0])), 'r-')
            plt.subplot(514)
            plt.plot(newAziChunks,newBGall, 'bo')
            plt.title('Background')
            plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot), np.array(bgfour[0])), 'r-')
            print ff.Fourier_expand(np.array(AziPlot), np.array(bgfour[0]))
            plt.subplot(515)
            plt.plot(newAziChunks,newProfileAll, 'bo')
            plt.title('Profile')
            plt.plot(AziPlot,ff.Fourier_expand(np.array(AziPlot)*symmetry, np.array(pfour[0])), 'r-')
            plt.show()
            plt.close()


        #put fits into data structure.
        #FIX ME: N.B. As written now this only works for a single peak.
        peaks = []
        peaks.append({"d-space":    dfour[0],
                        "height":   hfour[0],
                        "width":    wfour[0],
                        "profile":  pfour[0]})
        NewParams = []
        NewParams = {"background": bgfour[0], "peak":peaks}

        print NewParams
        ### Feed each d,h,w into versions of main equation to fit only or one at a time
        ### including the new paramters from the Fourier fits above as initial inputs
        ## Slice arrays for twothetarange

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

                newfour = ff.ParamFit(Parm, NewParams['peak'][peek][Parm], intens.flatten(),twotheta.flatten(),azimu.flatten()*sym,dspace.flatten(),d0s,heights,widths,profiles,wavelength,bg=bgexp)

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

            newfour = ff.ParamFit(Parm, NewParams[Parm], intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),d0s,heights,widths,profiles,wavelength,bg=bgexp)

            NewParams[Parm] = newfour[0]
            print 'New ', Parm, ' coefficients: ', NewParams[Parm]
            bgexp = ff.Fourier_expand(azimu.flatten(),NewParams[Parm])


            print  NewParams['peak'][peek]['profile']
            #print NewParams

            '''
            #newdfour = ff.dfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),dfour[0],heights,widths,wavelength,bg=backg)
            newdfour = ff.dfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),dfour[0],heights,widths,wavelength,bg=bgexp)
            print 'Old d coefficients: ', dfour[0]
            print 'New d coefficients: ', newdfour[0]

            ## Fit for h with new d and existing w

            newd0s = ff.Fourier_expand(azimu.flatten(),newdfour[0])
            # widths = ff.Fourier_expand(azimu.flatten(),wfour[0])
            #print hfour[0]
            #newhfour = ff.hfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newd0s,widths,hfour[0],wavelength,bg=backg)
            newhfour = ff.hfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newd0s,widths,hfour[0],wavelength,bg=bgexp)
            print 'Old h coefficients: ', hfour[0]
            print 'New h coefficients: ', newhfour[0]

            ## Fit for w with new d and h

            # newd0s     = ff.Fourier_expand(azimu.flatten(),newdfour[0])
            newheights = ff.Fourier_expand(azimu.flatten(),newhfour[0])
            #print wfour[0]
            #newwfour   = ff.wfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newd0s,newheights,wfour[0],wavelength,bg=backg)
            newwfour   = ff.wfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newd0s,newheights,wfour[0],wavelength,bg=bgexp)
            print 'Old w coefficients: ', wfour[0]
            print 'New w coefficients: ', newwfour[0]
            newwfour = wfour
            
            newbgfit = ff.bgfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),dspace.flatten(),d0s,heights,widths,wavelength,bg=bgexp)
            newbg = ff.update_backgrnd(newbgfit[0],lenbg,0)
            print 'Old bg coefficients: ', bgfour[0]
            print 'New bg coefficients: ', newbg

            backg = newbg
            '''
        # else:
        #     newdfour = dfour
        #     newhfour = hfour
        #     newwfour = wfour
        #     backg = bgfour

    else:
        NewParams = PreviousParams

        # newdfour = []
        # newdfour.append(NewParams['peak'][0]['d-space'])
        # newhfour = []
        # newhfour.append(NewParams['peak'][0]['height'])
        # newwfour = []
        # newwfour.append(NewParams['peak'][0]['width'])

        # backg = ff.update_backgrnd(NewParams['background'],lenbg,0)


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

    Finalparms = ff.Allfit(intens.flatten(),twotheta.flatten(),azimu.flatten(),newdfour[0],newhfour[0],newwfour[0],newpfour[0],wavelength,bg=backg, symm=symmetry, fix=fixed)

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
        inp = (twotheta.flatten(),azimu.flatten(),lenbg,np.array([len(fin_d),len(fin_h),len(fin_w)]), wavelength, symmetry, newpfour[0])
    else:
        inp = (twotheta.flatten(),azimu.flatten(),lenbg,np.array([len(fin_d),len(fin_h),len(fin_w),len(fin_p)]), wavelength, symmetry)
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
        ax1.set_ylim([0,360])
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
        ax2.set_ylim([0,360])
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        ax3 = plt.subplot(133)
        #plt.scatter(twotheta, azimu, s=1, c=np.log(intens-fullfit_intens), edgecolors='none', cmap=plt.cm.jet)
        plt.scatter(twotheta.flatten(), azimu.flatten(), s=dot_size, c=(intens.flatten()-fullfit_intens.flatten()), edgecolors='none', cmap=plt.cm.coolwarm, vmin=ValMin3, vmax=ValMax3)
        ax3.set_title('Residuals (data - model)')
        ax3.set_xlabel('2 Theta (deg)')
        ax3.set_ylabel('Azimuth (deg)')
        ax3.set_xlim([np.min(twotheta.flatten()),np.max(twotheta.flatten())])
        ax3.set_ylim([0,360])
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        plt.colorbar()

        plt.tight_layout()

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
    peaks.append({"d-space":     fin_d,
                "height":        fin_h,
                "width":         fin_w,
                "profile":       fin_p, 
                "profile_fixed": fin_p, #FIX ME: this is wrong here! it might not exist.
                "symmetry":      symmetry
    })
    collated = {"background": backg, "range": extent, "peak":peaks}

    return collated
    


