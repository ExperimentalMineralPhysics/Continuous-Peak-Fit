
import unittest
import os
import numpy as np
import matplotlib.pyplot as plt
import cpf
import cpf.series_functions as sf
import cpf.lmfit_model as lmm
from lmfit import Parameters
import cpf.output_formatters.convert_fit_to_crystallographic as cfc
from cpf.data_preprocess import image_adjust



#class TestCreateNewparams(unittest.TestCase):
#    def setUp(self):
#        pass#0.071*2
#
#    def test_FitFourierToDspacing(self):
"""
This is somewhere between a test case and a proove that the calcualtions are correct. 


example of data_prepare:

data_prepare = {"smooth": {"Gaussian": 2},
                "mask": Calib_mask,
                "order": {"smooth", "mask"}}
"""


#%% setup data for tests
os.chdir("../Example1-Fe/")
settings = cpf.XRD_FitPattern.initiate("BCC1_MultiPeak_input_Dioptas")

#get example image.
new_data = settings.data_class
new_data.fill_data(settings.image_list[0], settings=settings)





#%% test cosmic removal.
for x in range(2):
    if x==0:
        settings.data_prepare = {"cosmics": {"gain": 2.5, "sigclip": 3.0, "objlim": 30, "sigfrac": 0.3}}
    elif x==1:
        settings.data_prepare = {"cosmics": ""}
    
    # run removal
    cosmic_data = new_data.duplicate()
    cosmic_data = image_adjust(cosmic_data, settings.data_prepare)
    
    #plot original vs. processed
    fig, ax = plt.subplots(1,3)
    #ax.plot(azm, r, 'r')
    
    pc = new_data.plot_collected(fig_plot=fig, axis_plot=ax[0], show="intensity", limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[0], location='bottom')
    ax[0].set_title("original image")
    
    pc = cosmic_data.plot_collected(fig_plot=fig, axis_plot=ax[1], show="intensity", limits=[0, 99.9], location='bottom')
    #pc = ax[1].imshow((cosmic_data.intensity))
    #fig.colorbar(pc, shrink=0.6, ax=ax[1], location='bottom')
    ax[1].set_title("cosmic-removed image")
    
    #pc = ax[2].imshow(new_data.intensity-cosmic_data.intensity) 
    pc = cosmic_data.plot_collected(fig_plot=fig, axis_plot=ax[2], show="mask", limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[2], location='bottom')
    ax[2].set_title("differences")



#%% test image smoothing
for x in range(2):
    if x==0:
        settings.data_prepare = {"smooth": ""}
    elif x==1:
        settings.data_prepare = {"smooth": {"filter": "Gaussian", "kernel": 4}}
    
    # no need to modify the image for the test.
    
    # run removal
    smooth_data = new_data.duplicate()
    smooth_data = image_adjust(smooth_data, settings.data_prepare)
    
    #plot original vs. processed
    fig, ax = plt.subplots(1,3)
    #ax.plot(azm, r, 'r')
    pc = new_data.plot_collected(fig_plot=fig, axis_plot=ax[0], show="intensity", limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[0], location='bottom')
    ax[0].set_title("original image")
    
    pc = smooth_data.plot_collected(fig_plot=fig, axis_plot=ax[1], show="intensity", limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[1], location='bottom')
    ax[1].set_title("smoothed image")
    
    pc = smooth_data.plot_collected(fig_plot=fig, axis_plot=ax[2], show=new_data.intensity-smooth_data.intensity, limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[2], location='bottom')
    ax[2].set_title("differences")



#%% test rolling ball background. 
for x in range(2):
    if x==0:
        settings.data_prepare = {"background": ""}
    elif x==1:
        settings.data_prepare = {"background": {"kernel": 10}}

    # no need to modify the image for the test.
    
    # run removal
    rolling_bg_data = new_data.duplicate()
    rolling_bg_data = image_adjust(rolling_bg_data, settings.data_prepare)
    
    #plot original vs. processed
    fig, ax = plt.subplots(1,3)
    #ax.plot(azm, r, 'r')
    pc = new_data.plot_collected(fig_plot=fig, axis_plot=ax[0], show="intensity", limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[0], location='bottom')
    ax[0].set_title("original image")
    
    pc = rolling_bg_data.plot_collected(fig_plot=fig, axis_plot=ax[1], show="intensity", limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[1], location='bottom')
    ax[1].set_title("background removed image")
    
    pc = rolling_bg_data.plot_collected(fig_plot=fig, axis_plot=ax[2], show=new_data.intensity-rolling_bg_data.intensity, limits=[0, 99.9], location='bottom')
    #fig.colorbar(pc, shrink=0.6, ax=ax[2], location='bottom')
    ax[2].set_title("differences")


#%% test background removal. 






