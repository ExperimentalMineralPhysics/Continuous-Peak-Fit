#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import gridspec, cm, colors
from cpf.histograms import histogram2d


class _Plot_AngleDispersive:
    """
    This class cannot be imported as a stand alone class. Instead the methods 
    contained within it are created as methods in the angle dispersive diffraction
    classes (e.g. DioptasDetector and SERFlvpDetector).    
    """    
    
    def _dispersion_ticks(self, disp_ticks=None, unique=10, disp_lims = None):
        """
        Returns the labels for the dispersion axis/colour bars.
    
        :param disp_ticks: -- unused for maintained for compatibility with MED functions
        :param unique:
        :param disp_lims:
        :return new_tick_positions:
        """
                
        if disp_lims is not None:
            disp_lims = np.around(np.array(disp_lims) / self.azm_blocks) * self.azm_blocks
            disp_ticks = list(range(int(disp_lims[0]), int(disp_lims[1] + 1), int(self.azm_blocks)))
        elif self.azm_end is not None and self.azm_start is not None:
            disp_lims = np.around(np.array([self.azm_start, self.azm_end]) / self.azm_blocks) * self.azm_blocks
            disp_ticks = list(range(int(disp_lims[0]), int(disp_lims[1] + 1), int(self.azm_blocks)))
        elif len(np.unique(self.azm)) <= unique:
            # N.B. len(np.unique(self.azm)) is very slow if the data set is large.
            # as long as self.azm_start and self.azm_end are set we should never get here. 
            disp_ticks = np.unique(np.unique(self.azm))
        else: #len(np.unique(self.azm)) >= unique:
            # N.B. len(np.unique(self.azm)) is very slow if the data set is large.
            # as long as self.azm_start and self.azm_end are set we should never get here. 
            disp_lims = np.array(
                [np.min(self.azm.flatten()), np.max(self.azm.flatten())]
            )
            if disp_lims[1]-disp_lims[0] > 90:
                block = 45
            else:
                block = self.azm_blocks
            disp_lims = np.around(disp_lims / block) * block   
            disp_ticks = list(range(int(disp_lims[0]), int(disp_lims[1] + 1), int(block)))
            
        return disp_ticks
    
    
    
    def plot_masked(self, fig_plot=None):
        """
        Plot all the information needed to mask the data well.
        :param fig:
        :return:
        """
    
        x_plots = 3
        y_plots = 2
        spec = gridspec.GridSpec(
            ncols=x_plots,
            nrows=y_plots,
            width_ratios=[1, 1, 1],
            wspace=0.5,
            hspace=0.5,
            height_ratios=[2, 1],
        )
    
        ax1 = fig_plot.add_subplot(spec[0])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax1,
            show="unmasked_intensity",
            x_axis="default",
            limits=[0, 100],
        )
        ax1.set_title("All Data")
        ax2 = fig_plot.add_subplot(spec[1])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax2,
            show="mask",
            x_axis="default",
            limits=[0, 100],
        )
        ax2.set_title("Mask")
        ax3 = fig_plot.add_subplot(spec[2])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax3,
            show="intensity",
            x_axis="default",
            limits=[0, 100],
        )
        ax3.set_title("Masked Data")
    
        ax4 = fig_plot.add_subplot(spec[3])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax4,
            show="unmasked_intensity",
            x_axis="default",
            y_axis="intensity",
            limits=[0, 100],
        )
    
        ax5 = fig_plot.add_subplot(spec[4])
        # plot cdf of the intensities.
        # sort the data in ascending order
        x1 = np.sort(self.intensity.data)
        x2 = np.sort(self.intensity)
    
        # get the cdf values of y
        y1 = np.arange(np.size(x1)) / float(np.size(x1))
        y2 = np.arange(np.size(x2)) / float(ma.count(x2))
    
        # ax1 = fig_1.add_subplot(1, 1, 1)
        ax5.plot(
            x1,
            y1,
        )
        ax5.plot(
            x2,
            y2,
        )
        ax5.set_title("CDF of the intensities")
    
        ax6 = fig_plot.add_subplot(spec[5])
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax6,
            show="intensity",
            x_axis="default",
            y_axis="intensity",
            limits=[0, 100],
        )
    
    
    
    def plot_fitted(self, fig_plot=None, model=None, fit_centroid=None):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """
    
        # match max and min of colour scales
        limits = {
            "max": np.max([self.intensity.max(), model.max()]),
            "min": np.min([self.intensity.min(), model.min()]),
        }
    
        # plot data
        ax1 = fig_plot.add_subplot(1, 3, 1)
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax1,
            show="intensity",
            x_axis="default",
            limits=limits,
            colourmap="magma_r",
        )
        ax1.set_title("Data")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
        # plot model
        ax2 = fig_plot.add_subplot(1, 3, 2)
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax2,
            data=model,
            limits=limits,
            colourmap="magma_r",
        )
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax2.set_title("Model")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
    
        # plot residuals
        ax3 = fig_plot.add_subplot(1, 3, 3)
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax3,
            data=self.intensity - model,
            limits=[0, 100],
            colourmap="residuals-blanaced",
        )
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                plt.plot(fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5)
        ax3.set_title("Residuals")
        locs, labels = plt.xticks()
        plt.setp(labels, rotation=90)
    
        # tidy layout
        plt.tight_layout()
    
        
    
    def plot_collected(
        self, 
        fig_plot=None, 
        axis_plot=None, 
        show="intensity", 
        colourmap="jet",
        limits=[0, 99.9], 
        location='default',
    ):
        """
        Plots the collected data as images.
    
        Parameters
        ----------
        fig_plot : TYPE, optional
            DESCRIPTION. The default is None.
        axis_plot : TYPE, optional
            DESCRIPTION. The default is None.
        show : TYPE, optional
            DESCRIPTION. The default is "intensity".
        colourmap : TYPE, optional
            DESCRIPTION. The default is "jet".
        limits : TYPE, optional
            DESCRIPTION. The default is [0, 99.9].
        location : TYPE, optional
            DESCRIPTION. The default is 'bottom'.
    
        Returns
        -------
        TYPE
            DESCRIPTION.
    
        """
        
        if axis_plot==None and fig_plot == None:
            #make a figure
            fig_plot, axis_plot = plt.subplots()
            display=True
        
        elif axis_plot == None:
            #figure exists
            #assume will be shown outside of this function
            raise ValueError("Axes is needed to plot the data into")
            # FIX ME: the figure exists and so we whould be ableto find the axes
        
        elif fig_plot == None:
            # axis exists but figure not refrenced. 
            #assume will be shown outside of this function
            raise ValueError("A figure is needed to plot the data into")
            # FIX ME: the axes exists and so we whould be ableto find the figure
        
        if isinstance(show, str) and show == "unmasked_intensity":
            plot_i = self.intensity.data
        elif isinstance(show, str) and show == "mask":
            plot_i = np.array(ma.getmaskarray(self.intensity), dtype="uint8") + 1
            colourmap = "Greys"
            limits = [0,2]
        elif isinstance(show, np.ndarray) or ma.isMaskedArray(show):
            plot_i = show
        else:  # if show == "intensity"
            plot_i = self.intensity
    
        
        if np.ndim(plot_i) >2:
            # The parametrized images to be plotted
            def f(data, n):
                return data[n]
        else:
            # there is only 1 plane of data. 
            # (single detector rather than multiple detectors or positions.)
            def f(data, n):
                return data
    
        IMax = np.percentile(plot_i, limits[1])
        IMin = np.percentile(plot_i, limits[0])
        
        if limits[0] > 0 and limits[1] < 100:
            cb_extend = "both"
        elif limits[1] < 100:
            cb_extend = "max"
        elif limits[0] > 0:
            cb_extend = "min"
        else:
            cb_extend = "neither"
    
        im_num = 0
        the_plot = axis_plot.imshow(
            f(plot_i,im_num), vmin=IMin, vmax=IMax, cmap=colourmap
        )
        
        axis_plot.set_xlabel("x")
        axis_plot.set_ylabel("y")
        axis_plot.invert_yaxis()
    
        if location == "default":
            if plot_i.shape[-1] > plot_i.shape[-2]:
                location = "bottom"
                fraction = 0.046
            else:
                location = "right"
                fraction = 0.15
    
        cb = fig_plot.colorbar(mappable=the_plot, ax=axis_plot, extend=cb_extend, fraction=fraction, location=location)
    
        if np.ndim(plot_i) >2 and plot_i.shape[0]>1:
            print("why am I her)")
            # adjust the main plot to make room for the sliders
            fig_plot.subplots_adjust(bottom=0.25)
            
            # Make a horizontal slider to control the image being viewed.
            axim = fig_plot.add_axes([0.25, 0.1, 0.65, 0.03])
            im_slider = plt.Slider(
                ax=axim,
                label='Frame number',
                valmin=0,
                valmax=plot_i.shape[0],
                valinit=0,
            )
    
            
            # The function to be called anytime a slider's value changes
            def update(val):
                the_plot.imshow(f(plot_i, im_slider.val))
                fig_plot.canvas.draw_idle()
            # FIXME: we should also be able tp update the title and the colour scale.
            
            # register the update function with each slider
            im_slider.on_changed(update)
        
        #if display==True:
        #    the_plot.show()
    
    
    
    def plot_calibrated(
        self,
        fig_plot=None,
        axis_plot=None,
        show="default",
        x_axis="default",
        y_axis="default",
        data=None,
        limits=[1, 99.9],
        colourmap="jet",
        rastered=False,
        point_scale=3,
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """
     
        if axis_plot==None and fig_plot == None:
            #make a figure
            fig_plot, axis_plot = plt.subplots()
        
        elif axis_plot == None:
            #figure exists
            pass
        
        elif fig_plot == None:
            # axis exists but figure not refrenced. 
            pass    
        
        if self.intensity.size > 1000000: # was 50000 until fixed max/min functions 
            print(
                " Have patience. The plot(s) will appear but it can take its time to render."
            )
            rastered = True
            # print(type(axis_plot))
            # axis_plot = raster_axes.RasterAxes(axes=axis_plot)
            # print(type(axis_plot))
    
        if x_axis == "default":
            plot_x = self.tth
        else:
            plot_x = self.tth
        plot_y = self.azm
    
        if y_axis == "intensity":
            # plot y rather than azimuth on the y axis
            plot_y = self.intensity
            # organise colour scale as azimuth
            plot_i = self.azm
            label_y = "Intensity (a.u.)"
            y_ticks = False
        else:  # if y_axis is "default" or "azimuth"
            plot_y = self.azm
            plot_i = self.intensity
            label_y = "Azimuth (deg)"
    
            y_lims = [self.azm_start, self.azm_end]
            y_lims = [self.azm.min(), self.azm.max()]
            axis_plot.set_ylim(y_lims)
            # y_ticks = list(range(int(y_lims[0]),int(y_lims[1]+1),45))
    
            y_ticks = self.dispersion_ticks()
    
        # organise the data to plot
        if data is not None:
            plot_i = data
        elif show == "unmasked_intensity":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_i = plot_i.data
        elif show == "mask":
            plot_x = plot_x.data
            plot_y = plot_y.data
            plot_i = np.array(ma.getmaskarray(self.intensity), dtype="uint8") + 1
            colourmap = "Greys"
        else:  # if show == "intensity"
            plot_i = plot_i
    
        # set colour bar and colour maps.
        if colourmap == "Greys":
            IMax = 2.01
            IMin = 0
            cb_extend = "neither"
        elif isinstance(limits, dict):
            IMax = limits["max"]
            IMin = limits["min"]
            cb_extend = "neither"
        else:
            if limits[1] == 100:
                IMax = np.max(plot_i)
            else:
                IMax = np.percentile(plot_i, limits[1])
            if limits[0] == 0:
                IMin = np.min(plot_i)
            else:
                IMin = np.percentile(plot_i, limits[0])
            if limits[0] > 0 and limits[1] < 100:
                cb_extend = "both"
            elif limits[1] < 100:
                cb_extend = "max"
            elif limits[0] > 0:
                cb_extend = "min"
            else:
                cb_extend = "neither"
    
        # set colour map
        if colourmap == "residuals-blanaced":
            # colourmap = self.residuals_colour_scheme(
            colourmap = residuals_colour_scheme(
                IMax, IMin
            )
        else:
            colourmap = colourmap
    			
        # set axis limits
        x_lims = [plot_x.min(), plot_x.max()]
        print(x_lims)
        axis_plot.set_xlim(x_lims)
        if y_ticks:
            axis_plot.set_yticks(y_ticks)
    
        if rastered == False:
            the_plot = axis_plot.scatter(
                plot_x,
                plot_y,
                s=.5,
                c=plot_i,
                edgecolors="none",
                cmap=colourmap,
                vmin=IMin,
                vmax=IMax,
                rasterized=rastered,
            )
        else:
            # the_plot = self.raster_plot(
            the_plot = raster_plot(
                plot_i,
                plot_x,
                plot_y,
                fig_plot=fig_plot,
                axis_plot = axis_plot,
                vmin=IMin,
                vmax=IMax,
                colourmap=colourmap,
                pixels_per_bin = 2
            )
        axis_plot.set_xlabel(r"2$\theta$ (deg)")
        axis_plot.set_ylabel(label_y)
    
        fig_plot.colorbar(mappable=the_plot, extend=cb_extend)



def residuals_colour_scheme(maximum_value, minimum_value, **kwargs):
    # create custom colormap for residuals
    # ---------------
    # Need: a colour map that is white at 0 and the colours are equally scaled on each side. So it will match the intensity in black and white. Also one this is truncated so dont have lots of unused colour bar.
    # This can't be done using DivergingNorm(vcenter=0) or CenteredNorm(vcenter=0) so make new colourmap.
    #
    # create a colour map that truncates seismic so balanced around 0.
    # It is not perfect because the 0 point insn't necessarily perfectly white but it is close enough (I think).
    n_entries = 256
    all_colours = cm.seismic(np.arange(n_entries))

    if np.abs(maximum_value) > np.abs(minimum_value):
        n_cut = np.int_(
            (
                (2 * maximum_value - (maximum_value - np.abs(minimum_value)))
                / (2 * maximum_value)
            )
            * n_entries
        )
        keep = n_entries - n_cut
        all_colours = all_colours[keep:]
    else: #if np.abs(maximum_value) < np.abs(minimum_value):
        keep= np.int_(
            (
                (
                    2 * np.abs(minimum_value)
                    - (np.abs(minimum_value)-maximum_value)
                )
                / (2 * np.abs(minimum_value))
            )
            * n_entries
        )
        all_colours = all_colours[:keep]
    all_colours = colors.ListedColormap(
        all_colours, name="myColorMap", N=all_colours.shape[0]
    )

    return all_colours
   
    
    
def raster_plot(
    # self,
    data_plot, x_plot, y_plot,
    fig_plot=None,
    axis_plot=None,
    resample_shape = None,
    vmin=0,
    vmax=np.inf,
    colourmap="jet",
    pixels_per_bin=3,
):
    """
    Converts the diffraction data into an image and plots it. It is a more efficient replacement for 
    the default scatter plots. it should be called when the data sets are too large to plot in a 
    timely manner. 
    
    In essence it integrates the diffraction data into a plot with the pixel resolution of the screen.
    It calls the same numpy functions that pyFAI uses but without the same level of data cleaning. 

    Parameters
    ----------
    data_plot : masked array
        data to convert into the image to plot. 
    x_plot : masked array
        horizontal position of the data points (usually two theta)
    y_plot : masked array
        vertical position of the data points (usually azimuth)
    fig_plot : figure, optional
        Figure to add plot to. The default is None.
    axis_plot : axes, optional
        Axes to add plot to. The default is None.
    vmin : float, optional
        minimum of the plotted colour scale. The default is 0.
    vmax : float, optional
        minimum of the plotted colour scale. The default is np.inf, in effect the maximum value in data_plot
    colourmap : string, optional
        Colourmap for the plot. The default is "jet".
    pixels_per_bin : float, optional
        Scaler for the number of bins in the histogram  . The default is 3.
     : TYPE
        DESCRIPTION.

    Returns
    -------
    pl : axes
        filled set of axes.
    """
    # FIX ME: for massive data sets (e.g. ESFR) this is still slow. 
    # maybe it should be replaced by datashader or something else that is good for dynamicaly plotting massive data sets. 
    # an alternative is to use GIS image referencing packages. 
    import time
    start = time.time()
    
    
    if resample_shape != None:
        # use the size of the output set by resample_shape
        x_bins = resample_shape[0]
        y_bins = resample_shape[1]
    else:
        #base the sixe of the output on the number of pixels in the axes. 
        if axis_plot == None:
            if fig_plot==None:
                fig = plt.figure()
            axis_plot = fig.add_subplot(1, 1, 1)
        bbox = axis_plot.get_window_extent().transformed(fig_plot.dpi_scale_trans.inverted())
        width, height  = bbox.width, bbox.height
        width *= fig_plot.dpi
        height *= fig_plot.dpi
        x_bins = width/pixels_per_bin
        y_bins = height/pixels_per_bin
        
    #FIX ME: should I just replace this with a call to histogram2d_engine in pyFAI/engines/histogram_engine.py??
    # To have empyt pixels where there is no data use dummy = np.nan as an option. 
    # see comments in histogram2d for answer
    result, x_edges, y_edges, num_pix_per_bin  = histogram2d(data_plot, x_plot, y_plot,
                x_bins = x_bins,
                y_bins = y_bins)
    result = np.rot90(result)
    pl = axis_plot.imshow(result, extent=[x_edges[0],x_edges[-1],y_edges[0],y_edges[-1]], aspect="auto", vmin=vmin, vmax=vmax, cmap=colourmap)
    
    #axis_plot.invert_yaxis()
    
    end = time.time()
    print(end - start)
    
    return pl