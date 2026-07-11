#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib import cm, colors, gridspec, tri, colormaps

from cpf.histograms import histogram1d, histogram2d
from cpf.util.logging import get_logger

logger = get_logger("cpf.input_types._Plot_AngleDispersive")


class _Plot_AngleDispersive:
    """
    This class cannot be imported as a stand alone class. Instead the methods
    contained within it are created as methods in the angle dispersive diffraction
    classes (e.g. DioptasDetector and SERFlvpDetector).
    """
    
    default_colourmap = "magma_r"

    def _dispersion_ticks(self, disp_ticks=None, unique=10, disp_lims=None):
        """
        Returns the labels for the dispersion axis/colour bars.

        :param disp_ticks: -- unused for maintained for compatibility with MED functions
        :param unique:
        :param disp_lims:
        :return new_tick_positions:
        """

        if disp_lims is not None:
            disp_lims = (
                np.around(np.array(disp_lims) / self.azm_blocks) * self.azm_blocks
            )
            disp_ticks = list(
                range(int(disp_lims[0]), int(disp_lims[1] + 1), int(self.azm_blocks))
            )
        elif self.azm_end is not None and self.azm_start is not None:
            num_blocks = (self.azm_end - self.azm_start) / self.azm_blocks
            if num_blocks < 3:
                self.azm_blocks = np.around((self.azm_end - self.azm_start) / 3)
            disp_lims = (
                np.around(np.array([self.azm_start, self.azm_end]) / self.azm_blocks)
                * self.azm_blocks
            )
            disp_ticks = list(
                range(int(disp_lims[0]), int(disp_lims[1] + 1), int(self.azm_blocks))
            )
        elif len(np.unique(self.azm)) <= unique:
            # N.B. len(np.unique(self.azm)) is very slow if the data set is large.
            # as long as self.azm_start and self.azm_end are set we should never get here.
            disp_ticks = np.unique(np.unique(self.azm))
        else:  # len(np.unique(self.azm)) >= unique:
            # N.B. len(np.unique(self.azm)) is very slow if the data set is large.
            # as long as self.azm_start and self.azm_end are set we should never get here.
            disp_lims = np.array(
                [np.min(self.azm.flatten()), np.max(self.azm.flatten())]
            )
            if disp_lims[1] - disp_lims[0] > 90:
                block = 45
            else:
                block = self.azm_blocks
            disp_lims = np.around(disp_lims / block) * block
            disp_ticks = list(
                range(int(disp_lims[0]), int(disp_lims[1] + 1), int(block))
            )

        return disp_ticks


    def what_plot_type(self, plot_type=None):
        """
        Determines from the data which is the best sort of plot for the data. 
        
        The types and determining factors are:
        
        - surface : when < 2e4 observations
            Trianglular filling of the plot area. VERY slow for large data sets. 
            Only called if there are less than 10000 data points.
            
        - scatter
            sctter plot of the data
        
        - rastered : when > 1e6 observations. 
            Makes rastered image of the data points. Plots points within picels rather than 
            each point as a scatter. 
            
            
        If plot_type is set then it is sanity checked and a warning raised if 
        it is not sensible, but it is allowed to proceed. 

        Parameters
        ----------
        plot_type : str, bool, optional
            type of plot expected to be used. The default is None.

        Returns
        -------
        plot_type : str
            Best type of plot for the current data set.
        """
        
        recognised_plots = ["surface",
             "surf",
             "scatter",
             "raster",
             "rast",
             "im",
             "image"]
        
        surf_threshold = 5e3
        raster_threshold = 5e5
        
        """ if the azm values are nor rounded then numberical precision can make the number of unique values 
        greater than the actual number. Round to remove imprecision. Also correctly identifies the right number of 
        values then. see _plot_angleDispersive.plot_image()
        FIXME: rounding_precision is a fixed value. should be set from the precision of the data.
        """
        rounding_precision = 8
        
        """ if the underlying data is on an orthogonal grid, that is aligned with the tth and azm axes, then the data 
        is best plotted as images (i.e. it is fastest). Howver, cannot just test if:
           unique(azm).size * unique(tth).size == intensity.size
        because this is not true if the data has been compressed and masked values discarded.
        Therefore also test if the number of unique azm and tth is significnatly less than the intensity.size.
        if the data is on square grid then this is true. otherwise not. 
        
        image_fraction_unique_threshold is the arbitrary cutoff value to test for fraction of unique values. 
        In testing values are either way above or way below this threshold
        """
        image_fraction_unique_threshold = 1/4
        
        
        if plot_type not in recognised_plots:
            # then set it
                        
            # get unique azm and tth values. 
            unique_azm = np.unique(ma.MaskedArray(np.round(self.azm, rounding_precision)).data)
            unique_tth = np.unique(ma.MaskedArray(np.round(self.tth, rounding_precision)).data)
            
            if (unique_azm.size*unique_tth.size == ma.MaskedArray(self.intensity).data.size
                or unique_azm.size / ma.MaskedArray(self.azm).data.size <= image_fraction_unique_threshold
                or unique_tth.size / ma.MaskedArray(self.tth).data.size <= image_fraction_unique_threshold
            ):
                # the is orthonormal data
                plot_type = "image"
            elif ma.MaskedArray(self.intensity).compressed().size < surf_threshold:
                plot_type = "surface"
            elif ma.MaskedArray(self.intensity).compressed().size > raster_threshold:
                plot_type = "rastered"
            else:
                plot_type = "scatter"

        logger.effusive(f"plot_type is {plot_type}")
        if plot_type not in ["rast", "rastered"] and any(ma.MaskedArray(self.intensity).compressed() > raster_threshold):
            logger.moreinfo(" Have patience. The plot(s) will appear but it can take its time to render.")
        if plot_type not in ["surf", "surface"] and any(ma.MaskedArray(self.intensity).compressed() > surf_threshold):
            logger.moreinfo(" Have patience. The plot(s) will appear but it can take its time to render.")

        return plot_type        
        

    def plot_integrated(self, fig_plot=None, axis_plot=None, show=None):
        """
        Makes a plot of the integrated data.

        :param fig:
        :return:
        """

        if fig_plot == None:
            # make a figure
            fig_plot, axis_plot = plt.subplots()

        p, i, a = histogram1d(self.tth, self.intensity, self.azm)

        axis_plot.plot(p, i)

        axis_plot.set_xlabel(f"{self.Dispersionlabel} ({self.DispersionUnits})")
        axis_plot.set_ylabel(f"{self.Observationslabel} ({self.ObservationsUnits})")
        axis_plot.set_title("Integrated Data")


    def plot_masked(self, fig_plot=None, **kwargs):
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

    def plot_range(
        self,
        fig_plot=None,
        range_bounds=[-np.inf, np.inf],
        azm_bounds=[-np.inf, np.inf],
    ):
        """
        Plot data within given range.
        add data to axes.
        :param ax:
        :param show:
        :return:
        """

        self.set_limits(range_bounds=range_bounds, azm_bounds=range_bounds)

        # match max and min of colour scales
        limits = {
            "max": self.intensity.max(),
            "min": self.intensity.min(),
        }

        # plot data
        ax1 = fig_plot.add_subplot(1, 1, 1)
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

        # tidy layout
        plt.tight_layout()

    def plot_fitted(
        self,
        fig_plot=None,
        model=None,
        fit_centroid=None,
        plot_type="default",
        orientation="horizontal",
        plot_ColourRange=None,
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """

        #check plot type
        plot_type = self.what_plot_type(plot_type=plot_type)
        # match max and min of colour scales
        if plot_ColourRange:
            if not isinstance(plot_ColourRange, dict):
                limits = {
                    "max": plot_ColourRange[0],
                    "min": plot_ColourRange[1],
                }
            else:
                limits = plot_ColourRange
        else:
            limits = {
                "max": np.nanmax([self.intensity.max(), model.max()]),
                "min": np.nanmin([self.intensity.min(), model.min()]),
            }

        tight = False

        if orientation == "horizontal":
            up = 3
            across = 1
            # location = "right"
            tight = True
            loc = "left"
        else:
            up = 1
            across = 3
            # location = "right"
            loc = "center"
        # make axes
        # try:
        #     fig_plot.set_layout_engine(layout="constrained")
        # except:
        #     pass
        ax = []
        if tight == True:
            gs = gridspec.GridSpec(up, across, wspace=0.0, hspace=0.0)
            for i in range(3):
                ax.append(fig_plot.add_subplot(gs[i]))
        else:
            # =========================
            # FIXME: legends on Data/Model/Residual figures. 
            # It would be great if we could centre the colours bars under the axes 
            # in the 3 part figure. 
            # This is trivial using plt.subplots(... layout='constrained') (see example 1 below ) but 
            # cannot use used with fig.subplots() as easily. 
            # It is also possible with gridspec and a similar over lapping of the 
            # legends and figures. But the gridspec and mixed axes are imcompatible 
            # with 'tight_layout'. 
            # 
            # see: https://matplotlib.org/stable/users/explain/axes/colorbar_placement.html
            #      https://matplotlib.org/stable/gallery/subplots_axes_and_figures/gridspec_and_subplots.html
            #      https://matplotlib.org/stable/users/explain/axes/mosaic.html#mosaic
            #
            # This has to be soluble but is not trivial given that I am passing 
            # the figure round and reusing it. It needs some thought.
            # --------------------------
            # # # Example 1. 
            # fig, axs = plt.subplots(1, 3, layout='constrained')
            # for ax in axs.flat:
            #     pcm = ax.pcolormesh(np.random.random((20, 20)))
            # fig.colorbar(pcm, ax=axs[ :2], shrink=0.5, location='bottom')
            # fig.colorbar(pcm, ax=axs[ 1:], shrink=0.5, location='bottom')
            # =========================
            # here after are options that weere possible but not complete
            # the axes that the key is for a fed to self.plot_calibrated using 
            # cbar_axes... cbar_axes can be a list of axes which will then share a 
            # colour bar. 
            if 1:
                axs = fig_plot.subplots(nrows=up, ncols=across, sharey=True, sharex=True)
                ax.append(axs[0])
                ax.append(axs[1])
                ax.append(axs[2])
                # axs = fig_plot.subplots(up, across, sharey=True, sharex=True, layout='constrained')
            else:
                gs = gridspec.GridSpec(up, across, wspace=.08, hspace=0)
                for i in range(3):
                    ax.append(fig_plot.add_subplot(gs[i]))
                    if i > 0:
                        # ax[i].sharey(ax[0])
                        ax[i].set_yticklabels([])
                        # ax[i].yaxis.set_ticks_position('both')

        # plot data
        fig_plot = self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[0],
            show="intensity",
            limits=limits,
            colourmap="magma_r",
            # location=location,
            plot_type=plot_type,
            orientation=orientation,
            # cbar_axes=False
        )
        if orientation == "horizontal":
            ax[0].set_title("Data", loc=loc, y=0.75)
        else:
            ax[0].set_title("Data")
        locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)

        # plot model
        fig_plot = self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[1],
            data=model,
            y_label=None,
            limits=limits,
            colourmap="magma_r",
            # location=location,
            plot_type=plot_type,
            orientation=orientation,
            # cbar_axes=[ax[0],ax[1]]
        )
        
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                if orientation == "horizontal":
                    ax[1].plot(
                        fit_centroid[0], fit_centroid[1][i], "k--", linewidth=0.5
                    )
                else:
                    ax[1].plot(
                        fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5
                    )
        if orientation == "horizontal":
            ax[1].set_title("Model", loc=loc, y=0.75)
        else:
            ax[1].set_title("Model")
        locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)

        # plot residuals
        if "rmin" in limits:
            limits_resid = {"min": limits["rmin"], "max": limits["rmax"]}
            if "cb_extend" in limits:
                limits_resid["cb_extend"] = limits["cb_extend"]
        else:
            limits_resid = [0, 100]
        fig_plot = self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[2],
            data=self.intensity - model,
            y_label=None,
            limits=limits_resid,
            colourmap="residuals-blanaced",
            # location=location,
            plot_type=plot_type,
            orientation=orientation,
        )
        if fit_centroid is not None:
            for i in range(len(fit_centroid[1])):
                if orientation == "horizontal":
                    ax[2].plot(
                        fit_centroid[0], fit_centroid[1][i], "k--", linewidth=0.5
                    )
                else:
                    ax[2].plot(
                        fit_centroid[1][i], fit_centroid[0], "k--", linewidth=0.5
                    )
        if orientation == "horizontal":
            ax[2].set_title("Residuals", loc=loc, y=0.75)
        else:
            ax[2].set_title("Residuals")
        locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)

        # organise the axes and labelling.
        if tight == True:
            if orientation == "horizontal":
                bottom0, top0 = ax[0].get_ylim()
                bottom1, top1 = ax[1].get_ylim()
                bottom2, top2 = ax[2].get_ylim()
                top_max = np.max([top0, top1, top2])
                bottom_min = np.min([bottom0, bottom1, bottom2])
                ax[0].set(xticklabels=[])
                ax[1].set(xticklabels=[])
                if len(ax) > 1:
                    ax[0].set_title(ax[0].get_title(), y=0.8)
                    ax[1].set_title(ax[1].get_title(), y=0.7)
                    ax[2].set_title(ax[2].get_title(), y=0.9)
                # fig_plot.rcParams['axes.titley'] = 1.0    # y is in axes-relative coordinates.
                # fig_plot.rcParams['axes.titlepad'] = -14  # pad is in points...
                ax[2].set_xlabel(f"{self.Azimuthlabel} ({self.AzimuthUnits})")
            else:
                bottom0, top0 = ax[0].get_ylim()
                bottom1, top1 = ax[1].get_ylim()
                bottom2, top2 = ax[2].get_ylim()
                top_max = np.max([top0, top1, top2])
                bottom_min = np.min([bottom0, bottom1, bottom2])
                # ax[0].set_ylim(top=top_max, bottom=bottom_min)
                # ax[1].set_ylim(top=top_max, bottom=bottom_min)
                # ax[2].set_ylim(top=top_max, bottom=bottom_min)
                ax[1].set(yticklabels=[])
                ax[2].set(yticklabels=[])
                ax[1].set(ylabel=None)
                ax[2].set(ylabel=None)
                ax[0].yaxis.set_ticks_position("both")
                ax[1].yaxis.set_ticks_position("both")
                ax[2].yaxis.set_ticks_position("both")
                ax[0].xaxis.set_ticks_position("both")
                ax[1].xaxis.set_ticks_position("both")
                ax[2].xaxis.set_ticks_position("both")

        # tidy layout
        plt.tight_layout()
        
        return fig_plot

    def plot_collected(
        self,
        fig_plot=None,
        axis_plot=None,
        show="intensity",
        colourmap=default_colourmap,
        limits=[0.01, 99.9],
        location="default",
        cbar_axes = None,
        debug=False,
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
            DESCRIPTION. The default is the default_colourmap, magma_r.
        limits : TYPE, optional
            DESCRIPTION. The default is [0, 99.9].
        location : TYPE, optional
            DESCRIPTION. The default is 'bottom'.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """

        if axis_plot == None and fig_plot == None:
            # make a figure
            fig_plot, axis_plot = plt.subplots()
            display = True

        elif axis_plot == None:
            # figure exists
            # assume will be shown outside of this function
            raise ValueError("Axes is needed to plot the data into")
            # FIX ME: the figure exists and so we whould be ableto find the axes

        elif fig_plot == None:
            # axis exists but figure not refrenced.
            # assume will be shown outside of this function
            raise ValueError("A figure is needed to plot the data into")
            # FIX ME: the axes exists and so we whould be ableto find the figure

        if isinstance(show, str) and show == "unmasked_intensity":
            plot_i = self.intensity.data
        elif isinstance(show, str) and show == "mask":
            plot_i = np.array(ma.getmaskarray(self.intensity), dtype="uint8") + 1
            colourmap = "Greys"
            limits = [0, 2]
        elif isinstance(show, np.ndarray) or ma.isMaskedArray(show):
            plot_i = show
        else:  # if show == "intensity"
            plot_i = self.intensity

        if np.ndim(plot_i) > 2:
            # The parametrized images to be plotted
            def f(data, n):
                return data[n]
        else:
            # there is only 1 plane of data.
            # (single detector rather than multiple detectors or positions.)
            def f(data, n):
                return data

        if isinstance(limits, dict):
            IMax = limits["max"]
            IMin = limits["min"]
            if "cb_extend" in limits:
                cb_extend = limits["cb_extend"]
            else:
                cb_extend = "neither"
        else:
            if limits[1] == 100:
                IMax = np.max(plot_i)
            else:
                IMax = np.nanpercentile(ma.filled(plot_i, np.nan), limits[1])
            if limits[0] == 0:
                IMin = np.min(plot_i)
            else:
                IMin = np.nanpercentile(ma.filled(plot_i, np.nan), limits[0])

            if IMin > 0 and IMax < 100:
                cb_extend = "both"
            elif IMax < 100:
                cb_extend = "max"
            elif IMin > 0:
                cb_extend = "min"
            else:
                cb_extend = "neither"

        im_num = 0
        the_plot = axis_plot.imshow(
            f(plot_i, im_num), vmin=IMin, vmax=IMax, cmap=colourmap
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

        if cbar_axes is not False:
            if cbar_axes is None:
                cbar_axes = axis_plot
            try:
                shrink = 0.9 / len(cbar_axes)
            except:
                shrink = 0.5
            cb = fig_plot.colorbar(
                mappable=the_plot,
                ax=axis_plot,
                extend=cb_extend,
                fraction=fraction,
                location=location,
                shrink=shrink,
            )

        if np.ndim(plot_i) > 2 and plot_i.shape[0] > 1:
            # adjust the main plot to make room for the sliders
            fig_plot.subplots_adjust(bottom=0.25)

            # Make a horizontal slider to control the image being viewed.
            axim = fig_plot.add_axes([0.25, 0.1, 0.65, 0.03])
            im_slider = plt.Slider(
                ax=axim,
                label="Frame number",
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

        # if display==True:
        #    the_plot.show()

    def plot_calibrated(
        self,
        fig_plot=None,
        axis_plot=None,
        show="default",
        x_axis="default",
        y_axis="default",
        # x_label = r"2$\theta$ ($^\circ$)",
        y_label = True, #"Azimuth ($^\circ$)",
        data=None,
        limits=[1, 99.9],
        y_lims=None,
        colourmap= default_colourmap,
        plot_type=False,
        point_scale=2,
        resample_shape=None,
        location=None,
        orientation="vertical",
        cbar_axes = None,
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """

        y_ticks = None
        x_ticks = None

        if axis_plot == None and fig_plot == None:
            # make a figure
            fig_plot, axis_plot = plt.subplots()

        elif axis_plot == None:
            # figure exists
            pass

        elif fig_plot == None:
            # axis exists but figure not refrenced.
            pass

        plot_type = self.what_plot_type(plot_type = plot_type)
        # if self.intensity.size > 1e5:
        #     rastered = True
        # elif (
        #     rastered == False and self.intensity.size > 1e6
        # ):  # 100000000000:# 1000000: # was 50000 until fixed max/min functions
        #     logger.moreinfo(" Have patience. The plot(s) will appear but it can take its time to render.")
        #     rastered = True

        if x_axis == "azimuth":
            plot_x = self.azm
            x_lims = [self.azm_start, self.azm_end]
        else:  # if x_axis is "default" or "tth"
            plot_x = self.tth
        # plot_y = self.azm
        label_x = f"{self.Dispersionlabel} ({self.DispersionUnits})"

        if y_axis == "intensity":
            # plot y rather than azimuth on the y axis
            plot_y = self.intensity
            # organise colour scale as azimuth
            plot_i = self.azm
            label_y = f"{self.Observationslabel} ({self.ObservationsUnits})"
            y_ticks = False
            plot_type = False
            # sort the data in reverse order of azimuth
            o = plot_i.argsort()[::-1]
            plot_i = plot_i[o]
            plot_y = plot_y[o]
            plot_x = plot_x[o]
            
        else:  # if y_axis is "default" or "azimuth"
            plot_y = self.azm
            plot_i = self.intensity
            if y_label == None:
                label_y = None# label_y = y_label
            else:
                label_y = f"{self.Azimuthlabel} ({self.AzimuthUnits})"
            if y_lims == None:
                y_lims = [self.azm_start, self.azm_end]
                # y_lims = [self.azm.min(), self.azm.max()]
            # axis_plot.set_ylim(y_lims)
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

        # set axis limits
        x_lims = [plot_x.min(), plot_x.max()]

        # set colour bar and colour maps.
        if colourmap == "Greys":
            IMax = 2.01
            IMin = 0
            cb_extend = "neither"
        elif isinstance(limits, dict):
            IMax = limits["max"]
            IMin = limits["min"]
            if "cb_extend" in limits:
                cb_extend = limits["cb_extend"]
            else:
                cb_extend = "neither"
                
            if isinstance(IMax, str) and isinstance(IMin, str):
                IMax = limits["max"]
                IMin = limits["min"]
                cb_extend = "both"
            elif isinstance(IMin, str):
                cb_extend = "max"
            elif isinstance(IMax, str):
                cb_extend = "min"
            
            
        else:
            if limits[1] == 100:
                IMax = np.max(plot_i)
            else:
                IMax = np.nanpercentile(ma.filled(plot_i, np.nan), limits[1])
            if limits[0] == 0:
                IMin = np.min(plot_i)
            else:
                IMin = np.nanpercentile(ma.filled(plot_i, np.nan), limits[0])
            if limits[0] > 0 and limits[1] < 100:
                cb_extend = "both"
            elif limits[1] < 100:
                cb_extend = "max"
            elif limits[0] > 0:
                cb_extend = "min"
            else:
                cb_extend = "neither"
        if np.isnan(IMax) and np.isnan(IMin):
            IMax = 1
            IMin = -1
        elif IMax == IMin:
            IMax += 1
            IMin -= 1
        if location == None:
            if y_lims is None:
                location = "right"
            elif plot_y.ndim == 1:
                x_range = x_lims[1] - x_lims[0]
                y_range = y_lims[1] - y_lims[0]
                if x_range >= y_range:
                    location = "bottom"
                    fraction = 0.046
                else:
                    location = "right"
                    fraction = 0.15
            elif plot_i.squeeze().shape[-1] > plot_i.squeeze().shape[-2]:
                location = "bottom"
                fraction = 0.046
            else:
                location = "right"
                fraction = 0.15
        elif location == "default":
            location == "right"
        else:
            pass

        # set colour map
        if colourmap == "residuals-blanaced":
            # colourmap = self.residuals_colour_scheme(
            colourmap = residuals_colour_scheme(IMax, IMin)
        else:
            colourmap = colourmap

        # if horizontal swap everything around.
        if orientation.lower() == "horizontal":
            plot_y, plot_x = plot_x, plot_y
            x_lims, y_lims = y_lims, x_lims
            label_x, label_y = label_y, label_x

            location = "right"
        else:
            location = "bottom"

        if plot_type == False or plot_type == "scatter":
            the_plot = axis_plot.scatter(
                plot_x,
                plot_y,
                s=1.5,
                c=plot_i,
                edgecolors="none",
                cmap=colourmap,
                vmin=IMin,
                vmax=IMax,
            )
        elif plot_type == "surf" or plot_type == "surface":
            if hasattr(self, 'convert_tth_azm_to_x_y'):
                tth_azm2x_y = self.convert_tth_azm_to_x_y
            else:
                tth_azm2x_y = None
            the_plot = surface_plot(
                plot_i,
                plot_x,
                plot_y,
                fig_plot=fig_plot,
                axis_plot=axis_plot,
                vmin=IMin,
                vmax=IMax,
                colourmap=colourmap,
                triangle_cutoff = 99,
                tth_azm2x_y = tth_azm2x_y
            )
        elif plot_type == "image" or plot_type == "im":
            the_plot = image_plot(
                plot_i,
                plot_x,
                plot_y,
                show = show,
                fig_plot=fig_plot,
                axis_plot=axis_plot,
                vmin=IMin,
                vmax=IMax,
                colourmap=colourmap,
            )
        elif plot_type == True or plot_type == "rastered" or plot_type == "rast":
            the_plot = raster_plot(
                plot_i,
                plot_x,
                plot_y,
                fig_plot=fig_plot,
                axis_plot=axis_plot,
                vmin=IMin,
                vmax=IMax,
                colourmap=colourmap,
                pixels_per_bin=point_scale,
                resample_shape=resample_shape,
            )
        else:
            raise ValueError(f"The plot type '{plot_type}' is not recognised.")

        axis_plot.set_xlabel(label_x)
        axis_plot.set_ylabel(label_y)

        axis_plot.set_xlim(x_lims)
        axis_plot.set_ylim(y_lims)

        if y_axis != "intensity":
            if orientation == "horizontal":
                axis_plot.set_xticks(y_ticks)
            else:
                axis_plot.set_yticks(y_ticks)

        # fix colour bar. 
        # cbar_axes = False --> dont have colour bar
        # cbar_axes = None --> cbar for these axes (default)
        # cbar_axes = Axis --> make cbar for this/these axes. Used to make 
        # single colour bar for data and model in self.plot_fitted.
        if cbar_axes is not False:
            if cbar_axes is None:
                cbar_axes = axis_plot
            try:
                shrink = 0.9 / len(cbar_axes)
            except:
                shrink = 0.9
            cb = fig_plot.colorbar(
                mappable=the_plot,
                ax=cbar_axes,
                extend=cb_extend,
                location=location,
                shrink=shrink,
            )  # , pad=0.1, aspect=8)

        return fig_plot

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
    else:  # if np.abs(maximum_value) < np.abs(minimum_value):
        keep = np.int_(
            (
                (2 * np.abs(minimum_value) - (np.abs(minimum_value) - maximum_value))
                / (2 * np.abs(minimum_value))
            )
            * n_entries
        )
        all_colours = all_colours[:keep]
    all_colours = colors.ListedColormap(
        all_colours, name="myColorMap", N=all_colours.shape[0]
    )

    return all_colours


default_colourmap = _Plot_AngleDispersive.default_colourmap

def raster_plot(
    data_plot,
    x_plot,
    y_plot,
    fig_plot=None,
    axis_plot=None,
    resample_shape=None,
    vmin=0,
    vmax=np.inf,
    colourmap=default_colourmap,
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
        Colourmap for the plot. The default is the default_colourmap, magma_r.
    pixels_per_bin : float, optional
        Scaler for the number of bins in the histogram  . The default is 3.

    Returns
    -------
    pl : axes
        filled set of axes.
    """
    # FIX ME: for massive data sets (e.g. ESFR) this is still slow.
    # maybe it should be replaced by datashader or something else that is good for dynamicaly plotting massive data sets.
    # an alternative is to use GIS image referencing packages.

    if logger.is_below_level(level="DEBUG"):
        import time

        start = time.time()

    if resample_shape != None:
        # use the size of the output set by resample_shape
        x_bins = resample_shape[0]
        y_bins = resample_shape[1]
    else:
        # base the size of the output on the number of pixels in the axes.
        if axis_plot == None:
            if fig_plot == None:
                fig = plt.figure()
            axis_plot = fig.add_subplot(1, 1, 1)
        bbox = axis_plot.get_window_extent().transformed(
            fig_plot.dpi_scale_trans.inverted()
        )
        width, height = bbox.width, bbox.height
        width *= fig_plot.dpi
        height *= fig_plot.dpi
        x_bins = width / pixels_per_bin
        y_bins = height / pixels_per_bin

    # FIX ME: should I just replace this with a call to histogram2d_engine in pyFAI/engines/histogram_engine.py??
    # To have empyt pixels where there is no data use dummy = np.nan as an option.
    # see comments in histogram2d for answer
    result, x_edges, y_edges, num_pix_per_bin = histogram2d(
        data_plot, x_plot, y_plot, x_bins=x_bins, y_bins=y_bins
    )
    result = np.rot90(result)
    pl = axis_plot.imshow(
        result,
        extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
        cmap=colourmap,
    )

    # axis_plot.invert_yaxis()
    if logger.is_below_level(level="DEBUG"):
        end = time.time()
        logger.debug(" ".join(map(str, [(f"Time to make rastered plot {end-start}")])))

    return pl


def surface_plot(
    data_plot,
    x_plot,
    y_plot,
    fig_plot=None,
    axis_plot=None,
    resample_shape=None,
    vmin=0,
    vmax=np.inf,
    colourmap=default_colourmap,
    triangle_cutoff = 98,
    tth_azm2x_y = None
):
    """
    Plots the data on an irregular tripcolor gird.

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
        Colourmap for the plot. The default is the default_colourmap, magma_r.
    triangle_cutoff : float, optional
        Perceltile threshold for filtering the triangles. The default is 0.99

    Returns
    -------
    pl : axes
        filled set of axes.
    """

    if ma.is_masked(data_plot):
        data_plot = data_plot.compressed()
    elif not ma.is_masked(data_plot) and ma.is_masked(x_plot):
        # if the data is not masked while the axes are if causes promblems when plotting. 
        # so mask.
        data_plot = ma.array(data_plot, mask=x_plot.mask).compressed()
    if ma.is_masked(x_plot):
        x_plot = x_plot.compressed()
    if ma.is_masked(y_plot):
        y_plot = y_plot.compressed()
    
    if tth_azm2x_y:
        x_physical, y_physical = tth_azm2x_y(x_plot, y_plot)
    else:
        x_physical, y_physical  = x_plot, y_plot
        
    # rad = 1 * np.tan(np.deg2rad(x_plot.compressed()))
    # x_physical = rad * np.cos(np.deg2rad(ma.array(y_plot).compressed()))
    # y_physical = rad * np.sin(np.deg2rad(ma.array(y_plot).compressed()))
    # data_plot = ma.array(data_plot).compressed()
        
    triang = tri.Triangulation(x_physical, y_physical)
    corners = triang.triangles
    
    if triangle_cutoff != 100:
        areas = []
        x_range = []
        y_range = []

        import proglog
        progress_bar = proglog.default_bar_logger(
            "bar"
        )  # shorthand to generate a bar logger

        def PolyArea(x, y):
            return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

        # calculate size of triangles and discard ones taht are too large. 
        for i in range(len(corners)):
            areas.append(
                PolyArea(x_physical.flatten()[corners][i], y_physical.flatten()[corners][i])
            )
            x_range.append(
                x_physical.flatten()[triang.triangles][i].max()
                - x_physical.flatten()[triang.triangles][i].min()
            )
            y_range.append(
                y_physical.flatten()[triang.triangles][i].max()
                - y_physical.flatten()[triang.triangles][i].min()
            )
        
        areas = ma.array(areas)
        x_range = ma.array(x_range)
        y_range = ma.array(y_range)
        
        # Define an additional condition to mask elements greater than 80
        additional_condition1 = areas <= 0
        additional_condition2 = x_range <= 0
        additional_condition3 = y_range <= 0
        # Update the mask to include additional elements
        areas.mask = areas.mask | additional_condition1 | additional_condition2 | additional_condition3   
        x_range.mask = x_range.mask | additional_condition1 | additional_condition2 | additional_condition3   
        y_range.mask = y_range.mask | additional_condition1 | additional_condition2 | additional_condition3    
            
        if 0:
            fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)
            # We can set the number of bins with the *bins* keyword argument.
            axs[0].hist(np.log10(ma.array(areas).compressed()), bins=int(len(areas)/100))
            axs[1].hist(np.log10(ma.array(x_range).compressed()), bins= int(len(areas)/100))
            axs[2].hist(np.log10(ma.array(y_range).compressed()), bins= int(len(areas)/100))
            plt.show()
        
        if 0:
            cutoff_area = np.nanpercentile(ma.array(areas).compressed(), triangle_cutoff)
            cutoff_x = np.nanpercentile(ma.array(x_range).compressed(), triangle_cutoff)
            cutoff_y = np.nanpercentile(ma.array(y_range).compressed(), triangle_cutoff)
        else:
            multiples_of_median = 3
            cutoff_area = np.nanmedian(ma.array(areas).compressed())*multiples_of_median
            cutoff_x = np.nanmedian(ma.array(x_range).compressed())*multiples_of_median
            cutoff_y = np.nanmedian(ma.array(y_range).compressed())*multiples_of_median
        
        keep = []
        for i in progress_bar.iter_bar(FilterTriangles=range(len(corners))):
            # keep well behaved triangles
            if (np.abs(x_range[i]) <= cutoff_x
                and np.abs(y_range[i]) <= cutoff_y
                and np.abs(areas[i]) <= cutoff_area):
                keep.append(i)    

        triang.triangles = triang.triangles[keep]
        corners = corners[keep]

    # replce the coordinates of the triangles with coordinates to plot.
    triang.x = x_plot
    triang.y = y_plot

    pl = axis_plot.tripcolor(
        triang, data_plot, cmap=colourmap, vmin=vmin, vmax=vmax,
        shading='gouraud'
    )
    return pl   


def image_plot(
    data_plot,
    x_plot,
    y_plot,
    show = "intensity",
    fig_plot=None,
    axis_plot=None,
    resample_shape=None,
    vmin=0,
    vmax=np.inf,
    colourmap=default_colourmap
):
    """
    Plots the data as an image. Requires that the x_plot and y_plot are image coordinates.
    
    This is the fastest way to plot the data but is only representative if the underlying data is 
    on an orthonormal grid.

    Parameters
    ----------
    data_plot : masked array
        data to convert into the image to plot.
    x_plot : masked array
        horizontal position of the data points (usually two theta)
    y_plot : masked array
        vertical position of the data points (usually azimuth)
    show : bool, str 
        How to plot data_plot. The default is "intensity".
        If not 'intensity' or 'default' then any 'bad' data is shown as white. 
    fig_plot : figure, optional
        Figure to add plot to. The default is None.
    axis_plot : axes, optional
        Axes to add plot to. The default is None.
    vmin : float, optional
        minimum of the plotted colour scale. The default is 0.
    vmax : float, optional
        minimum of the plotted colour scale. The default is np.inf, in effect the maximum value in data_plot
    colourmap : string, optional
        Colourmap for the plot. The default is the default_colourmap, magma_r.

    Returns
    -------
    pl : axes
        filled set of axes.
    """

    # FIXME: should test that the data is convertable/plottable as an image. 
    # otherwise revert to another plotting type.

    # convert data to an image array.
    rounding_precision = 8
    # if the y values are nor rounded then numberical precision can make the number of unique values 
    # greater than the actual number. Round to remove imprecision. Also correctly identifies the right number of 
    # values then. 
    # FIXME: rounding_precision is a fixed value. should be set from the precision of the data.
    x_unique_vals, x_inverse = np.unique(np.round(x_plot.data,rounding_precision), return_inverse=True)
    y_unique_vals, y_inverse = np.unique(np.round(y_plot.data,rounding_precision), return_inverse=True)

    data_im = ma.zeros([y_unique_vals.shape[0], x_unique_vals.shape[0]])-1
    data_im[y_inverse, x_inverse] = data_plot
    
    if show == "intensity" or show == "default":
        data_im[data_im==-1] = np.nan
        data_im[data_im==-np.inf] = np.nan
        data_im[data_im==np.inf] = np.nan
    cmap = colormaps.get_cmap(colourmap)
    cmap.set_bad(color='white', alpha = 1.0)
    
    pl = axis_plot.imshow(
        data_im, cmap=cmap, vmin=vmin, vmax=vmax,
        origin='lower',
        aspect='auto',
        extent=[np.nanmin(x_unique_vals), np.nanmax(x_unique_vals), np.nanmin(y_unique_vals), np.nanmax(y_unique_vals)]
    )

    return pl
