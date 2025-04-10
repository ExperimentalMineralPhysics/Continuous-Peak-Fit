#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from matplotlib import cm, colors, gridspec, tri

from cpf.histograms import histogram1d, histogram2d
from cpf.logging import CPFLogger

logger = CPFLogger("cpf.input_types._Plot_AngleDispersive")



class _Plot_AngleDispersive:
    """
    This class cannot be imported as a stand alone class. Instead the methods
    contained within it are created as methods in the angle dispersive diffraction
    classes (e.g. DioptasDetector and SERFlvpDetector).
    """

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
    
    
    
    def plot_integrated(
        self,
        fig_plot=None,
        axis_plot=None, show=None):
        """
        Makes a plot of the integrated data. 
        
        :param fig:
        :return:
        """
    
        if fig_plot == None:
            #make a figure
            fig_plot, axis_plot = plt.subplots()
    
        p, i, a = histogram1d(self.tth, self.intensity, self.azm)

        axis_plot.plot(p,i)
                
        axis_plot.set_xlabel(f"{self.Dispersionlabel} ({self.DispersionUnits})")
        axis_plot.set_ylabel("Intensity (a.u.)")
        axis_plot.set_title("Integrated Data")
    
    
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
        plot_type="surface",
        orientation="horizontal",
        plot_ColourRange=None,
    ):
        """
        add data to axes.
        :param ax:
        :param show:
        :return:
        """

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
        ax = []
        if tight == True:
            gs = gridspec.GridSpec(up, across, wspace=0.0, hspace=0.0)
            for i in range(3):
                ax.append(fig_plot.add_subplot(gs[i]))
        else:
            axs = fig_plot.subplots(nrows=up, ncols=across, sharey=True, sharex=True)
            ax.append(axs[0])
            ax.append(axs[1])
            ax.append(axs[2])

        # plot data
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[0],
            show="intensity",
            limits=limits,
            colourmap="magma_r",
            # location=location,
            rastered=plot_type,
            orientation=orientation,
        )
        if orientation == "horizontal":
            ax[0].set_title("Data", loc=loc, y=0.75)
        else:
            ax[0].set_title("Data")
        locs, labels = plt.xticks()
        # plt.setp(labels, rotation=90)

        # plot model
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[1],
            data=model,
            y_label=None,
            limits=limits,
            colourmap="magma_r",
            # location=location,
            rastered=plot_type,
            orientation=orientation,
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
        else:
            limits_resid = [0, 100]
        self.plot_calibrated(
            fig_plot=fig_plot,
            axis_plot=ax[2],
            data=self.intensity - model,
            y_label=None,
            limits=limits_resid,
            colourmap="residuals-blanaced",
            # location=location,
            rastered=plot_type,
            orientation=orientation,
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
                ax[2].set_xlabel("Azimuth ($^\circ$)")
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

    def plot_collected(
        self,
        fig_plot=None,
        axis_plot=None,
        show="intensity",
        colourmap="jet",
        limits=[0, 99.9],
        location="default",
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
    
        IMax = np.percentile(plot_i.compressed(), limits[1])
        IMin = np.percentile(plot_i.compressed(), limits[0])
        
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

        cb = fig_plot.colorbar(
            mappable=the_plot,
            ax=axis_plot,
            extend=cb_extend,
            fraction=fraction,
            location=location,
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
        y_label = "Azimuth ($^\circ$)",
        data=None,
        limits=[1, 99.9],
        y_lims=None,
        colourmap="jet",
        rastered=False,
        point_scale=2,
        resample_shape=None,
        location=None,
        orientation="vertical",
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

        if (
            rastered == False and self.intensity.size > 1e6
        ):  # 100000000000:# 1000000: # was 50000 until fixed max/min functions
            logger.moreinfo(
                " ".join(
                    map(
                        str,
                        [
                            (
                                " Have patience. The plot(s) will appear but it can take its time to render."
                            )
                        ],
                    )
                )
            )
            rastered = True

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
            label_y = "Intensity (a.u.)"
            y_ticks = False
        else:  # if y_axis is "default" or "azimuth"
            plot_y = self.azm
            plot_i = self.intensity
            label_y = y_label
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
            cb_extend = "neither"
        else:
            if limits[1] == 100:
                IMax = np.max(plot_i)
            else:
                IMax = np.percentile(plot_i.compressed(), limits[1])
            if limits[0] == 0:
                IMin = np.min(plot_i)
            else:
                IMin = np.percentile(plot_i.compressed(), limits[0])
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
            IMax+=1
            IMin-=1
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

        if rastered == False or rastered == "scatter":
            the_plot = axis_plot.scatter(
                plot_x,
                plot_y,
                s=1.5,
                c=plot_i,
                edgecolors="none",
                cmap=colourmap,
                vmin=IMin,
                vmax=IMax,
                # rasterized=rastered,
            )
        elif rastered == "surf" or rastered == "surface":
            the_plot = surface_plot(
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

        axis_plot.set_xlabel(label_x)
        axis_plot.set_ylabel(label_y)

        axis_plot.set_xlim(x_lims)
        axis_plot.set_ylim(y_lims)

        if y_axis != "intensity":
            if orientation == "horizontal":
                axis_plot.set_xticks(y_ticks)
            else:
                axis_plot.set_yticks(y_ticks)

        # p2 = axis_plot.get_position().get_points().flatten()
        # ax_cbar1 = fig_plot.add_axes([p2[0], 0, p2[2]-p2[0], 0.05])

        cb = fig_plot.colorbar(
            mappable=the_plot,
            ax=axis_plot,
            extend=cb_extend,
            location=location,
            shrink=0.9,
        )  # , pad=0.1, aspect=8)


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


def raster_plot(
    # self,
    data_plot,
    x_plot,
    y_plot,
    fig_plot=None,
    axis_plot=None,
    resample_shape=None,
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

    if logger.is_below_level(level="DEBUG"):
        import time

        start = time.time()

    if resample_shape != None:
        # use the size of the output set by resample_shape
        x_bins = resample_shape[0]
        y_bins = resample_shape[1]
    else:
        # base the sixe of the output on the number of pixels in the axes.
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
    # self,
    data_plot,
    x_plot,
    y_plot,
    fig_plot=None,
    axis_plot=None,
    resample_shape=None,
    vmin=0,
    vmax=np.inf,
    colourmap="jet",
    pixels_per_bin=3,
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

    triang = tri.Triangulation(x_plot.flatten(), y_plot.flatten())
    corners = triang.triangles
    # FIXME: i might be better to make the triangles from the x,y, points rather than the tth+azm.
    # when made like this the plots spread the intensity in a broadway. -- the scatter plot looks better.

    triang.set_mask(tri.TriAnalyzer(triang).get_flat_tri_mask())
    mask = np.zeros(len(triang.triangles))
    for i in range(len(corners)):
        if any(x_plot.flatten()[corners[i, :]].mask == True):
            mask[i] = True
    mask = np.array(mask, dtype="bool")
    triang.set_mask((triang.mask == True) | (mask == True))

    if 0:
        areas = []
        x_range = []
        y_range = []
        import proglog

        progress_bar = proglog.default_bar_logger(
            "bar"
        )  # shorthand to generate a bar logger

        def PolyArea(x, y):
            return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))

        # for f in range(settings_class.image_number):
        for i in progress_bar.iter_bar(iteration=range(len(corners))):
            if triang.mask[i] == True:
                pass
            else:
                areas.append(
                    PolyArea(x_plot.flatten()[corners][i], y_plot.flatten()[corners][i])
                )
                x_range.append(
                    x_plot.flatten()[triang.triangles][i].max()
                    - x_plot.flatten()[triang.triangles][i].min()
                )
                y_range.append(
                    y_plot.flatten()[triang.triangles][i].max()
                    - y_plot.flatten()[triang.triangles][i].min()
                )

        fig, axs = plt.subplots(1, 3, sharey=True, tight_layout=True)

        # We can set the number of bins with the *bins* keyword argument.
        axs[0].hist(areas, bins=50)  # int(len(areas)/40))
        axs[1].hist(x_range, bins=50)  # int(len(areas)/40))
        axs[2].hist(y_range, bins=50)  # int(len(areas)/40))

        # FIXME: I should filter the triangles on their area or their shape.
        # mean_area = np.mean(areas)

    pl = axis_plot.tripcolor(
        triang, data_plot.flatten(), cmap=colourmap, vmin=vmin, vmax=vmax
    )

    return pl
