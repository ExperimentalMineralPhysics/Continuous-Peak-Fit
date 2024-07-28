#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__all__ = ["histogram1d"]


import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pyFAI.engines.histogram_engine as he

# from cpf.XRD_FitPattern import logger
from cpf.logger_functions import logger


def histogram1d(
    tth,
    intensity,
    azi=None,
    dspace=None,
    histogram_type="data",
    bin_n=None,
    debug=False,
):
    """
    Integration function for the two theta, intensity pixel data.
    It can make constant width bins or bins with constant number of data in them.


    Parameters
    ----------
    tth : array or list
        two theha position of pixel data.
    intensity : array or list
        intensities of pixel data.
    azi : array or list, optional
        azimuths of the data.
    bin_type : string, optional
        Switch for how to bin the data. The options are:
            "width" - constant width in two theta; number of data in each bin varies
            "data"  - each bin has the same number of data in it. variable width.
            The default is "data".
     bin_n : number
        the number of bins to make (if constant width) or the number of data
        to put in each bin (if constant data).
    debug : True/False
        Output control

    Returns
    -------
    bin_positions, intensities, azimuths

    """

    if azi.all != None:
        if tth.size != intensity.size:
            raise ValueError("the intensity and two theta arrays are not the same size")
    else:
        if tth.size != intensity.size & tth.size != azi.size:
            raise ValueError(
                "the intensity, two theta and azimuth arrays are not the same size"
            )

    if bin_n == None:
        # set number of bins according to Sturges' Rule
        # see: https://en.wikipedia.org/wiki/Sturges%27_Rule
        # or the number of data.
        # the multiply by 4 for Sturges' Rule is arbitrary but needed here becuase the data is wider than the distributions.
        bin_n = np.max([1 + 3.322 * np.log10(tth.size) * 4, tth.size / 50])
    if not isinstance(bin_n, int):
        bin_n = int(np.round(bin_n))

    if histogram_type == "data":
        # sort the data and then bin accordingly.
        order = np.argsort(tth)
        tth = tth[order]
        intensity = intensity[order]
        if azi.all != None:
            azi = azi[order]

        mean_per_bin = order.size / bin_n

        # could use a use pyFAI named tuple data holder.
        # histogram = namedtuple("Integrate1dtpl", "position intensity sigma signal variance normalization count std sem norm_sq", defaults=(None,) * 3)
        position = []
        intens = []
        signal = []
        count = []
        azm = []
        dspc = []

        for i in range(np.int(bin_n)):
            data_start = int(np.round(order.size / np.round(bin_n) * i))
            data_end = int(np.round(order.size / np.round(bin_n) * (i + 1)))

            position.append(np.mean(tth[data_start:data_end]))
            intens.append(np.mean(intensity[data_start:data_end]))
            signal.append(np.sum(intensity[data_start:data_end]))
            count.append(np.size(tth[data_start:data_end]))
            if azi.all != None:
                azm.append(np.mean(azi[data_start:data_end]))

        if debug == True:
            plt.plot(tth, intensity, ".", position, intens, "-", position, count, "-")
            plt.show()

    elif histogram_type == "width":
        # use pyFAI 1D integration to make equal width histogram.
        # pyFAI function is here :  pyFAI/src/pyFAI/containers.py
        histogram = he.histogram1d_engine(tth, bin_n, intensity)
        # returns a pyFAI named tuple data holder.
        # histogram = namedtuple("Integrate1dtpl", "position intensity sigma signal variance normalization count std sem norm_sq", defaults=(None,) * 3)
        # histogram.position or is tth position
        # histogram.intensity is intensity
        # histogram.singal is total counts?
        # histogram.normalization is normalisation factor
        # histogram.count is number pixels in range.
        position = histogram.position
        intens = histogram.intensity

        # get the mean azimuth of each bin.
        histogram = he.histogram1d_engine(tth, bin_n, azi)
        azm = histogram.intensity

        if debug == True:
            plt.plot(
                tth,
                intensity,
                ".",
                histogram.position,
                histogram.intensity,
                "-",
                histogram.position,
                histogram.normalization,
                "-",
            )
            plt.show()

    else:
        raise ValueError(
            "the intensity, two theta and azimuth arrays are not the same size"
        )

    return np.array(position), np.array(intens), np.array(azm)


def histogram2d(data, x, y, x_bins=500, y_bins=720):
    """
    Reduce the diffraction pixel data into regualar gridded data (in effect an image).
    This is basically pyFAI's integrate2d function without all the bells and whistles.

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    x_bins : TYPE, optional
        DESCRIPTION. The default is 500.
    y_bins : TYPE, optional
        DESCRIPTION. The default is 720.

    Returns
    -------
    result : TYPE
        DESCRIPTION.
    x_edges : TYPE
        DESCRIPTION.
    y_edges : TYPE
        DESCRIPTION.
    num_pix_per_bin : TYPE
        DESCRIPTION.

    """
    # FIX ME: should I just replace this with a call to histogram2d_engine in pyFAI/engines/histogram_engine.py??
    # To have empyt pixels where there is no data use dummy = np.nan as an option.
    # (April 2024) I think that it should stay as a function -- even if it eventually calls the pyFAI function
    # we can forace all the options we want here rather than having to set them everytime.

    if not ma.isMaskedArray(data):
        data = ma.array(data)
    if not ma.isMaskedArray(x):
        x = ma.array(x)
    if not ma.isMaskedArray(y):
        y = ma.array(y)

    x_edges = np.linspace(
        x[x.mask == False].min(), x[x.mask == False].max(), int(x_bins) + 1
    )
    y_edges = np.linspace(
        y[y.mask == False].min(), y[y.mask == False].max(), int(y_bins) + 1
    )

    num_pix_per_bin, _, _ = np.histogram2d(
        x[x.mask == False].flatten(),
        y[y.mask == False].flatten(),
        bins=[x_edges, y_edges],
    )
    nominator, _, _ = np.histogram2d(
        x[x.mask == False].flatten(),
        y[y.mask == False].flatten(),
        bins=[x_edges, y_edges],
        weights=data[data.mask == False].flatten(),
    )
    result = nominator / num_pix_per_bin.clip(1)

    i, j = np.where(num_pix_per_bin == 0)
    result[i, j] = np.nan  # make a white (empty) pixel
    num_pix_per_bin[i, j] = np.nan  # make a white (empty) pixel

    return result, x_edges, y_edges, num_pix_per_bin
