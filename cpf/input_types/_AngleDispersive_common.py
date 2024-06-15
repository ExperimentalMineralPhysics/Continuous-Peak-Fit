#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


class _AngleDispersive_common():
    """
    Common attributes and methods for the angle dispersive diffraction classes.
    They are in here to save replicating code.

    If the methods required by a class are different then they are included within the 
    class Functions file. 

    This class cannot be imported as a stand alone class. Instead the methods 
    contained within it are created as methods in the angle dispersive diffraction
    classes (e.g. DioptasDetector and ESRFlvpDetector).    
    """    

    def __init__(self, detector_class=None):
        """
        :param detector_class:
        """
        self = detector_class


    

    def bins(self, orders_class, cascade=False):
        """
        Determine bins to use in initial fitting.
        Assign each data to a chunk corresponding to its azimuth value
        Returns array with indices for each bin and array of bin centroids
        :param orders_class:
        :return chunks:
        :return bin_mean_azi:
        """

        # determine how to divide the data into bins and how many.
        if cascade:
            bt = orders_class.cascade_bin_type
            if bt == 1:
                b_num = orders_class.cascade_number_bins
            else:
                b_num = orders_class.cascade_per_bin
        else:
            bt = orders_class.fit_bin_type
            if bt == 1:
                b_num = orders_class.fit_number_bins
            else:
                b_num = orders_class.fit_per_bin

        # make the bins
        if bt == 0:
            # split the data into bins with an approximately constant number of data.
            # uses b_num to determine bin size
            num_bins = int(np.round(len(self.azm[self.azm.mask == False]) / b_num))
            bin_boundaries = self.equalObs(
                np.sort(self.azm[self.azm.mask == False]), num_bins
            )
        elif bt == 1:
            # split the data into a fixed number of bins
            # uses b_num to determine bin size
            lims = np.array(
                [
                    np.min(self.azm[self.azm.mask == False]),
                    np.max(self.azm[self.azm.mask == False] + 0.01),
                ]
            )
            lims = np.around(lims / 45) * 45
            bin_boundaries = np.linspace(lims[0], lims[1], num=b_num + 1)

        else:
            # issue an error
            err_str = "The bin type is not recognised. Check input file."
            raise ValueError(err_str)

        if 0:  # for debugging
            print(bin_boundaries)
            # print(np.sort(bin_boundaries))
            # create histogram with equal-frequency bins
            n, bins, patches = plt.hist(
                self.azm[self.azm.mask == False], bin_boundaries, edgecolor="black"
            )
            plt.show()
            # display bin boundaries and frequency per bin
            print("bins and occupancy", bins, n)
            print("expected number of data per bin", orders_class.cascade_per_bin)
            print("total data", np.sum(n))

        # fit the data to the bins
        chunks = []
        bin_mean_azi = []
        temp_azimuth = self.azm.flatten()
        for i in range(len(bin_boundaries) - 1):
            start = bin_boundaries[i]
            end = bin_boundaries[i + 1]
            azi_chunk = np.where((temp_azimuth > start) & (temp_azimuth <= end))
            chunks.append(azi_chunk)
            bin_mean_azi.append(np.mean(temp_azimuth[azi_chunk]))

        return chunks, bin_mean_azi




    def set_limits(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf]):
        """
        Cut the data to only data within range_bounds (two theta) and azi_bounds (azimuth).
        
        N.B. This is not a masking function. it removes all the data outside of the range.
        Data inside the range that is masked remains.

        Parameters
        ----------
        range_bounds : list, optional
            Two theta limits to apply to the data. The default is [-np.inf, np.inf].
        azi_bounds : list, optional
            Azimuth limits to apply to the data. . The default is [-np.inf, np.inf].

        Returns
        -------
        None.

        """
        local_mask = np.where(
            (self.tth >= range_bounds[0]) & 
            (self.tth <= range_bounds[1]) &
            (self.azm >= azi_bounds[0]) & 
            (self.azm <= azi_bounds[1])
        )
        self.intensity = self.intensity[local_mask]
        self.tth = self.tth[local_mask]
        self.azm = self.azm[local_mask]
        if "dspace" in self.__dict__.items():
            self.dspace = self.dspace[local_mask]
        
        if "x" in self.__dict__.items():
            self.x = self.x[local_mask]
        if "y" in self.__dict__.items():         
            self.y = self.y[local_mask]
        if "z" in self.__dict__.items():
            self.z = self.z[local_mask]
            
        self.azm_end = np.max(self.azm)
        self.azm_start = np.min(self.azm)
        
        
        
        # FIX ME This appears to be unused and has been removed. 
    # def test_azims(self, steps = 360):
    #     """
    #     Returns equally spaced set of aximuths within possible range.

    #     Parameters
    #     ----------
    #     steps : Int, optional
    #         DESCRIPTION. The default is 360.

    #     Returns
    #     -------
    #     array
    #         list of possible azimuths.

    #     """
        
    #     return np.linspace(self.azm_start,self.azm_end, steps+1)




def equalObs(x, nbin):
    """
    get equally populated bins for data set.
    copied from: https://www.statology.org/equal-frequency-binning-python/ on 26th May 2022.

    Parameters
    ----------
    x : TYPE
        data to disperse.
    nbin : TYPE
        number of bins.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    nlen = len(x)
    x = np.sort(x)
    return np.interp(np.linspace(0, nlen, nbin + 1), np.arange(nlen), np.sort(x))

