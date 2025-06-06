#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma

from cpf.util.logging import get_logger
from skimage.transform import rescale, resize, downscale_local_mean

logger = get_logger("cpf.input_types._AngelDispersive_common")


class _AngleDispersive_common:
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

    def _get_d_space(self, mask=None):
        """
        Give d-spacing value for detector x,y position; calibration info in data
        Get as twotheta and then convert into d-spacing
        :param mask:
        :return:
        """
        # FIX ME: i dont know that this method is needed. It could be wrapped into
        # fill_data. Check for calls and if this is the only one then remove it
        if mask is not None:
            return ma.array(
                self.conversion(self.tth, reverse=False),
                mask=mask,
            )
        else:
            return ma.array(self.conversion(self.tth, reverse=False))

    def conversion(self, tth_in, azm=None, reverse=False):
        """
        Convert two theta values into d-spacing.
        azm is needed to enable compatibility with the energy dispersive detectors
        :param tth_in:
        :param reverse:
        :param azm:
        :return:
        """
        tth_in = np.array(tth_in)

        if isinstance(tth_in, list):
            tth_in = np.array(tth_in)
            a = True
        else:
            a = False

        # wavelength = self.calibration.wavelength * 1e10
        if self.conversion_constant is None:
            # thnen the conversion has not been initated proerly somewhere
            raise ValueError("'Conversion_constant' cannot be None. Something has not been iniated properly.")
        else:
            wavelength = self.conversion_constant
        

        if self.conversion_constant == None:
            dspc_out = tth_in
        elif not reverse:
            # convert tth to d_spacing
            dspc_out = wavelength / 2 / np.sin(np.radians(tth_in / 2))
        else:
            # convert d-spacing to tth.
            # N.B. this is the reverse function so that labels tth and d_spacing are not correct.
            dspc_out = 2 * np.degrees(np.arcsin(wavelength / 2 / tth_in))

        if a == True:
            dspc_out = list(dspc_out)
        return np.squeeze(np.array(dspc_out))

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
            if bt == None:
                # force a default
                bt = 0
                b_num = 50
            elif bt == 1:
                b_num = orders_class.cascade_number_bins
            else:
                b_num = orders_class.cascade_per_bin
        else:
            bt = orders_class.fit_bin_type
            if bt == None:
                # force a default
                bt = 1
                b_num = 90
            elif bt == 1:
                b_num = orders_class.fit_number_bins
            else:
                b_num = orders_class.fit_per_bin

        # make the bins
        if bt == 0:
            # split the data into bins with an approximately constant number of data.
            # uses b_num to determine bin size
            num_bins = int(np.round(len(self.azm[self.azm.mask == False]) / b_num))
            bin_boundaries = equalObs(
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

        if 0:
            # create histogram with equal-frequency bins
            n, bins, patches = plt.hist(
                self.azm[self.azm.mask == False], bin_boundaries, edgecolor="black"
            )
            plt.show()
            logger.debug(" ".join(map(str, [("bins and occupancy", bins, n)])))
            logger.debug(" ".join(map(str, [("total data", np.sum(n))])))
        # display bin boundaries and frequency per bin
        logger.debug(" ".join(map(str, [("bin boundaries:", bin_boundaries)])))
        if bt == 1:
            logger.debug(
                " ".join(
                    map(
                        str,
                        [("expected number of chunks", b_num)],
                    )
                )
            )
        else:
            logger.debug(
                " ".join(
                    map(
                        str,
                        [("expected number of data per chunkbin", b_num)],
                    )
                )
            )

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
            (self.tth >= range_bounds[0])
            & (self.tth <= range_bounds[1])
            & (self.azm >= azi_bounds[0])
            & (self.azm <= azi_bounds[1])
        )
        self.intensity = self.intensity[local_mask]
        self.tth = self.tth[local_mask]
        self.azm = self.azm[local_mask]
        if "dspace" in dir(self):
            self.dspace = self.dspace[local_mask]

        if "x" in dir(self):
            if self.x is not None:
                self.x = self.x[local_mask]
        if "y" in dir(self):
            if self.y is not None:
                self.y = self.y[local_mask]
        if "z" in dir(self):
            if self.z is not None:
                self.z = self.z[local_mask]

        # self.azm_end = np.max(self.azm)
        # self.azm_start = np.min(self.azm)

    def test_azims(self, steps=360):
        """
        Returns equally spaced set of aximuths within possible range.

        Parameters
        ----------
        steps : Int, optional
            DESCRIPTION. The default is 360.

        Returns
        -------
        array
            list of possible azimuths.

        """

        return np.linspace(self.azm_start, self.azm_end, steps + 1)

    def GetDataType(self, rawData, numType="float", minimumPrecision=32):
        """
        A function to return the data type for the diffraction data.
        The data can be any format when it is read in, but should be a float
        for the data fitting because the lmfit model inherits the data type from
        the data.


        Parameters
        ----------
        rawData : array
            Diffraction data array.
        numType : string, optional
            Data type to make, without numerical precision. The default is "float".
        minimumPrecision : number, optional
            Minimum numerical precision to be applied to the data. If the numerical
            precision does not matter, set to False. The default is 32.

        Returns
        -------
        DataType : string
            A numpy data type to be applied to the data arrays.

        """
        precision = re.findall("\d+", rawData.dtype.name)[0]
        # force minimum precision
        if minimumPrecision != False:
            if precision < minimumPrecision:
                precision = minimumPrecision
        else:
            pass

        try:
            DataType = np.dtype(numType + precision)
        except:
            err_str = f"The datatype, {DataType}, is not a recognised data type"
            logger.critical(err_str)
            raise TypeError(err_str)

        return DataType
    
    
    def duplicate_without_detector(self, range_bounds=[-np.inf, np.inf], azi_bounds=[-np.inf, np.inf]):
        """
        Makes an independent copy of a Intensity, two theta and azimuth data. 
        If present the d-spacing and detector x,y,z are included as well 

        range_bounds and azi_bounds restrict the extent of the data if needed.
        The range resturictions should be applied upon copying (if required) for memory efficieny.
        
        Alaternatively run:
        region = data_class.get_SubRegion()
        region.set_limit2(range_bounds=[...], azi_bounds=[...])

        Although a Detector instance is returned the detector and calibration are empty. 

        Parameters
        ----------
        range_bounds : dict or array, optional
            Limits for the two theta range. The default is [-np.inf, np.inf].
        azi_bounds : dict or array, optional
            Limits for the azimuth range. The default is [-np.inf, np.inf].
        reduction : int, optional
            A fraction by which to reduce the data. The reduction for a Dioptas 
            detector instance in by binning the data. 

        Returns
        -------
        new : Detector Instance.
            Detector with only image related arrays.

        """  
        return self.duplicate(range_bounds=range_bounds, azi_bounds=azi_bounds, with_detector=False)

    
    def _reduce_array(self, data, reduce_by=None, method="downscale_local_mean", keep_FirstDim=False):
        """
        Function to reduce the size of the data arrays by factor of reduce_by.
        The reduction bins the pixel data into reduce_by x reduce_by bins. 

        Parameters
        ----------
        data : np.array
            Data to be reduced.
        reduce_by : int, optional
            DESCRIPTION. The default is None.
        sum : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        data : np.array
            DESCRIPTION.

        """
        
        if reduce_by is False or (reduce_by is None and self.reduce_by is None):
            # reduce_by = False is used by fill_data to make sure this function is passed
            # if both are none then there is nothing to do.
            return data
        
        elif reduce_by is not None:
            reduction = int(reduce_by)
        else:
            reduction = int(self.reduce_by)
        
        reduction = np.int_(np.ones([1,np.ndim(data)]) * reduction)[0]
        if keep_FirstDim == True:
            reduction[0] = 1
        reduction = tuple(reduction)
        # if np.ndim(data) == 3:
        #     reduction = (1, reduction, reduction)
        # else:
        #     reduction = (reduction, reduction)
                
        msk=None
        if ma.is_masked(data) == True:
            msk = data.mask
            if method == "rescale":
                raise NotImplementedError(f"'{method}' has not been implemented as a reuction method")
                msk = rescale(msk, 1/reduction, anti_aliasing=False)
            elif method == "resize":
                raise NotImplementedError(f"'{method}' has not been implemented as a reuction method")
                msk = resize(
                    msk, (msk.shape[0] // reduction, msk.shape[1] // reduction), anti_aliasing=True
                )
            elif method=="downscale_local_mean":
                msk = downscale_local_mean(msk, reduction)
            else:
                raise ValueError(f"Unrecognised data reduction method: {method}")
        else:
            msk=False
        
        if method == "rescale":
            raise NotImplementedError(f"'{method}' has not been implemented as a reuction method")
            data_out = ma.array(rescale(data, 1/reduction, anti_aliasing=False), mask=msk)
        elif method == "resize":
            raise NotImplementedError(f"'{method}' has not been implemented as a reuction method")
            data_out = ma.array(resize(
                data, (data.shape[0] // reduction, data.shape[1] // reduction), anti_aliasing=True
            ), mask=msk)
        elif method=="downscale_local_mean":
            data_out = ma.array(downscale_local_mean(data, reduction), mask=msk)
        else:
            raise ValueError(f"Unrecognised data reduction method: {method}")
                
        return data_out



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
