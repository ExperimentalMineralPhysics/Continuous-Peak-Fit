#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import re
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from copy import deepcopy
from skimage.util import shape

from cpf.util.logging import get_logger

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

    
    def _reduce_array(self, data, reduce_by=None, keep_FirstDim=False, polar=False):
        """
        Reduces the size of the data array by factor of reduce_by.
        The default returned values are the arithmetric mean of 
        'reduce_by' x 'reduce_by' blocks of data.
        
        Parameters
        ----------
        data : np.array
            Data array to be reduced.
        reduce_by : int, optional
            Factor to reduce the data by. The default is None.
        keep_FirstDim : bool, optional
            Do not reduce the first dimension of the data array if True.
            Used for compund detectors (e.g. MED and EXRFlvp type). The default is False.
        polar : bool, optional
            Account for wrapping of anglular data, where +/-180 degrees have the 
            same value and a straight mean of the values is incorrect.
            The default is False.

        Returns
        -------
        data_out : np.array
            Data array reduced by 'reduce_by'.
        """
        
        def do_reduction(data, reduction, mask=False, polar=False, rot=90, threshold=180, azm_min_max=[-180, 180]):
            """
            Do the resizing of the data array.
            Downscale using the local mean for each block
            
            :param data: array for be reduced
            :type data: ma.array
            :param reduction: Integer list of how much to reduce the data by
            :type reduction: list
            :param msk: DESCRIPTION, defaults to None
            :type msk: TYPE, optional
            :param method: methd to use for reduction, defaults to "downscale_local_mean"
            :type method: string, optional
            :raises NotImplementedError: Reduction method is not implemented
            :raises ValueError: Raised for unrecognised reduction method
            :return: reduced data array
            :rtype: ma.arrray

            """
            # Separated out from _reduce_array to prevent repeating of the same lines of code. 
            
            # Started from:
            # data_out = ma.array(downscale_local_mean(data, reduction), mask=mask)
            # but this doesn't work due to needing to check the angles for wrapping. 
            # downscale_local_mean calls other functions in skimage. call these as needed
            # shape.view_as_blocks() called by skimage.measure.block.block_reduce()
            
            # FIXME: something strange is going on here. if i dont round the 
            # FIXME: values then end up with np.int_(1) = 0 for BCC1 input and reduce_by = 3
            # FIXME: see following print statements for minimum working example. 
            if 0:
                print("data.shape:", data.shape)
                print("reduction:", reduction)
                print("run:   np.ceil(np.array(data.shape) / np.array(reduction)) , (np.array(data.shape) / np.array(reduction))")
                print(np.ceil(np.array(data.shape) / np.array(reduction)) , (np.array(data.shape) / np.array(reduction)))
                print("run:   np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction))")
                print(np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))
                print("run:   (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction))) * reduction")
                print( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction))) * reduction)
                print("run:   np.round(  (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction))) * reduction) " )
                print(  np.round(  (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction))) * reduction)  )
                print("First value of pevious array is:", np.round(  (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction))) * reduction)[0])
                print("THIS IS THE LINE I DO NOT UNDERSTAND")
                print("run:   1 - np.int32( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction )[0] " )
                print(1 - np.int32( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction )[0]  )
                print("run:   1 - ( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction )[0] " )
                print(1 - ( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction )[0]  )
                print("run:  np.int32( 1 - ( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction )[0] " )
                print(np.int32(1 - ( (np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction )[0]  ))
            
            #needs padding before can view as blocks. 
            #this round function is reuqired to deal with whatever mess is above in the numbers.
            pad = np.round((np.ceil(np.array(data.shape) / np.array(reduction)) - (np.array(data.shape) / np.array(reduction)))*reduction)
            # format pad as required by np.pad                
            pad = [(0, np.int_(pad[i])) for i in range(len(pad))]
            
            # pad and block the data
            data = np.pad(data, pad_width=pad, constant_values=np.nan)
            blocked = shape.view_as_blocks(data, reduction)

            #reduce the data. 
            data_out = np.nanmean(blocked, axis=tuple(range(data.ndim, blocked.ndim)))
                
            if polar == True:
                # the data could wrap around and we need to deal with that. 
                # the data wrapping is reason that cannot just call skimage.measure.block.block_reduce
                
                # see if the data actually wraps round 
                range1   = np.nanmax(blocked, axis=tuple(range(data.ndim, blocked.ndim))) - np.nanmin(blocked, axis=tuple(range(data.ndim, blocked.ndim)))
                wrapped = np.where(range1 >= threshold)
                
                if wrapped:
                    # some of the blocks have range greater than threshold
                    
                    # rotate the data so can get good values for average
                    # but need to make sure it all stays within expected ranges. 
                    rot = rot%360
                    rotated = (deepcopy(data)+rot)
                    rotated[rotated < azm_min_max[0]] = rotated[rotated < azm_min_max[0]] + 360
                    rotated[rotated > azm_min_max[1]] = rotated[rotated > azm_min_max[1]] - 360
                    # block and reduce the rotated data. 
                    blockedrot = shape.view_as_blocks(rotated, reduction)
                    data_outrot = np.nanmean(blockedrot, axis=tuple(range(data.ndim, blockedrot.ndim)))
                    
                    # if len(wrapped[0])!=0:
                    #     data_tmp = deepcopy(data_out)
                    #     data_wrapped_unclaned = data_tmp[wrapped]
                    
                    # replace the values that wrap around with good ones.                     
                    data_out[wrapped] = data_outrot[wrapped] - rot
                    
                    # if len(wrapped[0])!=0:
                    #     data_wrapped_cleaned = data_out[wrapped]
                    #     fig1, ax1 = plt.subplots(1,1)
                    #     ax1.scatter(data_wrapped_unclaned, data_wrapped_cleaned, 1, label="after")
                    #     ax1.legend()
                    #     ax1.set_xlabel("piror")
                    #     ax1.set_ylabel("after")

            return ma.MaskedArray(data_out, mask=mask)
        
        
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
                
        if ma.is_masked(data) == True:
            msk = do_reduction(data.mask, reduction, mask=False)
        else:
            msk=False
        
        data_out = do_reduction(data, reduction, mask=msk, polar=polar, azm_min_max=[self.azm_start, self.azm_end])
        
        if polar == True:
            data_out[data_out < self.azm_start] = data_out[data_out < self.azm_start] + 360
            data_out[data_out >= self.azm_end] = data_out[data_out >= self.azm_end] - 360
        
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
