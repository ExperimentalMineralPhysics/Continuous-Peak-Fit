#!/usr/bin/env python
# -*- coding: utf-8 -*-

__all__ = ["W2010_10element"]


"""
This file contains a list of the default energy dispersive detectors that continuous peak fit recognises.

New detectors can either be added to this list or added via a csv file. For details see ...
"""
# from cpf.XRD_FitPattern import logger
from cpf.logger_functions import logger


def W2010_10element():
    """
    Azimuthal detector positions for the 10 element detector of Weidner et al. (2010).

    Donald J. Weidner, Michael T. Vaughan, Liping Wang, Hongbo Long, Li Li, Nathaniel A. Dixon, and William B. Durham
    Precise stress measurements with white synchrotron x rays
    Review of Scientific Instruments 81, 013903 (2010); https://doi.org/10.1063/1.3263760

    """

    detector_azimuths = [
        0.0,  # 1, top element
        22.5,  # 2
        45.0,  # 3
        67.5,  # 4
        90.0,  # 5
        112.5,  # 6
        135.0,  # 7
        157.5,  # 8
        180.0,  # 9, bottom element
        270.0,
    ]  # 10, element on otherside (often obscured)

    return detector_azimuths
