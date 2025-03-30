# -*- coding: utf-8 -*-
"""
About
=====
Provides (or will provide) a number of options for throesholding or limiting the
pixel intensity in images.

Individual pixels in a diffraction ring are frequently much brighter than the
rest of the peak. This bright intensity skews the least squares fitting and
so the brightest spots are ether masked (or to be developed) eroded (keeping
total intentsity in the data the same)

Current useage of intensity limits:
    Each range has "imax" and "imin" as optional arguments
    These optional arguments can take the format:

    "imin": %i
    "imin": "[0<%f<100]%"
    "imin": {"percentile": "[0<%f<100]%",
             "minbelow": 20}
    "imax": {"limit": "[0<%f<100]%",
             "abovebelow": 200}
    e.g.
    "imin": 234
    "imin": "99.9%"
    "imax": {"limit": 300,
         "abovebelow": 200}
    "imax": {"limit": "99%",
         "abovebelow": 200}


    Where a numberical value masks the data above/below that value and the percentage
    as a string masks above or below that percentile.

    The combined condition only applies the limit condition if the maximum or minimum is
    above or below the "abovebelow" value.

Simon Hunt, March 2023
"""

__version__ = "0.1"

import numpy as np

from cpf.logging import CPFLogger

logger = CPFLogger("cpf.BrightSpots")


def SpotProcess(sub_data, settings_for_fit):
    imax = np.inf
    imin = -np.inf
    if "imax" in settings_for_fit.subfit_orders:
        if isinstance(settings_for_fit.subfit_orders["imax"], str):
            imax = np.percentile(
                sub_data.intensity.flatten(),
                float(settings_for_fit.subfit_orders["imax"].strip("%")),
            )
        elif isinstance(settings_for_fit.subfit_orders["imax"], dict):
            if (
                np.max(sub_data.intensity.flatten())
                > settings_for_fit.subfit_orders["imax"]["abovebelow"]
            ):
                if isinstance(settings_for_fit.subfit_orders["imax"]["limit"], str):
                    imax = np.percentile(
                        sub_data.intensity.flatten(),
                        float(
                            settings_for_fit.subfit_orders["imax"]["limit"].strip("%")
                        ),
                    )
                else:
                    imax = int(settings_for_fit.subfit_orders["imax"]["limit"])
        else:
            imax = int(settings_for_fit.subfit_orders["imax"])

    if "imin" in settings_for_fit.subfit_orders:
        if isinstance(settings_for_fit.subfit_orders["imin"], str):
            imin = np.percentile(
                sub_data.intensity.flatten(),
                float(settings_for_fit.subfit_orders["imin"].strip("%")),
            )
        elif isinstance(settings_for_fit.subfit_orders["imin"], dict):
            if (
                np.max(sub_data.intensity.flatten())
                > settings_for_fit.subfit_orders["imin"]["abovebelow"]
            ):
                if isinstance(settings_for_fit.subfit_orders["imin"]["limit"], str):
                    imax = np.percentile(
                        sub_data.intensity.flatten(),
                        float(
                            settings_for_fit.subfit_orders["imin"]["limit"].strip("%")
                        ),
                    )
                else:
                    imax = int(settings_for_fit.subfit_orders["imin"]["limit"])
        else:
            imin = int(settings_for_fit.subfit_orders["imin"])

    sub_data.set_mask(intensity_bounds=[imin, imax])  # (i_min=imin, i_max=imax)

    return sub_data
