# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 21:40:54 2022

@author: g05296ar
"""

import matplotlib.pyplot as pl


class matplotlib_qt:
    def __init__(self):
        super(matplotlib_qt, self).__init__()
        try:
            get_ipython().magic("matplotlib qt")
        except:
            pl.ion()
