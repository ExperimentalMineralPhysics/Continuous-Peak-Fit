# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 21:59:16 2022

@author: g05296ar
"""

import matplotlib.pyplot as pl


class matplotlib_inline():
    def __init__(self):
        super(matplotlib_inline, self).__init__()
        try:
            get_ipython().magic("matplotlib inline")
        except:
            pl.ion()
        import numpy as np
