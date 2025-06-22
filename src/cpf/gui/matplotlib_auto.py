###
import matplotlib.pyplot as pl


class matplotlib_inline:
    def __init__(self):
        super(matplotlib_inline, self).__init__()
        try:
            get_ipython().magic("matplotlib inline")
        except:
            pl.ion()
