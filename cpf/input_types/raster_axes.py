import numpy as np
import numpy.ma as ma

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import time

"""
edited after https://gist.github.com/astrofrog/1017095
downloaded: May 2022



if might be bettter to use 
matplot lib : pcolormesh (https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolormesh.html#differences-pcolor-pcolormesh) 
    Which makes plots with none rectangular grids. 
    e.g. https://stackoverflow.com/questions/32462548/plot-deformed-2d-mesh-with-python

The plot can be made using the corners of the pixels as the points.
The pixel corners can be got from pyFAI.Geometry.chi_corner or tth_corner etc...
    https://pyfai.readthedocs.io/en/master/api/pyFAI.html#module-pyFAI.geometry

"""

# def make_colormap(color):

#     r, g, b = colors.colorConverter.to_rgb(color)

#     cdict = {'red': [(0.0, r, r),
#                     (1.0, 1.0, 1.0)],

#              'green': [(0.0, g, g),
#                        (1.0, 1.0, 1.0)],

#              'blue':  [(0.0, b, b),
#                        (1.0, 1.0, 1.0)]}

#     return colors.LinearSegmentedColormap('custom', cdict)


class RasterAxes(plt.Axes):

    def __init__(self, *args, **kwargs):
        plt.Axes.__init__(self, *args, **kwargs)
        self.callbacks.connect('ylim_changed', self._update_all_scatter)
        self._scatter_objects = []
        self._scatter_images = {}

    def scatter(self, x, y, c='red'):
        scatter = (x, y, c)
        self._scatter_objects.append(scatter)
        self._update_scatter(scatter)

    def _update_all_scatter(self, event):
        for scatter in self._scatter_objects:
            self._update_scatter(scatter)

    def _update_scatter(self, scatter):

        x, y, c = scatter

        dpi = self.figure.get_dpi()

        autoscale = self.get_autoscale_on()

        if autoscale:
            self.set_autoscale_on(False)

        width = self.get_position().width * self.figure.get_figwidth()
        height = self.get_position().height * self.figure.get_figheight()

        nx = int(round(width * dpi)/7)
        ny = int(round(height * dpi)/7)

        nx, ny = GetRasterImageSizes(self)
        array_sum, array_num, extent = RasterImage(scatter, nx=nx, ny=ny)
        
        # xmin, xmax = self.get_xlim()
        # ymin, ymax = self.get_ylim()

        # ix = ((x - xmin) / (xmax - xmin) * float(nx)).astype(int)
        # iy = ((y - ymin) / (ymax - ymin) * float(ny)).astype(int)

        # array_num = ma.zeros((nx, ny), dtype=int)
        # array_sum = ma.zeros((nx, ny), dtype=float)

        # keep = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)

        # if len(x) == len(c):
        #     array_num[ix[keep], iy[keep]] +=1
        #     array_sum[ix[keep], iy[keep]] +=c[keep]
        # else:
        #     array_num[ix[keep], iy[keep]] +=1
        #     array_sum[ix[keep], iy[keep]] +=1   
        
        # array_num = array_num.transpose()
        # array_sum = array_sum.transpose()

        # array_num.mask = array_num == 0
        # array_sum.mask = array_sum == 0

        if id(x) in self._scatter_images:
            self._scatter_images[id(x)].set_data(array_num)
            # self._scatter_images[id(x)].set_extent([xmin, xmax, ymin, ymax])
        else:
            # self._scatter_images[id(x)] = self.imshow(array_sum, extent=[xmin, xmax, ymin, ymax], aspect='auto', interpolation='nearest')
            self._scatter_images[id(x)] = self.imshow(array_sum, aspect='auto', interpolation='nearest', extent=extent)

        if autoscale:
            self.set_autoscale_on(True)



def RasterImage(scatter, nx=200,ny=200):
    
    x, y, c = scatter

    xmin = np.min(x)
    xmax = np.max(x)
    ymin = np.min(y)
    ymax = np.max(y)
    # ymin, ymax = self.get_ylim()

    ix = ((x - xmin) / (xmax - xmin) * float(nx)).astype(int)
    iy = ((y - ymin) / (ymax - ymin) * float(ny)).astype(int)

    array_num = ma.zeros((nx, ny), dtype=int)
    array_sum = ma.zeros((nx, ny), dtype=float)

    keep = (ix >= 0) & (ix < nx) & (iy >= 0) & (iy < ny)

    if len(x) == len(c):
        array_num[ix[keep], iy[keep]] +=1
        array_sum[ix[keep], iy[keep]] +=c[keep]
    else:
        array_num[ix[keep], iy[keep]] +=1
        array_sum[ix[keep], iy[keep]] +=1   
    
    array_num = array_num.transpose()
    array_sum = array_sum.transpose()

    array_num.mask = array_num == 0
    array_sum.mask = array_sum == 0
    
    # print( array_sum)
    # print( array_num)
    
    return array_sum, array_num, [xmin, xmax, ymin, ymax]
    

def GetRasterImageSizes(ax, min_pix=5):
    
    dpi = ax.figure.get_dpi()

    autoscale = ax.get_autoscale_on()

    if autoscale:
        ax.set_autoscale_on(False)

    width = ax.get_position().width * ax.figure.get_figwidth()
    height = ax.get_position().height * ax.figure.get_figheight()

    nx = int(round(width * dpi)/min_pix)
    ny = int(round(height * dpi)/min_pix)
    
    if autoscale:
        ax.set_autoscale_on(True)
    
    return nx, ny


if __name__ == "__main__":

    
    n = 10000
    x = np.random.random(n)*6
    y = np.random.random(n)*2
    z = np.random.random(n)
    z = z**2 * 100
    
    

    # n = 100000
    # x = np.random.random(n)
    # y = np.random.random(n)
    t1 = time.time()
    fig = plt.figure()
    ax = RasterAxes(fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    ax.scatter(x, y, c=z)
    plt.show()
    t2 = time.time()
    print("Rastered image time: "+str(t2-t1))
    
    # # n = 100000
    # # x = np.random.random(n)
    # # y = np.random.random(n)
    t1 = time.time()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    fig.add_axes(ax)
    sc = (x,y,z)
    nx,ny = GetRasterImageSizes(ax, min_pix=5)
    im =RasterImage(sc, nx=nx,ny=ny)
    ax.imshow(im[0], extent=np.round(im[2]),aspect="auto")
    plt.show()
    t2 = time.time()
    print("Raster scatter time: "+str(t2-t1))
    
    t1 = time.time()
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)#fig, [0.1, 0.1, 0.8, 0.8])
    fig.add_axes(ax)
    ax.scatter(x, y, c=z, s=0.5)
    plt.show()
    t2 = time.time()
    print("Normal scatter time: "+str(t2-t1))




    # fig.canvas.draw()