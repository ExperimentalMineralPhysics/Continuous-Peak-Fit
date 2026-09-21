import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure


def mplfig_to_npimage(fig):
    """
    Converts a matplotlib figure to a RGB frame after updating the canvas.

    This is our own implementation of a function that was removed from MoviePy >=2,
    which we only make use of once in WriteOutput. This should fix the incompatibility
    with MoviePy >=2.

    Credits: https://github.com/Zulko/moviepy/issues/2297#issuecomment-2609419770
    """

    canvas = FigureCanvasAgg(fig)
    canvas.draw()  # update/draw the elements

    # get the width and the height to resize the matrix
    l, b, w, h = canvas.figure.bbox.bounds
    w, h = int(w), int(h)

    #  exports the canvas to a memory view and then to a numpy nd.array
    mem_view = canvas.buffer_rgba()  # Update to Matplotlib 3.8
    image = np.asarray(mem_view)
    return image[:, :, :3]  # Return only RGB, not alpha.


def figure_suptitle_space(figure: Figure, topmargin: float = 1):
    """increase figure size to make topmargin (in inches) space for
    titles, without changing the axes sizes.
    after: https://stackoverflow.com/questions/55767312/how-to-position-suptitle#55768955

            Acutally now does this by compresssing the axes away from the top of the figure.
    """

    axes = figure.axes
    pos = []
    for i in range(len(axes)):
        pos.append(axes[i].get_position().bounds)
    w, h = figure.get_size_inches()
    figh = h - topmargin  # - (1-s.y1)*h
    for i in range(len(axes)):
        al = pos[i][0]
        ab = pos[i][1] / h * figh
        aw = pos[i][2]
        ah = pos[i][3] / h * figh
        axes[i].set_position((al, ab, aw, ah))
