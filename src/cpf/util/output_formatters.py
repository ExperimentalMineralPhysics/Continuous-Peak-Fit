import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg


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
