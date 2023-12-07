# Import modules
import os

from pathlib import Path

# Import functions to test
from cpf.Cascade_Test import make_im_from_flts

# Set root working directory
# Root folder
root = Path(os.path.dirname(__file__))

# Data directory
data = "Example1-Fe"
data = root.joinpath(data)

# File name
file = r"BCC1_MultiPeak_input_Dioptas"

# Test function
os.chdir(data)
make_im_from_flts(file, debug=False)