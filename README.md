# Continuous-Peak-Fit

## Installation

Installation is available through git or via pip

From git use

`git clone https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit.git`

To install the required dependencies use, in the Continuous-Peak-Fit directory use 

`pip install requirements.txt`

# 

For pip installation use

`pip install continuous-peak-fit`

this should automatically install the required dependencies.


## Usage

Continuous-Peak-Fit can be run in a number of ways but requires an inputs file containing information about where your 
data files are stored as well as information of the number of peak and appropriate ranges for the fit. Example data and 
inputs file are available in the Example1-Fe directory from the git repository. 

If using a pip install an example inputs file can be generated using the 'CPF_generate_inputs' executable or from within
 python using

`import cpf`

`cpf.generate_inputs()`

Fitting can then be performed using 

`CPF_XRD_FitPattern <inputs_file>` from the command line

or from within python

`import cpf`

`cpf.XRD_FitPattern.execute(settings_file=<inputs_file>)`

