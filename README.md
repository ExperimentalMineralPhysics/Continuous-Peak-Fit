# Continuous-Peak-Fit

## Installation

### From Source
Source code can be obtained from this respository using `git`. `pip` can be used to build / install the dependencies.

Obtain the code using:

`git clone https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit.git`

To install the required dependencies, in the Continuous-Peak-Fit directory use:

`pip install requirements.txt` or `pip3 install -r requirements.txt` if using Python 3.

### From Pip

For pip installation use:

`pip install continuous-peak-fit` or `pip3 install continous-peak-fit` if using Python 3.

This should automatically install the required dependencies.

---

⚠️ **If this process fails to build the wheel for pycairo with "pkg-config not found" then take the following extra steps to resolve**

`brew install pkg-config` to install pkg-config using Homebrew
`brew install cairo` to isntall the cairo libraries which are linked to by pycairo.

---

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

