# Continuous-Peak-Fit

## Note on the state of the code and documentation

The documentation is undercooked and much of it has not been written. If you want to run the code it might be simplest to email me (Simon Hunt) at simon.hunt@manchester.ac.uk. Alternatively, please raise issues in the repository.

The development of the code is ongoing. Each new data set we process finds new limits to the code and I am happy to fix such issues as they arise.

## Installation

Installation is available through Git or via pip. In order for some of the dependencies to be compiled and installed, Microsoft Visual C++ 14.0 is required. The simplest way to acquire it is to download and install "Microsoft Visual C++ Build Tools" using the following link: https://visualstudio.microsoft.com/visual-cpp-build-tools/

### Using Git

You can clone the repository via Git using the following command:

`git clone https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit.git`

After cloning the repository, while having it as your current working directory, you can install it via pip as follows:

`pip install .`

### Using pip

Continuous-Peak-Fit can also be installed directly from pip using the following command:

`pip install continuous-peak-fit`

This should automatically install the required dependencies.

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
