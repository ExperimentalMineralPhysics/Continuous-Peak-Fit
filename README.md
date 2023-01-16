# Continuous-Peak-Fit

## Note on the state of the code and documentation
The documentation is undercooked and much of it has not been written. If you want to run the code it might be simplest to email me (Simon Hunt) at FirstName.LastName@manchester.ac.uk (replace as approprite). Or raise issues in the repository.

The development of the code is ongoing. Each new data set we process finds new limits to the code and I am happy to fix such issues as they arise.


## Installation
The following instructions are given for Python 3 on Ubuntu. The equivalent of the `apt` package manager on MacOS is to use Homebrew for package management via the `brew` command.

### Install Python and Pip
Ensure python and pip are installed on the system:

`sudo apt install python python3-pip python3-dev`

Clone the source code with

`git clone https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit.git`

### Install pre-requisites
Install packages that don't install properly with pip 'sudo apt install pkg-config libcairo2-dev ffmpeg libsm6 libxext6 python3-pyqt6 pyqt6-dev-tools`

### Configure environment using Pip
Install dependencies using `python3 -m pip install -r requirements.txt`.

## Docker
There is a docker image that can be built from the root `Continuous-Peak-Fit` directory using `docker build -t <name>:<tag> . -f docker/dockerfile`. You can then run the container using `docker run -it <name>:<tag> sh`.

## Usage

### GUI version

From the root Continuous-Peak-Fit directory run the GUI by running the entry point wrapper as `python3 main.py` or by running the `run.bat` script on Windows or the `run.sh` script on Linux/MacOS.

### Non-GUI version

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
