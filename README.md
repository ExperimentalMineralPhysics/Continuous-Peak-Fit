# Continuous-Peak-Fit

## Installation
The following instructions are given for Python 3 on Ubuntu. The equivalent of the `apt` command on MacOS is to use Homebrew for package management via the `brew` command.

Ensure python and pip are installed on the system:

`sudo apt install python` and `sudo apt install python3-pip`

You can check versions with:

```console
foo@bar:~$ python3 --version
Python 3.8.10
foo@bar:~$ pip3 --version
pip 20.0.2 from usr/lib/python3/dist-packages/pip (python 3.8)
```
Clone the source code with

`git clone https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit.git`

### Using Pipenv

You will need to first install `pkg-config` and `cairo` packages globally:

`sudo apt install libcairo2-dev pkg-config python3-dev`

A Python virtual environment can be used to install all required dependencies and the project has a Pipfile to facilitate this. First install `pipenv`:

`sudo apt install pipenv`

Create a virtual environment containing the dependencies listed in `Pipfile` by running the following inside the Continuous-Peak-Fit directory:

`pipenv install`

Then you can load the virtual environment with

`pipenv shell`

Typing `pipenv list` within the environment should confirm the dependencies have been installed correctly.

```console
Package           Version
----------------- ---------
asteval           0.9.28
contourpy         1.0.6
cycler            0.11.0
dill              0.3.6
fabio             0.14.0
fonttools         4.38.0
future            0.18.2
h5py              3.7.0
hdf5plugin        4.0.0
kiwisolver        1.4.4
lmfit             1.1.0
matplotlib        3.6.2
multiprocess      0.70.14
numexpr           2.8.4
numpy             1.23.5
packaging         21.3
pathos            0.3.0
Pillow            9.3.0
pip               20.0.2
pkg-resources     0.0.0
platformdirs      2.5.4
pox               0.3.2
ppft              1.7.6.6
pycairo           1.23.0
pyFAI             0.21.3
pyopencl          2022.3
pyparsing         3.0.9
python-dateutil   2.8.2
pytools           2022.1.13
scipy             1.9.3
setuptools        44.0.0
silx              1.1.0
six               1.16.0
typing-extensions 4.4.0
uncertainties     3.1.7
wheel             0.34.2
```

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

