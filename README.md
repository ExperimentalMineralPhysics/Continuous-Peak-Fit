# Continuous-Peak-Fit

## Note on the state of the code and documentation
The documentation is undercooked and much of it has not been written. If you want to run the code it might be simplest to email me (Simon Hunt) at FirstName.LastName@manchester.ac.uk (replace as approprite). Or raise issues in the repository.

The development of the code is ongoing. Each new data set we process finds new limits to the code and I am happy to fix such issues as they arise.


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

## Install using Pip

Install dependencies using `python3 -m pip install -r requirements.txt`.

You may need to first install `pkg-config` and `cairo` packages globally:

`sudo apt install libcairo2-dev pkg-config python3-dev`

## Usage

### GUI version

From the root Continuous-Peak-Fit directory run the GUI by running the entry point wrapper `python3 main.py`.

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

