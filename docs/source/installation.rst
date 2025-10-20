.. Continuous-Peak-Fit documentation master file, created by
   sphinx-quickstart on Fri Aug 28 15:27:46 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

======================
Download and Install
======================

.. _cpf:    https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit
.. _numpy:  http://numpy.org/
.. _Fabio:  https://github.com/silx-kit/fabio
.. _lmfit:  https://lmfit.github.io/lmfit-py/index.html

Prerequisits
------------

Continuous peak fit works with Python 3.7. 

Continuous peak fit requires the following Python packages:
   * `NumPy`_ version 1.10 or higher.
   * scipy version 1.6.1 or higher.
   * matplotlib version 3.3.4 or higher.
   * `lmfit`_ version 1.0.2 or higher.
   * pyFAI version 0.20.0 or higher.
   * `Fabio`_ version 0.11.0 or higher.
   * setuptools
   * pycairo 
   * h5py
   * pillow -- no longer required [I THINK]
   * pyopencl

These are avaliable on PyPI, and are installed automatically if installing with pypi.
Note that many of these packages are interdependent so not all of them will have to be installed as separate events. 


Download
---------

The latest release version of continuous peak fit is |release| and is avaliable from [ADD LINK].

However, it has been significantly updated since |release| was made.  The latest version of *continuous-peak-fit* is avaliable from GitHub at:

`https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit <https://github.com/ExperimentalMineralPhysics/Continuous-Peak-Fit>`_


Installation
------------

It is easiest to install cpf with:

 .. code-block:: python    

   pip install cpf


However, if `cpf` has been downloaded directly from GitHub, add the cpf directory to the python path or make the active directory and import the package:

 .. code-block:: python    

   Import cpf

For the import to complete, additional packages are needed are (see above). The import will fail if all the packages are not present and it will tell you which ones are missing. 



Installation with virtualenv
----------------------------
In your current environment run `pip install virtualenv`.

Pick/create a directory in which you will store your new Python environment. I'd typically recommend something memorable and simply like `~/envs`.

Run `virtualenv ~/envs/cpf` (or whatever folder name you want your new Python environment to be known by). If :code:`virtualenv` by itself fails, you might need to run `python -m virtualenv ~/envs/cpf` instead.

Activate it by running `source ~/envs/cpf/bin/activate`

Now change dir to where your CPF package is saved (cd /path/to/cpf/repo)

Install an editable version of the environment using `pip install -e .`. This will ensure the imports in your package are properly mapped internally, and that any edits to the repo get reflected when you run the environment.




Installation with conda
-----------------------

[This assumes that there is a python envoriment installed on the machine. If there is not then ...]

Step 1: Make a new python enviroment:

 .. code-block:: python    

   conda create --name *cpf*

where *cpf* is the name of the new enviroment. 


Step 2: Activate the enviroment

 .. code-block:: python    

   conda activate *cpf*


Step 3: Install cpf.

Either in a local directory: 

source
or install in enviroment's default location:

 .. code-block:: python    

   conda install 













