.. Continuous-Peak-Fit documentation master file


=====================================
Sytax for iterating over hdf5 files. 
=====================================

.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit


h5 files contain many images in a structure defined by keys. 

This iteration is designed to be as generic as possible while reflecting the structure the h5 files typical from synchrotron data sets.

We assume that the data set if organised in keys that increment numberically and that the key directs to an array (of either a single or multiple) diffraction patterns).



hdf5 key syntax
=====================================
The hdf5 file's data keys are pointed to via 'h5_datakey', which is a string. 
Asterixs in the key allow for iterations over key names, but the key must point to a dataset in the hdf5 file. 

How to iterate over the key names and the dataset are controlled by 'h5_iterate'. 


 .. code-block:: python
  h5_data    = '/*.1/measurement/p3'
  h5_iterate = [{start: 0, stop: 42, step: -1, index:"value", "label":"*"}
                    {"do":"iterate", 
                     "from": 0, 
                     "to": -1, 
                     "step": 1, 
                     "using":"position",
                     "label":['*','/*.1/instrument/positioners/dz1', '/*.1/instrument/positioners/dy']}]


**h5_data**: is a string with "*" to denote the position into which number(s) will be inserted to iterate over. 

**h5_iterate**: list of dictionaries. The length of the list has to be the same as, or 1 greater than, the number of "*" in 'h5_data_key'. 
The last entry in the list applies to the data set pointed to by h5_data string. 
Earlier entires in the list apply to the '*' in the h5_data key.  
The maximum length of h5_iterate is therefore 1 greater than the number of "*" in h5_data_key.

Possible values within each dictionary are: 

    * "do": Either "iterate" or "sum". Process to apply to the data array. If we are summing the data series then the label keys are averaged. This is used for example at Diamond Light Source where multiple exposures make up a single diffraction image. The detault is iterate. 

    * "from": where to start counting from. Does not have to be an integer. default is 0.

    * "to" -- where to stop. if set as -1 then goes to the end of the data set.

    * "step" - step in the list

    * NOT IMPLEMENTED: "list" this will allow a list of keys to be iterated over. 

    * "dim" -- which dimension of the dataset to iterate over. 

    * "using" This is a switch that determines if the indicies are either:

        * the "position" with in the key list or 

        * the numeric "value" of the h5 key to be used. 

    * "label": What to add as label for each data frame to give them unique file names.


Every argment in the dictionary can be optional. In extremis, when all the entries in the h5 key are to be processed, it is licit to have an empty list. This is equivalent to {"from": 0, "to": -1, step: 1, process:"iterate"}


for example to process every other image... step: 2



If the list is present then the list values are substituted into the key. If start/stop are used then the start^th to stop^th keys are iterated over. with step. if start = 0 then we start from the very begining and if stop = -1 we go to the very end. 
If step is <0 then the processing starts at the end and works backwards. 

For example. see above. 

h5_operation: a string which is either "iterate" or "sum". This instructs the code either to apply the fitting to each image separately or to sum the images and then apply the process.








original way (not recomended).
=====================================

 .. code-block:: python
  h5_key_list  = ['/', 'measurement/p3']
  h5_key_names = [["/",''], ["", 'instrument/positioners/dz1','instrument/positioners/dy']]
  #h5_key_names = ["",""]
  h5_key_start = [0, 0]
  h5_key_end   = [42,-1]
  h5_key_step  = [-1,1]
  h5_data      = 'iterate' # the other options would be "sum". "iterate" only applies to the bottom level.

but this cannot take lists for nonlineatrly increasing keys. 
Nor can it accept position values for the keys, rather than numerical values for the keys. 



