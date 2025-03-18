.. Continuous-Peak-Fit documentation master file


=====================================
Sytax for iterating over hdf5 files. 
=====================================

.. _cpf GitHub repository:   https://github.com/ExperimentalMineralPhysics/continuous-peak-fit


h5 files contain many images in a structure defined by keys. 

This iteration is designed to be as generic as possible while reflecting the structure the h5 files typical from synchrotron data sets.

We assume that the data set if organised in keys that increment numberically and that the ket directs to an array (of either a single or multiple) diffraction patterns).

There are two ways that the code accepts keys: a new way and the original way. 

original way.
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
Nor can it accept numberical values for the keys, rather than positions in a list of keys. 



New way
=====================================



 .. code-block:: python
  h5_data_key = "/{}/something/next_level/"
  h5_label_key = ["/{}/instrument/positioners/dz1","/{}/instrument/positioners/dy"]
  h5_iterate_key = [{start: 0, stop: -1, step: 1, index:"position"/"value" [default]}]
  h5_iterate_key = [{list: [5,6,9,10,24,15], index:"position"/"value" [default]}]
  h5_process = "iterate"



 .. code-block:: python
  h5_data_key = "/{}.1/something/next_level/"
  h5_label_key = ["/{}.2/instrument/positioners/dz1","/{}.2/instrument/positioners/dy"]
  h5_iterate_key = [{start: 0, stop: -1, step: 1, index:"position"/"value" [default], add: 0.1}]
  h5_iterate_key = [{list: [5,6,9,10,24,15], index:"position"/"value" [default]}, {process:"sum"}]


**h5_data_key**: is a f-string with "{}" to denote the position into which number(s) will be inserted to iterate over. 

**h5_label_key**: is an f-string or a list of f-strings, with "{}" to denote the position into which is being iterated over. The number of "{}" must be the same as that in h5_data_key. 

**h5_iterate_key**: list of dictionaries. The length of the list has to be the same as, or 1 greater than, the number of "{}" in 'h5_data_key'. 
If the list is 1 greater than the number of "{}" in h5_data_key, the last list applies to the data set pointed to by h5_data_key.

Possible values with in the list are: 

    * data set indicies. 

        * start / stop / step

        * list of values

    * index. This is a switch that determines if the indicies are either:

        * the "position" with in the key list or 

        * the numeric "value" of the h5 key to be used. 

    * "add". Value to add to the indicies. For example if all the indicies are n.1, then add = 0.1 allows use of start/stop/step. 

    * "process". Either "iterate" or "sum". Process to apply to the data array. If we are summing the data series then the label keys are averaged. This is used for example at Diamond where multiple exposures make up a single diffraction image. The detault is iterate. 

Every argment in the dictionary can be optional. In extremis, when all the entries in the h5 key are to be processed, it is licit to have an empty list. This is equivalent to {start: 0, stop: -1, step: 1, add: 0, process:"iterate"}



for example to process every other image... step: 2



If the list is present then the list values are substituted into the key. If start/stop are used then the start^th to stop^th keys are iterated over. with step. if start = 0 then we start from the very begining and if stop = -1 we go to the very end. 
If step is <0 then the processing starts at the end and works backwards. 

For example. see above. 

h5_operation: a string which is either "iterate" or "sum". This instructs the code either to apply the fitting to each image separately or to sum the images and then apply the process.



