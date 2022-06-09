#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from . import DioptasFunctions
# from . import GSASIIFunctions
# from . import MedFunctions
import os

"""
Loads all available input format modules

The output modules all require "Functions" in their filename.
This allows the addition of new output types by just adding them to the directory.
 
"""
module_list = []
for module in os.listdir(os.path.dirname(__file__)):
    if module == '__init__.py' or module[-3:] != '.py' or not 'Functions' in module:
        # do not list the file to be loaded
        pass
    else:
        module_list.append(module[:-3])
        
new_module = {}
for output_module in module_list:
    module = __import__("cpf.input_types." + output_module, fromlist=[None])
    new_module[output_module] = module

# del os
# del module
# del module_list
# del output_module