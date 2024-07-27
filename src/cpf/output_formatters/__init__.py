#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

"""
Loads all available output format modules

The output modules all require "Write" in their filename.
This allows the addition of new output types by just adding them to the directory.

Each output formatter must contain two modules called "Requirements" and "WriteOutput"
 
"""
module_list = []
for module in os.listdir(os.path.dirname(__file__)):
    if (
        module == "__init__.py"
        or module[-3:] != ".py"
        or module[:2] == "._"
        or not "Write" in module
    ):
        # do not list the file to be loaded
        pass
    else:
        module_list.append(module[:-3])

new_module = {}
for output_module in module_list:
    module = __import__("cpf.output_formatters." + output_module, fromlist=[None])
    new_module[output_module] = module

# del os
# del module
# del module_list
# del output_module
