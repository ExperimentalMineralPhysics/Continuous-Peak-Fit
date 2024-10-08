#!/usr/bin/env python
# -*- coding: utf-8 -*-

# from . import DioptasFunctions
# from . import GSASIIFunctions
# from . import MedFunctions
import os
from importlib import import_module
from types import ModuleType

"""
Loads all available input format modules

The output modules all require "Functions" in their filename.
This allows the addition of new output types by just adding them to the directory.

"""
module_list: list[str] = []
new_module: dict[str, object] = {}
for module_path in os.listdir(os.path.dirname(__file__)):
    if (
        module_path == "__init__.py"
        or module_path[-3:] != ".py"
        or not "Functions" in module_path
    ):
        # do not list the file to be loaded
        pass
    else:
        output_module = module_path[:-3]  # Remove ".py"
        module_list.append(output_module)
        module: ModuleType = import_module(f"cpf.input_types.{output_module}")
        new_module[output_module] = module
