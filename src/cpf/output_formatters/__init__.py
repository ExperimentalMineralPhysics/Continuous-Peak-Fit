#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from importlib import import_module
from types import ModuleType

"""
Loads all available output format modules

The output modules all require "Write" in their filename.
This allows the addition of new output types by just adding them to the directory.

Each output formatter must contain two modules called "Requirements" and "WriteOutput"

"""
module_list: list[str] = []
new_module: dict[str, ModuleType] = {}
for module_path in os.listdir(os.path.dirname(__file__)):
    if (
        module_path == "__init__.py"
        or module_path[-3:] != ".py"
        or module_path[:2] == "._"
        or "Write" not in module_path
    ):
        # do not list the file to be loaded
        pass
    else:
        output_module = module_path[:-3]  # Remove ".py"
        module_list.append(output_module)
        module: ModuleType = import_module(f"cpf.output_formatters.{output_module}")
        new_module[output_module] = module
