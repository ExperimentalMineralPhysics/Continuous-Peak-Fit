#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os as _os
from importlib import import_module as _import_module
from types import ModuleType as _ModuleType

"""
Loads all available output format modules

The output modules all require "Write" in their filename.
This allows the addition of new output types by just adding them to the directory.

Each output formatter must contain two modules called "Requirements" and "WriteOutput"

"""
module_list: list[str] = []
_new_module: dict[str, _ModuleType] = {}
for _module_path in _os.listdir(_os.path.dirname(__file__)):
    if (
        _module_path == "__init__.py"
        or _module_path[-3:] != ".py"
        or _module_path[:2] == "._"
        or "Write" not in _module_path
    ):
        # do not list the file to be loaded
        pass
    else:
        _output_module = _module_path[:-3]  # Remove ".py"
        module_list.append(_output_module)
        _module: _ModuleType = _import_module(f"cpf.output_formatters.{_output_module}")
        _new_module[_output_module] = _module