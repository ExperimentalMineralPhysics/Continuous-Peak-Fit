"""
Loads all available image processing modules

The modules all require "Write" in their filename.
This allows the addition of new modules by just adding them to the directory.

The modules loaded from each Write*.py file must contain a method called WriteOutput. 

"""

import os as _os
from importlib import import_module as _import_module
from types import ModuleType as _ModuleType

module_list: list[str] = []
_new_module: dict[str, _ModuleType] = {}
for _module_path in _os.listdir(_os.path.dirname(__file__)):
    if (
        _module_path == "__init__.py"
        or _module_path[-3:] != ".py"
        or _module_path[:2] == "._"
        or not "Write" in _module_path
    ):
        # do not list the file to be loaded
        pass
    else:
        _output_module = _module_path[:-3]  # Remove ".py"
        module_list.append(_output_module)
        _module: _ModuleType = _import_module(f"cpf.spot_outputs.{_output_module}")
        _new_module[_output_module] = _module