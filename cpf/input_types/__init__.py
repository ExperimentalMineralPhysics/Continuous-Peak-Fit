# import os
# for module in os.listdir(os.path.dirname(__file__)):
#     if module == '__init__.py' or module[-3:] != '.py':
#         continue
#     __import__(module[:-3], locals(), globals())
# del module

# for name in os.listdir("plugins"):
#     if name.endswith(".py"):
#           #strip the extension
#          module = name[:-3]
#          # set the module name in the current global name space:
#          globals()[module] = __import__(os.path.join("plugins", name)


# for module in os.listdir(os.path.dirname(__file__)):
#     if module == '__init__.py' or module[-3:] != '.py':
#         continue
#     __import__(module[:-3], locals(), globals())
# del module


# for name in os.listdir(os.path.dirname(__file__)):
#     if name.endswith(".py") and 'Functions' in name:
#           #strip the extension
#          module = name[:-3]
#          # set the module name in the current global name space:
#          # globals()[module] = __import__(os.path.join(os.path.dirname(__file__), name))
#          globals()[module] =  __import__(os.path.join(os.path.dirname(__file__),name))

# import os
# print(os.listdir(os.path.dirname(__file__)))
# for module in os.listdir(os.path.dirname(__file__)):
#     print(module)
#     print((module == '__init__.py' or module[-3:] != '.py') and not 'Functions' in module)
#     if module == '__init__.py' or module[-3:] != '.py' or not 'Functions' in module:
#         print("skip")
#     else:
#         print("load")
#         mdl = __import__(module[:-3], fromlist=[None])
# del module



# from os.path import dirname, basename, isfile
# import glob
# modules = glob.glob(dirname(__file__)+"/*.py")
# __all__ = [ basename(f)[:-3] for f in modules if isfile(f) and 'Functions' in basename(f)] # exclude __init__.py


# new_module = {}
# for output_module in __all__:
#     module = __import__(__all__, fromlist=[None])
#     new_module[output_module] = module
# return new_module

# for module in range(len(__all__)):
#     __import__(__all__[module], locals(), globals())
# del module

# import os
# for name in os.listdir("plugins"):
#     if name.endswith(".py"):
#           #strip the extension
#          module = name[:-3]
#          # set the module name in the current global name space:
#          globals()[module] = __import__(os.path.join("plugins", name))





from . import DioptasFunctions
from . import GSASIIFunctions
from . import MedFunctions




# import os
# print(os.listdir(os.path.dirname(__file__)))


# new_module = {}
# for module in os.listdir(os.path.dirname(__file__)):
#     if module.endswith(".py") and 'Functions' in module:
#         md = __import__(os.path.join(os.path.dirname(__file__),module), fromlist=[None])
#         new_module[module] = md
