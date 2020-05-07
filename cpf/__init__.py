__all__ = ['XRD_FitPattern', 'XRD_FitSubpattern', 'DioptasFunctions', 'GSASIIFunctions', 'MedFunctions',
           'PeakFunctions', 'Plot85Functions',  'WriteMultiFit', 'WritePolydefixED', 'XRD_FitPattern',
           'XRD_FitSubpattern', 'ImageMetaData', 'Mca', 'Med', 'generate_inputs']

from . import XRD_FitPattern
from . import XRD_FitSubpattern
from . import PeakFunctions
from . import DioptasFunctions
from . import GSASIIFunctions
from . import MedFunctions
from . import Plot85Functions
from . import WriteMultiFit
from . import WritePolydefixED
from . import Mca
from . import Med


def generate_inputs():

    try:
        fp = open('Example_MultiPeak_Inputs.py', 'r')
        fp.close()
        raise UserWarning('Example inputs already exist')
    except FileNotFoundError:
        pass

    try:
        fp = open('Example_MultiPeak_Inputs.py', 'w')

        fp.write("# Input parameters for BCC1_2GPa_100s_001.\n#\n\n\n")

        fp.write("# properties of the data files.\n")
        fp.write("datafile_directory = './Example1-Fe/'\n")
        fp.write("datafile_Basename  = 'BCC1_2GPa_10s_001_'\n")
        fp.write("datafile_Ending    = '.tif'\n")
        fp.write("datafile_StartNum  = 1\n")
        fp.write("datafile_EndNum    = 10\n")
        fp.write("datafile_NumDigit  = 5\n\n")

        fp.write("# Calibration and masking.\n")
        fp.write("Calib_type   = 'Dioptas'\n")
        fp.write("Calib_detector = 'Pilatus1M'\n")
        fp.write("Calib_data     = datafile_directory + 'CeO2_Pil207_E30_2Nov2016_001.tif'\n")
        fp.write("Calib_param    = datafile_directory + 'CeO2_cal_Dioptas.poni'\n")
        fp.write("Calib_mask     = datafile_directory + 'DiffractionMask_Dioptas.mask'\n")
        fp.write("Calib_pixels = 172\n\n")

        fp.write("# Fitting properties for peaks.\n")
        fp.write("fit_orders = [\n")
        fp.write("       {\n")
        fp.write("         'background': [0,0],\n")
        fp.write("         'peak': [{\n")
        fp.write("             'phase': 'Other',\n")
        fp.write("             'hkl': '000',\n")
        fp.write("             'd-space': 3,\n")
        fp.write("             'height': 1,\n")
        fp.write("             'profile': 0,\n")
        fp.write("             #'profile_fixed': 1,\n")
        fp.write("             'width': 0,\n")
        fp.write("             'symmetry': 2\n")
        fp.write("           },{\n")
        fp.write("             'phase': 'BCC-Fe',\n")
        fp.write("             'hkl': 110,\n")
        fp.write("             'd-space': 3,\n")
        fp.write("             'height': 16,\n")
        fp.write("             'profile': 0,\n")
        fp.write("             #'profile_fixed': 1,\n")
        fp.write("             'width': 0,\n")
        fp.write("             'symmetry': 2\n")
        fp.write("           }],\n")
        fp.write("         'range': [[10.8, 11.67]],\n")
        fp.write("         'PeakPositionSelection': [[1,-170.5,11.19],[1,-150.5,11.19],[1,-120.5,11.19],[1,-58.997,11.19],[1,59.289,11.19],[1,23.187,11.19],[1,23.212,11.19], [1,53.158,11.19], [1,73.246,11.19], [1,93.246,11.19], [1,123.246,11.19], [2,-170, 11.54], [2,-150, 11.54], [2,-110, 11.54], [2,-90, 11.54], [2,-45, 11.54], [2,-0, 11.54], [2,45, 11.54], [2,90, 11.54], [2,120, 11.54], [2,150, 11.54], [2,179, 11.54]]\n")
        fp.write("},\n")
        fp.write("       {\n")
        fp.write("         'background': [0,0],\n")
        fp.write("         'peak': [{\n")
        fp.write("             'phase': 'BCC-Fe',\n")
        fp.write("             'hkl': 200,\n")
        fp.write("             'd-space': 3,\n")
        fp.write("             'height': 8,\n")
        fp.write("             'profile': 0,\n")
        fp.write("             #'profile_fixed': 1,\n")
        fp.write("             'width': 0,\n")
        fp.write("             'symmetry': 2\n")
        fp.write("           }],\n")
        fp.write("         'range': [[16.1,16.6]]\n")
        fp.write("       },\n")
        fp.write("       {\n")
        fp.write("         'background': [0,0],\n")
        fp.write("         'peak': [{\n")
        fp.write("             'phase': 'BCC-Fe',\n")
        fp.write("             'hkl': 211,\n")
        fp.write("             'd-space': 3,\n")
        fp.write("             'height': 8,\n")
        fp.write("             'profile': 0,\n")
        fp.write("             #'profile_fixed': 1,\n")
        fp.write("             'width': 0,\n")
        fp.write("             'symmetry': 2\n")
        fp.write("           }],\n")
        fp.write("         'range': [[19.8, 20.27]]\n")
        fp.write("       },\n")
        fp.write("      {\n")
        fp.write("          'background': [1],\n")
        fp.write("          'peak': [{\n")
        fp.write("             'phase': 'BCC-Fe',\n")
        fp.write("             'hkl': 220,\n")
        fp.write("              'd-space': 3,\n")
        fp.write("              'height': 8,\n")
        fp.write("              'profile': 0,\n")
        fp.write("              'profile_fixed': 0.5,\n")
        fp.write("              'width': 0,\n")
        fp.write("              'symmetry': 2\n")
        fp.write("            }],\n")
        fp.write("#          'PeakPositionSelection': [[1,-120.5,23.226],[1,-58.997,23.186],[1,59.289,23.208],[1,23.187,120.893],[1,23.212,6.367], [1,23.158,90.025], [1,23.246,25.729]]\n")
        fp.write("          'range': [[23.1, 23.5]]\n")
        fp.write("        },\n")
        fp.write("       {\n")
        fp.write("         'background': [1,1],\n")
        fp.write("         'peak': [{\n")
        fp.write("             'phase': 'BCC-Fe',\n")
        fp.write("             'hkl': 310,\n")
        fp.write("             'd-space': 2,\n")
        fp.write("             'height': 1,\n")
        fp.write("             'profile': 0,\n")
        fp.write("             'profile_fixed': 0.5,\n")
        fp.write("             'width': 0,\n")
        fp.write("             'symmetry': 2\n")
        fp.write("           }],\n")
        fp.write("         'range': [[25.7, 26.4]]\n")
        fp.write("       },\n")
        fp.write("    ]\n\n")

        fp.write("#Number of bins for initial fitting.\n")
        fp.write("AziBins = 90\n\n")


        fp.write("#Output settings\n")
        fp.write("Output_directory   = './'\n")
        fp.write("Output_type        = 'MultiFit'\n")
        fp.write("Output_NumAziWrite = 90\n\n")

        fp.close()


    except PermissionError:
        print('Write permission required')