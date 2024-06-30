#!/usr/bin/env python

import argparse
from pathlib import Path


def run():
    # Create a parser object to take additional arguments for this function
    parser = argparse.ArgumentParser(
        # This is the help message that will be displayed when "-h" or "--help" is used
        description=(
            "Creates an input file for Example 1 to be processed by the cpf package"
        ),
    )

    # Add (optional) argument to decide where to write the example input file to
    parser.add_argument(
        "-p",
        "--path",
        default=".",  # Save in current working directory if not called
        type=str,
        # Will be displayed when "help" is called
        help="Specify the directory in which the example input file is generated",
    )

    args = parser.parse_args()

    input_file = Path(args.path) / "Example_MultiPeak_Inputs.py"
    if input_file.exists():
        print("Example input file already exists")
    else:
        try:
            with open(input_file, "w", encoding="utf-8") as file:
                file.write("# Input parameters for BCC1_2GPa_100s_001.\n#\n\n\n")

                file.write("# properties of the data files.\n")
                file.write("datafile_directory = './Example1-Fe/'\n")
                file.write("datafile_Basename  = 'BCC1_2GPa_10s_001_'\n")
                file.write("datafile_Ending    = '.tif'\n")
                file.write("datafile_StartNum  = 1\n")
                file.write("datafile_EndNum    = 10\n")
                file.write("datafile_NumDigit  = 5\n\n")

                file.write("# Calibration and masking.\n")
                file.write("Calib_type   = 'Dioptas'\n")
                file.write("Calib_detector = 'Pilatus1M'\n")
                file.write(
                    "Calib_data     = datafile_directory + 'CeO2_Pil207_E30_2Nov2016_001.tif'\n"
                )
                file.write(
                    "Calib_param    = datafile_directory + 'CeO2_cal_Dioptas.poni'\n"
                )
                file.write(
                    "Calib_mask     = datafile_directory + 'DiffractionMask_Dioptas.mask'\n"
                )
                file.write("Calib_pixels = 172\n\n")

                file.write("# Fitting properties for peaks.\n")
                file.write("fit_orders = [\n")
                file.write("       {\n")
                file.write("         'background': [0,0],\n")
                file.write("         'peak': [{\n")
                file.write("             'phase': 'Other',\n")
                file.write("             'hkl': '000',\n")
                file.write("             'd-space': 3,\n")
                file.write("             'height': 1,\n")
                file.write("             'profile': 0,\n")
                file.write("             #'profile_fixed': 1,\n")
                file.write("             'width': 0,\n")
                file.write("             'symmetry': 2\n")
                file.write("           },{\n")
                file.write("             'phase': 'BCC-Fe',\n")
                file.write("             'hkl': 110,\n")
                file.write("             'd-space': 3,\n")
                file.write("             'height': 16,\n")
                file.write("             'profile': 0,\n")
                file.write("             #'profile_fixed': 1,\n")
                file.write("             'width': 0,\n")
                file.write("             'symmetry': 2\n")
                file.write("           }],\n")
                file.write("         'range': [[10.8, 11.67]],\n")
                file.write(
                    "         'PeakPositionSelection': [[1,-170.5,11.19],[1,-150.5,11.19],[1,-120.5,11.19],[1,-58.997,11.19],[1,59.289,11.19],[1,23.187,11.19],[1,23.212,11.19], [1,53.158,11.19], [1,73.246,11.19], [1,93.246,11.19], [1,123.246,11.19], [2,-170, 11.54], [2,-150, 11.54], [2,-110, 11.54], [2,-90, 11.54], [2,-45, 11.54], [2,-0, 11.54], [2,45, 11.54], [2,90, 11.54], [2,120, 11.54], [2,150, 11.54], [2,179, 11.54]]\n"
                )
                file.write("},\n")
                file.write("       {\n")
                file.write("         'background': [0,0],\n")
                file.write("         'peak': [{\n")
                file.write("             'phase': 'BCC-Fe',\n")
                file.write("             'hkl': 200,\n")
                file.write("             'd-space': 3,\n")
                file.write("             'height': 8,\n")
                file.write("             'profile': 0,\n")
                file.write("             #'profile_fixed': 1,\n")
                file.write("             'width': 0,\n")
                file.write("             'symmetry': 2\n")
                file.write("           }],\n")
                file.write("         'range': [[16.1,16.6]]\n")
                file.write("       },\n")
                file.write("       {\n")
                file.write("         'background': [0,0],\n")
                file.write("         'peak': [{\n")
                file.write("             'phase': 'BCC-Fe',\n")
                file.write("             'hkl': 211,\n")
                file.write("             'd-space': 3,\n")
                file.write("             'height': 8,\n")
                file.write("             'profile': 0,\n")
                file.write("             #'profile_fixed': 1,\n")
                file.write("             'width': 0,\n")
                file.write("             'symmetry': 2\n")
                file.write("           }],\n")
                file.write("         'range': [[19.8, 20.27]]\n")
                file.write("       },\n")
                file.write("      {\n")
                file.write("          'background': [1],\n")
                file.write("          'peak': [{\n")
                file.write("             'phase': 'BCC-Fe',\n")
                file.write("             'hkl': 220,\n")
                file.write("              'd-space': 3,\n")
                file.write("              'height': 8,\n")
                file.write("              'profile': 0,\n")
                file.write("              'profile_fixed': 0.5,\n")
                file.write("              'width': 0,\n")
                file.write("              'symmetry': 2\n")
                file.write("            }],\n")
                file.write(
                    "#          'PeakPositionSelection': [[1,-120.5,23.226],[1,-58.997,23.186],[1,59.289,23.208],[1,23.187,120.893],[1,23.212,6.367], [1,23.158,90.025], [1,23.246,25.729]]\n"
                )
                file.write("          'range': [[23.1, 23.5]]\n")
                file.write("        },\n")
                file.write("       {\n")
                file.write("         'background': [1,1],\n")
                file.write("         'peak': [{\n")
                file.write("             'phase': 'BCC-Fe',\n")
                file.write("             'hkl': 310,\n")
                file.write("             'd-space': 2,\n")
                file.write("             'height': 1,\n")
                file.write("             'profile': 0,\n")
                file.write("             'profile_fixed': 0.5,\n")
                file.write("             'width': 0,\n")
                file.write("             'symmetry': 2\n")
                file.write("           }],\n")
                file.write("         'range': [[25.7, 26.4]]\n")
                file.write("       },\n")
                file.write("    ]\n\n")

                file.write("#Number of bins for initial fitting.\n")
                file.write("AziBins = 90\n\n")

                file.write("#Output settings\n")
                file.write("Output_directory   = './'\n")
                file.write("Output_type        = 'MultiFit'\n")
                file.write("Output_NumAziWrite = 90\n\n")

                print(f"Created example input file at {input_file.resolve()}")

        except PermissionError:
            print("Write permission required")


# Run script if called directly
if __name__ == "__main__":
    run()
