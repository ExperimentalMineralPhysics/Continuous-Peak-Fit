"""
Module to test the metadata functions present in XRD_FitPattern.
"""

import cpf
import os
from pathlib import Path

# Run the same test on different datasets and input files
execute_test_matrix = (
    # Dataset | Input file
    ("Example1-Fe", "BCC1_Dioptas_EmptyRanges_input.py", "default"),
    ("Example1-Fe", "BCC1_Dioptas_input.py", "default"),
    ("Example1-Fe", "BCC1_Dioptas_MultiPeak_input.py", "default"),
    ("Example1-Fe", "BCC1_Dioptas_SymmFixed_input.py", "default"),
    ("Example1-Fe", "BCC1_Dioptas_SeriesFunctions_input.py", "default"),
    ("Example1-Fe", "BCC1_Dioptas_EqualParams_input.py", "all"),
    ("Example2-MgO", "CoSi22_MgO_input.py", ['mean_start_time', 'mean_live_time', '6BMB_LVP:LVP_tc1_calcs.I', '6BMB_LVP:LVP_tc2_calcs.I', 'FILE_CREATION', 'time_start']),
    ("Example2-MgO", "CoSi22_MgO_Track_input.py", ['mean_start_time', 'mean_live_time', '6BMB_LVP:LVP_tc1_calcs.I', '6BMB_LVP:LVP_tc2_calcs.I', 'FILE_CREATION', 'time_start']),
    ("Example2-MgO", "CoSi22_MgO_Reverse_input.py", ['mean_start_time', 'mean_live_time', '6BMB_LVP:LVP_tc1_calcs.I', '6BMB_LVP:LVP_tc2_calcs.I', 'FILE_MODIFIED', 'time_start']),
    ("Example2-MgO", "CoSi22_MgO_DetectorPosition_input.py", 'all'),
)


# @fixture
# def reset_working_directory():
#     """
#     The CPF functions involve a lot of changes in working directories. These will not
#     get reset when running iterated tests, so this fixture funciton ensures that the
#     working directory is reset to its initial state at the start of every test run.
#     """

#     cwd = Path().cwd().absolute()  # Save current working directory
#     yield  # Test runs here
#     os.chdir(cwd)  # Reset working directory


for i in execute_test_matrix:

    cwd = Path().cwd().absolute()  # Save current working directory

    os.chdir("../" + i[0])

    # Unpack test params
    input_file = i[1]
    
    settings_class = cpf.XRD_FitPattern.initiate(input_file)
    settings_class.set_subpattern(0,0)
    
    #get data from settings class
    new_data = settings_class.data_class
    new_data.fill_data(settings=settings_class)
    
    for n, ii in enumerate(settings_class.image_list):
        meta = new_data.get_metadata()
        print(meta)
        
            
    for n, ii in enumerate(settings_class.image_list):
        meta = new_data.get_metadata(metadata_values=i[2])
        print(meta)
        
    os.chdir(cwd)  # Reset working directory