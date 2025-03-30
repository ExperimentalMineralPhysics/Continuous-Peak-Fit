"""
Module to test the functions present in XRD_FitPattern.
"""

import os
from pathlib import Path

from cpf.XRD_FitPattern import execute
from pytest import mark

# Run the same test on different datasets and input files
execute_test_matrix = (
    # Dataset | Input file
    ("Example1-Fe", "BCC1_MultiPeak_input_Dioptas.py"),
)


@mark.parametrize("test_params", execute_test_matrix)
def test_execute(test_params: tuple[str, str]):
    """
    Tests the top-level 'execute()' function to confirm that the workflow runs
    through to completion in its current state.
    """

    # Unpack test params
    dataset, input_file = test_params

    # Set things up for the function
    os.chdir(dataset)
    assert Path(input_file).exists()

    # Run the function
    execute(Path(input_file))

    # Execute doesn't return a helpful value at the moment, so we can assert True
    # at the end of the function to check that it ran through to completion.
    assert True
