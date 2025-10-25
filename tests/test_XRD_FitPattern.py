"""
Module to test the functions present in XRD_FitPattern.
"""

import os
from pathlib import Path

from cpf.XRD_FitPattern import execute
from pytest import fixture, mark


@fixture
def reset_working_directory():
    """
    The CPF functions involve a lot of changes in working directories. These will not
    get reset when running iterated tests, so this fixture funciton ensures that the
    working directory is reset to its initial state at the start of every test run.
    """

    cwd = Path().cwd().absolute()  # Save current working directory
    yield  # Test runs here
    os.chdir(cwd)  # Reset working directory


# Run the same test on different datasets and input files
execute_test_matrix = (
    # Dataset | Input file
    ("Example1-Fe", "BCC1_Dioptas_EmptyRanges_input.py"),
    ("Example1-Fe", "BCC1_Dioptas_EqualParams_input.py"),
    ("Example1-Fe", "BCC1_Dioptas_input.py"),
    ("Example1-Fe", "BCC1_Dioptas_MultiPeak_input.py"),
    # ("Example1-Fe", "BCC1_Dioptas_SeriesFunctions_input.py"),
    ("Example1-Fe", "BCC1_Dioptas_SymmFixed_input.py"),
    ("Example2-MgO", "CoSi22_MgO_DetectorPosition_input.py"),
    ("Example2-MgO", "CoSi22_MgO_input.py"),
    # ("Example2-MgO", "CoSi22_MgO_Reverse_input.py"),
    ("Example2-MgO", "CoSi22_MgO_Track_input.py"),
)


@mark.parametrize("test_params", execute_test_matrix)
def test_execute(
    test_params: tuple[str, str],
    reset_working_directory,
):
    """
    Tests the top-level 'execute()' function to confirm that the workflow runs
    through to completion in its current state.
    """

    # Unpack test params
    dataset, input_file = map(Path, test_params)

    # Set things up for the function
    os.chdir(dataset)
    assert input_file.exists()
    Path("results").mkdir(exist_ok=True)

    # Run the function
    execute(input_file)

    # 'execute()' doesn't return a helpful value at the moment, so we can assert True
    # at the end of the function for now to check that it ran to completion.
    assert True
