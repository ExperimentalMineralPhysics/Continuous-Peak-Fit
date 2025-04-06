"""
Module to test the functions present in XRD_FitPattern.
"""

import os
from pathlib import Path
from unittest.mock import patch

from cpf.IO_functions import make_outfile_name
from cpf.XRD_FitPattern import execute, initiate, write_output
from pytest import fixture, mark


@fixture
def reset_working_directory():
    """
    The CPF functions involve a lot of changes in working directories. These will not
    get reset when running iterated tests, so this fixture funciton ensures that the
    working directory is reset to its initial state at the start of every test run.
    """
    original_dir = os.getcwd()  # Save the current working directory
    yield  # Test runs here
    os.chdir(original_dir)  # After the test, reset the working directory


# Run the same test on different datasets and input files
# Q: How many of the Python input files are supposed to be compatible with 'execute()'?
execute_test_matrix = (
    # Dataset | Input file | Provide a Settings object? | Mode
    ("Example1-Fe", "BCC1_MultiPeak_input_Dioptas.py", False, "fit"),
    ("Example1-Fe", "BCC1_MultiPeak_input_Dioptas.py", False, "search"),
)


@mark.parametrize("test_params", execute_test_matrix)
@patch("cpf.XRD_FitPattern.write_output", wraps=write_output)  # spy_write_output
@patch(
    "cpf.XRD_FitPattern.make_outfile_name", wraps=make_outfile_name
)  # spy_make_outfile_name
@patch("cpf.XRD_FitPattern.initiate", wraps=initiate)  # spy_initiate
def test_execute(
    spy_initiate,
    spy_make_outfile_name,
    spy_write_output,
    test_params: tuple[str, str, bool, str],
    reset_working_directory,  # Will reset the working directory at the start of each run
):
    """
    Tests the top-level 'execute()' function to confirm that the workflow runs
    through to completion in its current state.
    """

    # Unpack test params
    dataset, input_file, use_settings_class, mode = (
        Path(test_params[0]),
        Path(test_params[1]),
        *test_params[2:],
    )

    # Set things up for the function
    os.chdir(dataset)
    (Path(".") / "results").mkdir(exist_ok=True)
    assert input_file.exists()

    # Run the function
    execute(
        settings_file=input_file,
        settings_class=None,
        mode=mode,
    )

    # 'initiate()' should be called if no settings file was provided
    if not use_settings_class:
        spy_initiate.assert_called_once_with(
            settings_file=input_file,
            inputs=None,
            report=False,
        )

    # Check that functions are called the right number of times in the different modes
    if mode == "fit":
        # This could be made more exact eventually
        assert spy_make_outfile_name.call_count > 1
        spy_write_output.assert_called_once()

    if mode == "search":
        # This could be made more exact eventually
        assert spy_make_outfile_name.call_count > 1

    # 'execute()' doesn't return a helpful value at the moment, so we 'assert True'
    # at the end of the function for now to check that it ran to completion.
    assert True
