from unittest.mock import MagicMock

import pytest
from cpf.util.output_formatters import figure_suptitle_space


@pytest.mark.parametrize(
    "test_params",
    (  # Fig size | List of axes bounds (list of (x0, y0, width, height)) | Top margin
        ((4, 3), [(0.1, 0.1, 0.8, 0.8)], 0.2),  # 1 plot
        (  # 2 vertical plots
            (4, 3),
            [
                (0.1, 0.1, 0.8, 0.4),
                (0.1, 0.55, 0.8, 0.4),
            ],
            0.2,
        ),
        (  # 3 vertical plots
            (4, 3),
            [
                (0.1, 0.1, 0.8, 0.25),
                (0.1, 0.4, 0.8, 0.25),
                (0.1, 0.7, 0.8, 0.25),
            ],
            0.2,
        ),
        (  # 2 horizontal plots
            (4, 3),
            [
                (0.1, 0.1, 0.4, 0.8),
                (0.55, 0.1, 0.4, 0.8),
            ],
            0.2,
        ),
        (  # 3 horizontal plots
            (4, 3),
            [
                (0.1, 0.1, 0.25, 0.8),
                (0.4, 0.1, 0.25, 0.8),
                (0.7, 0.1, 0.25, 0.8),
            ],
            0.2,
        ),
        (  # 2 x 2 plot
            (4, 3),
            [
                (0.1, 0.1, 0.4, 0.4),
                (0.1, 0.55, 0.4, 0.4),
                (0.55, 0.1, 0.4, 0.4),
                (0.55, 0.55, 0.4, 0.4),
            ],
            0.2,
        ),
    ),
)
def test_figure_suptitle_space(
    test_params: tuple[tuple[int, int], list[tuple[float, float, float, float]], float],
):
    # Unpack test params
    fig_size, bounds_list, margin_top = test_params

    # Create a list of mock Axes objects
    axes_list = []
    for bounds in bounds_list:
        mock_axis = MagicMock()
        mock_axis.get_position.return_value.bounds = bounds
        axes_list.append(mock_axis)

    # Create a mock Figure object
    mock_figure = MagicMock(axes=axes_list)
    mock_figure.get_size_inches.return_value = fig_size

    # Run the function with the mock Figure and the specified margin
    figure_suptitle_space(mock_figure, topmargin=margin_top)

    # Check that the expected calls were made
    mock_figure.get_size_inches.assert_called_once()

    # Check that the correct calculations were made for each set of axes
    height_new = fig_size[1] - margin_top
    for (al0, ab0, aw0, ah0), mock_axis in zip(bounds_list, axes_list):
        mock_axis.set_position.assert_called_once_with(
            (
                al0,  # No change to left
                ab0 * (height_new / fig_size[1]),  # Rescaled bottom
                aw0,  # No change to width
                ah0 * (height_new / fig_size[1]),  # Rescaled height
            )
        )
    pass
