"""
This module contains utility functions for the creation of grid coordinates
for the layowt.grids.grid.Grid object.
"""


from functools import lru_cache

import numpy as np


def _create_base_coords(
    n_rows: int,
    n_cols: int,
    row_step: float,
    col_step: float,
    row_offset: bool,
    col_offset: bool,
) -> np.ndarray:
    """_create_base_coords Creates a simple rectangular grid of points centered around (0, 0) in
    non-dimensional space.

    Parameters
    ----------
    n_rows : int
        Number of rows of points.
    n_cols : int
        Number of columns of points.
    row_step : float
        Distance between rows of points.
    col_step : float
        Distance between columns of points.
    row_offset : bool
        Offset every other row by half the row_step to create a diamond shape
        grid.
    col_offset : bool
        Offset every other columns by half the col_step to create a diamond
        shape grid.

    Returns
    -------
    numpy.ndarray
        Returns an ndarray of the [x,y] coordinate pairs of points.
    """

    x_max = col_step * (n_cols - 1) / 2
    x_list = np.arange(-x_max, x_max + col_step, col_step)

    y_max = row_step * (n_rows - 1) / 2
    y_list = np.arange(-y_max, y_max + row_step, row_step)

    x_grid, y_grid = np.meshgrid(x_list, y_list)

    if row_offset:
        x_grid[::2] += col_step / 2
    elif col_offset:
        y_grid[:, ::2] += row_step / 2

    return np.column_stack((x_grid.flatten(), y_grid.flatten()))


def _rotate_coords(coords: np.ndarray, angle: float) -> np.ndarray:
    """_rotate_coords Rotate a grid of points by an angle in degrees
    clockwise.

    Same convention as wind direction.

    Parameters
    ----------
    coords : np.ndarray
        Array of [x,y] coordinate pairs.
    angle : float
        Angle in degrees clockwise.

    Returns
    -------
    np.ndarray
        Rotated grid of points as an array of [x,y] coordinate pairs.
    """
    rotation = _rotation_matrix(-1 * angle)
    rotated = np.matmul(rotation, coords.T)
    rotated = rotated.T

    return rotated


@lru_cache
def _rotation_matrix(angle: float) -> np.ndarray:
    """_rotation_matrix Calculates the rotation matrix from an input angle in
    degrees. Positive anticlockwise.

    Parameters
    ----------
    angle : float
        Angle in degrees, positive anticlockwise.

    Returns
    -------
    np.ndarray
        Rotation matrix.
    """
    return np.array(
        [
            [np.cos(angle * np.pi / 180), -np.sin(angle * np.pi / 180)],
            [np.sin(angle * np.pi / 180), np.cos(angle * np.pi / 180)],
        ]
    )


def _shear_coords(coords: np.ndarray, x_shear: float, y_shear: float) -> np.ndarray:
    """_shear_coords Shears a grid of points based on horizontal and vertical
    shear angles.

    Parameters
    ----------
    coords : np.ndarray
        Array of [x,y] coordinate pairs to shear.
    x_shear : float
        Horizontal shear angle in degrees.
    y_shear : float
        Vertical shear angle in degrees.

    Returns
    -------
    np.ndarray
        Sheared grid of points as an array of [x,y] coordinate pairs.
    """
    shear = _shear_matrix(x_shear, y_shear)

    sheared = np.matmul(shear, coords.T)
    sheared = sheared.T

    return sheared


@lru_cache
def _shear_matrix(x_shear: float, y_shear: float) -> np.ndarray:
    """_shear_matrix Calculates the shear matrix based on horizontal and
    vertical shear angles.

    Parameters
    ----------
    x_shear : float
        Horizontal shear angle in degrees.
    y_shear : float
        Vertical shear angle in degrees.

    Returns
    -------
    np.ndarray
        Shear matrix.
    """
    return np.array(
        [
            [1, np.arctan(x_shear * np.pi / 180)],
            [np.arctan(y_shear * np.pi / 180), 1],  # type: ignore
        ]
    )


def create_coords(
    n_rows: int,
    n_cols: int,
    row_step: float,
    col_step: float,
    row_offset: bool = False,
    col_offset: bool = False,
    angle: float | None = None,
    x_shear: float | None = None,
    y_shear: float | None = None,
) -> np.ndarray:
    """create_coords Creates a grid of points centered around (0, 0) in non-dimensional space.

    Parameters
    ----------
    n_rows : int
        Number of rows of points in the grid.
    n_cols : int
        Number of columns of points in the grid.
    row_step : float
        Distance between rows in the grid.
    col_step : float
        Distance between columns in the grid.
    row_offset : bool, optional
        Shifts every other row by half the row_step to create a diamond shaped
        grid, by default False
    col_offset : bool, optional
        Shifts every other column by half the col_step to create a diamond
        shaped grid, by default False
    angle : float | None, optional
        Angle of orientation of the grid columns in degrees. Clockwise
        convention used to match wind direction convention, by default None
    x_shear : float | None, optional
        Horizontal shear angle of the grid in degrees. Defines the angle of
        the columns from vertical clockwise, by default None
    y_shear : float | None, optional
        Vertical shear angle of the grid in degrees. Defines the angle of the
        rows from the horizontal anticlockwise, by default None

    Returns
    -------
    numpy.ndarray
        Array of [x,y] coordinate pairs of all points in the grid.

    Raises
    ------
    ValueError
        When both row_offset and col_offset are True.
    """

    if row_offset and col_offset:
        raise ValueError("Both row_offset and col_offset cannot be True.")

    coords = _create_base_coords(
        n_rows=n_rows,
        n_cols=n_cols,
        row_step=row_step,
        col_step=col_step,
        row_offset=row_offset,
        col_offset=col_offset,
    )

    if x_shear is not None or y_shear is not None:
        coords = _shear_coords(
            coords, x_shear=float(x_shear or 0), y_shear=float(y_shear or 0)
        )

    if angle is not None:
        coords = _rotate_coords(coords, angle=angle)

    return coords
