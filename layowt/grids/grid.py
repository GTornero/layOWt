"""This module contains the Grid class which is the basic building block used
to generate layouts."""

import matplotlib.pyplot as plt
import numpy as np
from layowt.grids.utils import create_coords


class Grid:
    """Grid class.

    This class is responsible for holding all of the attributes required in
    order to define a grid of points and to generate the coordinate pairs of
    the points that make up the grid.

    Multiple numerical and graphical methods are available to manipulate the grid
    and to visualise the coordinates.

    """

    def __init__(
        self,
        n_rows: int = 10,
        n_cols: int = 10,
        row_step: float = 6,
        col_step: float = 6,
        row_offset: bool = False,
        col_offset: bool = False,
        angle: float = 0,
        x_shear: float | None = None,
        y_shear: float | None = None,
        origin: tuple = (0, 0),
        scale: float = 1,
    ) -> None:
        """__init__ Initialise the Grid instance.

        Parameters
        ----------
        n_rows : int, optional
            Number of rows of points in the Grid, by default 10
        n_cols : int, optional
            Number of columns of points in the Grid, by default 10
        row_step : float, optional
            Non-dimensional distance between rows, by default 6
        col_step : float, optional
            Non-dimensional distance between columns, by default 6
        row_offset : bool, optional
            Offset every other row by half row_step to create a diamond shaped
            grid, by default False
        col_offset : bool, optional
            Offset every other column by half col_step to create a diamond
            shaped grid, by default False
        angle : float, optional
            Angle of orientation of the grid columns in degrees. Clockwise
            convention used to match wind direction convention, by default 0
        x_shear : float | None, optional
            Horizontal shear angle of the grid in degrees. Defines the angle of
            the columns from vertical clockwise, by default None
        y_shear : float | None, optional
            Vertical shear angle of the grid in degrees. Defines the angle of the
            rows from the horizontal anticlockwise, by default None
        origin : tuple, optional
            (x,y) coordinate pair of the centroid of the grid, by default (0, 0)
        scale : float, optional
            Factor by which to convert the non-dimensional row_step and col_step
            distances into the real space, by default 1
        """
        self.n_rows = n_rows
        self.n_cols = n_cols
        self.row_step = row_step
        self.row_step = n_cols
        self.col_step = col_step
        self.row_offset = row_offset
        self.col_offset = col_offset
        self.angle = angle
        self.x_shear = x_shear
        self.y_shear = y_shear
        self.origin = origin
        self.scale = scale

    @property
    def coords(self) -> np.ndarray:
        """coords Coordinates of the points in the grid based on the grid attributes

        Returns
        -------
        np.ndarray
            Array of [x,y] coordinate pairs of the points of the grid.
        """
        grid = create_coords(
            n_rows=self.n_rows,
            n_cols=self.n_cols,
            row_step=self.row_step,
            col_step=self.col_step,
            row_offset=self.row_offset,
            col_offset=self.col_offset,
            angle=self.angle,
            x_shear=self.x_shear,
            y_shear=self.y_shear,
        )

        grid *= self.scale
        grid[:, 0] += self.origin[0]
        grid[:, 1] += self.origin[1]

        return grid

    @property
    def x(self) -> np.ndarray:
        """x 1-D array of the x coordinates of the points in the grid.

        Returns
        -------
        np.ndarray
            1-D array of the x coordinates of the points in the grid.
        """
        return self.coords[:, 0]

    @property
    def y(self) -> np.ndarray:
        """y 1-D array of the y coordinates of the points in the grid.

        Returns
        -------
        np.ndarray
            1-D array of the y coordinates of the points in the grid.
        """
        return self.coords[:, 1]

    def plot(
        self, ax: plt.Axes | None = None, figsize: tuple = (8, 6)
    ) -> tuple[plt.Figure, plt.Axes]:  # type: ignore
        """plot Plot the grid of points.

        Parameters
        ----------
        ax : plt.Axes | None, optional
            Axes on which the grid points are drawn. If None, a new Figure and
            Axes object are created, by default None
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6)

        Returns
        -------
        tuple[plt.Figure, plt.Axes]
            Tuple of Figure and Axes objects

        Raises
        ------
        TypeError
            The passed ax argument was not of type plt.Axes
        """

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            try:
                fig = ax.get_figure()
            except AttributeError:
                raise TypeError from AttributeError(
                    f"Invalid type in {type(ax)} for the 'ax' argument, "
                    f"must be Matplotlib Axes object."
                )

        ax.scatter(self.x, self.y)
        ax.set_aspect("equal", "box")

        return fig, ax

    def set_scale(self, rotor_diameter: float, inplace: bool = False) -> "Grid" | None:
        """set_scale Set the scaling factor attribute used to scale the non-dimensional
        row and column spacings into the real space.

        Parameters
        ----------
        rotor_diameter : float
            rotor diameter to use as the scaling factor
        inplace : tuple, optional
            Perform operation inplace, by default False

        Returns
        -------
        Grid or None
            Returns a new Grid object with the modified scale attribute or
            None if inplace=True.
        """
        if inplace:
            self.scale = rotor_diameter
        else:
            return self.__class__(
                n_rows=self.n_rows,
                n_cols=self.n_cols,
                row_step=self.row_step,
                col_step=self.col_step,
                row_offset=self.row_offset,
                col_offset=self.col_offset,
                angle=self.angle,
                x_shear=self.x_shear,
                y_shear=self.y_shear,
                origin=self.origin,
                scale=rotor_diameter,
            )

    def translate(self, x: float, y: float, inplace: bool = False) -> "Grid" | None:
        """translate Translate the grid by shifting the origin

        Parameters
        ----------
        x : float
            Horizontal distance by which to shift the origin
        y : float
            Vertical distance by which to shift the origin
        inplace : bool, optional
            Perform the operation inplace, by default False

        Returns
        -------
        Grid or None
            Return the translated grid or None if inplace=True.
        """
        if inplace:
            self.origin = (self.origin[0] + x, self.origin[1] + y)
        else:
            return self.__class__(
                n_rows=self.n_rows,
                n_cols=self.n_cols,
                row_step=self.row_step,
                col_step=self.col_step,
                row_offset=self.row_offset,
                col_offset=self.col_offset,
                angle=self.angle,
                x_shear=self.x_shear,
                y_shear=self.y_shear,
                origin=(x, y),
                scale=self.scale,
            )

    def rotate(self, angle: float, inplace: bool = False):
        """rotate Rotate the

        Parameters
        ----------
        angle : float
            Angle in degrees by which to rotate the grid cloclwise
        inplace : bool, optional
            Perform the operation inplace, by default False

        Returns
        -------
        Grid or None
            Return the rotated grid or None if inplace=True.
        """
        if inplace:
            self.angle += angle
        else:
            return self.__class__(
                n_rows=self.n_rows,
                n_cols=self.n_cols,
                row_step=self.row_step,
                col_step=self.col_step,
                row_offset=self.row_offset,
                col_offset=self.col_offset,
                angle=self.angle + angle,
                x_shear=self.x_shear,
                y_shear=self.y_shear,
                origin=self.origin,
                scale=self.scale,
            )
