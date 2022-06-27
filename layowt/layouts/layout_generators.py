""" This module contains generator classes in charge of procedurally generating layouts.
"""

from collections import OrderedDict
from itertools import product
from typing import Iterable

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from ..grids import Grid
from .layout import Layout


class GriddedLayoutGenerator:#
    """GriddedLayoutGenerator This class takes care of procedurally generating layouts by iterating through different combinations of Grid parameters. All resulting layouts will be generated from a uniform grid. Can also take care of clipping layouts to buildable areas as well as avoiding exclusions if passed.
    """
    def __init__(self,
                 n_rows: int = 10,
                 n_cols: int = 10,
                 row_steps: Iterable[float] = np.arange(5,10.5, 0.5),
                 col_steps: Iterable[float] = np.arange(5,10.5, 0.5),
                 row_offset: bool = False,
                 col_offset: bool = False,
                 angles: Iterable[float] = np.arange(0, 1),
                 x_shears: Iterable[float] = np.arange(0, 1),
                 y_shears: Iterable[float] = np.arange(0, 1),
                 origins: list[tuple[float, float]] = [(0, 0)],
                 scales: Iterable[float] = np.arange(250, 251),
                 n_wtg: int | None = None,
                 areas: list[Polygon | MultiPolygon] | None = None,
                 exclusions: list[Polygon | MultiPolygon] | None = None
                 ) -> None:
        """__init__ Initialises the GriddedLayoutGenerator instance.

        Parameters
        ----------
        n_rows : int, optional
            Number of grid rows, by default 10
        n_cols : int, optional
            Number of grid columns, by default 10
        row_steps : Iterable[float], optional
            An iterable with the non-dimensional row distances to iterate through when creating gridded layouts, by default np.arange(5,10.5, 0.5)
        col_steps : Iterable[float], optional
            An iterable with the non-dimensional column distances to iterate through when creating gridded layouts, by default np.arange(5,10.5, 0.5)
        row_offset : bool, optional
            Offset every other row by half row_step to create a diamond shaped grid. Cannot be True when col_offset is also True. By default False
        col_offset : bool, optional
            Offset every other column by half row_step to create a diamond shaped grid. Cannot be True when row_offset is also True. By default False
        angles : Iterable[float], optional
            An iterable of angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1)
        x_shears : Iterable[float], optional
            An iterable of vertical shear angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1)
        y_shears : Iterable[float], optional
            An iterable of horizontal shear angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1)
        origins : list[tuple[float, float]], optional
            A list of (x, y) tuple coordinate pairs to set the grid origin. Will iterate over all Grid parameters at for all grid origins, by default [(0, 0)]
        scales : Iterable[float], optional
            An iterable of factors by which to convert the non-dimensional Grid distances into the real space. Will iterate over all passed scales for all Grid parameters. By default np.arange(250, 251)
        n_wtg : int | None, optional
            Target number of valid layout positions. Used to filter only the resulting Layouts with the desired number of positions. By default None
        areas : list[Polygon  |  MultiPolygon], optional
            List of Shapely Polygons or MultiPolygons used to define the buildable areas. The generator wll take care to clip all layout positions to these areas. By default None
        exclusions : list[Polygon  |  MultiPolygon], optional
            List of Shapely Polygons or MultiPolygons used to define exclusion zones where layout positions cannot be placed. The generator will take care to remove any layout positions that fall within these exclusion areas. By default None
        """
        instructions = OrderedDict(
            {
            "n_rows": n_rows,
            "n_cols": n_cols,
            "row_steps": row_steps,
            "col_steps": col_steps,
            "row_offset": row_offset,
            "col_offset": col_offset,
            "angles": angles,
            "x_shears": x_shears,
            "y_shears": y_shears,
            "origins": origins,
            "scales": scales
        })
        
        self._grid_instructions = instructions
        self.n_wtg = n_wtg
        self.areas = areas
        self.exclusions = exclusions

    def _combinations(self):
        """_combinations Private helper method used to generate the different combinations of Grid parameters used to construct individual Grid objects.

        Yields
        ------
        tuple
            A tuple of all possible Grid parameter combinations.
        """
        instructions = []
        
        for k, v in self._grid_instructions.items():
            if k not in ["n_rows", "n_cols", "row_offset", "col_offset"]:
                instructions.append(v)
                
        for instruction in product(*instructions):
            yield instruction
    
    def _grid_generator(self):
        """_grid_generator Private helper method used to generate Grid objects from the Grid parameter combinations generated by the _combinations method.

        Yields
        ------
        Grid
            Returns a Grid object for the passed Grid parameters.
        """
        for instruction in self._combinations():
            yield Grid(
                n_rows = self._grid_instructions["n_rows"],
                n_cols = self._grid_instructions["n_cols"],
                row_step = instruction[0],
                col_step = instruction[1],
                row_offset = self._grid_instructions["row_offset"],
                col_offset = self._grid_instructions["col_offset"],
                angle = instruction[2],
                x_shear = instruction[3],
                y_shear = instruction[4],
                origin = instruction[5],
                scale = instruction[6]
            )
    
    def generate_layouts(self) -> list[Layout]:
        """generate_layouts This method returns a list of all the gridded Layouts generated by the GriddedLayoutGenerator object.

        Returns
        -------
        list[Layout]
            List of generated gridded Layout objects.
        
        Examples
        --------
        >>> layout_gen = GriddedLayoutGenerator(
                                    n_cols = 20,
                                    n_rows = 20,
                                    row_steps = np.arange(5, 6, 0.25),
                                    col_steps = np.arange(5, 6, 0.25),
                                    angles = np.arange(0, 90, 5),
                                    x_shears = np.arange(0, 45, 5),
                                    origins = [(0, 0) ,(2000, 5000)],
                                    scales = [250]
                                    )
        
        >>> layouts = layout_gen.generate_layouts()
        >>> len(layouts)
        5184
        
        This example shows how to create a GriddedLayoutGenerator instance and use it to automatically generate gridded layouts. This example will create layotus from a 20 by 20 grid, with row and column distances from 5 to 6 in steps of 0.25, angles from 0 to 90 in steps of 5 degrees, vertical shear angles from 0 to 45 in steps of 5 degrees, two different origins, and a scaling factor of 250.
        
        A total of 5184 layouts have been automatically generated.
        
        Note
        ----
        When passing parameters that results in hundreds of thousands (or more!) of individual grids, performance can suffer. This generator is a work in progress and still under heay development.
        
        See Also
        --------
        layowt.layouts.layout.Layout : Layout object.
        layowt.grids.grid.Grid : Grid object.
        
        Warning
        --------
        Under heavy developmnet.
        """
        layouts = []
        grid_gen = self._grid_generator()
        for grid in grid_gen:
            layout = Layout(grid, areas=self.areas, exclusions=self.exclusions)
            if self.n_wtg is None:
                layouts.append(layout)
            elif layout.n_wtg == self.n_wtg:
                layouts.append(layout)
        
        return layouts
    
