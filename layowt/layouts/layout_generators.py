""" This module contains generator classes in charge of procedurally generating layouts.
"""

from collections import OrderedDict
from itertools import product
from typing import Callable, Iterable

import numpy as np
import rasterio
from layowt.grids import Grid
from layowt.layouts import FloatingLayout, Layout
from layowt.layouts.utils import layouts_to_hdf5
from shapely.geometry import MultiPolygon, Polygon


class GriddedLayoutGenerator:
    """
    This class takes care of procedurally generating layouts by iterating through different combinations of Grid parameters. \
    All resulting layouts will be generated from a uniform grid. \
    Can also take care of clipping layouts to buildable areas as well as avoiding exclusions if passed.
    
    Parameters
    ----------
    n_rows : int, optional
        Number of grid rows, by default 10.
    n_cols : int, optional
        Number of grid columns, by default 10.
    row_steps : Iterable[float], optional
        An iterable with the non-dimensional row distances to iterate through when creating gridded layouts, by default np.arange(5,10.5, 0.5).
    col_steps : Iterable[float], optional
        An iterable with the non-dimensional column distances to iterate through when creating gridded layouts, by default np.arange(5,10.5, 0.5).
    row_offset : bool, optional
        Offset every other row by half row_step to create a diamond shaped grid. Cannot be True when col_offset is also True. By default False.
    col_offset : bool, optional
        Offset every other column by half row_step to create a diamond shaped grid. Cannot be True when row_offset is also True. By default False.
    angles : Iterable[float], optional
        An iterable of angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1).
    x_shears : Iterable[float], optional
        An iterable of vertical shear angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1).
    y_shears : Iterable[float], optional
        An iterable of horizontal shear angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1).
    origins : list[tuple[float, float]], optional
        A list of (x, y) tuple coordinate pairs to set the grid origin. Will iterate over all Grid parameters at for all grid origins, by default [(0, 0)].
    scales : Iterable[float], optional
        An iterable of factors by which to convert the non-dimensional Grid distances into the real space. Will iterate over all passed scales for all Grid parameters. By default np.arange(250, 251).
    n_wtg : tuple(int, int) | None, optional
        Tuple of minimum and maximum target number of valid layout positions. Used to filter only the resulting Layouts with the desired number of positions. By default None.
    epsg : int | None, optional
        EPSG code for the layout coordinates.
    areas : list[Polygon  |  MultiPolygon], optional
        List of Shapely Polygons or MultiPolygons used to define the buildable areas. The generator wll take care to clip all layout positions to these areas. By default None.
    exclusions : list[Polygon  |  MultiPolygon], optional
        List of Shapely Polygons or MultiPolygons used to define exclusion zones where layout positions cannot be placed. The generator will take care to remove any layout positions that fall within these exclusion areas. By default None.
    bathymetry_path : str
        filepath of the bathymetry raster to load. Must be a format supported by the rasterio.open function.
    bathymetry_sign : str, optional
        sign in which the bathymetry file represents water depths. By default '-' if bathymetry values are negative. Pass '+' if bathymetry values are positive.
    bathymetry_band : int, optional
        Band number of the raster file containing bathymetry data, by default 1. Following the GDAL convention, bands are indexed from 1.
    bathymetry_limits : tuple[float, float], optional
        Minimum and maximum bathymetry values, in the raster file units, for layout positions to be considered valid, by default (0., 60.).
    bathymetry_drop_na : bool, optional
        If True, consider layout positions with nan values as invalid. By default False.
    
    Attributes
    ----------
    _grid_instructions : dict
        Attributes: Ordered dictionary of the parameters used to create all grids.
    n_wtg : tuple(int, int)
        Attribute: The (minimum, maximum) target number of layout positions.
    epsg : int
        Attribute: EPSG code for the layout coordinates.
    areas : list[Polygon | MultiPolygon]
        Attribute: List of buildable area geometries of the layouts.
    exclusions : list[Polygon | MultiPolygon]
        Attribute: List of exclusion area geometries of the layouts.
    bathymetry_path : str
        Attribute: filepath of the bathymetry dataset.
    bathymetry_limits: tuple(float, float)
        Attribute: Tuple of minimum and maximum allowable bathymetry values for valid layout positions.
    bathymetry_sign : str
        Attribute: Sign ("+" or "-") in which bathymetry values are represented.
    bathymetry_drop_na : bool
        Attributes: If layout positions at locations of missing bathymetry have been removed.
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
                 n_wtg: tuple[int, int] | None = None,
                 epsg: int | None = None,
                 areas: list[Polygon | MultiPolygon] | None = None,
                 exclusions: list[Polygon | MultiPolygon] | None = None,
                 bathymetry_path: str | None = None,
                 bathymetry_sign: str = "-",
                 bathymetry_band: int = 1,
                 bathymetry_limits: tuple[float, float] = (0., 60.),
                 bathymetry_drop_na: bool = False
                 ) -> None:
        """__init__ Initialises the GriddedLayoutGenerator instance.
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
        self.epsg = epsg
        self.areas = areas
        self.exclusions = exclusions
        self.bathymetry_path = bathymetry_path
        self.bathymetry_sign = bathymetry_sign
        self.bathymetry_band = bathymetry_band
        self.bathymetry_limits = bathymetry_limits
        self.bathymetry_drop_na = bathymetry_drop_na

    def _combinations(self):
        """
        Private helper method used to generate the different combinations of Grid parameters used to construct individual Grid objects.

        Yields
        ------
        tuple
            A tuple of all possible Grid parameter combinations.
        """
        instructions = []
        
        for key, val in self._grid_instructions.items():
            if key not in ["n_rows", "n_cols", "row_offset", "col_offset"]:
                instructions.append(val)
                
        for instruction in product(*instructions):
            yield instruction
    

    def __len__(self):
        return len(list(self._combinations()))
    
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
    
    def generate_layouts(self, progress: bool = False, to_file: bool = False, filepath: str = "layouts.h5", chunksize: int = 5000) -> list[Layout] | None:
        """
        This method returns a list of all the gridded Layouts generated by the GriddedLayoutGenerator object. \
        Alternatively saves the layouts to an hdf5 file.
        
        Parameters
        ----------
        progress : bool, optional
            Print progress, by default False.
        to_file : bool, optional
            Save results to file rather than keep in memory. Recommended useful to avoid memory leaks when
            generating many thousands of layouts. Will store all of the layouts in an hdf5 file and return None. By default False.
        filepath: str, optional
            Filepath of the hdf5 file to store the layouts in. If the file already exists, the new layouts will be appended.
            By default "layouts.h5".
        chunksize: int, optional
            Number of layouts to create before saving to disk, by default 5000.
        
        Returns
        -------
        list[Layout] | None
            List of generated gridded Layout objects. Alternatively None if writing to file.
        
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
        When passing parameters that results in hundreds of thousands (or more!) of individual grids, performance can suffer.
        Please consider writting outputs to file to avoid memory leaks.
        This generator is a work in progress and still under heay development.
        
        See Also
        --------
        layowt.layouts.layout.Layout : Layout object.
        layowt.grids.grid.Grid : Grid object.
        layowt.layouts.utils.layouts_to_hdf5: Function to write layouts to an hdf5 file.
        
        Warning
        --------
        Under heavy development.
        """
        layouts = []
        total = len(self)
        grid_gen = self._grid_generator()
        hdf5_id_counter = 0
        # TODO work on creating a parallel/multithreaded compution option for this function.
        if self.bathymetry_path is not None:
            with rasterio.open(self.bathymetry_path, mode='r') as dataset:
                for counter, grid in enumerate(grid_gen, start=1):
                    layout = Layout(grid, areas=self.areas, exclusions=self.exclusions, epsg=self.epsg)
                    
                    layout.apply_bathymetry(dataset,
                                            sign=self.bathymetry_sign,
                                            band=self.bathymetry_band,
                                            limits=self.bathymetry_limits,
                                            drop_na=self.bathymetry_drop_na
                                            )
                        
                    if self.n_wtg is None:
                        layouts.append(layout)
                    elif self.n_wtg[0] <= layout.n_wtg <= self.n_wtg[1]:
                        layouts.append(layout)
                    
                    if progress:
                        print(f"Progress: {counter}/{total}", end="\r")
                    
                    if to_file and (len(layouts) != 0) and ((len(layouts) % chunksize == 0) or (counter == total)):
                        if progress:
                            print(f"\nSaving {len(layouts)} layouts to file.")
                        layouts_to_hdf5(layouts, filepath=filepath, start_id=hdf5_id_counter)
                        hdf5_id_counter += chunksize
                        layouts = []
        else:
            for counter, grid in enumerate(grid_gen, start=1):
                layout = Layout(grid, areas=self.areas, exclusions=self.exclusions, epsg=self.epsg)
                    
                if self.n_wtg is None:
                    layouts.append(layout)
                elif self.n_wtg[0] <= layout.n_wtg <= self.n_wtg[1]:
                    layouts.append(layout)
                
                if progress:
                    print(f"Progress: {counter}/{total}", end="\r")
                
                if to_file and (len(layouts) != 0) and ((len(layouts) % chunksize == 0) or (counter == total)):
                    if progress:
                        print(f"\nSaving {len(layouts)} layouts to file.")
                    layouts_to_hdf5(layouts, filepath=filepath, start_id=hdf5_id_counter)
                    hdf5_id_counter += chunksize
                    layouts = []
                
        if not to_file:
            return layouts

class FloatingGriddedLayoutGenerator(GriddedLayoutGenerator):
    """
    This class takes care of procedurally generating floating layouts by iterating through different combinations of Grid parameters. \
    All resulting layouts will be generated from a uniform grid. \
    Can also take care of clipping layouts and moorings to buildable areas as well as avoiding exclusions if passed.
    
    Parameters
    ----------
    n_rows : int, optional
        Number of grid rows, by default 10.
    n_cols : int, optional
        Number of grid columns, by default 10.
    row_steps : Iterable[float], optional
        An iterable with the non-dimensional row distances to iterate through when creating gridded layouts, by default np.arange(5,10.5, 0.5).
    col_steps : Iterable[float], optional
        An iterable with the non-dimensional column distances to iterate through when creating gridded layouts, by default np.arange(5,10.5, 0.5).
    row_offset : bool, optional
        Offset every other row by half row_step to create a diamond shaped grid. Cannot be True when col_offset is also True. By default False.
    col_offset : bool, optional
        Offset every other column by half row_step to create a diamond shaped grid. Cannot be True when row_offset is also True. By default False.
    angles : Iterable[float], optional
        An iterable of angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1).
    x_shears : Iterable[float], optional
        An iterable of vertical shear angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1).
    y_shears : Iterable[float], optional
        An iterable of horizontal shear angles in degrees to iterate through when creating gridded layouts, by default np.arange(0, 1).
    origins : list[tuple[float, float]], optional
        A list of (x, y) tuple coordinate pairs to set the grid origin. Will iterate over all Grid parameters at for all grid origins, by default [(0, 0)].
    scales : Iterable[float], optional
        An iterable of factors by which to convert the non-dimensional Grid distances into the real space. Will iterate over all passed scales for all Grid parameters. By default np.arange(250, 251).
    n_wtg : tuple(int, int) | None, optional
        Tuple of minimum and maximum target number of valid layout positions. Used to filter only the resulting Layouts with the desired number of positions. By default None.
    epsg : int | None, optional
        EPSG code for the layout coordinates.
    areas : list[Polygon  |  MultiPolygon], optional
        List of Shapely Polygons or MultiPolygons used to define the buildable areas. The generator wll take care to clip all layout positions to these areas. By default None.
    exclusions : list[Polygon  |  MultiPolygon], optional
        List of Shapely Polygons or MultiPolygons used to define exclusion zones where layout positions cannot be placed. The generator will take care to remove any layout positions that fall within these exclusion areas. By default None.
    bathymetry_path : str
        filepath of the bathymetry raster to load. Must be a format supported by the rasterio.open function.
    bathymetry_sign : str, optional
        sign in which the bathymetry file represents water depths. By default '-' if bathymetry values are negative. Pass '+' if bathymetry values are positive.
    bathymetry_band : int, optional
        Band number of the raster file containing bathymetry data, by default 1. Following the GDAL convention, bands are indexed from 1.
    bathymetry_limits : tuple[float, float], optional
        Minimum and maximum bathymetry values, in the raster file units, for layout positions to be considered valid, by default (0., 60.).
    bathymetry_drop_na : bool, optional
        If True, consider layout positions with nan values as invalid. By default False.
    n_moorings : int
        Number of mooring lines at each floating platfrom.
    mooring_radius_func : float | Callable[[float], float]
        Mooring radius function. This parameter defines how the mooring radius length is
        calculated for each turbine floating platform. If a float is passed, it is used as a
        constant mooring length for all turbine positions. If a callable is passed, it must
        take in a single argument of bathymetry for each turbine postion.
    platform_bearing : float | None, optional
        Cardinal bearing of the floating platform with respect to True North, by default None.
    mooring_bearings : tuple[float, ...] | None, optional
        Cardinal bearings of the mooring lines in the local reference frame of the floating
        platform. If None is passed, the moorings will be spread evenly around the 360°
        anti-clockwise, by default None.
    mooring_offsets : list[tuple[float, float]] | None, optional
        The coordinate pairs for the mooring fairleads with respect to the turbine position in
        the floating platform local reference frame. If None is passed, all mooring fairleads
        will be at the turbine location, by default None.
    contain : bool, optional
        If True, the mooring lines and associated turbine locations outwith the mooring_areas
        are removed from the layout. If set to True but no mooring_areas are passed, the
        normal turbine area constraint geometry will be used, by default False.
    mooring_areas : list[Polygon  |  MultiPolygon] | None, optional
        Polygon geometries to constrain the mooring lines to, used in conjunction with the
        contain argument, by default None. Has no effect if contain is set to False.
    avoid_exclusions : bool, optional
        If True, removes the mooring lines and associated turbine locations that overlap with
        the mooring_exclusions. If set to True but no mooring_exclusions are passed, the normal
        turbine exclusion geometry will be used, by default False.
    mooring_exclusions : list[Polygon  |  MultiPolygon] | None, optional
        Polygon geometries to use as mooring exclusion zones, used in conjunction with the
        avoid_exclusions argument, by default None. Has no effect if avoid_exclusions is set to False.
    
    Attributes
    ----------
    _grid_instructions : dict
        Attributes: Ordered dictionary of the parameters used to create all grids.
    n_wtg : tuple(int, int)
        Attribute: The (minimum, maximum) target number of layout positions.
    epsg : int
        Attribute: EPSG code for the layout coordinates.
    areas : list[Polygon | MultiPolygon]
        Attribute: List of buildable area geometries of the layouts.
    exclusions : list[Polygon | MultiPolygon]
        Attribute: List of exclusion area geometries of the layouts.
    bathymetry_path : str
        Attribute: filepath of the bathymetry dataset.
    bathymetry_limits: tuple(float, float)
        Attribute: Tuple of minimum and maximum allowable bathymetry values for valid layout positions.
    bathymetry_sign : str
        Attribute: Sign ("+" or "-") in which bathymetry values are represented.
    bathymetry_drop_na : bool
        Attributes: If layout positions at locations of missing bathymetry have been removed.
    n_moorings : int
        Attribute: Number of mooring lines per turbine position.
    mooring_radius_func : float | Callable
        Attribute: Function used to calculate the mooring radius at each turbine position. Can be
        a constant float or a Callable that takes local turbine bathymetry as its sole argument.
    platform_bearing: float
        Attribute: Cardinal bearing of the floating platforms in degrees from True North.
    mooring_bearings: list[float]
        Attribute: List of cardinal bearings of the platform mooring lines in the global frame
        of reference.
    mooring_offsets: list[tuple(float, float)]
        Attribute: List of tuple coordinate pairs for the mooring fairlead locations with respect
        to the turbine location in the local platform reference frame.
    mooring_area: Polygon | MultiPolygon
        Attribute: The polygon(s) geometry which must contain the mooring lines.
    mooring_exclusions: Polygon | MultiPolygon
        Attribute: The polygon(s) geometries which the mooring lines cannot fall into.
    
    See Also
    --------
    layowt.layouts.layout_generators.GriddedLayuotGenerator : Layout Generator Base Class.
    layowt.layouts.layout.FloatingLayout : FloatingLayout object.
    layowt.grids.grid.Grid : Grid object.
    layowt.layouts.utils.layouts_to_hdf5: Function to write layouts to an hdf5 file.
        
    Warning
    --------
    Under heavy development.
    
    Current performance is slow if all constraints to turbines and moorings are required, and if bathymetry constraints and 
    mooring radius functions are used.
    
    Currently, it is recommended to not apply mooring constraints within this generator and to apply them to a filtered subset of layouts
    to ease computation times. The recommended workflow is the following:
    
    1. Generate layouts only with area, exclusion, and bathymetry constraints
    2. Filter the output list to a (hopefully much smaller) subset of layous of interest discarding the rest
    3. Generate moorings for the remaining layouts
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
                 epsg: int | None = None,
                 areas: list[Polygon | MultiPolygon] | None = None,
                 exclusions: list[Polygon | MultiPolygon] | None = None,
                 bathymetry_path: str | None = None,
                 bathymetry_sign: str = "-",
                 bathymetry_band: int = 1,
                 bathymetry_limits: tuple[float, float] = (0., 60.),
                 bathymetry_drop_na: bool = False,
                 n_moorings: int | None = None,
                 mooring_radius_func: float | Callable[[float], float] | None = None,
                 platform_bearing: float | None = None,
                 mooring_bearings: tuple[float, ...] | None = None,
                 mooring_offsets: list[tuple[float, float]] | None = None,
                 contain: bool = False,
                 mooring_areas: list[Polygon | MultiPolygon] | None = None,
                 avoid_exclusions: bool = False,
                 mooring_exclusions: list[Polygon | MultiPolygon] | None = None
                 ) -> None:
        super().__init__(n_rows=n_rows,
                         n_cols=n_cols,
                         row_steps=row_steps,
                         col_steps=col_steps,
                         row_offset=row_offset,
                         col_offset=col_offset,
                         angles=angles,
                         x_shears=x_shears,
                         y_shears=y_shears,
                         origins=origins,
                         scales=scales,
                         n_wtg=n_wtg,
                         epsg=epsg,
                         areas=areas,
                         exclusions=exclusions,
                         bathymetry_path=bathymetry_path,
                         bathymetry_sign=bathymetry_sign,
                         bathymetry_band=bathymetry_band,
                         bathymetry_limits=bathymetry_limits,
                         bathymetry_drop_na=bathymetry_drop_na
                         )
        self.n_moorings = n_moorings
        self.mooring_radius_func = mooring_radius_func
        self.platform_bearing = platform_bearing
        self.mooring_bearings = mooring_bearings
        self.mooring_offsets = mooring_offsets
        self.contain = contain
        self.mooring_areas = mooring_areas
        self.avoid_exclusions = avoid_exclusions
        self.mooring_exclusions = mooring_exclusions
    
    def generate_layouts(self, progress: bool = False, to_file: bool = False, filepath: str = "layouts.h5", chunksize: int = 5000) -> list[FloatingLayout] | None:
        """
        This method returns a list of all the gridded FloatingLayouts generated by the FloatingGriddedLayoutGenerator object. \
        Alternatively saves the layouts to an hdf5 file.
        
        Parameters
        ----------
        progress : bool, optional
            Print progress, by default False.
        to_file : bool, optional
            Save results to file rather than keep in memory. Recommended useful to avoid memory leaks when
            generating many thousands of layouts. Will store all of the layouts in an hdf5 file and return None. By default False.
        filepath: str, optional
            Filepath of the hdf5 file to store the layouts in. If the file already exists, the new layouts will be appended.
            By default "layouts.h5".
        chunksize: int, optional
            Number of layouts to create before saving to disk, by default 5000.
        
        Returns
        -------
        list[FloatingLayout] | None
            List of generated gridded FloatingLayout objects. Alternatively None if writing to file.
        
        Note
        ----
        When passing parameters that results in hundreds of thousands (or more!) of individual grids, performance can suffer.
        Please consider writting outputs to file to avoid memory leaks.
        This generator is a work in progress and still under heay development.
        
        See Also
        --------
        layowt.layouts.floating_layout.FloatingLayout : FloatingLayout object.
        layowt.grids.grid.Grid : Grid object.
        layowt.layouts.utils.layouts_to_hdf5: Function to write layouts to an hdf5 file.
        
        Warning
        --------
        Under heavy development.
        """
        
        layouts = []
        total = len(self)
        grid_gen = self._grid_generator()
        moorings = False
        hdf5_id_counter = 0
        if self.n_moorings is not None and self.mooring_radius_func is not None:
            moorings = True
        
        if self.bathymetry_path is not None:
            with rasterio.open(self.bathymetry_path, mode='r') as dataset:
                for counter, grid in enumerate(grid_gen, start=1):
                    layout = FloatingLayout(grid, areas=self.areas, exclusions=self.exclusions, epsg=self.epsg)
                    
                    layout.apply_bathymetry(dataset,
                                            sign=self.bathymetry_sign,
                                            band=self.bathymetry_band,
                                            limits=self.bathymetry_limits,
                                            drop_na=self.bathymetry_drop_na
                                            )
                    
                    if moorings:
                        layout.generate_moorings(n_moorings=self.n_moorings,
                                                 radius=self.mooring_radius_func,
                                                 platform_bearing=self.platform_bearing,
                                                 mooring_bearings=self.mooring_bearings,
                                                 mooring_offsets=self.mooring_offsets,
                                                 contain=self.contain,
                                                 mooring_areas=self.mooring_areas,
                                                 avoid_exclusions=self.avoid_exclusions,
                                                 mooring_exclusions=self.mooring_exclusions
                                                 )
                        
                    if self.n_wtg is None:
                        layouts.append(layout)
                    elif self.n_wtg[0] <= layout.n_wtg <= self.n_wtg[1]:
                        layouts.append(layout)
                    
                    if progress:
                        print(f"Progress: {counter}/{total}", end="\r")
                    
                    if to_file and (len(layouts) != 0) and ((len(layouts) % chunksize == 0) or (counter == total)):
                        if progress:
                            print(f"\nSaving {len(layouts)} layouts to file.")
                        layouts_to_hdf5(layouts, filepath=filepath, start_id=hdf5_id_counter)
                        hdf5_id_counter += chunksize
                        layouts = []
        else:
            for counter, grid in enumerate(grid_gen, start=1):
                layout = FloatingLayout(grid, areas=self.areas, exclusions=self.exclusions, epsg=self.epsg)
                
                if moorings:
                    layout.generate_moorings(n_moorings=self.n_moorings,
                                             radius=self.mooring_radius_func,
                                             platform_bearing=self.platform_bearing,
                                             mooring_bearings=self.mooring_bearings,
                                             mooring_offsets=self.mooring_offsets,
                                             contain=self.contain,
                                             mooring_areas=self.mooring_areas,
                                             avoid_exclusions=self.avoid_exclusions,
                                             mooring_exclusions=self.mooring_exclusions
                                             )
                    
                if self.n_wtg is None:
                    layouts.append(layout)
                elif self.n_wtg[0] <= layout.n_wtg <= self.n_wtg[1]:
                    layouts.append(layout)
                
                if progress:
                    print(f"Progress: {counter}/{total}", end="\r")
                
                if to_file and (len(layouts) != 0) and ((len(layouts) % chunksize == 0) or (counter == total)):
                    if progress:
                        print(f"\nSaving {len(layouts)} layouts to file.")
                    layouts_to_hdf5(layouts, filepath=filepath, start_id=hdf5_id_counter)
                    hdf5_id_counter += chunksize
                    layouts = []
        
        return layouts
