"""This module contains the Layout class which is the basic data structure used
to define project scenario geometries."""

import copy

import fiona
import geopandas as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import rasterio
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from py_wake.wind_turbines import WindTurbine
from rasterio.plot import show
from shapely.geometry import (MultiPoint, MultiPolygon, Point, Polygon,
                              mapping, shape)
from shapely.ops import transform, unary_union
from sqlalchemy import create_engine

from ..grids import Grid


class Layout:
    """
    This class is responsible for holding all of the data structures required
    in order to define a project scenario and to generate the coordinate pairs
    of the points that make up the grid.

    Multiple numerical, spatial and graphical methods are available to
    manipulate the scenario o conform to constraint and set properties.
    
    Parameters
    ----------
    grid : Grid | None, optional
        Grid object used to define the Layout, by default None.
    areas : list[Polygon  |  MultiPolygon] | None, optional
        List of Shapely Polygons or MultiPolygons used to define the lease
        areas, by default None.
    exclusions : list[Polygon  |  MultiPolygon] | None, optional
        List of Shapely Polygons or MultiPolygons used to define exclusion
        zones where turbines cannot be placed, by default None.
    wtg : str
        Filepath for the .wtg file to load, by default None.
    
    See Also
    --------
    Polygon : Shapely Polygon object.
    MultiPolygon : Shapely MultiPolygon object.
    py_wake.wind_turbines.WindTurbines : Set of multiple py_wake wind turbines.
    """

    def __init__(
        self,
        grid: Grid | None = None,
        areas: list[Polygon | MultiPolygon] | None = None,
        exclusions: list[Polygon | MultiPolygon] | None = None,
        wtg: str | None = None,
    ) -> None:
        """__init__ Initialises the Layout class instance.
        """
        self.grid = grid
        self.bathymetry_path = None
        self.bathymetry_limits = None
        self.bathymetry_sign = None
        self.bathymetry_drop_na = None

        if grid is not None:
            self.geom = self.grid.to_multipoint()  # type: ignore
            self._raw_geom = copy.deepcopy(self.geom)
        else:
            self.geom = None
            self._raw_geom = None

        if (areas is None) and (exclusions is None):
            self._constrained = False
        else:
            self._constrained = True

        if areas is not None:
            self.area = unary_union(areas)
            self.geom = self.geom.intersection(self.area)  # type: ignore
        else:
            self.area = None

        if exclusions is not None:
            self.exclusion = unary_union(exclusions)
            self.geom = self.geom.difference(self.exclusion)  # type: ignore
        else:
            self.exclusion = None

        if wtg is not None:
            self.wtg = WindTurbine.from_WAsP_wtg(wtg)
        else:
            self.wtg = None

    @classmethod
    def from_shapefile(cls, filepath: str, target_epsg: int | None = None) -> "Layout":
        """from_shapefile Alternate constructor to read a shapefile into a Layout object.
        The shapefile must contain Point or MultiPoint geometry to be a valid Layout constructor. Can reproject geometries on the fly from the source CRS into the desired CRS defined by its EPSG code in the optional target_epsg argument.

        Parameters
        ----------
        filepath : str
            Valid string file path to the shapefile to read.
        target_epsg : int | None, optional
            EPSG code of the target projection for the geometries to be loaded in. By default, None.

        Returns
        -------
        Layout
            Layout object created from the shapefile.
        """
        with fiona.open(filepath) as src:
            geoms = [shape(rec["geometry"]) for rec in src]
            src_crs = pyproj.CRS(src.crs)
        
        if target_epsg is not None:
            target_crs = pyproj.CRS("EPSG:" + str(target_epsg))
            crs_transformer = pyproj.Transformer.from_crs(src_crs, target_crs, always_xy=True)
            transformed_geoms = []
            for geom in geoms:
                transformed_geom = transform(crs_transformer.transform, geom)
                transformed_geoms.append(transformed_geom)
            
            geoms = transformed_geoms

        if all(isinstance(x, (Point, MultiPoint)) for x in geoms):
            geoms = unary_union(geoms)
            layout = cls()
            layout._raw_geom = geoms
            layout.geom = geoms
            return layout
        else:
            raise TypeError("Geometry type must be Point or MultiPoint.")

    @classmethod
    def from_text(
        cls, filepath: str, x_header: str | int, y_header: str | int, source_epsg: int | None = None, target_epsg: int | None = None, **kwargs
    ) -> "Layout":
        """from_text Alternate constructor to read a text file into a Layout object. Text file must contain x and y coordinates to create Point

        Parameters
        ----------
        filepath : str
            Valid string path to the text file to read.
        x_header : str | int
            Column name or number to read the x-coordinates.
        y_header : str | int
            Column name or number to read the y-coordinates.
        source_epsg : int | None, optional
            EPSG code of the source projection for the points to be loaded in. Used alongside target_epsg to compute the coordinate transformation. By default, None.
        target_epsg : int | None, optional
            EPSG code of the target projection for the points to be loaded in. Used alongside source_epsg to compute the coordinate transformation. By default, None.
        **kwargs
            Additional keyword argumentts passed into pandas.read_csv.

        Returns
        -------
        Layout
            Layout object created from the text file.
        """
        data = pd.read_csv(filepath, **kwargs)
        geom = MultiPoint(list(zip(data[x_header], data[y_header])))  # type: ignore
        
        if (source_epsg is not None) and (target_epsg is not None):
            source_crs = pyproj.CRS("EPSG:" + str(source_epsg))
            target_crs = pyproj.CRS("EPSG:" + str(target_epsg))
            crs_transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
            geom = transform(crs_transformer.transform, geom)
        
        layout = cls()
        layout._raw_geom = geom
        layout.geom = geom
        return layout

    @classmethod
    def from_postgis(
        cls,
        username: str,
        password: str,
        schema: str,
        table: str,
        host: str = "ow-postgre.postgres.database.azure.com",
        db_name: str = "corp_ta_ea",
        geom_col: str = "geom",
        target_epsg: int | None = None,
        **kwargs,
    ) -> "Layout":
        """from_postgis Alternate constructor to read a PostGIS table into a Layout object. PostGIS table must contain Point or MultiPoint geometry to create a valid Layout. Can reproject geometries on the fly from the source CRS into the desired CRS defined by its EPSG code in the optional target_epsg argument.

        Parameters
        ----------
        username : str
            Username used to create the database connection.
        password : str
            Password used to create the database connection.
        schema : str
            Schema name where the target table is located within the database.
        table : str
            Table name used to create the Layout.
        host : str, optional
            Host url for the target database, by default "ow-postgre.postgres.database.azure.com"
        db_name : str, optional
            Database name, by default "corp_ta_ea"
        geom_col : str, optional
            Name of the table column name containing the geometry data, by default "geom"
        target_epsg : int | None, optional
            EPSG code of the target projection for the geometries to be loaded in. By default, None.
        **kwargs
            Additional keyword arguments passed into geopandas.GeoDataFrame.from_postgis.
        
        Returns
        -------
        Layout
            Returns a Layout object crated from the PostGIS table Point or MultiPoint geometries.

        Raises
        ------
        TypeError
            If the unary union of all the PostGIS geometries does not result in a MultiPoint geometry. I.e. all geometries within the PostGIS table must be Point or MultiPoint.
        """
        # TODO: Should change the geopandas method from_postgis to read_postgis
        db_string = f"postgresql://{username}:{password}@{host}/{db_name}"
        engine = create_engine(db_string)
        data = gp.GeoDataFrame.from_postgis(
            f'SELECT * from "{schema}"."{table}"',
            con=engine,
            geom_col=geom_col,
            **kwargs,
        )
        geom = unary_union(data[geom_col])
        if not isinstance(geom, MultiPoint):
            raise TypeError(
                f"Table geometry type must be Point or MultiPoint to construct a Layout, not {data[geom_col].geom_type.iloc[0]}."
            )
        
        if target_epsg is not None:
            src_crs = data.crs
            target_crs = pyproj.CRS("EPSG:" + str(target_epsg))
            crs_transformer = pyproj.Transformer.from_crs(src_crs, target_crs, always_xy=True)
            geom = transform(crs_transformer.transform, geom)
        
        layout = cls()
        layout._raw_geom = geom
        layout.geom = geom
        return layout

    def is_constrained(self) -> bool:
        """is_constrained Boolean test if the layout is constrained by areas, exclusions or bathymetry.

        Returns
        -------
        bool
            True if the layout is constrained by areas, exclusion or bathymetry. False otherwise.
        """
        return self._constrained

    def has_area(self) -> bool:
        """has_area Test if the layout has an area assigned.

        Returns
        -------
        bool
            Returns True if an area is assigned to the layout. False otherwise.
        """
        return self.area is not None

    def has_geom(self) -> bool:
        """has_geom Test if the layout has any geometry.

        Returns
        -------
        bool
            True if the layout has a geometry. False otherwise.
        """
        return self.geom is not None

    def has_exclusion(self) -> bool:
        """has_exclusion Test if the layout has any exclusions.

        Returns
        -------
        bool
            True if the layout has exclusions. False otherwise.
        """
        return self.exclusion is not None

    def has_bathymetry(self) -> bool:
        """has_bathymetry Test if the layout has a loaded bathymetry.

        Returns
        -------
        bool
            True if the layout has a loaded bathymetry. False otherwise.
        """
        return self.bathymetry_path is not None

    def clip_to_area(self, areas: list[Polygon | MultiPolygon], mode="a") -> "Layout":
        """clip_to_area Clips the layout geometry to the passed area geometries.

        Parameters
        ----------
        areas : list[Polygon  |  MultiPolygon]
            A list of Shapely Polygons or MultiPolygons. The layout geometry will be clipped to
            within these areas.
        mode : str, optional
            Speficies the mode in which the passed areas are added to the layout as constraints.
            Options are 'w' for write mode and 'a' for append mode (case insensitive). When calling
            this method in append mode, the passed geometries will be added to any existing areas as
            additional constraints. When calling this method in write mode, and existing area
            constraints will be overwritten by the passed geometries. By default "a".

        Returns
        -------
        Layout
            Layout with geometries clipped to the passed areas.

        Raises
        ------
        ValueError
            If mode is not 'a' or 'w'.
        ValueError
            If the layout has no geometries to clip.
        """
        if mode.lower() not in ["a", "w"]:
            raise ValueError(f"mode must be 'a' or 'w'. Value of {mode=} passed.")

        if not self.has_geom():
            raise ValueError("Layout has no geometry to clip.")

        area = unary_union(areas)

        if not self._constrained:
            self._constrained = True

        if self.has_area() and mode.lower() == "a":
            self.area = unary_union([self.area, area])
            self.geom = self._raw_geom
            if self.has_exclusion():
                self.avoid_exclusions(self.exclusion)  # type: ignore
        elif self.has_area() and mode.lower() == "w":
            self.reset_area()
            self.area = area
        else:
            self.area = area

        self.geom = self.geom.intersection(self.area)  # type: ignore

        return self

    def avoid_exclusions(
        self, exclusions: list[Polygon | MultiPolygon], mode="a"
    ) -> "Layout":
        """avoid_exclusions Removes all points that overlap with the passed geometries.

        Parameters
        ----------
        exclusions : list[Polygon  |  MultiPolygon]
            A list of Shapely Polygon or MultiPolygon geometries. Any layout geometries within these
            areas will be removed.
        mode : str, optional
            Speficies the mode in which the passed exclusions are added to the layout as constraints.
            Options are 'w' for write mode and 'a' for append mode (case insensitive). When calling
            this method in append mode, the passed geometries will be added to any existing
            exclusions as additional constraints. When calling this method in write mode, and
            existing exclusion constraints will be overwritten by the passed geometries.
            By default "a"

        Returns
        -------
        Layout
            Layout with geometries overlapping the passed exclusions removed.

        Raises
        ------
        ValueError
            If mode is not 'a' or 'w'.
        ValueError
            If the layout has no geometries to clip.
        """
        if mode.lower() not in ["a", "w"]:
            raise ValueError(f"mode must be 'a' or 'w'. Value of {mode=} passed.")

        if not self.has_geom():
            raise ValueError("Layout has no geometry to clip.")

        if not self._constrained:
            self._constrained = True
        exclusion = unary_union(exclusions)

        if self.has_exclusion() and mode.lower() == "a":
            self.exclusion = unary_union([self.exclusion, exclusion])
        elif self.has_exclusion() and mode.lower() == "w":
            self.reset_exclusion()
            self.exclusion = exclusion
        else:
            self.exclusion = exclusion

        self.geom = self.geom.difference(self.exclusion)  # type: ignore

        return self

    def load_bathymetry(self,
                        filepath: str,
                        sign: str = '-',
                        band: int = 1,
                        limits: tuple[float, float] = (0., 60.),
                        drop_na: bool = False
    ) -> "Layout":
        """load_bathymetry Loads bathymetry data, stores the path, limits, and sign, and removes any invalid positions.

        Parameters
        ----------
        filepath : str
            filepath of the bathymetry raster to load. Must be a format supported by the rasterio.open function.
        sign : str, optional
            sign in which the bathymetry file represents water depths. By default '-' if bathymetry values are negative. Pass '+' if bathymetry values are positive.
        band : int, optional
            Band number of the raster file containing bathymetry data, by default 1. Following the GDAL convention, bands are indexed from 1.
        limits : tuple[float, float], optional
            Minimum and maximum bathymetry values, in the raster file units, for layout positions to be considered valid, by default (0., 60.).
        drop_na : bool, optional
            If True, consider layout positions with nan values as invalid. By default False.

        Returns
        -------
        Layout
            Layout with invalid bathymetry positions removed.

        Raises
        ------
        ValueError
            If sign is not '-' or '+'.
        
        See Also
        --------
        rasterio.open : Rasterio function for opening datasets.
        """
        if sign not in ["-", "+"]:
            raise ValueError(f"sign must be '-' or '+'. Value of {sign=} passed.")
        
        if not self._constrained:
            self._constrained = True
        
        SIGN = {'-': -1, '+': 1}[sign]
        
        with rasterio.open(filepath, mode='r') as src:
            samples = [sample[0]*SIGN for sample in src.sample(self.coords, indexes=band)]
            
        valid_points = [point for i, point in enumerate(self.geom.geoms) if _valid_sample(samples[i], limits, drop_na)]
        
        self.geom = unary_union(valid_points)
        self.bathymetry_path = filepath
        self.bathymetry_sign = sign
        self.bathymetry_limits = limits
        self.bathymetry_drop_na = drop_na
            
        return self

    def apply_bathymetry(self,
                        dataset: rasterio.io.DatasetReader,
                        sign: str = '-',
                        band: int = 1,
                        limits: tuple[float, float] = (0., 60.),
                        drop_na: bool = False) -> "Layout":
        """apply_bathymetry Alternative method to apply minimum and maximum bathymetry limits to a Layout. This method takes in an open rasterio.IO.DatasetReader object. More efficient than Layout.load_bathymetry if applied within loops as this method avoids opening and closing the raster data file on every iteration.

        Parameters
        ----------
        dataset : rasterio.io.DatasetReader
            A rasterio.io.DatasetReader object in an open state.
        sign : str, optional
            sign in which the bathymetry file represents water depths. By default '-' if bathymetry values are negative. Pass '+' if bathymetry values are positive.
        band : int, optional
            Band number of the raster file containing bathymetry data, by default 1. Following the GDAL convention, bands are indexed from 1.
        limits : tuple[float, float], optional
            Minimum and maximum bathymetry values, in the raster file units, for layout positions to be considered valid, by default (0., 60.).
        drop_na : bool, optional
            If True, consider layout positions with nan values as invalid. By default False.

        Returns
        -------
        Layout
            Layout with invalid bathymetry positions removed.

        Raises
        ------
        IOError
            If dataset is in a closed state.
        ValueError
            If sign is not '-' or '+'.
        
        See Also
        --------
        rasterio.open : Rasterio function for opening datasets.
        rasterio.io.DatasetReader : A rasterio unbuffered data and metadata reader.
        
        Examples
        --------
        This examples shows how to open a raster file using rasterio, create a layout, and apply bathymetry constraints to the layout.
        
        >>> import rasterio
        >>> from layowt.grids import Grid
        >>> from layowt.layouts import Layout, geoms_from_shapefile
        >>> areas = geoms_from_shapefile('area.shp')
        >>> exclusions = geoms_from_shapefile('exclusions.shp')
        >>> grid = Grid(n_rows=40, n_cols=40, origin=(areas[0].centroid.x, areas[0].centroid.y), scale=236)
        >>> layout = Layout(grid, areas=areas, exclusions=exclusions)
        >>> with rasterio.open("bathymetry.tif") as dataset:
            >>> layout.apply_bathymetry(dataset, limits=(0, 56))
        >>> layout.plot(ax=ax, show_bathy=True)
        
        .. image:: ../../../../_static/Layout_apply_bathymetry_example1.png
        """
        if dataset.closed:
            raise IOError("Bathymetry dataset is closed.")
        
        if sign not in ["-", "+"]:
            raise ValueError(f"sign must be '-' or '+'. Value of {sign=} passed.")
        
        if not self._constrained:
            self._constrained = True
        
        SIGN = {'-': -1, '+': 1}[sign]
        
        samples = [sample[0]*SIGN for sample in dataset.sample(self.coords, indexes=band)]
            
        valid_points = [point for i, point in enumerate(self.geom.geoms) if _valid_sample(samples[i], limits, drop_na)]
        
        self.geom = unary_union(valid_points)
        self.bathymetry_sign = sign
        self.bathymetry_limits = limits
        self.bathymetry_drop_na = drop_na
        
        return self

    def reset_area(self) -> "Layout":
        """reset_area Removes any area constraints from the Layout object. Recovers any geometries that were lost due to the area constraint. Maintains any existing exclusion and bathymetry constraints.

        Returns
        -------
        Layout
            Layout without any area constraints.
        """
        if self.has_area():
            self.area = None
            if self.has_geom():
                self.geom = self._raw_geom
                if self.has_exclusion():
                    self.avoid_exclusions(self.exclusion)  # type: ignore
                if self.has_bathymetry():
                    self.load_bathymetry(
                        filepath=self.bathymetry_path,
                        limits=self.bathymetry_limits,
                        sign=self.bathymetry_sign,
                        drop_na=self.bathymetry_drop_na)
        return self

    def reset_exclusion(self) -> "Layout":
        """reset_exclusion Removes any exclusion constraints from the Layout object. Recovers any geometries that were lost due to the exclusion constraint. Maintaints any existing area and bathymetry constraints.

        Returns
        -------
        Layout
            Layout without any exclusion constraints.
        """
        if self.has_exclusion():
            self.exclusion = None
            self._constrained = False
            if self.has_geom():
                self.geom = self._raw_geom
                if self.has_area():
                    self.clip_to_area(self.area)  # type: ignore
                    self._constrained = True
                if self.has_bathymetry():
                    self.load_bathymetry(
                        filepath=self.bathymetry_path,
                        limits=self.bathymetry_limits,
                        sign=self.bathymetry_sign,
                        drop_na=self.bathymetry_drop_na)
                    self._constrained = True
        return self

    def reset_bathymetry(self) -> "Layout":
        """reset_bathymetry Removes any bathymetry constraints from the Layout object. Recovers any geometries that were lost due to the bathymetry constraint. Maintains any existing area and exclusion constraints.

        Returns
        -------
        Layout
            Layout without any bathymetry constraints.
        """
        if self.has_bathymetry():
            self.bathymetry_path = None
            self.bathymetry_limits = None
            self.bathymetry_sign = None
            self.bathymetry_drop_na = None
            self._constrained = False
        if self.has_geom():
            self.geom = self._raw_geom
            if self.has_area():
                self.clip_to_area(self.area)
                self._constrained = True
            if self.has_exclusion():
                self.avoid_exclusions(self.exclusion)
                self._constrained = True
        return self

    def reset_geom(self) -> "Layout":
        """reset_geom Resets the layout to its original constructor state, without any area or exclusion constraints on its original geometries.

        Returns
        -------
        Layout
            Layout with its original geometries and no area or exclusion constraints.
        """
        if self.has_geom():
            self.geom = self._raw_geom
            self._constrained = False
            self.area = None
            self.exclusion = None
            self.bathymetry_path = None
            self.bathymetry_limits = None
            self.bathymetry_sign = None
        return self

    @property
    def coords(self) -> np.ndarray:
        """coords x, y coordinates of the Point geometries of the Layout object. If

        Returns
        -------
        np.ndarray
            Array of [x, y] coordinate pairs of the points in the Layout object.
        """
        coord_array = []
        for geom in self.geom.geoms:  # type: ignore
            coord_array.append([geom.x, geom.y])
        return np.array(coord_array)

    @property
    def x(self) -> np.ndarray:
        """x 1-D array of the x coordinates of the points in the Layout.

        Returns
        -------
        np.ndarray
            1-D array of the x coordinates of the points in the Layout.
        """
        return self.coords[:, 0]

    @property
    def y(self) -> np.ndarray:
        """y 1-D array of the y coordinates of the points in the Layout.

        Returns
        -------
        np.ndarray
            1-D array of the y coordinates of the points in the Layout.
        """
        return self.coords[:, 1]

    @property
    def n_wtg(self) -> int:
        """n_wtg Number of valid points in the Layout.

        Returns
        -------
        int
            Number of valid points in the Layout.
        """
        return len(self.geom.geoms)  # type: ignore

    def plot(
        self, ax: Axes | None = None, figsize: tuple = (8, 6), show_bathy: bool = False
    ) -> tuple[Figure, Axes]:
        """plot Plots the Layout geometry, as well as areas and exclusions.

        Parameters
        ----------
        ax : matplotlib.axes.Axes | None, optional
            Axes on which the Layout is drawn. If None, a new Figure and Axes are created, by default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6)
        show_bathy : bool, optional
            If True, render the bathymetry onto the plot. Layout must have a loaded bathymetry filepath.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects

        Raises
        ------
        TypeError
            The passed ax argument was not of type Axes
        
        Caution
        ----
        Performance can be seriously reduced if show_bathy is set to true and high resolution bathymetry data is to be rendered.
        """
        # TODO: Update this method to use GeoPandas for easier plotting.
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
        
        if show_bathy and self.bathymetry_path is not None:
            with rasterio.open(self.bathymetry_path, mode='r') as src:
                rasterio.plot.show(src, ax=ax, alpha=0.5)
                plt.colorbar(mappable=ax.get_images()[0], ax=ax)
        
        # Plots the points of the valid turbines
        ax.scatter(self.x, self.y, color="g")

        # Plots the project area(s) and their interior holes if present
        if self.area is not None:
            if isinstance(self.area, MultiPolygon):
                for geom in self.area.geoms:  # type: ignore
                    ax.plot(*geom.exterior.xy, "r")  # type: ignore
                    if len(geom.interiors) != 0:  # type: ignore
                        for hole in geom.interiors:  # type: ignore
                            ax.plot(*hole.xy, "r")  # type: ignore
            if isinstance(self.area, Polygon):
                ax.plot(*self.area.exterior.xy, "r")  # type: ignore
                if len(self.area.interiors) != 0:  # type: ignore
                    for hole in self.area.interiors:  # type: ignore
                        ax.plot(*hole.xy, "r")

        # Plots the project exclusions(s) and their interior holes if present
        if self.exclusion is not None:
            if isinstance(self.exclusion, MultiPolygon):
                for geom in self.exclusion.geoms:  # type: ignore
                    ax.fill(*geom.exterior.xy, color="black", alpha=0.2)
                    if len(geom.interiors) != 0:
                        for hole in geom.interiors:
                            ax.fill(*hole.xy, color="green", hatch="/", alpha=0.2)
            if isinstance(self.exclusion, Polygon):
                ax.fill(*self.exclusion.exterior.xy, color="black", alpha=0.2)  # type: ignore
                if len(self.exclusion.interiors) != 0:  # type: ignore
                    for hole in self.exclusion.interiors:  # type: ignore
                        ax.fill(*hole.xy, color="green", hatch="/", alpha=0.2)

        ax.set_aspect("equal", "box")

        return fig, ax

    def load_wtg(self, file: str) -> None:
        """load_wtg Loads a WAsP .wtg file into the layout as a py_wake WindTurbine object.

        Parameters
        ----------
        file : str
            Filepath for the .wtg file to load.

        Returns
        -------
        None
        
        See Also
        --------
        py_wake.wind_turbines.WindTurbines : Set of multiple py_wake wind turbines.
        """
        self.wtg = WindTurbine.from_WAsP_wtg(file)
        return None

    def to_shapefile(self, filepath: str, epsg: int) -> None:
        """to_shapefile Saves the Layout Point geometries as a Shapefile.

        Parameters
        ----------
        filepath : str
            Filepath to save the Shapefile as.
        epsg : int
            Valid EPSG code for the projection of the Shapefile.
        """
        # TODO: Create a more detailed attribute schema (wakes, water depth, etc.).
        schema = {"geometry": "Point", "properties": {"id": "int"}}

        driver = "ESRI Shapefile"
        crs = pyproj.CRS.from_epsg(epsg).to_dict()

        with fiona.open(
            filepath, "w", crs=crs, driver=driver, schema=schema
        ) as new_file:
            for idx, geom in enumerate(self.geom.geoms, start=1):  # type: ignore
                new_file.write({"geometry": mapping(geom), "properties": {"id": idx}})

    def to_postgis(
        self,
        username: str,
        password: str,
        schema: str,
        table: str,
        host: str = "ow-postgre.postgres.database.azure.com",
        db_name: str = "corp_ta_ea",
        **kwargs,
    ):
        """to_postgis Not yet implemented.

        Parameters
        ----------
        username : str
            _description_
        password : str
            _description_
        schema : str
            _description_
        table : str
            _description_
        host : str, optional
            _description_, by default "ow-postgre.postgres.database.azure.com"
        db_name : str, optional
            _description_, by default "corp_ta_ea"

        Raises
        ------
        NotImplementedError
            Method not yet implemented.
        """
        # TODO write method to save layout as a PostGIS table.
        raise NotImplementedError("Method not yet implemented.")


def _valid_sample(sample: float, limits: tuple[float, float], drop_na: bool) -> bool:
    """_valid_sample helper function used to validate bathymetry sample points.

    Parameters
    ----------
    sample : float
        Sampled value from the bathymetry raster.
    limits : tuple[float, float]
        tuple of (min, max) acceptable bathymetry values as positive integers.
    drop_na : bool
        If True, regards nan samples as invalid. If False, keeps positions with non samples.

    Returns
    -------
    bool
        If a bathymetry sample is valid or not.
    """
    if drop_na:
        return (sample >= limits[0]) and (sample <= limits[1])
    else:
        return ((sample >= limits[0]) and (sample <= limits[1]) or sample is np.nan)
