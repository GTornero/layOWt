"""This module contains the Layout class which is the basic data structure used
to define project scenarios."""

import copy

import fiona
import geopandas as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
from layowt.grids import Grid
from matplotlib.figure import Figure
from py_wake.wind_turbines import WindTurbine
from shapely.geometry import MultiPoint, MultiPolygon, Point, Polygon, mapping, shape
from shapely.ops import unary_union
from sqlalchemy import create_engine


class Layout:
    """Layout class.

    This class is responsible for holding all of the data structures required
    in order to define a project scenario and to generate the coordinate pairs
    of the points that make up the grid.

    Multiple numerical, spatial and graphical methods are available to
    manipulate the scenario o conform to constraint and set properties.

    """

    def __init__(
        self,
        grid: Grid | None = None,
        areas: list[Polygon | MultiPolygon] | None = None,
        exclusions: list[Polygon | MultiPolygon] | None = None,
        wtg: str | None = None,
    ) -> None:
        """__init__ Initialises the Layout class instance.

        Parameters
        ----------
        grid : Grid | None, optional
            Grid object used to define the Layout, by default None
        areas : list[Polygon  |  MultiPolygon] | None, optional
            List of Shapely Polygons or MultiPolygons used to define the lease
            areas, by default None
        exclusions : list[Polygon  |  MultiPolygon] | None, optional
            List of Shapely Polygons or MultiPolygons used to define exclusion
            zones where turbines cannot be placed, by default None
        """
        self.grid = grid

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
    def from_shapefile(cls, filepath: str) -> "Layout":
        """from_shapefile Alternate constructor to read a shapefile into a Layout object.
        The shapefile must contain Point or 3DMultiPoint geometry to be a valid Layout constructor.

        Parameters
        ----------
        filepath : str
            Valid string file path to the shapefile to read.

        Returns
        -------
        Layout
            Layout object created from the shapefile.
        """
        with fiona.open(filepath) as src:
            geom = [shape(rec["geometry"]) for rec in src]

        if all(isinstance(x, (Point, MultiPoint)) for x in geom):
            geom = unary_union(geom)
            layout = cls()
            layout._raw_geom = geom
            layout.geom = geom
            return layout
        else:
            raise TypeError("Geometry type must be Point or MultiPoint.")

    @classmethod
    def from_text(
        cls, filepath: str, x_header: str | int, y_header: str | int, **kwargs
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

        Returns
        -------
        Layout
            Layout object created from the text file.
        """
        data = pd.read_csv(filepath, **kwargs)
        geom = MultiPoint(list(zip(data[x_header], data[y_header])))  # type: ignore
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
        geom_col="geom",
        **kwargs,
    ) -> "Layout":
        """from_postgis Alternate constructor to read a PostGIS table into a Layout object. PostGIS table must contain Point or MultiPoint geometry to create a valid Layout.

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

        Returns
        -------
        Layout
            Returns a Layout object crated from the PostGIS table Point or MultiPoint geometries.

        Raises
        ------
        TypeError
            If the unary union of all the PostGIS geometries does not result in a MultiPoint geometry. I.e. all geometries within the PostGIS table must be Point or MultiPoint.
        """
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
        layout = cls()
        layout._raw_geom = geom
        layout.geom = geom
        return layout

    def is_constrained(self) -> bool:
        """is_constrained Boolean test if the layout is constrained by areas or exclusions.

        Returns
        -------
        bool
            True if the layout is constrained by areas of exclusion geometries. False otherwise.
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

    def reset_area(self) -> "Layout":
        """reset_area Removes any area constraints from the Layout object. Recovers any geometries that were lost due to the area constraint. Maintains any existing exclusion constraints.

        Returns
        -------
        Layout
            Layout without any area constraints.
        """
        if self.has_area():
            self.area = None
            if self.has_geom():
                if self.has_exclusion():
                    self.geom = self._raw_geom.difference(self.exclusion)  # type: ignore
                else:
                    self.geom = self._raw_geom
        return self

    def reset_exclusion(self) -> "Layout":
        """reset_exclusion Removes any exclusion constraints from the Layout object. Recovers any geometries that were lost due to the exclusion constraint. Maintaints any existing area constraints.

        Returns
        -------
        Layout
            Layout without any exclusion constraints.
        """
        if self.has_exclusion():
            self.exclusion = None
            if self.has_geom():
                if self.has_area():
                    self.geom = self._raw_geom.intersection(self.area)  # type: ignore
                else:
                    self.geom = self._raw_geom
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
        self, ax: plt.Axes | None = None, figsize: tuple = (8, 6)
    ) -> tuple[Figure, plt.Axes]:
        """plot Plots the Layout geometry, as well as areas and exclusions.

        Parameters
        ----------
        ax : plt.Axes | None, optional
            Axes on which the Layout is drawn. If None, a new Figure and Axes are created, by default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6)

        Returns
        -------
        tuple[Figure, plt.Axes]
            Tuple of Figure and Axes objects

        Raises
        ------
        TypeError
            The passed ax argument was not of type plt.Axes
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
