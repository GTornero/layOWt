"""This module contains the Layout class which is the basic data structure used
to define project scenario geometries."""

import copy
import urllib.parse
from typing import Iterable, Optional

import contextily as ctx
import fiona
import geopandas as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import rasterio
import sqlalchemy as sa
import xyzservices
from geoalchemy2.shape import from_shape, to_shape
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from rasterio.errors import CRSError
from rasterio.plot import show
from scipy.spatial.distance import cdist
from shapely.geometry import (MultiPoint, MultiPolygon, Point, Polygon,
                              mapping, shape)
from shapely.ops import transform, unary_union

from layowt.errors import LayoutSchemaError
from layowt.grids import Grid


class Layout:
    """
    This class is responsible for holding all of the data structures required \
    in order to define a project scenario and to generate the coordinate pairs \
    of the points that make up the grid. \
    \
    Multiple numerical, spatial and graphical methods are available to \
    manipulate the scenario o conform to constraint and set properties. \
    
    Parameters
    ----------
    grid : Grid | None, optional
        Grid object used to define the Layout, by default None.
    epsg : int | None, optional
        EPSG code for the layout coordinates, by default None.
    areas : list[Polygon  |  MultiPolygon] | None, optional
        List of Shapely Polygons or MultiPolygons used to define the lease
        areas, by default None.
    exclusions : list[Polygon  |  MultiPolygon] | None, optional
        List of Shapely Polygons or MultiPolygons used to define exclusion
        zones where turbines cannot be placed, by default None.
    
    See Also
    --------
    Polygon : Shapely Polygon object.
    MultiPolygon : Shapely MultiPolygon object.
    
    Attributes
    ----------
    grid : Grid
        Attribute: Underlying grid to the layout. None for layouts not created from a grid.
    epsg : int
        Attribute: EPSG code of the layout coordinate system.
    area : list[Polygon | MultiPolygon]
        Attrbute: List of buildable area geometries of the layout.
    exclusion : list[Polygon | MultiPolygon]
        Attribute: List of exclusion area geometries of the layout.
    bathymetry_path : str
        Attribute: filepath of the layout's bathymetry dataset.
    bathymetry_limits: tuple(float, float)
        Attribute: Tuple of minimum and maximum allowable bathymetry values for valid layout positions.
    bathymetry_sign : str
        Attribute: Sign ("+" or "-") in which bathymetry values are represented.
    bathymetry_drop_na : bool
        Attributes: If layout positions at locations of missing bathymetry have been removed.
    bathymetry_epsg : int
        Attribute: EPSG code of the bathymetry raster file.
    bathymetry_band : int
        Attribute: Band number of the bathymetry raster file being used.
    samples : list[float]
        Attribute: List of sampled bathymetry values at the layout positions.
    geom : MultiPoint
        Attribute: MultiPoint geometry representation of the layout positions.
    """

    def __init__(
        self,
        grid: Grid | None = None,
        epsg: int | None = None,
        areas: list[Polygon | MultiPolygon] | None = None,
        exclusions: list[Polygon | MultiPolygon] | None = None
    ) -> None:
        """__init__ Initialises the Layout class instance.
        """
        self.grid = grid
        if _is_epsg(epsg) or epsg is None:
            self.epsg = epsg
        else:
            raise CRSError(f"{epsg=} not a valid EPSG code.")
        self.bathymetry_path = None
        self.bathymetry_limits = None
        self.bathymetry_sign = None
        self.bathymetry_drop_na = None
        self.bathymetry_epsg = None
        self.bathymetry_band = None
        self.samples = None

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

    @classmethod
    def from_shapefile(cls, filepath: str, target_epsg: int | None = None) -> "Layout":
        """from_shapefile Alternate constructor to read a shapefile into a Layout object.
        The shapefile must contain Point or MultiPoint geometry to be a valid Layout constructor. Can reproject geometries on the fly from the source CRS into the desired CRS defined by its EPSG code in the optional target_epsg argument.

        Parameters
        ----------
        filepath : str
            Valid string file path to the shapefile to read.
        target_epsg : int | None, optional
            `EPSG <https://epsg.io/>`_ code of the target projection for the geometries to be loaded in. By default, None.

        Returns
        -------
        Layout
            Layout object created from the shapefile.
        
        Raises
        ------
        CRSError
            If target_epsg is not a valid EPSG code.
        """
        if not _is_epsg(target_epsg) and target_epsg is not None:
            raise CRSError(f"{target_epsg=} not a valid EPSG code.")
        
        with fiona.open(filepath) as src:
            geoms = [shape(rec["geometry"]) for rec in src]
            src_crs = pyproj.CRS(src.crs_wkt)
            layout_epsg = src_crs.to_epsg()
        
        if target_epsg is not None:
            layout_epsg = target_epsg
            target_crs = pyproj.CRS.from_epsg(target_epsg)
            crs_transformer = pyproj.Transformer.from_crs(src_crs, target_crs, always_xy=True)
            transformed_geoms = []
            for geom in geoms:
                transformed_geom = transform(crs_transformer.transform, geom)
                transformed_geoms.append(transformed_geom)
            
            geoms = transformed_geoms

        if all(isinstance(x, (Point, MultiPoint)) for x in geoms):
            geoms = unary_union(geoms)
            layout = cls(epsg=layout_epsg)
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
            `EPSG <https://epsg.io/>`_ code of the source projection for the points to be loaded in. Used alongside target_epsg to compute the coordinate transformation. By default, None.
        target_epsg : int | None, optional
            `EPSG <https://epsg.io/>`_ code of the target projection for the points to be loaded in. Used alongside source_epsg to compute the coordinate transformation. By default, None.
        **kwargs
            Additional keyword argumentts passed into pandas.read_csv.

        Returns
        -------
        Layout
            Layout object created from the text file.
        
        Raises
        ------
        CRSError
            If source_epsg or target_epsg are not valid EPSG codes.
        """
        if not _is_epsg(source_epsg) and source_epsg is not None:
            raise CRSError(f"{source_epsg=} not a valid EPSG code.")
        elif not _is_epsg(target_epsg) and target_epsg is not None:
            raise CRSError(f"{target_epsg=} not a valid EPSG code.")
        
        data = pd.read_csv(filepath, **kwargs)
        geom = MultiPoint(list(zip(data[x_header], data[y_header])))  # type: ignore
        layout_epsg = None
        
        if (source_epsg is not None) and (target_epsg is not None):
            layout_epsg = target_epsg
            source_crs = pyproj.CRS("EPSG:" + str(source_epsg))
            target_crs = pyproj.CRS("EPSG:" + str(target_epsg))
            crs_transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
            geom = transform(crs_transformer.transform, geom)
        elif source_epsg is not None:
            layout_epsg = source_epsg
        
        layout = cls(epsg=layout_epsg)
        layout._raw_geom = geom
        layout.geom = geom
        return layout
    
    @classmethod
    def from_xy(cls, x: Iterable[float], y: Iterable[float], source_epsg: int | None = None, target_epsg: int | None = None) -> "Layout":
        """
        Alternate Layout constructor from xa nd y coordinate data. Pass x and y arrays directly from memory.

        Parameters
        ----------
        x : Iterable[float]
            Array of x-coordinates.
        y : Iterable[float]
            Array of y-coordinates.
        source_epsg : int | None, optional
            `EPSG <https://epsg.io/>`_ code of the source projection for the passed coordinates.
            Used alongside target_epsg to compute the coordinate transformation. By default, None.
        target_epsg : int | None, optional
            `EPSG <https://epsg.io/>`_ code of the target projection for the points to be Layout
            to be created in. Used alongside source_epsg to compute the coordinate transformation. By default, None.

        Returns
        -------
        Layout
            Layout object created from the x, y arrays.
        
        Raises
        ------
        ValueError
            If x and y arrays are of different length.
        CRSError
            If source_epsg or target_epsg are not valid EPSG codes.
        """
        if len(x) != len(y):
            raise ValueError(f"Length Mistmatch. x ({len(x)}) and y ({len(y)}) are not the same length.")
        
        if not _is_epsg(source_epsg) and source_epsg is not None:
            raise CRSError(f"{source_epsg=} not a valid EPSG code.")
        elif not _is_epsg(target_epsg) and target_epsg is not None:
            raise CRSError(f"{target_epsg=} not a valid EPSG code.")
        
        geom = MultiPoint(list(zip(x, y)))  # type: ignore
        layout_epsg = None
        
        if (source_epsg is not None) and (target_epsg is not None):
            layout_epsg = target_epsg
            source_crs = pyproj.CRS("EPSG:" + str(source_epsg))
            target_crs = pyproj.CRS("EPSG:" + str(target_epsg))
            crs_transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
            geom = transform(crs_transformer.transform, geom)
        elif source_epsg is not None:
            layout_epsg = source_epsg
        
        layout = cls(epsg=layout_epsg)
        layout._raw_geom = geom
        layout.geom = geom
        return layout

    @classmethod
    def from_postgis(
        cls,
        engine: sa.engine.Engine,
        schema: str,
        table: str,
        geom_col: str = "geom",
        target_epsg: int | None = None,
        **kwargs,
    ) -> "Layout":
        """from_postgis Alternate constructor to read a PostGIS table into a Layout object. PostGIS table must contain Point or MultiPoint
        geometry to create a valid Layout. Can reproject geometries on the fly from the source CRS into the desired CRS defined by its EPSG
        code in the optional target_epsg argument.

        Parameters
        ----------
        engine: sqlalchemy.engine.Engine
            SQLAlchemy database engine.
        schema : str
            Schema name where the target table is located within the database.
        table : str
            Table name used to create the Layout.
        geom_col : str, optional
            Name of the table column name containing the geometry data, by default "geom"
        target_epsg : int | None, optional
            `EPSG <https://epsg.io/>`_ code of the target projection for the geometries to be loaded in. By default, None.
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
        rasterio.errors.CRSError
            If target_epsg is not a valid EPSG code.
        """
        # TODO: Should change the geopandas method GeoDataFrame.from_postgis to geopandas.read_postgis
        new_engine = sa.create_engine(engine.url, echo=False, future=False)
        
        data = gp.read_postgis(
            f'SELECT * from "{schema}"."{table}"',
            con=new_engine,
            geom_col=geom_col,
            **kwargs,
        )
        geom = unary_union(data[geom_col])
        layout_epsg = data.crs.to_epsg()
        if not isinstance(geom, MultiPoint):
            raise TypeError(
                f"Table geometry type must be Point or MultiPoint to construct a Layout, not {data[geom_col].geom_type.iloc[0]}."
            )
        
        if target_epsg is not None:
            if not _is_epsg(target_epsg):
                raise CRSError(f"{target_epsg=} not a valid EPSG code.")
            layout_epsg = target_epsg
            src_crs = data.crs
            target_crs = pyproj.CRS.from_epsg(target_epsg)
            crs_transformer = pyproj.Transformer.from_crs(src_crs, target_crs, always_xy=True)
            geom = transform(crs_transformer.transform, geom)
        
        layout = cls(epsg=layout_epsg)
        layout._raw_geom = geom
        layout.geom = geom
        return layout
    
    @classmethod
    def from_database(
        cls,
        engine: sa.engine.Engine,
        schema: str,
        layout_id: int
    ) -> "Layout":
        """
        Alternate constructor that creates a Layout from the OW PostGIS Energy Assessment database \
        by specifying the schema and layout_id.
        
        The Layout instance will be initiated with the turbine, grid, area, and exclusion
        geometries from the database. The target schema must have been set up using the 
        :func:`layowt.layouts.utils.setup_database_schema` function.

        Parameters
        ----------
        engine : sqlalchemy.engine.Engine
            SQLAlchemy database engine.
        schema : str
            Name of the schema where the target layout is located within the database.
        layout_id : int
            layout_id of the Layout within the "ow__layowt_layout_info" table.

        Returns
        -------
        Layout
            Layout loaded from the database.
        
        Raises
        ------
        KeyError
            If layout_id is not found in the database.
        
        See Also
        --------
        layowt.layouts.utils.setup_database_schema : Function to setup a database schema with the necessary tables for layout upload.
        layowt.layouts.layout.Layout.database_engine : Static method to create a database engine.
        layowt.layouts.utils.database_engine : Function to create a database engine.
        """
        metadata = sa.MetaData()
        
        info_table = sa.Table("ow__layowt_layout_info", metadata, autoload_with=engine, schema=schema)
        turbine_table = sa.Table("ow__layowt_turbines", metadata, autoload_with=engine, schema=schema)
        boundary_table = sa.Table("ow__layowt_boundaries", metadata, autoload_with=engine, schema=schema)

        select_layout = sa.select(info_table).where(info_table.c.layout_id==layout_id)
        select_turbines = sa.select(sa.func.ST_Union(turbine_table.c.geom)).where(turbine_table.c.layout_id==layout_id)


        layout_epsg = _get_epsg_from_postgis(engine=engine, schema=schema, table="ow__layowt_turbines", geom_col="geom") 

        with engine.connect() as conn:
            layout = cls(epsg=layout_epsg) # change me to cls for the class method
            layout._constrained = True
            
            # 1 read the layout_info column
            layout_info = conn.execute(select_layout).fetchone()
            if layout_info is None:
                raise KeyError(f"Layout with layout_id {layout_id} not found.")
            # 2 read the turbine geometries
            layout.geom = to_shape(conn.execute(select_turbines).fetchone()[0])
            # 4 read the boundary geometry
            select_boundary = sa.select(boundary_table.c.geom).where(boundary_table.c.boundary_id==layout_info["boundary_id"])
            layout.area = to_shape(conn.execute(select_boundary).fetchone()[0])
            # 3 check the grid and load the parameters
            if layout_info["grid_id"] is not None:
                grid_table = sa.Table("ow__layowt_grids", metadata, autoload_with=engine, schema=schema)
                select_grid = sa.select(grid_table).where(grid_table.c.grid_id==layout_info["grid_id"])
                grid_info = conn.execute(select_grid).fetchone()
                layout.grid = Grid(n_rows=grid_info["n_rows"],
                                    n_cols=grid_info["n_cols"],
                                    row_step=float(grid_info["row_step"]),
                                    col_step=float(grid_info["col_step"]),
                                    row_offset=grid_info["row_offset"],
                                    col_offset=grid_info["col_offset"],
                                    angle=float(grid_info["angle"]),
                                    x_shear=float(grid_info["x_shear"]),
                                    y_shear=float(grid_info["y_shear"]),
                                    origin=(grid_info["origin_x"], grid_info["origin_y"]),
                                    scale=grid_info["scale"]
                                    )
                # 3.1 assign the raw_geom as the grid geom
                layout._raw_geom = to_shape(grid_info["geom"])
            else:
                # 3.2 if the layout has no grid, turbines geom is the raw_geom
                layout._raw_geom = layout.geom
            # 5 check if the layout has turbine constraints and read the geometry
            if layout_info["turbine_constraint_id"] is not None:
                constraint_table = sa.Table("ow__layowt_constraint_scenarios", metadata, autoload_with=engine, schema=schema)
                select_constraint = sa.select(constraint_table.c.geom).where(constraint_table.c.constraint_id==layout_info["turbine_constraint_id"])
                layout.exclusion = to_shape(conn.execute(select_constraint).fetchone()[0])
        
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
        
        if self.has_bathymetry():
            self.load_bathymetry(filepath=self.bathymetry_path,
                                 sign=self.bathymetry_sign,
                                 band=self.bathymetry_band,
                                 limits=self.bathymetry_limits,
                                 drop_na=self.bathymetry_drop_na)

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
        
        if self.has_bathymetry():
            self.load_bathymetry(filepath=self.bathymetry_path,
                                 sign=self.bathymetry_sign,
                                 band=self.bathymetry_band,
                                 limits=self.bathymetry_limits,
                                 drop_na=self.bathymetry_drop_na)

        return self

    def load_bathymetry(self,
                        filepath: str,
                        sign: str = '-',
                        band: int = 1,
                        limits: tuple[float, float] = (0., 60.),
                        drop_na: bool = False
    ) -> "Layout":
        """
        Loads bathymetry data, stores the path, limits, and sign, and removes any invalid positions.
        Bathymetry values are stored in the `samples` attribute of the Layout instance.
        
        If the raster data crs does not match the layout crs (and both the layout and raster data
        have valid epsg codes) the layout coordinates will be reprojected on-the-fly to match the
        raster projection and extract the sampled values. This reprojections does not change the layout instance.

        Parameters
        ----------
        filepath : str
            filepath of the bathymetry raster to load. Must be a format supported by the rasterio.open function.
        sign : str, optional
            sign in which the bathymetry file represents water depths. By default '-' if bathymetry values are negative. Pass '+' if bathymetry values are positive.
        band : int, optional
            Band number of the raster file containing bathymetry data, by default 1. Following the GDAL convention, bands are indexed from 1.
        limits : tuple[float, float], optional
            Minimum and maximum bathymetry values (as positive values), in the raster file units, for layout positions to be considered valid, by default (0., 60.).
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
        CRSError
            If the raster does not have a CRS with a valid ESPG code.
        
        See Also
        --------
        rasterio.open : Rasterio function for opening datasets.
        """
        if sign not in ["-", "+"]:
            raise ValueError(f"sign must be '-' or '+'. Value of {sign=} passed.")
        
        if not self._constrained:
            self._constrained = True
        
        if not self.geom.is_empty:
            SIGN = {'-': -1, '+': 1}[sign]
            
            with rasterio.open(filepath, mode='r') as src:
                if src.crs.is_epsg_code:
                    self.bathymetry_epsg = src.crs.to_epsg()

                if self.epsg != self.bathymetry_epsg and all((self.epsg, self.bathymetry_epsg)) and self.geom is not None:
                    if not src.crs.is_epsg_code:
                        raise CRSError("Raster CRS is not from a valid EPSG code.")
                    source_crs = pyproj.CRS("EPSG:" + str(self.epsg))
                    target_crs = pyproj.CRS("EPSG:" + str(self.bathymetry_epsg))
                    crs_transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
                    transformed_geom = transform(crs_transformer.transform, self.geom)
                    coord_list = []
                    for geom in transformed_geom.geoms:  # type: ignore
                        coord_list.append([geom.x, geom.y])
                else:
                    coord_list = self.coords
                
                samples = [sample[0]*SIGN for sample in src.sample(coord_list, indexes=band)]
            
            if self.geom is not None:    
                valid_points = [point for i, point in enumerate(self.geom.geoms) if _valid_sample(samples[i], limits, drop_na)]
                self.geom = unary_union(valid_points)
                
                self.samples = []
                for sample in samples:
                    if _valid_sample(sample, limits, drop_na):
                        self.samples.append(sample)

        self.bathymetry_path = filepath
        self.bathymetry_sign = sign
        self.bathymetry_limits = limits
        self.bathymetry_drop_na = drop_na
        self.bathymetry_band = band
            
        return self

    def apply_bathymetry(self,
                        dataset: rasterio.DatasetReader,
                        sign: str = '-',
                        band: int = 1,
                        limits: tuple[float, float] = (0., 60.),
                        drop_na: bool = False) -> "Layout":
        """
        Alternative method to apply minimum and maximum bathymetry limits to a Layout. This method takes
        in an open rasterio.IO.DatasetReader object. More efficient than :meth:`Layout.load_bathymetry` if applied within loops as this
        method avoids opening and closing the raster data file on every iteration.
        
        If the raster data crs does not match the layout crs (and both the layout and raster data
        have valid epsg codes) the layout coordinates will be reprojected on-the-fly to match the
        raster projection and extract the sampled values. This reprojections does not change the layout instance.

        Parameters
        ----------
        dataset : rasterio.io.DatasetReader
            An open rasterio.io.DatasetReader.
        sign : str, optional
            sign in which the bathymetry file represents water depths. By default '-' if bathymetry values are negative. Pass '+' if bathymetry values are positive.
        band : int, optional
            Band number of the raster file containing bathymetry data, by default 1. Following the GDAL convention, bands are indexed from 1.
        limits : tuple[float, float], optional
            Minimum and maximum bathymetry values (as positive values), in the raster file units, for layout positions to be considered valid, by default (0., 60.).
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
        CRSError
            If the raster does not have a CRS with a valid EPSG code.
        
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
        
        if not self.geom.is_empty:
            SIGN = {'-': -1, '+': 1}[sign]
            
            if dataset.crs.is_epsg_code:
                self.bathymetry_epsg = dataset.crs.to_epsg()
            
            if self.epsg != self.bathymetry_epsg and all((self.epsg, self.bathymetry_epsg)) and self.geom is not None:
                if not dataset.crs.is_epsg_code:
                    raise CRSError("Raster CRS is not from a valid EPSG code.")
                source_crs = pyproj.CRS("EPSG:" + str(self.epsg))
                target_crs = pyproj.CRS("EPSG:" + str(self.bathymetry_epsg))
                crs_transformer = pyproj.Transformer.from_crs(source_crs, target_crs, always_xy=True)
                transformed_geom = transform(crs_transformer.transform, self.geom)
                coord_list = []
                for geom in transformed_geom.geoms:  # type: ignore
                    coord_list.append([geom.x, geom.y])
            else:
                coord_list = self.coords
            
            samples = [sample[0]*SIGN for sample in dataset.sample(coord_list, indexes=band)]
            
            if self.geom is not None:
                valid_points = [point for i, point in enumerate(self.geom.geoms) if _valid_sample(samples[i], limits, drop_na)]
                self.geom = unary_union(valid_points)
                
                self.samples = []
                for sample in samples:
                    if _valid_sample(sample, limits, drop_na):
                        self.samples.append(sample)
        
        self.bathymetry_path = dataset.files[0]
        self.bathymetry_sign = sign
        self.bathymetry_limits = limits
        self.bathymetry_drop_na = drop_na
        self.bathymetry_band = band
        
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
            self.bathymetry_epsg = None
            self.samples = None
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
    def shortest_distance(self) -> float:
        """
        This property returns the shortest euclidiean distance between turbines in the layout projection units.
        
        If the Layout instance has a grid, consider calling this property from the Grid directly.
        
        Returns
        -------
        float
            Shortest euclidean distance between turbines in the layout projection units.
        
        Tip
        ---
        It is much faster to use the same property of the :class:`layowt.grids.grid.Grid` class
        than to call this directly from the Layout instance.
        
        Examples
        --------
        Calling the shortest_distance property from a Layout.
        
        >>> from layowt.grids import Grid
        >>> from layowt.layouts import Layout
        >>> import time
        >>> grid = Grid(col_step=5, row_step=10, scale=250)
        >>> layout = Layout(grid)
        >>> layout.shortest_distance
        1250.0
        
        Calling the shortest_distance property from the Layout's grid.
        
        >>> layout.grid.shortest_distance
        1250.0
        
        Checking the computation time of calling the shortest_distance property on a Layout
        
        >>> start = time.perf_counter()
        >>> layout.shortest_distance
        >>> end = time.perf_counter()
        >>> end - start
        0.01885549999860814
        
        Checking the computation time of calling the shortest_distance property on a Grid
        
        >>> start = time.perf_counter()
        >>> layout.grid.shortest_distance
        >>> end = time.perf_counter()
        >>> end - start
        0.0006460000004153699
        """
        distances = cdist(self.coords, self.coords)
        return distances[distances > 0].min()

    @property
    def n_wtg(self) -> int:
        """n_wtg Number of valid points in the Layout.

        Returns
        -------
        int
            Number of valid points in the Layout.
        """
        if self.geom.is_empty:
            return 0
        elif isinstance(self.geom, Point):
            return 1
        else:
            return len(self.geom.geoms)  # type: ignore

    def plot(
        self, ax: Axes | None = None, figsize: tuple = (8, 6), show_bathy: bool = False, basemap: bool = False, basemap_provider: str | xyzservices.TileProvider | None = None
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
        basemap : bool, optional
            If True, render a basemap onto the plot. Requires an internet connection.
        basemap_provider : str | xyzservices.TileProvider | None, optional
            The tile source: web tile provider, a valid input for a query of a
            :class:`xyzservices.TileProvider` by a name from ``xyzservices.providers`` or
            path to local file. The web tile provider can be in the form of a
            :class:`xyzservices.TileProvider` object or a URL. The placeholders for the XYZ
            in the URL need to be `{x}`, `{y}`, `{z}`, respectively. For local file paths,
            the file is read with `rasterio` and all bands are loaded into the basemap.
            The default is: Stamen Terrain web tiles.

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.

        Raises
        ------
        TypeError
            The passed ax argument was not of type Axes.
        
        Caution
        ----
        Performance can be seriously reduced if show_bathy is set to true and high resolution bathymetry data is to be rendered.
        
        See Also
        --------
        xyzservices.TileProvider : xyzservices TileProvider object
        contextly.providers : contextly wrapper of xzyservices TileProvider object
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
        
        if show_bathy and self.bathymetry_path is not None:
            with rasterio.open(self.bathymetry_path, mode='r') as src:
                rasterio.plot.show(src, ax=ax, alpha=0.5)
                plt.colorbar(mappable=ax.get_images()[0], ax=ax)
        
        # Plots the points of the valid turbines if there any any points left
        if not self.geom.is_empty:
            ax.scatter(self.x, self.y, color="g")

        # Plots the project area(s) and their interior holes if present
        if self.area is not None:
            gp.GeoSeries(self.area).plot(ax=ax, facecolor='none', edgecolor='r', lw=1.5)

        # Plots the project exclusions(s) and their interior holes if present
        if self.exclusion is not None:
            gp.GeoSeries(self.exclusion).plot(ax=ax, facecolor='black', alpha=0.2)
        
        if basemap:
            if self.epsg is None:
                raise AttributeError("Layout must have a valid EPSG code when plotting with basemaps.")
            else:
                ctx.add_basemap(ax=ax, crs=self.epsg, source=basemap_provider)

        ax.set_aspect("equal", "box")

        return fig, ax

    def to_shapefile(self, filepath: str, epsg: int = None) -> None:
        """to_shapefile Saves the Layout Point geometries as a Shapefile.

        Parameters
        ----------
        filepath : str
            Filepath to save the Shapefile as.
        epsg : int
            Valid `EPSG <https://epsg.io/>`_ code for the projection of the Shapefile. If passed,
            this code takes precedence over the layout epsg code (if defined). By default, None.
        
        Raises
        ------
        ValueError
            No EPSG code available for the shapefile CRS.
        """
        
        if epsg is None and self.epsg is not None:
            crs = pyproj.CRS.from_epsg(self.epsg).to_dict()
        elif epsg is None and self.epsg is None:
            raise ValueError("No EPSG code provided.")
        elif _is_epsg(epsg):
            crs = pyproj.CRS.from_epsg(epsg).to_dict()
            
        # TODO: Create a more detailed attribute schema (wakes, water depth, etc.).
        schema = {"geometry": "Point", "properties": {"id": "int"}}

        driver = "ESRI Shapefile"

        with fiona.open(
            filepath, "w", crs=crs, driver=driver, schema=schema
        ) as new_file:
            for idx, geom in enumerate(self.geom.geoms, start=1):  # type: ignore
                new_file.write({"geometry": mapping(geom), "properties": {"id": idx}})
    
    @staticmethod
    def database_engine(username: str, password: str, host: str = "ow-postgre.postgres.database.azure.com", db_name: str = "corp_ta_ea", dialect: str = "postgresql", echo: bool = False) -> sa.engine.Engine:
        """
        Create a SQLAlchemy database Engine.

        Parameters
        ----------
        username : str
            Username used to create the database connection.
        password : str
            Password used to create the database connection.
        host : str, optional
            Host url for the target database, by default "ow-postgre.postgres.database.azure.com"
        db_name : str, optional
            Database name, by default "corp_ta_ea"
        dialect : str, optional
            Database dialect, by default "postgresql"
        echo : bool, optional
            If True, the Engine will log all statements as well as a repr() of their parameter lists to
            the default log handler, which defaults to sys.stdout for output.

        Returns
        -------
        sqlalchemy.engine.Engine
            SQLAlchemy database engine.
        """
        password = urllib.parse.quote_plus(password)
        
        con_str = f"{dialect}://{username}:{password}@{host}/{db_name}"
        engine = sa.create_engine(con_str, echo=echo, future=True)
        
        return engine
    
    def to_database(
        self,
        engine: sa.engine.Engine,
        schema: str,
        boundary_id: Optional[int] = None,
        exclusion_id: Optional[int] = None,
        grid_id: Optional[int] = None,
        upload_boundary: bool = False,
        boundary_name: Optional[str] = None,
        boundary_type: Optional[str] = None,
        boundary_source: Optional[str] = None,
        boundary_comments: Optional[str] = None,
        upload_exclusion: bool = False,
        exclusion_name: Optional[str] = None,
        exclusion_source: Optional[str] = None,
        exclusion_comments: Optional[str] = None,
        upload_grid: bool = False,
        grid_comments: Optional[str] = None,
        layout_comments: Optional[str] = None,
        bathymetry_source: Optional[str] = None,
        bathymetry_datum: Optional[str] = None,
        bathymetry_units: Optional[str] = None,
    ) -> None:
        """
        Uploads the layout and associated data to the OW PostGIS database.
        
        The target schema must have been previously set-up with the necessary system tables for
        Layout upload. Can use the :func:`layowt.layouts.utils.setup_database_schema` function to
        do so. 
        
        Requires the following table system to exist:
            * ow__layowt_boundaries
            * ow__layowt_constraint_scenarios
            * ow__layowt_grids
            * ow__layowt_layout_info
            * ow__layowt_turbines

        Parameters
        ----------
        engine : sqlalchemy.engine.Engine
            SQLAlchemy database engine.
        schema : str
            Target database schema.
        boundary_id : Optional[int], optional
            boundary_id foreign key if not uploading the Layout's own area geometry as a new entry
            in the database and instead referencing an existing geometry in the database
            "ow__layowt_boundaries" table, by default None.
        exclusion_id : Optional[int], optional
            constraint_id foreign key if not uploading the Layout's own turbine exclusion geometry
            as a new entry in the database and instead referencing an existing geometry in the
            database "ow__layowt_constraint_scenarios" table, by default None.
        grid_id : Optional[int], optional
            grid_id foreign key if not uploading the Layout's own grid geometry as a new entry
            in the database and instead referencing an existing geometry in the database
            "ow__layowt_grids" table, by default None.
        upload_boundary : bool, optional
            True if uploading the Layout's area geometry as a new entry to the database
            "ow__layowt_boundaries" table, by default False.
        boundary_name : Optional[str], optional
            Descriptive name of the area geometry, by default None.
        boundary_type : Optional[str], optional
            Description of the type of area described by the Layout's area geometry, by default
            None.
        boundary_source : Optional[str], optional
            Description of the data source of area described by the Layout's area geometry, by
            default None.
        boundary_comments : Optional[str], optional
            General comments on the Layout's area, by default None.
        upload_exclusion : bool, optional
            True if uploading the Layout's turbine exclusion geometry as a new entry to the database
            "ow__layowt_constraint_scenarios" table, by default False.
        exclusion_name : Optional[str], optional
            Descriptive name of the turbine exclusion geometry, by default None.
        exclusion_source : Optional[str], optional
            Description of the data source of the turbine exclusion geometry, by default None.
        exclusion_comments : Optional[str], optional
            General comments on the Layout's turbine exclusions, by default None.
        upload_grid : bool, optional
            True if uploading the Layout's grid geometry as a new entry to the database
            "ow__layowt_grids" table, by default False.
        grid_comments : Optional[str], optional
            General comments on the Layout's grid, by default None.
        layout_comments : Optional[str], optional
            General comments on the Layout, by default None.
        bathymetry_source : Optional[str], optional
            Description of the data source of the bathymetry, by default None.
        bathymetry_datum : Optional[str], optional
            Vertical datum of the bathymetry, by default None.
        bathymetry_units : Optional[str], optional
            Units of the bathymetry, by default None.

        Raises
        ------
        KeyError
            If schema does not exist in the database.
        LayoutSchemaError
            If the schema exists but has not yet been setup correctly using the
            :func:`layowt.layouts.utils.setup_database_schema` function.
        rasterio.errors.CRSError
            If the EPSG code of the Layout does not match the EPSG of the geometry columns in the
            database.
        ValueError
            If boundary_id is None when upload_boundary=False. boundary_id is not a nullable field
            in the database.
        
        See Also
        --------
        layowt.layouts.utils.setup_database_schema : Function to setup a database schema with the necessary tables for layout upload.
        layowt.layouts.layout.Layout.database_engine : Static method to create a database engine.
        layowt.layouts.utils.database_engine : Function to create a database engine.
        """
        # Check schema exists in database 
        if schema not in sa.inspect(engine).get_schema_names():
           raise KeyError(f"Schema {schema} does not exist in the database.")
        
        if not _valid_schema(engine=engine, schema=schema, floating=False):
            err_msg = f"Schema {schema} is not yet set up with all the necessary system tables for layout upload."
            raise LayoutSchemaError(err_msg)
        
        # Check EPSG of layout is the same as the database EPSG
        get_srid = sa.select(sa.func.Find_SRID(schema, "ow__layowt_turbines", "geom"))
        
        with engine.connect() as conn:
            database_srid = conn.execute(get_srid).fetchone()[0]
        
        if self.epsg != database_srid:
            raise CRSError(f"Mismatch between Layout EPSG ({self.epsg}) and the database Geomatry SRID {database_srid}.")
        
        metadata = sa.MetaData()
        
        layout_table = sa.Table("ow__layowt_layout_info", metadata, autoload_with=engine, schema=schema)
        turbine_table = sa.Table("ow__layowt_turbines", metadata, autoload_with=engine, schema=schema)
        
        if boundary_id is None:
            if upload_boundary is True:
                if boundary_name is None:
                    boundary_name = input("Please enter a descriptive name for the Layout Boundary: \n")
                if boundary_type is None:
                    boundary_type = input("Please enter a descriptive type for the Layout Boundary (e.g. 'Lease Area'): \n")
                if boundary_source is None:
                    boundary_source = input("Please enter the source of the Layout Boundary (e.g. 'Lease Area'): \n")
            else:
                raise ValueError("boundary_id cannot be None when upload_boundary is False.")
        
        if exclusion_id is None:
            if upload_exclusion is True:
                if exclusion_name is None:
                    exclusion_name = input("Please enter a descriptive name for the Layout Turbine Exclusions: \n")
                if exclusion_source is None:
                    exclusion_source = input("Please enter the source of the Layout Turbine Exclusions (e.g. 'Subsea Cables'): \n")
        
        if self.has_bathymetry():
            if bathymetry_source is None:
                bathymetry_source = input("Please enter a description for the bathymetry data source: \n")
            if bathymetry_datum is None:
                bathymetry_datum = input("Please enter the bathymetry vertical datum: \n")
            if bathymetry_units is None:
                bathymetry_units = input("Please enter the bathymetry units: \n")
        
        with engine.connect() as conn:
            # Boundary Upload
            if (upload_boundary) and (boundary_id is None):
                if isinstance(self.area, MultiPolygon):
                    boundary_geom = from_shape(self.area, srid=self.epsg)
                else:
                    boundary_geom = from_shape(MultiPolygon([self.area]), srid=self.epsg)
                    
                boundary_table = sa.Table("ow__layowt_boundaries", metadata, autoload_with=engine, schema=schema)
                
                boundary_insert = sa.insert(boundary_table).values(
                    boundary_name=boundary_name,
                    type=boundary_type,
                    data_source=boundary_source,
                    comments=boundary_comments,
                    geom=boundary_geom
                    )
                
                result = conn.execute(boundary_insert)
                boundary_id = result.inserted_primary_key[0]

            # Exclusion Upload
            if (upload_exclusion) and (exclusion_id is None):
                if isinstance(self.exclusion, MultiPolygon):
                    constraint_geom = from_shape(self.exclusion, srid=self.epsg)
                else:
                    constraint_geom = from_shape(MultiPolygon([self.exclusion]), srid=self.epsg)
                
                constraint_table = sa.Table("ow__layowt_constraint_scenarios", metadata, autoload_with=engine, schema=schema)
                
                constraint_insert = sa.insert(constraint_table).values(
                    constraint_name=exclusion_name,
                    constraint_type="turbine",
                    data_source=exclusion_source,
                    comments=exclusion_comments,
                    geom=constraint_geom
                )
                
                result = conn.execute(constraint_insert)
                exclusion_id = result.inserted_primary_key[0]
            
            # Grid Upload
            if (upload_grid) and (grid_id is None):
                grid_table = sa.Table("ow__layowt_grids", metadata, autoload_with=engine, schema=schema)
                
                grid_geom = from_shape(self.grid.to_multipoint(), srid=self.epsg)
                
                grid_insert = sa.insert(grid_table).values(
                    row_step=self.grid.row_step,
                    col_step=self.grid.col_step,
                    n_rows=int(self.grid.n_rows),
                    n_cols=int(self.grid.n_cols),
                    angle=float(self.grid.angle),
                    row_offset=self.grid.row_offset,
                    col_offset=self.grid.col_offset,
                    x_shear=float(self.grid.x_shear or 0),
                    y_shear=float(self.grid.y_shear or 0),
                    origin_x=self.grid.origin[0],
                    origin_y=self.grid.origin[1],
                    scale= self.grid.scale,
                    comments=grid_comments,
                    geom=grid_geom
                )
                
                result = conn.execute(grid_insert)
                grid_id = result.inserted_primary_key[0]
                
            # Upload Layout info
            layout_insert = sa.insert(layout_table).values(
                n_wtg=self.n_wtg,
                grid_id=grid_id,
                boundary_id=boundary_id,
                turbine_constraint_id=exclusion_id,
                bathymetry_source=bathymetry_source,
                bathymetry_datum=bathymetry_datum,
                bathymetry_units=bathymetry_units,
                comments=layout_comments
            )
            
            result = conn.execute(layout_insert)
            layout_id = result.inserted_primary_key[0]
            
            # Upload Turbines
            turbines_to_insert = []
            for i, turbine in enumerate(self.geom.geoms):
                if self.has_bathymetry():
                    turbines_to_insert.append({"geom": from_shape(turbine, srid=self.epsg), "layout_id": layout_id, "bathymetry": self.samples[i]})
                else:
                    turbines_to_insert.append({"geom": from_shape(turbine, srid=self.epsg), "layout_id": layout_id})

            insert_turbines = sa.insert(turbine_table, turbines_to_insert)
            result = conn.execute(insert_turbines)
                
            conn.commit()
        

def _valid_schema(engine: sa.engine.Engine, schema: str, floating : bool = False) -> bool:
    """
    Check if the target schema already contains all the necessary system layout tables.

    Parameters
    ----------
    engine : sa.engine.Engine
        SQLAlchemy connection engine to the database.
    schema : str, optional
        Schema name, by default schema
    floating : bool
        If the schema should support floating layouts and have a "_moorings" table.

    Returns
    -------
    bool
        If True the schema already has all the required system tables to upload layouts.
    """
    if floating:
        TABLE_NAMES = ["ow__layowt_boundaries", "ow__layowt_constraint_scenarios", "ow__layowt_grids", "ow__layowt_layout_info", "ow__layowt_turbines", "ow__layowt_moorings"]
    else:
        TABLE_NAMES = ["ow__layowt_boundaries", "ow__layowt_constraint_scenarios", "ow__layowt_grids", "ow__layowt_layout_info", "ow__layowt_turbines"]
    
    valid_schema = True
    
    for table_name in TABLE_NAMES:
        if not sa.inspect(engine).has_table(table_name, schema=schema):
            valid_schema = False
    
    return valid_schema


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

def _is_epsg(epsg: int) -> bool:
    """
    Checks if number is a valid EPSG code against the pyproj database.

    Parameters
    ----------
    epsg : int
        number to check against EPSG code database.

    Returns
    -------
    bool
        Ture if valid EPSG number.
    """
    return str(epsg) in pyproj.database.get_codes("EPSG", pyproj.enums.PJType("CRS"))

def _get_epsg_from_postgis(engine, schema: str, table: str, geom_col: str) -> int:
    metadata = sa.MetaData()
    postgis_table = sa.Table(table, metadata, autoload_with=engine, schema=schema)
    geom_col = getattr(postgis_table.c, geom_col) 
    get_srid = sa.select(sa.func.ST_SRID(geom_col))
    
    with engine.connect() as conn:
        result = conn.execute(get_srid)
    
    return result.fetchone()[0]
