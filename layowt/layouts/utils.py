""" This module contains utility functions for the Layout class.
"""
import urllib.parse
from typing import Iterable

import fiona
import geopandas as gp
import h5py
import numpy as np
import pandas as pd
import pyproj
import rasterio
import sqlalchemy as sa
from geoalchemy2 import Geometry
from osgeo import gdal
from rasterio.crs import CRS
from rasterio.errors import CRSError
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import MultiPoint, shape
from shapely.ops import transform
from sqlalchemy.dialects.postgresql import DOUBLE_PRECISION

from layowt.grids.grid import Grid
from layowt.layouts.floating_layout import FloatingLayout
from layowt.layouts.layout import Layout


def geoms_from_shapefile(filepath: str, target_epsg: int | None = None) -> list:
    """
    Loads shapely geometry objects from a shapefile. Can reproject geometries on the fly from the \
    source CRS into the desired CRS defined by its EPSG code in the optional target_epsg argument.

    Parameters
    ----------
    filepath : str
        filepath of the shapefile to load.
    target_epsg : int | None, optional
        `EPSG <https://epsg.io/>`_ code of the target projection for the geometries to be loaded in. By default, None.

    Returns
    -------
    list[Geometry]
        list of shapely geometries contained in the shapefile.
    """
    with fiona.open(filepath) as src:
        geoms = [shape(rec["geometry"]) for rec in src]
        src_crs = pyproj.CRS(src.crs_wkt)
        
    if target_epsg is not None:
        target_crs = pyproj.CRS.from_epsg(target_epsg)
        crs_transformer = pyproj.Transformer.from_crs(src_crs, target_crs, always_xy=True)
        transformed_geoms = []
        for geom in geoms:
            transformed_geom = transform(crs_transformer.transform, geom)
            transformed_geoms.append(transformed_geom)
            
        geoms = transformed_geoms
        
    return geoms


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
    sa.engine.Engine
        SQLAlchemy database engine.
    """
    password = urllib.parse.quote_plus(password)
    
    con_str = f"{dialect}://{username}:{password}@{host}/{db_name}"
    engine = sa.create_engine(con_str, echo=echo, future=True)
    
    return engine

def _valid_schema(engine: sa.engine.Engine, schema: str) -> bool:
    """
    Check if the target schema already contains all the necessary system layout tables.

    Parameters
    ----------
    engine : sa.engine.Engine
        SQLAlchemy connection engine to the database.
    schema : str, optional
        Schema name, by default schema

    Returns
    -------
    bool
        If True the schema already has all the required system tables to upload layouts.
    """
    TABLE_NAMES = ["_boundaries", "_constraint_scenarios", "_grids", "_layout_info", "_turbines"]
    missing_tables = []
    valid_schema = True
    
    for table_name in TABLE_NAMES:
        if not sa.inspect(engine).has_table(table_name, schema=schema):
            valid_schema = False
    
    return valid_schema

def setup_database_schema(engine: sa.engine.Engine, schema: str, epsg: int) -> None:
        
    # # Check schema exists in database (once we upgrade to SQLAlchemy 2.0)
    # if not sa.inspect(engine).has_schema(schema):
    #    raise KeyError(f"Schema {schema} does not exist.")
    
    if _valid_schema(engine, schema):
        print(f'All the necessary layout system tables already exist in schema "{schema}".')
        return None
    
    if not _is_epsg(epsg=epsg):
        raise ValueError(f"Invalid EPSG code {epsg} passed.")
    
    metadata = sa.MetaData()
    
    boundary_table = _create_boundary_table(metadata=metadata, schema=schema, epsg=epsg)
    constraint_table = _create_constraint_table(metadata=metadata, schema=schema, epsg=epsg)
    grid_table = _create_grid_table(metadata=metadata, schema=schema, epsg=epsg)
    layout_table = _create_layout_table(metadata=metadata, schema=schema, epsg=epsg)
    turbine_table = _create_turbine_table(metadata=metadata, schema=schema, epsg=epsg)
    
    metadata.create_all(engine)
    
    print(f'Schema "{schema}" set up successfully!')
    
def _create_boundary_table(metadata: sa.schema.MetaData, schema: str, epsg: int) -> sa.Table:
    
    boundary_table = sa.Table(
        "_boundaries",
        metadata,
        sa.Column("boundary_id", sa.Integer, sa.Identity(always=True), primary_key=True),
        sa.Column("boundary_name", sa.String(50), nullable=False),
        sa.Column("geom", Geometry('MULTIPOLYGON', srid=epsg), nullable=False),
        sa.Column("type", sa.String(50), nullable=False),
        sa.Column("date_created", sa.Date, server_default=sa.func.current_date(), nullable=False),
        sa.Column("created_by", sa.String(50), server_default=sa.func.current_user(), nullable=False),
        sa.Column("comments", sa.Text),
        sa.Index(schema.lower() + "_boundary_id_idx", "boundary_id"),
        schema=schema
    )
    
    return boundary_table
    
def _create_constraint_table(metadata: sa.schema.MetaData, schema: str, epsg: int) -> sa.Table:
    
    constraint_table = sa.Table(
        "_constraint_scenarios",
        metadata,
        sa.Column("constraint_id", sa.Integer, sa.Identity(always=True), primary_key=True),
        sa.Column("constraint_name", sa.String(50), nullable=False),
        sa.Column("geom", Geometry('MULTIPOLYGON', srid=epsg), nullable=False),
        sa.Column("date_created", sa.Date, server_default=sa.func.current_date(), nullable=False),
        sa.Column("created_by", sa.String(50), server_default=sa.func.current_user(), nullable=False),
        sa.Column("comments", sa.Text),
        sa.Index(schema.lower() + "_constraint_id_idx", "constraint_id"),
        schema=schema
    )
    
    return constraint_table

def _create_grid_table(metadata: sa.schema.MetaData, schema: str, epsg: int) -> sa.Table:

    grid_table =sa.Table(
        "_grids",
        metadata,
        sa.Column('grid_id', sa.Integer, sa.Identity(always=True), primary_key=True),
        sa.Column("geom", Geometry('MULTIPOINT', srid=epsg), nullable=False),
        sa.Column("row_step", sa.Numeric(8, 4), nullable=False),
        sa.Column("col_step", sa.Numeric(8, 4), nullable=False),
        sa.Column("n_rows", sa.Integer, nullable=False),
        sa.Column("n_cols", sa.Integer, nullable=False),
        sa.Column("angle", sa.Numeric(6, 4), nullable=False),
        sa.Column("row_offset", sa.Boolean, nullable=False),
        sa.Column("col_offset", sa.Boolean, nullable=False),
        sa.Column("x_shear", sa.Numeric(6, 4), nullable=False),
        sa.Column("y_shear", sa.Numeric(6, 4), nullable=False),
        sa.Column("origin_x", DOUBLE_PRECISION, nullable=False),
        sa.Column("origin_y", DOUBLE_PRECISION, nullable=False),
        sa.Column("scale", sa.Float, nullable=False),
        sa.Column("date_created", sa.Date, server_default=sa.func.current_date(), nullable=False),
        sa.Column("created_by", sa.String(50), server_default=sa.func.current_user(), nullable=False),
        sa.Column("comments", sa.Text),
        sa.Index(schema.lower() + "_grid_id_idx", "grid_id"),
        sa.Index(schema.lower() + "_grid_angle_idx", "angle"),
        sa.Index(schema.lower() + "_grid_row_step_idx", "row_step"),
        sa.Index(schema.lower() + "_grid_col_step_idx", "col_step"),
        schema=schema
    )

    return grid_table

def _create_layout_table(metadata: sa.schema.MetaData, schema: str, epsg: int) -> sa.Table:
    
    layout_table = sa.Table(
        "_layout_info",
        metadata,
        sa.Column("layout_id", sa.Integer, sa.Identity(always=True), primary_key=True),
        sa.Column("n_wtg", sa.Integer, nullable=False),
        sa.Column("grid_id", sa.Integer, sa.ForeignKey(schema + "._grids.grid_id", onupdate="CASCADE", ondelete="SET NULL")),
        sa.Column("boundary_id", sa.Integer, sa.ForeignKey(schema + "._boundaries.boundary_id", onupdate="CASCADE", ondelete="SET NULL"), nullable=False),
        sa.Column("constraint_id", sa.Integer, sa.ForeignKey(schema + "._constraint_scenarios.constraint_id", onupdate="CASCADE", ondelete="SET NULL")),
        sa.Column("date_created", sa.Date, server_default=sa.func.current_date(), nullable=False),
        sa.Column("created_by", sa.String(50), server_default=sa.func.current_user(), nullable=False),
        sa.Column("base_case", sa.Boolean, server_default=sa.false()),
        sa.Column("comments", sa.Text),
        schema=schema
    )
    
    return layout_table

def _create_turbine_table(metadata: sa.schema.MetaData, schema: str, epsg: int) -> sa.Table:
    
    turbine_table = sa.Table(
        "_turbines",
        metadata,
        sa.Column("turbine_id", sa.Integer, sa.Identity(always=True), primary_key=True),
        sa.Column("layout_id", sa.Integer, sa.ForeignKey(schema + "._layout_info.layout_id", onupdate="CASCADE", ondelete="SET NULL"), nullable=False),
        sa.Column("geom", Geometry('POINT', srid=epsg), nullable=False),
        sa.Column("date_created", sa.Date, server_default=sa.func.current_date(), nullable=False),
        sa.Column("created_by", sa.String(50), server_default=sa.func.current_user(), nullable=False),
        sa.Index(schema.lower() + "_turbines_layout_id_idx", "layout_id"),
        schema=schema
    )
    
    return turbine_table

def geoms_from_postgis(
    engine: sa.engine.Engine,
    schema: str,
    table: str,
    geom_col: str = "geom",
    target_epsg: int | None = None,
    **kwargs,
) -> list:
    """
    Loads a list of shapely geometry objects from a PostGIS table.

    Parameters
    ----------
    engine: sa.engine.Engine
        SQLAlchemy database engine.
    schema : str
        name of the schema where the target table is located within the database.
    table : str
        name of the PostGIS table.
    geom_col : str, optional
        name of the geometry column within the table, by default "geom"
    target_epsg : int | None, optional
        `EPSG <https://epsg.io/>`_ code of the target projection for the geometries to be loaded in. By default, None.
    **kwargs : optional
        Additional keyword arguments passed into geopandas.from_postgis.

    Returns
    -------
    list[Geometry]
        list of shapely geometries contained in the PostGIS table.
    """
    # TODO: Should change the geopandas method from_postgis to read_postgis
    new_engine = sa.create_engine(engine.url, echo=False, future=False)
    
    data = gp.read_postgis(
        f'SELECT * from "{schema}"."{table}"',
        con=new_engine,
        geom_col=geom_col,
        **kwargs,
    )
    
    geoms = list(data[geom_col])
    
    if target_epsg is not None:
        if not _is_epsg(target_epsg):
            raise CRSError(f"{target_epsg=} not a valid EPSG code.")
        src_crs = data.crs
        target_crs = pyproj.CRS("EPSG:" + str(target_epsg))
        crs_transformer = pyproj.Transformer.from_crs(src_crs, target_crs, always_xy=True)
        transformed_geoms = []
        for geom in geoms:
            transformed_geom = transform(crs_transformer.transform, geom)
            transformed_geoms.append(transformed_geom)
        
        geoms = transformed_geoms
    
    return geoms


def layouts_to_legacy_csv(layouts: list[Layout], filepath: str = "layouts.csv") -> None:
    """layouts_to_legacy_csv Function for backwards compatibility with legacy multitech code. Exports a lists of layout coordiantes into a .csv file compatible with legacy style OW jupyter notebook codes.
    
    Will be removed in future versions.

    Parameters
    ----------
    layouts : list[Layout]
        A list of Layout objects.
    filepath : str
        Filepath of the csv to write, by default "layouts.csv".
    """
    layout_data = []
    layout_info = []
    for i, layout in enumerate(layouts):
        layout_data += list(zip(np.ones(layout.n_wtg)*i, layout.x, layout.y))
        layout_info.append([i, layout.grid.angle, layout.grid.row_step, layout.grid.col_step,
                            layout.grid.n_rows, layout.grid.n_cols,
                            layout.grid.row_offset, layout.grid.col_offset, layout.grid.x_shear,
                            layout.grid.y_shear, layout.grid.origin[0], layout.grid.origin[1],
                            layout.grid.scale, layout.n_wtg, layout.epsg])
        
    layout_df = pd.DataFrame(layout_data)
    layout_df[3] = 1
    layout_df.columns = ["id", "X", "Y", "center"]
    layout_df.to_csv(filepath, index=False)
    
    layout_info_df = pd.DataFrame(layout_info)
    layout_info_df.columns = ["id", "angle", "row", "col", "n_rows", "n_cols", "row_offset", "col_offset", "x_shear",
                              "y_shear", "origin_x", "origin_y", "scale", "n_wtg", "epsg"]
    layout_info_df.to_csv("INFO_" + filepath, index=False)
    
def legacy_csv_to_layouts(filepath: str, layout_type : str = "fixed") -> list[Layout]:
    """
    This function reads the legacy layouts .csv and the complementary INFO_ file, and rebuilds
    Layout objects from the data. 

    Parameters
    ----------
    filepath : str
        Filepath of the legacy layouts .csv file. Will find the complementary INFO_ file by adding
        "INFO_" as a preffix to this filepath.
    layout_type : str
        Type of layout class to create. Must be "fixed", or "floating". By default, "fixed".

    Returns
    -------
    list[Layout | FloatingLayout]
        List of layouts rebuilt from the .csv files.
    """
    layouts = []
    
    layout_df = pd.read_csv(filepath)
    layout_info = pd.read_csv("INFO_" + filepath)
    
    full_info = True
    grid = None
    
    if layout_type.lower() == "fixed":
        layout_class = Layout
    elif layout_type.lower() == "floating":
        layout_class = FloatingLayout
    else:
        raise ValueError(f"layout_type must be 'fixed' or 'floating'. '{layout_type}' passed.")
    
    for row in layout_info.itertuples():
        if full_info:
            try:
                grid = Grid(n_rows=row.n_rows, n_cols=row.n_cols, row_step=row.row, col_step=row.col,
                            row_offset=row.row_offset, col_offset=row.col_offset,
                            angle=row.angle, x_shear=row.x_shear, y_shear=row.y_shear, origin=(row.origin_x, row.origin_y),
                            scale=row.scale)
            except AttributeError:
                full_info = False
                grid = None
        
        wtg_coords = layout_df[layout_df["id"] == row.id]
        geom = MultiPoint(list(zip(wtg_coords["X"], wtg_coords["Y"])))
        layout = layout_class(grid=grid)
        layout._raw_geom = geom
        layout.geom = geom
        
        if full_info:
            layout.epsg = row.epsg
            
        layouts.append(layout)
    
    return layouts

def reproject_raster(filepath: str, output_path: str, target_epsg: int, resample_method: int = 0) -> None:
    """reproject_raster Utility function to reproject raster datasets into the CRS defined by the user input EPSG code. Can select from a variety of reampling methods.

    Parameters
    ----------
    filepath : str
        filepath of the raster to be reprojected. Must be a format supported by the rasterio.open function.
    output_path : str
        filepath of the reprojected raster dataset.
    target_epsg : int
        `EPSG <https://epsg.io/>`_ code of the target projection for the reprojected raster dataset.
    resample_method : int, optional
        Resampling algorithm. Integer value used by the :class:`rasterio.enums.Resampling` enumerator class to select the algorithm, by default 0.    
        The mapping of values to resampling algorithm is the following:
            * nearest = 0
            * bilinear = 1
            * cubic = 2
            * cubic_spline = 3
            * lanczos = 4
            * average = 5
            * mode = 6
            * gauss = 7
            * max = 8
            * min = 9
            * med = 10
            * q1 = 11
            * q3 = 12
            * sum = 13
            * rms = 14
    
    See Also
    --------
    rasterio.warp.reproject : Rasterio module for warping and reprojection of raster datasets.
    rasterio.enums.Resampling : Rasterio warp resampling algorithms.
    """
    dst_crs = CRS.from_epsg(target_epsg)
    
    with rasterio.open(filepath) as src:
        transform, width, height = calculate_default_transform(src.crs, dst_crs, src.width, src.height, *src.bounds)
        kwargs = src.meta.copy()
        kwargs.update({"crs": dst_crs,
                       "transform": transform,
                       "width": width,
                       "height": height
        })
    
        with rasterio.open(output_path, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling(resample_method)
                )

def raster_slope(filepath: str, output_path: str = "slope.tif", slope_format: str = "degree", scale: float = 1) -> None:
    """
    This function takes a raster elevation dataset and calculates the slope in degrees or percent.
    
    Parameters
    ----------
    filepath : str
        filepath of the raster to be processed. Must be a file format supported by GDAL.
    output_path : str, optional
        filepath of the slope raster dataset, by default "slope.tif".
    slope_format : str, optional
        Format of the output slope calculation. Either in degrees or percent, by default "degree".
    scale : float, optional
        Ratio of vertical units to horizontal, be default 1.
    """
    if slope_format.lower() not in ["degree", "percent"]:
        ValueError("slope_format but be either 'degree' or 'percent'.")
        
    input_dem = gdal.Open(filepath)
    gdal.DEMProcessing(output_path, input_dem, "slope", slopeFormat=slope_format, scale=scale)
    input_dem = None

def _get_grid_params(layout):
    return (layout.grid.angle, layout.grid.row_step, layout.grid.col_step,
            layout.grid.n_rows, layout.grid.n_cols, layout.grid.row_offset,
            layout.grid.col_offset, layout.grid.x_shear, layout.grid.y_shear,
            layout.grid.origin[0], layout.grid.origin[1], layout.grid.scale,
            layout.n_wtg, layout.epsg)

def _layouts_to_coords_array(layouts, start: int = 0):
    layout_dt = np.dtype({"names": ["id", "x", "y"],
                    'formats': [(float)]*4})
    
    layout_ids = [np.ones(layout.n_wtg)*idx for idx, layout in enumerate(layouts, start=start)]
    layout_xs = []
    layout_ys = []
    for layout in layouts:
        coords = layout.coords
        layout_xs.append(coords[:, 0])
        layout_ys.append(coords[:, 1])
    layout_ids = np.concatenate(layout_ids)
    layout_xs = np.concatenate(layout_xs)
    layout_ys = np.concatenate(layout_ys)
    
    layouts_arr = np.stack((layout_ids, layout_xs, layout_ys), axis=1)
    layouts_arr.dtype=layout_dt
    
    return layouts_arr

def _layouts_to_info_array(layouts, start: int = 0):
    vect_get_grid_params = np.vectorize(_get_grid_params)
    
    info_array = (np.arange(start, start + len(layouts)), *vect_get_grid_params(np.array(layouts)))
    info_array = np.stack(info_array, axis=1)
    
    names = ["id", "angle", "row_step", "col_step", "n_rows", "n_cols", "row_offset", "col_offset", "x_shear", "y_shear", "origin_x", "origin_y", "scale", "n_wtg", "epsg"]
    formats = [(float)]*15
    info_dt = np.dtype({"names": names,
                    'formats': formats})
    info_array.dtype = info_dt
    
    return info_array

def layouts_to_hdf5(layouts: list[Layout | FloatingLayout], filepath: str, mode: str = 'a', start_id: int = 0) -> None:
    """
    Function to save layouts to an hdf5 file.

    Parameters
    ----------
    layouts : list[Layout  |  FloatingLayout]
        List of layouts to save to disk.
    filepath : str
        Filepath of the hdf5 file to save layouts to.
    mode : str, optional
        File creation mode, by default "a". Valid modes are:
            * "r+": Read/write, file must exist
            * "w" : Create file, truncate if exists
            * "w-": Create file, fail if exists
            * "x" : Create file, fail if exists
            * "a" : Read/write if exists, create otherwise
    start_id : int, optional
        Number at which layout ids start counting from in the hdf5 file, by default 0.
    
    Raises
    ------
    ValueError
        If mode is not valid.
    """
    if mode not in ["r+", "w", "w-", "x", "a"]:
        raise ValueError(f"Invalid value for mode passed. {mode=} passed.")
    
    layouts_array = _layouts_to_coords_array(layouts, start=start_id)
    info_array = _layouts_to_info_array(layouts, start=start_id)
    
    with h5py.File(filepath, mode) as f:
        if "coordinates" not in f:
            f.create_dataset("coordinates", data=layouts_array, maxshape=(None, None))
        else:
            f["coordinates"].resize((f["coordinates"].shape[0] + layouts_array.shape[0]), axis = 0)
            f["coordinates"][-layouts_array.shape[0]:] = layouts_array
        
        if "info" not in f:
            f.create_dataset("info", data=info_array, maxshape=(None, None))
        else:
            f["info"].resize((f["info"].shape[0] + info_array.shape[0]), axis = 0)
            f["info"][-info_array.shape[0]:] = info_array

def layouts_from_hdf5(filepath: str, start_id: int = 0, end_id: int | None = None,
                      indexes: int | Iterable[int] | None = None, layout_type : str = "fixed") -> list[Layout | FloatingLayout]:
    """
    Function to read an hdf5 file and create layouts from the stored data. Can return Layout or \
    FloatingLayout classes.

    Parameters
    ----------
    filepath : str
        Filepath of the hdf5 file to read.
    start_id : int, optional
        Minimum bound of the layout id used to slice the hdf5 file. Uses the numpy convention of an inclusive lower limit. By default 0.
    end_id : int | None, optional
        Maximum bound of the layout id used to slice the hdf5 file. Uses the numpy convention of an non-inclusive upper limit. By default None.
    indexes : int | Iterable[int] | None, optional
        Integer or Iterable of integers of the layout indexes to read from the file. If a value is
        passed, the start_id and end_id argument are ignored. To extract a single layout it will be
        faster to set start_id and end_id to the same value rather than use the indexes argument.
        By default None.
    layout_type : str, optional
        Type of layout class to create. Must be "fixed", or "floating". By default, "fixed".

    Returns
    -------
    list[Layout | FloatingLayout]
        List of layouts rebuilt from the .csv files.

    Raises
    ------
    ValueError
        If layout_type is not "fixed" OR "floating".
    ValueError
        If start_id and/or end_id are out of range of the hdf5 file and if start_id and end_id are equal.
    """
    
    if layout_type.lower() == "fixed":
        layout_class = Layout
    elif layout_type.lower() == "floating":
        layout_class = FloatingLayout
    else:
        raise ValueError(f"layout_type must be 'fixed' or 'floating'. '{layout_type}' passed.")
    
    layouts = []
    
    with h5py.File(filepath, 'r') as f:
        if indexes is None:
            if end_id is None:
                end_id = f['info'].shape[0]

            if (start_id >= f['info'].shape[0]) or (end_id > f['info'].shape[0]) or (start_id == end_id):
                raise ValueError(f"Invalid start and/or end ids passed. {start_id=} and {end_id=} passed for file of shape {f['info'].shape}.")

            if end_id - start_id == 1:
                layout_info = f["info"][start_id:end_id][0]
            else:
                layout_info = np.squeeze(f["info"][start_id:end_id])

            if start_id == 0:
                coords_start = 0
            else:
                coords_start = int(f["info"][0:start_id]["n_wtg"].sum())

            if end_id is None:
                coords_end = f["coordinates"].shape[0]
            else:
                coords_end = int(coords_start + layout_info["n_wtg"].sum())

            coords_array = np.squeeze(f["coordinates"][coords_start:coords_end])
        else:
            if isinstance(indexes, int) or len(indexes) == 1:
                layout_info = f["info"][indexes]
            else:
                indexes = np.array(indexes).astype(int)
                layout_info = np.squeeze(f["info"][indexes])
            coords_array = np.squeeze(f["coordinates"][np.isin(f["coordinates"]["id"], indexes)])
    
    for info_row in layout_info:
        row_dict = {key: val for key, val in zip(layout_info.dtype.names, info_row)}
        row_dict.pop("id")
        row_dict.pop("epsg")
        row_dict.pop("n_wtg")
        row_dict["origin"] = (row_dict["origin_x"], row_dict["origin_y"])
        row_dict.pop("origin_x")
        row_dict.pop("origin_y")
        grid = Grid(**row_dict)
        
        geom = MultiPoint(list(zip(coords_array[coords_array["id"] == info_row["id"]]["x"], coords_array[coords_array["id"] == info_row["id"]]["y"])))
        
        layout = layout_class(grid=grid)
        layout.geom = geom
        layout.epsg = info_row["epsg"]
        layouts.append(layout)
    
    return layouts

def hdf5_to_dataframe(filepath: str, start_id: int = 0, end_id: int | None = None,
                      indexes: int | Iterable[int] | None = None, out: str = "both") -> tuple[pd.DataFrame, pd.DataFrame] | pd.DataFrame:
    """
    Function to read an hdf5 layouts file and create pandas DataFrames from the stored data.

    Parameters
    ----------
    filepath : str
        Filepath of the hdf5 file to read.
    start_id : int, optional
        Minimum bound of the layout id used to slice the hdf5 file. Uses the numpy convention of an inclusive lower limit. By default 0.
    end_id : int | None, optional
        Maximum bound of the layout id used to slice the hdf5 file. Uses the numpy convention of an non-inclusive upper limit. By default None.
    indexes : int | Iterable[int] | None, optional
        Integer or Iterable of integers of the layout indexes to read from the file. If a value is
        passed, the start_id and end_id argument are ignored. To extract a single layout it will be
        faster to set start_id and end_id to the same value rather than use the indexes argument.
        By default None.
    out : str, optional
        Select what outputs are returned. If "both" then both the coordinates and info dataframe are returned.
        If "info" or "coordinates", only the info or coordinates dataframes will be returned respectively.
        Case insensitive. By default, "both".

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple (coords_df, info_df) where the coords_df is the DataFrame of all turbine coordinates
        and info_df is the DataFrame with all the layout information.

    Raises
    ------
    ValueError
        If invalid value for out was passed.
    ValueError
        If start_id and/or end_id are out of range of the hdf5 file and if start_id and end_id are equal.
    """
    
    if out.lower() not in ["both", "info", "coordinates"]:
        raise ValueError(f"Invalid value for out. Must be one of 'both', 'info', or 'coordinates'. {out=} was passed.")
    
    with h5py.File(filepath, 'r') as f:
        if indexes is None:
            if end_id is None:
                end_id = f['info'].shape[0]

            if (start_id >= f['info'].shape[0]) or (end_id > f['info'].shape[0]) or (start_id == end_id):
                raise ValueError(f"Invalid start and/or end ids passed. {start_id=} and {end_id=} passed for file of shape {f['info'].shape}.")

            if end_id - start_id == 1:
                layout_info = f["info"][start_id:end_id][0]
            else:
                layout_info = np.squeeze(f["info"][start_id:end_id])
            
            if out.lower() == "info":
                info_df = pd.DataFrame(layout_info)
                info_df["row_offset"] = info_df["row_offset"].astype(bool)
                info_df["col_offset"] = info_df["col_offset"].astype(bool)
                return info_df
            else:
                if start_id == 0:
                    coords_start = 0
                else:
                    coords_start = int(f["info"][0:start_id]["n_wtg"].sum())

                if end_id is None:
                    coords_end = f["coordinates"].shape[0]
                else:
                    coords_end = int(coords_start + layout_info["n_wtg"].sum())

                coords_array = np.squeeze(f["coordinates"][coords_start:coords_end])
                coords_df = pd.DataFrame(coords_array)
                if out.lower() == "coordinates":
                    return coords_df
                else:
                    info_df = pd.DataFrame(layout_info)
                    info_df["row_offset"] = info_df["row_offset"].astype(bool)
                    info_df["col_offset"] = info_df["col_offset"].astype(bool)
                    return (coords_df, info_df)
        else:
            if out.lower() in ["both", "info"]:
                if isinstance(indexes, int) or len(indexes) == 1:
                    layout_info = f["info"][indexes]
                else:
                    indexes = np.array(indexes).astype(int)
                    layout_info = np.squeeze(f["info"][indexes])
                    
                info_df = pd.DataFrame(layout_info)
                info_df["row_offset"] = info_df["row_offset"].astype(bool)
                info_df["col_offset"] = info_df["col_offset"].astype(bool)
                
                if out.lower() == "info":
                    return info_df
                    
                coords_array = np.squeeze(f["coordinates"][np.isin(f["coordinates"]["id"], indexes)])
            
                coords_df = pd.DataFrame(coords_array)
                return (coords_df, info_df)
            if out.lower() in ["both", "coordinates"]:
                if not(isinstance(indexes, int) or len(indexes)) == 1:
                    indexes = np.array(indexes).astype(int)
                coords_array = np.squeeze(f["coordinates"][np.isin(f["coordinates"]["id"], indexes)])
                
                coords_df = pd.DataFrame(coords_array)
                
                if out.lower() == "coordinates":
                    return coords_df
                
                if isinstance(indexes, int) or len(indexes) == 1:
                    layout_info = f["info"][indexes]
                else:
                    layout_info = np.squeeze(f["info"][indexes])
                    
                info_df = pd.DataFrame(layout_info)
                info_df["row_offset"] = info_df["row_offset"].astype(bool)
                info_df["col_offset"] = info_df["col_offset"].astype(bool)
                
                return (coords_df, info_df)

def _is_epsg(epsg: int) -> bool:
    return str(epsg) in pyproj.database.get_codes("EPSG", pyproj.enums.PJType("CRS"))
