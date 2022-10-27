""" This module contains utility functions for the Layout class.
"""

import fiona
import geopandas as gp
import h5py
import numpy as np
import pandas as pd
import pyproj
import rasterio
from layowt.grids.grid import Grid
from layowt.layouts.floating_layout import FloatingLayout
from layowt.layouts.layout import Layout
from osgeo import gdal
from rasterio.crs import CRS
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import MultiPoint, shape
from shapely.ops import transform
from sqlalchemy import create_engine


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


def geoms_from_postgis(
    username: str,
    password: str,
    schema: str,
    table: str,
    host: str = "ow-postgre.postgres.database.azure.com",
    db_name: str = "corp_ta_ea",
    geom_col: str = "geom",
    target_epsg: int | None = None,
    **kwargs,
) -> list:
    """
    Loads a list of shapely geometry objects from a PostGIS table.

    Parameters
    ----------
    username : str
        username used to log into the PostGIS database connection.
    password : str
        password used to log into the PostGIS database connection.
    schema : str
        name of the schema where the target table is located within the database.
    table : str
        name of the PostGIS table.
    host : str, optional
        host used to connect to the PostGIS database, by default "ow-postgre.postgres.database.azure.com"
    db_name : str, optional
        name of the PostGIS database within the hose, by default "corp_ta_ea"
    geom_col : str, optional
        name of the geometry column within the table, by default "geom"
    target_epsg : int | None, optional
        `EPSG <https://epsg.io/>`_ code of the target projection for the geometries to be loaded in. By default, None.

    Returns
    -------
    list[Geometry]
        list of shapely geometries contained in the PostGIS table.
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
    
    geoms = list(data[geom_col])
    
    if target_epsg is not None:
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
        
        if "info" not in f:
            f.create_dataset("info", data=info_array, maxshape=(None, None))
        else:
            f["info"].resize((f["info"].shape[0] + info_array.shape[0]), axis = 0)

def layouts_from_hdf5(filepath: str, start_id: int = 0, end_id: int | None = None, layout_type : str = "fixed") -> list[Layout | FloatingLayout]:
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
        
        if (start_id >= f['info'].shape[0]) or (end_id > f['info'].shape[0]) or (start_id == end_id):
            raise ValueError(f"Invalid start and/or end ids passed. {start_id=} and {end_id=} passed for file of shape {f['info'].shape}.")
        
        layout_info = np.squeeze(f["info"][start_id:end_id])
                                 
        if start_id == 0:
            coords_start = 0
        else:
            coords_start = int(f["info"][0:start_id]["n_wtg"].sum())
        
        if end_id is None:
            coords_end = f["coordinates"].shape[0]
        else:
            coords_end = int(coords_start + layout_info["n_wtg"].sum() + 1)
        
        coords_array = np.squeeze(f["coordinates"][coords_start:coords_end])
    
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

def hdf5_to_dataframe(filepath: str, start_id: int = 0, end_id: int | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
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

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        Tuple (coords_df, info_df) where tehe coords_df is the DataFrame of all turbine coordinates
        and info_df is the DataFrame with all the layout information.

    Raises
    ------
    ValueError
        If start_id and/or end_id are out of range of the hdf5 file and if start_id and end_id are equal.
    """
    
    with h5py.File(filepath, 'r') as f:
        
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
        
        coords_df = pd.DataFrame(coords_array)
        info_df = pd.DataFrame(layout_info)
    
    return (coords_df, info_df)
