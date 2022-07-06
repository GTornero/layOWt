""" This module contains utility functions for the Layout class.
"""

import fiona
import geopandas as gp
import numpy as np
import pandas as pd
import pyproj
import rasterio
from rasterio.crs import CRS
from rasterio.warp import Resampling, calculate_default_transform, reproject
from shapely.geometry import shape
from shapely.ops import transform
from sqlalchemy import create_engine

from .layout import Layout


def geoms_from_shapefile(filepath: str, target_epsg: int | None = None) -> list:
    """geoms_from_shapefile loads shapely geometry objects from a shapefile. Can reproject geometries on the fly from the source CRS into the desired CRS defined by its EPSG code in the optional target_epsg argument.

    Parameters
    ----------
    filepath : str
        filepath of the shapefile to load.
    target_epsg : int | None, optional
        `EPSG <https://epsg.io/>`_ code of the target projection for the geometries to be loaded in. By default, None.

    Returns
    -------
    list
        list of shapely geometries contained in the shapefile.
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
    """geoms_from_postgis loads a list of shapely geometry objects from a PostGIS table.

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
    list
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
        layout_info.append([i, layout.grid.angle, layout.grid.row_step, layout.grid.col_step])
        
    layout_df = pd.DataFrame(layout_data)
    layout_df[3] = 1
    layout_df.columns = ["id", "X", "Y", "center"]
    layout_df.to_csv(filepath, index=False)
    
    layout_info_df = pd.DataFrame(layout_info)
    layout_info_df.columns = ["id", "angle", "row", "col"]
    layout_info_df.to_csv("INFO_" + filepath, index=False)

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
        Resampling algorithm. Integer value used by the rasterio.enums.Resampling enumerator class to select the algorithm, by default 0.    
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
