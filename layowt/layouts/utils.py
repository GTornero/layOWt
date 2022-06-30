""" This module contains utility functions for the Layout class.
"""

import fiona
import geopandas as gp
import numpy as np
import pandas as pd
from shapely.geometry import shape
from sqlalchemy import create_engine

from .layout import Layout


def geoms_from_shapefile(filepath: str) -> list:
    """geoms_from_shapefile loads shapely geometry objects from a shapefile.

    Parameters
    ----------
    filepath : str
        filepath of the shapefile to load.

    Returns
    -------
    list
        list of shapely geometries contained in the shapefile.
    """
    with fiona.open(filepath) as src:
        areas = [shape(rec["geometry"]) for rec in src]
    return areas


def geoms_from_postgis(
    username: str,
    password: str,
    schema: str,
    table: str,
    host: str = "ow-postgre.postgres.database.azure.com",
    db_name: str = "corp_ta_ea",
    geom_col="geom",
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
