""" This module contains utility functions for the Layout class.
"""

import fiona
import geopandas as gp
from shapely.geometry import shape
from sqlalchemy import create_engine


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
