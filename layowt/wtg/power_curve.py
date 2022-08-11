"""
This module contains functionality to load, create, and manipulate power curves from the internal
OW SQL Server database, text files, etc.

Can also be used to export power curves into common formats such as WAsP .wtg files, PyWake, etc.
"""

import os
import re
import urllib
import xml.etree.ElementTree as ET
from datetime import date
from typing import Iterable, Optional

import pandas as pd
import sqlalchemy as sa
from layowt import __version__


class PowerCurve:
    def __init__(self, wind_speed: Iterable[float], power: Iterable[float], ct: Iterable[float],
                 rotor_diameter: Optional[float]=None, hub_height: Optional[float]=None, air_density: Optional[float]=None,
                 cut_in: Optional[float]=None, rated_ws: Optional[float]=None, cut_out: Optional[float]=None,
                 curve_id: Optional[str]=None, manufacturer: Optional[str]=None):

        data = {"wind_speed": wind_speed, "power": power, "ct": ct}

        self.data = pd.DataFrame(data=data)
        self.data.set_index("wind_speed", inplace=True)
        self.attributes = {"rotor_diameter": rotor_diameter,
                           "rated_power": self.data["power"].max(),
                           "hub_height": hub_height,
                           "air_density": air_density,
                           "cut_in": cut_in,
                           "rated_ws": rated_ws,
                           "cut_out": cut_out,
                           "curve_id": curve_id,
                           "manufacturer": manufacturer,
                           "curve_date": date.today(),
                           "ref_document_id": "Unknown",
                           "revision": 0,
                           "curve_type": None}

    def set_hub_height(self, hub_height: float):
        if self.attributes["rotor_diameter"] is not None and hub_height < self.attributes["rotor_diameter"] / 2:
            raise ValueError(f"Hub Height cannot be less than half the rotor diameter. {hub_height=} passed for {self.attributes['rotor_diameter']=}.")
        self.attributes["hub_height"] = hub_height
    
    
    @staticmethod
    def database_engine(username: Optional[str]=None):
        
        if username is None:
            # gets the username from the OS, splits at capital letters using regex, and forms the OW email addres
            username = '.'.join([string for string in re.split('(?=[A-Z])', os.getlogin()) if string]) +"@oceanwinds.com"
        
        server = 'owsql.public.35b251897c24.database.windows.net,3342'
        database = 'CORP_WTG'
        
        params = urllib.parse.quote_plus('DRIVER={ODBC Driver 17 for SQL Server};SERVER='+server+';DATABASE='+database+';UID='+username+';Authentication=ActiveDirectoryInteractive')
        engine = sa.create_engine(f"mssql+pyodbc:///?odbc_connect={params}")
        
        return engine
    
    @classmethod
    def from_commercial_database(cls, curve_id: str, engine: Optional[sa.engine.Engine]=None, username: Optional[str]=None):
        
        if engine is None:
            engine = PowerCurve.database_engine(username=username)
        
        pc_df = pd.read_sql(f"SELECT * FROM power_curves.commercial_power_curves WHERE curve_id = '{curve_id}'", engine)
        
        pc_attrs = pd.read_sql(f"SELECT * FROM power_curves.commercial_turbines WHERE curve_id = '{curve_id}'", engine)
        
        pc = cls(wind_speed=pc_df["wind_speed"], power=pc_df["power_kw"], ct=pc_df["ct"], rotor_diameter=pc_attrs["rotor_diameter_m"][0],
                 air_density=pc_attrs["air_density"][0], cut_in=pc_attrs["cut_in_ws"][0], cut_out=pc_attrs["cut_out_ws"][0],
                 rated_ws=pc_attrs["rated_ws"][0], curve_id=curve_id, manufacturer=pc_attrs["manufacturer"][0])
        
        pc.attributes["curve_date"] = pc_attrs["curve_date"][0]
        pc.attributes["ref_document_id"] = pc_attrs["ref_document_id"][0]
        pc.attributes["revision"] = pc_attrs["revision"][0]
        pc.attributes["curve_type"] = "commercial"
        
        return pc
    
    @classmethod
    def from_internal_database(cls, curve_id: str, engine: Optional[sa.engine.Engine]=None, username: Optional[str]=None):
        # TODO: Implement method once the internal databse is set up.
        pass
    
    def to_wtg(self, filepath: str="power_curve.wtg"):
        
        if self.attributes["hub_height"] is None:
            raise ValueError("Must set hub height before creating .wtg file.")
        
        data = self.data[self.data["power"] > 0].copy()
        data["power"] = data["power"] * 1000    # converting from kW to W for WAsP format.
        
        # XML file structure
        wtg = ET.Element('WindTurbineGenerator')
        wtg.set('FormatVersion','1.02')
        wtg.set('LowestCompatibleFormatVersion','1.01')
        wtg.set('Description', self.attributes["curve_id"])
        wtg.set('ManufacturerName', self.attributes["manufacturer"])
        wtg.set('ReferenceURI','')
        wtg.set('RotorDiameter', str(self.attributes["rotor_diameter"]))
        wtg.set('ControlSystem','Pitch')
        wtg.set('RatedPower','%1.2e' % (1000 * self.attributes["rated_power"])) # kW to W (WAsP format) with power(10) format too.
        
        comment_1 = ET.SubElement(wtg, 'Comments')
        comment_1.text = 'This XML file has been generated automatically by LayOWt version' + __version__
        
        h_heights = ET.SubElement(wtg, 'SuggestedHeights')
        h_height = ET.SubElement(h_heights, 'Height')
        h_height.text = str(self.attributes["hub_height"])

        performance_table = ET.SubElement(wtg, 'PerformanceTable')
        performance_table.set('AirDensity', str(self.attributes["air_density"]))
        performance_table.set('MaximumNoiseLevel','0')
        performance_table.set('BladePitchAngle','-99999')
        performance_table.set('ReleaseDate', str(self.attributes["curve_date"]))
        performance_table.set('DataStatus', str(self.attributes["revision"]))
        performance_table.set('DataSource', self.attributes["ref_document_id"])
        performance_table.set('ReferenceURI','None')
        performance_table.set('StationaryThrustCoEfficient','%.3g' % data["ct"].iloc[-1]) # Ct of the last wind speed bin
        
        start_stop = ET.SubElement(performance_table, 'StartStopStrategy')
        start_stop.set('LowSpeedCutOut', str(self.attributes["cut_in"]))
        start_stop.set('LowSpeedCutIn', str(self.attributes["cut_in"]))
        start_stop.set('HighSpeedCutOut', str(self.attributes["cut_out"]))
        start_stop.set('HighSpeedCutIn', str(self.attributes["cut_out"]))

        comment_2 = ET.SubElement(performance_table, 'Comments')
        comment_2.text = f'Power curve from CORP_WTG.power_curves.{self.attributes["curve_type"]}_power_curves'
        
        data_table = ET.SubElement(performance_table, 'DataTable')
        for ws, power, ct in zip(data.index, data["power"], data["ct"]):
            data_point = ET.SubElement(data_table, 'DataPoint')
            data_point.set('WindSpeed','%.3g' % ws)
            data_point.set('PowerOutput','%1.3e' % power) # power(10) format
            data_point.set('ThrustCoEfficient','%.3g' % ct)
        
        context_information = ET.SubElement(wtg, 'ContextInformation')
        certificate_checksum = ET.SubElement(context_information, 'CertificateChecksum')

        # Write the XML file with the results
        mydata = ET.tostring(wtg,encoding='unicode',method='xml')
        with open(filepath, 'w', encoding="utf8") as myfile:
            myfile.write(mydata)
