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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sqlalchemy as sa
from layowt import __version__
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.interpolate import PchipInterpolator


class PowerCurve:
    """
    The PowerCurve class is used to represent the power curve performance of turbines and hold its attributes.
    The class also contains many useful utility methods and functions to connnect and query the internal OW
    corporate power curve database.

    Parameters
    ----------
    wind_speed : Iterable[float]
        Iterable of wind speed bins, must be the same length as `power` and `ct`.
    power : Iterable[float]
        Iterable of active power values in kW for each element in `wind_speed`, must be the same length as `wind_speed` and `ct`.
    ct : Iterable[float]
        Iterable of thrust coefficient values for each element in `wind_speed`, must be the same length as `wind_speed` and `power`.
    rotor_diameter : Optional[float], optional
        Rotor diameter in m, by default None.
    hub_height : Optional[float], optional
        Hub height in m, required to export to WAsP .wtg format, by default None.
    air_density : Optional[float], optional
        Air density for which the power curve is defined in kg/m3, by default None.
    cut_in : Optional[float], optional
        Cut-in wind speed in m/s, by default None.
    rated_ws : Optional[float], optional
        Lowest wind speed in m/s at which the power curve first reaches rated power, by default None.
    cut_out : Optional[float], optional
        Cut-out wind speed in m/s, by default None.
    curve_id : Optional[str], optional
        ID of the power curve for traceability purposes, by default None.
    manufacturer : Optional[str], optional
        Manufacturer of the power curve for traceability purposes, by default None.
    
    Attributes
    ----------
    data : pandas.DataFrame
        Attribute: DataFrame of the Power and Thrust curves.
    attributes : dict
        Attribute: Dictionary of the turbine attributes.
    """
    
    def __init__(self, wind_speed: Iterable[float], power: Iterable[float], ct: Iterable[float],
                 rotor_diameter: Optional[float]=None, hub_height: Optional[float]=None, air_density: Optional[float]=None,
                 cut_in: Optional[float]=None, rated_ws: Optional[float]=None, cut_out: Optional[float]=None,
                 curve_id: Optional[str]=None, manufacturer: Optional[str]=None) -> None:
        
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
        
        if cut_in is None:
            self.attributes["cut_in"] = cut_in = self.data[self.data["power"]!=0].index[0]
        if cut_out is None:
            self.attributes["cut_out"] = cut_in = self.data[self.data["power"]!=0].index[-1]

    def set_hub_height(self, hub_height: float) -> None:
        """
        This method sets the hub height for the power curve.
        If a rotor diameter is set for the power curve, the method checks that the hub height is
        not lower than half the rotor diameter.

        Parameters
        ----------
        hub_height : float
            Hub height in metres.

        Raises
        ------
        ValueError
            If the hub height is less than the rotor diameter if a rotor diameter is defined.
        """
        
        if self.attributes["rotor_diameter"] is not None and hub_height < self.attributes["rotor_diameter"] / 2:
            raise ValueError(f"Hub Height cannot be less than half the rotor diameter. {hub_height=} passed for {self.attributes['rotor_diameter']=}.")
        self.attributes["hub_height"] = hub_height
    
    
    @staticmethod
    def database_engine(username: Optional[str]=None) -> sa.engine.Engine:
        """
        Method used to create a database connection to the OW corporate power curve database using
        Active Directory authentication. When first called in a session, a authentication window
        will pop up asking for user credentials to connect to the database.
        
        The OW VPN must be active for a successful database connection.

        Parameters
        ----------
        username : Optional[str], optional
            OW corporate username in the @oceanwinds.com domain. If no username is provided, the
            username is derived friom the Windows user login name. By default None.

        Returns
        -------
        sqlalchemy.engine.Engine
            SqlAlchemy connection engine to the internal OW corporate power curve database.
        """
        
        if username is None:
            # gets the username from the OS, splits at capital letters using regex, and forms the OW email addres
            username = '.'.join([string for string in re.split('(?=[A-Z])', os.getlogin()) if string]) +"@oceanwinds.com"
        
        server = 'owsql.public.35b251897c24.database.windows.net,3342'
        database = 'CORP_WTG'
        
        params = urllib.parse.quote_plus('DRIVER={ODBC Driver 17 for SQL Server};SERVER='+server+';DATABASE='+database+';UID='+username+';Authentication=ActiveDirectoryInteractive')
        engine = sa.create_engine(f"mssql+pyodbc:///?odbc_connect={params}")
        
        return engine
    
    @classmethod
    def from_database(cls, curve_id: str, engine: Optional[sa.engine.Engine]=None, username: Optional[str]=None) -> "PowerCurve":
        """
        Alternative constructor method to create a :class:`PowerCurve` instance directly from the OW
        corporate power curve database by passing a unique curve_id.

        Parameters
        ----------
        curve_id : str
            Unique ID for a power curve in the OW corporate power curve database.
        engine : Optional[sqlalchemy.engine.Engine], optional
            SqlAlchemy connection engine to the OW corporate power curve database. If no open
            connection is passed, the static method :meth:`PowerCurve.database_engine()` will be
            called implicitly using `username` if passed, otherwise the username is derived from
            the Windows user login name. By default None.
        username : Optional[str], optional
            OW corporate username in the @oceanwinds.com domain. If no username is provided, the
            username is derived friom the Windows user login name. This argument is only required
            if no engine is passed. By default None.

        Returns
        -------
        PowerCurve
            PowerCurve instance built from the OW corporate power curve database query result.
        """
        
        curve_type = "commercial"
        
        if engine is None:
            engine = PowerCurve.database_engine(username=username)
        
        pc_attrs = pd.read_sql(f"SELECT * FROM power_curves.{curve_type}_turbines WHERE curve_id = '{curve_id}'", engine)
        
        if pc_attrs.empty:
            curve_type = "internal"
            pc_attrs = pd.read_sql(f"SELECT * FROM power_curves.{curve_type}_turbines WHERE curve_id = '{curve_id}'", engine)
        
            if pc_attrs.empty:
                raise ValueError(f"No {curve_id=} found in the database. Could not create a PowerCurve instance.")
        
        pc_df = pd.read_sql(f"SELECT * FROM power_curves.{curve_type}_power_curves WHERE curve_id = '{curve_id}'", engine)
        
        pc = cls(wind_speed=pc_df["wind_speed"], power=pc_df["power_kw"], ct=pc_df["ct"], rotor_diameter=pc_attrs["rotor_diameter_m"][0],
                 air_density=pc_attrs["air_density"][0], cut_in=pc_attrs["cut_in_ws"][0], cut_out=pc_attrs["cut_out_ws"][0],
                 rated_ws=pc_attrs["rated_ws"][0], curve_id=curve_id, manufacturer=pc_attrs["manufacturer"][0])
        
        pc.attributes["curve_date"] = pc_attrs["curve_date"][0]
        pc.attributes["ref_document_id"] = pc_attrs["ref_document_id"][0] # FIXME: this might throw an error if internal turbines dont have a reference doc.
        pc.attributes["revision"] = pc_attrs["revision"][0]
        pc.attributes["curve_type"] = curve_type.title()
        
        return pc
    
    def plot(self, ax: Axes | None = None, figsize: tuple = (8, 6)) -> tuple[Figure, Axes]:
        """
        Plot the :class:`PowerCurve` power and thrust curves onto a Matplotlib Figure.

        Parameters
        ----------
        ax : Axes | None, optional
            matplotlib.axes.Axes, by default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6).

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.

        Raises
        ------
        TypeError
            The passed ax argument was not of type Axes.
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
        
        ax.set_title(f'Power Curve: {self.attributes["curve_id"]}', fontweight="bold")
        
        ax.grid(True)
        ax.plot(self.data.index, self.data["power"], color="blue")
        ax.set_xlabel("Wind Speed [m/s]", fontweight="bold")
        ax.set_ylabel("Active Power [kW]", fontweight="bold")
        
        ax.spines['left'].set_color('blue')
        ax.yaxis.label.set_color('blue')
        ax.tick_params(axis='y', colors='blue')
        
        ax2 = ax.twinx()
        ax2.plot(self.data.index, self.data["ct"], color="red")
        ax2.set_ylabel("Thrust Coefficient", fontweight="bold")
        ax2.spines['right'].set_color('red')
        ax2.spines['left'].set_color('blue')
        ax2.yaxis.label.set_color('red')
        ax2.tick_params(axis='y', colors='red')
        
        return fig, ax
    
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
    
    @classmethod
    def from_wtg(cls, filepath: str) -> "PowerCurve":
        """
        Parses a WAsP .wtg file into a :class:`PowerCurve` instance.

        Parameters
        ----------
        filepath : str
            Filepath of the WAsP .wtg to parse.

        Returns
        -------
        PowerCurve
            PowerCurve instance from a WAsP .wtg file.
        """
        
        tree = ET.parse(filepath)
        root = tree.getroot()
        curve_id = root.attrib["Description"]
        manufacturer = root.attrib["ManufacturerName"]
        rotor_diameter = float(root.attrib["RotorDiameter"])
        hub_height = float(root[1][0].text)
        air_density = float(root[2].attrib["AirDensity"])
        curve_date = date.fromisoformat(root[2].attrib["ReleaseDate"])
        revision = int(root[2].attrib["DataStatus"])
        ref_document_id = root[2].attrib["DataSource"]
        cut_in = float(root[2][0].attrib["LowSpeedCutIn"])
        cut_out = float(root[2][0].attrib["HighSpeedCutOut"])
        
        wind_speed = []
        power = []
        ct = []
        for data_point in root[2][2]:
            wind_speed.append(float(data_point.attrib["WindSpeed"]))
            power.append(float(data_point.attrib["PowerOutput"]) / 1000)
            ct.append(float(data_point.attrib["ThrustCoEfficient"]))
        
        rated_ws = wind_speed[wind_speed.index(next(filter(lambda x: x!=0, wind_speed)))]
        
        pc = cls(wind_speed=wind_speed,
                 power=power,
                 ct=ct,
                 rotor_diameter=rotor_diameter,
                 hub_height=hub_height,
                 air_density=air_density,
                 cut_in=cut_in,
                 rated_ws=rated_ws,
                 cut_out=cut_out,
                 curve_id=curve_id,
                 manufacturer=manufacturer
                 )
        
        pc.attributes["curve_date"] = curve_date
        pc.attributes["ref_document_id"] = ref_document_id
        pc.attributes["revision"] = revision

        return pc
    
    @classmethod
    def compare_curves_from_database(cls, curve_ids: Iterable[str], engine: Optional[sa.engine.Engine]=None, username: Optional[str]=None, figsize: tuple=(18, 8)) -> tuple[list["PowerCurve"], Figure, Axes]:
        """
        Class method which returns and compares power curves by passing an iterable of `curve_id` to
        the OW corporate power curve database.
        
        Parameters
        ----------
        curve_ids : Iterable[str]
            Iterable of curve_id strings. Invalid curve_id strings that do not return a power curve
            from the database are ignored.
        engine : Optional[sqlalchemy.engine.Engine], optional
            SqlAlchemy connection engine to the OW corporate power curve database. If no open
            connection is passed, the static method :meth:`PowerCurve.database_engine()` will be
            called implicitly using `username` if passed, otherwise the username is derived from
            the Windows user login name. By default None.
        username : Optional[str], optional
            OW corporate username in the @oceanwinds.com domain. If no username is provided, the
            username is derived friom the Windows user login name. This argument is only required
            if no engine is passed. By default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (18, 8).

        Returns
        -------
        tuple[list[PowerCurve], matplotlib.figure.Figure, numpy.ndarray[matplotlib.axes.Axes]]
            Tuple made up of a list of :class:`PowerCurve` instances, and a Figure and Axes objects.
        """

        if engine is None:
            engine = cls.database_engine(username=username)
        
        power_curves = []
        
        for curve_id in curve_ids:
            try:
                power_curves.append(cls.from_database(curve_id=curve_id, engine=engine))
            except ValueError:
                continue
        
        fig, axs = plt.subplots(ncols=2, figsize=figsize)
        fig.suptitle('Power Curve Comparison', fontsize=16, fontweight="bold")
        
        axs[0].grid(True)
        axs[0].set_title("Power Curves", fontweight="bold")
        axs[0].set_xlabel("Wind Speed [m/s]", fontweight="bold")
        axs[0].set_ylabel("Active Power [kW]", fontweight="bold")
        axs[0].spines['left'].set_color('blue')
        axs[0].yaxis.label.set_color('blue')
        axs[0].tick_params(axis='y', colors='blue')
        
        
        axs[1].grid(True)
        axs[1].set_title("Thrust Curves", fontweight="bold")
        axs[1].set_xlabel("Wind Speed [m/s]", fontweight="bold")
        axs[1].set_ylabel("Thrust Coefficient", fontweight="bold")
        axs[1].spines['left'].set_color('red')
        axs[1].yaxis.label.set_color('red')
        axs[1].tick_params(axis='y', colors='red')
        
        for power_curve in power_curves:
            axs[0].plot(power_curve.data.index, power_curve.data["power"], label=power_curve.attributes["curve_id"])
            axs[1].plot(power_curve.data.index, power_curve.data["ct"], label=power_curve.attributes["curve_id"])
        
        axs[0].legend()
        axs[1].legend()

        return power_curves, fig, axs

    @classmethod
    def get_commercial_turbines(cls, engine: Optional[sa.engine.Engine]=None, username: Optional[str]=None) -> pd.DataFrame:
        """
        Class method that queries the OW corporate power curve database for commercial turbines and return a pandas DataFrame
        of all the available turbines and their attributes.

        Parameters
        ----------
        engine : Optional[sqlalchemy.engine.Engine], optional
            SqlAlchemy connection engine to the OW corporate power curve database. If no open
            connection is passed, the static method :meth:`PowerCurve.database_engine()` will be
            called implicitly using `username` if passed, otherwise the username is derived from
            the Windows user login name. By default None.
        username : Optional[str], optional
            OW corporate username in the @oceanwinds.com domain. If no username is provided, the
            username is derived friom the Windows user login name. This argument is only required
            if no engine is passed. By default None.

        Returns
        -------
        pd.DataFrame
            DataFrame of the available commercial turbines and their attributes.
        """
        if engine is None:
            engine = cls.database_engine(username=username)
        
        turbine_df= pd.read_sql(f"SELECT * FROM power_curves.commercial_turbines ORDER BY id ASC", engine)
        
        return turbine_df
    
    @classmethod
    def get_internal_turbines(cls, engine: Optional[sa.engine.Engine]=None, username: Optional[str]=None) -> pd.DataFrame:
        """
        Class method that queries the OW corporate power curve database for internal turbines and return a pandas DataFrame
        of all the available turbines and their attributes.

        Parameters
        ----------
        engine : Optional[sqlalchemy.engine.Engine], optional
            SqlAlchemy connection engine to the OW corporate power curve database. If no open
            connection is passed, the static method :meth:`PowerCurve.database_engine()` will be
            called implicitly using `username` if passed, otherwise the username is derived from
            the Windows user login name. By default None.
        username : Optional[str], optional
            OW corporate username in the @oceanwinds.com domain. If no username is provided, the
            username is derived friom the Windows user login name. This argument is only required
            if no engine is passed. By default None.

        Returns
        -------
        pd.DataFrame
            DataFrame of the available commercial turbines and their attributes.
        """
        if engine is None:
            engine = cls.database_engine(username=username)
        
        turbine_df= pd.read_sql(f"SELECT * FROM power_curves.internal_turbines ORDER BY id ASC", engine)
        
        return turbine_df

    def interpolate(self, ws_step: float) -> "PowerCurve":
        """
        Method that returns a copy of the :class:`PowerCurve` with interpolated power and thrust
        curves in wind speed steps of ws_step m/s.
        
        The interpolation of power and thrust values is carried out using a Piecewise Cubic Hermite
        Interpolating Polynomial as implemented in :class:`scipy.interpolate.PchipInterpolator`.

        Parameters
        ----------
        ws_step : float
            New wind speed step in m/s to interpolate the power and thurst curves.

        Returns
        -------
        PowerCurve
            A :class:`PowerCurve` object with interpolated power and thrust curves in steps of `ws_step` m/s.
        
        See Also
        --------
        scipy.interpolate.PchipInterpolator : Scipy PCHIP 1-D monotonic cubic interpolation.
        """
        power_interpolator = PchipInterpolator(self.data.index, self.data["power"])
        ct_interpolator = PchipInterpolator(self.data.index, self.data["ct"])
        
        new_ws_bins = np.arange(0, self.data.index.max() + ws_step, ws_step)
        interpolated_power = power_interpolator(new_ws_bins)
        interpolated_ct = ct_interpolator(new_ws_bins)
        
        interpolated_power[np.where((new_ws_bins<self.attributes["cut_in"]) | (new_ws_bins>self.attributes["cut_out"]))] = 0
        interpolated_ct[np.where((new_ws_bins<self.attributes["cut_in"]) | (new_ws_bins>self.attributes["cut_out"]))] = 0
        
        interpolated_power_curve = self.__class__(wind_speed=new_ws_bins, power=interpolated_power,
                                                  ct=interpolated_ct)
        
        interpolated_power_curve.attributes = self.attributes.copy()
        
        return interpolated_power_curve
