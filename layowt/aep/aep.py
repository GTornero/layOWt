""" aep module. This module contains functions to calculate aep using different methods and inputs.
"""
from collections.abc import Iterable

import numpy as np
import pandas as pd
import rasterio
from rasterio.errors import CRSError, DatasetIOShapeError
from scipy.integrate import simpson
from scipy.special import gamma
from scipy.stats import linregress

from layowt.stats.stats import weibull_mean
from layowt.wtg.power_curve import PowerCurve


def _weibull_dist(wind_speed: float | Iterable[float], weibull_a: float, weibull_k: float) -> float | Iterable[float]:
    return weibull_k / weibull_a * (wind_speed / weibull_a)**(weibull_k-1) * np.exp(-(wind_speed / weibull_a)**weibull_k)

def weibull_aep(pc: PowerCurve, weibull_a: float, weibull_k: float, units: str = "GWh", leap_year: bool=True) -> float:
    """
    This function calculates the gross AEP using a LayOWt PowerCurve and a Weibull distribution. \
    The method integrates the result of the dicrete function resulting from the multiplication of \
    the number of hours in a year, the PowerCurve power curve, and the 2-parameter Weibull \
    probability density function. The integration is carried out using Simpon's rule as implemented \
    in :func:`scipy.integrate.simpson`. The calculted AEP can be reported in Wh, kWh, MWh, and GWh. \
    

    .. math::
        aep = \\displaystyle [hours\_in\_year]\\int_{cut-in}^{cut-out} P(ws)PC(ws) dws

    Where the Weibull distribution Probability Density Function is given by:

    .. math::
        P(ws) = \\frac{k}{A}\\left(\\frac{ws}{A}\\right)^{k-1}e^{-(\\frac{ws}{A})^k}

    Parameters
    ----------
    pc : PowerCurve
        LayOWt PowerCurve object. This object contains all of the turbine specifications such as power and Ct vs wind speed curves, rotor size, etc.
    weibull_a : float
        The 2-parameter Weibull scale parameter in m/s.
    weibull_k : float
        The 2-parameter Weibull shape parameter.
    units : str, optional
        Units in which to report the calculated AEP. Must be one of "Wh", "kWh", "MWh", or "GWh" (case insensitive), by default "GWh".
    leap_year : bool, Optional
        If True, the number of hours in a year used in the calculation is 8766, otherwise 8760, by default True.

    Returns
    -------
    float
        The gross AEP for the passed PowerCurve assuming a theoretical 2-parameter Weibull distribution.

    Raises
    ------
    ValueError
        If weibull_a or weibull_k are non positive.
    ValueError
        If units is not in ["gwh", "mwh", "kwh", "wh"].

    Examples
    --------
    >>> from layowt.wtg import PowerCurve
    >>> from layowt.aep import weibull_aep
    >>> power_curve = PowerCurve.from_wtg("prototype_inputs/NGT260-21MW_Vestas_v3.wtg")
    >>> print(weibull_aep(power_curve, weibull_a=11, weibull_k=2))
    99.66977656017633

    See Also
    --------
    layowt.wtg.power_curve.PowerCurve : LayOWt PowerCurve object.

    scipy.integrate.simpson : Simpson's rule implemented by Scipy.
    """
    if (weibull_a <= 0) or (weibull_k <= 0):
        raise ValueError(f"a and k must be positive values. {weibull_a=}, {weibull_k=} passed.")

    if units.lower() not in ["gwh", "mwh", "kwh", "wh"]:
        raise ValueError(f'units must be in ["gwh", "mwh", "kwh", "wh"]. {units=} passed.')
    
    if leap_year:
        hours_in_year = 8766
    else:
        hours_in_year = 8760

    power_scale_factor = {"gwh": 1e6, "mwh": 1e3, "kwh": 1, "wh": 1e-3}[units.lower()]

    wind_speed_bins = pc.data.index
    power = pc.data["power"]
    weibull_dist = _weibull_dist(wind_speed_bins, weibull_a, weibull_k)
    aep = hours_in_year*simpson(weibull_dist*power, wind_speed_bins)

    return aep / power_scale_factor

def weibull_wpd(weibull_a: float, weibull_k: float, air_density: float = 1.225) -> float:
    """
    This function calculates the mean wind power density from a 2-parameter Weibull distrbution.
    The mean wind power density of a weibull distribution, in W/m\\ :sup:`2`\\  is defined as:
    
    .. math::
        wpd = \\frac{1}{2} \\rho A^3 \\Gamma(1 + \\frac{3}{k})
    
    Where:
        * :math:`\\rho = \\textrm{Air density}`
        * :math:`A = \\textrm{Weibull Scale parameter}`
        * :math:`k = \\textrm{Weibull Shape parameter}`
        * :math:`\\Gamma = \\textrm{Gamma function}`

    Parameters
    ----------
    weibull_a : float
        The 2-parameter Weibull scale parameter in m/s.
    weibull_k : float
        The 2-parameter Weibull shape parameter.
    air_density : float, optional
        Density of air in kg/m\\ :sup:`3`\\, by default 1.225.

    Return
    -------
    float
        Mean wind power density of the Weibull distribution in W/m\\ :sup:`2`\\.

    Raises
    ------
    ValueError
        If weibull_a or weibull_k are not positive values.
    
    See Also
    --------
    scipy.special.gamma : Gamma function implemented by Scipy.
    """
    if weibull_a <= 0 or weibull_k <= 0:
        raise ValueError(f"Invalid Weibull parameters passed. A and k must be positive values. {weibull_a=} {weibull_k=} passed.")
    
    wpd = 0.5 * air_density * (weibull_a**3) * gamma(1 + 3/weibull_k)
    return wpd

def weibull_aep_map(pc: PowerCurve, weibull_a: str, weibull_k: str, output_path: str = "aep.tif", units: str = "GWh", leap_year: bool=True) -> None:
    """
    This function takes in a PowerCurve and Weibull A (scale) and k (shape) parameter raster \
    datasets, and creates a gross annual energy production raster dataset. The :func:`weibull_aep` \
    function is used to calculate gross annual energy production. \
    
    The CRS, shape, and cell size of the Weibull A and k datasets must be the same. I.e. the rasters
    must have the same extent and cells must be aligned prior to being used with this function.

    Parameters
    ----------
    pc : PowerCurve
        LayOWt PowerCurve object. This object contains all of the turbine specifications such as power and Ct vs wind speed curves, rotor size, etc.
    weibull_a : str
        The filepath string for the 2-parameter Weibull scale parameter raster dataset in m/s.
    weibull_k : str
        The filepath string for the 2-parameter Weibull shape parameter raster dataset.
    output_path : str, optional
        Filepath of the output GAEP raster dataset, by default "aep.tif"
    units : str, optional
        Units in which to report the calculated AEP. Must be one of "Wh", "kWh", "MWh", or "GWh" (case insensitive), by default "GWh".
    leap_year : bool, optional
        If True, the number of hours in a year used in the calculation is 8766, otherwise 8760, by default True.

    Raises
    ------
    CRSError
        CRS of the Weibull A and k datasets is not the same.
    DatasetIOShapeError
        Size of the Weibull A and k datasets is not the same.
    ValueError
        Extant or cell size of the Weibull A and k datasets is not the same.
    
    See Also
    --------
    weibull_aep : Calculate Gross Annual Energy Production from a 2-Parameter Weibull Distribution.
    """
    weibull_aep_vector = np.vectorize(weibull_aep)
    
    with rasterio.open(weibull_a) as src_a:
        with rasterio.open(weibull_k) as src_k:
            if src_a.crs != src_k.crs:
                raise CRSError(f"A and k datasets must use the same CRS. A: {src_a.crs}, k: {src_k.crs}.")
            elif src_a.shape != src_k.shape:
                raise DatasetIOShapeError(f"A and k datasets must be of the same size. A: {src_a.shape}, k: {src_k.shape}.")
            elif src_a.bounds != src_k.bounds:
                raise ValueError("A and k datasets must have the same extent.")
            elif src_a.transform[0:5:4] != src_k.transform[0:5:4]:
                raise ValueError(f"A and k datasets must have the same cell size. A: {src_a.transform[0:5:4]}, k: {src_k.transform[0:5:4]}.")
            
            # clean A and k values below zero
            a_array = src_a.read(1)
            k_array = src_k.read(1)
            
            a_array[a_array <= 0] = np.nan
            k_array[k_array <= 0] = np.nan
            
            result = weibull_aep_vector(pc, a_array, k_array, units=units, leap_year=leap_year)

            with rasterio.open(
                output_path,
                'w',
                driver='GTiff',
                height=src_a.shape[0],
                width=src_a.shape[1],
                count=1,
                dtype=result.dtype,
                crs=src_a.crs,
                transform=src_a.transform
                ) as output:
                    output.write(result, 1)

def weibull_mean_map(weibull_a: str, weibull_k: str, output_path: str = "mean.tif") -> None:
    """
    This function takes in a Weibull A (scale) and k (shape) parameter raster \
    datasets, and creates a Weibull mean value raster dataset. The :func:`layowt.stats.stats.weibull_mean` \
    function is used to calculate the mean value of the 2-Parameter Weibull Distribution. \
    
    The CRS, shape, and cell size of the Weibull A and k datasets must be the same. I.e. the rasters
    must have the same extent and cells must be aligned prior to being used with this function.

    Parameters
    ----------
    weibull_a : str
        The filepath string for the 2-parameter Weibull scale parameter raster dataset.
    weibull_k : str
        The filepath string for the 2-parameter Weibull shape parameter raster dataset.
    output_path : str, optional
        Filepath of the output GAEP raster dataset, by default "mean.tif"

    Raises
    ------
    CRSError
        CRS of the Weibull A and k datasets is not the same.
    DatasetIOShapeError
        Size of the Weibull A and k datasets is not the same.
    ValueError
        Extant or cell size of the Weibull A and k datasets is not the same.
    
    See Also
    --------
    layowt.stats.stats.weibull_mean : Calculates the mean value of a 2-Parameter Weibull Distribution.
    """
    weibull_mean_vector = np.vectorize(weibull_mean)
    
    with rasterio.open(weibull_a) as src_a:
        with rasterio.open(weibull_k) as src_k:
            if src_a.crs != src_k.crs:
                raise CRSError(f"A and k datasets must use the same CRS. A: {src_a.crs}, k: {src_k.crs}.")
            elif src_a.shape != src_k.shape:
                raise DatasetIOShapeError(f"A and k datasets must be of the same size. A: {src_a.shape}, k: {src_k.shape}.")
            elif src_a.bounds != src_k.bounds:
                raise ValueError("A and k datasets must have the same extent.")
            elif src_a.transform[0:5:4] != src_k.transform[0:5:4]:
                raise ValueError(f"A and k datasets must have the same cell size. A: {src_a.transform[0:5:4]}, k: {src_k.transform[0:5:4]}.")
            
            # clean A and k values below zero
            a_array = src_a.read(1)
            k_array = src_k.read(1)
            
            a_array[a_array <= 0] = np.nan
            k_array[k_array <= 0] = np.nan
            
            result = weibull_mean_vector(a_array, k_array)
            
            result[result <= 0] = np.nan

            with rasterio.open(
                output_path,
                'w',
                driver='GTiff',
                height=src_a.shape[0],
                width=src_a.shape[1],
                count=1,
                dtype=result.dtype,
                crs=src_a.crs,
                transform=src_a.transform
                ) as output:
                    output.write(result, 1)

def weibull_ws_to_aep(pc: PowerCurve, weibull_a: float, weibull_k: float, ws_range: float = 5, ws_step: float = 0.1, how: str | None = None) -> float:
    """
    Function to calculate the approximate linear relationship between % changes in wind speed and \
    % changes in annual energy production for a given PowerCurve and mean 2-parameter Weibull wind \
    distribution.
    
    Annual energy production is calculated for a range of 2-parameter Weibull distributions with
    small % variations to their scale parameter A between **-ws_range** and **ws_range** in steps of
    **ws_step**. Annual energy production is calculated using the :func:`weibull_aep` function.
    
    A linear relationship is fit to the results of % change in wind speed vs % change in AEP using
    Scipy's implementation of linear regression using the :func:`scipy.stats.linregress` function.
    
    This method is applicable for small ranges of wind speed variations. For large deviations, the
    real effect on annual energy production deviates from a linear relationship.

    Parameters
    ----------
    pc : PowerCurve
        PowerCurve used to calculate AEP
    weibull_a : float
        The 2-parameter Weibull scale parameter in m/s.
    weibull_k : float
        The 2-parameter Weibull shape parameter.
    ws_range : float, optional
        Range of wind speed % deviation to consider in the analysis, by default 5.
    ws_step : float, optional
        wind speed % step increase within the range of % wind speed deviations, by default 0.1.
    how : str | None
        Determines the range of wind speed % deviations is calculated. If None or 'negative' is
        passed, only negative wind speed deviations are used in the anlysis, from -ws_range to 0.
        If 'positive', from 0 to +ws_range. If 'both' then the double-sided range is used, from
        -ws_range to +ws_range. By default None.

    Returns
    -------
    float
        The slope of the approximate linear relationship between % changes in wind speed and % changes in AEP.
    
    Raises
    ------
    ValueError
        If how not one of 'negative', 'positive', or 'both'.
    """
    if how is not None:
        if how.lower() not in ['negative', 'positive', 'both']:
            raise ValueError(f"Invalid value passed for {how=}. Must be one of 'negative', 'positive', or 'both'.")
    
    if (how is None) or (how.lower() == 'negative'):
        ws_deltas = 1 + (np.arange(-ws_range, ws_step, ws_step) / 100)
    elif how.lower() == "both":
        ws_deltas = 1 + (np.arange(-ws_range, ws_range + ws_step, ws_step) / 100)
    elif how.lower() == 'positive':
        ws_deltas = 1 + (np.arange(0, ws_range + ws_step, ws_step) / 100)
    
    aeps = np.empty(ws_deltas.shape)
    
    baseline = weibull_aep(pc=pc, weibull_a=weibull_a, weibull_k=weibull_k)
    
    for i, delta in enumerate(ws_deltas):
        aeps[i] = weibull_aep(pc=pc, weibull_a=weibull_a * delta, weibull_k=weibull_k)
    
    aep_deltas = aeps / baseline - 1
    
    r = linregress(ws_deltas, aep_deltas)
    
    return r.slope

def ts_ws_to_aep(ts: Iterable[float], pc: PowerCurve, ws_range: float = 5, ws_step: float = 0.1, how: str | None = None) -> float:
    """
    Function to calculate the approximate linear relationship between % changes in wind speed and \
    % changes in annual energy production for a given PowerCurve and wind speed timeseries.
    
    Annual energy production is calculated for a timeseries of wind speeds with small % variations
    in wind speed between **-ws_range** and **ws_range** in steps of **ws_step**. Annual energy
    production is calculated by looking up the wind speed values against the power curve.
    
    A linear relationship is fit to the results of % change in wind speed vs % change in AEP using
    Scipy's implementation of linear regression using the :func:`scipy.stats.linregress` function.
    
    This method is applicable for small ranges of wind speed variations. For large deviations, the
    real effect on annual energy production deviates from a linear relationship.

    Parameters
    ----------
    ts : Iterable
        Iterable containing wind speeds.
    pc : PowerCurve
        PowerCurve used to calculate AEP.
    ws_range : float, optional
        Range of wind speed % deviation to consider in the analysis, by default 5.
    ws_step : float, optional
        wind speed % step increase within the range of % wind speed deviations, by default 0.1.
    how : str | None
        Determines the range of wind speed % deviations is calculated. If None or 'negative' is
        passed, only negative wind speed deviations are used in the anlysis, from -ws_range to 0.
        If 'positive', from 0 to +ws_range. If 'both' then the double-sided range is used, from
        -ws_range to +ws_range. By default None.

    Returns
    -------
    float
        The slope of the approximate linear relationship between % changes in wind speed and % changes in AEP.
    
    Raises
    ------
    ValueError
        If how not one of 'negative', 'positive', or 'both'.
    """
    def _custom_round(x, base=0.5):
        return round(x*1/base)*base

    if how is not None:
        if how.lower() not in ['negative', 'positive', 'both']:
            raise ValueError(f"Invalid value passed for {how=}. Must be one of 'negative', 'positive', or 'both'.")
    
    if not isinstance(ts, pd.Series):
        ws_df = pd.Series(ts) 
    else:
        ws_df = ts.copy(deep=True)
    
    if (how is None) or (how.lower() == 'negative'):
        ws_deltas = 1 + (np.arange(-ws_range, ws_step, ws_step) / 100)
    elif how.lower() == "both":
        ws_deltas = 1 + (np.arange(-ws_range, ws_range + ws_step, ws_step) / 100)
    elif how.lower() == 'positive':
        ws_deltas = 1 + (np.arange(0, ws_range + ws_step, ws_step) / 100)
    
    aeps = np.empty(ws_deltas.shape)
    
    base = round(pd.Series(pc.data.index).diff(periods=1).mode().values[0], 3)
    
    round_ts = _custom_round(ws_df, base)
    
    baseline = pc.data['power'].loc[round_ts[(round_ts>=pc.attributes['cut_in']) & (round_ts <= pc.attributes['cut_out'])]].sum()
    
    for i, delta in enumerate(ws_deltas):
        delta_ts = _custom_round(round_ts * delta, base)
        aeps[i] =  pc.data['power'].loc[delta_ts[(delta_ts>=pc.attributes['cut_in']) & (delta_ts <= pc.attributes['cut_out'])]].sum()
    
    aep_deltas = aeps / baseline - 1
    
    r = linregress(ws_deltas, aep_deltas)

    return r.slope