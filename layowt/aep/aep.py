""" aep module. This module contains functions to calculate aep using different methods and inputs.
"""
from collections.abc import Iterable

import numpy as np
import rasterio
from layowt.wtg.power_curve import PowerCurve
from rasterio.errors import CRSError, DatasetIOShapeError
from scipy.integrate import simpson
from scipy.special import gamma


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
