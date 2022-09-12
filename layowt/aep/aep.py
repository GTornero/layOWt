""" aep module. This module contains functions to calculate aep using different methods and inputs.
"""
import numpy as np
from layowt.wtg.power_curve import PowerCurve
from scipy.integrate import simpson
from scipy.special import gamma


def weibull_aep(pc: PowerCurve, weibull_a: float, weibull_k: float, units: str = "GWh", leap_year: bool=True) -> float:
    """
    This function calculates the gross AEP using a LayOWt PowerCurve and a Weibull distribution.
    The method integrates the result of the dicrete function resulting from the multiplication of
    the number of hours in a year, the PowerCurve power curve, and the 2-parameter Weibull
    probability density function. The integration is carried out using Simpon's rule as implemented
    in :func:`scipy.integrate.simpson`. The calculted AEP can be reported in Wh, kWh, MWh, and GWh.
    

    .. math::
        aep = \\displaystyle [hours\_in\_year]\\int_{cut-in}^{cut-out} P(ws)PC(ws) dws

    Where the Weibull distribution Probability Density Function is given by:

    .. math::
        P(ws) = \\frac{k}{A}\\left(\\frac{ws}{A}\\right)^{k-1}e^{-(\\frac{ws}{A})^k}

    Parameters
    ----------
    wt : PowerCurve
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
    weibull_dist = weibull_k / weibull_a * (wind_speed_bins / weibull_a)**(weibull_k-1) * np.exp(-(wind_speed_bins/weibull_a)**weibull_k)
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

    Returns
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
