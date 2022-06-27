""" aep module. This modeule contains functions to calculate aep using different methods and inputs.
"""
import numpy as np
from py_wake.wind_turbines import WindTurbine
from scipy.integrate import simpson


def weibull_aep(wt: WindTurbine, weibull_a: float, weibull_k: float, units: str = "GWh"):
    """weibull_aep This function calculates the gross AEP using a py_wake WindTurbine and a Weibull distribution. The method integrates the result of dicrete function resulting from the multiplication of the WindTurbine power curve and the 2-parameter Weibull probability density function. The integration is carried out using Simpon's rule as implemented in Scipy.integrate.simpson. The calculted AEP can be reported in Wh, kWh, MWh, and GWh.

    Parameters
    ----------
    wt : WindTurbine
        py_wake WindTurbine object. This object contains all of the turbine specifications such as power and Ct vs wind speed curves, rotor size, etc.
    weibull_a : float
        The 2-parameter Weibull scale parameter. The
    weibull_k : float
        The 2-parameter Weibull shape parameter.
    units : str, optional
        Units in which to report the calculated AEP. Must be one of "Wh", "kWh", "MWh", or "GWh" (case insensitive), by default "GWh".

    Returns
    -------
    float
        The gross AEP for the passed WindTurbine assuming a theoretical 2-parameter Weibull distribution.

    Raises
    ------
    ValueError
        If weibull_a or weibull_k are non positive.
    ValueError
        If units is not in ["gwh", "mwh", "kwh", "wh"].
    
    Examples
    --------
    >>> from layowt.grids import Grid
    >>> from layowt.layouts import Layout
    >>> layout = Layout(Grid())
    >>> layout.load_wtg("prototype_inputs/NGT260-21MW_Vestas_v3.wtg")
    >>> weibull_aep(layout.wtg, weibull_a=11, weibull_k=2)
    99.66977656017633
    
    See Also
    --------
    layowt.layouts.layout.Layout.load_wtg : Loads a WAsP .wtg file into the layout as a py_wake WindTurbine object.
    
    py_wake.wind_turbines.WindTurbine : py_wake WindTurbine object.
    
    scipy.integrate.simpson : Simpson's rule implemented by Scipy.
    
    """
    if (weibull_a <= 0) or (weibull_k <= 0):
        raise ValueError(f"a and k must be positive values. {weibull_a=}, {weibull_k=} passed.")
    
    if units.lower() not in ["gwh", "mwh", "kwh", "wh"]:
        raise ValueError(f'units must be in ["gwh", "mwh", "kwh", "wh"]. {units=} passed.')
    
    power_scale_factor = {"gwh": 1e9, "mwh": 1e6, "kwh": 1e3, "wh": 1}[units.lower()]
    
    wind_speed_bins = wt.wt_data[0]["WindSpeed"]
    power = wt.wt_data[0]["PowerOutput"]
    weibull_dist = weibull_k/weibull_a*(wind_speed_bins/weibull_a)**(weibull_k-1) * np.exp(-(wind_speed_bins/weibull_a)**weibull_k)
    aep = 8766*simpson(weibull_dist*power, wind_speed_bins)
    
    return aep / power_scale_factor
