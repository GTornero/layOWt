""" aep module. This modeule contains functions to calculate aep using different methods and inputs.
"""
import numpy as np
from py_wake.wind_turbines import WindTurbine
from scipy.integrate import simps


def weibull_aep(wt: WindTurbine, weibull_a: float, weibull_k: float, units: str = "gwh"):
    
    if (weibull_a <= 0) or (weibull_k <= 0):
        raise ValueError(f"a and k must be positive values. {weibull_a=}, {weibull_k=} passed.")
    
    if units.lower() not in ["gwh", "mwh", "kwh", "wh"]:
        raise ValueError(f'units must be in ["gwh", "mwh", "kwh", "wh"]. {units=} passed.')
    
    power_scale_factor = {"gwh": 1e9, "mwh": 1e6, "kwh": 1e3, "wh": 1}[units.lower()]
    
    wind_speed_bins = wt.wt_data[0]["WindSpeed"]
    power = wt.wt_data[0]["PowerOutput"]
    weibull_dist = weibull_k/weibull_a*(wind_speed_bins/weibull_a)**(weibull_k-1) * np.exp(-(wind_speed_bins/weibull_a)**weibull_k)
    aep = 8766*simps(weibull_dist*power, wind_speed_bins)
    
    return aep / power_scale_factor
