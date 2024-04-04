import numpy as np


def saturation_vapor_pressure(temp):
    """
    Calculates the saturation vapor (partial) pressure from temperature in degrees.

    Parameters
    ----------
    temp : float
        Temperature in degrees Celcius.

    Returns
    -------
        Saturation water vapor (partial) pressure in Hectopascals [hPa].
    """
    return 6.112 * np.exp((17.67 * temp) / (temp + 243.5))

def vapor_pressure(rel_humidity, svp):
    """
    Calculates the vapor pressure from the relative humidity and saturation vapor pressure.

    Parameters
    ----------
    rel_humidity : float
        Relative humidty as a ratio (eg. 0.8)
    svp : float
        Saturation vapor (partial) pressure in Hectopascals [hPa].

    Returns
    -------
    float
        Vapor pressure in hectopascals [hPa].
    """
    return rel_humidity * svp

def dewpoint(vapor_pressure):
    """
    Calculates the dewpoint temperature in degrees celcius from the vapor pressure.

    Parameters
    ----------
    vapor_pressure : float
        Vapor pressure in Hectopascals [hPa].

    Returns
    -------
    float
        Dewpoint temperature in degrees Celcius.
    """
    val = np.log(vapor_pressure / 6.112)
    return (243.5 * val) / (17.67-val)

def dewpoint_from_rh(temp, rh):
    """
    Calculates the dewpoint temperature in degrees Celcius from air temperature and relative humidity.
    

    Parameters
    ----------
    temp : float
        Air temperature in degrees Celcius.
    rh : float
        Relative humidity expressed as a ratio in the range 0 < relative_humidity <= 1.

    Returns
    -------
    float
        Dewpoint temperature in degrees Celcius.
    """
    return dewpoint(vapor_pressure(rh, saturation_vapor_pressure(temp)))

