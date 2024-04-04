import numpy as np


def saturation_vapor_pressure(temp):
    """
    Calculates the saturation vapor (partial) pressure from temperature in degrees.
    
    Saturation vapor pressure is given by the formula from Bolton 1980 for temperature T in degrees Celcius:
    
    .. math:: svp = 6.112 e^\\frac{17.67T}{T + 243.5}

    Parameters
    ----------
    temp : float
        Temperature in degrees Celcius.

    Returns
    -------
        Saturation vapor (partial) pressure in Hectopascals [hPa].
    
    See Also
    --------
    layowt.meteo.thermo.vapor_pressure : Calculates the vapor pressure from the relative humidity and saturation vapor pressure.
    layowt.meteo.thermo.dewpoint : Calculates the dewpoint temperature in degrees celcius from the vapor pressure.
    
    References
    ----------
    .. [1] [Bolton1980] Bolton, D., 1980: The computation of equivalent potential temperature. Mon. Wea. Rev., 108, 1046-1053, doi:10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2.
        <https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml>
    """
    return 6.112 * np.exp((17.67 * temp) / (temp + 243.5))

def vapor_pressure(rel_humidity, svp):
    """
    Calculates the vapor (partial) pressure from the relative humidity and saturation vapor pressure.

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
    Calculates the ambient dewpoint temperature in degrees celcius from the vapor pressure.
    
    The dewpoint temperature is calculated by inverting the Bolton 1980 formula for saturation
    vapor pressure to calculate temperature:
    
    .. math:: T = \\frac{243.5 \\log(e / 6.112)}{17.67 - \\log(e / 6.112)}

    Parameters
    ----------
    vapor_pressure : float
        Vapor pressure in Hectopascals [hPa].

    Returns
    -------
    float
        Dewpoint temperature in degrees Celcius.
    
    See Also
    --------
    layowt.meteo.thermo.saturation_vapor_pressure : Calculates the saturation vapor (partial) pressure from temperature in degrees.
    layowt.meteo.thermo.vapor_pressure : Calculates the vapor pressure from the relative humidity and saturation vapor pressure.
    layowt.meteo.thermo.dewpoint_from_rh : Calculates the dewpoint temperature in degrees Celcius from air temperature and relative humidity.
    
    References
    ----------
    .. [1] [Bolton1980] Bolton, D., 1980: The computation of equivalent potential temperature. Mon. Wea. Rev., 108, 1046-1053, doi:10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2.
        <https://journals.ametsoc.org/view/journals/mwre/108/7/1520-0493_1980_108_1046_tcoept_2_0_co_2.xml>
    """
    val = np.log(vapor_pressure / 6.112)
    return (243.5 * val) / (17.67-val)

def dewpoint_from_rh(temp, rh):
    """
    Calculates the ambient dewpoint temperature in degrees Celcius from air temperature and relative humidity.
    

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
    
    See Also
    --------
    layowt.meteo.thermo.saturation_vapor_pressure : Calculates the saturation vapor (partial) pressure from temperature in degrees.
    layowt.meteo.thermo.vapor_pressure : Calculates the vapor pressure from the relative humidity and saturation vapor pressure.
    layowt.meteo.thermo.dewpoint : Calculates the ambient dewpoint temperature in degrees celcius from the vapor pressure.
    """
    return dewpoint(vapor_pressure(rh, saturation_vapor_pressure(temp)))

