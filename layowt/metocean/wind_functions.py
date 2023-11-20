from typing import Iterable

import numpy as np
from numpy.typing import NDArray
from scipy.special import gamma


def vref_ewts_exact(v: float, k: float, T: float=50, n: int=23037) -> float:
    """
    EWTS II Exact algorithm to estimate Vref from the mean wind speed and the Weibull k value.
    
    The algorithm appeared in the European Wind Turbine Standards II (EWTS II) Exact variation.
    
    The theoretical basis of the algorithm is that if a set of wind speeds conform to a Weibull
    distribution, then the Weibull k factor will determine the ratio of extreme wind speeds to the
    average wind speed.
    
    The algorithm uses the following equation:
    
    .. math::
        V_{ref} = \\frac{V_{ave}}{\\Gamma\\left(1 + 1/k\\right)}\\left[-\\ln\\left\\{1-\\exp\\left(\\frac{\\ln\\left(1 - 1/T\\right)}{n}\\right)\\right\\}\\right]^{1/k}
    
    Where:
        * :math:`V_{ref}` is the extreme wind speed
        * :math:`V_{ave}` is the long-term mean wind speed
        * :math:`\\Gamma` is the gamma function
        * :math:`k` is the Weibull k (shape) parameter
        * :math:`T` is the return period in years
        * :math:`n` is the number of independent events per year (23,037 for 10-min time steps in 1 year)

    Parameters
    ----------
    v : float
        Long term mean wind speed in m/s.
    k : float
        Weibull k (shape) parameter.
    T : float, optional
        Return period in years, by default 50.
    n : int, optional
        Number of independent events per year, by default 23,0037 for 10-min time steps. An
        explanation of this parameter is provided in <https://c2wind.com/f/content/ewts2new.pdf>.

    Returns
    -------
    float
        Extreme wind speed at the defined return period.
    
    References
    ----------
    .. [1] Dekker JWM, Pierik JTG (1998): 'European Wind Turbine Standards II', ECN-C-99-073, ECN Solar & Wind Energy, Netherlands.
        <https://c2wind.com/f/content/ewts2new.pdf>
    """
    a = 1. / gamma(1. + 1. / k)
    b = - np.log(1 - np.exp((np.log(1. - 1. / T)) / n))
    v_ref = v * a * (b ** (1. / k))
    
    return v_ref

def vref_ewts_gumbel(v: float, k: float, T: float=50, n: int=23037) -> float:
    """
    EWTS II Gumbel algorithm to estimate Vref from the mean wind speed and the Weibull k value.
    
    The algorithm appeared in the European Wind Turbine Standards II (EWTS II) Exact variation.
    
    The theoretical basis of the algorithm is that if a set of wind speeds conform to a Weibull
    distribution, then the Weibull k factor will determine the ratio of extreme wind speeds to the
    average wind speed.
    
    The algorithm uses the following equation:
    
    .. math::
        V_{ref} = \\frac{V_{ref}\\left(\\ln n\\right)^{\\frac{1}{k}-1}}{k\\Gamma\\left(1+1/k\\right)}\\left[k\\ln n - \\ln\\left\\{-\\ln\\left(1 - 1/T\\right)\\right\\}\\right]
    
    Where:
        * :math:`V_{ref}` is the extreme wind speed
        * :math:`V_{ave}` is the long-term mean wind speed
        * :math:`\\Gamma` is the gamma function
        * :math:`k` is the Weibull k (shape) parameter
        * :math:`T` is the return period in years
        * :math:`n` is the number of independent events per year (23,037 for 10-min time steps in 1 year)

    Parameters
    ----------
    v : float
        Long term mean wind speed in m/s.
    k : float
        Weibull k (shape) parameter.
    T : float, optional
        Return period in years, by default 50.
    n : int, optional
        Number of independent events per year, by default 23,0037 for 10-min time steps. An
        explanation of this parameter is provided in <https://c2wind.com/f/content/ewts2new.pdf>.

    Returns
    -------
    float
        Extreme wind speed at the defined return period.
    
    References
    ----------
    .. [1] Dekker JWM, Pierik JTG (1998): 'European Wind Turbine Standards II', ECN-C-99-073, ECN Solar & Wind Energy, Netherlands.
        <https://c2wind.com/f/content/ewts2new.pdf>
    """
    
    a = ((np.log(n))**(1/k -1)) / (k*gamma(1+1/k))
    b = k * np.log(n) - np.log(-np.log(1-1/T))
    vref = v * a * b
    return vref

def vref_ewts_davenport(v: float, k: float, T: float=50, n: int=23037) -> float:
    """
    EWTS II Davenport algorithm to estimate Vref from the mean wind speed and the Weibull k value.
    
    The algorithm appeared in the European Wind Turbine Standards II (EWTS II) Exact variation.
    
    The theoretical basis of the algorithm is that if a set of wind speeds conform to a Weibull
    distribution, then the Weibull k factor will determine the ratio of extreme wind speeds to the
    average wind speed.
    
    The algorithm uses the following equation:
    
    .. math::
        V_{ref} = \\frac{V_{ave}\\left(\\ln n \\right)^{\\frac{1}{k}-1}}{c_{1}k\\Gamma\\left(1 + 1/k\\right)}\\left[c_{1}c_{2}k\\ln n - \\ln \\left\\{-\\ln \\left(1 - 1/T\\right)\\right\\}\\right]
        
    Where:
        * :math:`V_{ref}` is the extreme wind speed
        * :math:`V_{ave}` is the long-term mean wind speed
        * :math:`\\Gamma` is the gamma function
        * :math:`k` is the Weibull k (shape) parameter
        * :math:`T` is the return period in years
        * :math:`n` is the number of independent events per year (23,037 for 10-min time steps in 1 year)
    
    And where the constants :math:`c_1` and :math:`c_2` are given by the following two equations:
    
    .. math::
        c_1 = 1 - \\frac{k-1}{k\\ln n}
        
    .. math::
        c_2 = 1 + \\frac{\\ln\\left[k\\Gamma\\left(\\frac{1}{k}+1\\right)\\left(\\ln n\\right)^{\\frac{k-1}{k}}\\right]}{k\\ln n - \\left(k-1\\right)}

    Parameters
    ----------
    v : float
        Long term mean wind speed in m/s.
    k : float
        Weibull k (shape) parameter.
    T : float, optional
        Return period in years, by default 50.
    n : int, optional
        Number of independent events per year, by default 23,0037 for 10-min time steps. An
        explanation of this parameter is provided in <https://c2wind.com/f/content/ewts2new.pdf>.

    Returns
    -------
    float
        Extreme wind speed at the defined return period.
    
    References
    ----------
    .. [1] Dekker JWM, Pierik JTG (1998): 'European Wind Turbine Standards II', ECN-C-99-073, ECN Solar & Wind Energy, Netherlands.
        <https://c2wind.com/f/content/ewts2new.pdf>
    """
    
    c1 = 1 - (k - 1) / (k * np.log(n))
    c2 = 1 + (np.log(k * gamma(1+1/k) * ((np.log(n))**(1-1/k)))) / (k * np.log(n) - k + 1)
    a = ((np.log(n))**(1/k-1)) / (c1 * k * gamma(1+1/k))
    b = c1 * c2 * k * np.log(n) - np.log(-np.log(1-1/T))
    vref = v * a * b
    return vref

def froya(u: float, T: int, T0: int, h_ref: float, z: float) -> float:
    """
    The Frøya wind profile can be used to convert the mean wind speed :math:`U_0` with averaging \
    period :math:`T_0` at height H above sea level to the mean wind speed U with averaging period T at \
    height z above sea level.
    
    The Frøya profile includes a vertical wind speed profile, and a gust factor which allows for
    conversion of mean wind speeds between different averaging periods.
    
    The Frøya profile provides an estimate to the mean wind speed :math:`U(T, z)` with averaging period of T seconds and at height z metres above mean sea level given by:
    
    .. math::
        U\\left(T, z\\right) = U_{0}\\left[1 + C \\ln\\left(\\frac{z}{H}\\right)\\right]\\left[1 - 0.41I_{U}\\left(z\\right)\\ln\\left(\\frac{T}{T_0}\\right)\\right]
    
    Where:
        * :math:`U_0` is the reference mean wind speed in metres per second
        * :math:`T` is the averaging period of :math:`U_0` in seconds
        * :math:`H` is the height of :math:`U_0` in metres above mean sea level
        * :math:`C = 0.0573\\sqrt{1 + 0.1483U_{0}}`
        * :math:`I_{U}\\left(z\\right) = 0.06\\left(1 + 0.043U_{0}\\right)\\left(\\frac{z}{H}\\right)^{-0.22}`
    
    The method is valid for :math:`T < T_0`.

    Parameters
    ----------
    u : float
        Reference mean wind seed with averaging period T0 seconds, and at height z metres above mean sea level.
    T : int
        Target averaging period in seconds. Must be less than or equal to T0.
    T0 : int
        Averaging period of the reference mean wind speed in seconds.
    h_ref : float
        Height of the reference mean wind speed in metres above sea level.
    z : float
        Target height in metre above mean sea level.

    Returns
    -------
    float
        Mean wind speed with averaging period of T seconds at height z metres above sea level from the Frøya wind profile.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    .. [2] Andersen, Odd & Løvseth, Jørgen. (2006). The Frøya database and maritime boundary layer wind description. Marine Structures. 19. 173-192. 10.1016/j.marstruc.2006.07.003.
        <https://www.researchgate.net/publication/229394889_The_Froya_database_and_maritime_boundary_layer_wind_description>
    """
    if not T <= T0:
        raise ValueError("T must smaller or equal to T0.")
    
    c = 0.0573 * np.sqrt(1 + 0.148 * u)
    I_u = 0.06 * (1 + 0.043 * u) * ((z / h_ref)**-0.22)
    
    a = 1 + (c * np.log(z / h_ref))
    b = 1 - (0.41 * I_u * np.log(T/T0))
    
    return u * a * b

def ntm(i_ref: float, ws: Iterable[float] | None=None, a: float = 0.75, b: float = 5.6) -> NDArray[np.float64]:
    """
    Normal Turbulence Model defined by the IEC.
    
    Calculates the Normal Turbulence profile vs wind speed for the standard IEC wind turbine classes.
    
    The normal turbulence model is given by:

    .. math::
        \\sigma_{1} = I_{ref}\\left(aV_{hub}+b\\right)
    
    Where:
        * :math:`\\sigma_{1}` is the reference representative value of the turbulence standard deviation, given by the 90% quantile for hte given hub height wind speed.
        * :math:`V_{hub}` is the hub height wind speed in m/s.
        * :math:`a` is equal to 0.75.
        * :math:`b` is equal to 5.6.

    Parameters
    ----------
    i_ref : float
        The reference value of representative ambient turbulence at 15m/s.
    ws : Iterable[float] | None, optional
        The wind speed values at which to calculate the normal turbulence model profile, by default None.
    a : float, optional
        Standard constant equal to 0.75 defined by the IEC, by default 0.75
    b : float, optional
        Standard constant equal to 5.6 defined by the IEC, by default 5.6.

    Returns
    -------
    NDArray[np.float64]
        Array of ambient turbulence values at the specified wind speeds.
    
    References
    ----------
    .. [1] IEC 61400-1: Wind Energy Generation Systems – Part 1: Design Requirements, edition 4.0, 2019-02.
        <https://webstore.iec.ch/publication/26423>
    """
    if ws is None:
        ws_array = np.arange(0.1, 35, 0.1, dtype=np.float64)
    else:
        ws_array = np.array(ws, dtype=np.float64)
    
    ti = i_ref * ((a * ws_array) + b) / ws_array
    
    return ti

def dnv_wind_profile(u: float, T: int, T0: int, h_ref: float, z: float) -> float:
    """
    The wind profile recommended by DNV-ST-N001 Marine Operations and Marine Warranty. \
    This wind profile can be used to convert the mean wind speed :math:`U_0` with averaging \
    period :math:`T_0` at height H above sea level to the mean wind speed U with averaging period \
    T at height z above sea level.
    
    The DNV profile includes a vertical wind speed profile, and a gust factor which allows for
    conversion of mean wind speeds between different averaging periods.
    
    The DNV profile provides an estimate to the mean wind speed :math:`U(T, z)` with averaging period of T seconds and at height z metres above mean sea level given by:

    .. math::
        U\\left(T, z\\right) = U_{0}\\left[1 + 0.137\\ln\\left(\\frac{z}{H}\\right)-0.047\\ln\\left(\\frac{T}{T_{0}}\\right)\\right]
    
    Where:
        * :math:`U_0` is the reference mean wind speed in metres per second
        * :math:`T` is the averaging period of :math:`U_0` in seconds
        * :math:`H` is the height of :math:`U_0` in metres above mean sea level

    Parameters
    ----------
    u : float
        Reference mean wind seed with averaging period T0 seconds, and at height z metres above mean sea level.
    T : int
        Target averaging period in seconds. Must be less than or equal to T0.
    T0 : int
        Averaging period of the reference mean wind speed in seconds.
    h_ref : float
        Height of the reference mean wind speed in metres above sea level.
    z : float
        Target height in metre above mean sea level.

    Returns
    -------
    float
        Mean wind speed with averaging period of T seconds at height z metres above sea level from the DNV wind profile.
    
    References
    ----------
    .. [1] DNVGL-ST-N001, September 2021 - Marine Operations and Marine Warranty.
        <https://www.dnv.com/oilgas/download/dnv-st-n001-marine-operations-and-marine-warranty.html>
    """
    correction_factor = 1 + 0.137 * np.log(z / h_ref) - 0.045 * np.log(T / T0)
    return u * correction_factor
