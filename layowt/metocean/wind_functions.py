import numpy as np
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

# TODO: Froya Equation 