import copy

import matplotlib.pyplot as plt
import numpy as np
import scipy


def deep_wavelength(wave_period: float):
    """
    Calculate the deep water wavelength from the wave period. Applicable for linear wave theory.
    
    The deep water wavelength is given by the following expression:
    
    .. math::
        \\lambda_{deep} = \\frac{gT}{2\\pi}
    
    Where:
        * :math:`\\lambda_{deep} = \\textrm{Deep water wavelength in metres.}`
        * :math:`g = \\textrm{Acceleration due to gravity in }ms^{-2}.`
        * :math:`T = \\textrm{Wave period in seconds.}`
    
    This formulation is recommended when:
    
    .. math::
        \\frac{d}{\\lambda} > \\frac{1}{2}
    
    Where:
        * :math:`d = \\textrm{Water depth in metres.}`

    Parameters
    ----------
    wave_period : float
        Wave period in seconds.

    Returns
    -------
    float
        Deep water wavelength in metres.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    return ((scipy.constants.g * wave_period**2) / (2 * np.pi))

def intermediate_wavelength_approx(depth: float, wave_period:  float):
    """
    Calculate the approximate intermediate water depth wavelength from the wave period and water depth. Applicable for linear wave theory.
    
    The deep water wavelength is given by the following expression:
    
    .. math::
        \\lambda_{int} = \\lambda_{deep}\\tanh\\left(\\varepsilon_0\\right)
    
    Where:
        * :math:`\\lambda_{deep} = \\frac{gT}{2\\pi} \\textrm{ (Deep water wavelength in metres).}`
        * :math:`g = \\textrm{Acceleration due to gravity in }ms^{-2}.`
        * :math:`T = \\textrm{Wave period in seconds.}`
        * :math:`\\varepsilon_0 = \\sqrt{k_{0}d}\\left(1 + \\frac{k_{0}d}{6} + \\frac{\\left(k_{0}d\\right)^2}{30}\\right)`
        * :math:`k_{0} = \\frac{2\\pi}{\\lambda_0}`
        * :math:`d = \\textrm{Water depth in metres.}`
    
    This approximate formulation is recommended when:
    
    .. math::
        \\frac{1}{20} <= \\frac{d}{\\lambda} <= \\frac{1}{2}
    
    Parameters
    ----------
    depth : float
        Water depth in metres.
    wave_period : float
        Wave period in seconds.

    Returns
    -------
    float
        Approximate intermediate water depth wavelength in metres.
    
    See Also
    --------
    layowt.metocean.wave_functions.deep_wavelength : Calculate the deep water wavelength for lienar wave theory.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    g = scipy.constants.g
    l_0 = deep_wavelength(wave_period=wave_period)
    
    k_0 = (2 * np.pi) / l_0
    eps = np.sqrt(k_0 * depth) * (1 + (k_0 * depth / 6) + (((k_0 * depth)**2)/30))
    
    return l_0 * np.tanh(eps)

def intermediate_wavelength(depth: float, wave_period: float, error: float=0.1, max_iter: int=100):
    """
    Calculate the intermediate water depth wavelength from the wave period and water depth iteratively. Applicable for linear wave theory.
    
    The deep water wavelength is given by the following expression:
    
    .. math::
        \\lambda_{int} = \\lambda_{deep}\\tanh\\left(\\frac{2\\pi d}{\\lambda_{int}}\\right)
    
    Where:
        * :math:`\\lambda_{deep} = \\frac{gT}{2\\pi} \\textrm{ (Deep water wavelength in metres).}`
        * :math:`g = \\textrm{Acceleration due to gravity in }ms^{-2}.`
        * :math:`T = \\textrm{Wave period in seconds.}`
        * :math:`d = \\textrm{Water depth in metres.}`
    
    This formulation is recommended when:
    
    .. math::
        \\frac{1}{20} <= \\frac{d}{\\lambda} <= \\frac{1}{2}
    
    Parameters
    ----------
    depth : float
        Water depth in metres.
    wave_period : float
        Wave period in seconds.
    error : float
        Error tolerance for result convergence.
    max_iter : int
        Maximum number of iterations allowable iterations.

    Returns
    -------
    float
        Intermediate water depth wavelength in metres.
    
    See Also
    --------
    layowt.metocean.wave_functions.deep_wavelength : Calculate the deep water wavelength for linear wave theory.
    layowt.metocean.wave_functions.intermediate_wavelength_approx : An apprximation to the linear wave theory intermediate depth wavelength.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    l_0 = deep_wavelength(wave_period=wave_period)
    l = l_0
    l_new = l_0 * np.tanh(2 * np.pi * depth / l)

    n = 0
    
    while (np.abs(l_new - l) > error):
        if n>max_iter:
            raise StopIteration(f"Not converged after {n-1} iterations.")
        
        l = l_new
        l_new = l_0 * np.tanh(2 * np.pi * depth / l)
        n+=1
    
    return l

def shallow_wave_length(depth: float, wave_period: float):
    """
    Calculate the shallow water depth wavelength from the wave period and water depth. Applicable for linear wave theory.
    
    The shallow water wavelength is given by the following expression:
    
    .. math::
        \\lambda_{shallow} = T\\sqrt{gd}
    
    Where:
        * :math:`g = \\textrm{Acceleration due to gravity in }ms^{-2}.`
        * :math:`T = \\textrm{Wave period in seconds.}`
        * :math:`d = \\textrm{Water depth in metres.}`
    
    This formulation is recommended when:
    
    .. math::
        \\frac{d}{\\lambda} <= \\frac{1}{20}
    
    Parameters
    ----------
    depth : float
        Water depth in metres.
    wave_period : float
        Wave period in seconds.

    Returns
    -------
    float
        Shallow water depth wavelength in metres.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    return wave_period * np.sqrt(scipy.constants.g * depth)

def breaking_height(depth: float, wavelength: float):
    """
    Calculates the breaking wave height for a regular wave in a plane sea bed.
    
    Breaking wave height is given by:
    
    .. math::
        H_b = 0.142\\lambda\\tanh\\frac{2\\pi d}{\\lambda}
    
    Where:
        * :math:`d = \\textrm{Water depth in metres.}`
        * :math:`\\lambda = \\textrm{Wavelength corresponding to water depth d in metres.}`
    
    In deep water the breaking wave limit corresponds to a maximum steepnes :math:`S_{max} = \\frac{H_b}{\\lambda} = \\frac{1}{7}`.

    Parameters
    ----------
    depth : float
        Water depth in metres.
    wavelength : float
        Wavelength in metres.

    Returns
    -------
    float
        Breaking wave height in metres.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    hb = 0.142 * wavelength * np.tanh(2 * np.pi * depth / wavelength)
    return hb