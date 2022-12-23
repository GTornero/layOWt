""" 
Circular Variable Statistics Module.
"""

from typing import Iterable

import numpy as np


def circ_mean(x: Iterable[float], units: str = 'degrees') -> float:
    """
    Calculates the mean value of circular data.
    
    Each observation is treated as a unit vector, or point in the unit circle. The resultant vector
    of the observations is found, and the direction of the resultant vector is returned in radians.
    
    The mean value of circular data is found using the following expression:
    
    .. math::
        \\bar{\\alpha} = arctan2\\left(\\frac{\\sum_{i=1}^n sin(\\alpha_{i})}{\\sum_{i=1}^n cos(\\alpha_{i})}\\right)
    
    Where:
        * :math:`\\alpha_{1}, \\alpha_{2}, ... \\alpha_{n} = \\textrm{Samples of circular data.}`
        * :math:`n = \\textrm{Number of samples.}`

    Parameters
    ----------
    x : Iterable[float]
        Circular data.
    units : str, optional
        Units of the circular data. Can be one of 'degrees' or 'radians', case insensitive.
        By default 'degrees'.

    Returns
    -------
    float
        Mean direction of circular data in radians.
    
    Raises
    ------
    ValueError
        If units is not one of 'degrees' or 'radians'.
    
    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
        Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
        Circular Statistics (2001)'". 2015.
        <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    
    if units.lower() not in ("degrees", "radians"):
        raise ValueError(f"Invalid value for units passed. Must be one of 'degrees' or 'radians'. Value of {units=} passed.")
    
    alpha = np.array(x)
    if units.lower() == "degrees":
        alpha = np.deg2rad(alpha)
    
    sinr = np.sum(np.sin(alpha))
    cosr = np.sum(np.cos(alpha))
    
    return np.arctan2(sinr, cosr)

def circ_corr(alpha: Iterable[float], beta: Iterable[float], units: str = 'degrees') -> float:
    """
    Computes the circular correlation coefficient between two arrays of direction data.
    
    The circular correlation coefficient is given by the following expression:
    
    .. math::
        \\rho_c = \\frac{\\sum_{i=1}^n \\sin(\\alpha_i - \\bar{\\alpha})\\sin(\\beta_i - \\bar{\\beta})}{\\sqrt{\\sum_{i=1}^n \\sin^{2}(\\alpha_i - \\bar{\\alpha})\\sum_{i=1}^n\\sin^{2}(\\beta_i - \\bar{\\beta})}}
    
    Where:
        * :math:`\\alpha_{1}, \\alpha_{2}, ... \\alpha_{n} = \\textrm{Samples of first circular variable.}`
        * :math:`\\beta_{1}, \\beta_{2}, ... \\beta_{n} = \\textrm{Samples of second circular variable.}`
        * :math:`\\bar{\\alpha} = \\textrm{Mean value of the samples of the first circular variable.}`
        * :math:`\\bar{\\beta} = \\textrm{Mean value of the samples of the second circular variable.}`
        * :math:`n = \\textrm{Number of samples.}`
    
    The circular correlation coefficient satisfies the following properties:
        * :math:`\\rho(\\alpha, \\beta)` does not depend on the zero direction used for either variable;
        * :math:`\\rho(\\alpha, \\beta) = \\rho(\\beta, \\alpha)`;
        * :math:`\\left|\\rho(\\alpha, \\beta)\\right| \\le 1`;
        * :math:`\\rho(\\alpha, \\beta) = 0` if the variables are independent;
        * :math:`\\rho(\\alpha, \\beta) = 1` if the variables are perfectly correlated.

    Parameters
    ----------
    alpha : Iterable[float]
        Array of circular (directional) data.
    beta : Iterable[float]
        Array of circular (directional) data.
    units : str, optional
        Units of the circular data. Can be one of 'degrees' or 'radians', case insensitive.
        By default 'degrees'.

    Returns
    -------
    float
        Circular Correlation Coefficient.

    Raises
    ------
    ValueError
        If units is not one of 'degrees' or 'radians'.
    
    References
    ----------
    .. [1] S. R. Jammalamadaka, A. SenGupta. "Topics in Circular Statistics".
        Series on Multivariate Analysis, Vol. 5, 2001.
    .. [2] Jammalamadaka, S. and Sarma, Y. (1988). A correlation coefficient for angular variables.
        Statistical Theory and Data Analysis 2. North Holland: New York.
    .. [3] C. Agostinelli, U. Lund. "Circular Statistics from 'Topics in
        Circular Statistics (2001)'". 2015.
        <https://cran.r-project.org/web/packages/CircStats/CircStats.pdf>
    """
    
    if units.lower() not in ("degrees", "radians"):
        raise ValueError(f"Invalid value for units passed. Must be one of 'degrees' or 'radians'. Value of {units=} passed.")
    
    a = np.array(alpha)
    b = np.array(beta)
    
    if units.lower() == "degrees":
        a = np.deg2rad(a)
        b = np.deg2rad(b)
    
    mu_a = circ_mean(a, units='radians')
    mu_b = circ_mean(b, units='radians')
    
    sin_1 = np.sin(a - mu_a)
    sin_2 = np.sin(b - mu_b)
    
    num = np.sum(sin_1 * sin_2)
    den = np.sqrt(np.sum(sin_1 ** 2) * np.sum(sin_2 ** 2))
    
    rho = num / den
    
    return rho