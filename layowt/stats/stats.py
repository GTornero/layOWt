""" 
Statistics Module.
"""

from scipy.special import gamma


def weibull_mean(weibull_a: float, weibull_k: float) -> float:
    """
    Calculates the mean value of a 2-parameter Weibull distribution.
    
    The mean value of a 2-parameter Weibull distribution is given by the following expression:
    
    .. math::
        E[x] = A \\Gamma\\left(1 + \\frac{1}{k}\\right)
    
    Where:
        * :math:`A = \\textrm{Weibull Scale parameter}`
        * :math:`k = \\textrm{Weibull Shape parameter}`
        * :math:`\\Gamma = \\textrm{Gamma function}`
    
    The Gamma function is implemented using :func:`scipy.special.gamma` from the Scipy library.

    Parameters
    ----------
    weibull_a : float
        The 2-parameter Weibull scale parameter.
    weibull_k : float
        The 2-parameter Weibull shape parameter.

    Returns
    -------
    float
        Mean value of the 2-parameter Weibull distribution.

    Raises
    ------
    ValueError
        If weibull_a or weibull_k are non positive.
    
    See Also
    --------
        scipy.special.gamma : Gamma function implemented by Scipy.
    """
    if (weibull_a <= 0) or (weibull_k <= 0):
        raise ValueError(f"a and k must be positive values. {weibull_a=}, {weibull_k=} passed.")
    
    return weibull_a * gamma(1 + 1 / weibull_k)