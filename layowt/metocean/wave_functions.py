import matplotlib.pyplot as plt
import numpy as np
import scipy
from coastlib.waves.support import wave_theories
from matplotlib.axes import Axes
from matplotlib.figure import Figure


def deep_wavelength(wave_period: float) -> float:
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

def intermediate_wavelength_approx(depth: float, wave_period:  float) -> float:
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
        \\frac{1}{20} \\leq \\frac{d}{\\lambda} \\leq \\frac{1}{2}
    
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

def intermediate_wavelength(depth: float, wave_period: float, error: float=0.1, max_iter: int=100) -> float:
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
        \\frac{1}{20} \\leq \\frac{d}{\\lambda} \\leq \\frac{1}{2}
    
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

def shallow_wavelength(depth: float, wave_period: float) -> float:
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
        \\frac{d}{\\lambda} \\leq \\frac{1}{20}
    
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

def breaking_height(depth: float, wavelength: float) -> float:
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

def iribarren(slope: float, breaking_height: float, wave_period: float) -> float:
    """
    Calculates the Iribarren number (also known as surf similarity parameter).
    
    The Iribarren number is a dimensionless parameter used to model several effects of breaking
    waves on beaches an coastal structures.
    
    The Iribarren number is given by:
    
    .. math::
        \\xi_{b} = \\frac{\\tan\\alpha}{\\sqrt{\\frac{H_{b}}{\\lambda_{deep}}}}
    
    Where:
        * :math:`\\alpha =  \\textrm{Seabed slope.}`
        * :math:`H_{b} =  \\textrm{Wave height at breaking point.}`
        * :math:`\\lambda_{deep} =  \\textrm{Deep water wavelength.}`

    The type of breakign wave - spilling, plunging, collapsing, or surging - depends on the
    Irabarren number. The approximate breaker type dependence  on the Irabarren number is given by:
    
    =====================  ============================
    Breaker Type                 :math:`\\xi_{b}`
    =====================  ============================
    Spilling               :math:`\\xi_{b} < 0.4`
    Plunging               :math:`0.4 < \\xi_{b} < 2.0`
    Surging or Collapsing  :math:`\\xi_{b} > 2.0`
    =====================  ============================
    
    .. image:: ../../../../_static/breaking_wave_types.png
        :width: 400
        :align: center

    Parameters
    ----------
    slope : float
        Angle of the seabed slope in degrees.
    breaking_height : float
        Wave height at breaking point in metres.
    wave_period : float
        Wave period in seconds.

    Returns
    -------
    float
        Iribarren number.
    
    See Also
    --------
    layowt.metocean.wave_functions.deep_wavelength : Calculate the deep water wavelength for linear wave theory.
    layowt.metocean.wave_functions.breaking_height : Calculates the breaking wave height for a regular wave in a plane sea bed.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    l_0 = deep_wavelength(wave_period=wave_period)
    
    return np.tan(np.deg2rad(slope)) / np.sqrt(breaking_height / l_0)

def breaker_type(iribarren: float) -> str:
    """
    Determine the type of breaking wave based on the Iribarren number.
    
    The type of breakign wave - spilling, plunging, collapsing, or surging - depends on the
    Irabarren number. The approximate breaker type dependence  on the Irabarren number is given by:
    
    =====================  ============================
    Breaker Type                 :math:`\\xi_{b}`
    =====================  ============================
    Spilling               :math:`\\xi_{b} < 0.4`
    Plunging               :math:`0.4 < \\xi_{b} < 2.0`
    Surging or Collapsing  :math:`\\xi_{b} > 2.0`
    =====================  ============================
    
    .. image:: ../../../../_static/breaking_wave_types.png
        :width: 400
        :align: center

    Parameters
    ----------
    iribarren : float
        Iribarren number.

    Returns
    -------
    str
        Type of breaking wave.
    
    See Also
    --------
    layowt.metocean.wave_functions.iribarren : Calculates the Iribarren number (also known as surf similarity parameter).
    """
    if iribarren < 0.4:
        return "Spilling"
    elif iribarren < 2:
        return "Pluging"
    else:
        return "Surging or Collapsing"

def second_order_crest(wave_height: float, wave_length: float) -> float:
    """
    Calculates the wave crest height of a Second-order Stokes wave in deep water.
    
    The wave crest height for a Second-order Stokes wave in deep water is given by:
    
    .. math::
        A_{C} = \\frac{H}{2}\\left(1 + \\frac{\\pi H}{2\\lambda}\\right)
    
    Where:
        * :math:`H =  \\textrm{Wave height.}`
        * :math:`\\lambda =  \\textrm{Wavelength.}`

    Parameters
    ----------
    wave_height : float
        Wave height in metres.
    wave_length : float
        Wavelength in metres.

    Returns
    -------
    float
        Second-order Stokes wave crest height.
    
    See Also
    --------
    layowt.metocean.wave_functions.second_order_trough : Calculates the wave trough height of a Second-order Stokes wave in deep water.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    return (wave_height / 2) * (1 + (np.pi * wave_height / (2 * wave_length)))

def second_order_trough(wave_height: float, wave_length: float) -> float:
    """
    Calculates the wave trough height of a Second-order Stokes wave in deep water.
    
    The wave trough height for a Second-order Stokes wave in deep water is given by:
    
    .. math::
        A_{T} = \\frac{H}{2}\\left(1 - \\frac{\\pi H}{2\\lambda}\\right)
    
    Where:
        * :math:`H =  \\textrm{Wave height.}`
        * :math:`\\lambda =  \\textrm{Wavelength.}`

    Parameters
    ----------
    wave_height : float
        Wave height in metres.
    wave_length : float
        Wavelength in metres.

    Returns
    -------
    float
        Second-order Stokes wave trough height.
    
    See Also
    --------
    layowt.metocean.wave_functions.second_order_crest : Calculates the wave crest height of a Second-order Stokes wave in deep water.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    return (wave_height / 2) * (1 - (np.pi * wave_height / (2 * wave_length)))

def pierson_moskowitz(f: float, hs: float, tp: float) -> float:
    """
    Calculates the spectral density of the Pierson-Moskowitz wave spectrum for a sea state at a given frequency.
    
    The spectral density of the Pierson-Moskowitz spectrum is given by:
    
    .. math::
        S_{PM}(\\omega) = \\frac{5}{16}H_{s}^{2}\\omega_{p}^{4}\\omega^{-5}e^{-\\frac{5}{4}\\left(\\frac{\\omega}{\\omega_{p}}\\right)^{-4}}
    
    Where:
        * :math:`H_s =  \\textrm{Significant wave height.}`
        * :math:`\\omega_p = \\frac{2\\pi}{T_p} \\textrm{, the angular spectral peak frequency.}`
        * :math:`\\omega = 2\\pi f \\textrm{, the angular frequency.}`
        

    Parameters
    ----------
    f : float
        Frequency in Hertz.
    hs : float
        Significant wave height of the sea state in metres.
    tp : float
        Spectral peak wave period of the sea state in seconds.

    Returns
    -------
    float
        The spectral density of the Pierson-Moskowitz spectrum of the sea state in :math:`m^2/Hz`.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    .. [2] Pierson, Willard J., Jr. and Moskowitz, Lionel A. Proposed Spectral Form for Fully Developed Wind Seas Based on the Similarity Theory of S. A. Kitaigorodskii, Journal of Geophysical Research, Vol. 69, p.5181-5190, 1964.
    """
    omega  = 2*np.pi*f
    omega_p = 2*np.pi/tp
    return (5 / 16) * (hs**2) * (omega_p**4) * (omega**(-5)) * np.exp(-(5 / 4) * ((omega / omega_p)**-4))

def plot_pierson_moskowitz(hs: float, tp: float, x_units: str = "hz", ax: Axes | None = None, figsize: tuple[float, float] = (10,10)) -> tuple[Figure, Axes]:
    """
    Plots the Pierson-Moskowitz wave spectrum for a given sea state.

    Parameters
    ----------
    hs : float
        Significant wave height in metres.
    tp : float
        Spectral peak wave period in seconds.
    x_units : str
        X-axis units, must be one of "Hs" or "rad". By default "Hz".
    ax : matplotlib.axes.Axes | None, optional
        matplotlib.axes.Axes, by default None
    figsize : tuple[float, float], optional
        Figure size in inches in (width, height) format, by default (10,10)

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        Tuple of Figure and Axes objects.

    Raises
    ------
    TypeError
        The passed ax argument was not of type Axes.
    ValueError
        If x_units is not "Hz" or "rad".
    
    See Also
    --------
    layowt.metocean.wave_functions.pierson_moskowitz : Calculates the spectral density of the Pierson-Moskowitz wave spectrum for a sea state at a given frequency.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        if x_units.lower() == "hz":
            ax.set_xlabel(r"$f [Hz]$")
        elif x_units.lower() == "rad":
            ax.set_xlabel(r"$\omega [rad/s]$")
        ax.set_ylabel(r"$S(\omega) [m^2/Hz]")
        ax.grid()
    else:
        try:
            fig = ax.get_figure()
        except AttributeError:
            raise TypeError from AttributeError(
                f"Invalid type in {type(ax)} for the 'ax' argument, "
                f"must be Matplotlib Axes object."
            )
    
    f =  np.arange(0.01, 1, 0.005)
    s = pierson_moskowitz(f=f, hs=hs, tp=tp)
    if x_units.lower() == "hz":
        ax.plot(f, s, label=f"Hs: {hs:.2f}m, Tp: {tp:.2f}s")
    elif x_units.lower() == "rad":
        ax.plot(2*np.pi*f, s, label=f"Hs: {hs:.2f}m, Tp: {tp:.2f}s")
    else:
        raise ValueError(f'Invalid value of X_units passed. Must be one of "Hz" or "rad". {x_units=} passed.')

    if ax.get_legend() is not None:
        ax.get_legend().remove()
        ax.legend()
    else:
        ax.legend()
    
    return fig, ax

def jonswap(f: float, hs: float, tp: float, gamma: float) -> float:
    """
    Calculates the spectral density of the JONSWAP wave spectrum for a sea state at a given frequency.
    
    The spectral density of the JONSWAP spectrum is given by:
    
    .. math::
        S_{J}(\\omega) = A_{\\gamma}S_{PM}(\\omega)\\gamma^{\\exp\\left(-0.5\\left(\\frac{\\omega - \\omega_p}{\\sigma \\omega_p}\\right)^2\\right)}
    
    Where:
        * :math:`S_{PM}(\\omega) =  \\textrm{Pierson-Moskowitz spectrum.}`
        * :math:`\\omega_p = \\frac{2\\pi}{T_p} \\textrm{, the angular spectral peak frequency.}`
        * :math:`\\omega = 2\\pi f \\textrm{, the angular frequency.}`
        * :math:`\\gamma =  \\textrm{JONSWAP peak enhancement factor.}`
        * :math:`\\sigma = 0.07 \\textrm{ for } \\omega \\leq \\omega_p, 0.09 \\textrm{ for } \\omega > \\omega_p`
        * :math:`A_{\\gamma} =  \\frac{0.2}{0.065\\gamma^{0.803} + 0.135}`
        
    Parameters
    ----------
    f : float
        Frequency in Hertz.
    hs : float
        Significant wave height of the sea state in metres.
    tp : float
        Spectral peak wave period of the sea state in seconds.
    gamma : float
        JONSWAP peak enhancement factor of the sea state.

    Returns
    -------
    float
        The spectral density of the JONSWAP spectrum of the sea state in :math:`m^2/Hz`.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    .. [2] Hasselmann K., T.P. Barnett, E. Bouws, H. Carlson, D.E. Cartwright, K. Enke, J.A. Ewing, H. Gienapp, D.E. Hasselmann, P. Kruseman, A. Meerburg, P. Mller, D.J. Olbers, K. Richter, W. Sell, and H. Walden. Measurements of wind-wave growth and swell decay during the Joint North Sea Wave Project (JONSWAP)' Ergnzungsheft zur Deutschen Hydrographischen Zeitschrift Reihe, A(8) (Nr. 12), p.95, 1973.
        <https://repository.tudelft.nl/islandora/object/uuid:f204e188-13b9-49d8-a6dc-4fb7c20562fc/datastream/OBJ/>
    
    See Also
    --------
    layowt.metocean.wave_functions.pierson_moskowitz : Calculates the spectral density of the Pierson Moskowitz wave spectrum for a sea state at a given frequency.
    """
    omega  = 2*np.pi*f
    omega_p = 2*np.pi/tp
    sigma = np.where(omega<=omega_p, 0.07, 0.09)
    
    a = 0.2 / ((0.065 * (gamma**0.803)) + 0.135)
    s_pm = pierson_moskowitz(f=f, hs=hs, tp=tp)
    
    s_j = a * s_pm * (gamma**np.exp(-0.5 * ((omega - omega_p)/sigma)**2))
    
    return s_j

def plot_jonswap(hs: float, tp: float, gamma: float, x_units: str = "hz", ax: Axes | None = None, figsize: tuple[float, float] = (10,10)) -> tuple[Figure, Axes]:
    """
    Plots the JONSWAP wave spectrum for a given sea state.

    Parameters
    ----------
    hs : float
        Significant wave height of the sea state in metres.
    tp : float
        Spectral peak wave period of the sea state in seconds.
    gamma : float
        JONSWAP peak enhancement factor of the sea state.
    x_units : str
        X-axis units, must be one of "Hs" or "rad". By default "Hz".
    ax : matplotlib.axes.Axes | None, optional
        matplotlib.axes.Axes, by default None
    figsize : tuple[float, float], optional
        Figure size in inches in (width, height) format, by default (10,10)

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        Tuple of Figure and Axes objects.

    Raises
    ------
    TypeError
        The passed ax argument was not of type Axes.
    ValueError
        If x_units is not "Hz" or "rad".
    
    See Also
    --------
    layowt.metocean.wave_functions.jonswap : Calculates the spectral density of the JONSWAP wave spectrum for a sea state at a given frequency.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
        if x_units.lower() == "hz":
            ax.set_xlabel(r"$f [Hz]$")
        elif x_units.lower() == "rad":
            ax.set_xlabel(r"$\omega [rad/s]$")
        ax.set_ylabel(r"$S(\omega) [m^2/Hz]$")
        ax.grid()
    else:
        try:
            fig = ax.get_figure()
        except AttributeError:
            raise TypeError from AttributeError(
                f"Invalid type in {type(ax)} for the 'ax' argument, "
                f"must be Matplotlib Axes object."
            )
    
    f =  np.arange(0.01, 1, 0.005)
    s = jonswap(f=f, hs=hs, tp=tp, gamma=gamma)
    if x_units.lower() == "hz":
        ax.plot(f, s, label=f"Hs: {hs:.2f}m, Tp: {tp:.2f}s, $\\gamma$: {gamma:.2f}")
    elif x_units.lower() == "rad":
        ax.plot(2*np.pi*f, s, label=f"Hs: {hs:.2f}m, Tp: {tp:.2f}s, $\\gamma$: {gamma:.2f}")
    else:
        raise ValueError(f'Invalid value of X_units passed. Must be one of "Hz" or "rad". {x_units=} passed.')
    
    if ax.get_legend() is not None:
        ax.get_legend().remove()
        ax.legend()
    else:
        ax.legend()
    
    return fig, ax

def ursell(wave_height: float, wave_length: float, depth: float) -> float:
    """
    Calculate the Ursell number.
    
    The Ursell number indicates the nonlinearity of waves on a fluid layer.
    
    The Ursell number is defined as:
    
    .. math::
        U_R = \\frac{H\\lambda^2}{h^3}
    
    Where:
        * :math:`H = \\textrm{Wave height in metres.}`
        * :math:`\\lambda = \\textrm{Wavelength in metres.}`
        * :math:`d = \\textrm{water depth in metres.}`

    Parameters
    ----------
    wave_height : float
        Wave height in metres.
    wave_length : float
        Wavelength in metres.
    depth : float
        Water depth in metres.

    Returns
    -------
    float
        Ursell number.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    .. [2] Ursell, F (1953). "The long-wave paradox in the theory of gravity waves". Proceedings of the Cambridge Philosophical Society. 49 (4): 685–694.
        <https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/longwave-paradox-in-the-theory-of-gravity-waves/5A178FB13BD7B9C3314A49A495A860DB>
    """
    return wave_height * (wave_length**2) / (depth**3)

def validate_theory(wave_height: float, wave_period: float, depth:float) -> tuple[Figure, Axes]:
    """
    Shows what wave theory is applicable for given wave parameters using the
    wave theories' figure by Le Mehaute per USACE CEM Part II Chap. 1 p.II-1-58
    
    .. image:: ../../../../_static/LeMehaute.png
        :align: center
    
    Parameters
    ----------
    wave_height : float
        Wave height in metres.
    wave_period : float
        Wave period in seconds.
    depth : float
        Water depth in metres

    Returns
    -------
    tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
        Tuple of Figure and Axes objects.
        
    References
    ----------
    .. [1] Part 2 of Coastal engineering manual, USA Corps of Engineers Volume 1110, Issues 2-1100 of Engineer manual, USA Corps of Engineers
        <https://www.plainwater.com/pubs/EM-1110-2-1100P2.pdf>
    """
    return wave_theories(wave_height=wave_height, wave_period=wave_period, depth=depth)

def _t02_over_t01(gamma: float , t02: float , t01: float) -> float:
    """
    Function of T02 / T01 based on JONSWAP period relationships given in DNVGL-RP-C205, 19th Edition, December 2020, section 3.5.5.4.
    Used to solve for a JONSWAP gamma given values ot T02 and T01.

    Parameters
    ----------
    gamma : float
        JONSWAP peak enhancement factor.
    t02 : float
        Mean wave period based on second sprectral moment.
    t01 : float
        Mean wave period based on first sprectral moment.

    Returns
    -------
    float
        Function output.
    """
    numerator = 0.6673 + (0.05037 * gamma) - (0.006230 * gamma**2) + (0.0003341 * gamma**3)
    denominator = 0.7303 + (0.04936 * gamma) - (0.006556 * gamma**2) + (0.0003610 * gamma**3)
    return (numerator / denominator) - (t02 / t01)

def gamma_from_mean_periods(t02: float, t01: float) -> float:
    """
    Estimate the JONSWAP peak enhancement factor :math:`\\gamma` from the spectral mean periods :math:`T_{m02}` and :math:`T_{m01}.`
    
    The spectal mean wave periods :math:`T_{m02}` and :math:`T_{m01}` can be related to the spectral peak wave period
    :math:`T_p` by the following approximate empirical relations (generally valid for :math:`1 \\leq \\gamma < 7`):
    
    .. math::
        \\frac{T_{m02}}{T_p} = 0.6673 + 0.05037\\gamma - 0.006230\\gamma^2 + 0.0003341\\gamma^3
        
        \\frac{T_{m01}}{T_p} = 0.7303 + 0.04936\\gamma - 0.006556\\gamma^2 + 0.0003610\\gamma^3
    
    These equations can be rearranged to eliminate :math:`T_p` and the root of the resulting equation gives us the value of gamma:
    
    .. math::
        \\frac{T_{m02}}{T_{m01}} - \\frac{0.6673 + 0.05037\\gamma - 0.006230\\gamma^2 + 0.0003341\\gamma^3}{0.7303 + 0.04936\\gamma - 0.006556\\gamma^2 + 0.0003610\\gamma^3} = 0

    Parameters
    ----------
    t02 : float
        Mean spectral wave period based on the second moment in seconds.
    t01 : float
        Mean spectral wave period based on the first moment in seconds.

    Returns
    -------
    float
        JONSWAP Peak Enhancement Factor, :math:`\\gamma`
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    return scipy.optimize.root(fun=_t02_over_t01, x0=[3.3], args=(t02, t01)).x[0]

def tp_from_mean_period(t: float, gamma: float, period: str = "t02") -> float:
    """
    Estimate the spectral peak wave period :math:`T_p` from the JONSWAP peak enhancement factor and a mean spectral period.
    
    The spectral peak wave period can be related to the JONSWAP peak enhancement factor and a mean spectral period
    (:math:`T_{m02}` or :math:`T_{m01}.`) using the following empirical relationships (generally valid for :math:`1 \\leq \\gamma < 7`):
    
    .. math::
        T_p = \\frac{T_{m02}}{0.6673 + 0.05037\\gamma - 0.006230\\gamma^2 + 0.0003341\\gamma^3}
        
        T_p = \\frac{T_{m01}}{0.7303 + 0.04936\\gamma - 0.006556\\gamma^2 + 0.0003610\\gamma^3} 

    Parameters
    ----------
    t : float
        Mean spectral wave period in seconds.
    gamma : float
        JONSWAP peak enhancement factor.
    period : str, optional
        Selection of which mean spectral period to use. Must be one of "t02" or "t01", by default "t02".

    Returns
    -------
    float
        Approximation of the spectral peak wave period.

    Raises
    ------
    ValueError
        If "t02" or "t01" was not passed to the "period" argument.
    
    References
    ----------
    .. [1] DNVGL-RP-C205, 19th Edition, December 2020 - Environmental conditions and environmental loads.
        <https://www.dnv.com/oilgas/download/dnv-rp-c205-environmental-conditions-and-environmental-loads.html>
    """
    if period.lower() == "t02":
        return  t / (0.6673 + (0.05037 * gamma) - (0.006230 * gamma**2) + (0.0003341 * gamma**3))
    elif period.lower() == "t01":
        return t / (0.7303 + (0.04936 * gamma) - (0.006556 * gamma**2) + (0.0003610 * gamma**3))
    else:
        raise ValueError(f'Invalid value passed for "period". Must be one of "t02" or "t01". {period=} was passed.')