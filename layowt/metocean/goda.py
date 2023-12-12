import numpy as np
import scipy
from scipy.optimize import fsolve

from .wave_functions import deep_wavelength


class Goda:
    """ 
    Goda (1975) model for maximum wave height.
    
    The maximum individual wave height is obtained by
    
    .. math::
        H_{max} = H_{0.1\\%} = \\begin{cases}
            1.86K_{s}K_{r}H_{s0}, for \\frac{d}{L_{0}} \\geq 0.2 \\\\
            min\\left\\{\\left(\\beta_{0}^{*}K_{r}H_{s0}+\\beta_{1}^{*}d\\right), \\beta_{max}^{*}K_{r}H_{s0}, 1.86K_{s}K_{r}H_{s0}\\right\\}, for \\frac{d}{L_{0}} < 0.2
            \\end{cases}
    
    The coefficients :math:`\\beta_{0}^{*}`, :math:`\\beta_{1}^{*}`, and :math:`\\beta_{max}^{*}` are formulated as follows:
    
    .. math::
        \\beta_{0}^{*} = 0.054\\left(\\frac{K_{r}H_{s0}}{L_{0}}\\right)^{-0.38}\\exp\\left(20\\tan^{1.5}\\theta\\right)
    
    .. math::
        \\beta_{1}^{*} = 0.63\\exp\\left(3.8\\tan{\\theta}\\right)
    
    .. math::
        \\beta_{max}^{*} = max\\left\\{1.7, 0.55\\left(\\frac{K_{r}H_{s0}}{L_{0}}\\right)^{-0.29}\\exp\\left(2.4\\tan{\\theta}\\right)\\right\\}
    
    Where :math:`\\theta` is the angle of the sea floor with respect to the horizontal in the direction of wave propagation,
    and :math:`L_{0}` is the wave length in deep water.
    
    :math:`K_{s}` is a shallow water coefficient, i.e. the shoaling coefficient, determined by:
    
    .. math::
        K_{s} = \\left[\\tanh{kd} + kd\\left(1 - \\tanh^{2}{kd}\\right)\\right]^{-\\frac{1}{2}} + 0.0015\\left(\\frac{d}{L_{0}}\\right)^{-2.87}\\left(\\frac{K_{r}H_{s0}}{L_{0}}\\right)^{1.27}
    
    .. math::
        k = \\frac{2\\pi}{L}
    
    where :math:`L` is the wavelength at water depth :math:`d`, and is calculated by:
    
    .. math::
        L = \\frac{2\\pi d}{\\sqrt{P\\left(P + \\frac{1}{Q}\\right)}}
    
    .. math::
        P = \\frac{\\left(2\\pi\\right)^{2}d}{gT_{D}^{2}}
    
    .. math::
        Q = 1 + P\\left[0.6522 + P\\left(0.4622 + P^{2}\\left(0.0864+0.0675P\\right)\\right)\\right]
    
    :math:`K_{r}` is the refraction coefficient, determined by:
    
    .. math::
        K_{r} = \\sqrt{\\frac{\\cos\\alpha_{0}}{\\cos\\alpha_{r}}}
    
    where :math:`\\alpha_{0}` is the wave direction, and :math:`\\alpha_{r}` is the refracted wave direction at water depth :math:`d`.
    
    .. math::
        \\alpha_{r} = \\sin^{-1}\\left(\\sin{\\alpha_{0}}\\tanh{\\frac{2\\pi d}{L}}\\right)
    
    The equivalent deep water significant wave height, :math:`H_{s0}` can be calculated by iteration to satisfy the following equiations:
    
    .. math::
        H_{s} = \\begin{cases}
            K_{s}K_{r}H_{s0}, for \\frac{d}{L_{0}} \\geq 0.2 \\\\
            min\\left\\{\\left(\\beta_{0}K_{r}H_{s0}+\\beta_{1}d\\right), \\beta_{max}K_{r}H_{s0}, K_{s}K_{r}H_{s0}\\right\\}, for \\frac{d}{L_{0}} < 0.2
            \\end{cases}
    
    .. math::
        \\beta_{0} = 0.028\\left(\\frac{K_{r}H_{s0}}{L_{0}}\\right)^{-0.38}\\exp\\left(20\\tan^{1.5}\\theta\\right)
    
    .. math::
        \\beta_{1} = 0.52\\exp\\left(4.2\\tan\\theta\\right)
    
    .. math::
        \\beta_{max} = max\\left\\{0.92, 0.32\\left(\\frac{K_{r}H_{s0}}{L_{0}}\\right)^{-0.29}\\exp\\left(2.4\\tan{\\theta}\\right)\\right\\}
    
    Parameters
    ----------
    Hm0 : float
        The spectral wave height [m].
    depth : float
        The water depth [m].
    tp : float
        The spectral peak period [s].
    slope : float or tuple
        The slope of the foreshore [rad], or as tuple formatted as (V,H).
    alpha : float
        The wave direction [rad].
    
    Attributes
    ----------
    Hm0 : float
        Spectral wave height [m].
    depth : float
        The water depth [m].
    tp : float
        The spectral peak period [s].
    slope : float or tuple
        The slope of the foreshore [rad], or as tuple formatted as (V,H).
    alpha : float
        The wave direction [rad].
    L0 : float
        The deep water wave length [m].
    L : float
        The wave length at water depth :math:`d` [m].
    alpha_r : float
        The refracted wave direction at water depth :math:`d` [rad].
    Kr : float
        The refraction coefficient.
    Ks : float
        The shallow water coefficient.
    Hs0 : float
        Equivalent deep water significant wave height [m].
    Hmax : float
        Maximum wave height [m].
    
    See Also
    --------
    layowt.metocean.wave_functions.deep_wavelength : Calculate the deep water wavelength for lienar wave theory.
    
    References
    ----------
    .. [1] Goda Y., 1975, Deformation of irregular waves due to depth-controlled wave breaking, Rept. Port and Harbour Res. Inst. 14(3), pp.59-106.
    .. [2] IEC 61400-3-1: Wind Energy Generation Systems – Part 3-1: Design Requirements for Fixed Offshore Wind Turbines, edition 4.0, 2019-02.
        <https://webstore.iec.ch/publication/26423>
    
    """
    def __init__(self, Hm0: float, depth : float, tp: float, slope: float | tuple[float, float], alpha: float):
        
        self.Hm0 = Hm0
        self.depth = depth
        self.tp = tp
        self.alpha = alpha
        
        # check if slope foreshore is given as a tuple (V, H)
        if isinstance(slope, tuple):
            # if so, compute the angle
            slope = np.arctan(slope[0]/slope[1])
        
        self.slope = slope
        self.L0 = deep_wavelength(wave_period=self.tp)
        self.P = (((2 * np.pi)**2) * self.depth) / (scipy.constants.g * (self.tp**2))
        self.Q = 1 + self.P * (0.6522 + self.P * (0.4622 + (self.P**2) * (0.0864 + 0.0675 * self.P)))
        self.L = 2 * np.pi * 21 / (np.sqrt(self.P * (self.P + 1 / self.Q)))
        self.alpha_r = np.arcsin(np.sin(self.alpha * np.pi / 180)*np.tanh(2 * np.pi * self.depth/self.L))
        self.Kr = np.sqrt(np.cos(self.alpha * np.pi /180) / np.cos(self.alpha_r))
        self.k = 2 * np.pi / self.L
        
        self.Hs0 = float(fsolve(self._solver, [self.Hm0]))
        self.Ks = self._Ks(self.Hs0)
        self.Hmax = self._Hmax
    
    def _beta_0(self, x):
        return 0.028 * (self.Kr * x / self.L0)**(-0.38) * np.exp(20 * (np.tan(self.slope))**(1.5))
    
    @property
    def _beta_1(self):
        return 0.52 * np.exp(4.2 * np.tan(self.slope))
    
    def _beta_max(self, x):
        return np.max([0.92, 0.32 * (self.Kr * x / self.L0)**(-0.29) * np.exp(2.4 * np.tan(self.slope))])
    
    def _Ks(self, x):
        return (np.tanh(self.k*self.depth) + (self.k*self.depth) * (1 - (np.tanh(self.k*self.depth))**2))**(-0.5) + 0.0015 * (self.depth / self.L0)**(-2.87)*((self.Kr * x / self.L0)**1.27)
    
    def _solver(self, x):
        
        if self.depth / self.L0 >= 0.2:
            eq_b13 = [(self._Ks(x[0]) * self.Kr * x[0]) - self.Hm0]
        else:
            a = (self._beta_0(x[0]) * self.Kr * x[0]) + (self._beta_1 * self.depth)
            b = self._beta_max(x[0]) * self.Kr * x[0]
            c = self._Ks(x[0]) * self.Kr * x[0]
            
            eq_b13 = [np.min([a, b, c]) - self.Hm0]
        
        return eq_b13
    
    @property
    def _beta_0_star(self):
        return 0.054 * ((self.Kr * self.Hs0 / self.L0)**(-0.38)) * np.exp(20 * (np.tan(self.alpha)**1.5))
    
    @property
    def _beta_1_star(self):
        return 0.63 * np.exp(3.8 * np.tan(self.alpha))
    
    @property
    def _beta_max_star(self):
        return np.max([1.7, 0.55 * ((self.Kr * self.Hs0 / self.L0)**(-0.29)) * np.exp(2.4 * np.tan(self.alpha))])
    
    @property
    def _Hmax(self):
        
        if self.depth / self.L0 >= 0.2:
            return 1.86 * self.Ks * self.Kr * self.Hs0
        else:
            a = (self._beta_0_star * self.Kr * self.Hs0) + (self._beta_1_star * self.depth)
            b = self._beta_max_star * self.Kr * self.Hs0
            c = 1.86 * self.Ks * self. Kr * self.Hs0
            
            return np.min([a, b, c])