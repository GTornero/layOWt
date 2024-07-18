import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sc
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from scipy.optimize import fsolve

from layowt.metocean.wave_functions import rayleigh_hp


class BattjesGroenendijk:
    """
    Formulas of Battjes and Groenendijk (2000) Wave Height Distributions on Shallow Foreshores. 

    The wave height distributions on shallow foreshores deviates from
    those in deep water due to the effects of the restricted
    depth-to-height ratio and of wave breaking. Battjes and Groenendijk
    (2000), therefore derived a generalised empirical parameterisations
    based on laboratory data to determine these effects. This allows for
    the computation of all statistical wave heights based on the spectral
    wave height, water depth and slope of the foreshore.
    
    The battjes and Groenendijk distribution is a Composite Weibull Distribution of wave heights
    matched at the transition wave height :math:`H_{tr}` given by:
    
    .. math::
        F(H) = Pr\\{x \\le H\\} = \\begin{cases}
            F_{1}(H) = 1 - \\exp\\left[-\\frac{H}{H_1}^{K_1}\\right], H \\le H_{tr} \\\\
            F_{2}(H) = 1 - \\exp\\left[-\\frac{H}{H_2}^{K_2}\\right], H \\geq H_{tr}
            \\end{cases}
    
    In order to obtain the continuity of the distribution function, the constraint
    :math:`F_{1}(H_{tr})=F_{2}(H_{tr})` is imposed. The shape parameters :math:`k_1` and :math:`k_1`
    are shape parameters of the distributions and have been set to 2 and 3.6 respectively
    from observations (these shape parameters can be modified).
    
    The transition wave height is determined from the depth and the seabed slope, given (from observations) by:
    
    .. math::
        H_{tr} = 0.35 + 5.8tan(\\alpha)d
    
    As per the IEC 61400-3-1, Appendix B.3.3 recommendations, the following conditions should be
    fulfilled if the Battjes and Groenendijk distribution is to be used:
    
    * The direction of wave propagation is within 30° of the direction of the steepest seabed slope.
    
    * Considering values of seabed slope averaged over distances of 1 and 3 deepwater wavelengths in \\
      the upwave direction from the turbine position, the greater of the two values should be \\
      positive and not greater than 0.05 rad.
    
    * Considering the seabed slope over a distance of 10 deepwater wavelengths in the upwave \\
      direction from the turbine position and the average seabed slope within each one-wavelength \\
      segment, each of the ten values of seabed slope should not differ by more than 0.02 rad \\
      from from the slope calculated above, and no value of slope should differ by more than \\
      0.005 rad from the slopes of adjacent segments.
    
    * The current velocity should be less than 1.5 m/s.

    Parameters
    ----------
    Hm0 : float
        The spectral wave height [m].
    d : float
        The water depth [m].
    slope : float or tuple
        The slope of the foreshore [rad], or as tuple formatted as (V,H).
    k1 : float, optional
        Shape parameter of the first part of the Composite Weibull Distribution (CWD) before Htr [-], by default 2.
    k2 : float, optional
        Shape parameter of the second part of the Composite Weibull Distribution (CWD) after Htr [-], by default 3.6.

    Attributes
    ----------
    Hm0 : float
        Spectral wave height [m].
    d : float
        Water depth [m].
    slope : float
        The slope of the foreshore [rad].
    Hrms : float
        The root-mean-square wave height [m]. Determined from Hm0 as :math:`(0.6725 + \\frac{0.2025H_{m0}}{d})H_{m0}`.
    Htr : float
        The transitional wave height [m].
    Htr_tilde : float
        The transitional wave height, made dimensionless with Hrms [-].
    k1 : float
        Shape parameter of the first part of the Composite Weibull Distribution (CWD) before Htr [-].
    k2 : float
        Shape parameter of the second part of the Composite Weibull Distribution (CWD) after Htr [-].
    H1 : float
        Scale parameter of the first part of the Composite Weibull Distribution (CWD) before Htr [m].
    H2 : float
        Scale parameter of the second part of the Composite Weibull Distribution (CWD) after Htr [m].
    H1_tilde : float
        Scale parameter of the first part of the Composite Weibull Distribution (CWD) before Htr, made dimensionless with Hrms [m].
    H2_tilde : float
        Scale parameter of the second part of the Composite Weibull Distribution (CWD) after Htr, made dimensionless with Hrms [m].
    
    References
    ----------
    .. [1] Battjes, J.A., and H.W. Groenendijk. 2000. Wave height distributions on shallow foreshores, Coastal Engineering, 40, 161-182.
        <https://www.researchgate.net/publication/228446847_Wave_height_distributions_on_shallow_foreshores>
    .. [2] IEC 61400-3-1: Wind Energy Generation Systems – Part 3-1: Design Requirements for Fixed Offshore Wind Turbines, edition 4.0, 2019-02.
        <https://webstore.iec.ch/publication/26423>
    """

    def __init__(self, Hm0: float, d: float, slope: float | tuple[float, float], k1: float = 2, k2: float = 3.6):
        """ See help(BattjesGroenendijk) for more info """
        # check if slope foreshore is given as a tuple (V, H)
        if isinstance(slope, tuple):
            # if so, compute the angle
            slope = np.arctan(slope[0]/slope[1])
        
        self.Hm0 = Hm0
        self.d = d
        self.slope = slope
        self.Htr = (0.35 + 5.8*np.tan(slope))*d
        self.Hrms = (0.6725 + 0.2025*(Hm0/d))*Hm0
        self.Htr_tilde = self.Htr/self.Hrms

        self.k1 = k1
        self.k2 = k2
        
        self.H1_tilde, self.H2_tilde = fsolve(self._solver, [1, 1])
        self.H1 = self.H1_tilde * self.Hrms
        self.H2 = self.H2_tilde * self.Hrms

    @staticmethod
    def gammainc_upper(a: float, x: float):
        """
        Upper incomplete gamma function.

        Defined as:

        .. math::
            \Gamma(a,x) = \int_x^\infty t^{a - 1}e^{-t} dt

        Implemented using :mod:`scipy.special` as gamma(a)*gammaincc(a,x).

        Parameters
        ----------
        a : float
            positive parameter
        x : float
            nonnegative argument

        Returns
        -------
        float
            Value of the upper incomplete gamma function.
        
        See Also
        --------
        scipy.special.gammaincc: Regularized upper incomplete gamma function.
        scipy.special.gamma : Gamma function.
        """
        if a >= 0:
            out = sc.gamma(a)*sc.gammaincc(a, x)
            return out
        else:
            return None

    @staticmethod
    def gammainc_lower(a, x):
        """
        Lower incomplete gamma function.

        Defined as:

        .. math::
            \Gamma(a,x) = \int_0^x t^{a - 1}e^{-t} dt

        Implemented using :mod:`scipy.special` as gamma(a)*gammainc(a,x).

        Parameters
        ----------
        a : float
            positive parameter
        x : float
            nonnegative argument

        Returns
        -------
        float
            value of the lower incomplete gamma function
        
        See Also
        --------
        scipy.special.gammainc : Regularized lower incomplete gamma function.
        scipy.special.gamma : gamma function.
        """
        if a >= 0:
            out = sc.gamma(a)*sc.gammainc(a, x)
            return out
        else:
            return None

    def _solver(self, x):
        """Computes H1~ and H2~ for a given Htr~

        Computes H1~ and H2~ by solving equation 7 from Battjes and
        Groenendijk (2000), and the continuity condition for the
        Composite Weibull Distribution. This continuity conditions is:
        .. math::
            \frac{H_tr}{H_1}^{k_1} = \frac{H_tr}{H_2}^{k_2}

        The values :math:`H_{tr}`, :math:`H_{1}` and :math:`H_{2}` are
        then made nondimensional by dividing them by :math:`H_{rms}`. A
        tilde is added to the names of the values to indicate that they
        are the nondimensional values.

        Parameters
        ----------
        x : list
            initial guess for the values of H1~ and H2~, [float, float]

        Returns
        -------
        list
            Containing the continuity equation and equation 7 of
            Battjes and Groenendijk (2000) used to find the values of
            H1~ and H2~ for a given Htr~
        """
        a1 = 2/self.k1 + 1
        a2 = 2/self.k2 + 1

        eq_continuity = (self.Htr_tilde
                         * (self.Htr_tilde/x[0])**(-self.k1/self.k2) - x[1])

        eq7 = (np.sqrt(x[0]**2
                       * self.gammainc_lower(a1,(self.Htr_tilde/x[0])**self.k1)
                       + x[1]**2
                       * self.gammainc_upper(a2,(self.Htr_tilde/x[1])**self.k2))
               - 1)

        return [eq_continuity, eq7]

    def get_Hp(self, P: float, converge_to_rayleigh: bool=False):
        """
        Computes :math:`H_{P}`.

        This method computes the wave height exceeded by P% of the
        waves.
        
        The wave height estimate can be made to converge to the Rayleigh estimate for deeper water
        depths by setting `converge_to_rayleigh` to True. The criteria for convergence is that the
        Battjes estimate must not exceed the Rayleigh estimate, and that Htr_tilde >= 2.75.

        Parameters
        ----------
        P : float
            Exceedance probability as a fractal [-].
        converge_to_rayleigh : bool
            Determines if the model shuold converge to the Rayleigh distribution for deeper water depths.

        Returns
        -------
        H_P : float
            Wave height exceeded by P\\% of the waves.

        Raises
        ------
        ValueError
            If P is entered as a percentage, must be entered as P%/100.
        
        References
        ----------
        .. [1] Caires, Sofia & Gent, Marcel. (2010). Wave height distribution in constant and finite depths. Coastal Engineering Proceedings. 1. 10.9753/icce.v33.waves.15.
            <https://www.researchgate.net/publication/259260987_Wave_height_distribution_in_constant_and_finite_depths>
        """
        if P >= 1:
            raise ValueError('Do not enter P as a percentage, but P%/100')

        H_CWD_1 = self.H1_tilde*(-np.log(P))**(1/self.k1)

        if H_CWD_1 <= self.Htr_tilde:
            hp = H_CWD_1*self.Hrms
        else:
            H_CWD_2 = self.H2_tilde*(-np.log(P))**(1/self.k2)
            hp = H_CWD_2*self.Hrms
        
        if converge_to_rayleigh and (hp < rayleigh_hp(Hm0=self.Hm0, P=P)) and (self.Htr_tilde >= 2.75):
            return rayleigh_hp(Hm0=self.Hm0, P=P)
        
        return hp

    def get_Hn(self, N):
        """
        Computes :math:`H_{1/N}`.

        This method computes the mean of the highest 1/N-part of the
        wave heights.

        Parameters
        ----------
        N : float
            Highest 1/N-part of the wave heights [-].

        Returns
        -------
        H_1/N : float
            Mean of the highest 1/N part of the waves.

        Raises
        ------
        ValueError
            If N is smaller than 1.
        """
        if N >= 1:

            a1 = 1/self.k1 + 1
            a2 = 1/self.k2 + 1

            H_N = self.H1_tilde*(np.log(N))**(1/self.k1)

            if self.Htr_tilde >= H_N:
                H_1 = (N*self.H1_tilde*(self.gammainc_upper(a1, np.log(N))
                       - self.gammainc_upper(a1, (self.Htr_tilde/self.H1_tilde)**self.k1))
                       + N*self.H2_tilde
                       * self.gammainc_upper(a2, (self.Htr_tilde/self.H2_tilde)**self.k2))

                return H_1*self.Hrms
            else:
                H_2 = N*self.H2_tilde*self.gammainc_upper(a2, np.log(N))

                return H_2*self.Hrms
        else:
            msg = f'N = {N} has no meaning, please enter a value N >= 1'
            raise ValueError(msg)
    
    def plot_poe(self, ax: Axes | None = None, figsize: tuple = (8, 6)) -> tuple[Figure, Axes]:
        """
        Plot the probability of exceedance of the distribution.

        Parameters
        ----------
        ax : matplotlib.axes.Axes | None, optional
            Axes on which the Layout is drawn. If None, a new Figure and Axes are created, by default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6)

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.

        Raises
        ------
        TypeError
            The passed ax argument was not of type Axes.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            ax.invert_xaxis()
            ax.set_xscale('log')
            ax.set_xticks([100, 40, 30, 20, 10, 5, 2, 1, 0.5, 0.1, 0.01])
            ax.set_xticklabels([100, 40, 30, 20, 10, 5, 2, 1, 0.5, 0.1, 0.01])
            ax.xaxis.grid(which='major')
            ax.yaxis.grid(which='major')
            ax.set_xlabel("Probability of Exceedance [%]")
            ax.set_ylabel("Wave Height [m]")
        else:
            try:
                fig = ax.get_figure()
            except AttributeError:
                raise TypeError from AttributeError(
                    f"Invalid type in {type(ax)} for the 'ax' argument, "
                    f"must be Matplotlib Axes object."
                )
        h001 = self.get_Hp(0.1/100)
        
        h_array = np.linspace(0, h001, 100) / self.Hrms
        
        probs = []
        for h in h_array:
            if h <= self.Htr_tilde:
                p = np.exp(-(h/self.H1_tilde)**self.k1)
                probs.append(p)
            elif h >= self.Htr_tilde:
                p = np.exp(-(h/self.H2_tilde)**self.k2)
                probs.append(p)
        
        probs = np.array(probs)
        
        ax.plot(probs*100, h_array*self.Hrms, label=f"Battjes Hm0: {self.Hm0}m, d:{self.d}m, slope: {self.slope:.2f}")
        ax.legend()
        return fig, ax
    
    def plot_pdf(self, normalize: bool = False,  ax: Axes | None = None, figsize: tuple = (8, 6)) -> tuple[Figure, Axes]:
        """
        Plot the probability density function of the distribution.

        Parameters
        ----------
        normalize : bool, optional
            When True, the wave heights are mormalized by Hrms, by default False.
        ax : matplotlib.axes.Axes | None, optional
            Axes on which the Layout is drawn. If None, a new Figure and Axes are created, by default None.
        figsize : tuple, optional
            Figure size in inches in (width, height) format, by default (8, 6).

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.

        Raises
        ------
        TypeError
            The passed ax argument was not of type Axes.
        """
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            ax.grid()
            if normalize:
                ax.set_xlabel("H/Hrms [-]")
                ax.set_ylabel("Probability Density Function f(H/Hrms)")
            else:
                ax.set_xlabel("H [m]")
                ax.set_ylabel("Probability Density Function f(H)")
        else:
            try:
                fig = ax.get_figure()
            except AttributeError:
                raise TypeError from AttributeError(
                    f"Invalid type in {type(ax)} for the 'ax' argument, "
                    f"must be Matplotlib Axes object."
                )
        if normalize:
            heights = np.linspace(0, 4, 100)
            probs = []
            
            for height in heights:
                if height <= self.Htr_tilde:
                    prob = ((self.k1 * height ) / (self.H1_tilde**self.k1)) * np.exp(-(height/self.H1_tilde)**self.k1)
                    probs.append(prob)
                elif height >= self.Htr_tilde:
                    prob = ((self.k2 * height ) / (self.H2_tilde**self.k2)) * np.exp(-(height/self.H2_tilde)**self.k2)
                    probs.append(prob)
        else:
            heights = np.linspace(0, 4*self.Hrms, 100)
            probs = []
            
            for height in heights:
                if height <= self.Htr:
                    prob = ((self.k1 * height ) / (self.H1**self.k1)) * np.exp(-(height/self.H1)**self.k1)
                    probs.append(prob)
                elif height >= self.Htr:
                    prob = ((self.k2 * height ) / (self.H2**self.k2)) * np.exp(-(height/self.H2)**self.k2)
                    probs.append(prob)
        
        probs = np.array(probs)
        if normalize:
            ax.plot(heights, probs, label=f"Hm0: {self.Hm0}m, d:{self.d}m, slope: {self.slope:.2f}, Htr: {self.Htr_tilde:.2f}")
        else:
            ax.plot(heights, probs, label=f"Hm0: {self.Hm0}m, d:{self.d}m, slope: {self.slope:.2f}, Htr: {self.Htr:.2f}m")
        ax.legend()
        
        return fig, ax
