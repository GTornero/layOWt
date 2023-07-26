from coastlib.waves import FentonWave
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .wave_functions import ursell, williams_breaking_criterion


class Wave(FentonWave):
    """
    This class provides interface to John Fenton's stream function wave theory program Fourier.
    See http://johndfenton.com/Steady-waves/Fourier.html
    
    .. image:: ../../../../_static/fentonwave.png
        :width: 800
        :align: center

    Parameters
    ----------
    wave_height : float
        Wave height in meters.
    depth : float
        Water depth in meters.
    wave_period : float, optional
        Wave period in seconds. Mandatory, if wave_length is not provided.
    wave_length : float, optional
        Wave length in meters. Mandatory, if wave_period is not provided.
    timeout : float, optional
        Timeout in seconds (default=5). Fourier.exe is terminated after timeout.
    kwargs
        g : float, optional
            Standard acceleration due to gravity in m/s^2 (default=9.80665).
        rho : float, optional
            Water density in kg/m^3 (default=1025).
        current_criterion : int, optional
            1 for Eulerian mean current, 2 for mass-transport velocity (default=1).
        current_velocity : float, optional
            Current velocity in meters per second (default=0).
        fourier_components : int, optional
            Number of Fourier components.
            10-20 OK for ordinary waves, more for highest (default=20).
        height_steps : int, optional
            Maximum number of iterations for each height step.
            10 OK for ordinary waves, 40 for highest (default=10).
        convergence : dict, optional
            Maximum number of iterations for each height step; 10 OK for ordinary waves, 40 for highest.
            Criterion for convergence, typically '1.e4' or '1.e5' for highest waves.
            default=dict(max_iter=40, crit_conv='1.e-5').
        points : dict, optional
            Number of points on free surface.
            Number of velocity/acceration profiles over half a wavelength, including 0 and lambda/2.
            Number of vertical points in each profile, including points at bottom and surface.
            default=dict(n_surface=10, n_profiles=100, n_vertical=50).

    Examples
    --------
    >>> wave = Wave(wave_height=2, wave_period=6, depth=20)
    >>> wave.wave_length
    55.70521
    >>> wave = Wave(wave_height=2, wave_length=55.70521, depth=20)
    >>> wave.wave_period
    6.0000000594928
    >>> fig, ax = wave.plot(what='ua', scale=1, nprof=4)
    
    References
    ----------
    .. [1] Use of the programs FOURIER, CNOIDAL and STOKES for steady waves
        <https://johndfenton.com/Steady-waves/Instructions.pdf>
    .. [2]  Fenton, J. D. & Rienecker, M. M. (1982) A Fourier method for solving nonlinear water wave problems: Application to solitary wave interactions, J. Fluid Mechanics 118, 411-443.
        <https://johndfenton.com/Papers/Fenton82c+Rienecker-A-Fourier-method-for-solving-nonlinear-water-wave-problems.pdf>
    .. [3] Rienecker, M. M. & Fenton, J. D. (1981) A Fourier approximation method for steady water waves, J. Fluid Mechanics 104, 119-137.
        <https://johndfenton.com/Papers/Rienecker81+Fenton-A-Fourier-approximation-method-for-steady-water-waves.pdf>
    """
    def __init__(self, wave_height, depth, wave_period=None, wave_length=None, timeout=5, **kwargs):
        super().__init__(wave_height=wave_height, depth=depth, wave_period=wave_period, wave_length=wave_length, timeout=timeout, **kwargs)
    
    def validate_theory(self) -> tuple[Figure, Axes]:
        """
        Shows what wave theory is applicable for given wave parameters using the
        wave theories' figure by Le Mehaute per USACE CEM Part II Chap. 1 p.II-1-58

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.
        
        See Also
        --------
        layowt.metocean.wave_functions.validate_theory : Shows what wave theory is applicable for given wave parameters using the wave theories' figure by Le Mehaute.
        
        References
        ----------
        .. [1] Part 2 of Coastal engineering manual, USA Corps of Engineers Volume 1110, Issues 2-1100 of Engineer manual, USA Corps of Engineers
            <https://www.plainwater.com/pubs/EM-1110-2-1100P2.pdf>
        """
        return super().validate()
    
    def williams_breaking_criterion(self, ax: Axes | None = None, figsize: tuple[float, float] = (10,10)) -> tuple[Figure, Axes]:
        """
        Plots the Williams (1981) analytical breaking wave criterion for the Wave instance.
        
        .. image:: ../../../../_static/williams_breaking_criterion.png
            :width: 800
            :align: center

        Parameters
        ----------
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
        
        References
        ----------
        .. [1] Use of the programs FOURIER, CNOIDAL and STOKES for steady waves
            <https://johndfenton.com/Steady-waves/Instructions.pdf>
        .. [2]  Fenton, J. D. & Rienecker, M. M. (1982) A Fourier method for solving nonlinear water wave problems: Application to solitary wave interactions, J. Fluid Mechanics 118, 411-443.
            <https://johndfenton.com/Papers/Fenton82c+Rienecker-A-Fourier-method-for-solving-nonlinear-water-wave-problems.pdf>
        .. [3] Rienecker, M. M. & Fenton, J. D. (1981) A Fourier approximation method for steady water waves, J. Fluid Mechanics 104, 119-137.
            <https://johndfenton.com/Papers/Rienecker81+Fenton-A-Fourier-approximation-method-for-steady-water-waves.pdf>
        .. [4] Williams, J.M. (1981) ”Limiting gravity waves in water of finite depth”, Phil. Trans. Roy. Soc. London A 302, 139-188.
        """
        return williams_breaking_criterion(wave_height=self.wave_height, wave_length=self.wave_length, depth=self.depth, ax=ax, figsize=figsize)
    
    @property
    def ursell(self) -> float:
        """
        Returns the URsell number of the wave.

        Returns
        -------
        float
            Ursell number.
        
        See Also
        --------
        layowt.metocean.wave_functions.ursell : Calculates the Ursell number.
        """
        return ursell(wave_height=self.wave_height, wave_length=self.wave_length, depth=self.depth)
    
    @property
    def crest_height(self) -> float:
        """
        Returns the crest height of the wave.

        Returns
        -------
        float
            Crest height above the still water level in metres.
        """
        return self.surface['eta (m)'].max() - self.depth
    
    @property
    def trough_height(self) -> float:
        """
        Returns the trough height of the wave.

        Returns
        -------
        float
            Trough height above the still water level in metres.
        """
        return self.surface['eta (m)'].min() - self.depth
    
    def plot(self, what: str = 'ua', scale: float = 1, nprof: int = 4) -> tuple[Figure, Axes]: 
        """
        Plots wave profile with multiple depth profiles for a given parameter.

        Parameters
        ----------
        what : str, optional
            Parameter to plot, by default 'ua'. \n
                'u': 'u (m/s)', \n
                'v': 'v (m/s)', \n
                'ua': 'du/dt (m/s^2)', \n
                'va': 'dv/dt (m/s^2)', \n
                'ux': 'du/dx (1/s)', \n
                'vx': 'du/dy (1/s)', \n
                'uy': 'du/dy (1/s)', \n
                'vy': 'du/dx (1/s)', \n
                'phi': 'dphi/dt (m^2/s^2)'
        scale : float, optional
            Parameter scale, by default 1
        nprof : int, optional
             Number of parameter profiles, by default 4

        Returns
        -------
        tuple[matplotlib.figure.Figure, matplotlib.axes.Axes]
            Tuple of Figure and Axes objects.
        """
        return super().plot(what=what, scale=scale, nprof=nprof)