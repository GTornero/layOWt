from dataclasses import asdict, dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress

from layowt.metocean.wind_functions import DNV_ntm, ntm


@dataclass
class BaseData:
    """
    Base class for all Data dataclasses.

    Expected child attributes:
        `data` : pd.DataFrame -- data
        `columns`: object -- dataclass with column names

    Raises
    ------
    ValueError
        If columns in `columns` and `data` dataframe don't match
    """
    data: pd.DataFrame
    columns: object

    def __post_init__(self):
        self._check_columns()
    
    def _check_columns(self):
        """
        Check if columns in dataframe match the columns in `columns` attribute

        Raises
        ------
        ValueError
            If columns don't match
        """
        data_columns = set(self.data.columns)
        schema_columns = {v for k,v in asdict(self.columns).items()}
        if not schema_columns.issubset(data_columns):
            raise ValueError("Missing Columns")
    
@dataclass
class WindDataFields:
    """
    data field names aka columns
    """
    wind_speed  : str = 'wind_speed'
    wind_dir    : str = 'wind_dir'
    std_dev     : str = 'std_dev'

@dataclass
class WindData(BaseData):
    data: pd.DataFrame
    columns: WindDataFields = WindDataFields()

class NTM:
    """
    Normal Turbulence Model Class used to analyze the turbulence intensity of wind dataset in \
    accordance with the definitions from IEC 61400-1 edition 4.0, 2019-02.

    Parameters
    ----------
    data : pd.DataFrame
        Input wind data. Must contain wind speed, wind direction, and standard deviation in columns \
        names 'wind_speed', 'wind_dir', and 'std_dev'.
    ws_bin_size : float, optional
        The size of the wind speed bins in m/s, by default 1.
    n_sectors : int, optional
        Number of direction sectors when carrying out directional analysis, by default 12.
    
    Attributes
    ----------
    data : pd.Dataframe
        Input wind data
    ws_bin_size : float
        The size of the wind speed bins in m/s.
    n_sectors : int
        Number of direction sectors.
    model : str
        The type of NTM model used (Offshore or Onshore).
    sigma_mode : str
        The type of standard deviation analysis (omnidirectional or directional).
    sigma_sigma_mode : str
        The type of standard deviation of wind speed standard deviation analysis (omnidirectional or directional).
    min_sigma_points : int
        Minimum number of points in a wind speed bin to use in the standard deviation model fitting.
    min_sigma_sigma_points : int
        Minimum number of points in a wind speed bin to use in the standard deviation of the standard deviation model fitting.
    sigma_samples : pd.DataFrame
        All of the mean standard deviation values per wind speed bin in the data.
    sigma_used_samples : pd.DataFrame
        The mean standard deviation values per wind speed bin used in the model fit.
    sigma_fit : list
        Coefficients of the standard deviation profiles of the Normal Turbulence Model.
    sigma_model : pd.DataFrame
        The standard deviation profile of the fitted Normal Turbulence Model.
    sigma_sigma_samples : pd.DataFrame
        All of the mean standard deviation of the standard deviation values per wind speed bin in the data.
    sigma_sigma_used_samples : pd.DataFrame
        The mean standard deviation of the standard deviation values per wind speed bin used in the model fit.
    sigma_sigma_fit : list
        Coefficients of the standard deviation of the standard deviation profiles of the Normal Turbulence Model.
    sigma_sigma_model : pd.DataFrame
        The standard deviation of the standard deviation profile of the fitted Normal Turbulence Model.
    overwrite_sigma : bool
        If the fitted mean standard deviation values of the Normal Turbulence Model have been overwritten by the measurements.
    """
    def __init__(self, data: pd.DataFrame, ws_bin_size: float = 1, n_sectors: int = 12) -> None:
        self.data = WindData(data)
        self.ws_bin_size = ws_bin_size
        self.n_sectors = n_sectors
        self._bin_wind_speed()
        self._deg_to_sectors()
        self.model = None
        self.sigma_mode = None
        self.sigma_sigma_mode = None
        self.min_sigma_points = None
        self.min_sigma_sigma_points = None
        self.sigma_samples = None
        self.sigma_used_samples = None
        self.sigma_fit = None
        self.sigma_model = None
        self.sigma_sigma_samples = None
        self.sigma_sigma_used_samples = None
        self.sigma_sigma_fit = None
        self.sigma_sigma_model = None
        self.overwrite_sigma = False
        self.ntm = None
    
    def _bin_wind_speed(self):
        """
        Adds the 'wind_speed_bin' column to the data using the defined wind speed bin size.

        Returns
        -------
        None
            None.
        """
        self.data.data['wind_speed_bin'] = np.round(self.data.data.wind_speed / self.ws_bin_size) * self.ws_bin_size
        return None
    
    def _deg_to_sectors(self):
        """
        Adds the 'sector' column to the data using the defined number of wind direction sectors.

        Returns
        -------
        None
            None
        """
        self.data.data['sector'] = np.floor((self.data.data.wind_dir + (360/self.n_sectors/2)) / (360 / self.n_sectors)) + 1
        self.data.data.loc[self.data.data.sector > self.n_sectors, 'sector'] = 1
        self.data.data.sector = self.data.data.sector.astype(int)
        return None
    
    def fit_model(
        self,
        sigma_mode: str = 'omni',
        sigma_sigma_mode: str = 'omni',
        model: str = 'offshore',
        min_sigma_points: int = 10,
        min_sigma_sigma_points: int = 50,
        overwrite_sigma: bool = True
        ):
        """
        Fit the Normal Turbulence Model to the data.
        
        [THEORY HERE]

        Parameters
        ----------
        sigma_mode : str, optional
            How to analyze the standard deviation, either 'omni' for and omnidirectional NTM model
            or 'directional' for a directional NTM model used the defined number of direction
            sectors, by default 'omni'.
        sigma_sigma_mode : str, optional
            How to analyze the standard deviation of the standard deviation, either 'omni' for and
            omnidirectional NTM model or 'directional' for a directional NTM model used the defined
            number of direction sectors. Directional analysis can only be done when sigma_mode is
            also directional. By default 'omni'.
        model : str, optional
            The type of NTM model to be used between the classic 'onshore' model as defined by the
            IEC, or the 'offshore' model which adds a second order effect to the standard deviation
            fit to account for the Charnock effect. By default 'offshore'
        min_sigma_points : int, optional
            Minimum number of points required to include the wind speed bin in the standard
            deviation model fit, by default 10.
        min_sigma_sigma_points : int, optional
            Minimum number of points required to include the wind speed bin in the standard
            deviation of the standard deviation model fit, by default 50.
        overwrite_sigma : bool, optional
            If the standard deviation profile values should be replaced by observations, by default True.

        Returns
        -------
        None
            None

        Raises
        ------
        ValueError
            _description_
        ValueError
            _description_
        ValueError
            _description_
        
        References
        ----------
        .. [1] IEC 61400-1: Wind Energy Generation Systems – Part 1: Design Requirements, edition 4.0, 2019-02.
            <https://webstore.iec.ch/publication/26423>
        """
        if (sigma_mode.lower() not in ('omni', 'directional')) or (sigma_sigma_mode.lower() not in ('omni', 'directional')):
            raise ValueError("Invalid mode passed. Must be one of 'omni' or 'directional'.")
        elif model.lower() not in ('offshore', 'onshore'):
            raise ValueError("Invalid model passed. Must be one of 'offshore' or 'onshore'.")
        
        if (sigma_mode.lower() == 'omni') and sigma_sigma_mode.lower() =='directional':
            raise ValueError("sigma_sigma analysis cannot be directional when sigma is omnidirectional.") 
        
        self._fit_sigma(mode=sigma_mode, model=model, min_points=min_sigma_points, overwrite=overwrite_sigma)
        self._fit_sigma_sigma(mode=sigma_sigma_mode, min_points=min_sigma_sigma_points)
        
        self.model = model
        self.sigma_mode = sigma_mode.lower()
        self.sigma_sigma_mode = sigma_sigma_mode.lower()
        self.min_sigma_points = min_sigma_points
        self.min_sigma_sigma_points = min_sigma_sigma_points
        self.overwrite_sigma = overwrite_sigma
        
        if self.sigma_mode.lower() =='omni':
            self.ntm = pd.DataFrame(self.sigma_model['std_dev'] + 1.28 * self.sigma_sigma_model['sigma_sigma'], columns=['ntm_sigma'])
        else:
            if self.sigma_sigma_mode.lower() == 'omni':
                self.ntm = self.sigma_model.add(1.28 * self.sigma_sigma_model['sigma_sigma'], axis='index')
            else:
                self.ntm = self.sigma_model.add(1.28*self.sigma_sigma_model, axis='columns')
        return None
    
    def _fit_sigma(self, mode: str = 'omni', model: str = 'offshore', min_points: int = 10, overwrite: bool = True):
        if mode.lower() == 'omni':
            self._fit_omni_sigma(model=model.lower(), min_points=min_points, overwrite=overwrite)
        else:
            self._fit_dir_sigma(model=model.lower(), min_points=min_points, overwrite=overwrite)
    
    def _fit_omni_sigma(self, model: str, min_points: int, overwrite: bool):
        sigma_profile = self.data.data.groupby('wind_speed_bin').mean()
        sigma_profile['count'] = self.data.data.groupby('wind_speed_bin').count()['wind_speed']
        used_sigma_samples = sigma_profile.loc[sigma_profile['count'] >= min_points]
        
        self.sigma_samples = sigma_profile.loc[:, 'std_dev']
        self.sigma_used_samples = used_sigma_samples.loc[:, 'std_dev']
        
        self.sigma_fit = self._fit(data=used_sigma_samples.loc[:, 'std_dev'], model=model)
        
        fit_x = np.arange(self.ws_bin_size, self.sigma_samples.index.max()+5, self.ws_bin_size)
         
        if len(self.sigma_fit) == 3:
            sigma_curve = self.sigma_fit[0] * fit_x**2 + self.sigma_fit[1] * fit_x + self.sigma_fit[2]
        else:
            sigma_curve = self.sigma_fit[0] * fit_x + self.sigma_fit[1]
        
        self.sigma_model = pd.DataFrame(data=sigma_curve, index=fit_x)
        self.sigma_model.index.name = 'wind_speed_bin'
        self.sigma_model.columns = ['std_dev']
        if overwrite:
            self.sigma_model.loc[self.sigma_used_samples.index, 'std_dev'] = self.sigma_used_samples.values
    
    def _fit_dir_sigma(self, model: str, min_points: int, overwrite: bool):
        sigma_samples = []
        used_sigma_samples = []
        fits = []
        
        sigma_profiles = pd.pivot_table(data=self.data.data, index='wind_speed_bin', columns='sector', values='std_dev', aggfunc='mean')
        sigma_counts = pd.pivot_table(data=self.data.data, index='wind_speed_bin', columns='sector', values='std_dev', aggfunc='count')
        for sector in range(1, self.n_sectors+1):
            sigma_samples.append(sigma_profiles.loc[:, sector])
            used_sigma_samples.append(sigma_profiles.loc[sigma_counts[sector]>=min_points, sector])
            
            fits.append(self._fit(data=used_sigma_samples[-1], model=model))
        
        self.sigma_samples = sigma_samples
        self.sigma_used_samples = used_sigma_samples
        self.sigma_fit = fits
        
        curves = []
        fit_x = np.arange(self.ws_bin_size, self.data.data['wind_speed_bin'].max()+5, self.ws_bin_size)
        for i, fit in enumerate(self.sigma_fit):
            if len(fit) == 3:
                sigma_curve = fit[0] * fit_x**2 + fit[1] * fit_x + fit[2]
            else:
                sigma_curve = fit[0] * fit_x + fit[1]
                
            sigma_curve = pd.Series(data=sigma_curve, index=fit_x, name=i)
            if overwrite:
                sigma_curve.loc[self.sigma_used_samples[i].index] = self.sigma_used_samples[i].values
            curves.append(sigma_curve)
        
        self.sigma_model = pd.concat(curves, axis=1)
        self.sigma_model.index.name = 'wind_speed_bin'
    
    @staticmethod
    def _fit(data, model: str):
        if model == 'offshore':
            fit = np.polyfit(x=data.index, y=data, deg=2)
            if fit[0] < 0:
                fit = np.polyfit(x=data.index, y=data, deg=1)
        else:
            fit = np.polyfit(x=data.index, y=data, deg=1)
        return fit
    
    def _fit_sigma_sigma(self, mode='omni', min_points: int = 50):
        if mode.lower() == 'omni':
            self._fit_omni_sigma_sigma(min_points=min_points)
        else:
            self._fit_dir_sigma_sigma(min_points=min_points)
    
    def  _fit_omni_sigma_sigma(self, min_points: int):
        sigma_profile = self.data.data.groupby('wind_speed_bin').std()
        sigma_profile['count'] = self.data.data.groupby('wind_speed_bin').count()['wind_speed']
        used_sigma_samples = sigma_profile.loc[sigma_profile['count'] >= min_points]
        
        self.sigma_sigma_samples = sigma_profile.loc[:, 'std_dev']
        self.sigma_sigma_used_samples = used_sigma_samples.loc[:, 'std_dev']
        
        self.sigma_sigma_fit = self._fit(data=used_sigma_samples.loc[:, 'std_dev'], model='onshore')
        
        fit_x = np.arange(self.ws_bin_size, self.sigma_sigma_samples.index.max()+5, self.ws_bin_size)
        sigma_curve = self.sigma_sigma_fit[0] * fit_x + self.sigma_sigma_fit[1]
        
        self.sigma_sigma_model = pd.DataFrame(data=sigma_curve, index=fit_x)
        self.sigma_sigma_model.index.name = 'wind_speed_bin'
        self.sigma_sigma_model.columns = ['sigma_sigma']
    
    def _fit_dir_sigma_sigma(self, min_points: int):
        sigma_samples = []
        used_sigma_samples = []
        fits = []
        
        sigma_profiles = pd.pivot_table(data=self.data.data, index='wind_speed_bin', columns='sector', values='std_dev', aggfunc='std')
        sigma_counts = pd.pivot_table(data=self.data.data, index='wind_speed_bin', columns='sector', values='std_dev', aggfunc='count')
        for sector in range(1, self.n_sectors+1):
            sigma_samples.append(sigma_profiles.loc[:, sector])
            used_sigma_samples.append(sigma_profiles.loc[sigma_counts[sector]>=min_points, sector])
            
            fits.append(self._fit(data=used_sigma_samples[-1], model='onshore'))
        
        self.sigma_sigma_samples = sigma_samples
        self.sigma_sigma_used_samples = used_sigma_samples
        self.sigma_sigma_fit = fits
        
        curves = []
        fit_x = np.arange(self.ws_bin_size, self.data.data['wind_speed_bin'].max()+5, self.ws_bin_size)
        for i, fit in enumerate(self.sigma_sigma_fit):
            sigma_curve = fit[0] * fit_x + fit[1]
            sigma_curve = pd.Series(data=sigma_curve, index=fit_x, name=i)
            curves.append(sigma_curve)
        
        self.sigma_sigma_model = pd.concat(curves, axis=1)
        self.sigma_sigma_model.index.name = 'wind_speed_bin'
    
    def plot_sigma(self, figsize=(8, 8)):
        
        if self.sigma_mode.lower() == 'omni':
            fig, ax = plt.subplots(figsize=figsize)
            
            ax.scatter(self.sigma_samples.index, self.sigma_samples, c='None', edgecolors='k', label="Measured")
            ax.scatter(self.sigma_used_samples.index, self.sigma_used_samples, c='k', label="Used")
            
            fit_x = np.arange(self.ws_bin_size, self.sigma_samples.index.max()+5, self.ws_bin_size)
            if len(self.sigma_fit) == 3:
                ax.plot(fit_x, self.sigma_fit[0] * fit_x**2 + self.sigma_fit[1] * fit_x + self.sigma_fit[2], label=f"Second-Order Fit: $y = {self.sigma_fit[0]:.3f}x^2 + {self.sigma_fit[1]:.3f}x + {self.sigma_fit[2]:.3f}$", c='k')
            else:
                ax.plot(fit_x, self.sigma_fit[0] * fit_x + self.sigma_fit[1], label=f"Linear Fit: $y = {self.sigma_fit[0]:.3f}x + {self.sigma_fit[1]:.3f}$", c='k')
            
            ax.plot(self.sigma_model.index, self.sigma_model['std_dev'], linestyle='--', c='r', label='NTM Model')
            ax.grid()
            ax.legend()
            ax.set_xlabel("Horizontal Wind Speed [m/s]", fontweight="bold")
            ax.set_ylabel("Horizontal Wind Speed Standard Deviation [m/s]", fontweight="bold")
            ax.set_title("Omnidirectional Mean Turbulence Model Fitting", fontweight="bold")
            
            return fig, ax
        else:
            layout = {12: (4, 3),
                      16: (4, 4),
                      32: (8, 4),
                      36: (6, 6)}
            
            sector_titles = {}
            for n in range(self.n_sectors):
                if n == 0:
                    sector_titles[n] = f"{360-360/self.n_sectors/2:.2f} - {360/self.n_sectors/2:.2f}"
                else:
                    sector_titles[n] = f"{360/self.n_sectors/2 + 360/self.n_sectors*n:.2f} - {360/self.n_sectors/2 + 360/self.n_sectors*(n+1):.2f}"
            
            fig, axs = plt.subplots(*layout[self.n_sectors], figsize=figsize)
            for sector, ax in enumerate(axs.flatten()):
                ax.scatter(self.sigma_samples[sector].index, self.sigma_samples[sector], c='None', edgecolors='k', label="Measured")
                ax.scatter(self.sigma_used_samples[sector].index, self.sigma_used_samples[sector], c='k', label="Used")
                fit_x = np.arange(self.ws_bin_size, self.sigma_samples[sector].index.max()+5, self.ws_bin_size)
                if len(self.sigma_fit[sector]) == 3:
                    ax.plot(fit_x, self.sigma_fit[sector][0] * fit_x**2 + self.sigma_fit[sector][1] * fit_x + self.sigma_fit[sector][2], label=f"Second-Order Fit: $y = {self.sigma_fit[sector][0]:.3f}x^2 + {self.sigma_fit[sector][1]:.3f}x + {self.sigma_fit[sector][2]:.3f}$", c='k')
                else:
                    ax.plot(fit_x, self.sigma_fit[sector][0] * fit_x + self.sigma_fit[sector][1], label=f"Linear Fit: $y = {self.sigma_fit[sector][0]:.3f}x + {self.sigma_fit[sector][1]:.3f}$", c='k')
                
                ax.plot(self.sigma_model[sector].index, self.sigma_model[sector], linestyle='--', c='r', label='NTM Model')
                ax.grid()
                ax.legend()
                ax.set_title(f"Sector {sector+1}: " + sector_titles[sector])
                # Set common labels
                fig.text(0.5, 0.085, "Horizontal Wind Speed [m/s]", ha='center', va='center', fontsize=16)
                fig.text(0.1, 0.5, "Horizontal Wind Speed Standard Deviation [m/s]", ha='center', va='center', rotation='vertical', fontsize=16)
            
            return fig, axs
    
    def plot_sigma_sigma(self, figsize=(8, 8)):
        if self.sigma_sigma_mode.lower() == 'omni':
            fig, ax = plt.subplots(figsize=figsize)
            
            ax.scatter(self.sigma_sigma_samples.index, self.sigma_sigma_samples, c='None', edgecolors='k', label="Measured")
            ax.scatter(self.sigma_sigma_used_samples.index, self.sigma_sigma_used_samples, c='k', label="Used")
            
            fit_x = np.arange(self.ws_bin_size, self.sigma_sigma_samples.index.max()+5, self.ws_bin_size)
            if len(self.sigma_sigma_fit) == 3:
                ax.plot(fit_x, self.sigma_sigma_fit[0] * fit_x**2 + self.sigma_sigma_fit[1] * fit_x + self.sigma_sigma_fit[2], label=f"Second-Order Fit: $y = {self.sigma_sigma_fit[0]:.3f}x^2 + {self.sigma_sigma_fit[1]:.3f}x + {self.sigma_sigma_fit[2]:.3f}$", c='k')
            else:
                ax.plot(fit_x, self.sigma_sigma_fit[0] * fit_x + self.sigma_sigma_fit[1], label=f"Linear Fit: $y = {self.sigma_sigma_fit[0]:.3f}x + {self.sigma_sigma_fit[1]:.3f}$", c='k')
            ax.grid()
            ax.legend()
            ax.set_xlabel("Horizontal Wind Speed [m/s]", fontweight="bold")
            ax.set_ylabel("Horizontal Wind Speed $\\sigma_{\\sigma}$ [m/s]", fontweight="bold")
            ax.set_title("Omnidirectional Mean $\\sigma_{\\sigma}$ Model Fitting", fontweight="bold")
            
            return fig, ax
        else:
            layout = {12: (4, 3),
                      16: (4, 4),
                      32: (8, 4),
                      36: (6, 6)}
            
            sector_titles = {}
            for n in range(self.n_sectors):
                if n == 0:
                    sector_titles[n] = f"{360-360/self.n_sectors/2:.2f} - {360/self.n_sectors/2:.2f}"
                else:
                    sector_titles[n] = f"{360/self.n_sectors/2 + 360/self.n_sectors*n:.2f} - {360/self.n_sectors/2 + 360/self.n_sectors*(n+1):.2f}"
            
            fig, axs = plt.subplots(*layout[self.n_sectors], figsize=figsize)
            for sector, ax in enumerate(axs.flatten()):
                ax.scatter(self.sigma_sigma_samples[sector].index, self.sigma_sigma_samples[sector], c='None', edgecolors='k', label="Measured")
                ax.scatter(self.sigma_sigma_used_samples[sector].index, self.sigma_sigma_used_samples[sector], c='k', label="Used")
                fit_x = np.arange(self.ws_bin_size, self.sigma_sigma_samples[sector].index.max()+5, self.ws_bin_size)
                ax.plot(fit_x, self.sigma_sigma_fit[sector][0] * fit_x + self.sigma_sigma_fit[sector][1], label=f"Linear Fit: $y = {self.sigma_sigma_fit[sector][0]:.3f}x + {self.sigma_sigma_fit[sector][1]:.3f}$", c='k')
                
                ax.grid()
                ax.legend()
                ax.set_title(f"Sector {sector+1}: " + sector_titles[sector])
                
            # Set common labels
            fig.text(0.5, 0.085, "Horizontal Wind Speed [m/s]", ha='center', va='center', fontsize=16)
            fig.text(0.1, 0.5, "Horizontal Wind Speed $\\sigma_{\\sigma}$ [m/s]", ha='center', va='center', rotation='vertical', fontsize=16)
            
            return fig, axs
    
    def plot_ntm(self, show_measurements: bool = False, show_iec: bool = False, show_dnv: bool = False, hub_height: float | None = None, figsize=(8, 8)):
        if self.sigma_mode.lower() == 'omni':
            fig, ax = plt.subplots(figsize=figsize)
            
            ax.plot(self.sigma_model.index, self.sigma_model.div(self.sigma_model.index, axis=0), c='k', linestyle='--', label='Mean Ambient TI [-]')
            ax.plot(self.sigma_model.index, (self.sigma_sigma_model*1.28 + self.sigma_model.values).div(self.sigma_model.index, axis=0), c='r', linestyle='--', label='Characteristic Ambient TI [-]')
            
            if show_measurements:
                ax.scatter(self.sigma_samples.index, self.sigma_samples / self.sigma_samples.index, c='k', marker='x', s=12, label='Mean TI Measurements')
                ax.scatter(self.sigma_samples.index, (self.sigma_samples + 1.28 * self.sigma_sigma_samples) / self.sigma_samples.index, c='r', marker='x', s=12, label='Characteristic TI Measurements')
            
            if show_iec:
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.18, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class A+')
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.16, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class A') 
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.14, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class B')
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.12, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class C')
            
            if show_dnv:
                if hub_height is None:
                    raise ValueError("hub_height cannot be None is show_dnv is True.")
                
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], DNV_ntm(i_ref=0.14, ws=self.sigma_model.index[self.sigma_model.index >= 2.5], hub_height=hub_height, ta=10, tb=0.556, charnock=0.018), label='Class OA') 
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], DNV_ntm(i_ref=0.12, ws=self.sigma_model.index[self.sigma_model.index >= 2.5], hub_height=hub_height, ta=10.5, tb=0.561, charnock=0.014), label='Class OB')
                ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], DNV_ntm(i_ref=0.10, ws=self.sigma_model.index[self.sigma_model.index >= 2.5], hub_height=hub_height, ta=11, tb=0.566, charnock=0.011), label='Class OC')
                
            
            ax.grid()
            ax.legend()
            ax.set_xlabel("Horizontal Wind Speed [m/s]", fontweight="bold")
            ax.set_ylabel("Ambient Turbulence [-]", fontweight="bold")
            ax.set_title("Omnidirectional Ambient Turbulence Profile", fontweight="bold")
            
            return fig, ax
        else:
            layout = {
                12: (4, 3),
                16: (4, 4),
                32: (8, 4),
                36: (6, 6)}

            sector_titles = {}
            for n in range(self.n_sectors):
                if n == 0:
                    sector_titles[n] = f"{360-360/self.n_sectors/2:.2f} - {360/self.n_sectors/2:.2f}"
                else:
                    sector_titles[n] = f"{360/self.n_sectors/2 + 360/self.n_sectors*n:.2f} - {360/self.n_sectors/2 + 360/self.n_sectors*(n+1):.2f}"
            
            fig, axs = plt.subplots(*layout[self.n_sectors], figsize=figsize)
            for sector, ax in enumerate(axs.flatten()):
                ax.plot(self.sigma_model.index, self.sigma_model[sector].div(self.sigma_model.index, axis=0), c='k', linestyle='--', label='Mean Ambient TI [-]')
                if self.sigma_sigma_mode.lower() == 'omni':
                    ax.plot(self.sigma_model.index, ((self.sigma_sigma_model*1.28).values.flatten() + self.sigma_model[sector].values) / self.sigma_model.index.values, c='r', linestyle='--', label='Characteristic Ambient TI [-]')
                else:
                    ax.plot(self.sigma_model.index, ((self.sigma_sigma_model[sector]*1.28).values + self.sigma_model[sector].values) / self.sigma_model.index.values, c='r', linestyle='--', label='Characteristic Ambient TI [-]')
                
                if show_measurements:
                    ax.scatter(self.sigma_samples[sector].index, self.sigma_samples[sector].values.flatten() / self.sigma_samples[sector].index, c='k', marker='x', s=12, label='Mean TI Measurements')
                    if self.sigma_sigma_mode.lower() == 'omni':
                        ax.scatter(self.sigma_samples[sector].index, (self.sigma_samples[sector].values.flatten() + 1.28 * self.sigma_sigma_samples) / self.sigma_samples[sector].index, c='r', marker='x', s=12, label='Characteristic TI Measurements')
                    else:
                        ax.scatter(self.sigma_sigma_samples[sector].index, (self.sigma_samples[sector].loc[self.sigma_sigma_samples[sector].index].values.flatten() + 1.28 * self.sigma_sigma_samples[sector].values.flatten()) / self.sigma_sigma_samples[sector].index, c='r', marker='x', s=12, label='Characteristic TI Measurements')
                
                if show_iec:
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.18, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class A+')
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.16, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class A') 
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.14, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class B')
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], ntm(i_ref=0.12, ws=self.sigma_model.index[self.sigma_model.index >= 2.5]), label='Class C')
                
                if show_dnv:
                    if hub_height is None:
                        raise ValueError("hub_height cannot be None is show_dnv is True.")
                    
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], DNV_ntm(i_ref=0.14, ws=self.sigma_model.index[self.sigma_model.index >= 2.5], hub_height=hub_height, ta=10, tb=0.556, charnock=0.018), label='Class OA') 
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], DNV_ntm(i_ref=0.12, ws=self.sigma_model.index[self.sigma_model.index >= 2.5], hub_height=hub_height, ta=10.5, tb=0.561, charnock=0.014), label='Class OB')
                    ax.plot(self.sigma_model.index[self.sigma_model.index >= 2.5], DNV_ntm(i_ref=0.10, ws=self.sigma_model.index[self.sigma_model.index >= 2.5], hub_height=hub_height, ta=11, tb=0.566, charnock=0.011), label='Class OC')
                
                ax.grid()
                ax.legend()
                ax.set_title(f"Sector {sector+1}: " + sector_titles[sector])
            
            # Set common labels
            fig.text(0.5, 0.085, "Horizontal Wind Speed [m/s]", ha='center', va='center', fontsize=16)
            fig.text(0.1, 0.5, "Ambient Turbulence [-]", ha='center', va='center', rotation='vertical', fontsize=16)
            
            return fig, axs