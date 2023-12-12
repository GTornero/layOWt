% (CHANGELOG template from https://gist.github.com/juampynr/4c18214a8eb554084e21d6e288a18a2c)
% (How to use Git to create a changelog https://www.freecodecamp.org/news/a-beginners-guide-to-git-what-is-a-changelog-and-how-to-generate-it/)

# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

=======================================================================
## [0.7.0] - [TBC]
 
LayOWt Version 0.7.0

Addition of the Goda maximum wave model and the Battjes & Groenendijk, Rayleigh, Glukhovski, and Forristal wave height distributions.
Minor additions/fixes to existing subpackages.
 
### Added
- {class}`layowt.wtg.power_curve.PowerCurve`:
    - {meth}`layowt.wtg.power_curve.PowerCurve.compare_curves` graphically compares PowerCurve instance with other PowerCurve objects.
- {mod}`layowt.metocean.wind_functions`: Metocean wind functions module:
    - {func}`layowt.metocean.wind_functions.vref_ewts_exact` EWTS II Exact algorithm to estimate Vref from the mean wind speed and the Weibull k value.
    - {func}`layowt.metocean.wind_functions.vref_ewts_gumbel` EWTS II Gumbel algorithm to estimate Vref from the mean wind speed and the Weibull k value.
    - {func}`layowt.metocean.wind_functions.vref_ewts_davenport` EWTS II Davenport algorithm to estimate Vref from the mean wind speed and the Weibull k value.
    - {func}`layowt.metocean.wind_functions.froya` The Frøya wind profile adjustment.
    - {func}`layowt.metocean.wind_functions.ntm` IEC Normal Turbulence Model based on Iref.
    - {func}`layowt.metocean.wind_functions.etm` IEC Extreme Turbulence Model based on Iref.
- {mod}`layowt.metocean.wave_functions`: Metocean wave functions module:
    - {func}`layowt.metocean.wave_functions.rayleigh_hp`: Computes the wave height exceeded by P% of the waves, :math:`H_{P}` of a Rayleigh distribution.
    - {func}`layowt.metocean.wave_functions.rayleigh_hmaxn`: Computes the maximum wave height out of N individual waves given a Rayleigh distribution.
    - {func}`layowt.metocean.wave_functions.rayleigh_hn`: Computes the mean of the highest 1/N-part of the wave heights for the given a Rayleigh distribution.
    - {func}`layowt.metocean.wave_functions.plot_rayleigh_pdf`: Plots the probability density function of a Rayleigh distribution.
    - {func}`layowt.metocean.wave_functions.plot_rayleigh_poe`: Plots the probability of exceedance of a Rayleigh distribution.
    - {func}`layowt.metocean.wave_functions.glukhovskiy_hp`: Computes the wave height exceeded by P% of the waves, :math:`H_{P}` of a Glukhovskiy distribution.
    - {func}`layowt.metocean.wave_functions.plot_glukhovskiy_pdf`: Plots the probability density function of a Glukhovskiy distribution.
    - {func}`layowt.metocean.wave_functions.plot_glukhovskiy_poe`: Plots the probability of exceedance of a Glukhovskiy distribution.
    - {func}`layowt.metocean.wave_functions.goda_breaking_height`: Estimates the breaking wave height for a regular wave in a sloped sea bed, Goda 2010.
    - {func}`layowt.metocean.wave_functions.forristal_hp`: Computes the wave height exceeded by P% of the waves, :math:`H_{P}` of a Forristal wave height distribution.
    - {func}`layowt.metocean.wave_functions.plot_forristal_poe`: Plots the probability of exceedance of a Forristal wave height distribution.
- {mod}`layowt.metocean.battjes`: Battjes and Groenendijk distrtibution module:
    - {class}`layowt.metocean.battjes.BattjesGroenendijk`: Battjes and Groenendijk distrtibution.
        - {meth}`layowt.metocean.battjes.BattjesGroenendijk.gammainc_upper`: Upper incomplete gamma function.
        - {meth}`layowt.metocean.battjes.BattjesGroenendijk.gammainc_lower`: Lower incomplete gamma function.
        - {meth}`layowt.metocean.battjes.BattjesGroenendijk.get_Hp`: Computes the wave height exceeded by P% of the waves.
        - {meth}`layowt.metocean.battjes.BattjesGroenendijk.get_Hn`: Computes the mean of the highest 1/N-part of the wave heights.
        - {meth}`layowt.metocean.battjes.BattjesGroenendijk.plot_poe`: Plots the probability of exceedance of the distribution.
        - {meth}`layowt.metocean.battjes.BattjesGroenendijk.plot_pdf`: Plots the probability density function of the distribution.
- {mod}`layowt.metocean.goda`: Goda model for maximum wave height module:
    - {class}`layowt.metocean.goda.Goda`: Goda model for maximum wave height.



 
### Changed
- {func}`layowt.metocean.wave_functions.dnv_breaking_height`: Renamed from `breaking_height` to differentiate from the new Goda 2010 implementation of the breaking wave height on sloped seabeds.



### Fixed
- {meth}`layowt.layouts.layout_generators.GriddedLayoutGenerator.generate_layouts`: Fixed error when trying to export layouts to file when EPSG is not specified.



=======================================================================
## [0.6.0] - [2023-07-21]
 
LayOWt Version 0.6.0

New Metocean Module!
 
### Added
- {mod}`layowt.metocean`: Metocean module:
    - {class}`layowt.metocean.wave.Wave` : Wave class using Fenton's stream function wave theory to estimate non-linear wave parameters using truncated Fourier series.
        - {meth}`layowt.metocean.wave.Wave.validate_theory` shows what wave theory is applicable for the Wave instance using the wave theories' figure by Le Mehaute.
        - {meth}`layowt.metocean.wave.Wave.ursell` calculates the Ursell number for the Wave instance.
        - {meth}`layowt.metocean.wave.Wave.williams_breaking_criterion` for plotting the Williams (1981) analytical breaking wave criterion for the Wave instance.
        - {meth}`layowt.metocean.wave.Wave.plot` for plotting the wave profile with multiple depth profiles for a given parameter.
    - {func}`layowt.metocean.wave_functions.deep_wavelength` for calculating the deep water wavelength from the wave period. Applicable for linear wave theory.
    - {func}`layowt.metocean.wave_functions.intermediate_wavelength_approx` for calculating the approximate intermediate water depth wavelength from the wave period and water depth. Applicable for linear wave theory.
    - {func}`layowt.metocean.wave_functions.intermediate_wavelength` for calculating the intermediate water depth wavelength from the wave period and water depth iteratively. Applicable for linear wave theory.
    - {func}`layowt.metocean.wave_functions.shallow_wavelength` for calculating the shallow water depth wavelength from the wave period and water depth. Applicable for linear wave theory.
    - {func}`layowt.metocean.wave_functions.breaking_height` for calculating the breaking wave height for a regular wave in a plane sea bed.
    - {func}`layowt.metocean.wave_functions.iribarren` for calculating the Iribarren number (also known as surf similarity parameter).
    - {func}`layowt.metocean.wave_functions.breaker_type` for determining the type of breaking wave based on the Iribarren number.
    - {func}`layowt.metocean.wave_functions.second_order_crest` for calculating the wave crest height of a Second-order Stokes wave in deep water.
    - {func}`layowt.metocean.wave_functions.second_order_trough` for calculating the wave trough height of a Second-order Stokes wave in deep water.
    - {func}`layowt.metocean.wave_functions.pierson_moskowitz` for calculating the spectral density of the Pierson-Moskowitz wave spectrum for a sea state at a given frequency.
    - {func}`layowt.metocean.wave_functions.plot_pierson_moskowitz` for plotting the Pierson-Moskowitz wave spectrum for a given sea state.
    - {func}`layowt.metocean.wave_functions.jonswap` for calculating the spectral density of the JONSWAP wave spectrum for a sea state at a given frequency.
    - {func}`layowt.metocean.wave_functions.plot_jonswap` for plotting the JONSWAP wave spectrum for a given sea state.
    - {func}`layowt.metocean.wave_functions.ursell` for calculating the Ursell number which indicates the nonlinearity of waves.
    - {func}`layowt.metocean.wave_functions.validate_theory` shows what wave theory is applicable for given wave parameters using the wave theories' figure by Le Mehaute.
    - {func}`layowt.metocean.wave_functions.gamma_from_mean_periods` for estimating the JONSWAP peak enhancement factor from the spectral mean periods.
    - {func}`layowt.metocean.wave_functions.tp_from_mean_period` for estimating the spectral peak wave period from the JONSWAP peak enhancement factor and a mean spectral period.
    - {func}`layowt.metocean.wave_functions.williams_breaking_criterion` for plotting the Williams (1981) analytical breaking wave criterion for the given wave parameters.
    - {func}`layowt.metocean.wave_functions.dingemans_gamma` for estimating the JONSWAP peak enhancement factor from the spectral mean periods using Dingemans' tabulated experimental results.
    - {func}`layowt.metocean.wave_functions.dingemans_tp` for estimating the spectral peak wave period from the spectral mean periods using Dingemans' tabulated experimental results.
 

 
### Changed
- {meth}`layowt.layouts.layout.Layout.plot` method now accepts the `basemap` and `basemap_provider` arguments to plot basemaps.




### Fixed
- Fixed data types of grid parameters when loading a layout from the database. Previously many parameters loaded as SQLAlchemy Decimal type. Now all are base python floats.
- FIxed erroneous read of rated wind speed when loading a {class}`layowt.wtg.power_curve.PowerCurve` using the {meth}`layowt.wtg.power_curve.PowerCurve.from_wtg()` constructor method.
- Fixed bug in AEP functions where the wind speed binning was not correct.




=======================================================================
## [0.5.0] - [2023-01-23]
 
LayOWt Version 0.5.0

Added the capability to upload layouts directly to the internal OW PostGIS Spatial Database. 
 
### Added
- {meth}`layowt.layouts.layout.Layout.database_engine` method to create an SQLAlchemy database engine.
- {func}`layowt.layouts.utils.database_engine` method to create an SQLAlchemy database engine.
- {func}`layowt.layouts.utils.setup_database_schema` method to automatically set up a database schema with all the necessary sytem tables in order to upload layouts.
- {meth}`layowt.layouts.layout.Layout.to_database()` method to upload a Layout to the OW postGIS database.
- {meth}`layowt.layouts.floating_layout.FloatingLayout.to_database()` method to upload a FloatingLayout to the OW postGIS database.
- {func}`layowt.layouts.utils.coordinate_transform` function to transform coordinates between coordinate reference systems using EPSG codes.
- {func}`layowt.layouts.utils.raster_from_postgis` function to load a raster dataset into memory from a PostGIS database.
- {func}`layowt.layouts.utils.geoms_from_database` function to load area or exclusion geometries from the OW Energy assessment PostGIS database table system.
- {func}`layowt.aep.aep.ts_aep` function to calculate the AEP of a PowerCurve using a timeseries of wind speeds.
- {func}`layowt.aep.aep.ts_power` function to calculate a power timeseries for a PowerCurve using a timeseries of horiztonal wind speeds.


 
### Changed
- {meth}`layowt.layouts.layout.Layout.from_postgis` no longer requires username, password, database name, and host as inputs. Takes in an SQLAlchemy engine directly.
- {func}`layowt.layouts.utils.geoms_from_postgis` no longer requires username, password, database name, and host as inputs. Takes in an SQLAlchemy engine directly.
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` now accepts a rasterio MemoryFile type for the `bathymetry_path` argument.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` now accepts a rasterio MemoryFile type for the `bathymetry_path` argument.



### Fixed




=======================================================================
## [0.4.0] - [2023-01-23]
 
LayOWt Version 0.4.0

Added the errors module for LayOWt specific Exceptions, statistics module, and minor additions to existing modules.
Can now import the package directly and access all subpackages with a single import.
 
### Added
- {mod}`layowt.errors`: Errors module with LayOWt specific Exception classes.
    - {class}`layowt.errors.errors.EmptyLayoutError` for raising errors when the layout geometry is empty.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` class:
    - Can now filter generated layouts but mooring line crossings and mooring buffer overlaps using the allow_mooring_crossings and allow_mooring_buffer_overlap arguments.
- {mod}`layowt.aep`: Module:
    - {func}`layowt.aep.aep.weibull_ws_to_aep` function to calculate the approximate linear relationship between wind speed changes and AEP changes for a 2-parameter Weibull distribution.
     {func}`layowt.aep.aep.ts_ws_to_aep` function to calculate the approximate linear relationship between wind speed changes and AEP changes given a wind speed timeseries.
    - {func}`layowt.aep.aep.weibull_mean_map` function to crate maps of mean 2-Parameter Weibull distribution mean values given shape and scale input raster datasets.
- {mod}`layowt.stats`: Statistics module:
    - {func}`layowt.stats.stats.weibull_mean` Function to calculate the mean value of a 2-parameter Weibull distribution.
    - {func}`layowt.stats.circ_stats.circ_mean` Function to compute the mean value of direction data.
    - {func}`layowt.stats.circ_stats.circ_corr` Function to compute the circular correlation coefficient between two arrays of direction data.
- {class}`layowt.layouts.layout.Layout` class:
    - {attr}`layowt.layouts.layout.Layout.shortest_distance` property added that returns the shortest euclidean distance between turbines in the layout projection units.
    - {meth}`layowt.layouts.layout.Layout.from_xy` constructor method to create layouts by passing arrays of x and y coordinates directly.
- {class}`layowt.grids.grid.Grid` class:
    - {attr}`layowt.grids.grid.Grid.shortest_distance` property added that returns the shortest euclidean distance between grid points in the grid projection units.
    - {meth}`layowt.grids.grid.Grid.to_shapefile` method to export a Grid object to a shapefile.
    - {meth}`layowt.grids.grid.Grid.to_multiline` method to return a Shapely MultiLineString representation of the Grid object.


 
### Changed
- Can now import the `layowt` package and access all of its subpackage namespaces directly.
- {func}`layowt.layouts.utils.hdf5_to_dataframe` now exports the row_offset and col_offset columns of the info dataframe as boolean.



### Fixed
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` class:
    - {meth}`layowt.layouts.layout_generators.GriddedLayoutGenerator.generate_layouts` method no longer crashes when checking the number of turbines of an empty layout.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` class:
    - {meth}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator.generate_layouts` method no longer crashes when checking the number of turbines of an empty layout.
- {class}`layowt.wtg.power_curve.PowerCurve` class:
    - {meth}`layowt.wtg.power_curve.PowerCurve.to_wtg()` method bugfix. No longer crashes when exporting to .wtg format with missing attributes.
- {class}`layowt.layouts.layout.Layout` class:
    - {attr}`layowt.layouts.layout.Layout.n_wtg` no longer raises an AttributeError when the Layout geometry is a single Point geometry.
    - {meth}`layowt.layouts.layout.Layout.from_postgis` method now correctly parses special characters such as @ and / in the user's password.
- {func}`layowt.layouts.utils.geoms_from_postgis` function now correctly parses special characters such as @ and / in the user's password.



=======================================================================
## [0.3.0] - [2022-11-07]
 
LayOWt Version 0.3.0

Added angular buffers to the {class}`layowt.layouts.floating_layout.FloatingLayout` class.
 
### Added
- {class}`layowt.layouts.floating_layout.FloatingLayout`:
    - {meth}`layowt.layouts.floating_layout.FloatingLayout.generate_mooring_buffers` method to generate angular buffers around mooring lines.
    - {attr}`layowt.layouts.floating_layout.FloatingLayout.has_mooring_buffers` property to identify if mooring buffers have been generated.
    - {attr}`layowt.layouts.floating_layout.FloatingLayout.mooring_buffer_overlap` property to indentify if mooring buffers overlap.


 
### Changed
- {class}`layowt.layouts.floating_layout.FloatingLayout`:
    - {meth}`layowt.layouts.floating_layout.FloatingLayout.plot` method now also plots mooring buffers if they have been generated.
    - {meth}`layowt.layouts.floating_layout.FloatingLayout.to_shapefile` method now also exports mooring buffers to a shapefile if they have been generated.
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` class:
    - Layout generators can now take [True, False] as an input for the row_offset and col_offset arguments to generated layouts with and without offests in a single pass.

### Fixed



=======================================================================
## [0.2.2] - [2022-10-28]
 
LayOWt Version 0.2.2

Bugfixes to {func}`layowt.layouts.utils.layouts_to_hdf5` function.
 
### Added
- {func}`layowt.layouts.utils.hdf5_to_dataframe` function has additional out and indexes arguments. The out argument lets the user
select if either or both of the layout info and coordinates dataframes are returned by the function. The indexes argument lets the user
select layouts by index rather than by continuous chunks as when using the start_id and end_id arguments.
- {func}`layowt.layouts.utils.layouts_from_hdf5` function has additional indexes argument. This new argument lets the user
select layouts by index rather than by continuous chunks as when using the start_id and end_id arguments.

 
### Changed


### Fixed
- Fixed bug where the {func}`layowt.layouts.utils.layouts_to_hdf5` function was not appending layout data to the hdf5 file from the second dump onwards.


=======================================================================
## [0.2.1] - [2022-10-26]
 
LayOWt Version 0.2.1

Fixing the layout generators memory leak problem and minor quality of life updates.
 
### Added
- {func}`layowt.layouts.utils.layouts_to_hdf5` function. This function allows to save a list of layouts to an hdf5 layouts file.
- {func}`layowt.layouts.utils.layouts_from_hdf5` function. This function allows to load a list of layouts from an hdf5 layouts file.
- {func}`layowt.layouts.utils.hdf5_to_dataframe` function. This function allows to load layouts DataFrames from an hdf5 layouts file.
 
### Changed
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` class:
    - The n_wtg argument has been changed from an int to a tuple(int, int). This now allows to select layouts between a minimum and maximum number of turbines.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` class:
    - The n_wtg argument has been changed from an int to a tuple(int, int). This now allows to select layouts between a minimum and maximum number of turbines.

### Fixed
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` class:
    - Memory leak bug has been addressed by adding the option to save results directly to file and clearing memory when using the {meth}`layowt.layouts.layout_generators.GriddedLayoutGenerator.generate_layouts()` method.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` class:
    - Memory leak bug has been addressed by adding the option to save results directly to file and clearing memory when using the {meth}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator.generate_layouts()` method.


=======================================================================
## [0.2.0] - [2022-10-17]
 
LayOWt Version 0.2.0

Added the FloatingLayout class and Moorings!
 
### Added
- {func}`layowt.aep.aep.weibull_wpd` function. This function calculates the wind power density of a weibull distribution.
- {func}`layowt.aep.aep.weibull_aep_map` function. This function allows for the creation of GAEP raster datasets.
- {class}`layowt.layouts.floating_layout.FloatingLayout` class: New Layout subclass with mooring design functionality!
- {class}`layowt.layouts.layout.Layout` class:
    - Now stores the attribute {attr}`layowt.layouts.layout.Layout.bathymetry_band`.
    - {meth}`layowt.layouts.layout.Layout.to_shapefile` method epsg argument now optional. If not passed, takes the layout epsg code if available.
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` class:
    - Now has an implemented \_\_len\_\_ method. Returns the number of layout combinations that will be created when running the {meth}`layowt.layouts.layout_generators.GriddedLayoutGenerator.generate_layouts()` method.
    - {meth}`layowt.layouts.layout_generators.GriddedLayoutGenerator.generate_layouts()` method has an optional `progress` argument to show a progress bar.
- {func}`layowt.layouts.utils.raster_slope` function. This function calculates a slope raster dataset from an elevation dataset.
- {func}`layowt.layouts.utils.legacy_csv_to_layouts` function. This function reads the legacy layout .csv files and rebuilds Layouts from the contents.
- Tutorials:
    - Added {class}`layowt.layouts.floating_layout.FloatingLayout` class tutorial.
 
### Changed
 - {func}`layowt.layouts.utils.layouts_to_legacy_csv` function now writes the number of turbines, epsg, and row and col numbers in the INFO_ .csv file.

 
### Fixed
- {class}`layowt.layouts.layout.Layout` class:
    - Now the {meth}`layowt.layouts.layout.Layout.clip_to_area()` and {meth}`layowt.layouts.layout.Layout.avoid_exclusions()` methods
    recalculate the bathymetry samples if they are applied to a layout with an existing bathymetry constraint. Previously the geometries
    would update correctly, but the list of bathymetry samples corresponded to the values for the layout prior to the new area or
    exclusion constraint being applied.
    - The attribute {attr}`layowt.layouts.layout.Layout.bathymetry_band` now stored as an attribute.  
- {class}`layowt.wtg.power_curve.PowerCurve` class:
    - Fixed errors when reading .wtg files without release dates or revision numbers using {meth}`layowt.wtg.power_curve.PowerCurve.from_wtg()`.
    - Fixed errors when writing .wtg files for PowerCurves without document references, revision, or release dates using {meth}`layowt.wtg.power_curve.PowerCurve.to_wtg()`.



=======================================================================
## [0.1.0] - [2022-09-05]
 
First official release of LayOWt!
 
### Added
- {class}`layowt.layouts.layout.Layout` class:
    - {meth}`layowt.layouts.layout.Layout.load_bathymetry()` and {meth}`layowt.layouts.layout.Layout.apply_bathymetry()` methods now store the bathymetry values in the {attr}`layowt.layouts.layout.Layout.samples` class attribute.
    - {meth}`layowt.layouts.layout.Layout.load_bathymetry()` and {meth}`layowt.layouts.layout.Layout.apply_bathymetry()` methods now reproject layout coordinates in the background if the raster CRS is different from the layout CRS.
    - Layout CRS EPSG code now stored in the {attr}`layowt.layouts.layout.Layout.epsg` attribute.
    - Layout bathymetry CRS EPSG code now stored in the {attr}`layowt.layouts.layout.Layout.bathymetry_epsg` attribute.
- New {mod}`layowt.wtg` module:
    - {class}`layowt.wtg.power_curve.PowerCurve` class with the capability to create power curves from the OW database directly.
        - {meth}`layowt.wtg.power_curve.PowerCurve.set_hub_height()` method to set the hub height on a {class}`layowt.wtg.power_curve.PowerCurve` instance.
        - {meth}`layowt.wtg.power_curve.PowerCurve.database_engine()` static method to create a {class}`sqlalchemy.engine.Engine` database connection to the corporate OW power curve database.
        - {meth}`layowt.wtg.power_curve.PowerCurve.from_database()` class method to create {class}`layowt.wtg.power_curve.PowerCurve` instances by `curve_id` from the corporate OW power curve database.
        - {meth}`layowt.wtg.power_curve.PowerCurve.plot()` method to plot the power and thurst curves.
        - {meth}`layowt.wtg.power_curve.PowerCurve.from_wtg()` class method to load WAsP .wtg files as an {class}`layowt.wtg.power_curve.PowerCurve` instance.
        - {meth}`layowt.wtg.power_curve.PowerCurve.to_wtg()` method to export a {class}`layowt.wtg.power_curve.PowerCurve` instance into a WAsP .wtg file.
        - {meth}`layowt.wtg.power_curve.PowerCurve.compare_curves_from_database()` class method to load multiple {class}`layowt.wtg.power_curve.PowerCurve` instances from the OW corporate power curve database, and compare them by potting their power and thurst curves.
        - {meth}`layowt.wtg.power_curve.PowerCurve.get_commercial_turbines()` class method that return a dataframe of all the commercial turbine attributes stored in the OW corporate power curve database.
        - {meth}`layowt.wtg.power_curve.PowerCurve.interpolate()` method to interpolate the {class}`layowt.wtg.power_curve.PowerCurve` instance with a new wind speed step.
        - {meth}`layowt.wtg.power_curve.PowerCurve.from_text()` method to load a {class}`layowt.wtg.power_curve.PowerCurve` instance using the Pandas {func}`pandas.read_csv` function.
- Documentation now has links to source code.

### Changed
- {class}`layowt.layouts.layout.Layout` class:
    - {meth}`layowt.layouts.layout.Layout.plot()` method now uses {class}`geopandas.GeoDataFrame` plotting functionality for simpler logic on multipolygons.
- {func}`layowt.aep.aep.weibull_aep()` function now takes in a LayOWt {class}`layowt.wtg.power_curve.PowerCurve` instead of a PyWake WindTurbine.
- {func}`layowt.layouts.utils.layouts_to_legacy_csv()` function now exports all grid parameters of the layout to the INFO file.

 
### Fixed
- {class}`layowt.layouts.layout.Layout` class:
    - {meth}`layowt.layouts.layout.Layout.plot()` method now does not raise an {exc}`AttributeError` if there are no turbines to plot. Plots an empty figure with no turbines.
    - All methods that take epsg codes now check for valid codes against the PROJ EPSG CRS database.
- Tutorials:
    - tutorial_inputs folder added to the repository with all example data used throughout the tutorials. Tutorial notebooks can now be executed locally with all the example data.
    - Updated the `Layout` class tutorial with latest CRS cahnges to the `Layout` class.



=======================================================================
## [0.1.0-beta] - [2022-07-06]
 
First pre-release of LayOWt!

Please test away and send your feedback!
 
### Added
- {mod}`layowt.grids` module
- {mod}`layowt.layouts` module
- {mod}`layowt.aep` module
 
### Changed
 

 
### Fixed




=======================================================================

