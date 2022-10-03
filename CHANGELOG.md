% (CHANGELOG template from https://gist.github.com/juampynr/4c18214a8eb554084e21d6e288a18a2c)
% (How to use Git to create a changelog https://www.freecodecamp.org/news/a-beginners-guide-to-git-what-is-a-changelog-and-how-to-generate-it/)

# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

=======================================================================
## [0.2.0] - [INSERT DATE HERE!]
 
LayOWt Version 0.2.0

Added the FloatingLayout class and Moorings!
 
### Added
- {func}`layowt.aep.aep.weibull_wpd` function.
- {class}`layowt.layouts.floating_layout.FloatingLayout` class: New Layout subclass with mooring design functionality!
- {class}`layowt.layouts.layout.Layout` class:
    - Now stores the attribute {attr}`layowt.layouts.layout.Layout.bathymetry_band`.
 
### Changed
 

 
### Fixed
- {class}`layowt.layouts.layout.Layout` class:
    - Now the {meth}`layowt.layouts.layout.Layout.clip_to_area()` and {meth}`layowt.layouts.layout.Layout.avoid_exclusions()` methods
    recalculate the bathymetry samples if they are applied to a layout with an existing bathymetry constraint. Previously the geometries
    would update correctly, but the list of bathymetry samples corresponded to the values for the layout prior to the new area or
    exclusion constraint being applied.
    - The attribute {attr}`layowt.layouts.layout.Layout.bathymetry_band` now stored as an attribute.
- {class}`layowt.wtg.power_curve.PowerCurve class:
    - Fixed errors when reading .wtg files without release dates or revision numbers using {meth}`layowt.wtg.power_curve.PowerCurve.from_wtg()`.
    - Fixed errors when writing .wtg files for PowerCurves without document references, revision, or release dates using meth}`layowt.wtg.power_curve.PowerCurve.to_wtg()`.



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

