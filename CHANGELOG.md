% (CHANGELOG template from https://gist.github.com/juampynr/4c18214a8eb554084e21d6e288a18a2c)
% (How to use Git to create a changelog https://www.freecodecamp.org/news/a-beginners-guide-to-git-what-is-a-changelog-and-how-to-generate-it/)

# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

=======================================================================
## [0.1.0] - [ADD DATE HERE!!]
 
First official release of LayOWt!
 
### Added
- {meth}`layowt.layouts.layout.Layout.load_bathymetry()` and {meth}`layowt.layouts.layout.Layout.apply_bathymetry()` methods now store the bathymetry values in the {attr}`layowt.layouts.layout.Layout.samples` class attribute.
- {mod}`layowt.wtg` module
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

### Changed
 - {meth}`layowt.layouts.layout.Layout.plot()` method now uses {class}`geopandas.GeoDataFrame` plotting functionality for simpler logic on multipolygons.
 - {func}`layowt.aep.aep.weibull_aep()` function now takes in a LayOWt {class}`layowt.wtg.power_curve.PowerCurve` instead of a PyWake WindTurbine.

 
### Fixed
- {meth}`layowt.layouts.layout.Layout.plot()` method now does not raise an {exc}`AttributeError` if there are no turbines to plot. Plots an empty figure with no turbines.



=======================================================================
## [0.1.0-beta] - 2022-07-06
 
First pre-release of LayOWt!

Please test away and send your feedback!
 
### Added
- {mod}`layowt.grids` module
- {mod}`layowt.layouts` module
- {mod}`layowt.aep` module
 
### Changed
 

 
### Fixed




=======================================================================

