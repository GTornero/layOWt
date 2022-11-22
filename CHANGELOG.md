% (CHANGELOG template from https://gist.github.com/juampynr/4c18214a8eb554084e21d6e288a18a2c)
% (How to use Git to create a changelog https://www.freecodecamp.org/news/a-beginners-guide-to-git-what-is-a-changelog-and-how-to-generate-it/)

# Change Log
All notable changes to this project will be documented in this file.
 
The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

=======================================================================
## [0.4.0] - [2022-11-22]
 
LayOWt Version 0.4.0

Added the errors module for LayOWt specific Exceptions and minor additions to existing modules.
 
### Added
- {mod}`layowt.errors`: Errors module with LayOWt specific Exception classes.
    - {class}`layowt.errors.errors.EmptyLayoutError` for raising errors when the layout geometry is empty.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` class:
    - Can now filter generated layouts but mooring line crossings and mooring buffer overlaps using the allow_mooring_crossings and allow_mooring_buffer_overlap arguments.
- {func}`layowt.aep.aep.weibull_ws_to_aep` function to calculate the approximate linear relationship between wind speed changes and AEP changes.
- {class}`layowt.layouts.layout.Layout` class:
    - {attr}`layowt.layouts.layout.Layout.shortest_distance` property added that returns the shortest euclidean distance between turbines in the layout projection units.
- {class}`layowt.grids.grid.Grid` class:
    - {attr}`layowt.grids.grid.Grid.shortest_distance` property added that returns the shortest euclidean distance between grid points in the grid projection units.


 
### Changed
- {func}`layowt.layouts.utils.hdf5_to_dataframe` now exports the row_offset and col_offset columns of the info dataframe as boolean.



### Fixed
- {class}`layowt.layouts.layout_generators.GriddedLayoutGenerator` class:
    - {meth}`layowt.layouts.layout_generators.GriddedLayoutGenerator.generate_layouts` method no longer crashes when checking the number of turbines of an empty layout.
- {class}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator` class:
    - {meth}`layowt.layouts.layout_generators.FloatingGriddedLayoutGenerator.generate_layouts` method no longer crashes when checking the number of turbines of an empty layout.



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

