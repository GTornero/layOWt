# LayOWt Tutorial

Welcome to LayOWt's tutorial page!

This page assumes that you are already familiar with the LayOWt installation process and how to use Jupyter Notebooks.

```{toctree}
:caption: 'Contents:'
:maxdepth: 2

Tutorial 1: Generating a Grid <notebooks/grid_tutorial>
Tutorial 2: The Layout Object <notebooks/layout_tutorial>
Tutorial 3: Procedural Generation of Gridded Layouts <notebooks/GriddedLayoutGenerator_tutorial>
Tutorial 4: The PowerCurve Object <notebooks/PowerCurve_tutorial>
```

## [Tutorial 1: Generating a Grid](notebooks/grid_tutorial)

This tutorial shows you how to generate unifrom grids of points using the custom LayOWt `Grid` class.

## [Tutorial 2: The Layout Object](notebooks/layout_tutorial)

This tutorial introduces you to the LayOWt `Layout` object. This is the building block to create project scenarios.
This section will cover the many ways to create a `Layout` and how to carry out geospatial analysis with it.

## [Tutorial 3: Procedural Generation of Gridded Layouts](notebooks/GriddedLayoutGenerator_tutorial)

This tutorial introduces you to LayOWt's layout generators, and how to use them to quickly produce many `Layout` configurations in only a few lines of code.
This will cover the specific `GriddedLayoutGenerator` objects that uses `Grid` objects to define layouts.

## [Tutorial 4: The PowerCurve Object](notebooks/PowerCurve_tutorial)

This tutorial introduces LayOWt's `PowerCurve` object. This object brings some handy functionality to easily interact with OW's internal power curve database
for both commercial and internally extrapolated turbines and their power curves. The `PowerCurve` object can be created from text data, from the database,
or from the commonly used DTU WAsP .wtg format. Can also be easily exported into these same formats.

Common day-to-day tasks can also be easliy carried out using its methods, taks such as power curve interpolation and AEP calculations.

