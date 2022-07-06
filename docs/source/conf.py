# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))

from datetime import date

from layowt import __release__, __version__

# -- Project information -----------------------------------------------------

project = 'LayOWt'
copyright = str(date.today().year) + ', Guillermo Tornero'
author = 'Guillermo Tornero'

# The short X.Y version
version = __version__

# The full version, including alpha/beta/rc tags
release = __release__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "myst_parser",
    "sphinx.ext.duration",
    "sphinx.ext.autosectionlabel",
    #"sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    #"sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "autoapi.extension",
    "nbsphinx", # can also use MyST-NB
    "hoverxref.extension",
    "sphinx_copybutton",
    "sphinx.ext.mathjax",
]

myst_enable_extensions = ["dollarmath", "amsmath"]

autosummary_generate = True

autoapi_dirs = ["../../layowt"]

autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
    'special-members',
    'imported-members',
    ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_static"]

# Intersphinx configuration
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy-1.8.1/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "py_wake": ("https://topfarm.pages.windenergy.dtu.dk/PyWake/", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable/", None),
    "rasterio": ("https://rasterio.readthedocs.io/en/latest/", None),
}

# Hoverxref Extension
hoverxref_auto_ref = True
hoverxref_mathjax = True
hoverxref_intersphinx = [
    "python",
    "numpy",
    "scipy",
    "matplotlib",
    "py_wake",
    "shapely",
    "rasterio",
]
hoverxref_domains = ["py"]
hoverxref_role_types = {
    "hoverxref": "modal",
    "ref": "modal",  # for hoverxref_auto_ref config
    "confval": "tooltip",  # for custom object
    "mod": "tooltip",  # for Python Sphinx Domain
    "class": "tooltip",  # for Python Sphinx Domain
    "meth": "tooltip",
    "obj": "tooltip",
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'

html_logo = '../../img/logo_small.png'

html_favicon = '../../img/favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
