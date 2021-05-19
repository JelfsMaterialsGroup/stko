# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a
# full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup -------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another
# directory,
# add these directories to sys.path here. If the directory is relative
# to the
# documentation root, use os.path.abspath to make it absolute, like
# shown here.

import sys
import os
# Appends main directory to the path to import the package
sys.path.insert(0, os.path.abspath('../../'))

# For mocking external dependencies when building.
from unittest.mock import MagicMock


class Mock(MagicMock):
    @classmethod
    def __getattr__(cls, name):
        return MagicMock()


MOCK_MODULES = [
    'rdkit',
    'rdkit.Chem',
    'rdkit.Chem.AllChem',
    'rdkit.Geometry',
    'rdkit.Geometry.Point3D'
    'numpy',
    'numpy.linalg',
    'scipy',
    'scipy.spatial',
    'scipy.spatial.transform',
    'scipy.spatial.distance',
    'scipy.constants',
    'scipy.optimize',
    'matplotlib',
    'matplotlib.pyplot',
    'pandas',
    'pathos',
    'seaborn',
    'stk',
]
sys.modules.update((mod_name, Mock()) for mod_name in MOCK_MODULES)


# -- Project information ----------------------------------------------

project = 'stko'
copyright = '2020, Steven Bennett, Andrew Tarzia'
author = 'Steven Bennett, Andrew Tarzia'

# The full version, including alpha/beta/rc tags
version = '0.0.1'


# -- General configuration --------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx_rtd_theme',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
]

autodoc_default_options = {
    'special-members': '__init__',
    'inherited-members': True,
    'show-inheritance': True,
    'ignore-module-all': True,
}


# add_module_names = False


# Add any paths that contain templates here, relative to this
# directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce
# nothing.
todo_include_todos = True


# -- Options for HTML output ------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation
# for a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"

html_theme_options = {
    'collapse_navigation': False,
}

# Theme options are theme-specific and customize the look and feel of
# a theme further.  For a list of options available for each theme,
# see the documentation.
# html_theme_options = {
#     'collapse_navigation': False,
# }
