# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# -- Project information ----------------------------------------------

project = "stko"
project_copyright = "2023, Steven Bennett, Andrew Tarzia, Lukas Turcani"
author = "Steven Bennett, Andrew Tarzia, Lukas Turcani"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.doctest",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx_copybutton",
    "moldoc",
]

autosummary_imported_members = True

autodoc_typehints = "description"
autodoc_member_order = "groupwise"
autoclass_content = "class"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "openff": (
        "https://docs.openforcefield.org/projects/toolkit/en/stable/",
        None,
    ),
    "openmm": ("http://docs.openmm.org/latest/api-python/", None),
}

templates_path = ["_templates"]
exclude_patterns: list[str] = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
