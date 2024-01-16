# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('../rocco'))
project = 'ROCCO'
author = 'Nolan H. Hamilton'


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.doctest',
    'myst_parser',
]


templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

latex_elements = {'classoptions': ',openany,oneside',     'babel' : '\\usepackage[english]{babel}',     'preamble': r'''
        \usepackage{charter}
        \usepackage[defaultsans]{lato}
        \usepackage{inconsolata}
    ''',}


html_theme = 'alabaster'
html_static_path = ['_static']

# -- Options for autodoc extension -------------------------------------------
autodoc_member_order = 'alphabetical'
autodoc_typehints = 'both'
autodoc_typehints_format = 'fully-qualified'
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}
add_module_names = True
# -- Options for napoleon extension ------------------------------------------
napoleon_sphinx_docstring = True
html_sidebars = { '**': ['searchbox.html'] }