# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys
sys.path.insert(0, os.path.abspath('..')) 

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



# -- Options for autodoc extension -------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#module-sphinx.ext.autodoc
autodoc_mock_imports = ["numpy","scipy", "pybedtools", "pandas", "pyBigWig", "deeptools", "ortools", "typing", "pysam", "myst_parser"]
autodoc_typehints = 'both'
autodoc_typehints_format = 'fully-qualified'
autodoc_default_options = {
    'members': True,
}

# -- Options for napoleon extension ------------------------------------------
napoleon_sphinx_docstring = True
html_theme_options = {
    'github_user': 'nolan-h-hamilton',
    'github_repo': 'ROCCO',
    'page_width': '1200px',
    'sidebar_width': '300px',
}
html_sidebars = { '*': ['searchbox.html', 'localtoc.html'] }