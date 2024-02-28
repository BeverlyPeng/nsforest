# Configuration file for the Sphinx documentation builder.

import os
import sys

sys.path.insert(0, os.path.abspath("../../nsforest"))

# -- Project information

project = 'NSForest'
authors = [
    "Renee Zhang <zhangy@jcvi.org>",
    "Richard Scheuermann <RScheuermann@jcvi.org>",
    "Brian Aevermann <baevermann@chanzuckerberg.com>", 
    "Angela Liu <aliu@jcvi.org>", 
    "Beverly Peng <bpeng@jcvi.org">
]

version = '4.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    'anndata': ("https://anndata.readthedocs.io/en/stable/", None),
    'pandas': ("https://pandas.pydata.org/pandas-docs/stable/", None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'
