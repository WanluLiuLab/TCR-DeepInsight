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
import os
import sys
sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath(os.path.join('..', '..')))
# sys.path.insert(0, os.path.abspath(os.path.join('..', '..', 't_deep_insight')))
print(sys.version)
# -- Project information -----------------------------------------------------

project = 'TCR-DeepInsight'
copyright = '2024, Ziwei Xue'
author = 'Ziwei Xue'

repository_url = "https://github.com/WanluLiuLab/TCR-DeepInsight"

# The full version, including alpha/beta/rc tags
release = 'latest'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx_rtd_theme',
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosectionlabel',
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.githubpages",
    # "sphinxext.opengraph",
    # "sphinxcontrib.jquery",
    # "scanpydoc",
    'hoverxref.extension',
    # "scanpydoc.elegant_typehints",
    # "scanpydoc.autosummary_generate_imported",
    "sphinx_design",
    "sphinx_copybutton",
    "nbsphinx",
    "sphinxcontrib.bibtex"
]

bibtex_bibfiles = ["references.bib"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'furo'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_css_files = [
    'css/custom.css',
]
html_js_files = [
    'js/custom.js',
]

html_show_sphinx = False
html_logo = "_static/logo.svg"
html_theme_options = {
    "repository_url": repository_url,
    'logo_only': False,
    'display_version': False,
    'dark_css_variables': {
        # Taken from: https://github.com/pradyunsg/furo/blob/c682d5d3502f3fa713c909eebbf9f3afa0f469d9/src/furo/assets/styles/variables/_colors.scss
        'color-problematic': '#b30000',

        # Base Colors
        'color-foreground-primary': 'black', # for main text and headings
        'color-foreground-secondary': '#5a5c63', # for secondary text
        'color-foreground-muted': '#646776', # for muted text
        'color-foreground-border': '#878787', # for content borders

        'color-background-primary': 'white', # for content
        'color-background-secondary': '#f8f9fb', # for navigation + ToC
        'color-background-hover': '#efeff4ff', # for navigation-item hover
        'color-background-hover--transparent': '#efeff400',
        'color-background-border': '#eeebee', # for UI borders
        'color-background-item': '#ccc', # for "background" items (eg: copybutton)

        # Announcements
        'color-announcement-background': '#000000dd',
        'color-announcement-text': '#eeebee',

        # Brand colors
        'color-brand-primary': '#2962ff',
        'color-brand-content': '#2a5adf',

        # Highlighted text (search)
        'color-highlighted-background': '#ddeeff',

        # GUI Labels
        'color-guilabel-background': '#ddeeff80',
        'color-guilabel-border': '#bedaf580',

        # API documentation
        'color-api-keyword': 'var(--color-foreground-secondary)',
        'color-highlight-on-target': '#ffffcc',

        # Admonitions
        'color-admonition-background': 'transparent',

        # Cards
        'color-card-border': 'var(--color-background-secondary)',
        'color-card-background': 'transparent',
        'color-card-marginals-background': 'var(--color-background-hover)',

        # Code blocks
        'color-code-foreground': 'black',
        'color-code-background': '#f8f9fb',
    }
}
# Force pygments style in dark mode back to the light variant
pygments_dark_style = 'tango'


# -- Config for hoverxref -------------------------------------------

hoverx_default_type = "tooltip"
hoverxref_domains = ["py"]
hoverxref_role_types = dict.fromkeys(
    ["ref", "class", "func", "meth", "attr", "exc", "data", "mod"],
    "tooltip",
)
hoverxref_intersphinx = [
    "python",
    "numpy",
    "scanpy",
    "anndata",
    "pytorch_lightning",
    "scipy",
    "pandas",
    "ml_collections",
    "ray",
]
# use proxied API endpoint on rtd to avoid CORS issues
if os.environ.get("READTHEDOCS"):
    hoverxref_api_host = "/_"