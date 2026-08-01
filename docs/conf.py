# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Sphinx configuration."""

import datetime

project = "SpectroscoPy"
author = "James N. Sturgis"
copyright = f"2023-{datetime.date.today().year}, {author}"   # noqa: A001

try:
    from importlib.metadata import version as _version
    release = _version("spectroscopy")
except Exception:                                            # noqa: BLE001
    release = "0.0.0"
version = release.split("+")[0]

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "matplotlib.sphinxext.plot_directive",
    "myst_nb",
    "sphinx_copybutton",
]

# -- Content ---------------------------------------------------------------

source_suffix = {".rst": "restructuredtext", ".md": "myst-nb", ".ipynb": "myst-nb"}
exclude_patterns = ["_build", "**.ipynb_checkpoints"]

myst_enable_extensions = ["colon_fence", "deflist", "substitution"]
#: Pages are executed at build time, so a code block that has gone stale
#: breaks the build instead of quietly lying to the reader.
nb_execution_mode = "auto"
nb_execution_timeout = 120
nb_execution_raise_on_error = True

autosummary_generate = True
autodoc_default_options = {"members": True, "undoc-members": False,
                           "show-inheritance": True}
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "sklearn": ("https://scikit-learn.org/stable/", None),
}

# -- Figures ---------------------------------------------------------------
#
# "Show the figure": the output of this library is pictures, so a reader should
# be able to scan a page for the picture that looks like what they want before
# reading a line of code. plot_directive runs the snippet and inserts the
# result, so a figure can never drift out of step with the code above it.

plot_include_source = True
plot_html_show_source_link = False
plot_html_show_formats = False
plot_formats = [("png", 110)]
plot_rcparams = {"figure.figsize": (7.0, 3.2), "figure.autolayout": True,
                 "font.size": 9, "savefig.bbox": "tight"}
plot_apply_rcparams = True
plot_pre_code = (
    "import numpy as np\n"
    "import matplotlib.pyplot as plt\n"
    "import spectroscopy as spc\n"
    "from spectroscopy import viz\n"
)

# -- HTML ------------------------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_title = f"SpectroscoPy {version}"
html_static_path = ["_static"]
html_theme_options = {
    "show_prev_next": True,
    "navigation_with_keys": True,
    "show_toc_level": 2,
    "icon_links": [{
        "name": "Source",
        "url": "https://github.com/jnsturgis/Spectroscopy1.0",
        "icon": "fa-brands fa-square-github",
        "type": "fontawesome",
    }],
    "announcement": (
        "SpectroscoPy is pre-1.0 &mdash; the API may still change between "
        "releases."
    ),
}
html_sidebars = {"index": []}
