# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Sphinx configuration."""

import datetime
from urllib.parse import quote

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

# -- Feedback (review 14.5) ------------------------------------------------
#
# The docs ask the reader to report missing questions and broken behaviour, so
# they have to say where to. The address is written down here once: the footer
# template (_templates/feedback.html) reads it from html_context, and the one
# inline link, at the end of "What do you want to know?", is built below from
# the same constant.

feedback_to = "james.sturgis@univ-amu.fr"

_question_body = (
    "The question I arrived with:\n\n"
    "The instrument or data I have:\n\n"
    "How I do this at the moment, if I do:\n\n\n"
    "-- do not delete: this says where the report came from --\n"
    "Page: what-do-you-want-to-know\n"
    f"SpectroscoPy: {release}\n"
)
_question_mailto = (
    f"mailto:{feedback_to}"
    f"?subject={quote('[SpectroscoPy] I want to know...')}"
    f"&body={quote(_question_body)}"
)

myst_substitutions = {
    "feedback_question_link": f"[tell me what it is](<{_question_mailto}>)",
}
#: Pages are executed at build time, so a code block that has gone stale
#: breaks the build instead of quietly lying to the reader.
nb_execution_mode = "auto"
nb_execution_timeout = 120
nb_execution_raise_on_error = True

autosummary_generate = True

# The build runs with -W so a genuine problem fails CI. Two classes of message
# are noise rather than problems and would otherwise mask the real ones:
#
#   mystnb   -- "Executing notebook ..." progress notes, emitted at warning
#               level; nothing is wrong.
#
# (Dataclass attributes were also being described twice, by autosummary's
# :recursive: member pages and by the module page. Dropped :recursive: in
# api.md rather than silencing it -- the per-member pages added nothing.)
suppress_warnings = ["mystnb"]
autodoc_default_options = {"members": True, "undoc-members": False,
                           "show-inheritance": True}
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False
# Render "Attributes" as :ivar: fields instead of separate object descriptions,
# which stops dataclass attributes being documented twice.
napoleon_use_ivar = True

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
templates_path = ["_templates"]
html_css_files = ["custom.css"]
#: Read by _templates/feedback.html.
html_context = {"feedback_to": feedback_to}
html_theme_options = {
    "show_prev_next": True,
    "navigation_with_keys": True,
    "show_toc_level": 2,
    "icon_links": [{
        "name": "Source",
        "url": "https://github.com/jnsturgis/SpectroscoPy",
        "icon": "fa-brands fa-square-github",
        "type": "fontawesome",
    }],
    # The feedback line sits in the middle of the footer, on every page.
    "footer_center": ["feedback"],
    "announcement": (
        "SpectroscoPy is pre-1.0 &mdash; the API may still change between "
        "releases."
    ),
}
html_sidebars = {"index": []}
