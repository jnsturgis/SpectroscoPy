# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Placeholder. This package reserves a name and does nothing else.

The library is SpectroscoPy -- one data model for Raman, FTIR, UV-Vis and
fluorescence spectra -- and it lives at https://github.com/jnsturgis/SpectroscoPy

If you reached this by installing ``pyspectroscopy`` and expecting a working
library, you have the placeholder rather than the thing itself. See the URL
above for what to install.
"""

from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _version

# Derived, not written down twice. This was hardcoded to "0.0.0" and stayed
# there when pyproject.toml went to 0.0.1, so the published 0.0.1 reports
# itself as 0.0.0 -- nothing checks the two against each other. The main
# package has always done it this way; the placeholder had drifted.
try:
    __version__ = _version("pyspectroscopy")
except _PackageNotFoundError:        # running from a source tree, not installed
    __version__ = "0.0.0+unknown"

__url__ = "https://github.com/jnsturgis/SpectroscoPy"
