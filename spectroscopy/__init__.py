# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
SpectroscoPy -- a common framework for handling and analysing spectra.

⚠️  API unstable pre-1.0: breaking changes are expected between 0.x releases.

Layout (see SpectroscoPy_Codebase_Review.md):
    spectroscopy.spectra      core data model
    spectroscopy.io           format readers/writers
    spectroscopy.processing   algorithms
    spectroscopy.lineshapes   gaussian / lorentzian / voigt-ish components
    spectroscopy.cli          command-line entry points
"""

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _version

from . import datasets  # noqa: F401
from .collection import SpectrumCollection  # noqa: F401
from .history import ProcessingStep  # noqa: F401
from .messages import *  # noqa: F401,F403
from .peaks import PeakTable  # noqa: F401
from .spectra import *  # noqa: F401,F403

try:
    __version__ = _version("spectroscopy")
except PackageNotFoundError:        # running from a source tree, not installed
    __version__ = "0.0.0+unknown"
