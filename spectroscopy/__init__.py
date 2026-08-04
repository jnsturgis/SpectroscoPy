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

from importlib.metadata import PackageNotFoundError as _PackageNotFoundError
from importlib.metadata import version as _version

from . import datasets, fitting, io, library, lineshapes, processing, units
from .collection import SpectrumCollection
from .fitting import FitResult
from .history import ProcessingStep
from .peaks import PeakTable
from .spectra import Spectrum

#: The public surface. Everything reachable as ``spectroscopy.<name>`` that is
#: not in this list is an implementation detail, and 1.0 does not promise to
#: keep it. Submodules stay importable either way (``import spectroscopy.io``
#: always works); what this list settles is what ``spc.<name>`` offers and what
#: the freeze covers. Roadmap section 14.2 has the audit that produced it.
__all__ = [
    # -- core objects --------------------------------------------------------
    'Spectrum',              # one spectrum: x, y, units, metadata, history
    'SpectrumCollection',    # many, with grouping and batch operations
    'PeakTable',             # what find_peaks() returns
    'FitResult',             # what fit_peaks() returns
    'ProcessingStep',        # what .history is made of
    # -- submodules ----------------------------------------------------------
    'datasets',              # bundled example data, so the docs run anywhere
    'io',                    # readers, writers, and the format registry
    'library',               # reference spectra and extinction coefficients
    'processing',            # common, ftir, multivariate, unmix
    'units',                 # unit names, conversion, axis quantities
    'lineshapes',            # gaussian / lorentzian components
    'viz',                   # plotting -- imported lazily, see below
]

#: ``__version__`` is public but deliberately not in __all__: ``from
#: spectroscopy import *`` would then bind it in the caller's namespace and
#: shadow their own.


def __getattr__(name):
    """Import :mod:`spectroscopy.viz` on first use rather than at import time.

    ``viz`` pulls in matplotlib, which costs about 350 ms -- most of a second
    added to every ``import spectroscopy``, including the ones that only read a
    file. Deferring it keeps the import cheap for people who never plot, while
    ``spc.viz.plot(...)`` still works without a separate import line.
    """
    if name == 'viz':
        # import_module, not ``from . import viz``: the latter goes through
        # _handle_fromlist, which calls getattr() on this package, which calls
        # this function again -- straight into infinite recursion.
        import importlib  # noqa: PLC0415
        return importlib.import_module('.viz', __name__)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


try:
    __version__ = _version("spectroscopy")
except _PackageNotFoundError:        # running from a source tree, not installed
    __version__ = "0.0.0+unknown"
