# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Deprecated alias for :mod:`spectroscopy.io`.

The format readers moved into the package as ``spectroscopy.io`` so that the
io layer sits below core rather than beside it. Importing ``formats`` still
works and returns the same module objects, but emits a DeprecationWarning and
will be removed in 0.2.

Old:  ``from formats import jcamp``  /  ``import formats.jcamp``
New:  ``from spectroscopy.io import jcamp``
"""

import sys
import warnings

from spectroscopy.io import csv, jcamp, spy

warnings.warn(
    "The top-level 'formats' package is deprecated; use 'spectroscopy.io' "
    "instead. 'formats' will be removed in 0.2.",
    DeprecationWarning,
    stacklevel=2,
)

# Make 'import formats.jcamp' resolve to the same module objects, so module
# state stays shared between the old and new names.
sys.modules[__name__ + '.jcamp'] = jcamp
sys.modules[__name__ + '.csv'] = csv
sys.modules[__name__ + '.spy'] = spy

__all__ = ['jcamp', 'csv', 'spy']
