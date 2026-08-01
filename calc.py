# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Deprecated alias for :mod:`spectroscopy.lineshapes`.

``calc`` moved into the package as ``spectroscopy.lineshapes``. Importing
``calc`` still works but emits a DeprecationWarning and will be removed in 0.2.

Old:  ``import calc``            ->  ``calc.gauss(...)``
New:  ``from spectroscopy import lineshapes``  ->  ``lineshapes.gauss(...)``
"""

import warnings

from spectroscopy.lineshapes import base, gauss, lorentz, spec_comp

warnings.warn(
    "The top-level 'calc' module is deprecated; use 'spectroscopy.lineshapes' "
    "instead. 'calc' will be removed in 0.2.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ['base', 'gauss', 'lorentz', 'spec_comp']
