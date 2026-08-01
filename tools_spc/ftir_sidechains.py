#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Deprecated alias for :mod:`spectroscopy.processing.ftir`.

The calculation moved to ``spectroscopy.processing.ftir`` and the command-line
front end to ``spectroscopy.cli.ftir_sidechains`` (installed as the
``spc-ftir-sidechains`` command). This shim keeps the old script path working;
it will be removed in 0.2.
"""

import warnings

from spectroscopy.cli.ftir_sidechains import main
from spectroscopy.processing.ftir import (
    SideChainData,
    calc_resid_spectrum,
    ftir_sidechain,
    get_composition,
)

warnings.warn(
    "tools_spc/ftir_sidechains.py is deprecated; use "
    "'spectroscopy.processing.ftir' (library) or the 'spc-ftir-sidechains' "
    "command (CLI). This shim will be removed in 0.2.",
    DeprecationWarning,
    stacklevel=2,
)

__all__ = ['SideChainData', 'calc_resid_spectrum', 'ftir_sidechain',
           'get_composition', 'main']


if __name__ == '__main__':
    raise SystemExit(main())
