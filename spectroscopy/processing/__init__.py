# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Processing and analysis algorithms.

Currently holds only technique-specific FTIR code (``ftir``). The
technique-agnostic algorithms (baseline, smoothing, normalisation, peak
detection and fitting) arrive in Phase 1/4 as ``processing.common``.
"""

__all__ = ['common', 'ftir', 'multivariate', 'structure']
