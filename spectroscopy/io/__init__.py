# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Readers and writers for spectroscopy file formats.

Each format module provides ``read(filehandle, my_spectrum, **kwargs)`` and
``write(filehandle, my_spectrum, **kwargs)``.

Note: this layer must not import ``spectroscopy.spectra`` at module scope --
the dependency runs core -> io, never the other way round. Where a reader has
to construct a Spectrum (JCAMP compound files) the import is made inside the
function.

This package was previously the top-level ``formats`` package; that name still
works but is deprecated.
"""

__all__ = ['jcamp', 'csv', 'spy']
