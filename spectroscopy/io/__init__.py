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

# Importing the format modules is what populates the registry -- each one
# registers itself with a decorator. Anything added here becomes visible to
# read_spectrum()/write_spectrum() and to Spectrum's file-type inference with
# no other change anywhere.
from spectroscopy.io import csv, dpt, jcamp, opus, spc, spy, table  # noqa: F401,E402
from spectroscopy.io.registry import (  # noqa: F401,E402
    describe_formats,
    detect_encoding,
    infer_file_type,
    known_extensions,
    known_types,
    read_spectra,
    read_spectrum,
    register_reader,
    register_writer,
    write_spectrum,
)

__all__ = [
    'jcamp', 'csv', 'dpt', 'opus', 'spc', 'spy', 'table',
    'read_spectrum', 'read_spectra', 'write_spectrum',
    'register_reader', 'register_writer',
    'known_types', 'known_extensions', 'infer_file_type', 'describe_formats',
    'detect_encoding',
]
