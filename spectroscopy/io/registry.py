# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Format registry -- readers and writers register themselves here.

Roadmap section 3 asks for "a registry of format readers, not a chain of
if/elif on file extension". Before this there were four places that had to be
kept in step by hand: ``FILE_EXTS``, ``KNOWNFILETYPES``, and the ``match``
statements in ``Spectrum.reload()`` and ``Spectrum.save()``. They had already
drifted apart -- that was defect D5, and it cost a silently truncated file.
Now there is one table and everything else is derived from it.

Adding a format is a decorator::

    @register_reader('dpt', extensions=['.dpt'], description='Bruker OPUS')
    def read(handle, spectrum, **kwargs):
        ...

Layering: this module must not import ``spectroscopy.spectra`` at module scope
(review C1). Where it needs to build a Spectrum it imports inside the function.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from dataclasses import dataclass, field

__all__ = [
    'FormatEntry', 'register_reader', 'register_writer',
    'read_spectrum', 'read_spectra', 'write_spectrum',
    'infer_file_type', 'known_types', 'known_extensions', 'describe_formats',
]

#: Byte-order marks that identify a text encoding, longest first so that
#: UTF-32's BOM is not mistaken for UTF-16's.
_BOMS = (
    (b'\xff\xfe\x00\x00', 'utf-32-le'),
    (b'\x00\x00\xfe\xff', 'utf-32-be'),
    (b'\xef\xbb\xbf',     'utf-8-sig'),
    (b'\xff\xfe',         'utf-16-le'),
    (b'\xfe\xff',         'utf-16-be'),
)

#: Tried in order when there is no BOM. latin-1 last because it maps every
#: byte and so can never itself fail -- it is the backstop, not a guess.
_FALLBACK_ENCODINGS = ('utf-8', 'latin-1')


@dataclass
class FormatEntry:
    """One registered format."""

    name: str
    extensions: tuple = ()
    reader: Callable | None = None
    writer: Callable | None = None
    description: str = ''
    #: True when the reader returns a list of Spectrum rather than filling one.
    multi: bool = False
    #: Opened in binary mode. Text formats sniff an encoding first; a
    #: binary one must not, and would be corrupted by the attempt.
    binary: bool = False
    #: Default keyword arguments handed to the reader/writer.
    defaults: dict = field(default_factory=dict)

    @property
    def readable(self) -> bool:
        return self.reader is not None

    @property
    def writable(self) -> bool:
        return self.writer is not None


#: name -> FormatEntry. The single source of truth.
REGISTRY: dict[str, FormatEntry] = {}


def _entry(name) -> FormatEntry:
    return REGISTRY.setdefault(name, FormatEntry(name=name))


def register_reader(name, extensions=(), description='', multi=False,
                    binary=False, **defaults):
    """Decorator registering a read function for a format."""
    def decorate(function):
        entry = _entry(name)
        entry.reader = function
        entry.multi = multi
        entry.binary = binary
        entry.extensions = tuple(dict.fromkeys(entry.extensions + tuple(extensions)))
        entry.description = description or entry.description
        entry.defaults = {**entry.defaults, **defaults}
        return function
    return decorate


def register_writer(name, extensions=(), description=''):
    """Decorator registering a write function for a format."""
    def decorate(function):
        entry = _entry(name)
        entry.writer = function
        entry.extensions = tuple(dict.fromkeys(entry.extensions + tuple(extensions)))
        entry.description = description or entry.description
        return function
    return decorate


def known_types(readable=None, writable=None):
    """Registered format names, optionally filtered by capability."""
    names = []
    for name, entry in sorted(REGISTRY.items()):
        if readable is not None and entry.readable != readable:
            continue
        if writable is not None and entry.writable != writable:
            continue
        names.append(name)
    return tuple(names)


def known_extensions():
    """Map of lower-case extension -> format name."""
    mapping = {}
    for entry in REGISTRY.values():
        for extension in entry.extensions:
            mapping.setdefault(extension.lower(), entry.name)
    return mapping


def infer_file_type(path):
    """Format name implied by a path's extension, or ``'unknown'``."""
    _, extension = os.path.splitext(os.fspath(path))
    return known_extensions().get(extension.lower(), 'unknown')


def describe_formats():
    """A human-readable table of what is supported -- handy in a notebook."""
    lines = [f"{'format':<8} {'r/w':<4} {'extensions':<22} description"]
    for name in known_types():
        entry = REGISTRY[name]
        flags = ('r' if entry.readable else '-') + ('w' if entry.writable else '-')
        lines.append(f"{name:<8} {flags:<4} {', '.join(entry.extensions):<22} "
                     f"{entry.description}")
    return "\n".join(lines)


def _lookup(file_type, path, need):
    """Resolve a format name, falling back to the path's extension."""
    name = file_type or infer_file_type(path)
    entry = REGISTRY.get(name)
    if entry is None:
        raise TypeError(
            f"Unknown filetype {name!r} for {os.path.basename(str(path))}; "
            f"known types are {', '.join(known_types())}"
        )
    if need == 'read' and not entry.readable:
        raise ValueError(f"Format {name!r} cannot be read (no reader registered)")
    if need == 'write' and not entry.writable:
        raise ValueError(
            f"Cannot write {name!r} files; writable types are "
            f"{', '.join(known_types(writable=True))}"
        )
    return entry


def detect_encoding(path):
    """
    Work out how to decode a text file.

    Checks for a byte-order mark first -- the AKTA chromatography exports are
    UTF-16-LE with a BOM, and reading them as UTF-8 produces nonsense rather
    than an error. Otherwise tries UTF-8 then latin-1, because vendor files
    from older instruments carry latin-1 bytes in comment fields and refusing
    a file whose numbers are perfectly readable helps nobody.
    """
    with open(path, 'rb') as raw:
        head = raw.read(4)
    for bom, encoding in _BOMS:
        if head.startswith(bom):
            return encoding

    for encoding in _FALLBACK_ENCODINGS:
        try:
            with open(path, encoding=encoding) as handle:
                handle.read(4096)
            return encoding
        except UnicodeDecodeError:
            continue
    return _FALLBACK_ENCODINGS[-1]


def _open(path, encoding=None):
    return open(path, encoding=encoding or detect_encoding(path))


def read_spectra(path, file_type=None, **kwargs):
    """
    Read a file and return a :class:`SpectrumCollection`.

    Every format goes through here, including the single-spectrum ones -- a
    file holding one spectrum is just a collection of length 1. That matters
    because several real formats hold many: Chloe's excitation-emission export
    is 173 columns, and a JCAMP compound file has children.
    """
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415
    from spectroscopy.spectra import Spectrum  # pylint: disable=C0415

    entry = _lookup(file_type, path, 'read')
    arguments = {**entry.defaults, **kwargs}
    encoding = arguments.pop('encoding', None)

    opener = (open(path, 'rb') if entry.binary
              else _open(path, encoding))
    with opener as handle:
        if entry.multi:
            spectra = list(entry.reader(handle, **arguments))
        else:
            spectrum = Spectrum()
            entry.reader(handle, spectrum, **arguments)
            spectra = [spectrum]

    directory, filename = os.path.split(os.fspath(path))
    for spectrum in spectra:
        spectrum.fileinfo = {'PATH': directory + os.sep if directory else '',
                             'NAME': filename, 'TYPE': entry.name}
        if spectrum.name in ('unnamed', '', None):
            spectrum.name = filename
    return SpectrumCollection(spectra, name=filename)


def read_spectrum(path, file_type=None, **kwargs):
    """
    Read a file expected to hold a single spectrum.

    Raises if the file turns out to hold several, naming
    :func:`read_spectra` -- silently returning the first would be the kind of
    quiet wrong this library exists to remove.
    """
    spectra = read_spectra(path, file_type, **kwargs)
    if len(spectra) != 1:
        raise ValueError(
            f"{os.path.basename(str(path))} holds {len(spectra)} spectra; "
            f"use read_spectra() to get them all"
        )
    return spectra[0]


def write_spectrum(spectrum, path, file_type=None, **kwargs):
    """
    Write a spectrum.

    The format is resolved **before** the file is opened. ``open(..., 'w')``
    truncates immediately, so validating afterwards would already have
    destroyed the target -- that was half of defect D5.
    """
    entry = _lookup(file_type, path, 'write')
    arguments = {**entry.defaults, **kwargs}
    arguments.pop('encoding', None)
    with open(path, 'w', encoding='utf-8') as handle:
        entry.writer(handle, spectrum, **arguments)
