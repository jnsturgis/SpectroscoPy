# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Created on Wed June 13 2025

@author: James STURGIS

The ``.spy`` native format: complete, lossless saving and restoring of a
Spectrum **or a whole collection**, including processing history.

Format 1.0
----------
A JSON header followed by tab-separated data, so the numbers stay greppable and
plottable with ordinary tools while the metadata is structured enough to be
read back exactly::

    # spy format 1.0
    # header
    {"name": "Glucose", "technique": "ATR-FTIR", "x_unit": "cm^-1", ...,
     "metadata": {...}, "history": [{"name": "crop", "params": {...}}, ...]}
    # data
    Wavenumber (cm^-1)	Absorbance
    900.00000	0.1234500000
    ...

The header is one JSON object, pretty-printed over as many lines as it needs
and terminated by the ``# data`` marker.

A collection: one format, not two
---------------------------------
A set of spectra is the same file with a ``# collection`` block in front and
the spectrum blocks repeated (James, 2026-08-13 -- one format that holds either,
rather than a second format beside it)::

    # spy format 1.0
    # collection
    {"name": "SMP180", "kind": "ReferenceSet",
     "info": {"source": "DichroWebGit SMP180", "licence": "MIT ...",
              "categories": [{"name": "helix", "states": ["G", "H", "I"]}]}}
    # spectrum
    {"name": "LACY", "technique": "CD", ..., "metadata": {...}}
    # data
    Wavelength (nm)	Delta epsilon
    240.00000	-0.1230000000
    # spectrum
    ...

``# spectrum`` and ``# header`` are the same marker, so **every file ever
written by this library still reads**, and a single spectrum is still written
exactly as it was before -- byte for byte. There is one function that writes a
spectrum block and both paths call it, so the two cannot drift apart.

**Why a set needs more than its spectra.** A collection has facts of its own
and there was nowhere to put them (ADR-0004 section 5). A reference set fetched
from DichroWebGit arrives under the MIT licence with a citation condition
attached; saved as 128 separate files the obligation is smeared across 128
copies of a string, or lost. A titration has a parameter that is *potential in
mV*, which is true of the series and not of any one spectrum. Those live in
``# collection``, once, and nothing is written twice -- a format that copied
``licence`` onto every spectrum would reintroduce exactly the disagreement
ADR-0004 removed.

No version bump. The version in the first line would not have protected anyone
anyway: the reader accepts any ``1.x`` and would have parsed a ``1.1``
collection as one spectrum with every block's numbers run together. What
protects a caller is the marker, which is checked.

Format 0.0 (legacy) is still read. It wrote name and axis labels but its reader
never parsed them back, so a round trip silently lost the spectrum's identity
-- review defect D4. Such files now load with whatever they can supply.

Writing always produces 1.0.
"""

# pylint: disable=W0718
import ast
import json
import warnings

import numpy as np

from spectroscopy.io.registry import (
    register_collection_writer,
    register_reader,
    register_writer,
)

CURRENT_VERSION = "1.0"

#: Simple attributes carried through the header, name -> default when absent.
_HEADER_FIELDS = {
    'name': 'unnamed',
    'technique': None,
    'x_quantity': 'Wavelength',
    'x_unit': 'nm',
    'y_quantity': 'Absorbance',
    'y_unit': 'absorbance',
    'x_label_override': None,
    'y_label_override': None,
}


def header_of(spectrum):
    """
    The JSON-ready header of one spectrum: everything but the numbers.

    Shared with the ``.spyc`` collection format, which embeds spectrum blocks
    of exactly this shape. One function writes a spectrum header and both
    formats call it, so the two cannot drift apart -- which is the failure the
    registry itself was built to stop (four tables kept in step by hand, defect
    D5).
    """
    header = {}
    for field, default in _HEADER_FIELDS.items():
        if field in ('x_label_override', 'y_label_override'):
            header[field] = getattr(spectrum, f'_{field}', None)
        else:
            header[field] = getattr(spectrum, field, default)
    header['metadata'] = spectrum.metadata
    header['history'] = [step.to_dict()
                         for step in getattr(spectrum, 'history', [])]
    return header


def apply_header(spectrum, header):
    """Set everything but the numbers on ``spectrum`` from a parsed header."""
    from spectroscopy.history import ProcessingStep  # pylint: disable=C0415

    for field, default in _HEADER_FIELDS.items():
        value = header.get(field, default)
        if field in ('x_label_override', 'y_label_override'):
            setattr(spectrum, f'_{field}', value)
        else:
            setattr(spectrum, field, value)
    spectrum.metadata = dict(header.get('metadata', {}))
    spectrum.history = [ProcessingStep.from_dict(step)
                        for step in header.get('history', [])]


def write_data(filehandle, spectrum):
    """The tab-separated block: a label row, then x and y. Shared with .spyc."""
    filehandle.write(f'{spectrum.x_label}\t{spectrum.y_label}\n')
    for x, y in zip(spectrum.x, spectrum.y):
        filehandle.write(f'{x:.5f}\t{y:.10f}\n')


def read_data(data_lines):
    """``(x, y)`` from a tab-separated block, skipping the label row."""
    xs, ys = [], []
    for line in data_lines:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            continue
        fields = stripped.split('\t')
        try:
            xs.append(float(fields[0]))
            ys.append(float(fields[1]))
        except (ValueError, IndexError):
            continue                                  # the column-label row
    return np.array(xs), np.array(ys)


def _detect_version(first_line):
    """Read the format version out of the first line, defaulting to 0.0."""
    text = first_line.strip().lstrip('#').strip().lower()
    if text.startswith('spy'):
        for part in text.split():
            if part and part[0].isdigit():
                return part
    return "0.0"


#: Markers that open a block. ``header`` and ``spectrum`` are the same thing
#: under two names: the first is what every file written before collections
#: existed uses, and dropping it would orphan them.
_SPECTRUM_MARKERS = ('header', 'spectrum')


def _blocks(lines):
    """
    Split a 1.x body into ``(marker, header_text, data_lines)`` triples.

    One pass, driven by the ``# collection`` / ``# spectrum`` / ``# header`` /
    ``# data`` markers. A block ends where the next begins.
    """
    blocks, marker, header, data, section = [], None, [], [], None
    for line in lines:
        stripped = line.strip().lower()
        if stripped.startswith('#'):
            tag = stripped.lstrip('#').strip()
            if tag in _SPECTRUM_MARKERS or tag == 'collection':
                if marker is not None:
                    blocks.append((marker, "".join(header), data))
                marker = 'collection' if tag == 'collection' else 'spectrum'
                header, data, section = [], [], 'header'
                continue
            if tag == 'data':
                section = 'data'
                continue
        if section == 'header':
            header.append(line)
        elif section == 'data':
            data.append(line)
    if marker is not None:
        blocks.append((marker, "".join(header), data))
    return blocks


def _read_v1(lines):
    """
    Parse a 1.x file into ``(spectra, collection_header)``.

    Returns a list because a file may hold one spectrum or a hundred and
    twenty-eight, and the caller decides which it wanted.
    """
    from spectroscopy.spectra import Spectrum  # pylint: disable=C0415

    spectra, collection = [], None
    for marker, header_text, data_lines in _blocks(lines[1:]):
        header = json.loads(header_text or "{}")
        if marker == 'collection':
            if collection is not None:
                raise ValueError(
                    "this file has two '# collection' blocks; a file describes "
                    "one set, and there is no way to tell which spectra belong "
                    "to which"
                )
            collection = header
            continue
        spectrum = Spectrum()
        apply_header(spectrum, header)
        spectrum.x, spectrum.y = read_data(data_lines)
        spectra.append(spectrum)
    return spectra, collection


def _read_v0(lines, my_spectrum):
    """
    Parse the legacy 0.0 layout: name, label row, '#', data, '#', metadata.

    The original reader skipped straight to the data and dropped the name and
    labels the writer had put there; they are recovered here.
    """
    xs, ys = [], []
    metadata_text = ""
    section = 'header'
    header_rows = []

    for line in lines[1:]:
        stripped = line.rstrip('\r\n')
        if section == 'header':
            if stripped.startswith('#'):
                section = 'body'
            else:
                header_rows.append(stripped)
        elif section == 'body':
            if stripped.startswith('#'):
                section = 'meta'
            else:
                fields = stripped.split('\t')
                if len(fields) >= 2:
                    xs.append(float(fields[0]))
                    ys.append(float(fields[1]))
        else:
            metadata_text += ' ' + stripped.strip()

    if header_rows:
        my_spectrum.name = header_rows[0]
    if len(header_rows) > 1 and '\t' in header_rows[1]:
        x_label, y_label = header_rows[1].split('\t')[:2]
        my_spectrum.x_label = x_label
        my_spectrum.y_label = y_label

    my_spectrum.x = np.array(xs)
    my_spectrum.y = np.array(ys)
    if metadata_text.strip():
        try:
            my_spectrum.metadata = ast.literal_eval(metadata_text.strip())
        except (ValueError, SyntaxError):
            my_spectrum.metadata = {}


# -- the collection block ----------------------------------------------------
#
# ``info`` is stored as JSON, for the same reason ``metadata`` is (ADR-0004
# section 2.5): what a format cannot represent, it degrades silently. The one
# value in practice that JSON has no form for is a list of Category, and it is
# also the one that must survive -- a category is a name *plus the DSSP states
# it claims*, and 'helix' alone does not say whether it covers 3-10 and pi.

#: Collection classes this format can restore, by the name written into
#: ``kind``. A small explicit table rather than an import by name: a file
#: should not be able to name an arbitrary class and have it constructed.
#:
#: A class earns an entry here when it **reads data the base class stores but
#: does not interpret**. ``ReferenceSet`` qualifies: its ``compositions`` reads
#: each spectrum's ``metadata['composition']`` against the set's
#: ``info['categories']``, and a plain collection does neither. A titration
#: does *not* qualify and is not a class -- it is a ``SpectrumCollection`` with
#: a parameter per spectrum and its name and unit in ``info``, both of which
#: the base class already reads. That rule is what stops this table growing by
#: habit.
COLLECTION_KINDS = ('SpectrumCollection', 'ReferenceSet')


def _resolve_kind(kind):
    """The class named by ``kind``, or the base class with a warning."""
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415

    if kind in (None, 'SpectrumCollection'):
        return SpectrumCollection
    if kind == 'ReferenceSet':
        from spectroscopy.library import ReferenceSet  # pylint: disable=C0415
        return ReferenceSet
    warnings.warn(
        f"this file was written by a version that knows a collection kind "
        f"{kind!r}; this one knows {list(COLLECTION_KINDS)}. The spectra and "
        f"the set-level info load as a plain SpectrumCollection, so nothing "
        f"is lost from the file -- but whatever {kind!r} added on top of them "
        f"is not restored.",
        stacklevel=4,
    )
    return SpectrumCollection


def _encode_info(info):
    encoded = {}
    for key, value in info.items():
        if key == 'categories':
            encoded[key] = [{'name': category.name,
                             'states': sorted(category.states),
                             'note': category.note}
                            for category in value]
        else:
            encoded[key] = value
    return encoded


def _decode_info(info):
    from spectroscopy.processing.structure import Category  # pylint: disable=C0415

    decoded = {}
    for key, value in info.items():
        if key == 'categories':
            decoded[key] = [Category(entry['name'],
                                     frozenset(entry.get('states', ())),
                                     entry.get('note'))
                            for entry in value]
        else:
            decoded[key] = value
    return decoded


# -- reading and writing -----------------------------------------------------

@register_reader('spy', extensions=['.spy'], multi=True,
                 description='native format, one spectrum or a whole collection')
def read(filehandle, **kwargs):
    """
    Read a .spy file, returning a collection.

    Every file is read as a set: one spectrum gives a set of one, which is the
    honest answer and is what lets a caller treat both the same way.
    :func:`~spectroscopy.io.registry.read_spectrum` is the singular, and keeps
    its own contract of refusing a file that holds several.

    The version is detected from the file itself. A ``format`` keyword is
    accepted for backwards compatibility but ignored -- trusting the caller
    over the file is how a 1.0 file could end up parsed as 0.0.
    """
    _ = kwargs
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415
    from spectroscopy.spectra import Spectrum  # pylint: disable=C0415

    lines = list(filehandle)
    if not lines:
        raise ValueError("empty .spy file")

    version = _detect_version(lines[0])
    if version.startswith("1."):
        spectra, collection = _read_v1(lines)
    elif version.startswith("0."):
        spectrum = Spectrum()
        _read_v0(lines, spectrum)
        spectra, collection = [spectrum], None
    else:
        raise ValueError(
            f"unsupported .spy format version {version!r}; this build writes "
            f"{CURRENT_VERSION} and reads 0.0 and 1.x"
        )

    if collection is None:
        return SpectrumCollection(spectra)
    return _resolve_kind(collection.get('kind'))(
        spectra,
        name=collection.get('name'),
        info=_decode_info(collection.get('info', {})),
    )


def _write_spectrum_block(filehandle, spectrum, marker):
    filehandle.write(f'# {marker}\n')
    filehandle.write(json.dumps(header_of(spectrum), indent=1,
                                default=str) + '\n')
    filehandle.write('# data\n')
    write_data(filehandle, spectrum)


@register_writer('spy')
def write(filehandle, my_spectrum, **kwargs):
    """
    Write ``my_spectrum`` as a .spy 1.0 file.

    Byte-for-byte what this has always written: ``# header`` rather than
    ``# spectrum``, and no collection block. Collections came later and cost
    single-spectrum files nothing.
    """
    _ = kwargs
    filehandle.write(f'# spy format {CURRENT_VERSION}\n')
    _write_spectrum_block(filehandle, my_spectrum, 'header')


@register_collection_writer('spy')
def write_collection(filehandle, collection, **kwargs):
    """
    Write a whole collection as one .spy 1.0 file.

    The set's own facts go in the ``# collection`` block, once. Per-spectrum
    data stays in the spectrum blocks. Nothing is written at both levels --
    duplicating ``licence`` onto every spectrum would put back exactly the
    disagreement ADR-0004 removed.
    """
    _ = kwargs
    filehandle.write(f'# spy format {CURRENT_VERSION}\n')
    filehandle.write('# collection\n')
    filehandle.write(json.dumps({
        'name': collection.name,
        'kind': type(collection).__name__,
        'n_spectra': len(collection),
        'info': _encode_info(getattr(collection, 'info', {})),
    }, indent=1, default=str) + '\n')

    for spectrum in collection:
        _write_spectrum_block(filehandle, spectrum, 'spectrum')


## ============================================================================

def main():
    """
    A main routine to do more or less nothing!
    """
    print("This file provides routines for reading and writing spy files")
    return True

if __name__ == '__main__':
    main()
