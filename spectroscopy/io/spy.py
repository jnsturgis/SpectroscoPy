# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Created on Wed June 13 2025

@author: James STURGIS

The ``.spy`` native format: complete, lossless saving and restoring of a
Spectrum, including its processing history.

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

Format 0.0 (legacy) is still read. It wrote name and axis labels but its reader
never parsed them back, so a round trip silently lost the spectrum's identity
-- review defect D4. Such files now load with whatever they can supply.

Writing always produces 1.0.
"""

# pylint: disable=W0718
import ast
import json

import numpy as np

from spectroscopy.io.registry import register_reader, register_writer

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


def _detect_version(first_line):
    """Read the format version out of the first line, defaulting to 0.0."""
    text = first_line.strip().lstrip('#').strip().lower()
    if text.startswith('spy'):
        for part in text.split():
            if part and part[0].isdigit():
                return part
    return "0.0"


def _read_v1(lines, my_spectrum):
    """Parse a 1.0 file: JSON header, then a label row, then x/y rows."""
    header_lines, data_lines = [], []
    section = None
    for line in lines[1:]:
        stripped = line.strip().lower()
        if stripped in ('# header', '#header'):
            section = 'header'
            continue
        if stripped in ('# data', '#data'):
            section = 'data'
            continue
        if section == 'header':
            header_lines.append(line)
        elif section == 'data':
            data_lines.append(line)

    header = json.loads("".join(header_lines) or "{}")

    for field, default in _HEADER_FIELDS.items():
        value = header.get(field, default)
        if field == 'x_label_override':
            my_spectrum._x_label_override = value    # pylint: disable=W0212
        elif field == 'y_label_override':
            my_spectrum._y_label_override = value    # pylint: disable=W0212
        else:
            setattr(my_spectrum, field, value)

    my_spectrum.metadata = dict(header.get('metadata', {}))

    from spectroscopy.history import ProcessingStep  # pylint: disable=C0415
    my_spectrum.history = [ProcessingStep.from_dict(step)
                           for step in header.get('history', [])]

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
    my_spectrum.x = np.array(xs)
    my_spectrum.y = np.array(ys)


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


@register_reader('spy', extensions=['.spy'], description='native format, carries units and processing history')
def read(filehandle, my_spectrum, **kwargs):
    """
    Read a .spy file into ``my_spectrum``.

    The version is detected from the file itself. A ``format`` keyword is
    accepted for backwards compatibility but ignored -- trusting the caller
    over the file is how a 1.0 file could end up parsed as 0.0.
    """
    _ = kwargs
    lines = list(filehandle)
    if not lines:
        raise ValueError("empty .spy file")

    version = _detect_version(lines[0])
    if version.startswith("1."):
        _read_v1(lines, my_spectrum)
    elif version.startswith("0."):
        _read_v0(lines, my_spectrum)
    else:
        raise ValueError(
            f"unsupported .spy format version {version!r}; this build writes "
            f"{CURRENT_VERSION} and reads 0.0 and 1.x"
        )


@register_writer('spy')
def write(filehandle, my_spectrum, **kwargs):
    """Write ``my_spectrum`` as a .spy 1.0 file."""
    _ = kwargs

    header = {}
    for field, default in _HEADER_FIELDS.items():
        if field == 'x_label_override':
            header[field] = getattr(my_spectrum, '_x_label_override', None)
        elif field == 'y_label_override':
            header[field] = getattr(my_spectrum, '_y_label_override', None)
        else:
            header[field] = getattr(my_spectrum, field, default)

    header['metadata'] = my_spectrum.metadata
    header['history'] = [step.to_dict()
                         for step in getattr(my_spectrum, 'history', [])]

    filehandle.write(f'# spy format {CURRENT_VERSION}\n')
    filehandle.write('# header\n')
    filehandle.write(json.dumps(header, indent=1, default=str) + '\n')
    filehandle.write('# data\n')
    filehandle.write(f'{my_spectrum.x_label}\t{my_spectrum.y_label}\n')
    for x, y in zip(my_spectrum.x, my_spectrum.y):
        filehandle.write(f'{x:.5f}\t{y:.10f}\n')


## ============================================================================

def main():
    """
    A main routine to do more or less nothing!
    """
    print("This file provides routines for reading and writing spy files")
    return True

if __name__ == '__main__':
    main()
