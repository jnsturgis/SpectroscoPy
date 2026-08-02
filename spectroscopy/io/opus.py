# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Bruker OPUS native files -- ``sample.0``, ``sample.1``, ...

OPUS is proprietary and undocumented; this reader was written against 41 real
files that also had ``.dpt`` exports of the same measurement, so every spectrum
it produces can be checked against what OPUS itself wrote out. That pairing is
the only reason it is trustworthy, and it is why roadmap section 15.2 asked for
the files before the code.

Structure, as far as is needed here::

    offset 0    magic 0xFEFE0A0A
    offset 4    program version (double)
    offset 12   pointer to the block directory
    offset 16   directory capacity, in entries
    offset 20   directory entries used

    directory   12 bytes per entry: (type, length in 4-byte words, offset)

A **data** block holds ``NPT`` float32 values. Its **parameter** block holds the
axis and the scaling, and is found by a rule that holds across every file
tested: its type is the data block's type with bit 4 set (``type | 0x10``).
Parameter blocks are ``NAM\\0`` + type + size + value, terminated by ``END``.

The y values are stored scaled: multiply by ``CSF``. The x axis is not stored
at all -- it is reconstructed from ``FXV``, ``LXV`` and ``NPT``, which is why a
file with a nonlinear axis would need more work than this.
"""

from __future__ import annotations

import struct

import numpy as np

from spectroscopy.io.registry import register_reader

__all__ = ['read_opus', 'read_opus_blocks', 'is_opus_file', 'OPUS_MAGIC']

#: First four bytes of every OPUS file.
OPUS_MAGIC = 0xFEFE0A0A

#: Bit that turns a data block's type into its parameter block's type.
_PARAMETER_BIT = 0x10

#: Data-block types worth returning, in the order they are preferred. The
#: names are OPUS's own. AB is what a user almost always wants; the single
#: channel spectra are kept because a reference is sometimes the point.
_BLOCK_NAMES = {
    0x0000100F: 'AB',        # absorbance
    0x0000140F: 'TR',        # transmittance
    0x0000180F: 'KM',        # Kubelka-Munk
    0x00001007: 'ScSm',      # single channel, sample
    0x00001008: 'IgSm',      # interferogram, sample
    0x00002007: 'ScRf',      # single channel, reference
    0x00002008: 'IgRf',      # interferogram, reference
    0x00000407: 'ScSm',
    0x0000040B: 'IgSm',
}

#: Parameter blocks describing the measurement rather than one spectrum.
_METADATA_BLOCKS = {
    0x00000020: 'instrument',
    0x00000040: 'acquisition',
    0x00000060: 'optics',
    0x000000A0: 'sample',
    0x00000030: 'fourier_transform',
    0x00000080: 'origin',
}

#: y units implied by the block name, so a reader does not have to guess.
_Y_UNITS = {
    'AB': ('Absorbance', 'absorbance'),
    'TR': ('Transmittance', 'transmittance'),
    'ScSm': ('Single channel', 'a.u.'),
    'ScRf': ('Single channel', 'a.u.'),
    'IgSm': ('Interferogram', 'a.u.'),
    'IgRf': ('Interferogram', 'a.u.'),
    'KM': ('Kubelka-Munk', 'a.u.'),
}


def is_opus_file(path) -> bool:
    """True if ``path`` starts with the OPUS magic number.

    Extension-based detection is not enough here: OPUS numbers its files by
    measurement (``sample.0``, ``sample.1``), so the extension carries no type
    information at all and a ``.0`` could be anything.
    """
    try:
        with open(path, 'rb') as handle:
            return struct.unpack('<I', handle.read(4))[0] == OPUS_MAGIC
    except (OSError, struct.error):
        return False


def _bytes_of(source):
    """Accept a path, an open binary handle, or the bytes themselves."""
    if isinstance(source, (bytes, bytearray)):
        return bytes(source)
    if hasattr(source, 'read'):
        data = source.read()
        if isinstance(data, str):
            raise ValueError(
                "OPUS is a binary format and this handle is in text mode; "
                "open it with mode='rb'."
            )
        return data
    with open(source, 'rb') as handle:
        return handle.read()


def _directory(raw):
    """(type, offset, length in bytes) for every block."""
    magic, = struct.unpack_from('<I', raw, 0)
    if magic != OPUS_MAGIC:
        raise ValueError(
            f"not an OPUS file: expected magic 0x{OPUS_MAGIC:08X}, "
            f"found 0x{magic:08X}"
        )
    pointer, _capacity, used = struct.unpack_from('<III', raw, 12)
    blocks = []
    for index in range(used):
        block_type, words, offset = struct.unpack_from('<III', raw,
                                                       pointer + 12 * index)
        blocks.append((block_type, offset, words * 4))
    return blocks


def _parameters(raw, offset, length):
    """Decode a parameter block into a plain dict."""
    values, position, end = {}, offset, offset + length
    while position < end - 8:
        name = raw[position:position + 3].decode('latin-1')
        kind, size = struct.unpack_from('<HH', raw, position + 4)
        position += 8
        if name == 'END':
            break
        payload = raw[position:position + size * 2]
        position += size * 2
        try:
            if kind == 0:
                values[name] = struct.unpack_from('<i', payload)[0]
            elif kind == 1:
                values[name] = struct.unpack_from('<d', payload)[0]
            else:
                values[name] = payload.split(b'\x00')[0].decode('latin-1')
        except struct.error:
            continue                       # a truncated entry ends the block
    return values


def read_opus_blocks(path):
    """
    Every spectrum in an OPUS file, as ``(name, x, y, parameters)``.

    The low-level entry point. :func:`read_opus` is the one the registry uses.
    """
    raw = _bytes_of(path)
    blocks = _directory(raw)
    by_type = {block_type: (offset, length) for block_type, offset, length in blocks}

    metadata = {}
    for block_type, label in _METADATA_BLOCKS.items():
        for candidate, (offset, length) in by_type.items():
            if candidate & 0x00FFFFFF == block_type:
                metadata.update({f"{label}.{key}": value for key, value
                                 in _parameters(raw, offset, length).items()})

    found = []
    for block_type, offset, length in blocks:
        status_type = block_type | _PARAMETER_BIT
        if status_type == block_type or status_type not in by_type:
            continue                       # not a data block
        status_offset, status_length = by_type[status_type]
        status = _parameters(raw, status_offset, status_length)
        points = status.get('NPT')
        if not points or points * 4 > length + 4:
            continue

        values = np.frombuffer(raw, dtype='<f4', count=points, offset=offset)
        values = values.astype(float) * float(status.get('CSF', 1.0) or 1.0)

        first, last = status.get('FXV'), status.get('LXV')
        if first is None or last is None:
            continue
        x = np.linspace(float(first), float(last), points)

        name = _BLOCK_NAMES.get(block_type & 0x00FFFFFF,
                                f"block_0x{block_type:08X}")
        found.append((name, x, values,
                      {**metadata, **status, 'opus_block_type': block_type}))
    return found


@register_reader('opus', extensions=['.0', '.1', '.2', '.3', '.4', '.5',
                                     '.6', '.7', '.8', '.9'],
                 description='Bruker OPUS native (binary)', multi=True,
                 binary=True)
def read_opus(path, block=None):
    """
    Read a native OPUS file into :class:`~spectroscopy.spectra.Spectrum`.

    Parameters
    ----------
    path : str or path-like
    block : str, optional
        Which spectrum to return -- ``'AB'``, ``'ScSm'``, ``'ScRf'`` ... The
        first absorbance-like block is used when omitted, falling back to
        whatever the file has.

    Returns
    -------
    Spectrum or list of Spectrum
        A list when the file holds several and none was named, which is the
        registry's convention for multi-spectrum files.

    Notes
    -----
    OPUS records what was done to a spectrum inside the file, and this reader
    keeps it: the acquisition parameters land in ``metadata`` with an
    ``acquisition.`` prefix, and anything OPUS did -- a subtraction, an
    atmospheric compensation -- is preserved verbatim under
    ``metadata['opus_history']`` rather than being silently lost. That matters
    more than it sounds: several of the files this reader was written against
    had already had a reference subtracted inside OPUS, which nothing in the
    exported ``.dpt`` reveals.
    """
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

    raw = _bytes_of(path)
    blocks = read_opus_blocks(raw)
    if not blocks:
        raise ValueError(
            f"{path}: no readable spectrum block. The file may hold only "
            "interferograms or parameters, or use a variant this reader does "
            "not handle -- please report it with the file."
        )

    history = _text_history(raw)

    spectra = []
    for name, x, y, parameters in blocks:
        quantity, unit = _Y_UNITS.get(name, ('Signal', 'a.u.'))
        ascending = x[0] < x[-1]
        spectrum = Spectrum(x if ascending else x[::-1],
                            y if ascending else y[::-1],
                            x_quantity='Wavenumber', x_unit='cm^-1',
                            y_quantity=quantity, y_unit=unit)
        spectrum.name = f"{parameters.get('sample.SNM', '')}".strip() or str(path)
        spectrum.name = f"{spectrum.name} [{name}]"
        spectrum.metadata.update({
            key: value for key, value in parameters.items()
            if key not in ('NPT', 'FXV', 'LXV', 'CSF', 'MXY', 'MNY', 'DPF')
        })
        spectrum.metadata['opus_block'] = name
        if history:
            spectrum.metadata['opus_history'] = history
        spectra.append(spectrum)

    # Always a list: this is a ``multi`` reader, and the registry wraps what it
    # returns into a SpectrumCollection. Spectrum.read() still hands back a
    # single spectrum when a file holds only one.
    names = [entry[0] for entry in blocks]
    if block is not None:
        chosen = [spectrum for name, spectrum in zip(names, spectra)
                  if name == block]
        if not chosen:
            available = ", ".join(sorted(set(names)))
            raise ValueError(f"{path}: no {block!r} block; it has {available}")
        return chosen

    # A file often holds the same quantity twice -- what was measured, and what
    # was left after something was done to it inside OPUS. Which one a .dpt
    # export contains depends on what the operator had selected, and is not
    # recoverable from the export. So return them all rather than guess: the
    # caller can look at opus_history and choose.
    for preferred in ('AB', 'TR', 'KM'):
        chosen = [spectrum for name, spectrum in zip(names, spectra)
                  if name == preferred]
        if chosen:
            return chosen
    return spectra


def _text_history(raw_or_path):
    """OPUS's own record of what it did, if the file carries one.

    Kept verbatim rather than parsed into ProcessingStep: it is another
    program's log, and inventing structure for it would be claiming to
    understand more of it than we do.
    """
    raw = _bytes_of(raw_or_path)
    for block_type, offset, length in _directory(raw):
        if block_type & 0x0000FF00 == 0x00003400:
            text = raw[offset:offset + length].decode('latin-1', errors='replace')
            return text.replace('\x00', ' ').strip() or None
    return None
