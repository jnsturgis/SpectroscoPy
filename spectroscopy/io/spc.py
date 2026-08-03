# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Galactic / Thermo ``.spc`` files.

Unlike OPUS (:mod:`spectroscopy.io.opus`), this format is **documented**:
Galactic published the "Universal Data Format Specification", which reproduces
the ``SPC.H`` header in full. This reader was written from that specification
and then checked against twelve files written by Galactic's own software.
``SPC_Format_Notes.md`` records the layout, the checking, and where the vendor
documentation is wrong -- its printed byte offsets do not add up, so the ones
used here were computed from the field sequence instead.

Structure::

    offset 0    512-byte main header (fversn 0x4B); 256-byte for old 0x4D
    offset 512  optional float32 x array, when TXVALS and not TXYXYS
                then, per subfile: 32-byte subheader + y payload

Three storage layouts, chosen by two bits of ``ftflgs``:

===========  ==========================================================
``gx-y``     no flags -- x generated from ``ffirst``/``flast``/``fnpts``
``x-y``      ``TXVALS`` -- one shared float32 x array follows the header
``-xy``      ``TXYXYS`` -- every subfile carries its own x and y
===========  ==========================================================

Two traps, both caught by real files rather than by reading:

* **y is usually not float.** Unless the exponent is ``0x80`` the values are
  block-scaled integers needing ``2**exp * y / 2**32``. Read naively as
  float32 a real mass spectrum here comes out around 8e15.
* **When ``TXYXYS`` is set, ``fnpts`` is not a point count** -- it is the byte
  offset of the subfile directory, and ``0`` means there is no directory at
  all, the subfiles simply follow the header.
"""

from __future__ import annotations

import struct

import numpy as np

from spectroscopy.io.registry import register_reader

__all__ = ['read_spc', 'read_spc_blocks', 'is_spc_file', 'SPC_VERSIONS']

#: ``fversn`` values and what they mean. Only 0x4B is implemented; the other
#: two are recognised so that they can be refused with a useful message.
SPC_VERSIONS = {
    0x4B: 'new format, LSB first',
    0x4C: 'new format, MSB first',
    0x4D: 'old LabCalc format',
}

_HEADER = "<cbbbiddicccci9s9sh32s130s30siicchf48sfifc187s"
_HEADER_FIELDS = (
    'ftflgs fversn fexper fexp fnpts ffirst flast fnsub fxtype fytype fztype '
    'fpost fdate fres fsource fpeakpt fspare fcmnt fcatxt flogoff fmods '
    'fprocs flevel fsampin ffactor fmethod fzinc fwplanes fwinc fwtype freserv'
).split()
_HEADER_SIZE = 512

_SUBHEADER = "<bbhfffiif4s"
_SUBHEADER_FIELDS = (
    'subflgs subexp subindx subtime subnext subnois subnpts subscan '
    'subwlevel subresv'
).split()
_SUBHEADER_SIZE = 32

_LOG_HEADER = "<iiiii44s"
_LOG_HEADER_SIZE = 64

#: ``ftflgs`` bits, from ``SPC.H``.
_TSPREC, _TCGRAM, _TMULTI = 0x01, 0x02, 0x04
_TRANDM, _TORDRD, _TALABS = 0x08, 0x10, 0x20
_TXYXYS, _TXVALS = 0x40, 0x80

#: The exponent that means "already IEEE float" -- 0x80 read as a signed char.
_FLOAT_EXPONENT = -128

#: ``fxtype``/``fztype``/``fwtype`` -> (quantity, unit). Units are the
#: project's own vocabulary (:mod:`spectroscopy.units`), so that a wavelength
#: axis is convertible and an arbitrary one is honestly labelled 'a.u.'.
_X_UNITS = {
    0:  ('Arbitrary', 'a.u.'),        16: ('Diode number', 'a.u.'),
    1:  ('Wavenumber', 'cm^-1'),      17: ('Channel', 'a.u.'),
    2:  ('Wavelength', 'um'),         18: ('Angle', 'degrees'),
    3:  ('Wavelength', 'nm'),         19: ('Temperature', 'F'),
    4:  ('Time', 's'),                20: ('Temperature', 'C'),
    5:  ('Time', 'min'),              21: ('Temperature', 'K'),
    6:  ('Frequency', 'Hz'),          22: ('Data points', 'a.u.'),
    7:  ('Frequency', 'kHz'),         23: ('Time', 'ms'),
    8:  ('Frequency', 'MHz'),         24: ('Time', 'us'),
    9:  ('Mass', 'm/z'),              25: ('Time', 'ns'),
    10: ('Chemical shift', 'ppm'),    26: ('Frequency', 'GHz'),
    11: ('Time', 'days'),             27: ('Distance', 'cm'),
    12: ('Time', 'years'),            28: ('Distance', 'm'),
    13: ('Raman shift', 'cm^-1'),     29: ('Distance', 'mm'),
    14: ('Photon energy', 'eV'),      30: ('Time', 'hours'),
    15: ('Label', 'a.u.'),            255: ('Double interferogram', 'a.u.'),
}

#: ``fytype`` -> (quantity, unit). Code 128 is '%T' rather than a 0-1
#: transmittance: the vendor's own IR sample runs to 105, which only makes
#: sense as percent.
_Y_UNITS = {
    0:  ('Intensity', 'a.u.'),        13: ('Relative intensity', 'a.u.'),
    1:  ('Interferogram', 'a.u.'),    14: ('Energy', 'a.u.'),
    2:  ('Absorbance', 'absorbance'), 16: ('Level', 'dB'),
    3:  ('Kubelka-Munk', 'a.u.'),     19: ('Temperature', 'F'),
    4:  ('Counts', 'counts'),         20: ('Temperature', 'C'),
    5:  ('Potential', 'V'),           21: ('Temperature', 'K'),
    6:  ('Angle', 'degrees'),         22: ('Refractive index', 'a.u.'),
    7:  ('Current', 'mA'),            23: ('Extinction coefficient', 'a.u.'),
    8:  ('Distance', 'mm'),           24: ('Real', 'a.u.'),
    9:  ('Potential', 'mV'),          25: ('Imaginary', 'a.u.'),
    10: ('log(1/R)', 'a.u.'),         26: ('Complex', 'a.u.'),
    11: ('Percent', '%'),             128: ('Transmittance', '%T'),
    12: ('Intensity', 'a.u.'),        129: ('Reflectance', 'reflectance'),
                                      130: ('Intensity', 'a.u.'),
    131: ('Emission', 'a.u.'),
}

#: ``fexper`` -> instrument technique. Note there is no code 6: the
#: enumeration in SPC.H skips it, and a reader that indexes a plain list here
#: mislabels everything from 7 upward.
_EXPERIMENTS = {
    0: 'General', 1: 'Gas chromatogram', 2: 'General chromatogram',
    3: 'HPLC chromatogram', 4: 'FT-IR/FT-NIR/FT-Raman', 5: 'NIR',
    7: 'UV-VIS', 8: 'X-ray diffraction', 9: 'Mass spectrum',
    10: 'NMR', 11: 'Raman', 12: 'Fluorescence', 13: 'Atomic',
    14: 'Chromatography diode array',
}

#: ``fmods`` bits -- what was done to the spectrum before it was saved. This
#: is the .spc counterpart of the OPUS history block, and worth keeping for
#: the same reason: section 21 turned up files that had already been
#: reference-subtracted, which nothing in the exported text revealed.
_MODIFICATIONS = {
    1: 'averaging', 2: 'baseline correction', 3: 'interferogram to spectrum',
    4: 'derivative', 6: 'resolution enhancement', 9: 'interpolation',
    14: 'noise reduction', 15: 'arithmetic', 19: 'spectral subtraction',
    20: 'truncation', 23: 'collection time modified', 24: 'x units conversion',
    25: 'y units conversion', 26: 'zap',
}

#: x-axis codes that identify a technique well enough to be worth recording.
#: fexper is not used for this -- see :func:`read_spc_blocks`.
_TECHNIQUE_FOR_X = {1: 'FTIR', 13: 'Raman'}


def is_spc_file(path) -> bool:
    """
    True if ``path`` looks like a Galactic SPC file.

    Sniffs rather than trusting the extension, because ``.spc`` is shared with
    Bruker EPR data -- which has no header at all, so anything that happens to
    sit in its second byte would otherwise be read as a version number.
    """
    try:
        with open(path, 'rb') as handle:
            head = handle.read(2)
        return len(head) == 2 and head[1] in SPC_VERSIONS
    except OSError:
        return False


def _bytes_of(source):
    """Accept a path, an open binary handle, or the bytes themselves."""
    if isinstance(source, (bytes, bytearray)):
        return bytes(source)
    if hasattr(source, 'read'):
        data = source.read()
        if isinstance(data, str):
            raise ValueError(
                "SPC is a binary format and this handle is in text mode; "
                "open it with mode='rb'."
            )
        return data
    with open(source, 'rb') as handle:
        return handle.read()


def _text(raw):
    """A NUL-terminated fixed-width field as a stripped string."""
    return raw.split(b'\x00')[0].decode('latin-1').strip()


def _header(raw):
    """Decode the 512-byte main header, or explain why we cannot."""
    if len(raw) < 2:
        raise ValueError(f"too short to be any kind of file: {len(raw)} bytes")

    # Identity before size. A truncated OPUS file is a broken SPC file, but a
    # Bruker EPR file is not an SPC file at all, and saying "too short" about
    # one sends the user off to check their disk instead of their format.
    version = raw[1]
    if version not in SPC_VERSIONS:
        raise ValueError(
            f"not a Galactic SPC file: byte 1 is 0x{version:02X}, and a "
            f"version byte must be one of "
            f"{', '.join(f'0x{key:02X}' for key in SPC_VERSIONS)}. "
            "Note that '.spc' is also used by Bruker for EPR data, which is a "
            "different format entirely -- a headerless block of floats with a "
            "separate .par descriptor. If this file has a .par beside it, "
            "that is what it is."
        )
    if version != 0x4B:
        raise NotImplementedError(
            f"{SPC_VERSIONS[version]} (0x{version:02X}) is not supported. "
            "Only the common little-endian format (0x4B) has been implemented, "
            "because no file of this kind was available to test against. "
            "Please report this file -- with it the support is small."
        )
    if len(raw) < _HEADER_SIZE:
        raise ValueError(
            f"too short to be an SPC file: {len(raw)} bytes, and the header "
            f"alone is {_HEADER_SIZE}"
        )

    values = dict(zip(_HEADER_FIELDS,
                      struct.unpack_from(_HEADER, raw, 0)))
    values['ftflgs'] = ord(values['ftflgs'])
    for key in ('fxtype', 'fytype', 'fztype', 'fpost', 'fprocs', 'flevel',
                'fwtype'):
        values[key] = ord(values[key])
    for key in ('fres', 'fsource', 'fmethod'):
        values[key] = _text(values[key])
    values['fcmnt'] = _text(values['fcmnt'])
    return values


def _axis_labels(header):
    """(x, y) as (quantity, unit), honouring custom text labels."""
    x_label = _X_UNITS.get(header['fxtype'], ('Arbitrary', 'a.u.'))
    y_label = _Y_UNITS.get(header['fytype'], ('Intensity', 'a.u.'))
    if header['ftflgs'] & _TALABS:
        # fcatxt holds NUL-separated x, y, z strings. Only the quantity is
        # overridden: a free-text label says what was measured, not in what
        # unit, and inventing a unit from it would defeat the point of
        # storing them separately.
        parts = header['fcatxt'].split(b'\x00')
        custom_x = parts[0].decode('latin-1').strip() if parts else ''
        custom_y = (parts[1].decode('latin-1').strip()
                    if len(parts) > 1 else '')
        if custom_x:
            x_label = (custom_x, x_label[1])
        if custom_y:
            y_label = (custom_y, y_label[1])
    return x_label, y_label


def _decode_y(raw, offset, points, exponent, sixteen_bit):
    """
    Y values, undoing the block scaling.

    ``exponent`` of ``0x80`` (-128 signed) means the payload is already IEEE
    float. Otherwise it is a fixed-point integer to be scaled by
    ``2**exponent / 2**bits`` -- get this wrong and the whole spectrum is off
    by a power of two, which looks like a units problem rather than a bug.
    """
    if exponent == _FLOAT_EXPONENT:
        return np.frombuffer(raw, dtype='<f4', count=points,
                             offset=offset).astype(float)
    dtype, bits = ('<i2', 16) if sixteen_bit else ('<i4', 32)
    values = np.frombuffer(raw, dtype=dtype, count=points, offset=offset)
    return values.astype(float) * (2.0 ** exponent) / (2.0 ** bits)


def _log(raw, header):
    """The ASCII log block as a dict, plus whatever would not split on '='."""
    offset = header['flogoff']
    if not offset or offset + _LOG_HEADER_SIZE > len(raw):
        return {}, []
    size, _memory, text_offset, _binary, _disk, _spare = struct.unpack_from(
        _LOG_HEADER, raw, offset)
    start = offset + text_offset
    text = raw[start:start + max(size - text_offset, 0)]
    entries, other = {}, []
    for record in text.replace(b'\r', b'').split(b'\n'):
        line = record.split(b'\x00')[0].decode('latin-1').strip()
        if not line:
            continue
        if '=' in line:
            key, value = line.split('=', 1)
            entries[key.strip()] = value.strip()
        else:
            other.append(line)
    return entries, other


def _modifications(fmods):
    """``fmods`` decoded into a list of names."""
    return [name for bit, name in sorted(_MODIFICATIONS.items())
            if fmods & (1 << bit)]


def read_spc_blocks(source):
    """
    Every spectrum in an SPC file, as ``(x, y, subheader, header)``.

    The low-level entry point; :func:`read_spc` is the one the registry uses.

    Notes
    -----
    The technique is taken from ``fxtype``, never from ``fexper``. All twelve
    vendor sample files report ``fexper = 0`` ("General") regardless of what
    they actually contain -- ``SPC.H`` explains that older software required
    ``TCGRAM`` to be set before ``fexper`` meant anything, and evidently most
    writers never bothered.
    """
    raw = _bytes_of(source)
    header = _header(raw)
    flags = header['ftflgs']
    multifile = bool(flags & _TMULTI)
    own_x = bool(flags & _TXYXYS)
    shared_x = bool(flags & _TXVALS)
    sixteen_bit = bool(flags & _TSPREC)
    points = header['fnpts']
    y_width = 2 if sixteen_bit else 4

    position = _HEADER_SIZE
    global_x = None
    if shared_x and not own_x:
        global_x = np.frombuffer(raw, dtype='<f4', count=points,
                                 offset=position).astype(float)
        position += 4 * points
    elif not shared_x:
        global_x = np.linspace(float(header['ffirst']),
                               float(header['flast']), points)

    starts = []
    if own_x and points:
        # For XYXY the directory lives at the byte offset that fnpts holds --
        # see the module docstring. fnpts == 0 means there is none, and the
        # subfiles simply follow the header, which is the branch below.
        for index in range(header['fnsub']):
            offset, _size, _time = struct.unpack_from(
                '<iif', raw, points + 12 * index)
            starts.append(offset)
    else:
        for _ in range(header['fnsub']):
            starts.append(position)
            block_points = points
            if own_x:
                block_points = struct.unpack_from('<i', raw,
                                                  position + 16)[0] or points
            position += _SUBHEADER_SIZE + block_points * (
                y_width + (4 if own_x else 0))

    found = []
    for start in starts:
        sub = dict(zip(_SUBHEADER_FIELDS,
                       struct.unpack_from(_SUBHEADER, raw, start)))
        block_points = sub['subnpts'] if own_x and sub['subnpts'] else points
        # The specification puts the per-subfile exponent in the subheader for
        # multifiles and in the main header otherwise. In all ten vendor files
        # tested the two agree, so this branch has never actually diverged.
        exponent = sub['subexp'] if multifile else header['fexp']
        data = start + _SUBHEADER_SIZE
        if own_x:
            x = np.frombuffer(raw, dtype='<f4', count=block_points,
                              offset=data).astype(float)
            data += 4 * block_points
        else:
            x = global_x
        y = _decode_y(raw, data, block_points, exponent, sixteen_bit)
        found.append((x, y, sub, header))
    return found


# '.cgm' is Galactic's extension for chromatograms, and the file is ordinary
# SPC -- the vendor's own gc_gasln.cgm and hplc.cgm are byte-for-byte the same
# format as the spectra. Note this extension collides with Computer Graphics
# Metafile, which is why the version byte is checked rather than the name.
@register_reader('spc', extensions=['.spc', '.cgm'],
                 description='Galactic/Thermo SPC (binary)', multi=True,
                 binary=True)
def read_spc(source):
    """
    Read a ``.spc`` file into :class:`~spectroscopy.spectra.Spectrum`.

    Returns
    -------
    list of Spectrum
        Always a list -- this is a ``multi`` reader and the registry wraps the
        result in a :class:`SpectrumCollection`. ``Spectrum.read()`` still
        hands back a single spectrum when the file holds only one.

    Notes
    -----
    What the instrument software did before saving is preserved rather than
    dropped: ``fmods`` is decoded into ``metadata['spc_modifications']`` and
    the ASCII log block, where acquisition settings live, lands in
    ``metadata`` under a ``log.`` prefix.
    """
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

    raw = _bytes_of(source)
    blocks = read_spc_blocks(raw)
    if not blocks:
        raise ValueError(
            "no spectrum found: the header declares 0 subfiles. Please report "
            "this file."
        )

    header = blocks[0][3]
    (x_quantity, x_unit), (y_quantity, y_unit) = _axis_labels(header)
    log_entries, log_other = _log(raw, header)
    modifications = _modifications(header['fmods'])
    packed = header['fdate']

    common = {f"log.{key}": value for key, value in log_entries.items()}
    common.update({
        'spc_experiment': _EXPERIMENTS.get(header['fexper'],
                                           f"unknown ({header['fexper']})"),
        'spc_version': SPC_VERSIONS[0x4B],
    })
    for key, value in (('resolution', header['fres']),
                       ('source', header['fsource']),
                       ('method', header['fmethod']),
                       ('comment', header['fcmnt'])):
        if value:
            common[key] = value
    if modifications:
        common['spc_modifications'] = modifications
    if log_other:
        common['log.text'] = "\n".join(log_other)
    if packed:
        common['collected'] = (
            f"{packed >> 20:04d}-{(packed >> 16) % 16:02d}-"
            f"{(packed >> 11) % 32:02d} {(packed >> 6) % 32:02d}:"
            f"{packed % 64:02d}"
        )

    technique = _TECHNIQUE_FOR_X.get(header['fxtype'])
    spectra = []
    for x, y, sub, _ in blocks:
        ascending = len(x) > 1 and x[0] < x[-1]
        spectrum = Spectrum(x if ascending else x[::-1],
                            y if ascending else y[::-1],
                            x_quantity=x_quantity, x_unit=x_unit,
                            y_quantity=y_quantity, y_unit=y_unit)
        spectrum.metadata.update(common)
        if technique:
            spectrum.metadata['technique_hint'] = technique
        if len(blocks) > 1:
            spectrum.name = f"subfile {sub['subindx']}"
            # The Z value is what distinguishes subfiles -- a time, a
            # temperature, a depth. Losing it makes a multifile a bag of
            # anonymous traces.
            spectrum.metadata['z_value'] = sub['subtime']
            spectrum.metadata['z_quantity'] = _X_UNITS.get(
                header['fztype'], ('Arbitrary', 'a.u.'))[0]
        if sub['subscan']:
            spectrum.metadata['scans'] = sub['subscan']
        spectra.append(spectrum)
    return spectra
