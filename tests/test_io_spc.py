# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Tests for the Galactic/Thermo ``.spc`` reader.

Two kinds of test here, and the distinction matters:

* Those that need no files at all -- sniffing, refusals, error messages.
  These always run.
* Those that read Galactic's own sample files. These **skip** unless the
  files have been fetched (``python scripts/fetch_spc_fixtures.py``), because
  the files carry no stated licence and are not committed.

The second kind is the point. A reader written from a specification and tested
only against files generated from the same specification proves nothing except
that the author was consistent -- roadmap section 15.2. These check it against
data written by the software that defined the format.
"""

import struct

import numpy as np
import pytest

from spectroscopy.io import read_spectra, read_spectrum
from spectroscopy.io.registry import infer_file_type
from spectroscopy.io.spc import (
    SPC_VERSIONS,
    is_spc_file,
    read_spc,
    read_spc_blocks,
)

# --------------------------------------------------------------------------
# No files needed
# --------------------------------------------------------------------------


def test_extension_is_registered():
    assert infer_file_type('sample.spc') == 'spc'
    assert infer_file_type('SAMPLE.SPC') == 'spc'


def test_sniffing_rejects_bruker_epr(tmp_path):
    """
    The confusion a user will actually hit.

    Bruker stores EPR data in ``.spc`` too, with no header at all -- a raw
    block of floats and a separate ``.par``. The first such file starts
    ``00 a0 9b be``, so byte 1 is 0xA0, not a version number.
    """
    path = tmp_path / 'epr.spc'
    path.write_bytes(bytes([0x00, 0xA0, 0x9B, 0xBE]) * 64)
    assert not is_spc_file(path)

    with pytest.raises(ValueError, match='Bruker'):
        read_spc(path)


def test_error_names_the_par_file(tmp_path):
    """The message has to be actionable, not just correct."""
    path = tmp_path / 'epr.spc'
    path.write_bytes(bytes([0x00, 0xA0, 0x9B, 0xBE]) * 64)
    with pytest.raises(ValueError) as caught:
        read_spc(path)
    assert '.par' in str(caught.value)


def test_old_and_msb_formats_are_refused_clearly(tmp_path):
    """
    Recognised but not implemented, and the difference is visible.

    No file of either kind was available to test against, so guessing at them
    would be exactly the untested-code problem the roadmap warns about.
    """
    for version in (0x4C, 0x4D):
        path = tmp_path / f'v{version:02x}.spc'
        path.write_bytes(bytes([0x00, version]) + b'\x00' * 600)
        with pytest.raises(NotImplementedError) as caught:
            read_spc(path)
        assert SPC_VERSIONS[version] in str(caught.value)
        assert 'report' in str(caught.value)


def test_truncated_file_says_so(tmp_path):
    path = tmp_path / 'short.spc'
    path.write_bytes(bytes([0x00, 0x4B]) + b'\x00' * 20)
    with pytest.raises(ValueError, match='too short'):
        read_spc(path)


def test_text_handle_is_refused(tmp_path):
    path = tmp_path / 'x.spc'
    path.write_text('not binary')
    with pytest.raises(ValueError, match="mode='rb'"):
        with open(path, encoding='utf-8') as handle:
            read_spc(handle)


def test_is_spc_file_on_missing_path(tmp_path):
    assert not is_spc_file(tmp_path / 'nope.spc')


# --------------------------------------------------------------------------
# Round trip through a file we build ourselves
# --------------------------------------------------------------------------


def _build_spc(y, first, last, *, fexp=-128, ftflgs=0, fxtype=1, fytype=2,
               x=None):
    """
    A minimal single-subfile SPC file.

    Deliberately built from the field sequence rather than by copying a real
    file, so that it exercises the writer side of our understanding of the
    layout.
    """
    points = len(y)
    header = struct.pack(
        "<cbbbiddicccci9s9sh32s130s30siicchf48sfifc187s",
        bytes([ftflgs]), 0x4B, 0, fexp, points, float(first), float(last), 1,
        bytes([fxtype]), bytes([fytype]), bytes([0]), bytes([0]),
        0, b'', b'', 0, b'', b'', b'', 0, 0, bytes([0]), bytes([0]), 0, 0.0,
        b'', 0.0, 0, 0.0, bytes([0]), b'',
    )
    subheader = struct.pack("<bbhfffiif4s", 0, fexp, 0, 0.0, 0.0, 0.0,
                            0, 0, 0.0, b'')
    body = b''
    if x is not None:
        body += np.asarray(x, dtype='<f4').tobytes()
    payload = (np.asarray(y, dtype='<f4').tobytes() if fexp == -128
               else np.asarray(y, dtype='<i4').tobytes())
    return header + body + subheader + payload


def test_round_trip_of_a_generated_file():
    y = np.linspace(0.0, 1.0, 64)
    raw = _build_spc(y, 400.0, 4000.0)
    spectra = read_spc(raw)
    assert len(spectra) == 1
    assert np.allclose(spectra[0].y, y)
    assert np.isclose(spectra[0].x[0], 400.0)
    assert np.isclose(spectra[0].x[-1], 4000.0)


def test_fixed_point_exponent_is_applied():
    """
    The trap from the notes: a non-0x80 exponent means integers.

    Read as float32 these bytes are denormal noise; scaled properly they are
    the counts that were stored.
    """
    counts = np.array([1000, 2000, 3000, 4000], dtype='<i4')
    raw = _build_spc(counts, 0.0, 3.0, fexp=15)
    spectrum = read_spc(raw)[0]
    expected = counts * (2.0 ** 15) / (2.0 ** 32)
    assert np.allclose(spectrum.y, expected)


def test_shared_x_array_is_used_not_generated():
    """TXVALS: the stored x wins over ffirst/flast interpolation."""
    x = np.array([100.0, 101.0, 105.0, 130.0], dtype='<f4')
    y = np.array([1.0, 2.0, 3.0, 4.0])
    raw = _build_spc(y, 100.0, 130.0, ftflgs=0x80, x=x)
    spectrum = read_spc(raw)[0]
    assert np.allclose(spectrum.x, x)
    # 105 is what distinguishes a stored axis from an evenly spaced one.
    assert not np.allclose(spectrum.x, np.linspace(100.0, 130.0, 4))


# --------------------------------------------------------------------------
# Galactic's own files
# --------------------------------------------------------------------------

#: (file, points, first x, last x, x unit, y quantity). Every value read off
#: the vendor's files, not off our reader.
EXPECTED = [
    ('ir-nh4.spc',   8289, 4398.10,  401.23, 'cm^-1', 'Transmittance'),
    ('raman.spc',    1260,   96.31, 1728.00, 'cm^-1', 'Counts'),
    ('uv-holm.spc',   901,  250.00,  700.00, 'nm',    'Absorbance'),
    ('nir-poly.spc',  700, 1100.00, 2498.00, 'nm',    'Absorbance'),
    ('vis-mirr.spc',  692,  890.00,  199.00, 'nm',    'Reflectance'),
    ('nmr-unk.spc',  8192,   14.83,   -5.24, 'ppm',   'Intensity'),
    ('fluor.spc',     491,  255.00,  500.00, 'nm',    'Intensity'),
    ('gc_gasln.cgm', 54000,   0.00,   30.00, 'min',   'Potential'),
]


@pytest.mark.parametrize('name,points,first,last,x_unit,y_quantity', EXPECTED)
def test_vendor_files_read_correctly(spc_sample, name, points, first, last,
                                     x_unit, y_quantity):
    spectrum = read_spectrum(spc_sample(name))
    assert len(spectrum.x) == points
    assert len(spectrum.y) == points
    # Stored high-to-low or low-to-high; the reader normalises to ascending.
    low, high = sorted((first, last))
    assert np.isclose(spectrum.x.min(), low, atol=0.01)
    assert np.isclose(spectrum.x.max(), high, atol=0.01)
    assert spectrum.x_unit == x_unit
    assert spectrum.y_quantity == y_quantity
    assert np.all(np.isfinite(spectrum.y))


def test_x_axis_is_ascending(spc_sample):
    """An IR file stored 4398 -> 401 comes back ascending."""
    spectrum = read_spectrum(spc_sample('ir-nh4.spc'))
    assert np.all(np.diff(spectrum.x) > 0)


def test_transmittance_is_percent_not_fraction(spc_sample):
    """
    Code 128 is %T.

    The vendor's own IR sample runs above 100, which only makes sense as
    percent -- a 0-1 transmittance would be a physical impossibility.
    """
    spectrum = read_spectrum(spc_sample('ir-nh4.spc'))
    assert spectrum.y_unit == '%T'
    assert spectrum.y.max() > 100
    assert 20 < spectrum.y.min() < 30


def test_counts_come_back_as_whole_numbers(spc_sample):
    """
    The signature of a correctly undone block scaling.

    These were stored as fixed-point integers; if the exponent or the divisor
    were wrong the values would still look plausible but would not be
    integral.
    """
    spectrum = read_spectrum(spc_sample('xraydiff.spc'))
    assert np.allclose(spectrum.y, np.round(spectrum.y))
    assert spectrum.y.max() > 4000


def test_xyxy_file_without_a_directory(spc_sample):
    """
    ``ms-barb.spc``: TXYXYS set, fnsub == 1, fnpts == 0.

    fnpts is a directory offset here, and 0 means there is no directory --
    the subheader follows the header directly. A reader that always seeks to
    a directory reads the main header as directory entries.
    """
    blocks = read_spc_blocks(spc_sample('ms-barb.spc'))
    assert len(blocks) == 1
    x, y, sub, header = blocks[0]
    assert header['fnpts'] == 0
    assert sub['subnpts'] == 37
    assert len(x) == len(y) == 37
    assert np.isclose(x[0], 27.05, atol=0.01)
    assert np.isclose(x[-1], 160.90, atol=0.01)
    assert np.all(np.diff(x) > 0)


def test_mass_spectrum_intensities_are_integers(spc_sample):
    """subexp = 15, so these are scaled integers -- and recover as integers."""
    _, y, _, _ = read_spc_blocks(spc_sample('ms-barb.spc'))[0]
    assert np.allclose(y, np.round(y))
    assert np.isclose(y.max(), 11798.0)


def test_base_peak_is_where_chemistry_says(spc_sample):
    """A barbiturate: base peak m/z 119, with tropylium at 91."""
    x, y, _, _ = read_spc_blocks(spc_sample('ms-barb.spc'))[0]
    assert np.isclose(x[int(np.argmax(y))], 119.15, atol=0.01)
    assert y[np.argmin(np.abs(x - 91.0))] > 0


def test_raman_is_recognised_from_the_axis_not_the_experiment_code(spc_sample):
    """
    fexper is 0 in every vendor file, including the Raman one.

    SPC.H explains why: older software needed TCGRAM set before fexper meant
    anything. So the technique comes from fxtype.
    """
    spectrum = read_spectrum(spc_sample('raman.spc'))
    assert spectrum.metadata['spc_experiment'] == 'General'
    assert spectrum.metadata['technique_hint'] == 'Raman'
    assert spectrum.x_quantity == 'Raman shift'


def test_collection_date_is_decoded(spc_sample):
    """The packed bitfield: year 12 bits, month 4, day 5, hour 5, minute 6."""
    spectrum = read_spectrum(spc_sample('ir-nh4.spc'))
    assert spectrum.metadata['collected'].startswith('1987-05-28')


def test_every_vendor_file_is_readable(spc_samples):
    """
    Nothing in the set falls over.

    Coverage across gx-y, x-y, XYXY, custom axis labels, fixed-point y and
    both chromatograms.
    """
    files = sorted(spc_samples.glob('*.spc')) + sorted(spc_samples.glob('*.cgm'))
    assert len(files) == 12
    for path in files:
        collection = read_spectra(path)
        assert len(collection) >= 1
        for spectrum in collection:
            assert len(spectrum.x) == len(spectrum.y) > 0
            assert np.all(np.isfinite(spectrum.y))
            assert np.all(np.isfinite(spectrum.x))


def test_size_predicted_from_the_header_matches_the_file(spc_samples):
    """
    The check that validated the layout in the first place.

    512 + (4*npts if TXVALS) + nsub*(32 + 4*npts) must equal the file length
    exactly for every evenly-spaced file. A single wrong offset anywhere in
    the header puts garbage in fnpts or fnsub and this stops closing.
    """
    checked = 0
    for path in sorted(spc_samples.iterdir()):
        raw = path.read_bytes()
        flags = raw[0]
        if flags & 0x40:          # XYXY is sized from its subheaders instead
            continue
        points, = struct.unpack_from('<i', raw, 4)
        subfiles, = struct.unpack_from('<i', raw, 24)
        width = 2 if flags & 0x01 else 4
        predicted = (512 + (4 * points if flags & 0x80 else 0)
                     + subfiles * (32 + width * points))
        assert predicted == len(raw), f"{path.name}: {predicted} != {len(raw)}"
        checked += 1
    assert checked == 11
