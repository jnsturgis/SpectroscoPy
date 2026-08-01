# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The JCAMP-DX reader.

Before Phase 2 this read **zero** of the 79 JCAMP files in data/. The reader
had been unreachable (FILE_EXTS listed a '.DX0' typo, fixed as part of D5), so
nothing exercised it and four bugs introduced when the code was adapted from
github.com/nzhagen/jcamp went unnoticed. These tests pin each of them.
"""

import glob
import os

import numpy as np
import pytest

from spectroscopy.io import jcamp
from spectroscopy.spectra import Spectrum

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
JCAMP_FILES = sorted(glob.glob(os.path.join(DATA, "**", "*.jdx"), recursive=True))


def _load(path):
    return Spectrum(os.path.dirname(path) + os.sep, os.path.basename(path))


@pytest.mark.skipif(not JCAMP_FILES, reason="sample data not present")
def test_every_sample_file_loads():
    """The headline: 0 of 79 used to load."""
    failures = []
    for path in JCAMP_FILES:
        try:
            _load(path)
        except Exception as error:                    # noqa: BLE001
            failures.append(f"{os.path.basename(path)}: {type(error).__name__}: {error}")
    assert not failures, "\n".join(failures)


@pytest.mark.skipif(not JCAMP_FILES, reason="sample data not present")
def test_most_files_carry_data():
    """
    A handful are compound (LINK) files whose parent legitimately holds no data
    of its own -- their spectra live in metadata['Children'].
    """
    with_data = sum(1 for path in JCAMP_FILES if len(_load(path)) > 1)
    assert with_data >= len(JCAMP_FILES) - 6


def test_the_xydata_header_is_not_parsed_as_data(tmp_path):
    """
    Bug 1: a `continue` was lost when the code was adapted, so the
    '##XYDATA=(X++(Y..Y))' line itself was fed to the data parser. Its 'X' and
    'Y' are DUP digits, so it failed inside the compression expansion -- which
    is what made every single file unreadable.
    """
    path = tmp_path / "tiny.jdx"
    path.write_text(
        "##TITLE=tiny\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n"
        "##XUNITS=1/CM\n##YUNITS=TRANSMITTANCE\n"
        "##FIRSTX=1000\n##LASTX=1003\n##NPOINTS=4\n"
        "##XFACTOR=1\n##YFACTOR=1\n"
        "##XYDATA=(X++(Y..Y))\n"
        "1000 0.1 0.2 0.3 0.4\n"
        "##END=\n"
    )
    spec = _load(path)
    assert len(spec) == 4
    assert np.allclose(spec.y, [0.1, 0.2, 0.3, 0.4])
    assert np.isclose(spec.x[0], 1000.0)


def test_asdf_detection_happens_on_the_first_data_line(tmp_path):
    """
    Bug 2: the test was inverted (`if len(y) > 0` instead of `if not len(y)`),
    so the flag was read before it was ever set -- UnboundLocalError on every
    file that got past bug 1.
    """
    path = tmp_path / "affn.jdx"
    path.write_text(
        "##TITLE=affn\n##JCAMP-DX=4.24\n##XUNITS=1/CM\n##YUNITS=ABSORBANCE\n"
        "##FIRSTX=1000\n##LASTX=1002\n##NPOINTS=3\n##XFACTOR=1\n##YFACTOR=1\n"
        "##XYDATA=(X++(Y..Y))\n1000 0.1 0.2 0.3\n##END=\n"
    )
    assert len(_load(path)) == 3


def test_dup_expansion_does_not_run_off_the_start_of_a_line():
    """
    Bug 3: expanding DUP compression walked backwards for a DIF digit and, when
    the previous value began with a SQZ digit instead, ran past index 0 into
    negative indices and raised IndexError.
    """
    # 'A' is SQZ +1, 'T' is DUP x2: the value 1, then a duplicate of it. The
    # previous value starts with a SQZ digit, which is the case upstream's
    # DIF-only search walked straight past.
    assert jcamp.jcamp_parse("AT") == [1.0, 1.0]
    assert jcamp.jcamp_parse("A1T") == [11.0, 11.0]     # SQZ +1 then digit 1
    # DIF-coded values keep working: +1, then +1 difference, then duplicate.
    assert jcamp.jcamp_parse("AJT") == [1.0, 2.0, 3.0]
    # A DUP as the very first character has nothing before it; must not crash.
    assert jcamp.jcamp_parse("T") == []


def test_units_are_not_all_called_wavenumber():
    """
    Bug 4: the x label was hardcoded to 'Wavenumber', mislabelling every
    UV/Vis, fluorescence and NMR file.
    """
    uvvis = os.path.join(DATA, "uvvis_spectra", "toluene.jdx")
    if os.path.exists(uvvis):
        assert "Wavelength" in _load(uvvis).x_label

    infrared = os.path.join(DATA, "infrared_spectra", "ethanol.jdx")
    if os.path.exists(infrared):
        assert "Wavenumber" in _load(infrared).x_label


def test_a_non_utf8_file_still_loads():
    """
    Older instruments write latin-1 bytes into comment fields. Refusing a file
    whose numbers are perfectly readable is not helpful.
    """
    path = os.path.join(DATA, "neutron_scattering_spectra", "emodine.jdx")
    if not os.path.exists(path):
        pytest.skip("sample file not present")
    assert len(_load(path)) > 1


def test_compound_files_expose_their_children():
    path = os.path.join(DATA, "infrared_spectra", "example_compound_file.jdx")
    if not os.path.exists(path):
        pytest.skip("sample file not present")
    spec = _load(path)
    children = spec.metadata.get('Children', [])
    assert len(children) > 1
    assert all(len(child) > 1 for child in children)
