# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Tests for the .dpt reader.

The cases here are taken from a survey of the 825 real .dpt files under
~/Documents/Research: 681 tab separated, 140 comma separated, 2 carrying a '#'
metadata header, 1 binary, and almost all of them CRLF.
"""

import numpy as np
import pytest

from spectroscopy.io import dpt
from spectroscopy.spectra import Spectrum

XS = [3997.63597, 3995.70754, 3993.77911]
YS = [-0.0096921418, -0.0096850535, -0.0096661681]


def _write(tmp_path, text, name="sample.dpt", newline=""):
    path = tmp_path / name
    path.write_text(text, newline=newline)
    return path


def _load(path):
    return Spectrum(str(path.parent) + "/", path.name, "dpt")


# --------------------------------------------------------------------------
# Separator variants -- the reason this reader exists
# --------------------------------------------------------------------------

@pytest.mark.parametrize("delimiter,label", [
    ("\t", "tab (681 of 825 real files)"),
    (",", "comma (140 of 825 -- OPUS follows the machine locale)"),
    (";", "semicolon"),
    ("   ", "runs of spaces"),
])
def test_separator_is_sniffed(tmp_path, delimiter, label):
    text = "".join(f"{x}{delimiter}{y}\r\n" for x, y in zip(XS, YS))
    spec = _load(_write(tmp_path, text))
    assert np.allclose(spec.x, XS), label
    assert np.allclose(spec.y, YS), label


def test_comma_separated_dpt_used_to_fail_outright(tmp_path):
    """
    Candice's whole FTIR_spectra directory is comma separated; reading it as
    'tsv' raised ValueError, which is why those notebooks used np.genfromtxt.
    """
    path = _write(tmp_path, "".join(f"{x},{y}\r\n" for x, y in zip(XS, YS)))
    with pytest.raises(ValueError):
        Spectrum(str(path.parent) + "/", path.name, "tsv")
    assert len(_load(path).x) == 3


def test_explicit_delimiter_overrides_sniffing(tmp_path):
    path = _write(tmp_path, "1.0,2.0\n3.0,4.0\n")
    spec = Spectrum()
    with open(path, encoding="utf-8") as handle:
        dpt.read(handle, spec, delimiter=",")
    assert np.allclose(spec.x, [1.0, 3.0])


# --------------------------------------------------------------------------
# Line endings and stray whitespace
# --------------------------------------------------------------------------

@pytest.mark.parametrize("ending", ["\r\n", "\n", "\r"])
def test_line_endings(tmp_path, ending):
    text = ending.join(f"{x}\t{y}" for x, y in zip(XS, YS)) + ending
    spec = _load(_write(tmp_path, text, newline=""))
    assert len(spec.x) == 3
    assert np.allclose(spec.y, YS)


def test_blank_lines_are_skipped_not_counted(tmp_path):
    text = "\n\n1000.0\t0.1\n\n1001.0\t0.2\n\n"
    spec = _load(_write(tmp_path, text))
    assert np.allclose(spec.x, [1000.0, 1001.0])


# --------------------------------------------------------------------------
# The '#' metadata header on James' own reference spectra
# --------------------------------------------------------------------------

def test_comment_header_is_preserved_verbatim(tmp_path):
    text = (
        "#\r\n"
        "#\tSample\t\t50mM Phosphate pH7.5\r\n"
        "#\tReference\tWater\r\n"
        "#\tOperator\tJames\r\n"
        "# Wavenumber\tAbsorption\r\n"
        "4000.00000\t-0.00120\r\n"
        "3999.03570\t-0.00128\r\n"
    )
    spec = _load(_write(tmp_path, text))

    assert len(spec.x) == 2                      # header rows are not data
    assert spec.x[0] == 4000.0
    header = spec.metadata["file_header"]
    assert any("50mM Phosphate pH7.5" in line for line in header)
    assert any("Operator" in line for line in header)


def test_no_header_means_no_metadata_key(tmp_path):
    spec = _load(_write(tmp_path, "1.0\t2.0\n"))
    assert "file_header" not in spec.metadata


# --------------------------------------------------------------------------
# Failure modes -- clear errors, not garbage
# --------------------------------------------------------------------------

def test_binary_file_raises_a_useful_error(tmp_path):
    """One real file, eau2.0.dpt, is binary despite the extension."""
    path = tmp_path / "binary.dpt"
    path.write_bytes(b"\xc3\xbe\xc3\xbe\x00\x00\x00\x00\\\x18,A\x18\x00\x00\x00")
    with pytest.raises((ValueError, UnicodeDecodeError)):
        _load(path)


def test_empty_file_raises(tmp_path):
    with pytest.raises(ValueError, match="no data points"):
        _load(_write(tmp_path, "\n\n"))


def test_ragged_line_reports_its_number(tmp_path):
    text = "1000.0\t0.1\n1001.0\t0.2\nGARBAGE HERE\n"
    with pytest.raises(ValueError, match="line 3"):
        _load(_write(tmp_path, text))


# --------------------------------------------------------------------------
# Round trip
# --------------------------------------------------------------------------

def test_round_trip(tmp_path):
    original = _load(_write(tmp_path, "".join(f"{x}\t{y}\r\n"
                                              for x, y in zip(XS, YS))))
    target = tmp_path / "out.dpt"
    original.save_as(str(target), "dpt")

    back = _load(target)
    assert np.allclose(back.x, original.x)
    assert np.allclose(back.y, original.y, atol=1e-9)
