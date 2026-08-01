# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Executable record of the defects found in the Phase 0 review.

These tests assert the *current, wrong* behaviour and are marked xfail(strict)
where a fix is planned. When Phase 0.5 fixes a defect the corresponding test
turns XPASS and fails the build -- that is the signal to invert the assertion
and drop the marker. Nothing here should be read as endorsing the behaviour.

See SpectroscoPy_Codebase_Review.md section 3.2 for D1-D7.
"""

import numpy as np
import pytest

from spectroscopy.spectra import Spectrum


@pytest.fixture
def dpt_file(tmp_path):
    """A Bruker-style .dpt export: tab separated, no header row."""
    path = tmp_path / "sample.dpt"
    xs = [1000.0, 1001.0, 1002.0, 1003.0]
    path.write_text("".join(f"{x}\t{x / 1000:.4f}\n" for x in xs))
    return path, xs


def _spectrum(xs, ys):
    spec = Spectrum()
    spec.x = np.asarray(xs, dtype=float)
    spec.y = np.asarray(ys, dtype=float)
    return spec


# --------------------------------------------------------------------------
# D1 -- loading a .dpt as 'tsv' eats the first data point and the axis labels
# --------------------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason="D1: fixed in Phase 0.5 by a real 'dpt' reader")
def test_dpt_keeps_all_points(dpt_file):
    path, xs = dpt_file
    spec = Spectrum(str(path.parent) + "/", path.name, "tsv")
    assert len(spec.x) == len(xs)
    assert spec.x[0] == xs[0]


def test_dpt_currently_loses_first_point(dpt_file):
    """Pin the damage precisely, so the Phase 0.5 fix can be verified."""
    path, xs = dpt_file
    spec = Spectrum(str(path.parent) + "/", path.name, "tsv")
    assert len(spec.x) == len(xs) - 1
    assert spec.x[0] == xs[1]
    # ... and the first data row was consumed as a header:
    assert spec.x_label == "1000.0"
    assert spec.y_label == "1.0000"


# --------------------------------------------------------------------------
# D2 -- arithmetic results alias the left operand's metadata dict
# --------------------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason="D2: fixed in Phase 0.5 by copying metadata")
def test_arithmetic_does_not_alias_metadata():
    left = _spectrum([1, 2, 3], [1, 1, 1])
    right = _spectrum([1, 2, 3], [2, 2, 2])
    left.set_sample("original")

    average = (left + right) / 2.0
    average.metadata["sample"] = "average"

    assert left.metadata["sample"] == "original"


def test_arithmetic_currently_aliases_metadata():
    left = _spectrum([1, 2, 3], [1, 1, 1])
    right = _spectrum([1, 2, 3], [2, 2, 2])
    left.set_sample("original")

    average = (left + right) / 2.0
    assert average.metadata is left.metadata
    average.metadata["sample"] = "average"
    assert left.metadata["sample"] == "average"      # <- the bug


# --------------------------------------------------------------------------
# D3 -- no axis compatibility checking
# --------------------------------------------------------------------------

def test_mismatched_axis_lengths_raise_a_numpy_error():
    short = _spectrum([1, 2, 3], [1, 1, 1])
    long_ = _spectrum([1, 2, 3, 4], [1, 1, 1, 1])
    with pytest.raises(ValueError, match="broadcast"):
        _ = short + long_


def test_differently_positioned_axes_are_silently_wrong():
    """Same length, different x -- arithmetic proceeds and the result is junk."""
    here = _spectrum([1000, 1001, 1002], [1, 1, 1])
    there = _spectrum([2000, 2001, 2002], [1, 1, 1])
    result = here + there
    assert np.allclose(result.y, 2.0)              # no complaint at all
    assert np.allclose(result.x, [1000, 1001, 1002])


# --------------------------------------------------------------------------
# D4 -- .spy does not round-trip name or axis labels
# --------------------------------------------------------------------------

@pytest.mark.xfail(strict=True, reason="D4: fixed in Phase 2 when .spy carries history")
def test_spy_round_trips_name_and_labels(tmp_path):
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    spec.name = "MySpectrum"
    spec.set_type("ATR-FTIR")

    target = tmp_path / "rt.spy"
    spec.save_as(str(target), "spy")
    back = Spectrum("", str(target), "spy")

    assert back.name == "MySpectrum"
    assert back.x_label == spec.x_label


def test_spy_currently_round_trips_data_but_not_identity(tmp_path):
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    spec.name = "MySpectrum"
    spec.set_type("ATR-FTIR")

    target = tmp_path / "rt.spy"
    spec.save_as(str(target), "spy")
    back = Spectrum("", str(target), "spy")

    # data and metadata survive ...
    assert np.allclose(back.x, spec.x)
    assert np.allclose(back.y, spec.y)
    assert back.metadata == spec.metadata
    # ... identity does not
    assert back.name != "MySpectrum"
    assert back.x_label == "Wavelength (nm)"        # the empty-ctor default


# --------------------------------------------------------------------------
# D5 -- format dispatch tables disagree
# --------------------------------------------------------------------------

def test_jcamp_extension_is_not_recognised(tmp_path):
    """
    spectroscopy.io.jcamp reads .dx perfectly well, but Spectrum() refuses the
    extension -- which is why the phage notebooks bypass Spectrum entirely.
    """
    path = tmp_path / "experiment.dx"
    path.write_text("##TITLE=nothing\n##END=\n")
    with pytest.raises(TypeError, match="Unknown filetype"):
        Spectrum(str(path))


def test_dpt_filetype_is_rejected_despite_having_a_reader_branch(tmp_path):
    """reload() has a `case 'dpt'` branch the constructor can never reach."""
    path = tmp_path / "sample.dpt"
    path.write_text("1000.0\t0.1\n")
    with pytest.raises(TypeError, match="Unknown filetype"):
        Spectrum(str(path.parent) + "/", path.name, "dpt")


def test_extension_matching_is_case_sensitive(tmp_path):
    path = tmp_path / "sample.CSV"
    path.write_text("x,y\n1.0,2.0\n")
    with pytest.raises(TypeError, match="Unknown filetype"):
        Spectrum(str(path))


@pytest.mark.parametrize("file_type", ["csv", "tsv"])
def test_supported_types_save_correctly(file_type, tmp_path):
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    target = tmp_path / f"out.{file_type}"
    spec.save_as(str(target), file_type)
    assert target.read_text().strip() != ""


@pytest.mark.xfail(strict=True,
                   reason="D5: save() must raise on an unhandled type, not truncate")
def test_saving_an_unhandled_type_raises_rather_than_truncating(tmp_path):
    """
    reload() knows 5 file types, save() only 4. An unhandled type falls through
    the match with the file already open for writing, silently producing an
    empty file. Latent today ('dpt' is unreachable via the constructor); a live
    data-loss bug as soon as Phase 0.5 adds a real 'dpt' type.
    """
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    target = tmp_path / "out.dpt"
    with pytest.raises((ValueError, TypeError)):
        spec.save_as(str(target), "dpt")


def test_unhandled_save_type_currently_truncates_silently(tmp_path):
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    target = tmp_path / "out.dpt"
    spec.save_as(str(target), "dpt")          # no error ...
    assert target.stat().st_size == 0          # ... and nothing written
