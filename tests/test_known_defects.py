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

import warnings

import numpy as np
import pytest

from spectroscopy.spectra import Spectrum, _infer_file_type


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
# D1 -- FIXED in Phase 0.5: .dpt has its own reader and keeps every point
# --------------------------------------------------------------------------

def test_dpt_keeps_all_points(dpt_file):
    """The defect: reading via 'tsv' consumed the first data row as a header."""
    path, xs = dpt_file
    spec = Spectrum(str(path.parent) + "/", path.name, "dpt")
    assert len(spec.x) == len(xs)
    assert spec.x[0] == xs[0]
    assert spec.x[-1] == xs[-1]


def test_dpt_is_inferred_from_the_extension(dpt_file):
    path, xs = dpt_file
    spec = Spectrum(str(path))
    assert len(spec.x) == len(xs)


def test_the_old_tsv_route_no_longer_loses_a_point_either(dpt_file):
    """
    The delimited-text reader now decides whether row 1 is a header by looking
    at it, rather than assuming one. So the 17 notebooks that say 'tsv' for a
    .dpt file get all their data even before they are updated to say 'dpt'.
    """
    path, xs = dpt_file
    spec = Spectrum(str(path.parent) + "/", path.name, "tsv")
    assert len(spec.x) == len(xs)
    assert spec.x[0] == xs[0]


def test_a_real_header_is_still_honoured(tmp_path):
    """Sniffing must not throw away a genuine header row."""
    path = tmp_path / "with_header.csv"
    path.write_text("Wavelength,Absorbance\n400.0,0.10\n401.0,0.20\n")
    spec = Spectrum(str(path))
    assert len(spec.x) == 2
    assert spec.x_label == "Wavelength"
    assert spec.y_label == "Absorbance"


def test_a_headerless_csv_keeps_every_point(tmp_path):
    """The .csv sibling of D1: data/uvvis_spectra/*.csv have no header row."""
    path = tmp_path / "no_header.csv"
    path.write_text("950,-0.1341\n949.5,-0.14331\n949,-0.15252\n")
    spec = Spectrum(str(path))
    assert len(spec.x) == 3
    assert spec.x[0] == 950.0


# --------------------------------------------------------------------------
# D2 -- FIXED in Phase 0.5: copies no longer alias the source metadata
# --------------------------------------------------------------------------

def test_arithmetic_does_not_alias_metadata():
    """The idiom from FTIR_sugars / Membrane_Analysis / Biofilm_CK."""
    left = _spectrum([1, 2, 3], [1, 1, 1])
    right = _spectrum([1, 2, 3], [2, 2, 2])
    left.set_sample("original")

    average = (left + right) / 2.0
    average.metadata["sample"] = "average"

    assert left.metadata["sample"] == "original"
    assert average.metadata is not left.metadata


def test_explicit_copy_does_not_alias_metadata():
    original = _spectrum([1, 2, 3], [1, 1, 1])
    original.set_sample("original")

    duplicate = Spectrum(original)
    duplicate.metadata["sample"] = "duplicate"

    assert original.metadata["sample"] == "original"


def test_copy_is_deep_so_nested_values_are_independent():
    """
    metadata values can be containers -- the .dpt reader stores the '#' header
    block as a list -- so a shallow copy would still share them.
    """
    original = _spectrum([1, 2, 3], [1, 1, 1])
    original.metadata["file_header"] = ["#\tSample\t\tbuffer"]

    duplicate = Spectrum(original)
    duplicate.metadata["file_header"].append("#\tOperator\tsomeone else")

    assert original.metadata["file_header"] == ["#\tSample\t\tbuffer"]


@pytest.mark.parametrize("operation", [
    lambda a, b: a + b,
    lambda a, b: a - b,
    lambda a, b: a * b,
    lambda a, b: a / b,
    lambda a, b: a + 1.0,
    lambda a, b: a * 2.0,
    lambda a, b: 1.0 - a,
    lambda a, b: -a,
])
def test_every_operator_returns_independent_metadata(operation):
    left = _spectrum([1, 2, 3], [1, 1, 1])
    right = _spectrum([1, 2, 3], [2, 2, 2])
    left.set_sample("original")

    result = operation(left, right)
    result.metadata["sample"] = "changed"

    assert left.metadata["sample"] == "original"


def test_accumulate_in_a_loop_keeps_sources_clean():
    """
    The averaging pattern from Biofilm_CK: `acc = a - a; acc += s/n`, where
    every step went through the copy constructor and so shared one dict.
    """
    sources = []
    for index in range(4):
        spec = _spectrum([1, 2, 3], [index, index, index])
        spec.set_sample(f"sample_{index}")
        sources.append(spec)

    accumulator = sources[0] - sources[0]
    for spec in sources:
        accumulator += spec / len(sources)
    accumulator.metadata["sample"] = "mean"

    assert [s.metadata["sample"] for s in sources] == \
        ["sample_0", "sample_1", "sample_2", "sample_3"]


# --------------------------------------------------------------------------
# D3 -- FIXED in Phase 1: axis compatibility is checked
# --------------------------------------------------------------------------

def test_mismatched_axis_lengths_raise_a_helpful_error():
    """Used to be a bare numpy 'could not be broadcast' with no way forward."""
    short = _spectrum([1, 2, 3], [1, 1, 1])
    long_ = _spectrum([1, 2, 3, 4], [1, 1, 1, 1])
    with pytest.raises(ValueError, match="different lengths"):
        _ = short + long_
    with pytest.raises(ValueError, match="resample"):
        _ = short + long_


def test_differently_positioned_axes_warn():
    """Same length, different x: used to proceed silently and produce junk."""
    here = _spectrum([1000, 1001, 1002], [1, 1, 1])
    there = _spectrum([2000, 2001, 2002], [1, 1, 1])
    with pytest.warns(RuntimeWarning, match="different x positions"):
        result = here + there
    assert np.allclose(result.x, [1000, 1001, 1002])


def test_float_noise_on_the_axis_does_not_warn():
    """
    Spectra off one instrument differ in the last decimal; refusing or nagging
    about those would break real workflows.
    """
    here = _spectrum([1000.0, 1001.0, 1002.0], [1, 1, 1])
    there = _spectrum([1000.0000001, 1001.0000001, 1002.0000001], [1, 1, 1])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _ = here + there


def test_resample_makes_incompatible_spectra_combinable():
    coarse = _spectrum([1000, 1002, 1004], [1, 1, 1])
    fine = _spectrum([1000, 1001, 1002, 1003, 1004], [2, 2, 2, 2, 2])
    combined = coarse + fine.resample(coarse.x)
    assert len(combined) == 3
    assert np.allclose(combined.y, 3.0)


# --------------------------------------------------------------------------
# D4 -- FIXED in Phase 2: .spy 1.0 round-trips identity, units and history
# --------------------------------------------------------------------------

def test_spy_round_trips_name_and_labels(tmp_path):
    """
    The defect: the 0.0 writer wrote name and labels, and the reader never
    parsed them back, so saving and reloading silently anonymised a spectrum.
    """
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    spec.name = "MySpectrum"
    spec.set_type("ATR-FTIR")

    target = tmp_path / "rt.spy"
    spec.save_as(str(target), "spy")
    back = Spectrum("", str(target), "spy")

    assert back.name == "MySpectrum"
    assert back.x_label == spec.x_label
    assert back.technique == "ATR-FTIR"
    assert back.x_unit == "cm^-1"


def test_spy_round_trips_data_and_metadata(tmp_path):
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    spec.name = "MySpectrum"
    spec.set_type("ATR-FTIR")
    spec.set_sample("buffer")

    target = tmp_path / "rt.spy"
    spec.save_as(str(target), "spy")
    back = Spectrum("", str(target), "spy")

    assert np.allclose(back.x, spec.x)
    assert np.allclose(back.y, spec.y)
    assert back.metadata == spec.metadata


# --------------------------------------------------------------------------
# D5 -- format dispatch tables disagree
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,expected", [
    ("experiment.dx", "jcamp"),      # was rejected: FILE_EXTS had a '.DX0' typo
    ("experiment.jdx", "jcamp"),
    ("experiment.jcamp", "jcamp"),
    ("sample.dpt", "dpt"),
    ("sample.DPT", "dpt"),
    ("data.csv", "csv"),
    ("data.CSV", "csv"),
    ("data.tsv", "tsv"),
    ("saved.spy", "spy"),
    ("mystery.qqq", "unknown"),
])
def test_extension_inference(name, expected):
    """
    D5 (fixed): io.jcamp always read .dx fine, but FILE_EXTS listed '.DX0' and
    matched case-sensitively, so Spectrum() refused the file -- which is why
    the phage notebooks bypass Spectrum and call formats.jcamp directly.
    """
    # pylint: disable=protected-access
    assert _infer_file_type(name) == expected


def test_extension_matching_is_case_insensitive(tmp_path):
    """D5 (fixed): .DPT / .CSV used to be rejected."""
    path = tmp_path / "sample.DPT"
    path.write_text("1000.0\t0.1\n1001.0\t0.2\n")
    spec = Spectrum(str(path))
    assert len(spec.x) == 2


@pytest.mark.parametrize("file_type", ["csv", "tsv"])
def test_supported_types_save_correctly(file_type, tmp_path):
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    target = tmp_path / f"out.{file_type}"
    spec.save_as(str(target), file_type)
    assert target.read_text().strip() != ""


def test_saving_an_unhandled_type_raises_without_touching_the_file(tmp_path):
    """
    D5 (fixed): reload() knew 5 types and save() only 4, so an unhandled type
    fell through the match with the file already open for writing and left a
    0-byte file. The type is now validated *before* open(..., 'w').
    """
    spec = _spectrum([1, 2, 3], [4, 5, 6])
    target = tmp_path / "precious.xyz"
    target.write_text("existing data that must not be destroyed\n")

    with pytest.raises(ValueError, match="Cannot write"):
        spec.save_as(str(target), "xyz")

    assert target.read_text().startswith("existing data")
