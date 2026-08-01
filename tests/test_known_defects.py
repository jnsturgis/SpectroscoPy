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


def test_dpt_reading_via_tsv_still_loses_a_point(dpt_file):
    """
    The old route is left working but is still wrong by construction: 'tsv'
    means "delimited text with a header row". Kept as a test so the difference
    between the two paths stays visible.
    """
    path, xs = dpt_file
    spec = Spectrum(str(path.parent) + "/", path.name, "tsv")
    assert len(spec.x) == len(xs) - 1


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
