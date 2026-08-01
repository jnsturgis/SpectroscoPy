# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The native .spy format.

.spy exists to be the lossless one (TODO.txt item 1). Version 0.0 was not: it
wrote name and axis labels that its own reader never read back. Version 1.0
carries identity, units, metadata and -- the point of it -- the processing
history, which is what makes provenance survive leaving memory.
"""

import json

import numpy as np
import pytest

from spectroscopy.io import spy
from spectroscopy.spectra import Spectrum


@pytest.fixture
def processed():
    """A spectrum that has actually been through a pipeline."""
    spec = Spectrum()
    spec.x = np.linspace(900.0, 1800.0, 91)
    spec.y = np.exp(-((spec.x - 1650.0) ** 2) / (2 * 40.0 ** 2)) + 0.1
    spec.set_type("ATR-FTIR")
    spec.set_sample("Glucose")
    spec.name = "Glucose average"

    reference = Spectrum(spec)
    reference.set_sample("H2O")

    return (spec
            .subtract_reference(reference, factor=0.7)
            .crop(1000, 1750)
            .baseline_correct('rubberband')
            .normalize('max', window=(1050, 1080)))


def _round_trip(spectrum, tmp_path, name="rt.spy"):
    target = tmp_path / name
    spectrum.save_as(str(target), "spy")
    return Spectrum("", str(target), "spy"), target


# ---------------------------------------------------------------------------
# 1.0 round trip
# ---------------------------------------------------------------------------

def test_data_survives(processed, tmp_path):
    back, _ = _round_trip(processed, tmp_path)
    assert np.allclose(back.x, processed.x)
    assert np.allclose(back.y, processed.y, atol=1e-9)


def test_identity_and_units_survive(processed, tmp_path):
    back, _ = _round_trip(processed, tmp_path)
    assert back.name == processed.name
    assert back.technique == processed.technique
    assert back.x_unit == processed.x_unit
    assert back.y_unit == processed.y_unit
    assert back.x_quantity == processed.x_quantity
    assert back.x_label == processed.x_label


def test_metadata_survives(processed, tmp_path):
    back, _ = _round_trip(processed, tmp_path)
    assert back.metadata == processed.metadata


def test_history_survives(processed, tmp_path):
    """The whole reason for a native format."""
    back, _ = _round_trip(processed, tmp_path)

    assert [s.name for s in back.history] == [s.name for s in processed.history]
    assert [s.params for s in back.history] == [s.params for s in processed.history]
    assert [s.timestamp for s in back.history] == \
        [s.timestamp for s in processed.history]


def test_the_hand_tuned_factor_survives(processed, tmp_path):
    """
    The scientific crux: the 0.7 water factor chosen by eye must come back out
    of the file, not just the numbers it produced.
    """
    back, _ = _round_trip(processed, tmp_path)
    step = next(s for s in back.history if s.name == 'subtract_reference')
    assert step.params['factor'] == 0.7
    assert step.params['reference']['sample'] == 'H2O'


def test_a_second_round_trip_changes_nothing(processed, tmp_path):
    once, _ = _round_trip(processed, tmp_path, "one.spy")
    twice, _ = _round_trip(once, tmp_path, "two.spy")
    assert np.allclose(twice.y, once.y)
    assert [s.to_dict() for s in twice.history] == [s.to_dict() for s in once.history]


def test_the_file_is_readable_by_eye(processed, tmp_path):
    """Data stays tab separated so ordinary tools still work on it."""
    _, target = _round_trip(processed, tmp_path)
    lines = target.read_text().splitlines()

    assert lines[0].startswith("# spy format 1.0")
    data_start = lines.index("# data")
    numbers = lines[data_start + 2].split("\t")
    assert len(numbers) == 2
    float(numbers[0]), float(numbers[1])

    header = json.loads("\n".join(lines[2:data_start]))
    assert header['name'] == processed.name
    assert len(header['history']) == len(processed.history)


def test_a_label_override_survives(tmp_path):
    """Labels read from a CSV header are an override, not a technique default."""
    spec = Spectrum()
    spec.x, spec.y = np.array([1.0, 2.0]), np.array([3.0, 4.0])
    spec.x_label = "Whatever the file said"
    back, _ = _round_trip(spec, tmp_path)
    assert back.x_label == "Whatever the file said"


# ---------------------------------------------------------------------------
# legacy 0.0
# ---------------------------------------------------------------------------

def test_legacy_files_still_load(tmp_path):
    legacy = tmp_path / "old.spy"
    legacy.write_text(
        "# spy format 0.0 file\n"
        "OldSpectrum\n"
        "Wavenumber\tAbsorbance\n"
        "# x,y data\n"
        "1000.000\t0.10000\n"
        "1001.000\t0.20000\n"
        "# metadata\n"
        "{'sample': 'legacy', 'spec_type': 'ATR-FTIR'}\n"
    )
    spec = Spectrum("", str(legacy), "spy")

    assert np.allclose(spec.x, [1000.0, 1001.0])
    assert np.allclose(spec.y, [0.1, 0.2])
    assert spec.metadata['sample'] == 'legacy'
    # ... and the name and labels the 0.0 reader used to drop:
    assert spec.name == "OldSpectrum"
    assert spec.x_label == "Wavenumber"


def test_version_is_taken_from_the_file_not_the_caller(tmp_path):
    """
    reload() used to pass format='0.0' unconditionally. Trusting the caller
    over the file is how a 1.0 file gets parsed as 0.0.
    """
    spec = Spectrum()
    spec.x, spec.y = np.array([1.0, 2.0]), np.array([3.0, 4.0])
    spec.name = "Modern"
    target = tmp_path / "modern.spy"
    spec.save_as(str(target), "spy")

    restored = Spectrum()
    with open(target, encoding="utf-8") as handle:
        spy.read(handle, restored, format='0.0')     # deliberately wrong
    assert restored.name == "Modern"


def test_an_unknown_version_is_refused(tmp_path):
    path = tmp_path / "future.spy"
    path.write_text("# spy format 9.9\n# header\n{}\n# data\n")
    with pytest.raises(ValueError, match="unsupported .spy format"):
        Spectrum("", str(path), "spy")


def test_an_empty_file_is_refused(tmp_path):
    path = tmp_path / "empty.spy"
    path.write_text("")
    with pytest.raises(ValueError, match="empty"):
        Spectrum("", str(path), "spy")
