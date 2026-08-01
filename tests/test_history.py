# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Processing provenance.

The requirement (roadmap section 2.3, settled) is that history records a
structured (name, params) pair rather than a sentence, so a chain built
interactively can later be replayed as a Pipeline. These tests pin that
property rather than the exact wording of anything.
"""

import json

import numpy as np
import pytest

from spectroscopy.history import ProcessingStep
from spectroscopy.spectra import Spectrum


@pytest.fixture
def spectrum():
    spec = Spectrum()
    spec.x = np.linspace(900.0, 1800.0, 451)
    spec.y = np.exp(-((spec.x - 1650.0) ** 2) / (2 * 25.0 ** 2)) + 0.1
    spec.set_type("ATR-FTIR")
    spec.set_sample("test")
    return spec


def test_a_fresh_spectrum_has_no_history(spectrum):
    assert spectrum.history == []


def test_operations_accumulate_in_order(spectrum):
    result = (spectrum
              .crop(1500, 1750)
              .baseline_correct('rubberband')
              .normalize('max')
              .smooth('savgol', window_length=11, polyorder=3))
    assert [step.name for step in result.history] == \
        ['crop', 'baseline_correct', 'normalize', 'smooth']


def test_the_original_is_left_alone(spectrum):
    before = len(spectrum)
    cropped = spectrum.crop(1500, 1750)
    assert len(spectrum) == before
    assert spectrum.history == []
    assert len(cropped) < before


def test_params_are_structured_not_prose(spectrum):
    """The distinction that makes Pipeline.from_history possible."""
    step = spectrum.smooth('SG', [2, 15]).history[-1]
    assert step.name == 'smooth'
    assert step.params == {'method': 'savgol', 'window_length': 15, 'polyorder': 2}


def test_legacy_positional_smooth_is_normalised(spectrum):
    """
    smooth('SG', [2, 15]) means polyorder=2, window=15 -- the reverse of
    scipy's order. It is recorded by keyword so the ambiguity cannot survive
    into history.
    """
    legacy = spectrum.smooth('SG', [2, 15])
    modern = spectrum.smooth('savgol', window_length=15, polyorder=2)
    assert np.allclose(legacy.y, modern.y)
    assert legacy.history[-1].params == modern.history[-1].params


def test_arithmetic_is_recorded(spectrum):
    """
    The scaled reference subtraction is the step that most needs recording --
    the hand-tuned factor is the scientific choice.
    """
    other = Spectrum(spectrum)
    other.set_sample("buffer")
    result = spectrum - 0.995 * other

    names = [step.name for step in result.history]
    assert 'subtract' in names
    step = result.history[-1]
    assert step.params['operand']['kind'] == 'spectrum'
    assert step.params['operand']['sample'] == 'buffer'


def test_scalar_operand_records_its_value(spectrum):
    step = (spectrum * 0.75).history[-1]
    assert step.name == 'multiply'
    assert step.params['operand'] == {'kind': 'scalar', 'value': 0.75}


def test_subtract_reference_records_the_factor(spectrum):
    reference = Spectrum(spectrum)
    reference.set_sample("water")
    result = spectrum.subtract_reference(reference, factor=0.93)

    step = result.history[-1]
    assert step.name == 'subtract_reference'
    assert step.params['factor'] == 0.93
    assert step.params['reference']['sample'] == 'water'


def test_baseline_correct_records_the_method_and_parameters(spectrum):
    step = spectrum.baseline_correct('als', lam=1e5, p=0.01).history[-1]
    assert step.name == 'baseline_correct'
    assert step.params['method'] == 'als'
    assert step.params['lam'] == 1e5
    assert step.params['p'] == 0.01


def test_guide_points_survive_into_history(spectrum):
    guides = [900, 1000, 1100, 1750, 1800]
    step = spectrum.baseline_correct('poly', degree=2, points=guides).history[-1]
    assert step.params['method'] == 'poly'
    assert step.params['points'] == guides


def test_history_is_json_serialisable(spectrum):
    """It has to survive a round trip through .spy or NetCDF."""
    result = (spectrum.crop(1500, 1750)
                      .baseline_correct('rubberband')
                      .normalize('max', window=(1600, 1700)))
    encoded = json.dumps([step.to_dict() for step in result.history])
    restored = [ProcessingStep.from_dict(d) for d in json.loads(encoded)]

    assert [s.name for s in restored] == [s.name for s in result.history]
    assert [s.params for s in restored] == [s.params for s in result.history]
    assert restored[0].timestamp == result.history[0].timestamp


def test_history_does_not_leak_backwards(spectrum):
    """Two branches from one spectrum must not see each other's steps."""
    base = spectrum.crop(1500, 1750)
    left = base.normalize('max')
    right = base.smooth('savgol')

    assert len(base.history) == 1
    assert [s.name for s in left.history] == ['crop', 'normalize']
    assert [s.name for s in right.history] == ['crop', 'smooth']


def test_steps_are_immutable():
    step = ProcessingStep("crop", {"x_min": 900})
    with pytest.raises(AttributeError):
        step.name = "something else"


def test_describe_history_reads_as_a_recipe(spectrum):
    text = spectrum.crop(1500, 1750).normalize('max').describe_history()
    assert text.startswith("1. crop(")
    assert "2. normalize(" in text


def test_baseline_does_not_claim_to_be_a_correction(spectrum):
    """
    baseline() builds a baseline; the correction is baseline_correct(). Doing
    `s - s.baseline()` by hand records an anonymous subtraction, which is why
    the correcting form exists.
    """
    assert spectrum.baseline('rubberband').history[-1].name == 'baseline'
    assert spectrum.baseline_correct('rubberband').history[-1].name == \
        'baseline_correct'

    by_hand = spectrum - spectrum.baseline('rubberband')
    assert by_hand.history[-1].name == 'subtract'
    assert np.allclose(by_hand.y, spectrum.baseline_correct('rubberband').y)
