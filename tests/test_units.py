# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Unit conversion.

Roadmap section 2.2: "spectrum.to('nm') should be a first-class, tested
operation", because nm-vs-cm^-1 and absorbance-vs-%T confusion is exactly the
kind of bug that erodes trust in a shared tool.
"""

import numpy as np
import pytest

from spectroscopy import units
from spectroscopy.spectra import Spectrum


def _uvvis():
    spec = Spectrum()
    spec.x = np.linspace(400.0, 700.0, 31)
    spec.y = np.linspace(0.1, 0.7, 31)
    spec.set_type("UV-Vis")
    spec.set_sample("sample")
    return spec


# ---------------------------------------------------------------------------
# x axis
# ---------------------------------------------------------------------------

def test_known_conversions_against_hand_values():
    """500 nm is 20000 cm^-1 and about 2.4797 eV."""
    converted, _ = units.convert_x(np.array([500.0]), 'nm', 'cm^-1')
    assert np.isclose(converted[0], 20000.0)

    converted, _ = units.convert_x(np.array([500.0]), 'nm', 'eV')
    assert np.isclose(converted[0], 2.4797, atol=1e-4)

    converted, _ = units.convert_x(np.array([1650.0]), 'cm^-1', 'nm')
    assert np.isclose(converted[0], 6060.6, atol=0.1)


@pytest.mark.parametrize("target", ["cm^-1", "eV", "um"])
def test_round_trip_is_exact(target):
    original = np.linspace(400.0, 700.0, 31)
    there, _ = units.convert_x(original, 'nm', target)
    back, _ = units.convert_x(there, target, 'nm')
    assert np.allclose(back, original)


def test_reciprocal_conversion_reverses_the_axis_and_reorders_y():
    """
    nm ascending becomes cm^-1 descending. If y were not reordered with it,
    every later interpolation, hull and crop would be working on scrambled
    data -- the quiet kind of wrong.
    """
    spec = _uvvis()
    converted = spec.to('cm^-1')

    assert np.all(np.diff(converted.x) > 0), "axis should be re-sorted ascending"
    # The point that was at 400 nm is now the largest wavenumber, and must
    # still carry the y value it had at 400 nm.
    assert np.isclose(converted.x[-1], 1e7 / 400.0)
    assert np.isclose(converted.y[-1], spec.y[0])
    assert np.isclose(converted.y[0], spec.y[-1])


def test_conversion_relabels_the_quantity():
    converted = _uvvis().to('cm^-1')
    assert converted.x_unit == 'cm^-1'
    assert converted.x_quantity == 'Wavenumber'      # not 'Wavelength (cm^-1)'
    assert 'cm$^{-1}$' in converted.x_label


def test_a_zero_on_the_axis_is_refused():
    spec = Spectrum()
    spec.x = np.array([0.0, 1.0, 2.0])
    spec.y = np.array([1.0, 1.0, 1.0])
    with pytest.raises(ValueError, match="reciprocal"):
        spec.to('cm^-1')


def test_unknown_unit_is_refused():
    with pytest.raises(ValueError, match="cannot convert x unit"):
        units.convert_x(np.array([1.0]), 'nm', 'furlongs')


# ---------------------------------------------------------------------------
# y axis
# ---------------------------------------------------------------------------

def test_absorbance_to_transmittance_and_back():
    absorbance = np.array([0.0, 1.0, 2.0])
    transmittance = units.convert_y(absorbance, 'absorbance', 'transmittance')
    assert np.allclose(transmittance, [1.0, 0.1, 0.01])
    assert np.allclose(
        units.convert_y(transmittance, 'transmittance', 'absorbance'), absorbance)


def test_absorbance_to_percent_transmittance():
    percent = units.convert_y(np.array([0.0, 1.0, 2.0]), 'absorbance', '%T')
    assert np.allclose(percent, [100.0, 10.0, 1.0])


def test_non_positive_transmittance_cannot_become_absorbance():
    """log10 of zero is not a number the user wants silently in their data."""
    with pytest.raises(ValueError, match="non-positive transmittance"):
        units.convert_y(np.array([0.5, 0.0, -0.1]), 'transmittance', 'absorbance')


def test_arbitrary_units_are_not_convertible():
    """
    Fluorescence counts have no defined conversion; saying so beats inventing
    a factor.
    """
    with pytest.raises(ValueError, match="instrument-specific"):
        units.convert_y(np.array([1.0]), 'a.u.', 'absorbance')


# ---------------------------------------------------------------------------
# the Spectrum method
# ---------------------------------------------------------------------------

def test_to_records_history():
    step = _uvvis().to('cm^-1', '%T').history[-1]
    assert step.name == 'to'
    assert step.params == {'x_unit': 'cm^-1', 'from_x_unit': 'nm',
                           'y_unit': '%T', 'from_y_unit': 'absorbance'}


def test_to_leaves_the_original_alone():
    spec = _uvvis()
    before = spec.x.copy()
    spec.to('cm^-1')
    assert np.allclose(spec.x, before)
    assert spec.x_unit == 'nm'


def test_converting_to_the_same_unit_is_a_no_op():
    spec = _uvvis()
    same = spec.to('nm')
    assert np.allclose(same.x, spec.x)
    assert same.history == spec.history          # nothing worth recording


def test_to_needs_something_to_do():
    with pytest.raises(ValueError, match="needs an x_unit"):
        _uvvis().to()


def test_conversion_survives_the_rest_of_the_pipeline():
    """A converted spectrum must still crop, baseline and peak-pick normally."""
    result = (_uvvis().to('cm^-1')
              .crop(15000, 24000)
              .baseline_correct('rubberband'))
    assert len(result) > 0
    assert np.all(np.diff(result.x) > 0)
    assert [s.name for s in result.history] == ['to', 'crop', 'baseline_correct']
