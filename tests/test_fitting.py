# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Band decomposition, checked against synthetic ground truth.

A fitter is one of the few things that can be tested properly without real
data: build a spectrum from components whose positions, widths and areas are
known exactly, and see whether they come back. Real protein spectra check
something different -- that the recipe around the fitter is right -- and are
still wanted (roadmap section 15.2).
"""

import numpy as np
import pytest

from spectroscopy.fitting import FitResult, fit_components
from spectroscopy.lineshapes import gauss, lorentz, spec_comp
from spectroscopy.spectra import Spectrum

#: (position, fwhm, amplitude) -- roughly an amide I envelope.
TRUTH = [(1628.0, 14.0, 0.55), (1645.0, 16.0, 0.30),
         (1655.0, 12.0, 0.90), (1672.0, 14.0, 0.35)]


@pytest.fixture
def amide_i():
    x = np.linspace(1600.0, 1700.0, 401)
    y = sum(gauss(x, position, width, height) for position, width, height in TRUTH)
    return x, y


# ---------------------------------------------------------------------------
# recovering what was put in
# ---------------------------------------------------------------------------

def test_it_recovers_known_components(amide_i):
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)

    np.testing.assert_allclose(fit.position, [p for p, _, _ in TRUTH], atol=1e-3)
    np.testing.assert_allclose(fit.fwhm, [w for _, w, _ in TRUTH], atol=1e-3)
    np.testing.assert_allclose(fit.amplitude, [a for _, _, a in TRUTH], atol=1e-4)
    assert fit.r_squared > 0.9999


def test_areas_are_analytic_not_sampled(amide_i):
    """Area must not depend on how finely the axis is sampled."""
    x, y = amide_i
    coarse_x = np.linspace(1600.0, 1700.0, 101)
    coarse_y = sum(gauss(coarse_x, p, w, a) for p, w, a in TRUTH)

    fine = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    coarse = fit_components(coarse_x, coarse_y, [p for p, _, _ in TRUTH],
                            position_tolerance=4)
    np.testing.assert_allclose(fine.area, coarse.area, rtol=1e-4)


def test_area_matches_the_closed_form(amide_i):
    """Gaussian area is amplitude * fwhm * sqrt(pi / 4 ln2)."""
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    expected = [a * w * np.sqrt(np.pi / (4 * np.log(2))) for _, w, a in TRUTH]
    np.testing.assert_allclose(fit.area, expected, rtol=1e-3)


def test_fractions_sum_to_one(amide_i):
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    assert fit.fractions().sum() == pytest.approx(1.0)


def test_well_determined_components_survive_noise(amide_i):
    x, y = amide_i
    noisy = y + np.random.default_rng(0).normal(0, 0.005, len(y))
    fit = fit_components(x, noisy, [p for p, _, _ in TRUTH], position_tolerance=5)

    # The three that sit on their own maxima come back.
    for index, expected in [(0, 1628.0), (2, 1655.0), (3, 1672.0)]:
        assert fit.position[index] == pytest.approx(expected, abs=0.6)
    assert fit.r_squared > 0.999


def test_a_poorly_determined_component_says_so(amide_i):
    """The uncertainty is the point, not a formality.

    The 1645 component is weak (0.30) and sits between two strong ones, so
    with 0.5% noise it drifts several wavenumbers -- and under a tolerance of
    5 it will happily sit on the bound. Fitting cannot rescue information the
    data does not contain; what it can do is report that the number is soft,
    which is why stderr exists and why anything reading `position` should read
    `stderr['position']` too.
    """
    x, y = amide_i
    noisy = y + np.random.default_rng(0).normal(0, 0.005, len(y))
    fit = fit_components(x, noisy, [p for p, _, _ in TRUTH], position_tolerance=5)

    assert fit.stderr is not None
    soft = fit.stderr['position'][1]
    firm = np.delete(fit.stderr['position'], 1)
    assert soft == fit.stderr['position'].max()
    assert soft > 5 * np.median(firm), (
        "the weak shoulder should carry a visibly larger uncertainty than the "
        f"resolved bands; got {soft:.3f} against a median of {np.median(firm):.3f}"
    )


def test_the_curve_reproduces_the_data(amide_i):
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    np.testing.assert_allclose(fit.curve(), y, atol=1e-6)
    assert fit.components().shape == (4, len(x))


def test_it_can_be_evaluated_on_a_different_axis(amide_i):
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    assert fit.curve(np.linspace(1600, 1700, 51)).shape == (51,)


# ---------------------------------------------------------------------------
# models
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('model, shape', [
    ('gaussian', gauss),
    ('lorentzian', lorentz),
])
def test_each_model_recovers_its_own_shape(model, shape):
    x = np.linspace(1600.0, 1700.0, 401)
    y = shape(x, 1650.0, 15.0, 1.0)
    fit = fit_components(x, y, [1650.0], model=model)
    assert fit.position[0] == pytest.approx(1650.0, abs=1e-3)
    assert fit.fwhm[0] == pytest.approx(15.0, abs=1e-3)
    assert fit.r_squared > 0.9999


def test_voigt_recovers_its_mixing_fraction():
    x = np.linspace(1600.0, 1700.0, 401)
    y = spec_comp(x, 1650.0, 15.0, 1.0, 0.3)
    fit = fit_components(x, y, [1650.0], model='voigt')
    assert fit.eta[0] == pytest.approx(0.3, abs=0.02)
    assert fit.r_squared > 0.9999


def test_pseudo_voigt_peak_height_is_the_amplitude():
    """The lineshapes bug fixed in 0.1.1.

    spec_comp multiplied the Lorentzian half by ``ext`` twice, so a mixed
    component peaked at ``fg * ext + (1 - fg) * ext**2``. With ext = 1 -- the
    obvious thing to try -- the error is invisible.
    """
    for fraction in (0.0, 0.25, 0.5, 1.0):
        peak = spec_comp(np.array([1650.0]), 1650.0, 15.0, 2.5, fraction)
        assert peak[0] == pytest.approx(2.5), f"wrong height at fg={fraction}"


# ---------------------------------------------------------------------------
# constraints
# ---------------------------------------------------------------------------

def test_position_tolerance_holds_components_in_place():
    """Two components on one band collapse together unless constrained.

    This is the failure that makes an unconstrained amide I fit report
    structure that is not there.
    """
    x = np.linspace(1600.0, 1700.0, 401)
    y = gauss(x, 1650.0, 15.0, 1.0)

    free = fit_components(x, y, [1640.0, 1660.0])
    held = fit_components(x, y, [1640.0, 1660.0], position_tolerance=2.0)

    assert abs(free.position[0] - free.position[1]) < 6.0
    assert abs(held.position[0] - held.position[1]) > 15.0


def test_non_negative_is_the_default(amide_i):
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    assert np.all(fit.amplitude >= 0)


def test_negative_amplitudes_are_allowed_when_asked_for():
    """A difference spectrum has real negative bands."""
    x = np.linspace(1600.0, 1700.0, 401)
    y = gauss(x, 1655.0, 12.0, 0.8) - gauss(x, 1630.0, 12.0, 0.5)
    fit = fit_components(x, y, [1630.0, 1655.0], position_tolerance=3,
                         non_negative=False)
    assert fit.amplitude.min() < 0
    assert fit.r_squared > 0.999


# ---------------------------------------------------------------------------
# assignment
# ---------------------------------------------------------------------------

def test_assign_groups_components_by_range(amide_i):
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    structure = fit.assign({'beta-sheet': (1620, 1640), 'random': (1640, 1650),
                            'alpha-helix': (1650, 1660), 'turns': (1660, 1685)})
    assert sum(structure.values()) == pytest.approx(1.0)
    assert structure['alpha-helix'] > structure['beta-sheet']
    assert 'unassigned' not in structure


def test_assign_never_silently_drops_a_component(amide_i):
    """A band that matches nothing is the interesting one."""
    x, y = amide_i
    fit = fit_components(x, y, [p for p, _, _ in TRUTH], position_tolerance=4)
    structure = fit.assign({'alpha-helix': (1650, 1660)})
    assert 'unassigned' in structure
    assert sum(structure.values()) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# refusals
# ---------------------------------------------------------------------------

def test_positions_outside_the_data_are_refused(amide_i):
    x, y = amide_i
    with pytest.raises(ValueError, match='outside the data'):
        fit_components(x, y, [1200.0])


def test_no_positions_is_refused(amide_i):
    x, y = amide_i
    with pytest.raises(ValueError, match='at least one'):
        fit_components(x, y, [])


def test_unknown_model_is_refused(amide_i):
    x, y = amide_i
    with pytest.raises(ValueError, match='unknown model'):
        fit_components(x, y, [1650.0], model='trapezoid')


# ---------------------------------------------------------------------------
# through Spectrum
# ---------------------------------------------------------------------------

def test_spectrum_fit_peaks(amide_i):
    x, y = amide_i
    fit = Spectrum(x, y, technique='ATR-FTIR', name='synthetic').fit_peaks(
        [p for p, _, _ in TRUTH], position_tolerance=4)
    assert isinstance(fit, FitResult)
    assert fit.x_unit == 'cm^-1'
    assert fit.source == 'synthetic'
    assert 'synthetic' in repr(fit)


def test_fit_peaks_finds_its_own_starting_positions():
    """Well-separated bands need no help."""
    x = np.linspace(1600.0, 1700.0, 401)
    y = gauss(x, 1630.0, 8.0, 0.6) + gauss(x, 1670.0, 8.0, 0.9)
    fit = Spectrum(x, y).fit_peaks(n_peaks=2, position_tolerance=6)
    np.testing.assert_allclose(fit.position, [1630.0, 1670.0], atol=1.0)


def test_str_is_a_readable_table(amide_i):
    x, y = amide_i
    text = str(Spectrum(x, y, technique='ATR-FTIR').fit_peaks(
        [p for p, _, _ in TRUTH], position_tolerance=4))
    assert 'position' in text and 'fraction' in text
    assert '1628' in text


# ---------------------------------------------------------------------------
# derivative weighting and lineshape -- the two things that decide whether a
# composition is right. Both measured rather than assumed; see the numbers in
# the fit_components docstring.
# ---------------------------------------------------------------------------

def _overlapping_pair():
    """Two bands 10 cm^-1 apart with fwhm 16 -- one maximum, one shoulder."""
    truth = [(1643.0, 16.0, 0.55), (1653.0, 16.0, 0.85)]
    x = np.linspace(1600.0, 1700.0, 401)
    y = sum(gauss(x, p, w, a) for p, w, a in truth)
    return x, y, truth


def test_the_pair_really_does_overlap():
    """The premise: if the envelope had two maxima none of this would matter."""
    _, y, _ = _overlapping_pair()
    maxima = (np.diff(np.sign(np.diff(y))) < 0).sum()
    assert maxima == 1


def test_derivative_weighting_rescues_a_fit_from_a_residual_background():
    """The reason to fit d2A/dv2 as well as A.

    A second derivative annihilates a smooth background; the absorbance
    envelope absorbs it into the component areas instead. This is why the
    weighting matters on real data, where water subtraction is never perfect,
    and why it does nothing on a synthetic band with no background at all.
    """
    x, bands, truth = _overlapping_pair()
    curvature = (x - x.mean()) / 50.0
    y = bands + 0.10 * curvature ** 2 + 0.04 * curvature

    plain = fit_components(x, y, [1640.0, 1656.0], model='gaussian',
                           position_tolerance=8)
    weighted = fit_components(x, y, [1640.0, 1656.0], model='gaussian',
                              position_tolerance=8, derivative_weight=2.0)

    expected = np.array([p for p, _, _ in truth])
    plain_error = np.abs(plain.position - expected).max()
    weighted_error = np.abs(weighted.position - expected).max()
    assert weighted_error < plain_error / 2, (
        f"derivative weighting should more than halve the position error on a "
        f"residual background; got {plain_error:.2f} -> {weighted_error:.2f}"
    )


def test_an_uneven_axis_is_refused_for_derivative_weighting():
    x = np.concatenate([np.linspace(1600, 1650, 200), np.linspace(1651, 1700, 50)])
    y = gauss(x, 1650.0, 15.0, 1.0)
    with pytest.raises(ValueError, match='uniformly spaced'):
        fit_components(x, y, [1650.0], derivative_weight=1.0)


def test_a_good_r_squared_does_not_validate_the_lineshape():
    """Why the default is pseudo-Voigt.

    Fit intermediate-shaped bands with pure Gaussians and the fit looks
    excellent while the composition is wrong by tens of points. Nothing in the
    residual tells the user; only using the wrong shape does.
    """
    x = np.linspace(1600.0, 1700.0, 401)
    truth = [(1628.0, 14.0, 0.55), (1645.0, 16.0, 0.30),
             (1655.0, 12.0, 0.90), (1672.0, 14.0, 0.35)]
    y = sum(spec_comp(x, p, w, a, 0.5) for p, w, a in truth)

    positions = [p for p, _, _ in truth]
    wrong = fit_components(x, y, positions, model='gaussian',
                           position_tolerance=4, derivative_weight=2.0)
    right = fit_components(x, y, positions, model='voigt',
                           position_tolerance=4, derivative_weight=2.0)

    assert wrong.r_squared > 0.97, "the wrong-shape fit should still look good"

    true_areas = np.array([
        (0.5 * a * w * np.sqrt(np.pi / (4 * np.log(2))) + 0.5 * a * w * np.pi / 2)
        for _, w, a in truth])
    true_fractions = true_areas / true_areas.sum()

    wrong_error = np.abs(wrong.fractions() - true_fractions).max()
    right_error = np.abs(right.fractions() - true_fractions).max()
    assert wrong_error > 0.15, "expected the pure-Gaussian fit to be badly wrong"
    assert right_error < 0.02
    assert right_error < wrong_error / 5


def test_voigt_costs_nothing_when_the_bands_really_are_gaussian():
    """The other half of the default: floating the mixing is not a risk."""
    x = np.linspace(1600.0, 1700.0, 401)
    y = gauss(x, 1650.0, 15.0, 1.0)
    fit = fit_components(x, y, [1650.0], model='voigt')
    assert fit.eta[0] == pytest.approx(1.0, abs=0.02)
    assert fit.position[0] == pytest.approx(1650.0, abs=1e-2)


def test_voigt_is_the_default():
    x = np.linspace(1600.0, 1700.0, 401)
    y = gauss(x, 1650.0, 15.0, 1.0)
    assert fit_components(x, y, [1650.0]).model == 'voigt'
    assert Spectrum(x, y).fit_peaks([1650.0]).model == 'voigt'
