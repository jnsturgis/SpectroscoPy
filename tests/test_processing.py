# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Processing algorithms, tested against synthetic ground truth.

Roadmap section 6: "synthetic spectra with known ground truth (e.g. a synthetic
Gaussian peak on a known baseline) to verify baseline correction/peak fitting
recover the right answer within tolerance."
"""

import numpy as np
import pytest

from spectroscopy.lineshapes import gauss
from spectroscopy.processing import common
from spectroscopy.spectra import Spectrum


@pytest.fixture
def synthetic():
    """Two Gaussians on a known sloping baseline, in FTIR-like units."""
    x = np.linspace(900.0, 1800.0, 901)
    peaks = gauss(x, 1650.0, 40.0, 1.0) + gauss(x, 1540.0, 30.0, 0.6)
    baseline = 0.15 + 2.0e-4 * (x - 900.0)
    spec = Spectrum()
    spec.x, spec.y = x, peaks + baseline
    spec.set_type("ATR-FTIR")
    spec.set_sample("synthetic")
    return spec, x, peaks, baseline


# ---------------------------------------------------------------------------
# baselines
# ---------------------------------------------------------------------------

def test_rubberband_recovers_a_linear_baseline(synthetic):
    _, x, peaks, baseline = synthetic
    estimate = common.rubberband_baseline(x, peaks + baseline)
    # Hull-based, so it sits at or below the true baseline; it should be close
    # and never above it by more than numerical slack.
    assert np.max(estimate - baseline) < 1e-6
    assert np.mean(np.abs(estimate - baseline)) < 0.02


def test_rubberband_returns_usable_anchor_points(synthetic):
    _, x, peaks, baseline = synthetic
    base, anchor_x, anchor_y = common.rubberband_baseline(
        x, peaks + baseline, return_points=True)
    assert len(anchor_x) >= 2
    assert len(anchor_x) == len(anchor_y)
    assert anchor_x.min() >= x.min() and anchor_x.max() <= x.max()
    # The anchors are exactly what a polynomial baseline can be fitted through.
    refit = common.poly_baseline(x, peaks + baseline, degree=1, points=anchor_x)
    assert np.mean(np.abs(refit - baseline)) < 0.05
    assert len(base) == len(x)


def test_poly_baseline_from_guide_points(synthetic):
    """The SMALP-notebook form: anchor positions, not coefficients."""
    _, x, peaks, baseline = synthetic
    guides = [900, 950, 1000, 1150, 1750, 1800]      # away from both bands
    estimate = common.poly_baseline(x, peaks + baseline, degree=1, points=guides)
    assert np.allclose(estimate, baseline, atol=0.01)


def test_poly_baseline_from_coefficients(synthetic):
    _, x, _, _ = synthetic
    estimate = common.poly_baseline(x, np.zeros_like(x), coefficients=[2.0e-4, 0.03])
    assert np.allclose(estimate, 2.0e-4 * x + 0.03)


def test_poly_baseline_needs_enough_guide_points(synthetic):
    _, x, peaks, baseline = synthetic
    with pytest.raises(ValueError, match="more than degree"):
        common.poly_baseline(x, peaks + baseline, degree=3, points=[900, 1800])


def test_poly_baseline_needs_points_or_coefficients(synthetic):
    _, x, peaks, baseline = synthetic
    with pytest.raises(ValueError, match="points.*or.*coefficients"):
        common.poly_baseline(x, peaks + baseline)


def test_als_baseline_runs_and_sits_under_the_peaks(synthetic):
    """
    The notebooks' copy of this never executed -- it used `sparse`/`spsolve`
    without importing them and was never called.
    """
    _, x, peaks, baseline = synthetic
    estimate = common.als_baseline(peaks + baseline, lam=1e5, p=0.001, niter=10)
    assert len(estimate) == len(x)
    at_peak = np.argmin(np.abs(x - 1650.0))
    assert estimate[at_peak] < (peaks + baseline)[at_peak]
    off_band = (x < 1200)
    assert np.mean(np.abs(estimate[off_band] - baseline[off_band])) < 0.05


def test_unknown_baseline_method_is_rejected(synthetic):
    _, x, peaks, baseline = synthetic
    with pytest.raises(ValueError, match="Unknown baseline method"):
        common.baseline(x, peaks + baseline, method="magic")


# ---------------------------------------------------------------------------
# smoothing, derivatives, normalisation
# ---------------------------------------------------------------------------

def test_smoothing_reduces_noise_but_keeps_the_peak(synthetic):
    _, x, peaks, _ = synthetic
    rng = np.random.default_rng(0)
    noisy = peaks + rng.normal(0.0, 0.02, len(x))
    smoothed = common.smooth(noisy, window_length=15, polyorder=3)
    assert np.std(smoothed - peaks) < np.std(noisy - peaks)
    assert abs(x[np.argmax(smoothed)] - 1650.0) < 5.0


def test_even_savgol_window_is_corrected():
    """
    Every notebook calls savgol_filter(y, 10, 3, deriv=2) -- an even window,
    almost certainly a slip for 11, propagated by copy-paste.
    """
    y = np.sin(np.linspace(0, 10, 200))
    assert len(common.smooth(y, window_length=10, polyorder=3)) == len(y)
    assert len(common.derivative(y, window_length=10, polyorder=3)) == len(y)


def test_second_derivative_is_negative_at_a_maximum(synthetic):
    _, x, peaks, _ = synthetic
    second = common.derivative(peaks, order=2, window_length=11, polyorder=3)
    assert second[np.argmin(np.abs(x - 1650.0))] < 0


@pytest.mark.parametrize("method", ["max", "area", "minmax", "vector", "snv"])
def test_normalisation_methods_run(synthetic, method):
    _, x, peaks, _ = synthetic
    assert len(common.normalize(x, peaks, method=method)) == len(x)


def test_normalise_to_max_in_a_window(synthetic):
    _, x, peaks, _ = synthetic
    scaled = common.normalize(x, peaks, method='max', window=(1520, 1560))
    assert np.isclose(np.max(scaled[(x >= 1520) & (x <= 1560)]), 1.0)


def test_normalisation_window_is_order_insensitive(synthetic):
    """FTIR ranges get quoted high-to-low as often as low-to-high."""
    _, x, peaks, _ = synthetic
    assert np.allclose(common.normalize(x, peaks, window=(1800, 900)),
                       common.normalize(x, peaks, window=(900, 1800)))


def test_empty_normalisation_window_is_rejected(synthetic):
    _, x, peaks, _ = synthetic
    with pytest.raises(ValueError, match="no points"):
        common.normalize(x, peaks, window=(10, 20))


# ---------------------------------------------------------------------------
# peak detection
# ---------------------------------------------------------------------------

def test_second_derivative_finds_both_bands(synthetic):
    spec, _, _, _ = synthetic
    corrected = spec.baseline_correct('rubberband')
    table = corrected.find_peaks(prominence=1e-5, distance=5)
    found = sorted(table.position)
    assert any(abs(p - 1540.0) < 10 for p in found), found
    assert any(abs(p - 1650.0) < 10 for p in found), found


def test_second_derivative_finds_a_shoulder_that_direct_detection_misses():
    """The reason every notebook picks peaks off the second derivative."""
    # A genuine shoulder: weak and close enough that the sum has only ONE
    # local maximum, so direct detection cannot see the second band at all.
    x = np.linspace(1500.0, 1750.0, 501)
    y = gauss(x, 1650.0, 50.0, 1.0) + gauss(x, 1615.0, 30.0, 0.20)
    spec = Spectrum()
    spec.x, spec.y = x, y

    direct = spec.find_peaks(method='direct')
    second = spec.find_peaks(method='second_derivative', prominence=1e-7)

    assert len(direct) == 1, f"expected a single maximum, got {direct.position}"
    assert not any(abs(p - 1615.0) < 12 for p in direct.position)
    assert any(abs(p - 1615.0) < 12 for p in second.position), second.position


def test_troughs(synthetic):
    spec, _, _, _ = synthetic
    inverted = spec * -1.0
    table = inverted.find_peaks(troughs=True, prominence=1e-5)
    assert table.kind == 'trough'
    assert len(table) > 0


def test_rubberband_survives_a_degenerate_spectrum():
    """
    A flat or perfectly linear spectrum has no 2-D convex hull. Without a
    fallback scipy raises QhullError with a page of diagnostics.
    """
    x = np.linspace(900.0, 1800.0, 91)

    flat = common.rubberband_baseline(x, np.full_like(x, 3.0))
    assert np.allclose(flat, 3.0)

    line = 0.5 + 1e-3 * x
    assert np.allclose(common.rubberband_baseline(x, line), line)

    base, anchor_x, anchor_y = common.rubberband_baseline(
        x, np.full_like(x, 3.0), return_points=True)
    assert len(anchor_x) == len(anchor_y) == 2
    assert np.allclose(base, 3.0)


def test_even_savgol_window_biases_peak_positions_by_half_a_sample():
    """
    Why the even window is corrected rather than preserved.

    Every notebook calls savgol_filter(y, 10, 3, deriv=2). An even-length
    Savitzky-Golay window is asymmetric about the centre point, which gives the
    filter a half-sample group delay and biases every recovered peak position
    low by half the sample spacing. On the real .dpt sampling (~1.93 cm-1) that
    is about -1 cm-1 on every reported band.
    """
    from scipy.signal import find_peaks as scipy_find_peaks
    from scipy.signal import savgol_filter as scipy_savgol

    x = np.arange(900.0, 1800.0, 1.928)
    spacing = float(np.diff(x)[0])

    errors = {10: [], 11: []}
    for true_position in [1000.0, 1250.5, 1400.0, 1550.7]:
        y = gauss(x, true_position, 30.0, 1.0)
        for window in (10, 11):
            second = -scipy_savgol(y, window, 3, deriv=2)
            found = scipy_find_peaks(second)[0]
            nearest = x[found][np.argmin(np.abs(x[found] - true_position))]
            errors[window].append(nearest - true_position)

    # The even window is biased low by roughly half a sample ...
    assert np.mean(errors[10]) < -0.3 * spacing
    # ... the odd one is not, and is more accurate overall.
    assert abs(np.mean(errors[11])) < 0.3 * spacing
    assert np.sqrt(np.mean(np.square(errors[11]))) < \
        np.sqrt(np.mean(np.square(errors[10])))

    # Which is what _as_odd_at_least guarantees for us.
    assert common._as_odd_at_least(10, 3) == 11
    assert common._as_odd_at_least(11, 3) == 11
    assert common._as_odd_at_least(2, 3) == 5


def test_peak_thresholds_apply_to_the_detection_signal(synthetic):
    """
    A trap worth pinning: with the default method, height/prominence apply to
    the *second derivative*, which is orders of magnitude smaller than the
    spectrum. A value that looks sensible in spectrum units finds nothing.
    This is where the unexplained small constants in the notebooks come from.
    """
    spec, _, _, _ = synthetic
    corrected = spec.baseline_correct('rubberband').normalize('max')

    assert len(corrected.find_peaks(prominence=0.02)) == 0        # spectrum units
    assert len(corrected.find_peaks(prominence=1e-4)) > 0         # d2 units


def test_relative_thresholds_are_scale_free(synthetic):
    """relative=True means 'a fraction of the detection signal's range'."""
    spec, _, _, _ = synthetic
    corrected = spec.baseline_correct('rubberband').normalize('max')

    assert len(corrected.find_peaks(prominence=0.02, relative=True)) > 0

    # Fewer peaks as the bar rises, and the same answer whatever the spectrum
    # is scaled by -- which is the point of it being relative.
    loose = len(corrected.find_peaks(prominence=0.02, relative=True))
    strict = len(corrected.find_peaks(prominence=0.30, relative=True))
    assert strict < loose

    scaled = corrected * 1000.0
    assert len(scaled.find_peaks(prominence=0.02, relative=True)) == loose


def test_relative_works_for_direct_detection_too(synthetic):
    spec, _, _, _ = synthetic
    found = spec.find_peaks(method='direct', prominence=0.1, relative=True)
    assert len(found) > 0
