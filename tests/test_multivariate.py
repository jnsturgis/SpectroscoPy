# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Multivariate decomposition, against mixtures whose components are known.

This is the CK-notebook analysis (PCA/NMF/ICA plus stability testing) moved
into the library. Testing it means building spectra from known component
spectra in known proportions and checking the decomposition recovers them.
"""

import numpy as np
import pytest

from spectroscopy import Spectrum, SpectrumCollection
from spectroscopy.lineshapes import gauss

sklearn = pytest.importorskip("sklearn", reason="needs the [multivariate] extra")

from spectroscopy.processing import multivariate as mv  # noqa: E402

X_AXIS = np.linspace(900.0, 1800.0, 451)
#: An amide-I-like band and a glycan-like pair -- deliberately non-overlapping
#: so a successful separation is unambiguous.
COMPONENT_A = gauss(X_AXIS, 1650.0, 60.0, 1.0)
COMPONENT_B = gauss(X_AXIS, 1080.0, 40.0, 1.0) + gauss(X_AXIS, 1230.0, 30.0, 0.5)


def _mixtures(weights, noise=0.0, seed=0):
    """Build a collection of mixtures of the two known components."""
    generator = np.random.default_rng(seed)
    spectra = []
    for index, (first, second) in enumerate(weights):
        spectrum = Spectrum()
        values = first * COMPONENT_A + second * COMPONENT_B
        if noise:
            values = values + generator.normal(0.0, noise, len(X_AXIS))
        spectrum.x, spectrum.y = X_AXIS.copy(), values
        spectrum.set_type("ATR-FTIR")
        spectrum.set_sample("A-rich" if first > second else "B-rich")
        spectrum.name = f"mix{index}"
        spectra.append(spectrum)
    return SpectrumCollection(spectra)


#: Weights that do NOT sum to a constant. Closed mixtures (fractions summing
#: to 1) are rank-1 once centred, so asking any method for two components is
#: ill-posed -- which is a real trap when composing test data.
WEIGHTS = np.array([
    [1.00, 0.10], [0.80, 0.35], [0.55, 0.60],
    [0.30, 0.85], [0.12, 1.10], [0.90, 0.25],
    [0.45, 0.70], [0.65, 0.15],
])


@pytest.fixture
def mixtures():
    return _mixtures(WEIGHTS)


# ---------------------------------------------------------------------------
# decomposition
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["nmf", "pca", "ica"])
def test_two_components_reconstruct_a_two_component_mixture(mixtures, method):
    result = mv.decompose(mixtures, method=method, n_components=2)
    assert result.explained_variance > 0.999
    assert np.nanmin(result.r2_per_sample) > 0.999
    assert result.reconstruction.shape == (len(mixtures), len(X_AXIS))


def test_nmf_recovers_the_two_component_shapes(mixtures):
    """The components are the point, not just the reconstruction."""
    result = mv.decompose(mixtures, 'nmf', 2)
    peaks = sorted(spectrum.find_peaks(prominence=1e-4).strongest(1).position[0]
                   for spectrum in result.components)
    assert any(abs(p - 1080.0) < 15 for p in peaks), peaks
    assert any(abs(p - 1650.0) < 15 for p in peaks), peaks


def test_components_are_spectra_with_the_collections_axis(mixtures):
    result = mv.decompose(mixtures, 'nmf', 2)
    assert len(result.components) == 2
    for component in result.components:
        assert isinstance(component, Spectrum)
        assert np.allclose(component.x, X_AXIS)
        assert component.x_unit == 'cm^-1'          # inherited from the input
        assert component.history[-1].name == 'decompose'
        # ... so the ordinary API works on them
        assert len(component.crop(1000, 1700)) < len(component)


def test_contribution_fractions_track_the_true_proportions(mixtures):
    """
    The table the CK notebooks build. Absolute values depend on how NMF splits
    the scaling between W and H, so what must hold is the *ordering*: samples
    with relatively more of component A must rank higher on it.
    """
    result = mv.decompose(mixtures, 'nmf', 2)
    fractions = result.contributions()

    assert fractions.shape == (len(mixtures), 2)
    assert np.allclose(fractions.sum(axis=1), 1.0)

    true_ratio = WEIGHTS[:, 0] / WEIGHTS.sum(axis=1)
    which = np.argmax([abs(np.corrcoef(true_ratio, fractions[:, i])[0, 1])
                       for i in range(2)])
    assert abs(np.corrcoef(true_ratio, fractions[:, which])[0, 1]) > 0.98


def test_unnormalised_contributions_scale_with_amount():
    """Doubling a sample's intensity doubles its unnormalised contribution."""
    collection = _mixtures(np.array([[1.0, 0.2], [2.0, 0.4], [0.5, 0.9]]))
    result = mv.decompose(collection, 'nmf', 2)
    absolute = result.contributions(normalize=False)
    assert absolute[1].sum() > absolute[0].sum() > absolute[2].sum() * 0.5


def test_reconstruction_and_residual_come_back_as_spectra(mixtures):
    result = mv.decompose(mixtures, 'nmf', 2)
    rebuilt = result.reconstructed(0)
    residual = result.residual(0)
    assert np.allclose(rebuilt.x, X_AXIS)
    assert np.allclose(rebuilt.y + residual.y, mixtures[0].y, atol=1e-8)


def test_accepts_a_bare_matrix(mixtures):
    x, matrix = mixtures.to_matrix()
    result = mv.decompose((x, matrix), 'pca', 2)
    assert result.n_components == 2


def test_unknown_method_is_refused(mixtures):
    with pytest.raises(ValueError, match="Unknown method"):
        mv.decompose(mixtures, 'magic', 2)


def test_nmf_on_negative_data_says_what_to_do():
    """
    Baseline-corrected spectra routinely dip below zero, so this is the first
    thing that goes wrong on real data. sklearn's own message does not help.
    """
    collection = _mixtures(WEIGHTS[:3])
    shifted = collection.map(lambda s: s - 0.5)
    with pytest.raises(ValueError, match="non-negative"):
        mv.decompose(shifted, 'nmf', 2)
    with pytest.raises(ValueError, match="clip|Baseline"):
        mv.decompose(shifted, 'nmf', 2)


def test_pca_is_happy_with_negative_data():
    collection = _mixtures(WEIGHTS[:4]).map(lambda s: s - 0.5)
    assert mv.decompose(collection, 'pca', 2).explained_variance > 0.99


# ---------------------------------------------------------------------------
# stability
# ---------------------------------------------------------------------------

def test_a_well_posed_fit_is_stable_under_resampling(mixtures):
    """bootstrap is the default because it is the probe that discriminates."""
    result = mv.stability(mixtures, 'nmf', 2, n_trials=8)
    assert result.mode == 'bootstrap'
    assert result.overall > 0.9
    assert result.correlations.shape == (8, 2)
    assert "component 1" in str(result)


def test_reseeding_mode_runs_but_is_weak_for_nmf(mixtures):
    """
    Kept available and documented rather than removed. NMF's initialisation is
    deterministic and coordinate descent converges to the same optimum, so
    reseeding usually returns 1.0 whatever k is -- on the real biofilm data it
    does so for every k from 2 to 8. It is not a useful criterion by itself.
    """
    result = mv.stability(mixtures, 'nmf', 2, mode='runs', n_trials=5)
    assert result.mode == 'runs'
    assert 0.0 <= result.overall <= 1.0


def test_too_many_components_is_less_stable():
    """
    The reason stability is worth measuring: explained variance rises
    monotonically with k and so never tells you to stop, but components that
    are not really there move about between fits.
    """
    # Noise takes the baseline slightly negative, as it does on real
    # baseline-corrected data; clip, which is the remedy the NMF guard names.
    collection = _mixtures(WEIGHTS, noise=0.02, seed=1).map(
        lambda s: s._derive(y=np.clip(s.y, 0.0, None)))    # noqa: SLF001

    honest = mv.stability(collection, 'nmf', 2, n_trials=10)
    excessive = mv.stability(collection, 'nmf', 6, n_trials=10)
    assert honest.overall > excessive.overall


def test_bad_mode_is_refused(mixtures):
    with pytest.raises(ValueError, match="mode must be"):
        mv.stability(mixtures, 'nmf', 2, mode='sideways')


def test_match_components_pairs_and_survives_reordering():
    reference = np.vstack([COMPONENT_A, COMPONENT_B])
    shuffled = np.vstack([COMPONENT_B, COMPONENT_A])
    _, columns, correlations = mv.match_components(reference, shuffled)
    assert list(columns) == [1, 0]
    assert np.all(np.abs(correlations) > 0.99)


def test_match_components_keeps_the_sign():
    """A flipped component is still the same component; the sign says so."""
    reference = np.vstack([COMPONENT_A, COMPONENT_B])
    flipped = np.vstack([-COMPONENT_A, COMPONENT_B])
    _, _, correlations = mv.match_components(reference, flipped)
    assert correlations[0] < -0.99
    assert correlations[1] > 0.99


# ---------------------------------------------------------------------------
# choosing k
# ---------------------------------------------------------------------------

def test_scan_components_reports_a_table(mixtures):
    table = mv.scan_components(mixtures, 'nmf', candidates=range(1, 4),
                               stability_trials=3)
    assert list(table['n_components']) == [1, 2, 3]
    assert table['explained_variance'][1] > table['explained_variance'][0]
    assert 'stability' in table and len(table['stability']) == 3


def test_explained_variance_never_falls_as_k_rises(mixtures):
    """Which is exactly why it cannot be the criterion for choosing k."""
    values = mv.scan_components(mixtures, 'pca',
                                candidates=range(1, 5))['explained_variance']
    assert np.all(np.diff(values) > -1e-9)


# ---------------------------------------------------------------------------
# interop
# ---------------------------------------------------------------------------

def test_to_dataframe_is_the_pandas_escape_hatch(mixtures):
    pandas = pytest.importorskip("pandas")
    frame = mv.decompose(mixtures, 'nmf', 2).to_dataframe()
    assert isinstance(frame, pandas.DataFrame)
    assert list(frame.columns) == ['Comp_1', 'Comp_2']
    assert len(frame) == len(mixtures)


def test_the_core_result_needs_no_pandas(mixtures):
    """Everything above to_dataframe() is numpy, per review section 5.6."""
    result = mv.decompose(mixtures, 'nmf', 2)
    assert isinstance(result.contributions(), np.ndarray)
    assert isinstance(result.weights, np.ndarray)
    assert isinstance(result.r2_per_sample, np.ndarray)
