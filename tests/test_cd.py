# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Circular dichroism: structure from spectrum shape, and thermal melts.

Everything here is tested against synthetic data built from a known basis or a
known transition, so that the recovered answer can be compared with a truth.
That proves the arithmetic. It does not prove that any real basis describes a
real protein -- see CD_Branch_Plan.md WP5, which is what real reference
proteins are for.
"""

import numpy as np
import pytest

import spectroscopy as spc
from spectroscopy import units
from spectroscopy.processing import melting, structure
from spectroscopy.processing.structure import Category, from_cd

X = np.linspace(190.0, 250.0, 121)
GAS_CONSTANT = 8.314462618


def _band(x, centre, width, height):
    return height * np.exp(-((x - centre) ** 2) / (2 * width ** 2))


#: Textbook far-UV shapes. Not reference data -- shapes with the right signs
#: and rough positions, enough to test that a mixture comes back apart.
HELIX = Category('helix', frozenset({'G', 'H', 'I'}))
SHEET = Category('sheet', frozenset({'E', 'B'}))
COIL = Category('other', frozenset({'S', '-'}))

SHAPES = {
    HELIX: _band(X, 193, 7, 75) + _band(X, 208, 7, -37) + _band(X, 222, 9, -37),
    SHEET: _band(X, 195, 8, 32) + _band(X, 217, 9, -18),
    COIL: _band(X, 198, 7, -40) + _band(X, 220, 10, 3),
}


@pytest.fixture
def basis():
    spectra = []
    for category, y in SHAPES.items():
        spectrum = spc.Spectrum(X, y, technique='CD', name=category.name)
        spectrum.metadata['category'] = category
        spectra.append(spectrum)
    return spectra


def _mixture(fractions, name='synthetic'):
    y = sum(fractions[c.name] * shape for c, shape in SHAPES.items())
    return spc.Spectrum(X, y, technique='CD', name=name)


# ---------------------------------------------------------------------------
# CD as a technique
# ---------------------------------------------------------------------------

def test_cd_has_its_own_axis_conventions():
    spectrum = spc.Spectrum(X, np.zeros_like(X), technique='CD')
    assert spectrum.x_unit == 'nm'
    assert spectrum.y_unit == 'mdeg'
    assert spectrum.y_quantity == 'Ellipticity'
    # Unlike infrared, CD is plotted short-to-long wavelength.
    assert not spectrum.reversed_x


@pytest.mark.parametrize('unit', ['mdeg', 'deg', 'deg cm^2 dmol^-1',
                                  'delta epsilon'])
def test_every_cd_unit_is_bipolar(unit):
    """A helix is negative at 222 and positive at 193, in one spectrum."""
    assert units.band_direction(unit) == 'both'


def test_peaks_of_a_cd_spectrum_are_found_both_ways():
    spectrum = _mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0})
    peaks = spectrum.find_peaks(method='direct', prominence=2)
    assert peaks.kind == 'both'
    assert len(peaks.maxima()) and len(peaks.minima())
    assert np.any(np.isclose(peaks.maxima().position, 193.0, atol=3.0))


# ---------------------------------------------------------------------------
# structure from the shape
# ---------------------------------------------------------------------------

def test_a_known_mixture_comes_back_apart(basis):
    truth = {'helix': 0.55, 'sheet': 0.20, 'other': 0.25}
    result = from_cd(_mixture(truth), 'basis-spectra', basis=basis)

    for category, fraction in result.fractions.items():
        assert fraction == pytest.approx(truth[category.name], abs=1e-3)
    assert result.technique == 'CD'
    assert result.method == 'basis-spectra'


def test_the_fit_is_scale_free(basis):
    """
    Shape, not size: an unknown concentration must not prevent an answer.
    Getting to mean residue ellipticity needs three sample facts, and
    requiring them before any structure could be estimated would block the
    common case for a reason that does not apply to it.
    """
    truth = {'helix': 0.55, 'sheet': 0.20, 'other': 0.25}
    normal = from_cd(_mixture(truth), 'basis-spectra', basis=basis)
    scaled = _mixture(truth)
    scaled.y = scaled.y * 37.0
    louder = from_cd(scaled, 'basis-spectra', basis=basis)

    for category in normal.fractions:
        assert louder.fractions[category] == pytest.approx(
            normal.fractions[category], abs=1e-6)


def test_fractions_sum_to_one_and_are_never_negative(basis):
    """A negative fraction of helix is not a small number, it is a wrong fit."""
    result = from_cd(_mixture({'helix': 0.9, 'sheet': 0.05, 'other': 0.05}),
                     'basis-spectra', basis=basis)
    assert sum(result.fractions.values()) == pytest.approx(1.0, abs=1e-6)
    assert all(value >= 0.0 for value in result.fractions.values())


def test_reference_proteins_carry_their_own_compositions():
    """The other kind of standard: fit proteins, then mix their structures."""
    first = {'helix': 0.80, 'sheet': 0.05, 'other': 0.15}
    second = {'helix': 0.10, 'sheet': 0.60, 'other': 0.30}
    proteins, compositions = [], []
    for index, truth in enumerate((first, second)):
        proteins.append(_mixture(truth, f'ref{index}'))
        compositions.append(structure.Composition(
            fractions={c: truth[c.name] for c in SHAPES},
            method='dssp', technique='X-ray'))

    unknown = spc.Spectrum(X, 0.25 * proteins[0].y + 0.75 * proteins[1].y,
                           technique='CD', name='unknown')
    result = from_cd(unknown, 'reference-proteins', basis=proteins,
                     compositions=compositions)

    for category in SHAPES:
        expected = 0.25 * first[category.name] + 0.75 * second[category.name]
        assert result.fractions[category] == pytest.approx(expected, abs=1e-3)


def test_the_method_must_be_named(basis):
    """The two kinds of standard answer differently; the result must say."""
    with pytest.raises(ValueError, match='method must be one of'):
        from_cd(_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}),
                basis=basis)


def test_no_basis_says_where_to_get_one():
    with pytest.raises(ValueError, match='no reference spectra ship'):
        from_cd(_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}),
                'basis-spectra')


def test_a_basis_spectrum_must_declare_its_category():
    plain = spc.Spectrum(X, SHAPES[HELIX], technique='CD', name='helix')
    with pytest.raises(ValueError, match="metadata\\['category'\\]"):
        from_cd(_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}),
                'basis-spectra', basis=[plain])


def test_a_basis_that_cannot_describe_the_spectrum_shows_in_the_rmsd(basis):
    """
    The number to look at before believing any fraction. A fit against a basis
    missing the dominant component still returns fractions summing to one.
    """
    good = from_cd(_mixture({'helix': 0.55, 'sheet': 0.20, 'other': 0.25}),
                   'basis-spectra', basis=basis)
    unusual = spc.Spectrum(X, _band(X, 230, 6, 40) - _band(X, 205, 5, 25),
                           technique='CD', name='not in the basis')
    poor = from_cd(unusual, 'basis-spectra', basis=basis)

    assert sum(poor.fractions.values()) == pytest.approx(1.0, abs=1e-6)
    assert poor.quality['rmsd_relative'] > 50 * good.quality['rmsd_relative']


# ---------------------------------------------------------------------------
# melting curves
# ---------------------------------------------------------------------------

TEMPERATURES = np.arange(20.0, 90.1, 2.5)
MELT_X = np.linspace(200.0, 250.0, 101)
_FOLDED = _band(MELT_X, 208, 7, -37) + _band(MELT_X, 222, 9, -37)
_INTERMEDIATE = _band(MELT_X, 212, 8, -28) + _band(MELT_X, 228, 9, -14)
_UNFOLDED = _band(MELT_X, 203, 7, -25) + _band(MELT_X, 222, 10, -1.5)


def _melt(three_state=False, tm=52.0, enthalpy=250.0, noise=0.15, seed=0):
    generator = np.random.default_rng(seed)
    spectra = []
    for temperature in TEMPERATURES:
        kelvin = temperature + 273.15
        if three_state:                      # F -> I at 42 C, I -> U at 62 C
            first = np.exp(-(220e3 / GAS_CONSTANT)
                           * (1 / kelvin - 1 / (42 + 273.15)))
            second = np.exp(-(220e3 / GAS_CONSTANT)
                            * (1 / kelvin - 1 / (62 + 273.15)))
            y = ((_FOLDED + first * _INTERMEDIATE + first * second * _UNFOLDED)
                 / (1 + first + first * second))
        else:
            equilibrium = np.exp(-(enthalpy * 1e3 / GAS_CONSTANT)
                                 * (1 / kelvin - 1 / (tm + 273.15)))
            y = (_FOLDED + equilibrium * _UNFOLDED) / (1 + equilibrium)
        spectra.append(spc.Spectrum(
            MELT_X, y + noise * generator.normal(size=MELT_X.size),
            technique='CD', name=f'{temperature:g} C'))
    return spc.SpectrumCollection(spectra, name='melt').with_parameters(
        TEMPERATURES, name='temperature', unit='C')


def test_tm_and_enthalpy_come_back():
    result = melting.from_collection(_melt())
    assert result.tm == pytest.approx(52.0, abs=0.5)
    assert result.enthalpy == pytest.approx(250.0, rel=0.05)


def test_a_single_wavelength_agrees_with_the_whole_spectrum():
    whole = melting.from_collection(_melt())
    at_222 = melting.from_collection(_melt(), 222.0)
    assert at_222.tm == pytest.approx(whole.tm, abs=0.5)


def test_sloping_baselines_are_fitted_not_assumed_flat():
    """
    A state's own signal drifts with temperature. Forced flat, that drift is
    pushed into the transition and moves Tm.
    """
    result = melting.from_collection(_melt())
    assert len(result.folded) == 2 and len(result.unfolded) == 2
    assert np.isfinite(result.folded[1]) and np.isfinite(result.unfolded[1])


def test_fraction_unfolded_is_a_half_at_tm():
    result = melting.from_collection(_melt())
    assert result.fraction_unfolded(result.tm) == pytest.approx(0.5, abs=1e-9)
    assert result.fraction_unfolded(0.0) < 0.01
    assert result.fraction_unfolded(100.0) > 0.99


def test_kelvin_is_refused_rather_than_fitted():
    """It would otherwise succeed and report a Tm around 325 C."""
    with pytest.raises(ValueError, match='Kelvin rather than Celsius'):
        melting.two_state(TEMPERATURES + 273.15,
                          np.linspace(0.0, 1.0, TEMPERATURES.size))


def test_a_transition_the_measurement_never_reached_warns():
    """
    Stopping at 25 C when Tm is 60 leaves a curve that never turned over. The
    baselines can still fit it, and Tm can still land inside the range, so
    the honest test is how much of the transition actually happened.
    """
    cold = np.arange(5.0, 25.1, 1.0)
    kelvin = cold + 273.15
    equilibrium = np.exp(-(250e3 / GAS_CONSTANT)
                         * (1 / kelvin - 1 / (60 + 273.15)))
    fraction = equilibrium / (1 + equilibrium)
    signal = 10.0 - 0.02 * cold + fraction * (2.0 - 10.0)

    with pytest.warns(UserWarning, match='of the transition happened'):
        melting.two_state(cold, signal)


def test_too_few_points_for_six_parameters():
    with pytest.raises(ValueError, match='six parameters'):
        melting.two_state([20.0, 40.0, 60.0], [1.0, 0.5, 0.0])


# -- the two-state checks, which are the point of the module ----------------

def test_a_two_state_melt_passes_both_checks():
    series = _melt()
    assert melting.isodichroic_point(series)['is_tight']
    assert melting.two_state_rank(series)['is_two_state']


def test_a_three_state_melt_fails_both_checks():
    series = _melt(three_state=True)
    assert not melting.isodichroic_point(series)['is_tight']
    assert not melting.two_state_rank(series)['is_two_state']


def test_the_three_state_melt_still_fits_acceptably():
    """
    Why the checks exist. A two-state curve describes three-state data at one
    wavelength with an unremarkable residual, and reports a single Tm for a
    protein with two transitions: it finds the one at 62 C, and the one at
    42 C leaves no trace in Tm at all.

    The fitted dH does carry a signature -- far below the 220 kJ/mol that
    built each step, which is the classical sign of a transition that is not
    two-state. But reading it needs an expectation to compare against, which
    an unknown protein does not come with. The two structural checks need no
    such expectation, and that is the difference.
    """
    series = _melt(three_state=True)
    result = melting.from_collection(series, 222.0)
    span = float(np.ptp(series.to_matrix()[1]))

    assert result.residual_rms / span < 0.02        # residual looks fine
    assert abs(result.tm - 42.0) > 15.0             # the first is invisible
    assert result.enthalpy < 0.5 * 220.0            # the dH signature
    assert melting.two_state_rank(series)['third_over_noise'] > 5.0


def test_the_rank_check_needs_enough_spectra():
    few = _melt()[:4]
    with pytest.raises(ValueError, match='at least five'):
        melting.two_state_rank(few)


def test_free_energy_and_entropy_are_consistent():
    result = melting.from_collection(_melt())
    assert result.free_energy(result.tm) == pytest.approx(0.0, abs=1e-9)
    assert result.entropy == pytest.approx(
        result.enthalpy * 1000.0 / (result.tm + 273.15))
