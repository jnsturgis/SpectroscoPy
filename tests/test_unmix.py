# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Supervised unmixing and the reference library.

Tested against synthetic ground truth, per roadmap section 6: mixtures are
built from known amounts of known references, and the test is whether the
amounts come back. That is the only way to check an unmixing -- goodness of
fit cannot do it, as several tests here demonstrate.
"""

import numpy as np
import pytest

from spectroscopy import Spectrum, SpectrumCollection
from spectroscopy import library as lib
from spectroscopy.library import Library, Reference, from_series
from spectroscopy.lineshapes import gauss
from spectroscopy.processing.unmix import (
    absorbance_ratio,
    best_wavelengths,
    nucleic_acid_and_protein,
    unmix,
)

X = np.linspace(220.0, 340.0, 241)
DNA_EPS = 20.0 * gauss(X, 260, 22, 1.0) + 8.0 * gauss(X, 230, 18, 1.0)
PROTEIN_EPS = 1.0 * gauss(X, 280, 14, 1.0) + 12.0 * gauss(X, 225, 12, 1.0)
HAEM_EPS = 30.0 * gauss(X, 300, 10, 1.0)


def _spectrum(y, name='s'):
    return Spectrum(X, y, x_quantity='Wavelength', x_unit='nm',
                    y_quantity='Absorbance', y_unit='absorbance', name=name)


def _reference(name, epsilon):
    return Reference(name, _spectrum(epsilon, name),
                     unit='(ug/mL)^-1 cm^-1', source='synthetic, for tests')


@pytest.fixture
def references():
    return Library([_reference('dsDNA', DNA_EPS),
                    _reference('protein', PROTEIN_EPS)], name='test')


def test_known_amounts_come_back(references):
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS, 'mixture')
    result = unmix(mixture, references)

    assert result['dsDNA'] == pytest.approx(0.03, rel=1e-6)
    assert result['protein'] == pytest.approx(0.40, rel=1e-6)
    assert result.r_squared > 0.9999
    assert np.max(np.abs(result.residual.y)) < 1e-9


def test_amounts_are_recovered_under_noise(references):
    rng = np.random.default_rng(0)
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS
                        + 0.005 * rng.normal(size=X.size))
    result = unmix(mixture, references)

    assert result['dsDNA'] == pytest.approx(0.03, abs=0.002)
    assert result['protein'] == pytest.approx(0.40, abs=0.01)
    assert np.all(result.stderr > 0)


def test_a_missing_component_shows_in_the_residual_not_in_r_squared(references):
    """
    The reason the residual is the documented diagnostic.

    A species the fit was never offered leaves R^2 looking excellent while
    putting a band-shaped feature in the residual, exactly where it absorbs.
    """
    contaminated = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS
                             + 0.02 * HAEM_EPS, 'contaminated')
    result = unmix(contaminated, references)

    assert result.r_squared > 0.98          # would pass any casual check
    peak = result.residual.x[int(np.argmax(np.abs(result.residual.y)))]
    assert peak == pytest.approx(300.0, abs=3.0)
    assert np.max(np.abs(result.residual.y)) > 0.5


def test_amounts_are_never_negative(references):
    """A negative concentration is not a small number, it is a wrong model."""
    protein_only = _spectrum(0.40 * PROTEIN_EPS - 0.05 * HAEM_EPS)
    result = unmix(protein_only, references)
    assert np.all(result.amounts >= 0.0)


def test_the_constraint_can_be_lifted_for_a_difference_spectrum(references):
    """There a negative coefficient means a component was lost, not gained."""
    difference = _spectrum(0.40 * PROTEIN_EPS - 0.02 * DNA_EPS)
    result = unmix(difference, references, non_negative=False)
    assert result['dsDNA'] == pytest.approx(-0.02, abs=1e-6)


def test_references_are_resampled_onto_the_sample_axis(references):
    """A reference from another instrument is never on the same grid."""
    coarse = np.linspace(225.0, 335.0, 60)
    mixture = Spectrum(coarse,
                       np.interp(coarse, X, 0.03 * DNA_EPS + 0.40 * PROTEIN_EPS),
                       x_unit='nm', y_unit='absorbance')
    result = unmix(mixture, references)
    assert result['dsDNA'] == pytest.approx(0.03, abs=1e-3)
    assert len(result.residual.x) == len(coarse)


def test_an_underdetermined_fit_is_refused(references):
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    with pytest.raises(ValueError, match='underdetermined'):
        unmix(mixture, references, wavelengths=[260.0])


def test_fitting_at_chosen_wavelengths_still_reconstructs_everywhere(references):
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    result = unmix(mixture, references, wavelengths=[230.0, 260.0, 280.0])
    assert result['dsDNA'] == pytest.approx(0.03, abs=1e-6)
    assert len(result.residual.x) == len(X)      # not just the three points


def test_best_wavelengths_choose_separability_over_peak_position(references):
    """
    Both components absorb at 280, so 280 is a poor place to tell them apart.
    Choosing by "where does each component peak" gets this wrong; choosing by
    the conditioning of the extinction matrix does not.
    """
    chosen = best_wavelengths(references)
    assert len(chosen) == 2
    # It should not pick the pair a naive reading would.
    assert not (np.any(np.isclose(chosen, 260.0, atol=2))
                and np.any(np.isclose(chosen, 280.0, atol=2)))
    # And the pair it picks must be better conditioned than 260/280.
    def conditioning(wavelengths):
        return np.linalg.cond(references.on(np.asarray(wavelengths)).T)
    assert conditioning(chosen) < conditioning([260.0, 280.0])


def test_too_few_wavelengths_for_the_components_is_refused(references):
    with pytest.raises(ValueError, match='cannot separate'):
        best_wavelengths(references, count=1)


# ---------------------------------------------------------------------------
# A260/A280
# ---------------------------------------------------------------------------

def test_the_ratio_is_computed_at_the_nearest_points():
    spectrum = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    at_260 = spectrum.y[int(np.abs(X - 260).argmin())]
    at_280 = spectrum.y[int(np.abs(X - 280).argmin())]
    assert absorbance_ratio(spectrum) == pytest.approx(at_260 / at_280)


def test_the_ratio_cannot_distinguish_what_unmixing_can():
    """
    The argument for preferring the whole spectrum, made as a test.

    A clean mixture and one contaminated with a third species can share an
    A260/A280 to within a few percent, while unmixing separates them plainly
    because it has a residual to show.
    """
    references = Library([_reference('dsDNA', DNA_EPS),
                          _reference('protein', PROTEIN_EPS)])
    clean = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    contaminated = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS
                             + 0.02 * HAEM_EPS)

    assert absorbance_ratio(clean) == pytest.approx(
        absorbance_ratio(contaminated), rel=0.05)

    clean_residual = np.max(np.abs(unmix(clean, references).residual.y))
    dirty_residual = np.max(np.abs(unmix(contaminated, references).residual.y))
    assert dirty_residual > 100 * max(clean_residual, 1e-12)


def test_nucleic_acid_and_protein_is_the_two_component_case():
    dna, protein = _reference('dsDNA', DNA_EPS), _reference('protein', PROTEIN_EPS)
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    result = nucleic_acid_and_protein(mixture, dna, protein)
    assert result.names == ('dsDNA', 'protein')
    assert result['dsDNA'] == pytest.approx(0.03, rel=1e-6)


# ---------------------------------------------------------------------------
# Library and calibration
# ---------------------------------------------------------------------------

def test_a_library_refuses_to_shadow_a_name(references):
    with pytest.raises(KeyError, match='already in this library'):
        references.add(_reference('dsDNA', DNA_EPS))


def test_an_unknown_reference_lists_what_there_is(references):
    with pytest.raises(KeyError, match='dsDNA, protein'):
        references['haem']


def test_select_keeps_the_order_asked_for(references):
    assert references.select('protein', 'dsDNA').names == ('protein', 'dsDNA')


def test_extinction_spectrum_is_recovered_from_a_concentration_series():
    """Beer-Lambert the other way round: known concentrations in, epsilon out."""
    rng = np.random.default_rng(0)
    concentrations = [2.0, 5.0, 10.0, 20.0, 40.0]
    series = SpectrumCollection(
        [_spectrum(c * DNA_EPS + 0.002 * rng.normal(size=X.size))
         for c in concentrations])

    reference = from_series(series, concentrations, 'dsDNA',
                            unit='(ug/mL)^-1 cm^-1')

    assert np.allclose(reference.spectrum.y, DNA_EPS, atol=1e-3)
    assert reference.is_absolute
    assert reference.uncertainty is not None
    # The stated uncertainty must actually describe the error.
    error = np.abs(reference.spectrum.y - DNA_EPS)
    assert np.mean(error < 3.0 * reference.uncertainty) > 0.9


def test_a_calibration_reference_can_be_used_to_unmix():
    """The loop closes: measure a standard, then use it on an unknown."""
    concentrations = [5.0, 10.0, 20.0]
    dna_series = SpectrumCollection([_spectrum(c * DNA_EPS)
                                     for c in concentrations])
    protein_series = SpectrumCollection([_spectrum(c * PROTEIN_EPS)
                                         for c in concentrations])
    measured = Library([
        from_series(dna_series, concentrations, 'dsDNA'),
        from_series(protein_series, concentrations, 'protein'),
    ])

    unknown = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    result = unmix(unknown, measured)
    assert result['dsDNA'] == pytest.approx(0.03, rel=1e-6)
    assert result['protein'] == pytest.approx(0.40, rel=1e-6)


def test_a_series_needs_matching_concentrations():
    series = SpectrumCollection([_spectrum(DNA_EPS), _spectrum(2 * DNA_EPS)])
    with pytest.raises(ValueError, match='one to one'):
        from_series(series, [1.0], 'dsDNA')


def test_a_single_concentration_is_refused():
    series = SpectrumCollection([_spectrum(DNA_EPS)])
    with pytest.raises(ValueError, match='at least two'):
        from_series(series, [1.0], 'dsDNA')


# ---------------------------------------------------------------------------
# Published coefficients
# ---------------------------------------------------------------------------

def test_published_coefficients_convert_absorbance_to_concentration():
    """An A260 of 1.0 is 50 ug/mL of double-stranded DNA."""
    assert lib.concentration_from_absorbance(1.0, 'dsDNA', 260) == \
        pytest.approx(50.0)
    assert lib.concentration_from_absorbance(1.0, 'ssDNA', 260) == \
        pytest.approx(33.0)
    assert lib.concentration_from_absorbance(0.5, 'RNA', 260) == \
        pytest.approx(20.0)


def test_path_length_is_applied():
    assert lib.coefficient('dsDNA', 260).concentration(1.0, path_length=0.1) == \
        pytest.approx(500.0)


def test_every_shipped_coefficient_says_where_it_came_from():
    """A number without a provenance is not usable in a paper."""
    for key, entry in lib.COEFFICIENTS.items():
        assert entry.source, f"{key} has no source"


def test_an_unknown_species_lists_the_known_ones():
    with pytest.raises(KeyError, match='dsDNA'):
        lib.coefficient('unobtainium', 260)


def test_protein_epsilon_from_sequence():
    """Pace et al.: 5500 per Trp, 1490 per Tyr, 125 per cystine."""
    assert lib.protein_epsilon_280(1, 0, 0) == pytest.approx(5500.0)
    assert lib.protein_epsilon_280(0, 1, 0) == pytest.approx(1490.0)
    # Lysozyme: 6 Trp, 3 Tyr, 4 cystines -> about 38 000, the accepted figure.
    assert lib.protein_epsilon_280(6, 3, 4) == pytest.approx(37970.0)


# ---------------------------------------------------------------------------
# Path length
# ---------------------------------------------------------------------------

def test_path_length_scales_the_concentrations(references):
    """A = eps*c*l, so the fit gives c*l and the concentration needs l."""
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)

    standard = unmix(mixture, references)
    short = unmix(mixture, references, path_length=0.1)

    assert standard.path_length == 1.0
    assert short['dsDNA'] == pytest.approx(10.0 * standard['dsDNA'])
    assert short.stderr[0] == pytest.approx(10.0 * standard.stderr[0])


def test_path_length_is_taken_from_metadata_when_not_given(references):
    """So a reader that knows the cuvette can record it once."""
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    mixture.metadata['path_length'] = 0.5

    result = unmix(mixture, references)
    assert result.path_length == 0.5
    assert result['dsDNA'] == pytest.approx(0.06, rel=1e-6)
    # An explicit argument still wins over the recorded value.
    assert unmix(mixture, references, path_length=1.0)['dsDNA'] == \
        pytest.approx(0.03, rel=1e-6)


def test_path_length_does_nothing_to_relative_references():
    """With no absolute reference there is no concentration to get wrong."""
    relative = Library([Reference('a', _spectrum(DNA_EPS)),
                        Reference('b', _spectrum(PROTEIN_EPS))])
    mixture = _spectrum(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    assert np.allclose(unmix(mixture, relative, path_length=0.1).amounts,
                       unmix(mixture, relative).amounts)


def test_a_nonsense_path_length_is_refused(references):
    mixture = _spectrum(DNA_EPS)
    with pytest.raises(ValueError, match='must be positive'):
        unmix(mixture, references, path_length=0.0)


def test_extinction_units_are_recognised_but_not_convertible():
    """
    Epsilon is absorbance-shaped, not another spelling of absorbance.

    Bands point up, so peak finding works; but eps -> A needs a concentration
    and a path length that are not properties of the spectrum, so it is
    deliberately absent from the convertible table.
    """
    from spectroscopy import units

    for unit in ('M^-1 cm^-1', 'mM^-1 cm^-1', '(ug/mL)^-1 cm^-1'):
        assert units.is_extinction(unit)
        assert units.band_direction(unit) == 'up'
        assert not units.can_convert_y(unit)

    # cm^-1 is wavenumber -- an x unit, and emphatically not an extinction.
    assert not units.is_extinction('cm^-1')
    assert not units.is_extinction('absorbance')
    assert not units.is_extinction('')
