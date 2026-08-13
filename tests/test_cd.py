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
from spectroscopy import library as lib
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
    """
    A structural basis, as a ReferenceSet: one spectrum per pure structure,
    each carrying the category it is of.

    The category is stored as a **name**, not as a ``Category`` object, and the
    objects are declared once on the set. That is what survives ``.spy``, which
    serialises metadata as JSON and would otherwise hand back a bare string
    with the DSSP states silently gone (ADR-0004 section 2.5).
    """
    spectra = []
    for category, y in SHAPES.items():
        spectrum = spc.Spectrum(X, y, technique='CD', name=category.name)
        spectrum.metadata['category'] = category.name
        spectrum.metadata['known_from'] = 'synthetic, built in this test'
        spectra.append(spectrum)
    return lib.ReferenceSet(spectra, name='synthetic basis',
                            info={'categories': list(SHAPES)})


@pytest.fixture
def proteins():
    """The other kind of standard: whole proteins of known composition."""
    truths = ({'helix': 0.80, 'sheet': 0.05, 'other': 0.15},
              {'helix': 0.10, 'sheet': 0.60, 'other': 0.30})
    spectra = [_mixture(truth, f'ref{index}')
               for index, truth in enumerate(truths)]
    compositions = [structure.Composition(
        fractions={c: truth[c.name] for c in SHAPES},
        method='dssp', technique='X-ray') for truth in truths]
    return lib.ReferenceSet.from_compositions(spectra, compositions,
                                              name='synthetic proteins')


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
    result = from_cd(_mixture(truth), 'basis-spectra', references=basis)

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
    normal = from_cd(_mixture(truth), 'basis-spectra', references=basis)
    scaled = _mixture(truth)
    scaled.y = scaled.y * 37.0
    louder = from_cd(scaled, 'basis-spectra', references=basis)

    for category in normal.fractions:
        assert louder.fractions[category] == pytest.approx(
            normal.fractions[category], abs=1e-6)


def test_fractions_sum_to_one_and_are_never_negative(basis):
    """A negative fraction of helix is not a small number, it is a wrong fit."""
    result = from_cd(_mixture({'helix': 0.9, 'sheet': 0.05, 'other': 0.05}),
                     'basis-spectra', references=basis)
    assert sum(result.fractions.values()) == pytest.approx(1.0, abs=1e-6)
    assert all(value >= 0.0 for value in result.fractions.values())


def test_reference_proteins_carry_their_own_compositions(proteins):
    """The other kind of standard: fit proteins, then mix their structures."""
    first = {'helix': 0.80, 'sheet': 0.05, 'other': 0.15}
    second = {'helix': 0.10, 'sheet': 0.60, 'other': 0.30}

    unknown = spc.Spectrum(X, 0.25 * proteins[0].y + 0.75 * proteins[1].y,
                           technique='CD', name='unknown')
    result = from_cd(unknown, 'reference-proteins', references=proteins)

    for category in SHAPES:
        expected = 0.25 * first[category.name] + 0.75 * second[category.name]
        assert result.fractions[category] == pytest.approx(expected, abs=1e-3)


def test_the_method_must_match_the_kind_of_standard(basis, proteins):
    """
    A structural basis and a set of reference proteins give different answers
    from the same spectrum, so the method names which was used -- and naming
    the wrong one is now caught rather than quietly answered.
    """
    sample = _mixture({'helix': 0.55, 'sheet': 0.20, 'other': 0.25})
    with pytest.raises(ValueError, match='does not match the standards'):
        from_cd(sample, 'reference-proteins', references=basis)
    with pytest.raises(ValueError, match='does not match the standards'):
        from_cd(sample, 'basis-spectra', references=proteins)


def test_the_method_must_be_named(basis):
    """The two kinds of standard answer differently; the result must say."""
    with pytest.raises(ValueError, match='method must be one of'):
        from_cd(_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}),
                references=basis)


def test_no_basis_says_where_to_get_one():
    with pytest.raises(ValueError, match='no reference spectra ship'):
        from_cd(_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}),
                'basis-spectra')


def test_a_reference_must_say_what_it_is_of():
    """A set of spectra with no known structure cannot answer the question."""
    plain = spc.Spectrum(X, SHAPES[HELIX], technique='CD', name='helix')
    with pytest.raises(ValueError, match='no known structure'):
        from_cd(_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}),
                'basis-spectra', references=lib.ReferenceSet([plain]))


def test_two_parallel_lists_are_refused_at_the_door(basis):
    """
    ADR-0004's whole point. Spectra and their compositions as two arguments
    could fall out of step, and the case that did not raise -- same length,
    wrong order -- moved a real helix estimate from 0.464 to 0.307 in silence.
    The pairing now happens once, in from_compositions, where it is named.
    """
    sample = _mixture({'helix': 0.55, 'sheet': 0.20, 'other': 0.25})
    with pytest.raises(TypeError, match='must be a library.ReferenceSet'):
        from_cd(sample, 'basis-spectra', references=list(basis))


def test_a_reversed_pairing_can_no_longer_be_expressed(basis):
    """
    There is one list, so reversing the spectra reverses their structures with
    them. The composition of the reversed set is the same set of answers in a
    different order -- not a different set of answers.
    """
    forwards = {s.name: c.get('helix')
                for s, c in zip(basis, basis.compositions)}
    backwards = {s.name: c.get('helix')
                 for s, c in zip(basis[::-1], basis[::-1].compositions)}
    assert forwards == backwards


def test_a_basis_that_cannot_describe_the_spectrum_shows_in_the_rmsd(basis):
    """
    The number to look at before believing any fraction. A fit against a basis
    missing the dominant component still returns fractions summing to one.
    """
    good = from_cd(_mixture({'helix': 0.55, 'sheet': 0.20, 'other': 0.25}),
                   'basis-spectra', references=basis)
    unusual = spc.Spectrum(X, _band(X, 230, 6, 40) - _band(X, 205, 5, 25),
                           technique='CD', name='not in the basis')
    poor = from_cd(unusual, 'basis-spectra', references=basis)

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


# ---------------------------------------------------------------------------
# mean residue ellipticity, and the single-wavelength helix estimate
# ---------------------------------------------------------------------------

def test_mre_conversion_matches_the_worked_example():
    """
    AqpZ-W14A: -65.21 mdeg at 222 nm, 20 uM, 0.5 mm cell, 243 residues
    (231 + a 12-residue N-terminal extension). Checked by hand so the factors
    of ten in the formula are pinned rather than trusted.
    """
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-65.21]),
                            technique='CD')
    mre = spectrum.to_mean_residue_ellipticity(
        concentration=20e-6, path_length=0.05, residues=243)

    assert mre.y_unit == 'deg cm^2 dmol^-1'
    assert mre.y[0] == pytest.approx(-65.21 / (10 * 0.05 * 20e-6 * 243))
    assert mre.y[0] == pytest.approx(-26836, abs=1)


def test_the_three_sample_facts_are_never_defaulted():
    """Each scales the answer linearly and leaves it looking like a protein."""
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-65.21]),
                            technique='CD')
    with pytest.raises(ValueError, match='path_length'):
        spectrum.to_mean_residue_ellipticity(concentration=20e-6, residues=243)
    with pytest.raises(ValueError, match='concentration'):
        spectrum.to_mean_residue_ellipticity(path_length=0.05, residues=243)


def test_the_conversion_reads_the_agreed_metadata_keys():
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-65.21]),
                            technique='CD')
    spectrum.metadata.update({'concentration': 20e-6, 'path_length': 0.05,
                              'n_residues': 243})
    assert spectrum.to_mean_residue_ellipticity().y[0] == pytest.approx(
        -26836, abs=1)


def test_theta222_refuses_millidegrees():
    """A guessed path length rescales helix without changing how it looks."""
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-65.21]),
                            technique='CD')
    with pytest.raises(ValueError, match='mean residue ellipticity'):
        structure.helix_from_theta222(spectrum)


def test_theta222_fills_one_category_and_leaves_the_rest_none():
    """
    ADR-0002: it is not a decomposition. Every other category must be absent
    rather than zero -- zero would be a claim about sheet that one wavelength
    cannot support.
    """
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-26836.0]),
                            technique='CD')
    spectrum.y_unit = 'deg cm^2 dmol^-1'
    result = structure.helix_from_theta222(spectrum, residues=243)

    assert result.method == 'theta-222'
    assert list(result.fractions) == [Category('helix', frozenset({'G', 'H', 'I'}))]
    assert next(iter(result.fractions.values())) == pytest.approx(0.687, abs=0.005)
    assert result.quality['single_wavelength'] is True


def test_the_chain_length_correction_is_applied_and_recorded():
    """A helix has two ends that make no hydrogen bonds; short chains signal
    less per residue, and ignoring it overestimates helix."""
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-26836.0]),
                            technique='CD')
    spectrum.y_unit = 'deg cm^2 dmol^-1'
    short = structure.helix_from_theta222(spectrum, residues=20)
    long = structure.helix_from_theta222(spectrum, residues=1000)

    assert next(iter(short.fractions.values())) > next(iter(long.fractions.values()))
    assert short.quality['chain_length_corrected']


def test_an_impossible_helix_fraction_warns():
    """Out of 0-1 means the conversion inputs were wrong, not the protein."""
    spectrum = spc.Spectrum(np.array([222.0]), np.array([-67090.0]),
                            technique='CD')
    spectrum.y_unit = 'deg cm^2 dmol^-1'
    with pytest.warns(UserWarning, match='outside 0-1'):
        structure.helix_from_theta222(spectrum, residues=243)


# ---------------------------------------------------------------------------
# shape without amplitude
# ---------------------------------------------------------------------------

def _helical_shape(x):
    return (_band(x, 208, 7, -37) + _band(x, 222, 9, -37)
            + _band(x, 193, 7, 75))


def test_every_descriptor_is_blind_to_amplitude():
    """
    The point of the module. Concentration, path length, residue count and a
    pipetting slip all multiply a spectrum by a constant; none of them should
    move a single number here.
    """
    x = np.linspace(190.0, 250.0, 241)
    quiet = spc.Spectrum(x, _helical_shape(x), technique='CD')
    loud = spc.Spectrum(x, 17.3 * _helical_shape(x), technique='CD')

    first = structure.cd_shape_descriptors(quiet)
    second = structure.cd_shape_descriptors(loud)
    for key, value in first.items():
        if isinstance(value, float) and np.isfinite(value):
            assert second[key] == pytest.approx(value, rel=1e-9), key


def test_the_zero_crossing_separates_helix_from_coil():
    """
    A position rather than a size, so it carries structural information with
    no amplitude in it. Helix crosses near 200 nm on the way up from its
    positive 193 band; a coil, which has no positive band there, crosses much
    further red on the way up to its weak positive near 218.
    """
    x = np.linspace(190.0, 250.0, 241)
    helix = structure.cd_shape_descriptors(
        spc.Spectrum(x, _helical_shape(x), technique='CD'))['zero_crossing']
    coil = structure.cd_shape_descriptors(spc.Spectrum(
        x, _band(x, 198, 7, -40) + _band(x, 218, 11, 4),
        technique='CD'))['zero_crossing']

    assert helix == pytest.approx(202.0, abs=2.0)
    assert coil > helix + 5.0


def test_a_minimum_at_the_edge_is_called_out():
    """
    The real AqpZ case: cropped at the detector's limit, the most negative
    point was the crop itself. Reading it as a band position reads the
    instrument's failure as a property of the protein.
    """
    x = np.linspace(210.0, 250.0, 161)          # cut off above the band
    spectrum = spc.Spectrum(x, _band(x, 205, 8, -40), technique='CD')
    with pytest.warns(UserWarning, match='where the data stops'):
        result = structure.cd_shape_descriptors(spectrum)
    assert result['minimum_at_edge']
    assert result['minimum'] == pytest.approx(210.0, abs=0.5)


def test_a_real_minimum_is_not_called_an_edge():
    x = np.linspace(190.0, 250.0, 241)
    result = structure.cd_shape_descriptors(
        spc.Spectrum(x, _helical_shape(x), technique='CD'))
    assert not result['minimum_at_edge']


def test_from_cd_needs_no_amplitude_either():
    """
    Stated as a test because it is the reason a composition can be had from a
    sample of unknown concentration: the fit is on shape alone.
    """
    x = np.linspace(190.0, 250.0, 241)
    spectra, categories = [], []
    for name, y in (('helix', _helical_shape(x)),
                    ('other', _band(x, 198, 7, -40) + _band(x, 220, 10, 3))):
        reference = spc.Spectrum(x, y, technique='CD', name=name)
        reference.metadata['category'] = name
        spectra.append(reference)
        categories.append(Category(
            name, frozenset({'H'} if name == 'helix' else {'-'})))
    basis = lib.ReferenceSet(spectra, info={'categories': categories})

    truth = 0.7 * basis[0].y + 0.3 * basis[1].y
    quiet = from_cd(spc.Spectrum(x, truth, technique='CD'),
                    'basis-spectra', references=basis)
    loud = from_cd(spc.Spectrum(x, 250.0 * truth, technique='CD'),
                   'basis-spectra', references=basis)

    for category in quiet.fractions:
        assert loud.fractions[category] == pytest.approx(
            quiet.fractions[category], abs=1e-6)


# ---------------------------------------------------------------------------
# loading a basis somebody else measured
# ---------------------------------------------------------------------------

def _write_basis_files(tmp_path, shapes):
    for name, y in shapes.items():
        (tmp_path / f'{name}.csv').write_text(
            'wavelength,cd\n' + '\n'.join(f'{a},{b}' for a, b in zip(X, y)))


def test_a_structural_basis_loads_from_a_manifest(tmp_path):
    _write_basis_files(tmp_path, {c.name: y for c, y in SHAPES.items()})
    (tmp_path / 'basis.csv').write_text(
        'file,category,source,citation\n'
        'helix.csv,helix,measured here,doi:10.0/x\n'
        'sheet.csv,sheet,measured here,doi:10.0/x\n'
        'other.csv,other,measured here,doi:10.0/x\n')

    basis = lib.load_basis(tmp_path / 'basis.csv')

    assert isinstance(basis, lib.ReferenceSet)
    assert basis.is_structural
    assert len(basis) == 3
    assert {s.metadata['category'] for s in basis} == {
        'helix', 'sheet', 'other'}
    assert basis[0].metadata['reference_citation'] == 'doi:10.0/x'
    # A value every row repeats describes the set, and a citation condition
    # needs somewhere to live as an obligation of the whole set.
    assert basis.info['citation'] == 'doi:10.0/x'
    assert basis.info['source'] == 'measured here'
    # And each spectrum's one category reads back as a one-hot composition,
    # so nothing downstream has to branch on which kind of set this is.
    helix = basis.compositions[0]
    assert helix.get('helix') == 1.0 and helix.get('sheet') == 0.0


def test_set_level_provenance_survives_selection(tmp_path):
    """
    A subset of SP175 is still SP175 data and still carries its citation
    condition. ``info`` travels with ``select``, ``crop`` and slicing exactly
    as ``name`` does -- and the subset is still a ReferenceSet, not a plain
    collection that has quietly lost its provenance.
    """
    _write_basis_files(tmp_path, {c.name: y for c, y in SHAPES.items()})
    (tmp_path / 'basis.csv').write_text(
        'file,category,citation\n'
        'helix.csv,helix,doi:10.0/x\n'
        'sheet.csv,sheet,doi:10.0/x\n'
        'other.csv,other,doi:10.0/x\n')
    basis = lib.load_basis(tmp_path / 'basis.csv')

    subset = basis.select(lambda s: s.name != 'sheet')
    assert isinstance(subset, lib.ReferenceSet)
    assert len(subset) == 2
    assert subset.info['citation'] == 'doi:10.0/x'
    assert [c.get('helix') for c in subset.compositions] == [1.0, 0.0]


def test_a_loaded_basis_goes_straight_into_from_cd(tmp_path):
    """The whole point: obtain a basis, load it, get a composition."""
    _write_basis_files(tmp_path, {c.name: y for c, y in SHAPES.items()})
    (tmp_path / 'basis.csv').write_text(
        'file,category\nhelix.csv,helix\nsheet.csv,sheet\nother.csv,other\n')
    basis = lib.load_basis(tmp_path / 'basis.csv')

    truth = {'helix': 0.55, 'sheet': 0.20, 'other': 0.25}
    result = from_cd(_mixture(truth), 'basis-spectra', references=basis)
    for category, fraction in result.fractions.items():
        assert fraction == pytest.approx(truth[category.name], abs=2e-3)


def test_reference_proteins_load_with_their_compositions(tmp_path):
    first = {'helix': 0.80, 'sheet': 0.05, 'other': 0.15}
    second = {'helix': 0.10, 'sheet': 0.60, 'other': 0.30}
    _write_basis_files(tmp_path, {
        'p1': sum(first[c.name] * y for c, y in SHAPES.items()),
        'p2': sum(second[c.name] * y for c, y in SHAPES.items())})
    (tmp_path / 'refs.csv').write_text(
        'file,name,helix,sheet,other\n'
        'p1.csv,protein one,0.80,0.05,0.15\n'
        'p2.csv,protein two,0.10,0.60,0.30\n')

    basis = lib.load_basis(tmp_path / 'refs.csv')
    assert not basis.is_structural
    assert len(basis.compositions) == 2
    assert basis[0].name == 'protein one'

    unknown = spc.Spectrum(X, 0.25 * basis[0].y + 0.75 * basis[1].y,
                           technique='CD')
    result = from_cd(unknown, 'reference-proteins', references=basis)
    for category in result.fractions:
        expected = (0.25 * first[category.name] + 0.75 * second[category.name])
        assert result.fractions[category] == pytest.approx(expected, abs=3e-3)


def test_a_manifest_must_choose_one_kind_of_basis(tmp_path):
    _write_basis_files(tmp_path, {'helix': SHAPES[HELIX]})
    (tmp_path / 'both.csv').write_text(
        'file,category,helix\nhelix.csv,helix,1.0\n')
    with pytest.raises(ValueError, match='not both and not neither'):
        lib.load_basis(tmp_path / 'both.csv')


def test_percentages_instead_of_fractions_are_caught(tmp_path):
    """80/5/15 sums to 100, not 1 -- a composition it is not."""
    _write_basis_files(tmp_path, {'p1': SHAPES[HELIX]})
    (tmp_path / 'refs.csv').write_text(
        'file,helix,sheet,other\np1.csv,80,5,15\n')
    with pytest.raises(ValueError, match='not a composition'):
        lib.load_basis(tmp_path / 'refs.csv')


def test_a_missing_reference_file_says_why_it_is_missing(tmp_path):
    (tmp_path / 'basis.csv').write_text('file,category\nabsent.csv,helix\n')
    with pytest.raises(FileNotFoundError, match='not shipped with this package'):
        lib.load_basis(tmp_path / 'basis.csv')


def test_convergence_at_a_dead_end_is_not_an_isodichroic_point():
    """
    Found on the real AqpZ melt, which reported a "tight crossing" at 238.8 nm
    -- the red end, where every spectrum has decayed towards zero and so they
    trivially agree. Spectra agreeing because there is no signal is the
    absence of information, not evidence of two states.
    """
    x = np.linspace(200.0, 260.0, 121)
    decay = np.exp(-((x - 205.0) ** 2) / (2 * 12.0 ** 2))
    spectra = [spc.Spectrum(x, -(40.0 - 3.0 * step) * decay, technique='CD',
                            name=f'{step}')
               for step in range(8)]
    series = spc.SpectrumCollection(spectra).with_parameters(
        np.arange(30.0, 30.0 + 8 * 5, 5.0), name='temperature', unit='C')

    result = melting.isodichroic_point(series)
    # The spectra converge towards zero at the red end but never cross: every
    # one stays the same side of the others, so there is no inversion.
    assert not result['is_tight']
    assert not result['inverts_here'] or result['wavelength'] < 240.0


def test_a_real_crossing_is_still_found():
    x = np.linspace(200.0, 260.0, 121)
    folded = _band(x, 208, 9, -40) + _band(x, 235, 9, 12)
    unfolded = _band(x, 202, 8, -30) + _band(x, 235, 9, -6)
    fractions = np.linspace(0.0, 1.0, 9)
    spectra = [spc.Spectrum(x, (1 - f) * folded + f * unfolded,
                            technique='CD', name=f'{f}') for f in fractions]
    series = spc.SpectrumCollection(spectra).with_parameters(
        np.linspace(30.0, 90.0, 9), name='temperature', unit='C')

    result = melting.isodichroic_point(series)
    assert result['is_tight']
    assert result['inverts_here']


def test_seven_points_cannot_determine_six_parameters():
    """
    The real 4 uM AqpZ series has seven temperatures. It fitted, and returned
    a Tm with a standard error of 1e17 -- the fit saying it has no idea while
    still printing a number.
    """
    with pytest.raises(ValueError, match='Below ten'):
        melting.two_state(np.linspace(30.0, 90.0, 7), np.linspace(1.0, 0.0, 7))


def test_a_degenerate_fit_says_so_rather_than_quoting_a_number():
    generator = np.random.default_rng(3)
    temperature = np.linspace(30.0, 90.0, 12)
    noise = 1.0 + 0.01 * generator.normal(size=temperature.size)
    with pytest.warns(UserWarning):
        melting.two_state(temperature, noise)


def test_nearest_references_is_amplitude_free(basis):
    """Shape only: both spectra are unit-normalised before comparison."""
    target = _mixture({'helix': 0.9, 'sheet': 0.05, 'other': 0.05})
    quiet = structure.nearest_references(target, basis, count=2,
                                         region=(190.0, 250.0))
    louder = spc.Spectrum(X, 31.0 * target.y, technique='CD')
    loud = structure.nearest_references(louder, basis, count=2,
                                        region=(190.0, 250.0))
    assert ([n[1] for n in quiet['neighbours']]
            == [n[1] for n in loud['neighbours']])
    assert quiet['neighbours'][0][0] == pytest.approx(
        loud['neighbours'][0][0], rel=1e-9)


def test_nearest_references_finds_the_right_shape(basis):
    helical = _mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0})
    result = structure.nearest_references(helical, basis, count=1,
                                          region=(190.0, 250.0))
    assert result['neighbours'][0][1] == 'helix'


def test_more_references_than_the_data_can_determine_is_refused():
    """
    Rows are not information. Resampling a 1 nm basis onto a 0.1 nm spectrum
    made 439 equations of rank 46 look like a well-posed fit for 71 unknowns;
    it returned a composition 34 points from the crystal at an rmsd of 0.6%.
    """
    x = np.linspace(200.0, 240.0, 401)                # 0.1 nm sample grid
    coarse = np.arange(200.0, 240.5, 1.0)             # 1 nm references
    spectra = []
    for index in range(60):
        y = _band(coarse, 205 + 0.4 * index, 8, -30)
        reference = spc.Spectrum(coarse, y, technique='CD', name=f'r{index}')
        reference.metadata['category'] = 'helix'
        spectra.append(reference)
    basis = lib.ReferenceSet(
        spectra, info={'categories': [Category('helix', frozenset({'H'}))]})

    sample = spc.Spectrum(x, _band(x, 210, 8, -30), technique='CD')
    with pytest.raises(ValueError, match='cannot be determined|cannot determine'):
        from_cd(sample, 'basis-spectra', references=basis,
                region=(200.0, 240.0))


# ---------------------------------------------------------------------------
# the method registry and held-out validation
# ---------------------------------------------------------------------------

from spectroscopy.processing import cd as cdm  # noqa: E402


def _synthetic_set(count=24, seed=0, noise=0.3):
    x = np.linspace(190.0, 240.0, 51)
    pure = {HELIX: _band(x, 208, 7, -37) + _band(x, 222, 9, -37)
                   + _band(x, 193, 7, 60),
            SHEET: _band(x, 217, 9, -18) + _band(x, 195, 8, 32),
            COIL: _band(x, 198, 7, -40) + _band(x, 220, 10, 3)}
    generator = np.random.default_rng(seed)
    spectra, compositions = [], []
    for index in range(count):
        weights = generator.dirichlet([1.4, 1.0, 1.0])
        y = sum(w * pure[c] for w, c in zip(weights, pure))
        spectra.append(spc.Spectrum(x, y + noise * generator.normal(size=x.size),
                                    technique='CD', name=f'ref{index}'))
        compositions.append(structure.Composition(
            fractions=dict(zip(pure, weights)), method='known',
            technique='X-ray'))
    return lib.ReferenceSet.from_compositions(spectra, compositions,
                                              name='synthetic set')


@pytest.mark.parametrize('method', sorted(cdm.METHODS))
def test_every_method_recovers_a_known_mixture(method):
    # Three references for all-references, which correctly refuses anything
    # the data cannot determine -- these spectra span a three-dimensional
    # space, so three is exactly what it can take.
    count = 3 if method == 'all-references' else 12
    references = _synthetic_set(count=count, noise=0.0)
    truth = {HELIX: 0.6, SHEET: 0.25, COIL: 0.15}
    x = references[0].x
    pure = {HELIX: _band(x, 208, 7, -37) + _band(x, 222, 9, -37)
                   + _band(x, 193, 7, 60),
            SHEET: _band(x, 217, 9, -18) + _band(x, 195, 8, 32),
            COIL: _band(x, 198, 7, -40) + _band(x, 220, 10, 3)}
    sample = spc.Spectrum(x, sum(truth[c] * pure[c] for c in pure),
                          technique='CD', name='unknown')

    result = cdm.estimate(sample, method, references,
                          region=(190.0, 240.0),
                          **({'draws': 300} if method == 'subset-average' else {}))
    assert result.method == method
    # Loose: this is asking that each method is in the right area, not that it
    # is accurate. Accuracy is what benchmark() measures.
    assert result.fractions[HELIX] > result.fractions[SHEET]


@pytest.mark.parametrize('method', ['nearest-shapes', 'ridge'])
def test_the_shape_only_methods_ignore_amplitude(method):
    """
    The constraint that decides which methods can be used at all: a reference
    set is in delta epsilon and a measurement is in millidegrees until three
    sample facts are supplied.
    """
    references = _synthetic_set(count=12, noise=0.0)
    x = references[0].x
    y = _band(x, 208, 7, -30) + _band(x, 222, 9, -30)
    quiet = cdm.estimate(spc.Spectrum(x, y, technique='CD'), method, references)
    loud = cdm.estimate(spc.Spectrum(x, 47.0 * y, technique='CD'), method,
                        references)
    for category in quiet.fractions:
        assert loud.fractions[category] == pytest.approx(
            quiet.fractions[category], abs=1e-6)


def test_sum_to_one_needs_matched_units_and_says_so():
    """
    The classical self-consistency test compares an amplitude. On a spectrum
    in millidegrees against a basis in delta epsilon, nothing is accepted.
    """
    references = _synthetic_set(count=12, noise=0.0)
    x = references[0].x
    mismatched = spc.Spectrum(x, 900.0 * (_band(x, 208, 7, -30)
                                          + _band(x, 222, 9, -30)),
                              technique='CD')
    with pytest.raises(ValueError, match='sum to about 1'):
        cdm.estimate(mismatched, 'subset-average', references,
                     draws=300, sum_to_one=True)


def test_all_references_refuses_when_it_cannot_determine_them():
    references = _synthetic_set(count=80, noise=0.1)
    held_out = references[0]
    with pytest.raises(ValueError, match='cannot be determined'):
        cdm.estimate(held_out, 'all-references', references[1:])


def test_benchmark_reports_bias_as_well_as_error():
    """
    Bias is the number that catches a method pulling towards the reference
    set's mean, which is the failure that grows with how unusual a protein is.
    """
    references = _synthetic_set(count=18, seed=1)
    scores = cdm.benchmark(references,
                           ['nearest-shapes', 'ridge'], folds=3,
                           options={'subset-average': {'draws': 100}})
    for method in ('nearest-shapes', 'ridge'):
        entry = scores[method]
        assert entry['n'] == 18 and entry['failures'] == 0
        for category in ('helix', 'sheet', 'other'):
            assert 0.0 <= entry[category]['rmse'] < 1.0
            assert abs(entry[category]['bias']) <= entry[category]['rmse'] + 1e-9


def test_benchmark_counts_failures_rather_than_hiding_them():
    # Noiseless, so the twelve references span only the three shapes they were
    # built from: eleven references, rank three, every fold refused.
    references = _synthetic_set(count=12, noise=0.0)
    scores = cdm.benchmark(references, ['all-references'], folds=3)
    # 11 references cannot be determined by this data, so every fold fails.
    assert scores['all-references']['failures'] == 12
    assert scores['all-references']['n'] == 0


def test_the_design_matrix_uses_the_coarser_grid():
    fine = np.linspace(200.0, 240.0, 401)
    coarse = np.arange(200.0, 240.5, 1.0)
    reference = spc.Spectrum(coarse, _band(coarse, 220, 8, -20), technique='CD')
    sample = spc.Spectrum(fine, _band(fine, 220, 8, -20), technique='CD')
    grid, measured, design = cdm.design_matrix(sample, [reference],
                                               region=(200.0, 240.0))
    assert len(grid) == pytest.approx(41, abs=1)
    assert design.shape == (len(grid), 1)


# ---------------------------------------------------------------------------
# SELCON, the class signature, and measured noise
# ---------------------------------------------------------------------------

def test_selcon_is_in_the_registry_and_runs():
    references = _synthetic_set(count=20, seed=2, noise=0.2)
    x = references[0].x
    pure = {HELIX: _band(x, 208, 7, -37) + _band(x, 222, 9, -37)
                   + _band(x, 193, 7, 60),
            SHEET: _band(x, 217, 9, -18) + _band(x, 195, 8, 32),
            COIL: _band(x, 198, 7, -40) + _band(x, 220, 10, 3)}
    truth = {HELIX: 0.55, SHEET: 0.25, COIL: 0.20}
    sample = spc.Spectrum(x, sum(truth[c] * pure[c] for c in pure),
                          technique='CD', name='unknown')

    result = cdm.estimate(sample, 'selcon', references,
                          ignore_sum_rule=True, max_references=12)
    assert result.method == 'selcon'
    assert result.quality['n_solutions'] > 0
    assert sum(result.fractions.values()) == pytest.approx(1.0, abs=1e-6)


def test_selcon_needs_the_amplitude_and_says_so():
    """
    The self-consistent step puts the query into the basis beside the
    references, so its magnitude is part of the model. A wrong-scale spectrum
    is refused at some amplitudes -- and, worse, quietly answered at others,
    which is why the docstring quotes the measured drift instead of promising
    scale-freedom.
    """
    references = _synthetic_set(count=16, seed=3, noise=0.1)
    x = references[0].x
    y = _band(x, 208, 7, -30) + _band(x, 222, 9, -30)

    with pytest.raises(ValueError, match='amplitude'):
        cdm.estimate(spc.Spectrum(x, 0.1 * y, technique='CD'), 'selcon',
                     references, max_references=10)

    with pytest.warns(UserWarning, match='not scale-free'):
        relaxed = cdm.estimate(spc.Spectrum(x, 0.1 * y, technique='CD'),
                               'selcon', references,
                               ignore_sum_rule=True, max_references=10)
    assert relaxed.quality['ignore_sum_rule'] is True


def test_ignore_sum_rule_does_not_make_selcon_scale_free():
    """
    The flag was originally called scale_free. It is not: renormalising the
    fractions does not remove the query's amplitude from the model.
    """
    references = _synthetic_set(count=16, seed=3, noise=0.1)
    x = references[0].x
    y = _band(x, 208, 7, -30) + _band(x, 222, 9, -30)
    with pytest.warns(UserWarning):
        quiet = cdm.estimate(spc.Spectrum(x, 0.1 * y, technique='CD'),
                             'selcon', references,
                             ignore_sum_rule=True, max_references=10)
        loud = cdm.estimate(spc.Spectrum(x, y, technique='CD'), 'selcon',
                            references, ignore_sum_rule=True,
                            max_references=10)
    assert quiet.fractions[HELIX] != pytest.approx(loud.fractions[HELIX],
                                                   abs=0.01)


def test_measured_noise_weights_the_fit():
    """
    Real CD noise is strongly wavelength-dependent -- on the AqpZ scans the
    standard error at 215 nm is three times that at 240. An unweighted fit
    treats both alike.
    """
    references = _synthetic_set(count=16, seed=4, noise=0.1)
    x = references[0].x
    y = _band(x, 208, 7, -30) + _band(x, 222, 9, -30)
    sample = spc.Spectrum(x, y, technique='CD')

    # Noise huge at the blue end, tiny at the red.
    sigma = np.linspace(20.0, 0.05, x.size)
    plain = cdm.estimate(sample, 'ridge', references)
    weighted = cdm.estimate(sample, 'ridge', references, sigma=sigma)

    assert plain.quality['weighted'] is False
    assert weighted.quality['weighted'] is True
    assert any(weighted.fractions[c] != plain.fractions[c]
               for c in plain.fractions)


def test_resampling_reports_how_far_the_answer_moves_under_its_own_noise():
    references = _synthetic_set(count=16, seed=5, noise=0.1)
    x = references[0].x
    sample = spc.Spectrum(x, _band(x, 208, 7, -30) + _band(x, 222, 9, -30),
                          technique='CD')
    result = cdm.estimate(sample, 'ridge', references,
                          sigma=np.full(x.size, 1.5), resamples=25)
    assert result.quality['n_resamples'] > 0
    assert set(result.quality['uncertainty']) == {'helix', 'sheet', 'other'}
    assert all(v >= 0.0 for v in result.quality['uncertainty'].values())


def test_uncertainty_from_replicates_is_the_standard_error():
    x = np.linspace(200.0, 240.0, 41)
    generator = np.random.default_rng(0)
    scans = [spc.Spectrum(x, _band(x, 222, 9, -30)
                          + 0.5 * generator.normal(size=x.size), technique='CD')
             for _ in range(5)]
    sem = cdm.uncertainty_from_replicates(spc.SpectrumCollection(scans))
    assert len(sem) == len(x)
    assert 0.05 < float(np.median(sem.y)) < 1.0


def test_known_truth_survives_a_spy_round_trip(tmp_path, basis):
    """
    ADR-0004 section 2.5. ``.spy`` serialises metadata as JSON and silently
    degrades anything else -- a ``Category`` written into metadata came back a
    bare ``str``, with the frozenset of DSSP states gone and nothing saying so.

    So the stored form is JSON-native: a category is its name, a composition is
    ``{name: fraction}``, and the objects are rebuilt from the set's own
    declaration. This is what makes the round trip lossless rather than
    nearly-lossless.
    """
    for index, spectrum in enumerate(basis):
        spectrum.save_as(str(tmp_path / f'ref{index}.spy'))
    reloaded = lib.ReferenceSet(
        [spc.read(tmp_path / f'ref{index}.spy') for index in range(len(basis))],
        info={'categories': list(SHAPES)})

    for before, after in zip(basis.compositions, reloaded.compositions):
        assert {c.name: f for c, f in before.fractions.items()} == \
               {c.name: f for c, f in after.fractions.items()}
    # And the states came back, which is the part JSON cannot carry by itself.
    assert all(c.states for c in reloaded.compositions[0].fractions)


def test_a_category_the_set_does_not_declare_is_refused():
    """
    A bare name does not say which DSSP states it covers, so a set cannot
    silently invent the vocabulary its numbers are in.
    """
    spectrum = spc.Spectrum(X, SHAPES[HELIX], technique='CD', name='r')
    spectrum.metadata['composition'] = {'helix': 0.8, 'polyproline': 0.2}
    references = lib.ReferenceSet([spectrum], info={'categories': [HELIX]})
    with pytest.raises(ValueError, match='which this set does not declare'):
        _ = references.compositions


def test_a_set_cannot_hold_two_meanings_of_one_category():
    """Two Category objects of one name and different states are two claims."""
    strict = Category('helix', frozenset({'H'}))
    loose = Category('helix', frozenset({'G', 'H', 'I'}))
    spectra = [_mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}, 'a'),
               _mixture({'helix': 0.5, 'sheet': 0.5, 'other': 0.0}, 'b')]
    compositions = [
        structure.Composition(fractions={strict: 1.0}, method='dssp',
                              technique='X-ray'),
        structure.Composition(fractions={loose: 0.5}, method='dssp',
                              technique='X-ray'),
    ]
    with pytest.raises(ValueError, match='disagree about what'):
        lib.ReferenceSet.from_compositions(spectra, compositions)


def test_from_compositions_checks_the_counts_and_says_why(basis):
    with pytest.raises(ValueError, match='matched by position'):
        lib.ReferenceSet.from_compositions(list(basis), basis.compositions[:2])


def test_from_compositions_does_not_touch_the_callers_spectra():
    spectrum = _mixture({'helix': 1.0, 'sheet': 0.0, 'other': 0.0}, 'r')
    composition = structure.Composition(fractions={HELIX: 1.0}, method='dssp',
                                        technique='X-ray')
    lib.ReferenceSet.from_compositions([spectrum], [composition])
    assert 'composition' not in spectrum.metadata


def test_a_result_records_which_reference_set_it_came_from(basis):
    """
    The same spectrum against SP175 and against SMP180 is two results, not one,
    and a composition that cannot say which is not reproducible.
    """
    basis.info['source'] = 'DichroWebGit SMP180'
    result = cdm.estimate(_mixture({'helix': 0.6, 'sheet': 0.2, 'other': 0.2}),
                          'ridge', basis, region=(190.0, 250.0))
    assert result.quality['references'] == 'DichroWebGit SMP180'


# ---------------------------------------------------------------------------
# the reference sets that ship
#
# SP175 and SMP180 are somebody else's data, and they ship because their terms
# allow it -- MIT, via github.com/pcddb/DichroWebGit -- rather than because
# they are useful. ADR-0002 section 9 forbade shipping reference data until the
# terms had been checked; they have been, and these tests are what keep the
# conditions attached to the data rather than to a memory of a conversation.
# ---------------------------------------------------------------------------

def test_the_reference_sets_ship_and_load():
    for name, count in (('sp175', 71), ('smp180', 128)):
        references = spc.datasets.reference_set(name)
        assert isinstance(references, lib.ReferenceSet)
        assert len(references) == count
        assert references.has_truth and not references.is_structural


def test_the_licence_notice_ships_beside_the_data():
    """
    MIT permits redistribution **provided the notice travels with the copy**.
    Shipping the numbers without this file is the one way this arrangement
    becomes a licence violation, and it is a file rather than code -- exactly
    what a build backend drops silently.
    """
    import pathlib

    import spectroscopy
    notice = (pathlib.Path(spectroscopy.__file__).parent / 'data'
              / 'cd_reference' / 'LICENSE.DichroWebGit')
    assert notice.is_file(), "the DichroWebGit MIT notice is not shipped"
    text = notice.read_text()
    assert 'MIT License' in text
    assert 'Andy Miles' in text


def test_a_shipped_set_carries_its_licence_and_citation():
    """
    Citation is a *condition of use* of these data, so it has to travel with
    the set rather than live in a docstring. A result can then say what it
    owes without the caller having looked it up.
    """
    for name in ('sp175', 'smp180'):
        references = spc.datasets.reference_set(name)
        assert 'MIT' in references.info['licence']
        assert 'Andy Miles' in references.info['licence']
        assert references.info['citation']
        assert references.info['unit'] == 'delta epsilon'


def test_the_aqpz_example_is_cropped_to_where_it_was_measured():
    """
    The shipped scan stops at 197 nm because the J-815's photomultiplier
    voltage passes 600 V below that, and starved-detector output is not
    measurement. Shipping the full 180 nm scan would hand a reader 17 nm of
    detector noise that looks exactly like a band.
    """
    protein = spc.datasets.load('aqpz')
    assert protein.technique == 'CD'
    assert protein.y_unit == 'mdeg'
    assert protein.x.min() >= 196.0
    assert protein.metadata['n_residues'] == 254
    assert protein.metadata['sample'] == 'AqpZ-W14A'


def test_the_shipped_example_gives_the_documented_disagreement():
    """
    The tutorial's point, pinned: three methods, three answers, and the one
    with the best held-out score furthest from the crystal structure. If this
    ever stops being true the tutorial's argument has to change with it.
    """
    protein = spc.datasets.load('aqpz')
    references = spc.datasets.reference_set('smp180')
    crystal_helix = 178 / (231 + 23)

    answers = {method: cdm.estimate(protein, method, references)
               for method in ('nearest-shapes', 'ridge')}
    helix = {m: c.get('helix') for m, c in answers.items()}

    assert helix['nearest-shapes'] == pytest.approx(0.633, abs=0.01)
    assert helix['ridge'] == pytest.approx(0.449, abs=0.01)
    assert abs(helix['nearest-shapes'] - helix['ridge']) > 0.15
    assert (abs(helix['ridge'] - crystal_helix)
            > abs(helix['nearest-shapes'] - crystal_helix))


def test_the_answer_does_not_depend_on_where_the_scan_starts():
    """
    The wavelengths a fit uses come from the reference set, not from the
    measurement. A scan starting at 196.2 nm and the same scan resampled to
    whole nanometres describe the same protein, so they must give the same
    answer -- and once did not: 0.705 against 0.633 in helix, which is a large
    reply to a question about where a file begins.
    """
    references = spc.datasets.reference_set('smp180')
    shipped = spc.datasets.load('aqpz')

    # the same spectrum on a finer grid, offset from whole nanometres
    offset_x = np.arange(196.2, 280.0, 0.1)
    offset = spc.Spectrum(offset_x, np.interp(offset_x, shipped.x, shipped.y),
                          technique='CD', name='finer, offset')

    grids = [cdm.design_matrix(s, references, (190.0, 240.0))[0]
             for s in (shipped, offset)]
    assert np.allclose(grids[0], grids[1])
    assert np.allclose(grids[0], np.arange(197.0, 240.5, 1.0))

    for method in ('nearest-shapes', 'ridge'):
        answers = [cdm.estimate(s, method, references).get('helix')
                   for s in (shipped, offset)]
        assert answers[0] == pytest.approx(answers[1], abs=1e-9), method


def test_a_basis_that_disagrees_with_itself_falls_back_to_a_common_grid():
    """
    A basis assembled from several sources has no wavelengths they all hold,
    so there is nothing to anchor on and everybody gets interpolated. That is
    the compromise, and it should still produce a usable grid rather than
    refusing.
    """
    a = np.arange(190.0, 250.5, 1.0)
    b = np.arange(190.0, 250.25, 0.5)
    spectra = [spc.Spectrum(a, _band(a, 208, 7, -30), technique='CD', name='a'),
               spc.Spectrum(b, _band(b, 222, 9, -30), technique='CD', name='b')]
    references = lib.ReferenceSet(
        spectra, info={'categories': [HELIX, SHEET]})
    for spectrum, category in zip(references, ('helix', 'sheet')):
        spectrum.metadata['category'] = category

    sample = spc.Spectrum(a, _band(a, 210, 8, -30), technique='CD')
    grid, measured, design = cdm.design_matrix(sample, references,
                                               (190.0, 240.0))
    assert len(grid) > 5
    assert design.shape == (len(grid), 2)
