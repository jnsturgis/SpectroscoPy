# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Secondary structure from FTIR, and the vocabulary it reports in.

The estimate is checked against synthetic amide I envelopes whose composition
is known exactly. That tests the arithmetic, not the biology: whether the band
assignments are right for real proteins needs real proteins, and one with a
published structure is still wanted (roadmap section 15.2).

What is tested here just as hard is the comparison machinery, because that is
the contract the CD branch has to satisfy (roadmap section 17.2).
"""

import numpy as np
import pytest

from spectroscopy.lineshapes import spec_comp
from spectroscopy.processing.structure import AMIDE_I_BANDS, Category, Composition, from_ftir
from spectroscopy.spectra import Spectrum

#: (position, fwhm, amplitude) chosen so each falls in a different band.
TRUTH = [(1630.0, 16.0, 0.50), (1645.0, 14.0, 0.25),
         (1653.0, 13.0, 0.85), (1670.0, 15.0, 0.30)]
ETA = 0.6


def _expected_fractions():
    """Area fractions implied by TRUTH, from the closed forms."""
    gauss = np.sqrt(np.pi / (4 * np.log(2)))
    areas = np.array([(ETA * gauss + (1 - ETA) * np.pi / 2) * w * a
                      for _, w, a in TRUTH])
    return areas / areas.sum()


@pytest.fixture
def protein():
    x = np.linspace(1580.0, 1720.0, 561)
    y = sum(spec_comp(x, p, w, a, ETA) for p, w, a in TRUTH)
    return Spectrum(x, y, technique='ATR-FTIR', name='synthetic protein')


# ---------------------------------------------------------------------------
# the estimate
# ---------------------------------------------------------------------------

def test_it_recovers_a_known_composition(protein):
    composition = from_ftir(protein, method='amide-i-curve-fit',
                            positions=[p for p, _, _ in TRUTH])
    sheet, other, helix, turn = _expected_fractions()

    assert composition.get('sheet') == pytest.approx(sheet, abs=0.01)
    assert composition.get('other') == pytest.approx(other, abs=0.01)
    assert composition.get('helix') == pytest.approx(helix, abs=0.01)
    assert composition.get('turn') == pytest.approx(turn, abs=0.01)


def test_it_survives_a_residual_background(protein):
    """Water subtraction is never perfect; the derivative weighting is why."""
    x = protein.x
    sloped = Spectrum(x, protein.y + 0.03 * ((x - 1650) / 100) ** 2,
                      technique='ATR-FTIR')
    composition = from_ftir(sloped, method='amide-i-curve-fit',
                            positions=[p for p, _, _ in TRUTH])
    helix = _expected_fractions()[2]
    assert composition.get('helix') == pytest.approx(helix, abs=0.03)


def test_fractions_sum_to_one(protein):
    composition = from_ftir(protein, method='amide-i-curve-fit',
                            positions=[p for p, _, _ in TRUTH])
    assert sum(composition.estimated.values()) == pytest.approx(1.0)


def test_quality_is_populated(protein):
    composition = from_ftir(protein, method='amide-i-curve-fit',
                            positions=[p for p, _, _ in TRUTH])
    assert composition.quality['r_squared'] > 0.999
    assert composition.quality['components'] == 4
    assert composition.quality['model'] == 'voigt'
    assert len(composition.quality['position_stderr']) == 4


def test_it_finds_its_own_bands(protein):
    """The second derivative is the standard way in; check it works unaided."""
    composition = from_ftir(protein, method='amide-i-curve-fit')
    assert composition.get('helix') > composition.get('turn')


def test_the_method_must_be_named(protein):
    with pytest.raises(ValueError, match='needs an explicit method'):
        from_ftir(protein)


def test_an_unknown_method_is_refused(protein):
    with pytest.raises(ValueError, match='unknown method'):
        from_ftir(protein, method='eyeballing-it')


# ---------------------------------------------------------------------------
# aggregation: not a structure
# ---------------------------------------------------------------------------

def test_aggregation_is_reported_but_not_as_structure():
    """A failed preparation must not come back as beta-sheet."""
    x = np.linspace(1580.0, 1720.0, 561)
    y = (spec_comp(x, 1615.0, 14.0, 0.60, ETA)       # intermolecular
         + spec_comp(x, 1653.0, 13.0, 0.60, ETA))    # helix
    spectrum = Spectrum(x, y, technique='ATR-FTIR', name='aggregated')

    composition = from_ftir(spectrum, method='amide-i-curve-fit',
                            positions=[1615.0, 1653.0])
    aggregated = next(c for c in composition.fractions if c.name == 'aggregated')

    assert composition.get('aggregated') > 0.4
    assert composition.get('sheet') == 0.0
    assert aggregated.is_structural is False
    assert 'DSSP' in aggregated.note


def test_categories_know_whether_they_are_structural():
    assert Category('helix', frozenset({'H'})).is_structural
    assert not Category('aggregated', frozenset()).is_structural


def test_the_band_table_is_contiguous():
    """No gaps: a fitted component must not fall between two bands."""
    edges = [span for _, span in AMIDE_I_BANDS]
    for (_, high), (low, _) in zip(edges, edges[1:]):
        assert high == low, f"gap between {high} and {low}"


# ---------------------------------------------------------------------------
# comparison -- the contract the CD branch must satisfy
# ---------------------------------------------------------------------------

def _cd_style():
    """A CDSSTR-shaped estimate: positional splits, and unordered
    swallowing G, I and B."""
    return Composition(
        fractions={
            Category('regular-helix', frozenset({'H'}), note='positional'): 0.30,
            Category('distorted-helix', frozenset({'H'}), note='positional'): 0.12,
            Category('regular-sheet', frozenset({'E'})): 0.16,
            Category('distorted-sheet', frozenset({'E'})): 0.09,
            Category('turns', frozenset({'T'})): 0.14,
            Category('unordered', frozenset({'S', '-', 'B', 'I', 'G'})): 0.19,
        },
        method='cdsstr', technique='CD', quality={'nrmsd': 0.03})


def test_comparison_offers_the_by_name_view(protein):
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    result = ftir.compare(_cd_style())

    assert set(result.nominal) == {'helix', 'sheet', 'turn'}
    # regular + distorted must be summed into the parent
    assert result.nominal['helix'][1] == pytest.approx(0.42)
    assert result.nominal['sheet'][1] == pytest.approx(0.25)


def test_comparison_names_the_states_that_make_it_approximate(protein):
    """The actionable half: not 'incomparable' but 'and here is why'."""
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    caveats = " ".join(ftir.compare(_cd_style()).caveats)

    assert 'G, I' in caveats, "the 3-10/pi mismatch should be named"
    assert 'B' in caveats
    assert 'helix' in caveats


def test_the_strict_grouping_refuses_to_separate_what_overlaps(protein):
    """CD's unordered claims G, I and B, so it touches helix and sheet.

    The rigorous answer is that only turns are separately comparable. It is a
    disappointing answer and it is the true one, which is why the by-name view
    exists beside it rather than instead of it.
    """
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    result = ftir.compare(_cd_style())

    assert any('turn' in name for name in result.groups)
    assert any(len(name.split('+')) > 3 for name in result.groups)


def test_matching_vocabularies_compare_cleanly(protein):
    """Two amide I estimates share a vocabulary, so nothing merges."""
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    result = ftir.compare(ftir)

    assert result.rmsd == pytest.approx(0.0)
    assert all(len(name.split('+')) == 1 for name in result.groups)
    assert not result.caveats


def test_non_structural_categories_are_excluded_from_comparison():
    x = np.linspace(1580.0, 1720.0, 561)
    y = (spec_comp(x, 1615.0, 14.0, 0.5, ETA)
         + spec_comp(x, 1653.0, 13.0, 0.8, ETA))
    ftir = from_ftir(Spectrum(x, y, technique='ATR-FTIR'),
                     method='amide-i-curve-fit', positions=[1615.0, 1653.0])
    result = ftir.compare(_cd_style())

    assert any('aggregated' in key for key in result.excluded)


def test_unestimated_categories_are_excluded_not_zeroed(protein):
    """theta-222 estimates helix and nothing else; that is not 'no sheet'."""
    partial = Composition(
        fractions={Category('helix', frozenset({'H', 'G', 'I'})): 0.40,
                   Category('sheet', frozenset({'E', 'B'})): None},
        method='theta-222', technique='CD')
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    result = ftir.compare(partial)

    assert any('not estimated' in reason for reason in result.excluded.values())
    assert 'sheet' not in result.nominal


def test_largest_disagreement_points_at_something(protein):
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    name, difference = ftir.compare(_cd_style()).largest_disagreement
    assert name
    assert isinstance(difference, float)


def test_str_shows_both_views(protein):
    ftir = from_ftir(protein, method='amide-i-curve-fit',
                     positions=[p for p, _, _ in TRUTH])
    text = str(ftir.compare(_cd_style()))
    assert 'by name' in text and 'strictly comparable' in text and 'RMSD' in text


def test_composition_str_marks_the_non_dssp_category(protein):
    composition = from_ftir(protein, method='amide-i-curve-fit',
                            positions=[p for p, _, _ in TRUTH])
    assert 'not a DSSP state' in str(composition)
