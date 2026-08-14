# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Scattering correction, against synthetic ground truth.

The test that matters is not "does the background look plausible" but "does
removing it recover the concentrations", so most of these end in an unmix.
"""

import warnings

import numpy as np
import pytest

from spectroscopy import Spectrum
from spectroscopy.library import Library, Reference
from spectroscopy.lineshapes import gauss
from spectroscopy.processing.scattering import (
    correct_scattering,
    from_references,
    scatter_baseline,
)
from spectroscopy.processing.unmix import unmix

X = np.linspace(220.0, 400.0, 361)
DNA_EPS = 20.0 * gauss(X, 260, 22, 1.0) + 8.0 * gauss(X, 230, 18, 1.0)
PROTEIN_EPS = 1.0 * gauss(X, 280, 14, 1.0) + 12.0 * gauss(X, 225, 12, 1.0)

#: Deliberately not a single power law: two exponents, neither of them 4.
SCATTER = 0.6 * (X / 400.0) ** -3.2 + 0.25 * (X / 400.0) ** -1.5


def _uv(y, name='s'):
    return Spectrum(X, y, x_quantity='Wavelength', x_unit='nm',
                    y_quantity='Absorbance', y_unit='absorbance', name=name)


@pytest.fixture
def references():
    return Library([
        Reference('dsDNA', _uv(DNA_EPS), unit='(ug/mL)^-1 cm^-1'),
        Reference('protein', _uv(PROTEIN_EPS), unit='(ug/mL)^-1 cm^-1'),
    ])


@pytest.fixture
def turbid():
    return _uv(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS + SCATTER, 'turbid')


def test_scattering_wrecks_an_uncorrected_unmix(references, turbid):
    """
    Why this belongs before unmixing rather than after.

    The background is not one of the components, so the fit distributes it
    across the ones it has -- here inflating the nucleic acid eightfold.
    """
    wrong = unmix(turbid, references)
    assert wrong['dsDNA'] > 5 * 0.03
    assert wrong['protein'] > 1.5 * 0.40


def test_correction_recovers_the_concentrations(references, turbid):
    corrected = correct_scattering(turbid)
    result = unmix(corrected, references)
    assert result['dsDNA'] == pytest.approx(0.03, abs=0.002)
    assert result['protein'] == pytest.approx(0.40, abs=0.01)


def test_the_basis_beats_a_fixed_rayleigh_exponent(turbid):
    """
    A single exponent is a guess about particle size. The sample here is a
    mixture of two, which is what a polydisperse suspension looks like.
    """
    basis = np.max(np.abs(scatter_baseline(X, turbid.y) - SCATTER))
    rayleigh = np.max(np.abs(scatter_baseline(X, turbid.y, exponents=[4.0])
                             - SCATTER))
    assert basis < rayleigh / 5


def test_the_background_is_returned_for_inspection(turbid):
    corrected, background = correct_scattering(turbid, return_background=True)
    assert np.allclose(background.y, SCATTER, atol=0.02)
    assert np.allclose(corrected.y + background.y, turbid.y)
    assert 'scattering' in background.name


def test_the_correction_is_recorded(turbid):
    corrected = correct_scattering(turbid)
    step = corrected.history[-1]
    assert 'correct_scattering' in str(step)
    assert '320' in str(step)


def test_the_background_is_never_negative(turbid):
    """Scattering adds signal; a negative contribution is unphysical."""
    background = scatter_baseline(X, turbid.y)
    assert np.all(background >= -1e-12)


def test_fitting_through_a_band_is_warned_about():
    """
    The failure this cannot silently survive: a fit region that is not free of
    absorption lets the background eat the band.
    """
    absorbing_everywhere = _uv(2.0 * gauss(X, 360, 40, 1.0) + SCATTER)
    with pytest.warns(UserWarning, match='absorption-free'):
        correct_scattering(absorbing_everywhere, region=(340.0, 380.0))


def test_too_narrow_a_region_is_refused(turbid):
    with pytest.raises(ValueError, match='not enough'):
        correct_scattering(turbid, region=(399.0, 400.0))


def test_measured_backgrounds_can_be_used_instead(references, turbid):
    """The honest reference where the blanks exist."""
    blanks = [_uv(0.5 * (X / 400.0) ** -3.2), _uv(0.5 * (X / 400.0) ** -1.5)]
    corrected = from_references(turbid, blanks)

    result = unmix(corrected, references)
    assert result['dsDNA'] == pytest.approx(0.03, abs=1e-3)
    assert result['protein'] == pytest.approx(0.40, abs=1e-3)
    assert 'measured backgrounds' in str(corrected.history[-1])


def test_measured_backgrounds_are_resampled(references, turbid):
    coarse = np.linspace(220.0, 400.0, 90)
    blanks = [Spectrum(coarse, np.interp(coarse, X, 0.5 * (X / 400.0) ** -3.2),
                       x_unit='nm', y_unit='absorbance'),
              Spectrum(coarse, np.interp(coarse, X, 0.5 * (X / 400.0) ** -1.5),
                       x_unit='nm', y_unit='absorbance')]
    corrected = from_references(turbid, blanks)
    assert len(corrected.x) == len(X)
    assert unmix(corrected, references)['dsDNA'] == pytest.approx(0.03, abs=5e-3)


def test_no_backgrounds_is_refused(turbid):
    with pytest.raises(ValueError, match='no scattering backgrounds'):
        from_references(turbid, [])


def test_a_clean_spectrum_is_left_almost_alone(references):
    """Correcting something that does not scatter should do nearly nothing."""
    clean = _uv(0.03 * DNA_EPS + 0.40 * PROTEIN_EPS)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        corrected = correct_scattering(clean)
    assert np.max(np.abs(corrected.y - clean.y)) < 0.02


# ---------------------------------------------------------------------------
# a blank says nothing outside where it was measured
# ---------------------------------------------------------------------------

def test_a_blank_is_not_followed_past_its_own_range():
    """
    Scattering only ever adds signal, so a negative background is not a small
    error but a contradiction. Extrapolating a blank measured 250-350 nm out
    to 200 gave -4.04 on a sample spanning 0.5 to 0.501 -- four absorbance
    units of invented band, waiting to be subtracted.
    """
    sample_x = np.linspace(200.0, 400.0, 201)
    sample = Spectrum(sample_x, 0.5 + 1e6 * sample_x ** -4.0,
                      technique='UV-Vis', name='turbid')
    blank_x = np.linspace(250.0, 350.0, 101)
    blank = Spectrum(blank_x, 1e6 * blank_x ** -4.0, technique='UV-Vis',
                     name='blank')

    with pytest.warns(UserWarning, match='returned uncorrected'):
        corrected, background = from_references(
            sample, [blank], region=(300.0, 350.0), return_background=True)

    assert background.y.min() >= 0.0


def test_the_spectrum_keeps_its_full_length_and_is_left_alone_outside():
    """
    No cropping: what cannot be corrected is returned as it arrived, rather
    than removed. Declining to correct is honest; inventing a background to
    subtract is not.
    """
    sample_x = np.linspace(200.0, 400.0, 201)
    sample = Spectrum(sample_x, 0.5 + 1e6 * sample_x ** -4.0,
                      technique='UV-Vis', name='turbid')
    blank_x = np.linspace(250.0, 350.0, 101)
    blank = Spectrum(blank_x, 1e6 * blank_x ** -4.0, technique='UV-Vis')

    with pytest.warns(UserWarning):
        corrected = from_references(sample, [blank],
                                               region=(300.0, 350.0))

    assert len(corrected) == len(sample)
    outside = (corrected.x < 250.0) | (corrected.x > 350.0)
    assert np.allclose(corrected.y[outside], sample.y[outside])
    assert corrected.history[-1].params['corrected_over'] == (250.0, 350.0)


def test_backgrounds_that_share_no_wavelengths_are_refused():
    x = np.linspace(200.0, 400.0, 201)
    sample = Spectrum(x, np.ones_like(x), technique='UV-Vis')
    low_x, high_x = np.linspace(200.0, 240.0, 41), np.linspace(300.0, 400.0, 101)
    with pytest.raises(ValueError, match='no wavelengths in common'):
        from_references(
            sample, [Spectrum(low_x, np.ones_like(low_x), technique='UV-Vis'),
                     Spectrum(high_x, np.ones_like(high_x), technique='UV-Vis')])
