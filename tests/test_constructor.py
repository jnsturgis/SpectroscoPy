# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Building a Spectrum: from data, from a file, from another Spectrum.

Until 0.1.1 there was no way to construct a spectrum from arrays. Anything
computed rather than read -- a simulation, a model, a difference worked out
elsewhere -- had to build an empty object and assign ``.x`` and ``.y``, which
the test fixtures in this suite still show being done. ADR-0001 section 7.1
called that a blocker on the 1.0 freeze, because freezing it would have made
empty-then-assign the API.
"""

import pathlib

import numpy as np
import pytest

from spectroscopy import datasets
from spectroscopy.spectra import Spectrum


@pytest.fixture
def band():
    x = np.linspace(900.0, 1800.0, 901)
    return x, np.exp(-((x - 1650.0) / 40.0) ** 2)


# ---------------------------------------------------------------------------
# from data
# ---------------------------------------------------------------------------

def test_arrays_make_a_spectrum(band):
    x, y = band
    spectrum = Spectrum(x, y)
    assert len(spectrum) == 901
    np.testing.assert_array_equal(spectrum.x, x)
    np.testing.assert_array_equal(spectrum.y, y)


def test_lists_are_accepted_and_become_floats():
    spectrum = Spectrum([1, 2, 3], [4, 5, 6])
    assert spectrum.x.dtype == np.float64
    assert spectrum.y.dtype == np.float64


def test_technique_sets_the_conventional_axes(band):
    spectrum = Spectrum(*band, technique='ATR-FTIR')
    assert spectrum.technique == 'ATR-FTIR'
    assert spectrum.x_unit == 'cm^-1'
    assert spectrum.metadata['spec_type'] == 'ATR-FTIR'
    assert spectrum.reversed_x is True


def test_explicit_units_win_over_the_technique(band):
    """A fluorescence excitation scan in nm is not the technique's default."""
    spectrum = Spectrum(*band, technique='Raman', x_unit='nm',
                        x_quantity='Wavelength')
    assert spectrum.x_unit == 'nm'
    assert spectrum.technique == 'Raman'


def test_explicit_units_survive_a_later_set_type(band):
    """Passing units marks them authoritative, like a file declaring them.

    Otherwise set_type() would treat them as defaults and relabel the axis,
    which is the same failure that mislabels a JCAMP transmittance file.
    """
    spectrum = Spectrum(*band, y_unit='%T', y_quantity='Transmittance')
    spectrum.set_type('FTIR')
    assert spectrum.y_unit == '%T'


def test_metadata_is_copied_not_referenced(band):
    original = {'sample': 'PG_coli', 'nested': {'a': 1}}
    spectrum = Spectrum(*band, metadata=original)
    spectrum.metadata['nested']['a'] = 2
    assert original['nested']['a'] == 1


def test_building_from_data_records_no_history(band):
    """Creating data is not processing it."""
    assert Spectrum(*band).history == []


def test_history_can_be_supplied(band):
    from spectroscopy.history import ProcessingStep
    step = ProcessingStep(name='smooth', params={'method': 'savgol'})
    assert Spectrum(*band, history=[step]).history == [step]


def test_a_constructed_spectrum_chains_like_any_other(band):
    result = Spectrum(*band, technique='ATR-FTIR').crop(1500, 1700).normalize()
    assert result.y.max() == pytest.approx(1.0)
    assert [s.name for s in result.history] == ['crop', 'normalize']


def test_it_has_no_path(band):
    assert Spectrum(*band).path is None


# ---------------------------------------------------------------------------
# rejections
# ---------------------------------------------------------------------------

def test_mismatched_lengths_are_refused(band):
    x, y = band
    with pytest.raises(ValueError, match='same length'):
        Spectrum(x, y[:5])


def test_two_dimensional_input_is_refused_by_name():
    """The error should say where a 2-D array belongs."""
    grid = np.zeros((2, 3))
    with pytest.raises(ValueError, match='SpectrumCollection'):
        Spectrum(grid, grid)


def test_keywords_are_refused_on_the_file_form():
    with pytest.raises(TypeError, match='only accepted by the data form'):
        Spectrum('sample.dpt', technique='Raman')


def test_unknown_initialiser_still_raises():
    with pytest.raises(TypeError):
        Spectrum(42)


# ---------------------------------------------------------------------------
# read()
# ---------------------------------------------------------------------------

def test_read_infers_the_format():
    assert len(Spectrum.read(datasets.path('ethanol'))) > 0


def test_read_takes_the_format_as_a_named_argument():
    """The call the old __init__ docstring wrongly promised.

    ``Spectrum("ethanol.jdx", "jcamp")`` raised TypeError, because two
    positional strings mean (directory, filename). ADR-0001 section 7.2.
    """
    assert len(Spectrum.read(datasets.path('ethanol'), 'jcamp')) > 0


def test_read_and_the_positional_form_agree():
    path = datasets.path('ethanol')
    np.testing.assert_array_equal(Spectrum.read(path).y, Spectrum(str(path)).y)


def test_read_accepts_a_path_object():
    assert len(Spectrum.read(pathlib.Path(datasets.path('ethanol')))) > 0


def test_the_constructor_accepts_a_path_object():
    assert len(Spectrum(pathlib.Path(datasets.path('ethanol')))) > 0


# ---------------------------------------------------------------------------
# the older forms still work
# ---------------------------------------------------------------------------

def test_empty_is_still_empty():
    spectrum = Spectrum()
    assert spectrum.name == 'unnamed'
    assert spectrum.history == []


def test_copy_is_still_a_deep_copy(band):
    original = Spectrum(*band, metadata={'sample': 'a'})
    duplicate = Spectrum(original)
    duplicate.metadata['sample'] = 'b'
    assert original.metadata['sample'] == 'a'


def test_directory_and_filename_form_still_works():
    path = pathlib.Path(datasets.path('ethanol'))
    assert len(Spectrum(str(path.parent) + '/', path.name)) > 0
