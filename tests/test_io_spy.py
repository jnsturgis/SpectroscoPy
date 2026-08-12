# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The native .spy format.

.spy exists to be the lossless one (TODO.txt item 1). Version 0.0 was not: it
wrote name and axis labels that its own reader never read back. Version 1.0
carries identity, units, metadata and -- the point of it -- the processing
history, which is what makes provenance survive leaving memory.
"""

import json

import numpy as np
import pytest

import spectroscopy as spc
from spectroscopy.io import spy
from spectroscopy.spectra import Spectrum


@pytest.fixture
def processed():
    """A spectrum that has actually been through a pipeline."""
    spec = Spectrum()
    spec.x = np.linspace(900.0, 1800.0, 91)
    spec.y = np.exp(-((spec.x - 1650.0) ** 2) / (2 * 40.0 ** 2)) + 0.1
    spec.set_type("ATR-FTIR")
    spec.set_sample("Glucose")
    spec.name = "Glucose average"

    reference = Spectrum(spec)
    reference.set_sample("H2O")

    return (spec
            .subtract_reference(reference, factor=0.7)
            .crop(1000, 1750)
            .baseline_correct('rubberband')
            .normalize('max', window=(1050, 1080)))


def _round_trip(spectrum, tmp_path, name="rt.spy"):
    target = tmp_path / name
    spectrum.save_as(str(target), "spy")
    return Spectrum("", str(target), "spy"), target


# ---------------------------------------------------------------------------
# 1.0 round trip
# ---------------------------------------------------------------------------

def test_data_survives(processed, tmp_path):
    back, _ = _round_trip(processed, tmp_path)
    assert np.allclose(back.x, processed.x)
    assert np.allclose(back.y, processed.y, atol=1e-9)


def test_identity_and_units_survive(processed, tmp_path):
    back, _ = _round_trip(processed, tmp_path)
    assert back.name == processed.name
    assert back.technique == processed.technique
    assert back.x_unit == processed.x_unit
    assert back.y_unit == processed.y_unit
    assert back.x_quantity == processed.x_quantity
    assert back.x_label == processed.x_label


def test_metadata_survives(processed, tmp_path):
    back, _ = _round_trip(processed, tmp_path)
    assert back.metadata == processed.metadata


def test_history_survives(processed, tmp_path):
    """The whole reason for a native format."""
    back, _ = _round_trip(processed, tmp_path)

    assert [s.name for s in back.history] == [s.name for s in processed.history]
    assert [s.params for s in back.history] == [s.params for s in processed.history]
    assert [s.timestamp for s in back.history] == \
        [s.timestamp for s in processed.history]


def test_the_hand_tuned_factor_survives(processed, tmp_path):
    """
    The scientific crux: the 0.7 water factor chosen by eye must come back out
    of the file, not just the numbers it produced.
    """
    back, _ = _round_trip(processed, tmp_path)
    step = next(s for s in back.history if s.name == 'subtract_reference')
    assert step.params['factor'] == 0.7
    assert step.params['reference']['sample'] == 'H2O'


def test_a_second_round_trip_changes_nothing(processed, tmp_path):
    once, _ = _round_trip(processed, tmp_path, "one.spy")
    twice, _ = _round_trip(once, tmp_path, "two.spy")
    assert np.allclose(twice.y, once.y)
    assert [s.to_dict() for s in twice.history] == [s.to_dict() for s in once.history]


def test_the_file_is_readable_by_eye(processed, tmp_path):
    """Data stays tab separated so ordinary tools still work on it."""
    _, target = _round_trip(processed, tmp_path)
    lines = target.read_text().splitlines()

    assert lines[0].startswith("# spy format 1.0")
    data_start = lines.index("# data")
    numbers = lines[data_start + 2].split("\t")
    assert len(numbers) == 2
    float(numbers[0]), float(numbers[1])

    header = json.loads("\n".join(lines[2:data_start]))
    assert header['name'] == processed.name
    assert len(header['history']) == len(processed.history)


def test_a_label_override_survives(tmp_path):
    """Labels read from a CSV header are an override, not a technique default."""
    spec = Spectrum()
    spec.x, spec.y = np.array([1.0, 2.0]), np.array([3.0, 4.0])
    spec.x_label = "Whatever the file said"
    back, _ = _round_trip(spec, tmp_path)
    assert back.x_label == "Whatever the file said"


# ---------------------------------------------------------------------------
# legacy 0.0
# ---------------------------------------------------------------------------

def test_legacy_files_still_load(tmp_path):
    legacy = tmp_path / "old.spy"
    legacy.write_text(
        "# spy format 0.0 file\n"
        "OldSpectrum\n"
        "Wavenumber\tAbsorbance\n"
        "# x,y data\n"
        "1000.000\t0.10000\n"
        "1001.000\t0.20000\n"
        "# metadata\n"
        "{'sample': 'legacy', 'spec_type': 'ATR-FTIR'}\n"
    )
    spec = Spectrum("", str(legacy), "spy")

    assert np.allclose(spec.x, [1000.0, 1001.0])
    assert np.allclose(spec.y, [0.1, 0.2])
    assert spec.metadata['sample'] == 'legacy'
    # ... and the name and labels the 0.0 reader used to drop:
    assert spec.name == "OldSpectrum"
    assert spec.x_label == "Wavenumber"


def test_version_is_taken_from_the_file_not_the_caller(tmp_path):
    """
    reload() used to pass format='0.0' unconditionally. Trusting the caller
    over the file is how a 1.0 file gets parsed as 0.0.
    """
    spec = Spectrum()
    spec.x, spec.y = np.array([1.0, 2.0]), np.array([3.0, 4.0])
    spec.name = "Modern"
    target = tmp_path / "modern.spy"
    spec.save_as(str(target), "spy")

    with open(target, encoding="utf-8") as handle:
        restored = spy.read(handle, format='0.0')    # deliberately wrong
    assert restored[0].name == "Modern"


def test_an_unknown_version_is_refused(tmp_path):
    path = tmp_path / "future.spy"
    path.write_text("# spy format 9.9\n# header\n{}\n# data\n")
    with pytest.raises(ValueError, match="unsupported .spy format"):
        Spectrum("", str(path), "spy")


def test_an_empty_file_is_refused(tmp_path):
    path = tmp_path / "empty.spy"
    path.write_text("")
    with pytest.raises(ValueError, match="empty"):
        Spectrum("", str(path), "spy")


# ---------------------------------------------------------------------------
# a whole collection in one .spy file
#
# ADR-0004 section 5 left this open and it is on the 1.0 critical path, because
# the native format freezes in November. The decision (James, 2026-08-13) is
# one format that holds either a spectrum or a set, not a second format beside
# it. What these tests are really about is that a set has facts of its own --
# where it came from, what may be done with it, what its numbers are in -- and
# that saving the spectra one at a time loses every one of them.
# ---------------------------------------------------------------------------

def _melt():
    x = np.linspace(200.0, 250.0, 11)
    spectra = []
    for index in range(3):
        spectrum = Spectrum(x, np.sin(x / 10.0) + index, technique='CD',
                            name=f'scan {index}')
        spectrum.set_parameter(30.0 + 10.0 * index)
        spectra.append(spectrum)
    return spc.SpectrumCollection(spectra, name='AqpZ melt',
                                  info={'parameter_name': 'temperature',
                                        'parameter_unit': 'C'})


def _reference_set():
    from spectroscopy import library as lib
    from spectroscopy.processing.structure import Category, Composition

    x = np.linspace(190.0, 240.0, 11)
    helix = Category('helix', frozenset({'G', 'H', 'I'}))
    sheet = Category('sheet', frozenset({'E', 'B'}))
    spectra = [Spectrum(x, np.cos(x / 8.0) + index, technique='CD',
                        name=f'ref{index}') for index in range(3)]
    compositions = [
        Composition(fractions={helix: 0.8 - 0.3 * index,
                               sheet: 0.2 + 0.3 * index},
                    method='DSSP on a crystal structure', technique='X-ray')
        for index in range(3)]
    return lib.ReferenceSet.from_compositions(
        spectra, compositions, name='SMP180',
        info={'source': 'DichroWebGit SMP180',
              'licence': 'MIT (c) 2023 Andy Miles',
              'citation': 'Abdul-Gader et al. 2011',
              'unit': 'delta epsilon'})


def test_a_collection_round_trips(tmp_path):
    original = _melt()
    original.save_as(tmp_path / 'melt.spy')
    back = spc.io.read_spectra(tmp_path / 'melt.spy')

    assert len(back) == 3
    assert back.name == 'AqpZ melt'
    assert [s.name for s in back] == [s.name for s in original]
    for before, after in zip(original, back):
        assert np.allclose(after.x, before.x)
        assert np.allclose(after.y, before.y)
        assert after.technique == before.technique
    assert list(back.parameters) == [30.0, 40.0, 50.0]


def test_the_set_level_facts_survive(tmp_path):
    """
    The whole reason a set needed a file of its own. Saved as three separate
    spectra, 'temperature in C' has nowhere to live and is simply gone.
    """
    _melt().save_as(tmp_path / 'melt.spy')
    back = spc.io.read_spectra(tmp_path / 'melt.spy')

    assert back.info['parameter_name'] == 'temperature'
    assert back.parameter_unit == 'C'


def test_a_single_spectrum_file_is_written_exactly_as_before(tmp_path):
    """
    Collections cost single-spectrum files nothing: same '# header' marker,
    same layout, no collection block. Every file ever written still reads, and
    every file written from now on is still readable by what came before.
    """
    spectrum = Spectrum(np.array([1.0, 2.0]), np.array([3.0, 4.0]),
                        name='Glucose')
    spectrum.save_as(str(tmp_path / 'one.spy'))
    text = (tmp_path / 'one.spy').read_text()

    assert text.startswith('# spy format 1.0\n# header\n')
    assert '# collection' not in text
    assert '# spectrum' not in text


def test_a_single_spectrum_file_reads_as_a_set_of_one(tmp_path):
    """
    Not an error: a set of one is a set. It is what lets a caller treat both
    kinds of file the same way.
    """
    spectrum = Spectrum(np.array([1.0, 2.0]), np.array([3.0, 4.0]),
                        name='Glucose')
    spectrum.save_as(str(tmp_path / 'one.spy'))
    loaded = spc.io.read_spectra(tmp_path / 'one.spy')

    assert len(loaded) == 1
    assert loaded[0].name == 'Glucose'


def test_reading_one_spectrum_from_a_set_of_several_refuses(tmp_path):
    """The existing contract of read_spectrum, unchanged by collections."""
    _melt().save_as(tmp_path / 'melt.spy')
    with pytest.raises(ValueError, match='holds 3 spectra'):
        spc.read(tmp_path / 'melt.spy')


def test_reading_one_spectrum_from_a_set_of_one_just_works(tmp_path):
    """
    No special rule for a collection file that holds one spectrum: you asked
    for a spectrum, there is exactly one, you get it.
    """
    single = spc.SpectrumCollection([_melt()[0]], name='just the one',
                                    info={'source': 'somewhere'})
    single.save_as(tmp_path / 'one_set.spy')

    assert spc.read(tmp_path / 'one_set.spy').name == 'scan 0'


def test_writing_a_set_to_a_one_spectrum_format_refuses(tmp_path):
    """
    Better than writing the first of three, and better than writing three files
    that no longer know they belong together.
    """
    with pytest.raises(ValueError, match='stores one spectrum'):
        _melt().save_as(tmp_path / 'melt.csv', 'csv')


# -- what makes a reference set a reference set ------------------------------

def test_a_reference_set_comes_back_a_reference_set(tmp_path):
    """
    ``kind`` is recorded, so the set does not degrade to a plain collection
    that has quietly lost what it was -- and lost ``.compositions`` with it,
    while the truth sits unread in each spectrum's metadata.
    """
    from spectroscopy import library as lib

    _reference_set().save_as(tmp_path / 'SMP180.spy')
    back = spc.io.read_spectra(tmp_path / 'SMP180.spy')

    assert isinstance(back, lib.ReferenceSet)
    assert back.has_truth and not back.is_structural
    assert [c.get('helix') for c in back.compositions] == pytest.approx(
        [0.8, 0.5, 0.2])


def test_the_licence_travels_with_the_data(tmp_path):
    """
    A citation condition is an obligation of the whole set. Smeared across 128
    spectra it is no record of one; dropped on save it is no record at all.
    """
    _reference_set().save_as(tmp_path / 'SMP180.spy')
    back = spc.io.read_spectra(tmp_path / 'SMP180.spy')

    assert back.info['licence'] == 'MIT (c) 2023 Andy Miles'
    assert back.info['citation'] == 'Abdul-Gader et al. 2011'
    assert back.info['source'] == 'DichroWebGit SMP180'


def test_categories_come_back_with_their_dssp_states(tmp_path):
    """
    JSON has no form for a Category, and a bare name is not a vocabulary --
    'helix' does not say whether it covers 3-10 and pi. The states are written
    out and read back, which is what makes this lossless rather than nearly so.
    """
    _reference_set().save_as(tmp_path / 'SMP180.spy')
    back = spc.io.read_spectra(tmp_path / 'SMP180.spy')

    helix = next(c for c in back.categories if c.name == 'helix')
    assert helix.states == frozenset({'G', 'H', 'I'})
    assert all(category.states for category in back.compositions[0].fractions)


def test_known_truth_survives_and_still_names_its_source(tmp_path):
    _reference_set().save_as(tmp_path / 'SMP180.spy')
    back = spc.io.read_spectra(tmp_path / 'SMP180.spy')

    assert back[0].metadata['known_from'] == 'DSSP on a crystal structure'
    assert back.compositions[0].method == 'DSSP on a crystal structure'


def test_an_unknown_kind_loads_with_a_warning(tmp_path):
    """
    A file from a later version is still worth reading for its spectra and its
    provenance; what it must not do is pretend the class it named was restored.
    """
    _reference_set().save_as(tmp_path / 'SMP180.spy')
    path = tmp_path / 'SMP180.spy'
    path.write_text(path.read_text().replace('"ReferenceSet"', '"FutureSet"', 1))

    with pytest.warns(UserWarning, match='FutureSet'):
        back = spc.io.read_spectra(path)
    assert type(back) is spc.SpectrumCollection
    assert len(back) == 3
    assert back.info['licence'] == 'MIT (c) 2023 Andy Miles'


# -- the file itself ---------------------------------------------------------

def test_nothing_is_written_at_both_levels(tmp_path):
    """
    A format that copied the licence onto every spectrum would put back exactly
    the disagreement ADR-0004 removed: two copies that can differ, with no way
    to tell which was meant.
    """
    _reference_set().save_as(tmp_path / 'SMP180.spy')
    text = (tmp_path / 'SMP180.spy').read_text()

    assert text.count('MIT (c) 2023 Andy Miles') == 1
    assert text.count('# spectrum') == 3


def test_the_collection_file_stays_greppable(tmp_path):
    """Same bargain as always: structured to read back, plain to look at."""
    _reference_set().save_as(tmp_path / 'SMP180.spy')
    text = (tmp_path / 'SMP180.spy').read_text()

    assert text.startswith('# spy format 1.0')
    assert 'DichroWebGit SMP180' in text
    assert '\t' in text                     # the numbers are still columns


def test_two_collection_blocks_are_refused(tmp_path):
    """Two sets concatenated into one file: there is no answer, so it says so."""
    _melt().save_as(tmp_path / 'melt.spy')
    doubled = tmp_path / 'doubled.spy'
    text = (tmp_path / 'melt.spy').read_text()
    doubled.write_text(text + text.split('\n', 1)[1])

    with pytest.raises(ValueError, match='two .# collection. blocks'):
        spc.io.read_spectra(doubled)


def test_an_empty_collection_round_trips(tmp_path):
    """A set of nothing is still a set, and still has a provenance."""
    empty = spc.SpectrumCollection([], name='nothing yet',
                                   info={'source': 'a plan'})
    empty.save_as(tmp_path / 'empty.spy')
    back = spc.io.read_spectra(tmp_path / 'empty.spy')

    assert len(back) == 0
    assert back.info['source'] == 'a plan'
