# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
SpectrumCollection -- the answer to friction item #1, replicate averaging by
hard-coded index.
"""

import numpy as np
import pytest

from spectroscopy.collection import SpectrumCollection
from spectroscopy.spectra import Spectrum

X = np.linspace(900.0, 1800.0, 91)


def _spectrum(sample, level, flat=True):
    """
    A replicate at intensity ``level``.

    Flat by default because it makes the averaging assertions exact and easy to
    read; pass ``flat=False`` for a spectrum with an actual band, needed by
    anything that baselines or normalises (a flat spectrum baseline-corrects to
    all zeros, which cannot then be scaled -- correctly).
    """
    spec = Spectrum()
    if flat:
        spec.y = np.full_like(X, float(level))
    else:
        spec.y = float(level) * np.exp(-((X - 1650.0) ** 2) / (2 * 40.0 ** 2)) \
            + 0.05 * float(level)
    spec.x = X.copy()
    spec.set_type("ATR-FTIR")
    spec.set_sample(sample)
    return spec


@pytest.fixture
def collection():
    """Three replicates each of two samples, as the notebooks always have."""
    return SpectrumCollection([
        _spectrum("PG_coli", 1), _spectrum("PG_coli", 2), _spectrum("PG_coli", 3),
        _spectrum("PG_myxo", 10), _spectrum("PG_myxo", 20), _spectrum("PG_myxo", 30),
    ])


@pytest.fixture
def dpt_tree(tmp_path):
    """A directory of .dpt replicates named the way OPUS writes them."""
    for sample, count, level in [("PG_coli", 3, 1.0), ("PG_myxo", 2, 5.0)]:
        for index in range(count):
            path = tmp_path / f"{sample}.{index}.dpt"
            path.write_text("".join(f"{x}\t{level + index}\n" for x in X))
    return tmp_path


# ---------------------------------------------------------------------------
# the point of the class
# ---------------------------------------------------------------------------

def test_group_by_sample_then_mean(collection):
    """
    Replaces `(spectra[0] + spectra[1] + spectra[2]) / 3.0`, which breaks
    silently the moment a file is added.
    """
    averages = {name: group.mean()
                for name, group in collection.group_by('sample').items()}

    assert set(averages) == {"PG_coli", "PG_myxo"}
    assert np.allclose(averages["PG_coli"].y, 2.0)
    assert np.allclose(averages["PG_myxo"].y, 20.0)


def test_the_average_matches_the_hand_written_index_arithmetic(collection):
    by_hand = (collection[0] + collection[1] + collection[2]) / 3.0
    by_group = collection.group_by('sample')["PG_coli"].mean()
    assert np.allclose(by_hand.y, by_group.y)


def test_adding_a_replicate_changes_the_average(collection):
    """The failure mode of index arithmetic: it would not have noticed."""
    before = collection.group_by('sample')["PG_coli"].mean()
    extended = collection + _spectrum("PG_coli", 6)
    after = extended.group_by('sample')["PG_coli"].mean()
    assert np.allclose(before.y, 2.0)
    assert np.allclose(after.y, 3.0)


def test_mean_records_how_many_spectra_went_in(collection):
    mean = collection.group_by('sample')["PG_coli"].mean()
    step = mean.history[-1]
    assert step.name == 'mean'
    assert step.params['n_spectra'] == 3


def test_spread_statistics(collection):
    group = collection.group_by('sample')["PG_myxo"]
    assert np.allclose(group.mean().y, 20.0)
    assert np.allclose(group.median().y, 20.0)
    assert np.allclose(group.std().y, 10.0)
    assert np.allclose(group.sem().y, 10.0 / np.sqrt(3))


# ---------------------------------------------------------------------------
# loading
# ---------------------------------------------------------------------------

def test_from_files_globs_and_names_samples(dpt_tree):
    collection = SpectrumCollection.from_files(
        str(dpt_tree / "*.dpt"), technique="ATR-FTIR")

    assert len(collection) == 5
    assert sorted(set(collection.samples)) == ["PG_coli", "PG_myxo"]
    assert all(s.technique == "ATR-FTIR" for s in collection)
    assert all(len(s) == len(X) for s in collection)


def test_from_files_then_group_is_the_whole_load_step(dpt_tree):
    """What the notebooks spend 40 lines on."""
    averages = {name: group.mean() for name, group in
                SpectrumCollection.from_files(str(dpt_tree / "*.dpt"),
                                              technique="ATR-FTIR")
                .group_by('sample').items()}
    assert np.allclose(averages["PG_coli"].y, 2.0)      # levels 1, 2, 3
    assert np.allclose(averages["PG_myxo"].y, 5.5)      # levels 5, 6


def test_from_files_reports_a_bad_pattern(tmp_path):
    with pytest.raises(FileNotFoundError, match="no files matched"):
        SpectrumCollection.from_files(str(tmp_path / "nothing-here-*.dpt"))


# ---------------------------------------------------------------------------
# batch operations and interop
# ---------------------------------------------------------------------------

def test_batch_operations_return_a_collection():
    """The whole per-sample pipeline, applied across a collection at once."""
    banded = SpectrumCollection([
        _spectrum("a", 1, flat=False), _spectrum("a", 2, flat=False),
        _spectrum("b", 5, flat=False),
    ])
    processed = (banded
                 .crop(1000, 1700)
                 .baseline_correct('rubberband')
                 .normalize('max'))

    assert isinstance(processed, SpectrumCollection)
    assert len(processed) == len(banded)
    assert all(len(s) < len(X) for s in processed)
    # Normalisation removes the intensity difference between replicates.
    for spectrum in processed:
        assert np.isclose(spectrum.y.max(), 1.0)


def test_batch_operations_record_history_on_each(collection):
    processed = collection.crop(1000, 1700).normalize('max')
    assert all([step.name for step in s.history] == ['crop', 'normalize']
               for s in processed)


def test_to_matrix_is_the_handover_to_sklearn(collection):
    x, matrix = collection.to_matrix()
    assert matrix.shape == (len(collection), len(X))
    assert np.allclose(x, X)
    assert np.allclose(matrix[:, 0], [1, 2, 3, 10, 20, 30])


def test_to_matrix_refuses_ragged_spectra(collection):
    ragged = SpectrumCollection(list(collection) + [collection[0].crop(1000, 1500)])
    with pytest.raises(ValueError, match="different lengths"):
        ragged.to_matrix()


def test_resample_makes_a_ragged_collection_stackable(collection):
    ragged = SpectrumCollection(list(collection) + [collection[0].crop(1000, 1500)])
    common_axis = np.linspace(1000.0, 1500.0, 51)
    _, matrix = ragged.resample(common_axis).to_matrix()
    assert matrix.shape == (7, 51)


def test_select_filters(collection):
    coli = collection.select(lambda s: s.metadata['sample'] == 'PG_coli')
    assert len(coli) == 3


def test_slicing_gives_a_collection(collection):
    assert isinstance(collection[:2], SpectrumCollection)
    assert isinstance(collection[0], Spectrum)


def test_rejects_non_spectra():
    with pytest.raises(TypeError, match="Spectrum objects"):
        SpectrumCollection([1, 2, 3])


def test_empty_collection_reduction_is_an_error():
    with pytest.raises(ValueError, match="empty collection"):
        SpectrumCollection().mean()
