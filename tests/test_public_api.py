# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""The public surface, pinned.

1.0 promises not to break what it exports, so what it exports has to be a
decision rather than an accident. Before this test the package re-exported
whatever ``spectra.py`` happened to import -- ``spectroscopy.os``,
``spectroscopy.np``, ``spectroscopy.copy``, forty names in all, seventeen of
them modules. Freezing that would have promised to keep ``spc.np`` forever.

These tests fail when a name appears or disappears without the list below being
edited, which is the point: adding to the public API should take a deliberate
line in a diff. Roadmap section 14.2.
"""

import types

import pytest

import spectroscopy as spc

#: Exactly what ``from spectroscopy import *`` should bind.
EXPECTED_ALL = {
    'Spectrum',
    'SpectrumCollection',
    'PeakTable',
    'FitResult',
    'ProcessingStep',
    'datasets',
    'io',
    'library',
    'metadata',
    'processing',
    'units',
    'lineshapes',
    'viz',
}

#: Submodules bound on the package as a side effect of the imports in
#: ``__init__`` and below it. They cannot be hidden -- ``import
#: spectroscopy.spectra`` works whatever we do -- but they are not in __all__
#: and 1.0 does not promise them.
TOLERATED_SUBMODULES = {
    'collection', 'fitting', 'history', 'messages', 'peaks', 'spectra',
}


def _fresh_public_names():
    """What a plain ``import spectroscopy`` binds, in a clean interpreter.

    Not ``dir(spc)`` in this process: importing any submodule anywhere binds it
    as an attribute of the package, so by the time the suite reaches this file
    another test has pulled in ``spectroscopy.cli`` and the answer depends on
    test ordering. A subprocess is what a user actually sees.
    """
    import subprocess
    import sys

    code = ('import spectroscopy; '
            "print(' '.join(n for n in dir(spectroscopy) "
            "if not n.startswith('_')))")
    out = subprocess.run([sys.executable, '-c', code],
                         capture_output=True, text=True, check=True)
    return set(out.stdout.split())


def test_all_is_exactly_the_expected_surface():
    assert set(spc.__all__) == EXPECTED_ALL


def test_every_exported_name_resolves():
    missing = [n for n in spc.__all__ if not hasattr(spc, n)]
    assert not missing, f"__all__ promises names that do not exist: {missing}"


def test_star_import_binds_only_the_public_surface():
    namespace = {}
    exec('from spectroscopy import *', namespace)          # noqa: S102
    bound = {n for n in namespace if not n.startswith('__')}
    assert bound == EXPECTED_ALL


def test_no_foreign_modules_are_reachable():
    """``spectroscopy.np`` and friends must stay gone.

    This is the specific regression: a bare ``from .spectra import *`` with no
    ``__all__`` in spectra.py re-exports numpy, os, copy, operator, warnings and
    scipy's CubicSpline.
    """
    foreign = [
        name for name in dir(spc)
        if not name.startswith('_')
        and isinstance(getattr(spc, name), types.ModuleType)
        and not getattr(spc, name).__name__.startswith('spectroscopy')
    ]
    assert not foreign, f"third-party or stdlib modules leaked into the API: {foreign}"


@pytest.mark.parametrize('name', ['np', 'os', 'copy', 'operator', 'warnings',
                                  'CubicSpline', 'PackageNotFoundError'])
def test_known_leaks_stay_fixed(name):
    assert not hasattr(spc, name), f"{name} is back in the public namespace"


def test_nothing_unexpected_is_public():
    public = _fresh_public_names()
    # viz is lazy, so a fresh import does not bind it until it is asked for.
    unexpected = public - EXPECTED_ALL - TOLERATED_SUBMODULES - {'viz'}
    assert not unexpected, (
        f"new public names appeared: {sorted(unexpected)}. Add them to "
        "EXPECTED_ALL deliberately, or make them private."
    )


def test_the_surface_is_no_larger_than_it_looks():
    """A blunt count, so a slow drift upward is visible in a diff.

    It was 40 names, 17 of them modules, before roadmap section 14.2.
    Raised to 18 on 2026-08-05 for ``metadata``, which documents the keys that
    freeze with the ``.spy`` format (roadmap D2).
    """
    assert len(_fresh_public_names()) <= 18


def test_viz_is_not_imported_eagerly():
    """matplotlib costs ~350 ms; importing spectroscopy should not pay it.

    Checked through a fresh interpreter, because by the time this test runs
    another test has almost certainly imported viz already.
    """
    import subprocess
    import sys

    code = (
        'import sys, spectroscopy; '
        "print('matplotlib.pyplot' in sys.modules)"
    )
    out = subprocess.run([sys.executable, '-c', code],
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False', "importing spectroscopy pulled in pyplot"


def test_viz_is_still_reachable_without_a_separate_import():
    assert spc.viz.plot is not None


# ---------------------------------------------------------------------------
# metadata keys freeze with the .spy format, whether or not anyone says so:
# the writer serialises the dictionary verbatim. Roadmap D2.
# ---------------------------------------------------------------------------

def test_the_metadata_keys_the_library_reads_are_pinned():
    """
    These names are in every .spy file ever written. Renaming one silently
    orphans the data in existing files -- the key is still there, and nothing
    reads it any more.
    """
    from spectroscopy import metadata

    assert set(metadata.KNOWN_KEYS) == {
        # sample conditions
        'path_length', 'concentration', 'mass_concentration',
        'n_residues', 'mean_residue_weight', 'temperature', 'pH',
        'reference_electrode',
        # identification
        'sample', 'reference', 'spec_type',
        'parameter', 'parameter_name', 'parameter_unit',
        # acquisition
        'excitation_nm', 'z_value', 'z_quantity', 'scans',
    }


def test_the_keys_the_code_actually_reads_are_in_the_schema():
    """
    The schema is only worth having if it describes the real thing. This
    catches a module that starts reading a key nobody wrote down -- which is
    how path_length arrived in the first place.
    """
    import pathlib
    import re

    root = pathlib.Path(__file__).resolve().parent.parent / 'spectroscopy'
    from spectroscopy import metadata

    pattern = re.compile(r"metadata(?:\[|\.get\()'([a-zA-Z_]+)'")
    found = set()
    for source in root.rglob('*.py'):
        if source.name in ('metadata.py', 'jcamp.py'):
            continue                       # the schema itself; vendored parser
        found |= set(pattern.findall(source.read_text(encoding='utf-8')))

    undocumented = {key for key in found
                    if key not in metadata.KNOWN_KEYS
                    and not key.startswith(metadata.PROVENANCE_PREFIXES)}
    assert not undocumented, (
        f"these metadata keys are read by the code but are not in the "
        f"schema: {sorted(undocumented)}. Add them to spectroscopy/metadata.py "
        f"or give them a provenance prefix."
    )


def test_unknown_keys_finds_the_near_miss():
    """'pathlength' is ignored by every consumer and looks exactly like
    having forgotten to set it."""
    from spectroscopy import metadata

    assert metadata.unknown_keys({'path_length': 1.0}) == []
    assert metadata.unknown_keys({'pathlength': 1.0}) == ['pathlength']
    assert metadata.unknown_keys({'opus_history': '...'}) == []


# ---------------------------------------------------------------------------
# Guessability: the wrong guess should work, or fail loudly -- never work and
# be wrong. SpectroscoPy_API_Guessability.md, applied 2026-08-05.
# ---------------------------------------------------------------------------

def test_assigning_technique_does_what_set_type_does():
    """
    A1, the one that was silently wrong: `.technique` was a plain attribute,
    so the guessable line left an infrared spectrum labelled in nm and reading
    it back reported exactly what had been asked for.
    """
    import numpy as np

    def fresh():
        return spc.Spectrum(np.linspace(900, 1800, 50), np.ones(50))

    assigned, called = fresh(), fresh()
    assigned.technique = 'ATR-FTIR'
    called.set_type('ATR-FTIR')

    for attribute in ('technique', 'x_unit', 'x_quantity', 'y_unit'):
        assert getattr(assigned, attribute) == getattr(called, attribute)
    assert assigned.metadata['spec_type'] == 'ATR-FTIR'
    assert assigned.x_unit == 'cm^-1'


def test_assigning_an_unknown_technique_is_refused():
    with pytest.raises(TypeError, match='Unknown spectrum type'):
        spc.Spectrum().technique = 'NMR'


def test_technique_can_still_be_cleared():
    spectrum = spc.Spectrum()
    spectrum.set_type('UV-Vis')
    spectrum.technique = None
    assert spectrum.technique is None
    assert 'spec_type' not in spectrum.metadata


def test_copying_keeps_the_axes_the_original_had():
    """
    The property setter applies the technique's default axes, so the copy
    constructor must not go through it -- a UV-Vis spectrum deliberately held
    in cm^-1 would come back in nm.
    """
    import numpy as np

    original = spc.Spectrum(np.linspace(400, 700, 10), np.ones(10),
                            technique='UV-Vis')
    original.x_unit = 'cm^-1'
    assert spc.Spectrum(original).x_unit == 'cm^-1'


@pytest.mark.parametrize('alias, canonical', [
    ('get_info', 'describe'),
    ('normalise', 'normalize'),
    ('write', 'save_as'),
])
def test_spectrum_aliases_exist(alias, canonical):
    assert hasattr(spc.Spectrum, alias) and hasattr(spc.Spectrum, canonical)


@pytest.mark.parametrize('alias, canonical', [
    ('groupby', 'group_by'),
    ('filter', 'select'),
    ('normalise', 'normalize'),
])
def test_collection_aliases_exist(alias, canonical):
    assert (hasattr(spc.SpectrumCollection, alias)
            and hasattr(spc.SpectrumCollection, canonical))


def test_aliases_return_what_the_canonical_names_do():
    import numpy as np

    x = np.linspace(900, 1800, 50)
    spectra = [spc.Spectrum(x, np.full_like(x, level), technique='ATR-FTIR')
               for level in (1.0, 2.0)]
    for index, spectrum in enumerate(spectra):
        spectrum.set_sample(f"s{index}")
    collection = spc.SpectrumCollection(spectra)

    assert set(collection.groupby('sample')) == set(
        collection.group_by('sample'))
    assert len(collection.filter(lambda s: True)) == len(collection)
    assert spectra[0].describe() == spectra[0].get_info()


def test_the_mutating_and_copying_prefixes_are_kept_apart():
    """
    A2: `set_` mutates and returns None on a Spectrum; the collection method
    returns a new collection, so it is `with_parameters` and not `set_`.
    """
    assert hasattr(spc.SpectrumCollection, 'with_parameters')
    assert not hasattr(spc.SpectrumCollection, 'set_parameters')
    assert spc.Spectrum().set_sample('x') is None


def test_dropped_aliases_stay_dropped():
    """
    Both borrowed a pandas name and then broke its contract, which is worse
    than not having the alias: pandas' to_numpy returns an ndarray and ours
    returned a 2-tuple; pandas' nlargest ranks by value and ours returned
    position order. A guess that works and returns the wrong shape or order
    is the failure this whole audit is about. Added and removed 2026-08-05.
    """
    assert not hasattr(spc.SpectrumCollection, 'to_numpy')
    assert not hasattr(spc.PeakTable, 'nlargest')


def test_strongest_returns_them_strongest_first():
    """It returned position order until 2026-08-05, despite the name."""
    import numpy as np

    table = spc.PeakTable(position=np.array([100.0, 200.0, 300.0]),
                          height=np.array([0.5, 0.9, 0.7]),
                          index=np.arange(3))
    assert list(table.strongest(3).height) == [0.9, 0.7, 0.5]
    assert list(table.strongest(2).position) == [200.0, 300.0]
    # and the tutorials' explicit re-sort is meaningful again
    assert list(table.strongest(3).sorted_by_position().position) == [
        100.0, 200.0, 300.0]


def test_a_dpt_file_knows_it_is_infrared():
    """
    899-3998 cm^-1 labelled 'Wavelength (nm)' plots as a mislabelled mirror
    image, and near-infrared nanometres are plausible enough that nothing
    looks wrong.
    """
    import pathlib as _pathlib

    root = _pathlib.Path(__file__).resolve().parent.parent
    spectrum = spc.Spectrum.read(
        root / 'spectroscopy/data/ftir_replicates/Glucose.1.dpt')
    assert spectrum.technique == 'FTIR'
    assert spectrum.x_unit == 'cm^-1'
    assert spectrum.x_quantity == 'Wavenumber'
    assert spectrum.reversed_x
    # and it stays refinable to the sampling accessory the file cannot know
    spectrum.set_type('ATR-FTIR')
    assert spectrum.x_unit == 'cm^-1'


def test_a_coefficient_says_which_way_round_it_is():
    """.value is the reciprocal of the number people quote."""
    dsdna = spc.library.coefficient('dsDNA', 260)
    assert dsdna.value == pytest.approx(0.02)
    assert dsdna.quoted_as == pytest.approx(50.0)
    assert '50 ug/mL' in str(dsdna)
