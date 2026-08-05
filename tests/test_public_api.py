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
        found |= set(pattern.findall(source.read_text()))

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
