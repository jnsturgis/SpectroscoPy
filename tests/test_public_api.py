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
    'ProcessingStep',
    'datasets',
    'io',
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
    'collection', 'history', 'messages', 'peaks', 'spectra',
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
    """
    assert len(_fresh_public_names()) <= 15


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
