# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Phase 0 checkpoint tests: the package imports cleanly from anywhere, the layer
boundaries hold, and the deprecated top-level names are gone.
"""

import pathlib
import subprocess
import sys

import pytest


def test_package_imports_without_path_hacks():
    """The whole point of Phase 0: no sys.path.append needed."""
    result = subprocess.run(
        [sys.executable, "-c",
         "import spectroscopy; import spectroscopy.io; "
         "import spectroscopy.processing.ftir; "
         "print(spectroscopy.Spectrum.__name__)"],
        capture_output=True, text=True, cwd="/", check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "Spectrum" in result.stdout


def test_version_is_exposed():
    import spectroscopy
    assert spectroscopy.__version__
    assert spectroscopy.__version__ != "0.0.0+unknown"


def test_importing_package_is_silent():
    """Regression for C5: formats/__init__.py used to print on import."""
    result = subprocess.run(
        [sys.executable, "-c", "import spectroscopy"],
        capture_output=True, text=True, cwd="/", check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout == "", f"import wrote to stdout: {result.stdout!r}"


def test_io_layer_does_not_import_core_at_module_scope():
    """
    Regression for C1: the io layer must not depend on core at import time.

    Checked statically rather than by importing, because importing any
    submodule necessarily runs spectroscopy/__init__.py (which does pull in
    core) -- so sys.modules cannot distinguish the two cases. Walking the AST
    for *module-scope* imports asks the question we actually care about, and
    keeps working as more readers are added.
    """
    import ast
    import pathlib

    import spectroscopy.io

    io_dir = pathlib.Path(spectroscopy.io.__file__).parent
    offenders = []

    for source in sorted(io_dir.glob("*.py")):
        tree = ast.parse(source.read_text(), filename=str(source))
        for node in tree.body:                      # module scope only
            names = []
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                names = [node.module or ""]
            for name in names:
                if name == "spectroscopy" or name.startswith("spectroscopy.spectra"):
                    offenders.append(f"{source.name}: {name}")

    assert not offenders, (
        "io modules import core at module scope, recreating the cycle: "
        + ", ".join(offenders)
    )


@pytest.mark.parametrize("name", ["calc", "formats", "tools_spc"])
def test_the_deprecated_shims_are_gone(name):
    """
    Removed 2026-08-13, roadmap section 14.2 blocker 5.

    They promised removal "in 0.2" and there is no 0.2 -- the next release is
    1.0, and a version whose whole point is keeping its promises cannot be the
    first one to break this one.

    What they were for -- keeping existing notebooks importing -- they had
    stopped doing well before they were removed. The three notebooks that
    import ``formats.jcamp`` call ``jcamp.readfile()``, which is the upstream
    nzhagen API and has never existed in this package: the shim aliased
    ``formats.jcamp`` to ``spectroscopy.io.jcamp``, so the import succeeded and
    the *next* cell raised ``AttributeError``. A shim that makes a failure
    arrive one cell later is worse than no shim, because it looks like support.

    They also squatted on three plausible top-level names in every environment
    that installed this package. ``import formats`` is a thing somebody else's
    project may well want to be.

    Asserted against the source tree rather than against ``import``, because a
    stale copy outlives the repo: ``calc.py`` was ``force-include``-d, so the
    editable install put a real *file* in site-packages and ``import calc``
    kept working after the repo file was gone. That is the squatting problem
    demonstrating itself, and it is an environment to clean rather than a
    reason to weaken the test. What the repo can promise is that it no longer
    ships them.
    """
    root = pathlib.Path(__file__).resolve().parent.parent
    for candidate in (root / f"{name}.py", root / name):
        assert not candidate.exists(), (
            f"{candidate} is back. It was removed deliberately; returning it "
            f"needs a decision rather than a reappearance."
        )

    packaging = (root / "pyproject.toml").read_text()
    packages_line = next(line for line in packaging.splitlines()
                         if line.startswith("packages = "))
    assert name not in packages_line, f"{name!r} is shipped again"


def test_cli_entry_point_is_importable():
    from spectroscopy.cli.ftir_sidechains import main
    assert callable(main)


def test_py_typed_marker_is_present_and_shipped():
    """
    PEP 561: without this file no type checker reads the annotations, however
    many there are, so it is what makes the signatures visible to a user's
    tooling. It is data rather than code, which is exactly the kind of file a
    build backend drops silently -- hence the test.
    """
    import spectroscopy
    marker = pathlib.Path(spectroscopy.__file__).parent / 'py.typed'
    assert marker.is_file(), "spectroscopy/py.typed is missing"
