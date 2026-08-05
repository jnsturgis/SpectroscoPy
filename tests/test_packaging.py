# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Phase 0 checkpoint tests: the package imports cleanly from anywhere, the layer
boundaries hold, and the deprecated top-level names still work for existing
notebooks.
"""

import pathlib
import subprocess
import sys
import warnings

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


@pytest.mark.parametrize("statement", [
    "import calc; calc.gauss",
    "import formats; formats.jcamp",
    "from formats import jcamp; jcamp.read",
    "import formats.jcamp; formats.jcamp.read",
])
def test_deprecated_shims_still_work(statement):
    """Existing notebooks must keep running through the 0.1 series."""
    result = subprocess.run([sys.executable, "-c", statement],
                            capture_output=True, text=True, cwd="/", check=False)
    assert result.returncode == 0, result.stderr


def test_deprecated_shims_warn():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        import importlib

        import formats
        importlib.reload(formats)
    assert any(issubclass(w.category, DeprecationWarning) for w in caught)


def test_formats_and_io_share_module_objects():
    """The shim must alias, not duplicate, the reader modules."""
    import formats
    import spectroscopy.io
    assert formats.jcamp is spectroscopy.io.jcamp
    assert formats.csv is spectroscopy.io.csv
    assert formats.spy is spectroscopy.io.spy


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
