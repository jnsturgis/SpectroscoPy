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


def _module_scope_import_graph():
    """
    Which module imports which, counting *module-scope* imports only.

    Checked statically rather than by importing, because importing any
    submodule necessarily runs spectroscopy/__init__.py (which pulls in most
    of the package) -- so sys.modules cannot distinguish an import-time
    dependency from a deferred one. Walking the AST asks the question we
    actually care about, and keeps working as more modules are added.

    Imports inside a function body are deliberately not edges: deferring one
    is the standard way to break a cycle, and the point of this graph is to
    find the cycles nobody has broken.
    """
    import ast

    import spectroscopy

    root = pathlib.Path(spectroscopy.__file__).parent
    sources, packages = {}, {}
    for source in sorted(root.rglob("*.py")):
        parts = source.relative_to(root.parent).with_suffix("").parts
        # A package's __init__.py *is* the package; anything else is a module
        # inside its parent package. The distinction decides what a relative
        # import counts from, so `from . import x` resolves differently in
        # io/__init__.py than in io/registry.py.
        if parts[-1] == "__init__":
            parts = parts[:-1]
            packages[".".join(parts)] = ".".join(parts)
        else:
            packages[".".join(parts)] = ".".join(parts[:-1])
        sources[".".join(parts)] = source

    def owning_module(target):
        while target and target not in sources:
            target = target.rsplit(".", 1)[0] if "." in target else None
        return target

    graph = {}
    for name, source in sources.items():
        tree = ast.parse(source.read_text(encoding='utf-8'), filename=str(source))
        targets = set()
        for node in tree.body:                      # module scope only
            if isinstance(node, ast.Import):
                targets.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                if node.level:                      # `from . import x`
                    base = packages[name].rsplit(".", node.level - 1)[0]
                    module = f"{base}.{node.module}" if node.module else base
                else:
                    module = node.module or ""
                targets.add(module)
                targets.update(f"{module}.{a.name}" for a in node.names)
        graph[name] = {
            owner for target in targets
            if target.startswith("spectroscopy")
            and (owner := owning_module(target)) and owner != name
        }
    return graph


def test_no_module_scope_import_cycles():
    """
    Regression for C1, generalised.

    C1 was `core` reaching down into `io`: spectra.py imported the readers and
    dispatched on file type. The registry inverted that edge, so the package
    now runs one way -- io imports core, core does not import io -- and
    io/registry.py names Spectrum at its top like any other dependency.

    The guard this replaced forbade io from importing core at module scope,
    which was the *avoidance mechanism* of the era when the edge still pointed
    both ways, not the invariant. It would have passed a cycle routed through
    spectroscopy.collection, which it did not check for. Asking directly
    whether a cycle exists is both stricter and less likely to outlive its
    reason -- and a cycle here does not merely offend a diagram, it makes the
    package's importability depend on which module is named first.
    """
    graph = _module_scope_import_graph()

    cycles, visiting, done = [], [], set()

    def walk(name, path):
        path.append(name)
        visiting.append(name)
        for target in sorted(graph.get(name, ())):
            if target in visiting:
                cycles.append(path[path.index(target):] + [target])
            elif target not in done:
                walk(target, path)
        visiting.pop()
        path.pop()
        done.add(name)

    for name in sorted(graph):
        if name not in done:
            walk(name, [])

    assert not cycles, "module-scope import cycles: " + "; ".join(
        " -> ".join(cycle) for cycle in cycles
    )


@pytest.mark.parametrize("first", [
    "spectroscopy", "spectroscopy.spectra", "spectroscopy.collection",
    "spectroscopy.io.registry", "spectroscopy.io.spy", "spectroscopy.library",
    "spectroscopy.processing.structure",
])
def test_any_module_can_be_imported_first(first):
    """
    A cycle broken by a deferred import still fails if the deferral is undone,
    and the failure looks like 'partially initialized module' in whichever
    module the user happened to name first. Importing each one into a fresh
    interpreter is the check a static graph cannot make.
    """
    result = subprocess.run(
        [sys.executable, "-c", f"import {first}"],
        capture_output=True, text=True, cwd="/", check=False,
    )
    assert result.returncode == 0, result.stderr


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

    packaging = (root / "pyproject.toml").read_text(encoding="utf-8")
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
