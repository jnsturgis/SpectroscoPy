#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
One-shot migration for the Phase 1.1 rename:  .x_data -> .x  and  .y_data -> .y

Handles .py files and Jupyter .ipynb notebooks (rewriting the ``source`` of code
cells only, leaving outputs and markdown untouched).

IMPORTANT -- ordering
    The library and the notebooks must flip together. Running this on notebooks
    while spectroscopy.Spectrum still exposes x_data/y_data will break them, and
    vice versa. Run it over the repo and the notebooks in a single pass.

Usage
    python scripts/migrate_rename_xy.py PATH [PATH ...]           # dry run
    python scripts/migrate_rename_xy.py PATH [PATH ...] --apply   # rewrite
    python scripts/migrate_rename_xy.py PATH --apply --no-backup

By default each rewritten file is backed up alongside as ``<name>.bak``.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# Attribute access only: a preceding '.' and a word boundary after. This will
# not match a bare identifier called x_data, nor a string 'x_data'.
PATTERNS = [
    (re.compile(r"\.x_data\b"), ".x"),
    (re.compile(r"\.y_data\b"), ".y"),
]

# Reported but not rewritten: these are almost certainly Spectrum attributes,
# but a human should confirm rather than have a regex guess.
SUSPICIOUS = re.compile(r"""(?<![.\w])['"](x_data|y_data)['"]""")


def rewrite_text(text: str) -> tuple[str, int]:
    """Apply the renames to a blob of source, returning (new_text, n_changes)."""
    total = 0
    for pattern, replacement in PATTERNS:
        text, count = pattern.subn(replacement, text)
        total += count
    return text, total


def process_python(path: Path) -> tuple[str | None, int, list[str]]:
    original = path.read_text(encoding="utf-8")
    updated, count = rewrite_text(original)
    notes = [f"{path.name}: string literal {m.group(0)}"
             for m in SUSPICIOUS.finditer(original)]
    return (updated if count else None), count, notes


def process_notebook(path: Path) -> tuple[str | None, int, list[str]]:
    notebook = json.loads(path.read_text(encoding="utf-8"))
    total = 0
    notes: list[str] = []

    for index, cell in enumerate(notebook.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = cell.get("source", [])
        joined = "".join(source)
        updated, count = rewrite_text(joined)
        if count:
            # Preserve the list-of-lines shape Jupyter writes.
            cell["source"] = updated.splitlines(keepends=True)
            total += count
        notes += [f"{path.name} cell {index}: string literal {m.group(0)}"
                  for m in SUSPICIOUS.finditer(joined)]

    if not total:
        return None, 0, notes
    # indent=1 + trailing newline matches nbformat's own output, keeping the
    # git diff to the lines that actually changed.
    return json.dumps(notebook, indent=1, ensure_ascii=False) + "\n", total, notes


def collect(paths: list[str]) -> list[Path]:
    found: list[Path] = []
    for raw in paths:
        path = Path(raw).expanduser()
        if path.is_dir():
            for pattern in ("**/*.py", "**/*.ipynb"):
                found += [p for p in path.glob(pattern)
                          if ".ipynb_checkpoints" not in p.parts
                          and "site-packages" not in p.parts
                          and ".git" not in p.parts]
        elif path.suffix in {".py", ".ipynb"}:
            found.append(path)
    return sorted(set(found))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("paths", nargs="+", help="Files or directories to migrate.")
    parser.add_argument("--apply", action="store_true",
                        help="Actually rewrite files (default is a dry run).")
    parser.add_argument("--no-backup", action="store_true",
                        help="Do not write <name>.bak copies.")
    args = parser.parse_args(argv)

    targets = collect(args.paths)
    if not targets:
        print("No .py or .ipynb files found.", file=sys.stderr)
        return 1

    grand_total = 0
    touched = 0
    all_notes: list[str] = []

    for path in targets:
        try:
            handler = process_notebook if path.suffix == ".ipynb" else process_python
            new_text, count, notes = handler(path)
        except (OSError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            print(f"  SKIP  {path}: {type(exc).__name__}: {exc}", file=sys.stderr)
            continue

        all_notes += notes
        if not count:
            continue

        touched += 1
        grand_total += count
        print(f"  {count:5d}  {path}")

        if args.apply:
            if not args.no_backup:
                path.with_suffix(path.suffix + ".bak").write_text(
                    path.read_text(encoding="utf-8"), encoding="utf-8")
            path.write_text(new_text, encoding="utf-8")

    verb = "rewrote" if args.apply else "would rewrite"
    print(f"\n{verb} {grand_total} references across {touched} file(s) "
          f"(scanned {len(targets)}).")

    if all_notes:
        print("\nNeeds a human -- 'x_data'/'y_data' appearing as string literals,\n"
              "left untouched because they may be dict keys or column names:")
        for note in all_notes:
            print(f"  {note}")

    if not args.apply:
        print("\nDry run. Re-run with --apply to write changes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
