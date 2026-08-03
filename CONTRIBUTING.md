# Contributing to SpectroscoPy

Bug reports from people using this on real data are worth more than almost any
patch, so start there if that is what you have.

## Reporting something

**You do not need a GitHub account.** Every page of the
[documentation](https://jnsturgis.github.io/SpectroscoPy/) carries a feedback
link that opens an email to the maintainer, and that is a perfectly good way
to report a problem. Use it if it is easier — a report that arrives is better
than a report you did not send.

If you do have an account, the
[issue tracker](https://github.com/jnsturgis/SpectroscoPy/issues) is better,
because other people can find it.

Either way, the useful things to include are:

- What you ran and what happened, ideally with the traceback.
- What you expected instead.
- **The file, if a file is involved.** Most of the hard bugs in this project
  have been in file readers, and almost none of them were reproducible without
  the actual file. A file that fails to read is a genuinely valuable
  contribution — the OPUS and SPC readers both exist because somebody supplied
  files.
- Your operating system and the output of `python -c "import spectroscopy;
  print(spectroscopy.__version__)"`.

If the data are unpublished or sensitive, say so and send a truncated or
anonymised version — the first few hundred points of a spectrum are usually
enough to diagnose a reader.

## Setting up to work on the code

```bash
git clone https://github.com/jnsturgis/SpectroscoPy.git
cd SpectroscoPy
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

A virtual environment is worth the two extra lines: many Linux distributions
now ship a Python that refuses `pip install` into the system directories
(PEP 668), and `--user` installs will fight with your package manager sooner
or later.

Then:

```bash
pytest -q                 # the test suite
ruff check .              # lint, same invocation CI uses
python -m sphinx -b html -W --keep-going docs docs/_build/html
```

All three run in CI, and the documentation build uses `-W`, so a warning is a
failure there even though it is not locally.

### Optional test fixtures

Some tests need files that are not in the repository because their licence
does not permit redistribution. They skip cleanly when absent, so the suite
passes either way. To run them:

```bash
python scripts/fetch_spc_fixtures.py
```

That fetches Galactic's own sample `.spc` files from the Internet Archive into
`tests/data/spc/`, which is gitignored. Without them you get 10 passes and 18
skips in `tests/test_io_spc.py`; with them, 28 passes.

## What a change should look like

**Explain why, not what.** The code says what it does. Comments and commit
messages exist for the reasoning that is not recoverable from reading it —
why this threshold, why this order, what went wrong when it was done the
obvious way. Much of this codebase's commentary records a measurement or a
failed approach, and that is deliberate.

**Add a test that would have failed before.** For a bug fix, write the test
first and watch it fail; a test that passes before and after is not testing
what you think.

**Do not tune a number until the answer looks right.** If a method gives a
wrong answer on real data, the finding is that it gives a wrong answer. There
is a worked example of this in roadmap §19–20: an FTIR secondary-structure
estimator with excellent fit quality and a twenty-percentage-point error,
which was recorded as not usable rather than adjusted until it agreed with
expectation. Fit quality is not validation.

**Match the surrounding code.** Naming, comment density, docstring style. New
code should be hard to pick out.

### Decisions that are hard to reverse

Anything that changes the public API, the on-disk format, or a default that
alters numerical results wants an
[ADR](https://jnsturgis.github.io/SpectroscoPy/adr/) — a short record of the
decision and, more importantly, of the alternatives and why they lost. See
`docs/adr/` for the two that exist. Open an issue before writing code for
these; it is much cheaper to disagree about an approach than about a
pull request.

### Adding a file format

The most likely useful contribution, and the path is short:

1. A module in `spectroscopy/io/` that registers itself with
   `@register_reader(...)`. There is no dispatch table to update — the
   registry is the single source of truth, and `describe_formats()` picks it
   up automatically.
2. Tests, including at least one real file if you can supply one.
3. A row in the format tables in `README.md`, `docs/index.md` and
   `docs/guide/reading-files.md`.

Two conventions worth knowing. Binary formats register with `binary=True`, so
the registry does not try to sniff a text encoding and corrupt them. And a
format whose extension is ambiguous should check the file's own magic value
rather than trusting the name — `.spc` is used by both Galactic and Bruker
EPR, and `.0` says nothing at all.

**We do not write undocumented formats.** Reading a reverse-engineered format
is fine; writing one is not, because nothing can establish that the file we
produced is one the vendor's software would accept. `spectroscopy/io/README.md`
has the reasoning.

## Pull requests

Branch from `main`, keep the change focused, and say in the description what
problem it solves. If it is your first contribution, add yourself to the
authors list in `CITATION.cff` in the same PR.

CI must be green: lint, tests on the supported Python versions, and the
documentation build.

## Licence

SpectroscoPy is [MPL-2.0](LICENSE). Contributions are accepted under the same
licence, and by opening a pull request you agree to that. There is no CLA.

MPL-2.0 is file-level copyleft: it can be combined with proprietary code, but
modifications to these files stay open. Note this means **GPL-licensed code
cannot be copied into this project** — the compatibility runs the other way.
If you are porting an algorithm, port it from a specification or from
permissively licensed code, and say in the PR where it came from.

## Conduct

By taking part you agree to the [Code of Conduct](CODE_OF_CONDUCT.md).
