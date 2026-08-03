# Contributing and feedback

## Telling us something is wrong

**You do not need a GitHub account.** Every page here has a feedback link at
the bottom that opens an email. Use it — a report that arrives beats a report
you did not send because the tooling was in the way.

If you do have an account, the
[issue tracker](https://github.com/jnsturgis/SpectroscoPy/issues) is better,
because other people can find it.

What helps most:

- What you ran, what happened, and what you expected instead.
- The traceback, if there was one.
- **The file, if a file is involved.** Most of the hard bugs in this project
  have been in file readers, and almost none were reproducible without the
  actual file. A file that will not read is a genuinely valuable contribution
  — the OPUS and SPC readers both exist because somebody supplied files. If
  the data are unpublished, a truncated or anonymised version is usually
  enough.
- Your operating system, and:

```python
import spectroscopy
print(spectroscopy.__version__)
```

Unclear error messages count as bugs. They stop people using a tool more
effectively than missing features do.

## Working on the code

The full guide lives in
[CONTRIBUTING.md](https://github.com/jnsturgis/SpectroscoPy/blob/main/CONTRIBUTING.md)
— setup, the test suite, optional fixtures, how to add a file format, and the
licensing constraints on where code may be copied from.

The short version:

```bash
git clone https://github.com/jnsturgis/SpectroscoPy.git
cd SpectroscoPy
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

pytest -q
ruff check .
```

Three things about this codebase that are worth knowing before you send a
patch:

**Explain why, not what.** The code says what it does. Comments exist for the
reasoning that cannot be recovered by reading it — why this threshold, why
this order, what went wrong when it was done the obvious way.

**Do not tune a number until the answer looks right.** If a method gives a
wrong answer on real data, that is the finding. There is a worked example in
the roadmap: an FTIR secondary-structure estimator with excellent fit quality
and a twenty-percentage-point error, recorded as not usable rather than
adjusted until it agreed with expectation. Fit quality is not validation.

**Decisions that are hard to reverse get written down.** Anything touching the
public API, the on-disk format, or a default that changes numerical results
wants an [architecture decision record](adr/index.md) first — the alternatives
and why they lost matter more than the choice.

## Conduct

Participation is covered by the
[Code of Conduct](https://github.com/jnsturgis/SpectroscoPy/blob/main/CODE_OF_CONDUCT.md)
(Contributor Covenant 2.1). One point follows from this being a scientific
project: reporting that a method gives a wrong answer, including one written
by somebody here, is a contribution and not an attack.

## Citing

See
[CITATION.cff](https://github.com/jnsturgis/SpectroscoPy/blob/main/CITATION.cff),
or GitHub's "Cite this repository" button. Cite the repository URL rather than
a `pip` name for the moment — the PyPI name is not yet settled, and the URL is
stable whichever way that goes.
