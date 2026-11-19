# Releasing to PyPI

Publishing is done by `.github/workflows/release.yml` when a version tag is
pushed. Nobody uploads from a laptop, and no token exists to leak.

## Once, before the first release

On PyPI, under the project's *Manage → Publishing*, add a trusted publisher:

| Field | Value |
|---|---|
| Owner | `jnsturgis` |
| Repository | `SpectroscoPy` |
| Workflow | `release.yml` |
| Environment | `pypi` |

This needs the Owner role on the PyPI project. On GitHub, create the environment `pypi` under
*Settings → Environments*. Adding yourself as a required reviewer there gives
a final "approve" click before anything is uploaded.

## Each release

Numbering starts at **0.20.0**, never below. Joey Knightbrook's `0.10` stays on
the project, as promised. Versions compare part by part, so 0.2 to 0.9 all sort
below 0.10, and `pip install spectroscopy` would install his package instead.
The release workflow refuses any tag below 0.20.0. **Never yank or delete 0.10.**

1. `version` and `date-released` in `CITATION.cff` updated, committed.
2. `main` is green in CI and pushed.
3. `git tag 0.20.0 && git push origin 0.20.0`. There is no `v` prefix, to match `0.1.0`.
4. The workflow builds the sdist and wheel and refuses to continue if the
   version is not exactly the tag. It then runs the test suite against the
   **installed wheel** on Linux and Windows, at Python 3.10 and 3.13, and
   publishes only if all of that passes.

A version on PyPI can never be reused, even after deletion. For a trial run,
tag a release candidate (`0.20.0rc1`) first.
