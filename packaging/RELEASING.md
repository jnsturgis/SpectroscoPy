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

This needs owner rights on the PyPI project, so it waits on the name (see
`name-reservation/`). On GitHub, create the environment `pypi` under
*Settings → Environments*. Adding yourself as a required reviewer there gives
a final "approve" click before anything is uploaded.

**If the release goes out as `pyspectroscopy` instead**, change `name` in
`pyproject.toml`. The import name stays `spectroscopy`; only the
`pip install` name changes. Update the README install line to match.

## Each release

1. `main` is green in CI and pushed.
2. `git tag 0.2.0 && git push origin 0.2.0`. There is no `v` prefix, to match `0.1.0`.
3. The workflow builds the sdist and wheel and refuses to continue if the
   version is not exactly the tag. It then runs the test suite against the
   **installed wheel** on Linux and Windows, at Python 3.10 and 3.13, and
   publishes only if all of that passes.

A version on PyPI can never be reused, even after deletion. For a trial run,
tag a release candidate (`0.2.0rc1`) first.
