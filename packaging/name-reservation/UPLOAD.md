# Reserving `pyspectroscopy` — how to upload

This is a placeholder: version 0.0.0, one module containing a docstring that
says where the real library is, "Development Status :: 1 - Planning", and a
README that tells anyone who installs it not to depend on it.

## Steps

0. **Build first.** `dist/` is not committed — it is in `.gitignore`, like
   every other build artefact — so on a fresh checkout it does not exist and
   `twine upload dist/*` fails with

       InvalidDistribution: Cannot find file (or expand pattern): 'dist/*'

   That error means "nothing has been built yet", not "something is wrong
   with your account or token". From this directory:

       pip install --user build twine
       python -m build
       python -m twine check dist/*

   which should leave you with

       pyspectroscopy-0.0.0-py3-none-any.whl
       pyspectroscopy-0.0.0.tar.gz

   ⚠️ Build **here**, in `packaging/name-reservation/`, not at the repository
   root. The root builds the real library under the name `spectroscopy`,
   which belongs to somebody else on PyPI — uploading it would be rejected,
   and it is not what you want to publish yet in any case.

1. **Register on PyPI** if you have not: https://pypi.org/account/register/
   Enable 2FA — PyPI requires it for uploading.

2. **Create an API token**: https://pypi.org/manage/account/token/
   Scope it to "Entire account" for the first upload (there is no project to
   scope to yet). Afterwards, create a project-scoped token and delete the
   account-wide one.

3. **Look inside the sdist before uploading it.** Not optional — this went
   wrong once already:

       tar tzf dist/*.tar.gz

   You should see exactly four entries: the `pyspectroscopy/` module,
   `README.md`, `pyproject.toml`, `PKG-INFO` (and a harmless `.gitignore`
   hatchling insists on). **Anything else, stop.**

   `twine check` does *not* do this. It validates that the metadata renders;
   it has no opinion about a credentials file in the archive.

   Version 0.0.0 shipped `PyPI-Recovery-Codes-….txt` because the file happened
   to be sitting in this directory when the build ran, and hatchling's default
   is to include everything. `pyproject.toml` now carries an `only-include`
   whitelist so that new files are excluded unless named — but check anyway,
   because the whole failure mode is a file nobody thought about.

4. **Upload**, from this directory:

       python -m twine upload dist/*

   Username: `__token__` — literally that, underscores included. Not your
   PyPI username; that is the commonest way this step fails.
   Password: the token, `pypi-…` prefix included.

   Keep the token out of the shell history and out of this repository. If you
   want it stored, `~/.pypirc` with mode 600 is the usual place:

       [pypi]
         username = __token__
         password = pypi-…

       chmod 600 ~/.pypirc

   Or test the mechanics first against TestPyPI, which is a separate registry
   and does not reserve the real name:

       python -m twine upload --repository testpypi dist/*

5. **Check** https://pypi.org/project/pyspectroscopy/ renders sensibly, and
   that the file list holds only what you expect.

## If something private ever gets published again

In order, and the order is the point:

1. **Rotate the secret first.** Regenerate the recovery codes, revoke the
   token, change the password — whatever it was. PyPI files are mirrored,
   cached and scraped within minutes; removing the file does not unpublish it,
   and treating deletion as the fix leaves a live credential in the wild.
2. **Delete the affected *release*** — Manage → Releases → the version →
   Options → Delete. Yanking is *not* enough: a yanked file is still
   downloadable, it merely stops resolvers picking it by default.
3. ⚠️ **Never delete the project.** That releases the *name*, which PyPI will
   not let you re-register. Deleting `pyspectroscopy` would destroy the very
   thing this placeholder exists to hold.
4. **Bump the version and re-upload.** A deleted version number can never be
   reused, so go forward: 0.0.0 → 0.0.1.

## What this does and does not buy

It stops the fallback name being taken while the `spectroscopy` request sits,
which is the whole point. It does **not** commit you to the name: if the PEP 541
request succeeds you simply publish the real library as `spectroscopy` and leave
this at 0.0.0, or delete it.

One caveat worth knowing before you press the button: **PyPI does not let a
deleted project name be re-registered.** Deleting `pyspectroscopy` later would
burn the name rather than free it. That is an argument for leaving the
placeholder in place indefinitely rather than tidying it up — it costs nothing.

## Do not publish the real library under two names

If `spectroscopy` comes through, publish there and leave this as a placeholder.
Two live packages containing the same library is how you get a citation pointing
at one and users installing the other — precisely the problem the name decision
exists to avoid.
