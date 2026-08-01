# Reserving `pyspectroscopy` — how to upload

Built and ready in `dist/`:

    pyspectroscopy-0.0.0-py3-none-any.whl
    pyspectroscopy-0.0.0.tar.gz

This is a placeholder: version 0.0.0, one module containing a docstring that
says where the real library is, "Development Status :: 1 - Planning", and a
README that tells anyone who installs it not to depend on it.

## Steps

1. **Register on PyPI** if you have not: https://pypi.org/account/register/
   Enable 2FA — PyPI requires it for uploading.

2. **Create an API token**: https://pypi.org/manage/account/token/
   Scope it to "Entire account" for the first upload (there is no project to
   scope to yet). Afterwards, create a project-scoped token and delete the
   account-wide one.

3. **Upload**, from this directory:

       pip install --user twine
       python -m twine upload dist/*

   Username: `__token__`
   Password: the token, `pypi-…` included.

   Or test the mechanics first against TestPyPI, which is a separate registry
   and does not reserve the real name:

       python -m twine upload --repository testpypi dist/*

4. **Check** https://pypi.org/project/pyspectroscopy/ renders sensibly.

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
