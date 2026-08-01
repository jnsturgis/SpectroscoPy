# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Command-line entry points.

These are thin wrappers: argument parsing, file reading and plotting only.
All computation lives in ``spectroscopy.processing``, so that a CLI can do
nothing a notebook user could not also do through the library API.
"""
