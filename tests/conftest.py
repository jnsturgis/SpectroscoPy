# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Shared fixtures."""

from pathlib import Path

import pytest

#: Vendor sample files that are deliberately not committed -- see
#: scripts/fetch_spc_fixtures.py for why, and for how to get them.
SPC_SAMPLE_DIR = Path(__file__).parent / 'data' / 'spc'

_MISSING = (
    "Galactic sample .spc files are not present. They are not committed "
    "(unstated licence); fetch them with:\n"
    "    python scripts/fetch_spc_fixtures.py"
)


@pytest.fixture(scope='session')
def spc_samples():
    """Directory of Galactic sample files, skipping the test when absent."""
    if not SPC_SAMPLE_DIR.is_dir() or not any(SPC_SAMPLE_DIR.iterdir()):
        pytest.skip(_MISSING)
    return SPC_SAMPLE_DIR


@pytest.fixture(scope='session')
def spc_sample(spc_samples):
    """Look one sample up by name, skipping if that particular file is absent."""
    def get(name):
        path = spc_samples / name
        if not path.exists():
            pytest.skip(f"{name} not present. {_MISSING}")
        return path
    return get
