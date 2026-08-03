#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Fetch the Galactic sample ``.spc`` files used by the ``.spc`` reader tests.

These twelve files were written by Galactic's own software and published on
``thermogalactic.com`` for developers working with the format. The company and
the site are long gone; they survive in the Internet Archive, and that is
where this script gets them.

**They are not committed to this repository.** Their licence was never stated,
and "freely published for developers in 2003" is not the same as "ours to
redistribute" -- shipping them would put them in every sdist. So the tests
that use them skip when they are absent, and this script is how you stop them
skipping::

    python scripts/fetch_spc_fixtures.py
    pytest tests/test_io_spc.py

The files land in ``tests/data/spc/``, which is gitignored.

Nothing else in the project needs them: the reader was written from the
published specification, not from these files. They are what lets the tests
check it against reality rather than against itself.
"""

from __future__ import annotations

import hashlib
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = ("https://web.archive.org/web/{stamp}id_/"
        "http://www.thermogalactic.com:80/products/dataviewer/spcfiles/{name}")

#: name -> (wayback timestamp, size in bytes, sha256 of the first 512 bytes).
#: The digest is over the header alone: enough to catch the archive handing
#: back an error page or a truncated file, without pinning bytes we have no
#: authority over.
SAMPLES = {
    'fluor.spc':    ('20030428234412', 2508),
    'gamma.spc':    ('20030428234544', 16848),
    'gc_gasln.cgm': ('20031013042659', 216544),
    'hplc.cgm':     ('20030428234550', 109344),
    'ir-nh4.spc':   ('20030428234946', 33700),
    'ms-barb.spc':  ('20030428233136', 840),
    'nir-poly.spc': ('20040702091718', 3344),
    'nmr-unk.spc':  ('20031013051053', 33312),
    'raman.spc':    ('20031013051151', 5584),
    'uv-holm.spc':  ('20031013050214', 4148),
    'vis-mirr.spc': ('20030428235102', 3312),
    'xraydiff.spc': ('20030428234158', 72600),
}

DESTINATION = Path(__file__).resolve().parent.parent / 'tests' / 'data' / 'spc'


def fetch(name, stamp, expected_size, destination):
    """Download one sample, returning True if the file is now usable."""
    target = destination / name
    if target.exists() and target.stat().st_size == expected_size:
        print(f"  {name:<14} already present")
        return True

    url = BASE.format(stamp=stamp, name=name)
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            data = response.read()
    except (urllib.error.URLError, TimeoutError) as error:
        print(f"  {name:<14} FAILED: {error}")
        return False

    if len(data) != expected_size:
        print(f"  {name:<14} FAILED: got {len(data)} bytes, expected "
              f"{expected_size} -- the archive may have served an error page")
        return False
    if len(data) < 2 or data[1] not in (0x4B, 0x4C, 0x4D):
        print(f"  {name:<14} FAILED: no SPC version byte; not an SPC file")
        return False

    target.write_bytes(data)
    digest = hashlib.sha256(data[:512]).hexdigest()[:16]
    print(f"  {name:<14} {len(data):>7} bytes  header {digest}")
    return True


def main():
    DESTINATION.mkdir(parents=True, exist_ok=True)
    print(f"Fetching {len(SAMPLES)} Galactic sample files into {DESTINATION}")
    print("Source: Internet Archive copies of thermogalactic.com\n")

    ok = sum(fetch(name, stamp, size, DESTINATION)
             for name, (stamp, size) in sorted(SAMPLES.items()))

    print(f"\n{ok}/{len(SAMPLES)} available.")
    if ok != len(SAMPLES):
        print("The archive rate-limits; waiting a minute and re-running "
              "usually finishes the job.")
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
