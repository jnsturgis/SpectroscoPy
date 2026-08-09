#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Fetch a CD reference set from the PCDDB into a local directory.

**Nothing is redistributed by SpectroscoPy.** The PCDDB's terms grant free
*access* to its data; they do not grant redistribution, and copyright is
asserted only over the site's design and implementation. So this script
downloads to your machine, under the terms you accept from the PCDDB, and the
files stay out of this repository (``data/cd_reference/`` is gitignored).

The same reasoning as ``fetch_spc_fixtures.py``: freely downloadable is not
the same as ours to ship.

    python scripts/fetch_pcddb_basis.py --set sp175

**A condition of use of PCDDB data is citation.** Any publication using it
must cite both the original reference for each spectrum and the PCDDB itself.
The manifest this writes carries the accession of every entry so that is
possible; keep it with the data.

.. warning::

   **Untested against real files.** The PCDDB was unreachable throughout
   2026-08-09 and 2026-08-10 -- DNS resolves, both ports time out -- so the
   download URL and the file layout below are inferred from the site's own
   metadata and have never been run against a live server. Under the working
   agreement (roadmap section 15.2) that makes this provisional: the ``.dpt``
   precedent is a reader written against a specification where the format
   turned out not to be what it looked like.

   It is written to **fail loudly** rather than quietly write rubbish: every
   download is checked for being a plausible CD file before it is kept. When
   the site returns, run it and fix what it reports. Do not assume it worked
   because it did not crash.
"""

from __future__ import annotations

import argparse
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = "https://pcddb.cryst.bbk.ac.uk"

#: The published reference sets and their accession ranges.
#:
#: SP175: Lees, Miles, Wien & Wallace (2006) Bioinformatics 22, 1955.
#: SMP180 adds 30 membrane proteins to SP175's soluble ones:
#: Abdul-Gader, Miles & Wallace (2011) Bioinformatics 27, 1630.
SETS = {
    'sp175': [f"CD{n:010d}" for n in range(1000, 71001, 1000)],
    'mp180': [f"CD{n:011d}" for n in range(99000, 128001, 1000)],
}
SETS['smp180'] = SETS['sp175'] + SETS['mp180']

#: Minimum plausible size of a real entry, in bytes. A server error page,
#: a redirect stub or a "temporarily unavailable" notice is far smaller.
MINIMUM_BYTES = 2000

#: Something every PCDDB data file should contain. Checked case-insensitively.
#: If none of these appears, what came back is not a spectrum.
EXPECTED_MARKERS = ('pcddb', 'wavelength', 'cd')


def looks_like_a_spectrum(text):
    """Cheap sanity check, so an error page is never saved as data."""
    lowered = text.lower()
    if not any(marker in lowered for marker in EXPECTED_MARKERS):
        return False, "none of the expected markers appear"
    if lowered.lstrip().startswith('<!doctype') or '<html' in lowered[:400]:
        return False, "the server returned an HTML page, not a data file"
    numbers = sum(1 for line in text.splitlines()
                  if line[:1].isdigit() or line[:1] == '-')
    if numbers < 50:
        return False, f"only {numbers} lines start with a number"
    return True, ""


def fetch(accession, destination, timeout=30):
    """One entry. Returns the path written, or raises."""
    url = f"{BASE}/deposit/{accession}"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as response:
            raw = response.read()
    except (urllib.error.URLError, TimeoutError, OSError) as error:
        raise RuntimeError(f"{accession}: {url} -- {error}") from error

    if len(raw) < MINIMUM_BYTES:
        raise RuntimeError(
            f"{accession}: {len(raw)} bytes, below the {MINIMUM_BYTES} a real "
            f"entry should have. Probably an error page.")
    text = raw.decode('utf-8', errors='replace')
    ok, why = looks_like_a_spectrum(text)
    if not ok:
        raise RuntimeError(f"{accession}: {why}. First line: "
                           f"{text.splitlines()[0][:80]!r}")

    path = destination / f"{accession}.txt"
    path.write_text(text, encoding='utf-8')
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[1])
    parser.add_argument('--set', dest='which', default='sp175',
                        choices=sorted(SETS),
                        help="which published reference set to fetch")
    parser.add_argument('--into', default='data/cd_reference',
                        help="destination directory (gitignored)")
    parser.add_argument('--limit', type=int, default=None,
                        help="stop after N entries, for a trial run")
    arguments = parser.parse_args(argv)

    accessions = SETS[arguments.which][:arguments.limit]
    destination = Path(arguments.into) / arguments.which
    destination.mkdir(parents=True, exist_ok=True)

    print(f"{len(accessions)} entries -> {destination}")
    print("PCDDB grants access, not redistribution: these stay on your "
          "machine.\nCiting both the original reference and the PCDDB is a "
          "condition of use.\n")

    written, failures = [], []
    for accession in accessions:
        try:
            path = fetch(accession, destination)
        except RuntimeError as error:
            failures.append(str(error))
            print(f"  FAILED  {error}")
        else:
            written.append((accession, path))
            print(f"  ok      {accession}  {path.stat().st_size:>8} bytes")

    if written:
        manifest = destination / 'basis.csv'
        with manifest.open('w', encoding='utf-8') as handle:
            handle.write('file,name,accession,source,category\n')
            for accession, path in written:
                handle.write(
                    f"{path.name},{accession},{accession},"
                    f"PCDDB {arguments.which.upper()},\n")
        print(f"\nwrote {manifest}")
        print("The 'category' column is deliberately empty: only you can say "
              "what each\nspectrum is of. For a reference-protein set, replace "
              "it with helix/sheet/\nturn/other fraction columns. See "
              "library.load_basis.")

    if failures:
        print(f"\n{len(failures)} of {len(accessions)} failed. This script "
              f"has never run against a live PCDDB -- if every entry failed "
              f"the same way, the URL pattern or the layout is wrong rather "
              f"than the server.", file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
