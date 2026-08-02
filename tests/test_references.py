# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""The bibliography and the references page must not drift apart.

Two files describe the same thing -- ``docs/references.bib`` for machines and
the JOSS paper, ``docs/references.md`` for readers -- which is exactly the
duplication that review section 14.4 says should be caught structurally rather
than by review. Without sphinxcontrib-bibtex there is no automatic link between
them, so these tests are the link.
"""

import re
from pathlib import Path

import pytest

DOCS = Path(__file__).resolve().parent.parent / 'docs'
BIB = DOCS / 'references.bib'
PAGE = DOCS / 'references.md'


def _entries():
    """(key, body) for every entry in the .bib."""
    text = BIB.read_text(encoding='utf-8')
    return [(match.group(2), match.group(0))
            for match in re.finditer(r'@(\w+)\{([^,]+),(.*?)\n\}',
                                     text, re.S)]


def test_the_bibliography_parses():
    assert _entries(), "no entries found -- has the file format changed?"


@pytest.mark.parametrize('field', ['author', 'title', 'year'])
def test_every_entry_has_the_minimum_fields(field):
    missing = [key for key, body in _entries()
               if not re.search(rf'\n\s*{field}\s*=', body)]
    assert not missing, f"entries with no {field}: {missing}"


def test_every_entry_is_either_verified_or_marked_unverified():
    """No entry may be silently unchecked.

    An unverified citation that looks like a verified one is worse than no
    citation: it will be copied into a methods section by someone who assumed
    somebody had looked.
    """
    unmarked = [
        key for key, body in _entries()
        # \s* because the fields are column-aligned: 'verified  = {...}'.
        if not re.search(r'\bverified\s*=', body) and 'UNVERIFIED' not in body
    ]
    assert not unmarked, (
        f"entries neither verified nor marked UNVERIFIED: {unmarked}. "
        "Add verified = {YYYY-MM-DD} after checking against the publisher, "
        "or note = {UNVERIFIED -- ...} until someone does."
    )


def test_the_page_marks_exactly_the_unverified_entries():
    """The warning triangle on the page must match the .bib, both ways."""
    page = PAGE.read_text(encoding='utf-8')
    unverified_authors = {
        re.search(r'author\s*=\s*\{([^,}]+)', body).group(1).split()[-1].strip('{}')
        for key, body in _entries() if 'UNVERIFIED' in body
    }
    for surname in unverified_authors:
        pattern = rf'⚠[^\n]*(?:\n[^\n]*){{0,4}}{re.escape(surname)}'
        assert re.search(pattern, page), (
            f"{surname} is UNVERIFIED in references.bib but is not marked ⚠ "
            "in references.md"
        )


def test_verified_entries_are_not_marked_on_the_page():
    page = PAGE.read_text(encoding='utf-8')
    for key, body in _entries():
        if not re.search(r'\bverified\s*=', body):
            continue
        surname = re.search(r'author\s*=\s*\{([^,}]+)',
                            body).group(1).split()[-1].strip('{}')
        for line in page.splitlines():
            if surname in line and line.strip().startswith('- ⚠'):
                pytest.fail(f"{surname} is verified but still marked ⚠ on the page")


def test_the_page_is_linked_from_the_documentation_index():
    index = (DOCS / 'index.md').read_text(encoding='utf-8')
    assert 'references' in index, "references.md is not in the toctree"
