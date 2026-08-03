# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Generic delimited-text reader: many spectra from one wide file.

Roadmap section 3 calls a configurable delimited reader "high value early,
since it unblocks 'weird format my colleague sent me' cases without a
dedicated parser". Three real files in the notebooks need exactly this, and
each defeats a plain two-column reader in a different way:

===============================  =========================================
``J_Peri.csv`` (Chloe)           173 columns as (x, y) **pairs** -- an
                                 excitation-emission series, one spectrum
                                 per excitation wavelength, with the series
                                 names on a sparse first header row.
``sfGFP.csv`` (GFP binding)      one **shared** wavelength column then many
                                 named sample columns.
``Sepharose ... .csv`` (AKTA)    **UTF-16-LE** with a BOM, two header rows,
                                 tab separated, paired columns.
===============================  =========================================

So the reader supports paired and shared-x layouts, sniffs the separator and
the number of header rows, and takes column names from the header. Encoding is
handled upstream in :mod:`spectroscopy.io.registry`.
"""

from __future__ import annotations

import numpy as np

from spectroscopy.io.registry import register_reader

__all__ = ['read_table', 'sniff_delimiter', 'sniff_header_rows',
           'sniff_format', 'parse_number']

CANDIDATE_DELIMITERS = ('\t', ',', ';', None)

#: Decimal separators worth trying. A spreadsheet exported under a French,
#: German, Spanish or Italian locale writes ``0,1234``, and because the comma
#: is then taken by the decimal it exports ``;`` as the field separator -- so a
#: file called ``.csv`` is very often neither comma separated nor dot decimal.
#: Nothing in the file says which; it has to be worked out from the numbers.
CANDIDATE_DECIMALS = ('.', ',')

#: Digit-group separators to discard before parsing: ASCII space, non-breaking
#: space and narrow non-breaking space are all used as thousands separators,
#: and Excel emits the non-breaking ones.
_GROUPING = (' ', ' ', ' ')


def parse_number(text, decimal='.'):
    """
    Parse one field, honouring a comma decimal separator.

    Handles the awkward combination too: ``1.234,56`` is one thousand two
    hundred and thirty-four, and the dots have to go before the comma becomes
    a point. The test is positional -- a dot appearing *before* the last comma
    is grouping, not a decimal.
    """
    text = text.strip()
    for space in _GROUPING:
        text = text.replace(space, '')
    if decimal == ',':
        if '.' in text and text.rfind('.') < text.rfind(','):
            text = text.replace('.', '')
        text = text.replace(',', '.')
    return float(text)


def _is_number(text, decimal='.'):
    try:
        parse_number(text, decimal)
        return True
    except (TypeError, ValueError):
        return False


def _split(line, delimiter):
    return line.split(delimiter) if delimiter is not None else line.split()


def sniff_format(lines):
    """
    Work out ``(delimiter, decimal)`` together.

    They cannot be decided separately. Given ``400,5;0,1234``, splitting on
    ``,`` finds two things that look like numbers -- ``400`` and ``1234`` --
    and splitting on ``;`` finds none at all unless you already believe the
    comma is a decimal point. Scoring every combination and taking the best is
    the only way round it, and it is cheap: four separators by two decimals.

    Ties go to the dot, and then to the earlier separator in
    :data:`CANDIDATE_DELIMITERS`, so an unambiguous file reads exactly as it
    always did.
    """
    best, best_score = (',', '.'), -1e9
    for delimiter in CANDIDATE_DELIMITERS:
        for decimal in CANDIDATE_DECIMALS:
            if delimiter == decimal:
                continue                  # a comma cannot be both at once
            score = 0
            for line in lines:
                fields = [f.strip() for f in _split(line.strip(), delimiter)]
                if len(fields) < 2:
                    continue
                numeric = sum(1 for f in fields if _is_number(f, decimal))
                if not numeric:
                    continue              # a header row; it judges nothing
                # Leftovers count against. Splitting "400,5;0,1234" on the
                # comma yields 400 and 1234, which look like two fine numbers
                # until you notice the "5;0" stranded between them -- counting
                # hits alone scores that as highly as the right answer.
                score += numeric - (len(fields) - numeric)
            if score > best_score:
                best, best_score = (delimiter, decimal), score
    return best


def sniff_delimiter(lines):
    """Pick the separator that yields the most numeric fields."""
    return sniff_format(lines)[0]


def sniff_decimal(lines, delimiter=None):
    """
    Pick the decimal separator, optionally for a delimiter already decided.

    Worth calling directly when the separator is known -- a ``.tsv`` say --
    but the locale of whoever exported it is not.
    """
    if delimiter is None:
        return sniff_format(lines)[1]
    best, best_score = '.', -1
    for decimal in CANDIDATE_DECIMALS:
        if delimiter == decimal:
            continue
        score = sum(1 for line in lines
                    for f in _split(line.strip(), delimiter)
                    if _is_number(f.strip(), decimal))
        if score > best_score:
            best, best_score = decimal, score
    return best


def sniff_header_rows(lines, delimiter, decimal='.'):
    """Count leading rows that are not data (i.e. hold no parsable numbers)."""
    count = 0
    for line in lines:
        if not line.strip():
            count += 1
            continue
        fields = [f.strip() for f in _split(line.strip(), delimiter)]
        if any(_is_number(f, decimal) for f in fields):
            break
        count += 1
    return count


def _carry_forward(header, column, fallback=''):
    """
    Value of ``header`` at ``column``, carrying the last non-empty one forward.

    The naming row of a paired export is sparse: the series name sits over the
    x column of each pair and the partner cell is blank.
    """
    if not header:
        return fallback
    value = ''
    for index in range(min(column + 1, len(header))):
        if header[index].strip():
            value = header[index].strip()
    return value or fallback


def choose_name_row(headers, columns):
    """
    Pick which header row names the series, when there is more than one.

    Two real files disagree about where the name lives, so this decides by
    looking rather than by position:

    * Chloe's fluorimeter export -- row 0 is the sample names, row 1 repeats
      "Wavelength (nm)" / "Intensity (a.u.)" for every pair.
    * The AKTA export -- row 0 repeats the run name "Chrom.1", row 1 holds the
      trace names UV / Conductivity / Injection.

    The naming row is the one whose values actually *distinguish* the series,
    i.e. the one with the most distinct entries at the series positions.
    """
    if not headers:
        return None
    best, best_score = 0, -1
    for index, header in enumerate(headers):
        values = {_carry_forward(header, column) for column in columns}
        values.discard('')
        if len(values) > best_score:
            best, best_score = index, len(values)
    return best


@register_reader('table', extensions=(), multi=True,
                 description='generic delimited text, wide or paired columns')
def read_table(handle, *, x_col=0, y_cols=None, paired=False, delimiter=None,
               decimal=None, header_rows=None, names=None, comments='#',
               max_rows=None):
    """
    Read a delimited text file into a list of Spectrum.

    Parameters
    ----------
    paired : bool
        Columns are ``(x0, y0, x1, y1, ...)``, each pair one spectrum. This is
        how spectrofluorimeters and chromatography systems export a series.
    x_col : int
        Shared x column when ``paired`` is False.
    y_cols : sequence of int, optional
        Which columns are data. Defaults to every column except ``x_col``.
    header_rows : int, optional
        Leading non-data rows. Sniffed when omitted.
    names : sequence of str, optional
        Names for the resulting spectra; taken from the header otherwise.
    """
    from spectroscopy.spectra import Spectrum  # pylint: disable=C0415

    raw = [line.rstrip('\r\n') for line in handle]
    body = [line for line in raw
            if line.strip() and not line.lstrip().startswith(comments)]
    if not body:
        raise ValueError("no data found in table")

    # Separator and decimal are sniffed together -- see sniff_format. Either
    # can be pinned; whichever is left is worked out around it.
    if delimiter is None and decimal is None:
        delimiter, decimal = sniff_format(body[:20])
    elif decimal is None:
        decimal = sniff_decimal(body[:20], delimiter)
    elif delimiter is None:
        delimiter = sniff_delimiter(body[:20])
    if header_rows is None:
        header_rows = sniff_header_rows(body[:10], delimiter, decimal)

    # ﻿ is a byte-order mark left at the start of a decoded UTF-16 file;
    # without stripping it the first series comes out named '﻿Chrom.1'.
    headers = [[f.strip().lstrip('﻿') for f in _split(line, delimiter)]
               for line in body[:header_rows]]
    data_lines = body[header_rows:]
    if max_rows is not None:
        data_lines = data_lines[:max_rows]
    if not data_lines:
        raise ValueError("table has header rows but no data rows")

    width = max(len(_split(line, delimiter)) for line in data_lines)
    table = np.full((len(data_lines), width), np.nan)
    for row, line in enumerate(data_lines):
        for column, text in enumerate(_split(line, delimiter)):
            text = text.strip()
            if text and _is_number(text, decimal):
                table[row, column] = parse_number(text, decimal)

    if paired:
        pairs = [(index, index + 1) for index in range(0, width - 1, 2)]
    else:
        if y_cols is None:
            y_cols = [index for index in range(width) if index != x_col]
        pairs = [(x_col, index) for index in y_cols]

    # In a paired layout the series name sits over the pair's x column; with a
    # shared x column it is the y column that identifies the series.
    name_columns = [x for x, _ in pairs] if paired else [y for _, y in pairs]
    name_row_index = choose_name_row(headers, name_columns)
    name_row = headers[name_row_index] if name_row_index is not None else []

    spectra = []
    for number, (x_index, y_index) in enumerate(pairs):
        xs, ys = table[:, x_index], table[:, y_index]
        keep = ~(np.isnan(xs) | np.isnan(ys))
        if not keep.any():
            continue

        spectrum = Spectrum()
        spectrum.x, spectrum.y = xs[keep], ys[keep]

        if names is not None and number < len(names):
            spectrum.name = names[number]
        else:
            spectrum.name = _carry_forward(
                name_row, name_columns[number], f"column {y_index}")

        # Axis labels come from a *different* header row, and only when its x
        # and y entries differ -- that is what tells a real label row
        # ("Wavelength (nm)" / "Intensity (a.u.)") from a run name repeated
        # across every column ("Chrom.1"), which would be a misleading label.
        for index, header in enumerate(headers):
            if index == name_row_index or len(header) <= max(x_index, y_index):
                continue
            if header[x_index] and header[y_index] and \
                    header[x_index] != header[y_index]:
                spectrum.x_label = header[x_index]
                spectrum.y_label = header[y_index]

        spectrum.metadata['column'] = y_index
        spectrum.set_sample(spectrum.name)
        spectra.append(spectrum)

    if not spectra:
        raise ValueError("no usable (x, y) column pairs found in table")
    return spectra
