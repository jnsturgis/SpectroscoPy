# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Reader/writer for ``.dpt`` data point table files.

``.dpt`` is the plain two-column export from Bruker OPUS (and the format most
of the ATR-FTIR data in this project arrives in). There is no header: every
line is ``x<sep>y``.

Why this is not just ``csv`` with different options
---------------------------------------------------
A survey of the 825 ``.dpt`` files in ~/Documents/Research found the format is
not as uniform as it looks:

    tab separated     681
    comma separated   140      <- OPUS honours the machine's locale
    other               4      <- 2 with a '#' metadata header, 1 binary

So the separator is sniffed per file rather than assumed. Reading these as
``'tsv'`` -- which is what every notebook did -- silently dropped the first
data point of every file (the csv reader's ``skiprows=1`` default consumed it
as a header) and failed outright on the comma-separated ones. That is defect D1
in SpectroscoPy_Codebase_Review.md.

Files are also almost universally CRLF, being written on Windows.

Comment lines beginning with '#' are kept verbatim in
``metadata['file_header']`` rather than being interpreted; some of James' own
reference spectra carry a Sample/Reference/Operator/Date/Machine block there.
"""

import numpy as np

#: Separators tried, in order, when sniffing. ``None`` means "any whitespace".
CANDIDATE_DELIMITERS = ('\t', ',', ';', None)


def split_fields(line, delimiter):
    """Split ``line``; ``delimiter=None`` means any run of whitespace."""
    return line.split(delimiter) if delimiter is not None else line.split()


def sniff_delimiter(line, usecols=(0, 1)):
    """
    Work out which separator a data line uses.

    Returns
    -------
    (found, delimiter) : tuple[bool, str | None]
        ``found`` is False when no candidate yields parsable numbers -- callers
        should treat that as "this is not a data line". ``delimiter`` is None
        when the fields turned out to be whitespace separated, which is why
        this returns a pair rather than using None to signal failure.
    """
    needed = max(usecols) + 1
    for delimiter in CANDIDATE_DELIMITERS:
        fields = split_fields(line, delimiter)
        if len(fields) < needed:
            continue
        try:
            for index in usecols:
                float(fields[index])
        except ValueError:
            continue
        return True, delimiter
    return False, None


def read(filehandle, my_spectrum, **kwargs):
    """
    Read a .dpt file and populate ``my_spectrum``.

    Keyword arguments
    -----------------
    comments : str
        Prefix marking comment lines, default ``'#'``.
    delimiter : str or None
        Force a separator instead of sniffing. ``None`` (the default) sniffs.
    usecols : tuple[int, int]
        Which columns hold x and y, default ``(0, 1)``.
    """
    options = {
        'comments': '#',
        'delimiter': None,
        'usecols': (0, 1),
    }
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    usecols = options['usecols']
    delimiter = options['delimiter']
    # 'delimiter' arriving as None means "sniff it"; once sniffed, None is a
    # legitimate value meaning whitespace, hence the separate flag.
    sniffed = 'delimiter' in kwargs and kwargs['delimiter'] is not None

    xs, ys = [], []
    header = []

    for number, raw in enumerate(filehandle, start=1):
        line = raw.rstrip('\r\n').strip()

        if not line:
            continue
        if line.startswith(options['comments']):
            header.append(line)
            continue

        if not sniffed:
            found, delimiter = sniff_delimiter(line, usecols)
            if not found:
                raise ValueError(
                    f"{getattr(filehandle, 'name', '<dpt>')}: line {number} is not a "
                    f"pair of numbers and no separator could be determined: {line[:40]!r}"
                )
            sniffed = True

        fields = split_fields(line, delimiter)
        try:
            xs.append(float(fields[usecols[0]]))
            ys.append(float(fields[usecols[1]]))
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"{getattr(filehandle, 'name', '<dpt>')}: cannot parse line {number}: "
                f"{line[:40]!r}"
            ) from exc

    if not xs:
        raise ValueError(
            f"{getattr(filehandle, 'name', '<dpt>')}: no data points found"
        )

    my_spectrum.x = np.array(xs)
    my_spectrum.y = np.array(ys)
    if header:
        my_spectrum.metadata['file_header'] = header


def write(filehandle, my_spectrum, **kwargs):
    """
    Write ``my_spectrum`` as a .dpt file: two columns, no header.

    Keyword arguments
    -----------------
    delimiter : str
        Separator to write, default tab (what OPUS emits most often).
    """
    options = {'delimiter': '\t'}
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    delimiter = options['delimiter']
    for x, y in zip(my_spectrum.x, my_spectrum.y):
        filehandle.write(f'{x:.5f}{delimiter}{y:.10f}\n')
