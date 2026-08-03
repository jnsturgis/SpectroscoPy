# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Created on Wed Feb 19 2025

@author: James STURGIS

This file implements import and export of spectra with CSV format files
for my spectroscopy programmes and utilities.

"""

# pylint: disable=W0718
import numpy as np

from spectroscopy.io.registry import register_reader, register_writer
from spectroscopy.io.table import (
    parse_number,
    sniff_decimal,
    sniff_delimiter,
    sniff_format,
)


def looks_like_a_header(line, delimiter, usecols, decimal='.'):
    """
    True when ``line`` cannot be read as a pair of numbers.

    Used to decide whether the first row is column labels or data. Assuming a
    header unconditionally is what silently ate the first data point of every
    headerless file -- the .csv sibling of defect D1.
    """
    fields = line.split(delimiter)
    try:
        for index in usecols:
            parse_number(fields[index], decimal)
    except (ValueError, IndexError):
        return True
    return False


@register_reader('csv', extensions=['.csv'],
                 description='separator and decimal comma sniffed')
def read(  filehandle, my_spectrum, **kwargs):
    """
    Update the contents of my_spectrum based on the file.

    Keyword arguments
    -----------------
    delimiter : str or None
        Field separator. ``None`` (the default for ``.csv``) sniffs it, which
        matters because a file named ``.csv`` frequently is not comma
        separated -- see below.
    decimal : str or None
        ``'.'`` or ``','``. ``None`` sniffs it.
    skiprows : int or None
        Number of leading rows to treat as a header. ``None`` (the default)
        decides by looking: a first row that does not parse as numbers is a
        header, anything else is data. Pass an integer to force it.

    Notes
    -----
    A spreadsheet exported under a French, German, Spanish or Italian locale
    writes ``0,1234``, and because the comma is then taken by the decimal it
    uses ``;`` for the field separator. So a file called ``.csv`` from a
    European colleague is routinely neither comma separated nor dot decimal,
    and nothing in the file says so. Both are worked out from the numbers,
    together, because neither can be decided without the other.

    Pass either explicitly to override the sniffing. ``.tsv`` and ``.txt``
    come through here with the separator already pinned to a tab, and still
    sniff the decimal -- the locale problem is independent of the separator.
    """
    options = {
        'comments': '#',
        'delimiter': None,
        'decimal': None,
        'skiprows': None,
        'usecols' :(0,1),
    }
    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    usecols = options['usecols']
    skiprows = options['skiprows']

    # Read it all: the separator can only be chosen by looking at several
    # lines at once, and these are small text files.
    lines = [line.rstrip('\r\n') for line in filehandle]
    body = [line for line in lines
            if line.strip() and not line.startswith(options['comments'])]

    delimiter, decimal = options['delimiter'], options['decimal']
    if delimiter is None and decimal is None:
        delimiter, decimal = sniff_format(body[:20])
    elif decimal is None:
        decimal = sniff_decimal(body[:20], delimiter)
    elif delimiter is None:
        delimiter = sniff_delimiter(body[:20])

    if skiprows is None:
        skiprows = (1 if body and looks_like_a_header(body[0], delimiter,
                                                      usecols, decimal)
                    else 0)

    x = []
    y = []
    for index, line in enumerate(body):
        cols = line.split(delimiter)
        if index < skiprows:
            if index == skiprows - 1:      # the row nearest the data names it
                try:
                    my_spectrum.x_label = cols[usecols[0]]
                    my_spectrum.y_label = cols[usecols[1]]
                except IndexError:
                    pass
            continue
        x.append(parse_number(cols[usecols[0]], decimal))
        y.append(parse_number(cols[usecols[1]], decimal))

    my_spectrum.x = np.array(x)
    my_spectrum.y = np.array(y)

@register_writer('csv')
def write( filehandle, my_spectrum, **kwargs):
    """
    Write my_spectrum to a file.
    """
    options = {
        'delimiter': ',',
    }
    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})
    delimiter = options['delimiter']
    filehandle.write(f'{my_spectrum.x_label}{delimiter}{my_spectrum.y_label}\n')
    for x,y in zip(my_spectrum.x, my_spectrum.y):
        filehandle.write(f'{x:.3f}{delimiter}{y:.5f}\n')

## ============================================================================

def main():
    """
    A main routine to do more or less nothing!
    """
    print("This file provides routines for reading and writing csv files")
    return True

if __name__ == '__main__':
    main()


# 'tsv' is the same reader with a tab separator, registered under its own name
# because the notebooks name it explicitly. '.txt' lands here too: it is what
# np.savetxt writes, and those files travel between notebooks.
register_reader('tsv', extensions=['.tsv', '.txt'], delimiter='\t',
                description='tab separated, header row sniffed')(read)
register_writer('tsv', extensions=['.tsv', '.txt'])(
    lambda handle, spectrum, **kwargs: write(handle, spectrum,
                                             **{'delimiter': '\t', **kwargs}))
