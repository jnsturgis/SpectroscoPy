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


def looks_like_a_header(line, delimiter, usecols):
    """
    True when ``line`` cannot be read as a pair of numbers.

    Used to decide whether the first row is column labels or data. Assuming a
    header unconditionally is what silently ate the first data point of every
    headerless file -- the .csv sibling of defect D1.
    """
    fields = line.split(delimiter)
    try:
        for index in usecols:
            float(fields[index])
    except (ValueError, IndexError):
        return True
    return False


@register_reader('csv', extensions=['.csv'], description='comma separated, header row sniffed')
def read(  filehandle, my_spectrum, **kwargs):
    """
    Update the contents of my_spectrum based on the file.

    Keyword arguments
    -----------------
    skiprows : int or None
        Number of leading rows to treat as a header. ``None`` (the default)
        decides by looking: a first row that does not parse as numbers is a
        header, anything else is data. Pass an integer to force it.
    """
    options = {
        'comments': '#',
        'delimiter': ',',
        'skiprows': None,
        'usecols' :(0,1),
    }
    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    delimiter = options['delimiter']
    usecols = options['usecols']
    skiprows = options['skiprows']
    sniffing = skiprows is None

    x = []
    y = []
    for line in filehandle:
        line = line.rstrip('\r\n')
        if not line:
            continue
        if line[0] == options['comments']:
            continue

        if sniffing:
            # Decide once, on the first non-blank, non-comment row.
            skiprows = 1 if looks_like_a_header(line, delimiter, usecols) else 0
            sniffing = False

        if skiprows == 1:
            cols = line.split(delimiter)
            try:
                my_spectrum.x_label = cols[usecols[0]]
                my_spectrum.y_label = cols[usecols[1]]
            except IndexError:
                pass
        elif skiprows <= 0:
            cols = line.split(delimiter)
            x.append(float(cols[usecols[0]]))
            y.append(float(cols[usecols[1]]))
        skiprows -= 1

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
