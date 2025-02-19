"""
Created on Wed Feb 19 2025

@author: James STURGIS

This file implements import and export of spectra with CSV format files
for my spectroscopy programmes and utilities.

"""

# pylint: disable=W0718, W0611
import numpy as np
import spectroscopy as spc

def read(  filehandle, my_spectrum, **kwargs):
    """
    Update the contents of my_spectrum based on the file.
    """
    options = {
        'comments': '#',
        'delimiter': ',',
        'skiprows': 1,
        'usecols' :(0,1),
    }
    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    skiprows = options['skiprows']
    x = []
    y = []
    for line in filehandle:
        if line[-1] == '\n':
            line = line[:-1]
        if len(line) > 0:
            if line[0]==options['comments']:
                continue
            if skiprows == 1:
                try:
                    cols = line.split(options['delimiter'])
                    my_spectrum.x_label = cols[options['usecols'][0]]
                    my_spectrum.y_label = cols[options['usecols'][1]]
                except Exception:
                    pass
            if skiprows <= 0:
                cols = line.split(options['delimiter'])
                x.append(float(cols[options['usecols'][0]]))
                y.append(float(cols[options['usecols'][1]]))
        skiprows -= 1
    my_spectrum.x_data = np.array(x)
    my_spectrum.y_data = np.array(y)

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
    for x,y in zip(my_spectrum.x_data, my_spectrum.y_data):
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
