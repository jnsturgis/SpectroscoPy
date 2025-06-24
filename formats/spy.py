"""
Created on Wed June 13 2025

@author: James STURGIS

This file implements import and export of spectra with spy format files
for my spectroscopy programmes and utilities.

This file format is designed to allow "complete" saving and restoral of
spectra objects.
"""

# pylint: disable=W0718, W0611
import ast
import numpy as np
import spectroscopy as spc

def read(  filehandle, my_spectrum, **kwargs):
    """
    Update the contents of my_spectrum based on the file.
    """
    options = {
    }

    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    x = []
    y = []
    the_metadata = ""
    version = None

    for line in filehandle:
        if version is None:         # Parse the first line and set version.
            version = 0.0
            header = True
        elif header:                # Handle header info (name and columns)
            if line[0] == '#':
                header = False
                body = True
        elif body:                  # Read x,y values
            if line[0] == '#':
                body = False
                meta = True
            else:
                words = line.strip().split('\t')
                x.append(float(words[0]))
                y.append(float(words[1]))
        elif meta:                  # Read metadata
            the_metadata += ' ' + line.strip()

    my_spectrum.x_data = np.array(x)
    my_spectrum.y_data = np.array(y)
    my_spectrum.metadata = ast.literal_eval(the_metadata)

def write( filehandle, my_spectrum, **kwargs):
    """
    Write my_spectrum to a file.
    """
    options = {
    }
    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    filehandle.write('# spy format 0.0 file\n')
    filehandle.write(f'{my_spectrum.name}\n')
    filehandle.write(f'{my_spectrum.x_label}\t{my_spectrum.y_label}\n')
    filehandle.write('# x,y data\n')
    for x,y in zip(my_spectrum.x_data, my_spectrum.y_data):
        filehandle.write(f'{x:.3f}\t{y:.5f}\n')
    filehandle.write('# metadata\n')
    filehandle.write(f'{repr(my_spectrum.metadata)}\n')

## ============================================================================

def main():
    """
    A main routine to do more or less nothing!
    """
    print("This file provides routines for reading and writing spy files")
    return True

if __name__ == '__main__':
    main()
