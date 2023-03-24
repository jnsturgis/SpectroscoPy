""" Spectra structure to hold a spectrum and associated metadata

The spectrum structure
The metadata associated with the spectrum and default values:

name	 a name for the spectrum initially derived from the filename.
filename the name of the file the spectrum came from.
x_label  the x axis label (defaults to Wavelength)
x_units  the x axis units (defaults to nm)
y_label  the y axis label (defaults to OD)
y_units  the y axis units (defaults to NONE)

data     a  numpy array of 2 columns with x, y pairs, with the x values in 
         either ascending or descending order

Acquisition metadata - for original untreated data it is interesting to have
information on the acquisition parameters.

Treatment metadata - is it interesting to have the history? Certainly not 
for the moment.

Then there are derived values that can either be stored or calculated on 
as they are needed.

x_min    minimum value of x data
x_max    maximum value of y data
y_min    minimum value of y data
y_max    maximum value of y data
n_points the number of points in the spectrum
even_p   are the x values evenly spaced - T/F

Functions to load and save spectra - these functions should be able to
handle multiple pertinent formats (JCAMP, CSV w/wo header, etc)

export   write a spectrum in a given format
import   read a spectrum from a given format
import2  read a spectrum guessing the format
save     save a spectrum
save_as  save a spectrum with a new name

Then there are functions and operators that produce a new spectrum with
modified data. TODO think about what to do to any acquisition metadata
in these cases. As these modify spectra they are esentially given as functions
that return a new spectrum derived from an original spectrum or more.

x_shift  shift the x scale by a fixed amount
x_mod    modify the x scale (for example change from Wavelength to Frequency)
x_even   resample to get an even x_axis with the requested number of points
y_shift  shift the specrtum on the y_axis
y_mod    modify the y scale (for example Absorption to Transmission)

+-/*     simple maths operators with real constants or compatible spectra



"""

import numpy as np

class Spectrum:
    def __init__(self):
        self.name = ''
        self.filename = ''
        self.x_label = 'Wavelength'
        self.x_units = 'nm'
        self.y_label = 'OD'
        self.y_units = 'NONE'
        self.data = np.empty(shape=[0,2])

    def x_min(self):
        return np.min(self.data[:,0])

    def x_max(self):
        return np.max(self.data[:,0])

    def y_min(self):
        return np.min(self.data[:,1])

    def y_max(self):
        return np.max(self.data[:,1])

    def n_points(self):
        return np.shape(self.data)[0]

    def even_p(self):
        n = self.n_points
        if n <= 2 return True                    # Always true for small
        x_data = self.data[:,0]
        first = x_data[0]
        gap = x_data[1] - first
        if gap == 0:                             # Constant x values
            return all(x == first for x in lst)
        r = range( first, first + gap * n, gap)
        return all( x == y for x,y in zip( x_data, r ))
        
    def y_string(self):
        return self.y_label + 
               ( "" if self.y_units == 'NONE' else (" ("+self.y_units+")") )
        
    def x_string(self):
        return self.x_label + 
               ( "" if self.x_units == 'NONE' else (" ("+self.x_units+")") )
        
    def __str__(self):
        return f"Spectrum :{self.name} {self.y_string} vs {self.x_string}"

    def read_csv(self, filename: str):
        self.data = np.genfromtxt(filename, delimiter=',')
        self.filename = filename

    def write_csv(self, filename: str):
        numpy.savetxt(filename, self.data, delimiter=",")
    
    def export(self, filename: str, file_fmt: str):
        return True
    
def smooth(a_spectrum, parameters: list):
    return a_spectrum
