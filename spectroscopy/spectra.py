""" Spectra structure to hold a spectrum and associated metadata

The spectrum structure
The metadata associated with the spectrum and default values:

name	 a name for the spectrum initially derived from the filename.
filename the name of the file the spectrum came from.
x_labels a list of 2 tuples giving (label, units) for each x dimension.
y_label  a 2 tuple describing the y data (name, units)

x_data   a numpy array of x values ordered along all axes
y_data   a numpy array of y values at positions given by the x values

Acquisition metadata - for original untreated data it is interesting to have
information on the acquisition parameters.

Treatment metadata - is it interesting to have the history? Certainly not
for the moment.

TODO:
    Should define functions for jcamp required info and use it.
    Should make sanity checks and die gracefully

"""

# pylint: disable=W0511, W0107
import os

import numpy as np

from formats import jcamp, csv

import spectroscopy.messages

FILE_EXTS = {
    '.csv':'csv','.tsv':'tsv',
    '.jcamp':'jcamp','.DX0':'jcamp'
}

KNOWNFILETYPES = ('csv','jcamp')

def _infer_file_type( name ):
    """
    Work out if possible from the file extension the possible file types.
    """
    name, ext = os.path.splitext(name)
    if ext in FILE_EXTS:
        return FILE_EXTS[ext]
    return 'unknown'

class Spectrum:
    """
    A class for spectra
    """
    def __init__(self, *args):
        """ Initialize a Spectrum object

        Simulate overloading the different possitilities for initializing
        a spectrum:

        Parameters
        ----------
        args a set of parameters that allow for different initiation routines
        depending on the parameters.


        Spectrum()                     : no arguments an empty spectrum.
        Spectrum( filename )           : a file
        Spectrum( filename, filetype ) : a file where filetype can be from
                 ('jcamp','csv')
        Spectrum( a_spectrum )         : a deepcopy of the original

        """

        if len(args) == 0 :
            #     This is the empty initializer
            self.name        = 'unnamed'
            self.fileinfo    = {'PATH':'','NAME':'','TYPE':'csv'}
            self.x_label     = 'Wavelength (nm)'
            self.y_label     = 'Absorbance'
            self.x_data      = np.empty(1)
            self.y_data      = np.empty(1)
            self.metadata    = {}

        elif isinstance(args[0], Spectrum ) :
            # This is the copy method - initiate from another Spectrum()
            other = args[0]
            self.name        = str(other.name)
            self.fileinfo    = {'PATH':'','NAME':'','TYPE':'csv'}
            self.x_label     = other.x_label
            self.y_label     = other.y_label
            self.x_data      = np.copy(other.x_data)
            self.y_data      = np.copy(other.y_data)
            self.metadata    = other.metadata

        elif isinstance(args[0], str ) :
            # This is the open method - initialize from a file
            # First parse the args to get filetype and full filename
            # Len 1: filename - infer type from 1
            # Len 2: path, filename - infer type  from 2
            # Len 3: path, filename, type
            if len(args) == 3:
                self.fileinfo['PATH'] = args[0]
                self.fileinfo['NAME'] = args[1]
                self.fileinfo['TYPE'] = args[2]
            elif len(args) == 2:
                self.fileinfo['PATH'] = args[0]
                self.fileinfo['NAME'] = args[1]
                self.fileinfo['TYPE'] = _infer_file_type(args[0]+args[1])
            elif len(args) == 1:
                self.fileinfo['NAME'] = args[0]
                self.fileinfo['TYPE'] = _infer_file_type(args[0])

            # Now read the file according to type
            if self.fileinfo['TYPE'] in KNOWNFILETYPES:
                self.reload()
            else:
                raise TypeError(f"Unknown filetype {self.fileinfo['TYPE']}")
        else :
            raise TypeError("Unknown spectrum initializer")
        # TODO Should finish by checking that we have a well formed Spectrum()
        pass

##=============================================================================
#
#   Some Magic Dunder methods
#
#   Not comparison operations (does not make sense)
#   Not floordiv, mod or power operators
#
#   TODO add an iterator method for point in spectrum.
#
##=============================================================================

    def __str__(self) -> str:       # Users string version of object
        return f"Spectrum :{self.name} {self.y_label()} vs {self.x_label()}"

    def __repr__(self) -> str:      # Developper string version of object
        # TODO
        pass

    def __len__(self) -> int:       # Length of object
        return len(self.x_data)

    def __add__(self, other):
        # Check that self and other are compatible x and y ranges and scales
        # And copy almost all of self to new
        new_spectrum = Spectrum(self)
        if isinstance(other, Spectrum):
            # Resample data if necessary
            new_spectrum.y_data = self.y_data + other.y_data
        elif isinstance(other, int | float):
            new_spectrum.y_data = self.y_data + other
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        return new_spectrum

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other ):
        # Check that self and other are compatible x and y ranges and scales
        # And copy almost all of self to new
        new_spectrum = Spectrum(self)
        if isinstance(other, Spectrum):
            new_spectrum.y_data = self.y_data - other.y_data
        elif isinstance(other, int | float):
            new_spectrum.y_data = self.y_data - other
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        return new_spectrum

    def __rsub__(self, other ):
        new_spectrum = Spectrum(self)
        if isinstance(other, int | float):
            new_spectrum.y_data = other - self.y_data
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        return new_spectrum

    def __neg__(self):
        return 0.0 - self

    def __pos__(self):
        return self + 0.0

    def __mul__(self, other):
        # Check that self and other are compatible x and y ranges and scales
        # And copy almost all of self to new
        new_spectrum = Spectrum(self)
        if isinstance(other, Spectrum):
            new_spectrum.y_data = self.y_data * other.y_data
        elif isinstance(other, int | float):
            new_spectrum.y_data = self.y_data * other
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        return new_spectrum

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        # Check that self and other are compatible x and y ranges and scales
        # And copy almost all of self to new
        new_spectrum = Spectrum(self)
        if isinstance(other, Spectrum):
            new_spectrum.y_data = self.y_data / other.y_data
        elif isinstance(other, int | float):
            new_spectrum.y_data = self.y_data / other
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        return new_spectrum

    def __rtruediv__(self, other):
        new_spectrum = Spectrum(self)
        if isinstance(other, int | float):
            new_spectrum.y_data = other / self.y_data
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        return new_spectrum


##=============================================================================
#
#   Some functions to implement for current applications
#
##=============================================================================

    def plot(self, ax, **kwargs) -> None:
        """
        A routine to plot a spectrum on some matplotlib axes.
        """
        ax.plot(self.x_data, self.y_data, label=self.name, **kwargs)

##=============================================================================
#
#   Some functions to implement for current applications
#
##=============================================================================

#   smooth( method_name, parameteres)
#   derivative( method_name, parameters )
#   clip( x_range )
#   convert( new_label )
#   peaks( method, parameters )
    def resample( self, x_values, **kwargs ):
        """
        Change the position of the spectrum points

        Use the new series of x_values and estimate the y_values at these
        positions using one of several methods provided in the arg list.

        This is necessary for doing arithmatic on spectra if they are not
        sampled at the same positions, due to differences in spectrometer
        settings or coming from different machines for example.


        Parameters
        ----------
        x_values : numpy array of float.
            New x values at which to estimate the spectrum points.
        **kwargs : dictionary of optional arguments
            Information for the different possible methods to use for the
            estimation.

        Returns
        -------
        Spectrum
            A new spectrum containing the provided x_values and calculated
            y_values.
        """
        result = Spectrum()
        result.name    = self.name + ' baseline'
        result.x_label = self.x_label
        result.y_label = self.y_label
        result.x_data  = x_values
        result.y_data = np.zeros(len(x_values))

        return result

    def baseline(self, **kwargs ):
        """
        Create a baseline spectrum

        Create a baseline spectrum using a method and a set of control
        points possibly using the spectrum to find them.

        Parameters
        ----------
        **kwargs : a dictionary of optional parameters.

        Returns
        -------
        Spectrum
            returns a baseline spectrum constructed using the defined method
            and control points at x_values defined by the current spectrum.
        """
        result = Spectrum()
        result.name    = self.name + ' baseline'
        result.x_label = self.x_label
        result.y_label = self.y_label
        result.x_data  = self.x_data

        return result

##=============================================================================
#
#   Some gui functions to implement
#
##=============================================================================

#   save_as( new_filename: str)
#   save()

#   edit_metadata()
#   fit_components()
#   adjust_components()
#   display()

##=============================================================================
#
#   Here are the routines to interface to different file formats and returned
#   data structures.
#
#   TODO: Populate the data structure properly from the dictionary
#
##=============================================================================
    def reload(self) -> None:
        """
        Reload the spectrum from the file, or load a first time after setting fileinfo
        """
        filename = self.fileinfo['PATH']+self.fileinfo['NAME']
        with open( filename, 'r', encoding="utf-8") as f:
            match self.fileinfo['TYPE']:
                case 'jcamp':
                    jcamp.read(f, self)
                case 'csv':
                    csv.read(f, self )
                case 'tsv':
                    csv.read(f, self, delimiter='\t')

    def save(self) -> None:
        """
        Write the spectrum to the file described by file info.
        """
        filename = self.fileinfo['PATH']+self.fileinfo['NAME']
        with open(filename, 'w', encoding="utf-8") as f:
            match self.fileinfo['TYPE']:
                case 'jcamp':
                    jcamp.write(f, self)
                case 'csv':
                    csv.write(f, self )
                case 'tsv':
                    csv.write(f, self, delimiter='\t')
