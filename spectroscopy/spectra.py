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

import numpy as np

import formats.jcamp
import spectroscopy.messages

class Spectrum:
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
            self.name     = 'unnamed'
            self.filename = ''
            self.filetype = 'csv'
            self.x_label  = ('Wavelength','nm')
            self.y_label  = ('OD','None')
            self.x_data   = np.empty(1)
            self.y_data   = np.empty(1)
            self.acquisition = {}
            self.treatment   = {}
            self.children    = {}
        elif isinstance(args[0], Spectrum ) :
            # This is the copy method - initiate from another Spectrum()
            other = args[0]
            self.name        = str(other.name)
            self.filename    = ''
            self.filetype    = 'csv'
            self.x_label     = (str(other.x_label[0]),str(other.x_label[1]))
            self.y_label     = (str(other.y_label[0]),str(other.y_label[1]))
            self.acquisition = {} # other.acquisition.copy()
            self.treatment   = {} #other.treatment.copy()
            self.children    = {} #other.children.copy()
            self.x_data      = np.copy(other.x_data)
            self.y_data      = np.copy(other.y_data)
        elif isinstance(args[0], str ) :
            # This is the open method - initialize from a file.
            #
            # TODO implement file initializer
            if len(args) > 1:
                self.filetype = args[1]
            else:
                # TODO detect the file type!!
                self.filetype  = "jcamp"
            self.filename  = args[0]
            match self.filetype :
                case "csv":
                    self.__reread_csv()
                case "jcamp":
                    self.__reread_jcamp()
                case _:
                    raise TypeError(f"Unknown filetype {self.filetype}")
        else :
            raise TypeError("Unknown spectrum initializer")
        # TODO Should finish by checking that we have a well formed Spectrum()
        pass
            
    def y_string(self) -> str:
        # TODO Get this to work
        return self.y_label[0] + (
            "" if self.y_label[1] == 'None' 
            else (" ("+self.y_label[1]+")") )

    def x_string(self) -> str:
        # TODO Get this to work
        return self.x_label[0] + (
            "" if self.x_label[1] == 'None' 
            else (" ("+self.x_label[1]+")") )

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
        return f"Spectrum :{self.name} {self.y_string()} vs {self.x_string()}"

    def __repr__(self) -> str:      # Developper string version of object
        # TODO
        pass
    
    def __len__(self) -> int:       # Length of object
        return self.npts()
    
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
# Some usefull private functions.
#
##=============================================================================


    
##=============================================================================
#
# Simple information about the spectrum.
#
##=============================================================================

    def xmax(self) -> float:
        return np.amax(self.x_data)

    def xmin(self) -> float:
        return np.amin(self.x_data)

    def ymax(self) -> float:
        return np.amax(self.y_data)

    def ymin(self) -> float:
        return np.amin(self.y_data)

    def npts(self) -> int:
        return len(self.x_data)

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
    def resample( self, x_values, *args ):
        """Change the position of the spectrum points
        
        Use the new series of x_values and estimate the y_values at these
        positions using one of several methods provided in the arg list.
        
        This is necessary for doing arithmatic on spectra if they are not
        sampled at the same positions, due to differences in spectrometer
        settings or coming from different machines for example.
        

        Parameters
        ----------
        x_values : numpy array of float.
            New x values at which to estimate the spectrum points.
        *args : array of optional arguments
            Information for the different possible methods to use for the
            estimation.

        Returns
        -------
        Spectrum
            A new spectrum containing the provided x_values and calculated
            y_values.
            name     = original + 'resampled'
            filename = ''
            filetype = 'csv'
            x_label  = original
            y_label  = original
            acquisition = {}
            treatment   = {}
            children    = {}
        """
        result = Spectrum()
        result.name    = self.name + ' baseline'
        result.x_label = self.x_label
        result.y_label = self.y_label
        result.x_data  = x_values
        result.y_data = np.zeros(len(x_values))
        
        # Version 0.0 
        # Linear inter/extra polation using 2 closest points.
        # TODO more advanced version with polynomial fits to series of points
        
        n = len(args)
        if n > 0 :
            if args[0] == "linear":
                control = 0
            elif args[0] == "polynomial":
                if n < 2 :
                    raise ValueError("Spectrum resample 'polynomial' requires a control value")
                control = int(args[1])
            pass
        else:
            # Use default method
            method = "linear"
            control = 0
            pass

        # TODO implement choice of methods using args
        # get closest points from self
        for i in range(len(x_values)):
            x = x_values[i]
            d = abs(self.x_data[0] - x) + 1 # Ensure first test true
            for j in range(len(self.x_data)):
                if d > abs(self.x_data[j] - x):
                    d = abs(self.x_data[j] - x)
                    index0 = j
                else :
                    break
            
            # index0 is the index of the closest x value in self.x_data.
            
            if index0 == 0:
                index1 = 1
            elif index0 == len(self.x_data) - 1:
                index1 = index0
                index0 = index1 - 1
            elif abs(self.x_data[index0+1] - x) > abs(self.x_data[index0-1] - x):
                index1 = index0
                index0 = index0 - 1
            else:
                index1 = index0 + 1

            assert index0 < index1
            assert index0 >= 0
            assert index1 < len(self.x_data)
            
            r = (x - self.x_data[index0])/(self.x_data[index1]-self.x_data[index1])
            delta = self.y_data[index0] - self.y_data[index1]
            result.y_data[i] = self.y_data[index0] + r * delta
            
        # TODO Check we have a well formed spectrum
        return result
    
    def baseline(self, method, *args ):
        """Create a baseline spectrum
        
        Create a baseline spectrum using a method and a set of control
        points possibly using the spectrum to find them.

        Parameters
        ----------
        method : a tuple of string and parameters
            Describes how to construct the baseline methods include
            1. "polynomial", order - use a polynomial of order n to fit the
            points given the args array(s).
            2. "spline", order - use a spline function to interpolate the 
            control points given the args array(s).
        *args : a list of one or two lists of values.
            If there is one list it is interpreted as the set of x_values
            for the control points, the y values are extracted from the
            current spectrum.
            If there are two lists, they are interpreted as a set of x_values
            and a set of y_values for the control points.

        Returns
        -------
        Spectrum
            returns a baseline spectrum constructed using the defined method
            and control points at x_values defined by the current spectrum.
            name     = original + 'baseline'
            filename = ''
            filetype = 'csv'
            x_label  = original
            y_label  = original
            acquisition = {}
            treatment   = {}
            children    = {}


        """
        result = Spectrum()
        result.name    = self.name + ' baseline'
        result.x_label = self.x_label
        result.y_label = self.y_label
        result.x_data  = self.x_data
        
        if len(args) == 1:
            x_values = args[0]
            y_values = np.zeros(len(x_values))
            
            # Version 0.0
            # get closest points from self
            for i in range(len(x_values)):
                x = x_values[i]
                d = abs(self.x_data[0] - x) + 1 # Ensure first test true
                for j in range(len(self.x_data)):
                    if d > abs(self.x_data[j] - x):
                        d = abs(self.x_data[j] - x)
                        xv = self.x_data[j]
                        yv = self.y_data[j]
                    else :
                        break
                x_values[i] = xv
                y_values[i] = yv
        elif len(args) == 2:
            x_values = args[0]
            y_values = args[1]
        else :
            # Should never arrive here throw an error
            raise ValueError("Spectrum baseline requires control value array")
        
        # Have the control points
        if method[0] == "polynomial" :
            order = method[1]
            assert order < len(x_values)
            p_coef = np.polyfit(x_values, y_values, order)
            result.y_data = np.polyval(p_coef, result.x_data)
            pass
        elif method[0] == "spline" :
            order = method[1]
            raise NotImplementedError("Spline interpolation is not yet implemented")
            pass
        else:
            # Should never arrive here
            raise ValueError("Spectrum baseline requires a method tuple")
            pass
        
        # TODO Check we have a well formed spectrum
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

    def __reread_jcamp(self) -> None:

        with open(self.filename, 'rb') as filehandle:
            datadict = formats.jcamp.jcamp_read(filehandle)    
        self.x_data   = datadict['x']
        self.y_data   = datadict['y']
        self.name     = datadict['title']
        self.x_label  = ( "Wavenumber", datadict['xunits'].lower() )
        self.y_label  = ( datadict['yunits'].capitalize(), "" )
    
def write_jcamp_file(filename: str, my_spectrum: Spectrum ) -> None:

    linewidth=75

    # Pack the data into a jcamp dictionary
    
    jcamp_dict = {}
    jcamp_dict['filename'] = my_spectrum.filename
    jcamp_dict['x']        = my_spectrum.x_data
    jcamp_dict['y']        = my_spectrum.y_data

    with open(filename, 'w') as filehandle:
        jcamp_str = formats.jcamp.jcamp_write(jcamp_dict, linewidth=linewidth)
        filehandle.write(jcamp_str)
    
    return

##=============================================================================

    def __reread_csv(self):
        data = np.genfromtxt(self.filename, delimiter=',')
        self.x_data   = data[:,0]
        self.y_data   = data[:,1]

def write_csv_file(filename: str, my_spectrum: Spectrum) -> None:
    np.savetxt( filename, 
               np.column_stack((my_spectrum.x_data, my_spectrum.y_data)),
               fmt = '%.6f',
               delimiter = ',')
    
    return

##=============================================================================

