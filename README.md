# SpectroscoPy

## Introduction

This project is born out of frustration as a scientist that uses a variety of
spectroscopy instruments but frequently runs into problems with idiosyncratic
and low quality anaysis software, proprietary file formats. The objectives of
this project are therefor:
* Provide a simple to use spectroscopy analysis package that is able to handle
and perform standard analyses for a variety of different spectroscopies.
* Provide an interface for producing high quality scientific figures using
matplotlib.
* Provide handles to allow data acquisition interfacing in cases where
instrument manufacturers make this possible.
* Be able to handle multidimensional data, in particular spectra vs time, or
emission vs excitation spectra and even in the future spectral imaging data.
* Permit scripting, using python

## Main things to do

Need to create an initial user interface that is able to handle "spectrum"
objects and multidimensional super-objects.
Then there needs to be the code for the "spectrum" object and its interface.

## User manual

### Installation

How to install the program in various different platforms.

The program can be used as a library and scripted from the command line after
importing, for example, as:
import spectroscopy as spec
into an interactive python session, an Jupyter notebook or a python script.
Alternatively the program be run and used via a GUI.

A simple python script to load a pair of UV-vis spectra and generate a
figure showing the spectra and their difference could be:

    # Setup the program
    import spectroscopy as spc
    # Simple import instruction
    spectrum_one = spc.import("File1.txt")
    spectrum_two = spc.import("File2.txt")
    # Basic maths operations
    difference   = spectrum_one - spectrum_two
    # Make a plot using matplotlib

TODO: Need to imagine this    
1. mise en place, File1.txt and File2.txt two spectra in txt files.
2. imagine code for plotting simply a nice graphic (using internal info)
3. code import (and of course, export, load and save) methods
4. code simple math operations on spectra (+ - / *) with constants and other
   spectra

### Loading data from different sources

TODO: The program should be able to load data from different sources.
Different open source spectroscopy file types and any others that I can find
how to read.

TODO: In the future it would be good if it could load directly from different
machines through an integrated interface (but that is not for tomorrow)

TODO: Editing spectrum metadata

### Spectrum Calculations

TODO: The program should be able to do a number of simple, and less simple,
calculations on spectra, these should be easily scriptable and produce as
output a new spectrum. These include:
1. Simple arithmetic operations (+,-,*,/) between spectra or with a number.
2. Resampling the x-axis, shifting x-axis etc.
3. Smoothing with different algorithms (notable Savitsky-Golay and Fourier)
4. Baseline correction using several algorithms (compare chromatography)
5. Calculation of derivatives and integrals
6. Unit conversions Abs/T/1-T, Fluorescence/corrected (Io or Sens), CD ?

### Spectrum Analysis

TODO: Analysis programs that can produce different types of information:
1. Spectral decomposition into predetermined components, spectra or extinctions
outputing [component, weight] pairs
2. Fitting with different model curves outputing [parameters, weight] list

### Output and figures

TODO: Output of data.
1. Internal format that is as versatile as possible
2. Export as standard open-source formats as complete as possible
3. Visual exploration scripted fro figures or GUI for exploration

### Program Status

TODO: Enquiries about current status of the program
1. List of loaded spectra and their information
2. List of views and what they show
3. About, Help and manual

## API - what can you do with a spectrum object
