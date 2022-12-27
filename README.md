# SpectrscoPy

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

### Spectrum Calculations

### Spectrum Analysis

### Output and figures

## API - what can you do with a spectrum object
