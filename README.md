# SpectroscoPy1.0

## Introduction

This project is born out of frustration as a scientist that uses a variety of
spectroscopy instruments but frequently runs into problems with idiosyncratic
and low quality analysis software, proprietary file formats etc. 

The objective of this project is therefor to make my life simpler providing
a common framework for handling and analysing spectra from different sources.
I my work I use mostly UV-visible, absorption fluorescence and circular dichroism 
spectra and ATR-FTIR absorption spectra so these spectroscopies will probably
be over represented.
I envisage a three tier system in the organization of the project.
1. Stand alone tools to do specific tasks, working esentially like most unix/linux
tools.
2. A library of functions, covering the same tools, that can be used for scripting
analyses and working on them for example in a jupyter notebook.
3. A user friendly programme driven via a GUI.

## Main things to do

### Managing data from different sources

TODO: The program should be able to load data from different sources.
Different open source spectroscopy file types and any others that I can find
how to read.

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
3. Calculating different sorts of synthetic spectrum.

### Output and figures

TODO: Output of data.
1. Internal format that is as versatile as possible
2. Export as standard open-source formats as complete as possible
3. Visual exploration scripted fro figures or GUI for exploration
## User manual

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
    difference.plot()

### Program Status

TODO: Enquiries about current status of the program
1. List of loaded spectra and their information
2. List of views and what they show
3. About, Help and manual

## API - what can you do with a spectrum object
