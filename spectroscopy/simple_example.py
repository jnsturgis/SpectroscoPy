#! /usr/bin/python3
""" A simple example script using the most basic operability

This simple script should be almost the Hello World example for the API.
"""

# Setup the program
import numpy as np
import matplotlib.pyplot as plt
import spectra as spectra

# Simple import instruction
spectrum_one = spectra.load("../data/Spectrum1.csv")
spectrum_two = spectra.load("../data/Spectrum2.csv")

spectrum_one.write_csv("../data/toto.csv")

# Basic maths operations
difference   = spectrum_one - spectrum_two
result = spectra.smooth(difference, ['savitsky-golay', 3, 9])
decomposition = spectra.decompose(result, [rules])

# Make a plot using matplotlib
spectra.display([[spectrum_one, spectrum_two ], difference ])
print(result)
