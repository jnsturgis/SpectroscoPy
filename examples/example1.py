""" A simple example script using the most basic operability

This simple script should be almost the Hello World example for the API.
"""

# Find the files
import sys
import os
sys.path.append(os.path.abspath("/home/james/src/SpectrscoPy/SpectroscoPy"))

# Setup the program
import numpy as np
import matplotlib.pyplot as plt
import spectra

# Simple import instruction
spectrum_one = spectra.import("../data/File1.txt")
spectrum_two = spectra.import("../data/File2.txt")

# Basic maths operations
difference   = spectrum_one - spectrum_two

# Make a plot using matplotlib
spactra.display([[spectrum_one, spectrum_two ], difference ])
