""" A simple example script using the most basic operability

This simple script should be almost the Hello World example for the API.
"""

# import sys
#
# sys.path.append(r'//home/james/src/SpectroscoPy/src')

# Setup the program
import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc

# Simple import instruction
spectrum_one = spc.import("../data/File1.txt")
spectrum_two = spc.import("../data/File2.txt")

# Basic maths operations
difference   = spectrum_one - spectrum_two

# Make a plot using matplotlib
spc.display([[spectrum_one, spectrum_two ], difference ])
