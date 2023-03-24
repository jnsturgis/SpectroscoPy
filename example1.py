#! /usr/bin/python3
""" A simple example script using the most basic operability

This simple script should be almost the Hello World example for the API.
"""

# Setup the program
import spectroscopy as spec

# Simple import instruction
spectrum_one = spec.load("data/File1.txt")
spectrum_two = spec.load("data/File2.txt")

# Basic maths operations
difference   = spectrum_one - spectrum_two

# Make a plot using matplotlib
spec.display([[spectrum_one, spectrum_two ], difference ])

