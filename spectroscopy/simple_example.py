#! /usr/bin/python3
""" A simple example script using the most basic operability

This simple script should be almost the Hello World example for the API.
"""

# Setup the program
import spectra

# Simple import instruction
spectrum_one = spectra.Spectrum("../data/uvvis_spectra/Spectrum1.csv", "csv")
spectrum_two = spectra.Spectrum("../data/uvvis_spectra/Spectrum2.csv", "xml")

spectra.write_csv_file("../data/toto.csv",spectrum_one)

# Basic maths operations
difference   = spectrum_one - spectrum_two
result = difference.smooth('savitsky-golay', 3, 9)

# Make a plot using matplotlib
spectra.display([[spectrum_one, spectrum_two ], difference ])
print(result)
