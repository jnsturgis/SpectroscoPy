""" A simple example script using the most basic operability

This simple script should be almost the Hello World example for the API.
"""

    # Setup the program
    import spectroscopy as spc

    # Simple import instruction
    spectrum_one = spc.import("File1.txt")
    spectrum_two = spc.import("File2.txt")

    # Basic maths operations
    difference   = spectrum_one - spectrum_two

    # Make a plot using matplotlib
