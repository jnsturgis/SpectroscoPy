
""" File to allow context independent execution of tests

To use in individual test files start with the line:
'from .context import spectra'
This should work independent of the installation method see:
https://docs.python-guide.org/writing/structure/
"""

import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                '../SpectroscoPy')))

import spectra
