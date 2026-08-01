# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
""" Module containing the different messages given by the program

"""

ABOUT = ''.join(("                   SpectroscoPy                           \n",
                 "         Copyright (C) 2023-2026 James N. Sturgis         \n",
                 "      This program comes with ABSOLUTELY NO WARRANTY;     \n",
                 "This is free software, and you are welcome to redistribute\n",
                 " it under the terms of the Mozilla Public License v2.0.   "))

# Error messages from spectra.py module

SPEC_MATH_ERR = "Unsupported type for spectrum maths."

# Error messages from processing.ftir

WN_RANGE_ERR  = "wn_range must be a 2-tuple (start_freq, end_freq)"
PH_FLOAT_ERR  = "pH must be a float"
RES_FLOAT_ERR = "res (spectrum resolution) must be a float"
D2O_BOOL_ERR  = "D2O must be a boolean"
ADDN_BOOL_ERR = "addN must be a boolean"
ADDC_BOOL_ERR = "addC must be a boolean"
