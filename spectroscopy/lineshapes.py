# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
calc.py this module contains various useful calculations for spectroscopy

"""

import math

import numpy as np


def base( pka, ph ):
    """
    Use Henderson Hasselbach equation to calculate the fraction of base form.
    """
    return 10**(ph-pka)/(1+10**(ph-pka))

def gauss( x_values, posn: float, fwhh: float, ext: float ):
    """
    Calculate a gaussian at the x_values centered at posn and with width fwhh.
    """
    return ext * np.exp(-(math.log(2.))*(2.*(x_values-posn)/fwhh)**2)

def lorentz( x_values, posn: float, fwhh: float, ext: float ):
    """
    Calculate a lorentzian at the x_values centered at posn and with width fwhh.
    """
    return ext * fwhh**2/(fwhh**2+4.0*(x_values-posn)**2)

def spec_comp( x_values, posn:float, fwhh:float, ext: float, fg: float ):
    """
    A pseudo-Voigt component: a Gaussian and a Lorentzian of the same position
    and width, mixed in the ratio ``fg`` : ``1 - fg``.

    At ``x = posn`` the value is ``ext`` for any ``fg``, which is what makes
    ``ext`` mean "peak height" whatever the mixing is.

    .. note::

       Until 0.1.1 the Lorentzian half was multiplied by ``ext`` twice --
       :func:`lorentz` already scales by it -- so the peak height was
       ``fg * ext + (1 - fg) * ext**2``. That is only ``ext`` when ``fg`` is 1
       or ``ext`` is 1, so every mixed component was wrong by a factor of
       ``ext``, and the error vanished in exactly the test case (a unit-height
       Gaussian) most likely to be tried first.
    """
    return fg * gauss(x_values, posn, fwhh, ext) + \
        (1.0 - fg) * lorentz( x_values, posn, fwhh, ext )

def main():
    """
    The main routine
    """
    print('This module is part of the spectroscopy programmes and should be imported not run.')
if __name__ == '__main__':
    main()
