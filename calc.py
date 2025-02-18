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

def gauss( x_values, posn, fwhh ):
    """
    Calculate a gaussian at the x_values centered at posn and with width fwhh.
    """
    return np.exp(-(math.log(2.))*(2.*(x_values-posn)/fwhh)**2)

def lorentz( x_values, posn, fwhh ):
    """
    Calculate a lorentzian at the x_values centered at posn and with width fwhh.
    """
    return fwhh**2/(fwhh**2+4.0*(x_values-posn)**2)

def spec_comp( x_values, posn, fwhh, fg ):
    """
    Calculate a compont
    """
    return fg * gauss(x_values, posn, fwhh) + (1.0 - fg) * lorentz( x_values, posn, fwhh )

def main():
    """
    The main routine
    """

if __name__ == '__main__':
    main()
