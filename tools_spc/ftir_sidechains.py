"""
sidechain - this module caclulates the ir-spectra of protein sidechains.

The code is based on that of Joëlle De Meutter and Eric Goormaghtigh 2021 (Eur
Biophys J 50: 641–651) and the parameters in that article which in turn are based
on the spectral parameters from Venyaminov S. Yu. and Kalnin N.N. 1990

Note that the pKa values used by Goormaghtigh are sometimes a bit strange, and
perhaps the ionization of Cyteine and Lysine should be taken into account. Though
Cys runs into the problem of oxidised or not.

Amino acid    Code  pKa1    pKa2    pKa3     pI     pKa3(here)
Cysteine  	    C 	1.96 	10.28 	8.18 	5.07 	 N.A !!!
Aspartic acid  	D 	1.88 	9.60 	3.65 	2.77 	 4.25
Glutamic acid  	E 	2.19 	9.67 	4.25 	3.22 	 3.65
Lysine      	K 	2.18 	8.95 	10.53 	9.74 	 N.A.
Arginine  	    R 	2.17 	9.04 	12.48 	10.76 	 N.A.
Histidine    	H 	1.82 	9.17 	6.00 	7.59 	 8.97
Tyrosine        Y   2.20    9.11                    10.00

Data from: Carey and Giuliano (2011) Amino acids, peptides and proteins.
Organic Chemistry 8 th Edition, 25, 1126 McGraw Hill, ISBN-13: 978-0077354770

AA             Value        Here
Asp          3.5 ± 1.2      4.25
Glu          4.2 ± 0.9      3.65
His          6.6 ± 1.0      8.97
Cys          6.8 ± 2.7      N.A. !!!
Tyr         10.3 ± 1.2     10.00
Lys         10.5 ± 1.1      N.A.
C-terminus   3.3 ± 0.8      5.5
N-terminus   7.7 ± 0.5      9.5

Data from Grimsley, G. R., Scholtz, J. M. & Pace, C. N. A summary of the
measured pK values of the ionizable groups in folded proteins.
Protein Sci 18, 247–251, doi: 10.1002/pro.19 (2009).

================================================================================

Proposed solution.
Add a switch to use original values but default to pKa's from Grimsley et al.
as experimental from folded proteins.

================================================================================

When called as an independant program the user can specify using the command line
parameters, the protein sequence and pH, the wavenumber range of interest, if
the protein is in light or heavy water (H2O or D2O), and what output to produce
form the options, a plot, a saved plot as a pdf file and the spectrum in a csv
file.

This programme usage:
sidechains.py [-h] [--wn_range X Y] [--pH PH] [--res RES] [--D2O]
            [--addN] [--addC] [--plot] [--output OUTPUT] filename

Process command-line arguments.

Positional argument:
  filename         File name for sequence read as fasta format

options:
  -h, --help       show this help message and exit
  --wn_range X Y   An optional 2-tuple giving the wavenumber range, default 1300 2000.
  --pH PH          Optional pH value, default 7.0
  --res RES        Optional spectral resolution, default 1.0
  --D2O            Optional switch if D2O used, default False
  --addN           Optional switch to include N-term contribution, default True
  --addC           Optional switch to include C-term contribution, default True
  --plot           Optional switch to plot spectrum for sequence and individual components.
  --output OUTPUT  An optional filename for saving the spectrum in csv format.


Desired API sidechain function
spec.ir.sidechain( sequence, (wavenumbers), **kwargs)

    sequence      - is a fast format string representing the protein sequence.
    (wavenumbers) - is an optional 2 tuple of fist and last wavenumber and
                    defaults to (1300.0, 2000.0).
    kwargs
    pH            - (Float) the pH to use for the calculations defaults to 7.0
    d2o           - (Boolean) calculate for deuterated protein defaults to False
    addN          - (Boolean) add the N-teminus primary amine contribution.
    addC          - (Boolean) add the C-teminus free carboxyl contribution.
    res           - (Float) resolution defaults to 1 cm-1

Return a spectrum object containing the calculated side chain spectrum of the
protein given the sequence.
"""

# TODO Use Grimsley pKa values as default, switch for Goormaghtigh

# pylint: disable=W0511

import argparse
import numpy as np
import matplotlib.pyplot as plt

import calc
import spectroscopy as spc

SideChainData = {
  'D': [
    4.25,               # pKa to use in calculation
    [5, 3, 1, 1],       # number of ir components [high pH, low pH, high pD, low pD]
                        # Spectral components an array of tuples (freq, fwhh, ext, fg)
    [   (1598, 44, 175, 0.3), (1570, 44, 402, 0.3), (1472, 44, 70, 0.3),
        (1421, 44, 186, 0.3), (1395, 44, 351, 0.3),
        (1729, 44, 282, 0.6), (1456, 44,  71, 0.8), (1410, 44, 131, 0.8),
        (1584, 44, 820, 0.8), (1713, 44, 290, 0.8) ]],
  'E': [3.65,[3,3,1,1],
    [   (1570, 48, 546, 0.9), (1451, 48,  63, 0.9), (1404, 48, 290, 0.9),
        (1728, 56, 219, 0.5), (1454, 56,  23, 0.5), (1417, 56,  33, 0.5),
        (1567, 34, 830, 0.4), (1706, 45, 280, 0.4) ]],
  'Y': [10.0,[4,4,2,2],
    [   (1601, 14, 320, 0.5), (1560, 14, 319, 0.5), (1499, 14, 468, 0.5), (1444, 14,  63, 0.5),
        (1617, 14, 222, 0.5), (1599, 14, 115, 0.5), (1514, 10, 241, 0.5), (1455, 14, 104, 0.5),
        (1603, 14, 350, 0.0), (1500, 14, 650, 0.4),
        (1615,  9, 160, 0.0), (1516,  7, 500, 0.4) ]],
  'H': [8.97,[5,3,1,1],
    [   (1591, 14,  10, 0.4), (1568, 14,  97, 0.4), (1498, 14,  74, 0.4),
        (1465, 14,  10, 0.4), (1439, 14,  30, 0.4),
        (1603, 14,  97, 0.4), (1526, 14,  74, 0.4), (1438, 14,  30, 0.4),
        (1596, 14,  70, 0.4), (1596, 14,  70, 0.4) ]],
  'F': [0.00,[4,1],
    [   (1606,  6,  66, 0.2), (1499,  6,  55, 0.2), (1457,  6,  35, 0.2), (1446,  6,  80, 0.2),
        (1497,  6,  80, 0.2) ]],
  'Q': [0.00,[5,1],
    [   (1672, 32, 360, 0.8), (1610, 44, 275, 0.0), (1523, 44,  79, 0.0),
        (1452, 44,  72, 0.0), (1411, 44, 149, 0.0), (1635, 36, 560, 0.6) ]],
  'N': [0.00,[5,1],
    [   (1681, 32, 274, 0.8), (1618, 44, 187, 0.0), (1502, 44,  56, 0.0),
        (1421, 44,  99, 0.0), (1404, 44, 103, 0.0), (1648, 31, 570, 0.6) ]],
  'R': [0.00,[6,2],
    [   (1673, 40, 235, 0.9), (1633, 40, 219, 0.5), (1598, 40, 151, 0.5),
        (1522, 40,  75, 0.5), (1475, 40,  18, 0.5), (1454, 40,  28, 0.5),
        (1608, 21, 500, 0.4), (1586, 22, 460, 0.4) ]],
  'K': [0.00,[7,0],
    [   (1651, 46, 108, 0.5), (1634, 48, 242, 0.7), (1608, 48, 152, 0.7), (1521, 48, 138, 0.7),
        (1476, 48,  36, 0.7), (1462, 48,  35, 0.7), (1445, 48,  28, 0.7) ]],
  '-': [5.50,[1,1,1,1],
    [   (1582, 47, 575, 0.4), (1740, 50, 170, 1.0), (1592, 32, 830, 0.4), (1720, 45, 230, 0.8) ]],
  '+': [9.50,[1,2,0,0],
    [   (1560, 46, 450, 0.0), (1630, 54, 330, 0.2), (1515, 60, 200, 0.0) ]]
}

PKA   =  0
NCOMP =  1
COMP  =  2
FREQ  =  0
FWHH  =  1
INTE  =  2
FG    =  3

# pKa values from Grimsley
Grimsley = {'D':3.5,'E':4.2,'H':6.6,'Y':10.3,'K':10.5,'-':3.3,'+':7.7}

# Default values
LOW_FREQ  = 1300
HIGH_FREQ = 2000
DEF_PH    = 7.0
DEF_RES   = 1.0

def get_composition( sequence: str ):
    """
    From a sequence string of single-letter amino-acid codes calculate the composition.

    The single letter codes are augmented by + for an unblocked N-terminal
    and - for an unblocked C terminal. The composition is returned as an array
    of lists in the form ['Code', count ].
    The routine can also handle fasta files so all lines starting with '>'
    are ignored.
    """
    composition = []
    ignore = False
    for res in sequence:
        if res == '>':
            ignore = True
        elif res== '\n':
            ignore = False
        elif res.isspace():
            pass
        elif not ignore:
            found = False
            for item in composition:
                if item[0] == res:
                    item[1] = item[1] + 1
                    found = True
            if not found:
                composition.append([res,1])

    return composition

def calc_resid_spectrum(resid, ph:float, d2o: bool, wn_info ):
    """
    Calculate the expected spectrum of an amino-acid residue at a given pH.

    Parameters:
        resid   is a one letter amino acid code (or '+' or '-' for N- and C- terms).
        ph      the pH of the solution used to calculate acid and base fractions
        d2o     is this in light or heavy water
        wn_info is a 3 tuple containing:
            start   beginning of wavenumber range in cm-1
            stop    end of wavenumber range in cm-1
            res     resolution for the spectrum in cm-1
    Returns:
        a spectrum
    """
    my_spectrum = spc.Spectrum()
    count = int(abs(wn_info[1]-wn_info[0])/wn_info[2]+1)
    my_spectrum.x_data = np.linspace(wn_info[0], wn_info[1], count)
    my_spectrum.y_data = np.zeros(count)
    my_spectrum.name = resid
    my_spectrum.x_label = 'Wavenumber (cm$^{-1}$)'
    my_spectrum.y_label = 'Absorbance'

    try:
        data = SideChainData[resid]
    except KeyError:
        return my_spectrum

    if data[PKA] != 0:
        base_form = calc.base(data[PKA],ph)
        h_weights = base_form * np.ones(data[NCOMP][0])
        temp      = (1.0-base_form) * np.ones(data[NCOMP][1])
        h_weights = np.concatenate((h_weights,temp))
        d_weights = base_form * np.ones(data[NCOMP][2])
        temp      = (1.0-base_form) * np.ones(data[NCOMP][3])
        d_weights = np.concatenate((d_weights,temp))
    else:
        h_weights = np.ones(data[NCOMP][0])
        d_weights = np.ones(data[NCOMP][1])

    weights = np.concatenate(( h_weights * int(not d2o), d_weights * int(d2o)))
    for weight, component in zip(weights, data[COMP]):
        my_spectrum.y_data = my_spectrum.y_data + \
                    weight * calc.spec_comp(my_spectrum.x_data, component[FREQ],
                        component[FWHH], component[INTE], component[FG])
    return my_spectrum

def ftir_sidechain(sequence, wn_range=(LOW_FREQ, HIGH_FREQ), **kwargs ):
    """
    Calculate the side chain spectrum for a protein sequence at the given pH in h2o or d2o.

    See API definition above.
    """

    if not isinstance(wn_range, tuple) or len(wn_range) != 2:
        raise ValueError(spc.WN_RANGE_ERR)

    # Set defaults
    options = {
        "pH":   DEF_PH,
        "D2O":  False,
        "addN": True,
        "addC": True,
        "res":  DEF_RES
    }
    # Parse kwargs
    if kwargs:
        options.update({k: v for k, v in kwargs.items() if k in options})

    if not isinstance(options.get("pH"), float ):
        raise ValueError(spc.PH_FLOAT_ERR)
    if not isinstance(options.get("res"), float ):
        raise ValueError(spc.RES_FLOAT_ERR)
    if not isinstance(options.get("D2O"), bool ):
        raise ValueError(spc.D2O_BOOL_ERR)
    if not isinstance(options.get("addN"), bool ):
        raise ValueError(spc.ADDN_BOOL_ERR)
    if not isinstance(options.get("addC"), bool ):
        raise ValueError(spc.ADDC_BOOL_ERR)

    ph    = options.get("pH")
    d2o   = options.get("D2O")
    res   = options.get("res")

    if options.get("addN"):
        sequence = sequence + '+'
    if options.get("addC"):
        sequence = sequence + '-'

    composition = get_composition(sequence)
    sum_spectrum = spc.Spectrum()
    count = int(abs(wn_range[1]-wn_range[0])/res+1)
    sum_spectrum.x_data  = np.linspace(wn_range[0], wn_range[1], count)
    sum_spectrum.y_data  = np.zeros(count)
    sum_spectrum.sum     = 'Sum'
    sum_spectrum.x_label = 'Wavenumber (cm$^{-1}$)'
    sum_spectrum.y_label = 'Absorbance'

    for residue in composition:
        resid_spectrum = calc_resid_spectrum(residue[0], ph, d2o,
            (wn_range[0], wn_range[1], res)) * residue[1]
        sum_spectrum = sum_spectrum + resid_spectrum

    return sum_spectrum

# These routines are private and belong with main only to get data from the
# command line.

def main():
    """
    If run as a program use the command line to determine the parameters,
    for the pH, d2o/h2o and a sequence also wavenumber range and to output
    a figure and/or a csv file.
    """

    parser = argparse.ArgumentParser(description="Process command-line arguments.")

    parser.add_argument("filename", type=str,
        help="File name for sequence or '-' for stdin.")
    parser.add_argument("--wn_range", type=int, nargs=2,
        default=(LOW_FREQ, HIGH_FREQ),
        metavar=('X', 'Y'), help="An optional 2-tuple giving the wavenumber range.")
    parser.add_argument("--pH", type=float, help="Optional pH value, default 7.0")
    parser.add_argument("--res", type=float,
        help="Optional spectral resolution, default 1.0")
    parser.add_argument("--D2O", action="store_true",
        help="Optional switch if D2O used, default False")
    parser.add_argument("--addN", action="store_false",
        help="Optional switch to include N-term contribution, default True")
    parser.add_argument("--addC", action="store_false",
        help="Optional switch to include C-term contribution, default True")
    parser.add_argument("--plot", action="store_true",
        help="Optional switch to plot spectrum for sequence and individual components.")
    parser.add_argument("--output", type=str,
        help="An optional filename for saving the spectrum in csv format.")
    args = parser.parse_args()
    # Convert argparse.Namespace to a dictionary for **kwargs
    kwargs = {k: v for k, v in vars(args).items()
        if k not in ["filename", "wn_range"] and v is not None}

    # Call the function with parsed arguments
    with open(args.filename, encoding="UTF-8") as f:
        sequence = f.read()

    sum_spectrum = ftir_sidechain(sequence, tuple(args.wn_range), **kwargs)

    if args.plot:
        options = {
            "pH":   DEF_PH,
            "D2O":  False,
            "addN": True,
            "addC": True,
            "res":  DEF_RES
        }
        # Parse kwargs
        options.update({k: v for k, v in kwargs.items() if k in options})

        if options.get("addN"):
            sequence = sequence + '+'
        if options.get("addC"):
            sequence = sequence + '-'

        composition = get_composition(sequence)
        components = []

        for residue in composition:
            resid_spectrum = calc_resid_spectrum(residue[0],
                options.get("pH"), options.get("D2O"),
                (args.wn_range[0], args.wn_range[1], options.get("res"))) * residue[1]
            if np.any(resid_spectrum.y_data):
                components.append(resid_spectrum)

        # Plot the spectra
        fig , ax = plt.subplots()
        for element in components:
            element.plot(ax)

        sum_spectrum.plot(ax)
        ax.set_xlabel(sum_spectrum.x_label)
        ax.set_ylabel(sum_spectrum.y_label)
        fig.legend()
        plt.show()

    if args.output:
        np.savetxt(args.output, np.swapaxes(np.array(sum_spectrum),0,1), delimiter = ',' )

if __name__ == '__main__':
    main()
