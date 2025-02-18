"""
sidechain - this module caclulates the ir-spectra of protein sidechains.

The code is based on that of Joëlle De Meutter and Eric Goormaghtigh (Eur Biophys J
2021 50: 641–651) and the parameters in that article.

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

# pylint: disable=W0511, C0200, C0103

# TODO: Refactor data array
# TODO: Get default True flags to work.

import argparse
import numpy as np
import matplotlib.pyplot as plt

import calc
import spectroscopy as spc

SideChainData = {
  'D': [
    4.25,               # pKa to use in calculation
    [5, 3, 1, 1],       # number of ir components [high pH, low pH, high pD, low pD]
                        # frequencies for all components in order
    [1598,1570,1472,1421,1395,1729,1456,1410,1584,1713],
                        # FWHH for all components in order
    [44,44,44,44,44,44,44,44,44,44],
                        # intensities for all components in order
    [349/2,402,70,186,351,282,71,131,820,290],      # 349/2 as a shoulder
                        # fraction gaussian for all components in order
    [.3,.3,.3,.3,.3,.6,.8,.8,.8,.8]
  ],
  'E': [3.65,[3,3,1,1],[1570,1451,1404,1728,1454,1417,1567,1706],
    [48,48,48,56,56,56,34,45],[546,63,290,219,23,33,830,280],[.9,.9,.9,.5,.5,.5,.4,.4]],
  'Y': [
    10,[4,4,2,2],
    [1601,1560,1499,1444,1617,1599,1514,1455,1603,1500,1615,1516],
    [14,14,14,14,14,14,10,14,14,14,9,7],
    [320,319,468,63,222,115,241,104,350,650,160,500],
    [.5,.5,.5,.5,.5,.5,.5,.5,.5,.5,0,.4,0,.4]
  ],
  'H': [8.97,[5,3,1,1],[1591,1568,1498,1465,1439,1603,1526,1438,1596,1596],
    [14,14,14,14,14,14,14,14,14,14],[10,97,74,10,30,97,74,30,70,70],
    [.4,.4,.4,.4,.4,.4,.4,.4,.4,.4]],
  'F': [0,[4,1],[1606,1499,1457,1446,1494],[6,6,6,6,6],
    [66,55,35,20,80],[.2,.2,.2,.2,.2]],
  'Q': [0,[5,1],[1672,1610,1523,1452,1411,1635],[32,44,44,44,44,36],
    [360,275,79,72,149,560],[.8,0,0,0,0,.6]],
  'N': [0,[5,1],[1681,1618,1502,1421,1404,1648],[32,44,44,44,44,31],
    [274,187,56,99,103,570],[.8,0,0,0,0,.6]],
  'R': [0,[6,2],[1673,1633,1598,1522,1475,1454,1608,1586],[40,40,40,40,40,40,40,21,22],
    [235,219,151,75,18,28,500,460],[.9,.5,.5,.5,.5,.5,.5,.4,.4]],
  'K': [0,[7,0],[1651,1634,1608,1521,1476,1462,1445],[46,48,48,48,48,48,48],
    [108,242,152,138,36,35,28],[.5,.7,.7,.7,.7,.7,.7]],
  '-': [5.5,[1,1,1,1],[1582,1740,1592,1720],[47,50,32,45],[575,170,830,230],[.4,1,.4,.8]],
  '+': [9.5,[1,2,0,0],[1560,1630,1515],[46,54,60],[450,330,200],[0,.2,0]]
}

PKA   =  0
NCOMP =  1
FREQ  =  2
FWHH  =  3
INTE  =  4
FG    =  5

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
    Calculate the spectrum of an amino-acid residue at a given pH.

    Parameters:
        resid is a one letter amino acid code (or '+' or '-').
        ph    the pH of the solution used to calculate acid and base fractions
        d2o   is this in light or heavy water
        start beginning of wavenumber range
        stop  end of wavenumber range
        res   resolution for the spectrum
    Returns:
        a spectrum
    """
    my_spectrum = spectrum_init(wn_info[0], wn_info[1], wn_info[2] )

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
    weights = weights *  data[INTE]
    for weight, freq, fwhh, fg in zip(
        weights, data[FREQ], data[FWHH], data[FG] ):
        my_spectrum[1] = my_spectrum[1] + weight * calc.spec_comp(my_spectrum[0], freq, fwhh, fg)

    return my_spectrum

def sidechain(sequence, wn_range=(LOW_FREQ, HIGH_FREQ), **kwargs ):
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

    ph = options.get("pH")
    d2o = options.get("D2O")
    res = options.get("res")
    addN = options.get("addN")
    addC = options.get("addC")

    if addN:
        sequence = sequence + '+'
    if addC:
        sequence = sequence + '-'

    composition = get_composition(sequence)
    sum_spectrum = spectrum_init(wn_range[0], wn_range[1], res)

    for residue in composition:
        resid_spectrum = calc_resid_spectrum(residue[0], ph, d2o,
            (wn_range[0], wn_range[1], res)) * residue[1]
        sum_spectrum = spectrum_add(sum_spectrum, resid_spectrum )

    return sum_spectrum

# These routines manipulate spectra and should probably be in another module where
# spectra and their manipulation are described.

def spectrum_init(start, end, step):
    """
    Create a new empty spectrum
    """
    if step is None:
        raise ValueError("Step size can not be None")

    count = int(abs(end-start)/step+1)
    return [np.linspace(start, end, count), np.zeros(count)]

def spectrum_add(spectrum1, spectrum2):
    """
    Add 2 spectra and return the sum
    """
    assert (spectrum1[0]==spectrum2[0]).all()

    new_spectrum = [np.copy(spectrum1[0]),spectrum1[1]+spectrum2[1]]
    return new_spectrum

def spectrum_zero(a_spectrum):
    """
    Is the spectrum a zero baseline?
    """
    return np.max(a_spectrum[1])==0.0 and np.min(a_spectrum[1]) == 0.0

def spectrum_wavenumber(a_spectrum):
    """
    return the list of wavenumbers for the spectrum
    """
    return a_spectrum[0]

def spectrum_absorbance(a_spectrum):
    """
    return the list of absorbance values for the spectrum
    """
    return a_spectrum[1]

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

    sum_spectrum = sidechain(sequence, tuple(args.wn_range), **kwargs)

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

        print(f"Calculating individual components at pH {options.get("pH")}")

        composition = get_composition(sequence)
        components = []

        for residue in composition:
            resid_spectrum = calc_resid_spectrum(residue[0],
                options.get("pH"), options.get("D2O"),
                (args.wn_range[0], args.wn_range[1], options.get("res"))) * residue[1]
            if not spectrum_zero(resid_spectrum):
                components.append([resid_spectrum,residue[0]])

        # Plot the spectra
        fig , ax = plt.subplots()
        for element in components:
            ax.plot(spectrum_wavenumber(element[0]),spectrum_absorbance(element[0]),
                label=element[1])
        ax.plot(spectrum_wavenumber(sum_spectrum),spectrum_absorbance(sum_spectrum),label='Sum')
        ax.set_xlabel('Wavenumber (cm$^{-1}$)')
        ax.set_ylabel('Absorbance')
        fig.legend()
        plt.show()

    if args.output:
        np.savetxt(args.output, np.swapaxes(np.array(sum_spectrum),0,1), delimiter = ',' )

if __name__ == '__main__':
    main()
