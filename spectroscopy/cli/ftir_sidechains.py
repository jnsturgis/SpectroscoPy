#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Command-line front end for the protein side-chain IR spectrum calculation.

Usage:
    spc-ftir-sidechains [-h] [--wn_range X Y] [--pH PH] [--res RES] [--D2O]
                        [--addN] [--addC] [--plot] [--output OUTPUT] filename

All the science is in :mod:`spectroscopy.processing.ftir`; this module only
parses arguments, reads the sequence file, and optionally plots or saves.
"""

import argparse

import numpy as np

from spectroscopy.processing.ftir import (
    DEF_PH,
    DEF_RES,
    HIGH_FREQ,
    LOW_FREQ,
    calc_resid_spectrum,
    ftir_sidechain,
    get_composition,
)


def _parse_args(argv=None):
    """Build the argument parser and parse ``argv`` (defaults to sys.argv)."""
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
    return parser.parse_args(argv)


def _plot_components(sequence, args, kwargs, sum_spectrum):
    """Plot the summed side-chain spectrum together with its components."""
    import matplotlib.pyplot as plt  # pylint: disable=C0415

    options = {
        "pH":   DEF_PH,
        "D2O":  False,
        "addN": True,
        "addC": True,
        "res":  DEF_RES,
    }
    options.update({k: v for k, v in kwargs.items() if k in options})

    if options.get("addN"):
        sequence = sequence + '+'
    if options.get("addC"):
        sequence = sequence + '-'

    components = []
    for residue in get_composition(sequence):
        resid_spectrum = calc_resid_spectrum(residue[0],
            options.get("pH"), options.get("D2O"),
            (args.wn_range[0], args.wn_range[1], options.get("res"))) * residue[1]
        if np.any(resid_spectrum.y):
            components.append(resid_spectrum)

    fig, ax = plt.subplots()
    for element in components:
        element.plot(ax)

    sum_spectrum.plot(ax)
    ax.set_xlabel(sum_spectrum.x_label)
    ax.set_ylabel(sum_spectrum.y_label)
    fig.legend()
    plt.show()


def main(argv=None):
    """Entry point for the ``spc-ftir-sidechains`` command."""
    args = _parse_args(argv)

    kwargs = {k: v for k, v in vars(args).items()
        if k not in ["filename", "wn_range"] and v is not None}

    with open(args.filename, encoding="UTF-8") as handle:
        sequence = handle.read()

    sum_spectrum = ftir_sidechain(sequence, tuple(args.wn_range), **kwargs)

    if args.plot:
        _plot_components(sequence, args, kwargs, sum_spectrum)

    if args.output:
        sum_spectrum.fileinfo['NAME'] = args.output
        sum_spectrum.save()

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
