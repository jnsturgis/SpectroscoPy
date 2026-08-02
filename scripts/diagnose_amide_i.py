#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Diagnose an amide I secondary-structure analysis, one page per sample.

Roadmap section 20. Section 19 found that the same protein at four
concentrations gives compositions spread over 20 percentage points while every
fit reports R^2 > 0.998. This script exists so that the afternoon spent finding
out why is spent looking at spectra rather than writing plumbing.

It answers section 20.1's questions in order, and prints the numbers as well as
drawing them, so the summary is readable without opening the PDF.

    1. Is there enough signal?          -> panel A, and the printed table
    2. Do the samples differ physically? -> the water-factor column. If it
       varies across a concentration series, the samples are at different
       stages of drying and are not repeats of one measurement.
    3. Is the subtraction recoverable?   -> panel B, and the factor scan in C
    4. How many components are supported? -> panel D, second derivative
    5. Are the assignments right?        -> panel E, fit and residual

Usage
-----
    python scripts/diagnose_amide_i.py                 # both default series
    python scripts/diagnose_amide_i.py --series bsa
    python scripts/diagnose_amide_i.py --data DIR --samples "X-1:2,X-2:2" \\
        --reference "H2O:4" --label "my series"

Output goes to diagnostics/ beside the data, one PDF per series.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use('Agg')                                   # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

import spectroscopy as spc  # noqa: E402
from spectroscopy.processing import structure  # noqa: E402

#: The two series from roadmap section 18.2, so the script runs with no
#: arguments. ``stem: replicate count``.
DEFAULT_DATA = Path(
    "/hdd/james/Documents/Research/Students/Candice Gomez/FTIR_spectra")
SERIES = {
    'bsa': {
        'label': 'BSA (dry film)',
        'samples': {'BSA-1': 2, 'BSA-2': 2, 'BSA-5': 2, 'BSA-10': 2},
        'reference': None,          # no water band to remove
    },
    'lysozyme': {
        'label': 'lysozyme (solution)',
        'samples': {'lysosyme-1': 3, 'lysosyme-2': 2,
                    'lysosyme-5': 2, 'lysosyme-10': 2},
        'reference': ('H2O', 4),
    },
}

#: Where water absorbs and protein does not. The anchor for the subtraction
#: factor: the OH stretch near 3300 overlaps amide A and is not usable.
WATER_WINDOW = (1950.0, 2250.0)
WATER_BAND = 2130.0

#: A conventional amide I set, plus amide II anchors so the baseline is
#: constrained by both bands (which is why the fit region is wide).
AMIDE_II_POSITIONS = [1515.0, 1545.0, 1570.0]
AMIDE_I_POSITIONS = [1620.0, 1631.0, 1640.0, 1651.0, 1658.0, 1668.0, 1680.0]

FIT_REGION = (1480.0, 1720.0)
BASELINE_REGION = (1400.0, 1800.0)


def average(directory, stem, count):
    """Mean of a sample's replicates."""
    paths = [f"{directory}/{stem}.{index}.dpt" for index in range(count)]
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        raise FileNotFoundError(f"missing replicate(s): {missing}")
    return spc.SpectrumCollection.from_files(paths, technique='ATR-FTIR').mean()


def factor_scan(sample, reference, grid=None):
    """Residual in the water window against subtraction factor."""
    if grid is None:
        grid = np.linspace(0.0, 1.3, 261)
    scores = np.array([
        float(np.sum(sample.subtract_reference(reference, factor=float(f))
                     .crop(*WATER_WINDOW).y ** 2))
        for f in grid
    ])
    return grid, scores


def value_at(spectrum, wavenumber):
    ascending = spectrum.x[0] < spectrum.x[-1]
    xs = spectrum.x if ascending else spectrum.x[::-1]
    ys = spectrum.y if ascending else spectrum.y[::-1]
    return float(np.interp(wavenumber, xs, ys))


def diagnose(directory, samples, reference_spec, label, positions,
             derivative_weight, derivative_span, output):
    """One page per sample, plus a summary page. Returns the printed rows."""
    rows = []
    with PdfPages(output) as pdf:
        for stem, count in samples.items():
            raw = average(directory, stem, count)

            factor, grid, scores = None, None, None
            if reference_spec is not None:
                grid, scores = factor_scan(raw, reference_spec)
                factor = float(grid[int(np.argmin(scores))])
                subtracted = raw.subtract_reference(reference_spec, factor=factor)
            else:
                subtracted = raw

            prepared = (subtracted.crop(*BASELINE_REGION)
                                  .baseline_correct('rubberband'))
            prepared.name = stem
            band = prepared.crop(1600.0, 1700.0)

            detected = band.find_peaks(method='second_derivative')
            noise = float(np.std(prepared.crop(1750.0, 1800.0).y))
            signal = float(np.max(band.y))

            composition, fit = None, None
            try:
                composition = structure.from_ftir(
                    prepared, method='amide-i-curve-fit', region=FIT_REGION,
                    positions=positions, position_tolerance=6.0,
                    derivative_weight=derivative_weight,
                    derivative_span=derivative_span)
                fit = prepared.crop(*FIT_REGION).fit_peaks(
                    positions, position_tolerance=6.0,
                    derivative_weight=derivative_weight,
                    derivative_span=derivative_span)
            except Exception as error:                  # noqa: BLE001
                print(f"  {stem}: fit failed -- {error}", file=sys.stderr)

            rows.append({
                'stem': stem, 'factor': factor, 'signal': signal,
                'noise': noise, 'detected': len(detected),
                'composition': composition,
            })

            # -- the page ---------------------------------------------------
            figure = plt.figure(figsize=(11.7, 8.3))
            figure.suptitle(f"{label} -- {stem}", fontsize=13)
            grid_spec = figure.add_gridspec(3, 2, hspace=0.45, wspace=0.25)

            axis = figure.add_subplot(grid_spec[0, 0])
            axis.plot(raw.x, raw.y, lw=0.8, label='sample')
            if reference_spec is not None:
                axis.plot(reference_spec.x, factor * reference_spec.y, lw=0.8,
                          label=f'reference x {factor:.3f}')
            axis.axvspan(*WATER_WINDOW, color='0.9', zorder=0)
            axis.set_xlim(3800, 900)
            axis.set_title('A. raw, with the water window shaded', fontsize=9)
            axis.legend(fontsize=7, frameon=False)

            axis = figure.add_subplot(grid_spec[0, 1])
            axis.plot(subtracted.x, subtracted.y, lw=0.8)
            axis.axvline(WATER_BAND, color='crimson', lw=0.8, ls=':')
            axis.axhline(0.0, color='0.6', lw=0.5)
            axis.set_xlim(2400, 900)
            axis.set_title('B. after subtraction (2130 marked): flat there?',
                           fontsize=9)

            axis = figure.add_subplot(grid_spec[1, 0])
            if grid is not None:
                axis.plot(grid, scores, lw=0.9)
                axis.axvline(factor, color='crimson', lw=0.8)
                axis.set_yscale('log')
                axis.set_xlabel('subtraction factor')
                axis.set_title(f'C. factor scan -- minimum {factor:.3f}',
                               fontsize=9)
                if factor <= grid[1] or factor >= grid[-2]:
                    axis.set_title(f'C. factor at the edge ({factor:.3f}) -- '
                                   'widen the grid', fontsize=9, color='crimson')
            else:
                axis.text(0.5, 0.5, 'no reference subtracted',
                          ha='center', va='center')
                axis.set_title('C. factor scan', fontsize=9)

            axis = figure.add_subplot(grid_spec[1, 1])
            spacing = abs(float(np.mean(np.diff(band.x))))
            window = max(int(round(derivative_span / spacing)) | 1, 5)
            from spectroscopy.processing import common
            second = common.derivative(band.y, order=2, window_length=window,
                                       polyorder=3, delta=spacing)
            axis.plot(band.x, -second, lw=0.9)
            for position in detected.position:
                axis.axvline(position, color='crimson', lw=0.6, ls=':')
            axis.set_xlim(1700, 1600)
            axis.set_title(f'D. inverted 2nd derivative -- {len(detected)} '
                           f'bands found (amide I wants 5-7)', fontsize=9)

            axis = figure.add_subplot(grid_spec[2, :])
            fitted = prepared.crop(*FIT_REGION)
            axis.plot(fitted.x, fitted.y, lw=1.0, color='k', label='data')
            if fit is not None:
                axis.plot(fit.x, fit.curve(), lw=0.9, color='crimson',
                          label=f'fit R^2={fit.r_squared:.4f}')
                for component in fit.components():
                    axis.plot(fit.x, component, lw=0.5, color='0.5')
                axis.plot(fit.x, fit.residual - 0.15 * np.max(fitted.y),
                          lw=0.7, color='tab:blue', label='residual (offset)')
            axis.set_xlim(FIT_REGION[1], FIT_REGION[0])
            axis.set_xlabel('wavenumber (cm$^{-1}$)')
            axis.legend(fontsize=7, frameon=False, ncol=3)
            if composition is not None:
                summary = "  ".join(
                    f"{category.name} {100 * value:.0f}%"
                    for category, value in composition.fractions.items()
                    if value
                )
                axis.set_title(f'E. fit and components -- {summary}', fontsize=9)

            pdf.savefig(figure)
            plt.close(figure)

        # -- summary page ---------------------------------------------------
        figure, axes = plt.subplots(1, 2, figsize=(11.7, 5.0))
        figure.suptitle(f"{label} -- across the series", fontsize=13)

        names = [row['stem'] for row in rows]
        categories = ['helix', 'sheet', 'turn', 'other', 'aggregated']
        for category in categories:
            values = [100 * (row['composition'].get(category) or 0.0)
                      if row['composition'] else 0.0 for row in rows]
            axes[0].plot(names, values, marker='o', label=category)
        axes[0].set_ylabel('percent')
        axes[0].set_title('composition across the series\n'
                          '(a flat line is what correctness looks like)',
                          fontsize=9)
        axes[0].legend(fontsize=7, frameon=False)
        axes[0].tick_params(axis='x', rotation=30)

        factors = [row['factor'] for row in rows]
        if any(f is not None for f in factors):
            axes[1].plot(names, factors, marker='s', color='crimson')
            axes[1].set_title('water subtraction factor\n'
                              '(varying means the samples differ physically)',
                              fontsize=9)
            axes[1].tick_params(axis='x', rotation=30)
        else:
            axes[1].axis('off')
        figure.tight_layout()
        pdf.savefig(figure)
        plt.close(figure)

    return rows


def report(rows, label):
    """The same information as text, because a table beats a figure for this."""
    print(f"\n=== {label} ===")
    header = (f"{'sample':<16}{'factor':>8}{'amide I':>9}{'noise':>10}"
              f"{'S/N':>8}{'bands':>7}   composition")
    print(header)
    print("-" * len(header))
    for row in rows:
        composition = row['composition']
        factor = '   --   ' if row['factor'] is None else f"{row['factor']:8.3f}"
        ratio = row['signal'] / row['noise'] if row['noise'] else float('inf')
        text = "fit failed"
        if composition is not None:
            text = "  ".join(
                f"{name} {100 * (composition.get(name) or 0):.0f}%"
                for name in ('helix', 'sheet', 'turn', 'other', 'aggregated'))
        print(f"{row['stem']:<16}{factor}{row['signal']:9.4f}"
              f"{row['noise']:10.5f}{ratio:8.0f}{row['detected']:7d}   {text}")

    spreads = {}
    for name in ('helix', 'sheet'):
        values = [100 * (row['composition'].get(name) or 0.0)
                  for row in rows if row['composition']]
        if values:
            spreads[name] = max(values) - min(values)
    if spreads:
        print("\nspread across the series (this is the number to minimise):")
        for name, spread in spreads.items():
            verdict = "usable" if spread < 5 else "NOT usable"
            print(f"  {name:<10}{spread:6.1f} percentage points   {verdict}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--series', choices=[*SERIES, 'all'], default='all')
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--samples',
                        help='override, as "stem:count,stem:count"')
    parser.add_argument('--reference', help='override, as "stem:count"')
    parser.add_argument('--label', default='series')
    parser.add_argument('--out', type=Path, default=None,
                        help='output directory (default: diagnostics/ beside the data)')
    parser.add_argument('--derivative-weight', type=float, default=0.5)
    parser.add_argument('--derivative-span', type=float, default=10.0)
    parser.add_argument('--auto-positions', action='store_true',
                        help='let the second derivative choose the components, '
                             'instead of the conventional set')
    arguments = parser.parse_args(argv)

    if not arguments.data.exists():
        parser.error(f"data directory not found: {arguments.data}")

    output_directory = arguments.out or (arguments.data.parent / 'diagnostics')
    output_directory.mkdir(parents=True, exist_ok=True)

    positions = None if arguments.auto_positions else (
        AMIDE_II_POSITIONS + AMIDE_I_POSITIONS)

    if arguments.samples:
        jobs = {arguments.label: {
            'label': arguments.label,
            'samples': {part.split(':')[0]: int(part.split(':')[1])
                        for part in arguments.samples.split(',')},
            'reference': ((arguments.reference.split(':')[0],
                           int(arguments.reference.split(':')[1]))
                          if arguments.reference else None),
        }}
    elif arguments.series == 'all':
        jobs = SERIES
    else:
        jobs = {arguments.series: SERIES[arguments.series]}

    for key, job in jobs.items():
        reference_spec = None
        if job['reference'] is not None:
            stem, count = job['reference']
            reference_spec = average(arguments.data, stem, count)
        output = output_directory / f"amide_i_{key}.pdf"
        rows = diagnose(arguments.data, job['samples'], reference_spec,
                        job['label'], positions,
                        arguments.derivative_weight, arguments.derivative_span,
                        output)
        report(rows, job['label'])
        print(f"\nfigures: {output}")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
