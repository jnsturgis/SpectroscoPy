#!/usr/bin/env python3
"""
Phase 1 validation: re-implement the FTIR_sugars pipeline against the new API
and check it gives the same numbers as the original notebook code.

The original (Notebook/2025/10/FTIR_sugars.ipynb) does, per sample:
    load N replicates -> average by hard-coded index -> crop 900-1800
    -> subtract water * hand factor -> rubberband baseline -> normalise to the
    max in 1050-1080 -> second-derivative peak pick
"""
import warnings
import numpy as np

warnings.simplefilter("ignore", DeprecationWarning)

from scipy.signal import find_peaks, savgol_filter        # noqa: E402
from scipy.spatial import ConvexHull                      # noqa: E402

import spectroscopy as spc                                # noqa: E402
from spectroscopy import SpectrumCollection               # noqa: E402

ROOT = "/home/james/Documents/Research/Notebook/2025/10/"


# --------------------------------------------------------------------------
# ORIGINAL notebook code, verbatim
# --------------------------------------------------------------------------

def rubberband(x, y):
    x = np.asarray(x); y = np.asarray(y)
    pts = np.column_stack((x, y))
    hull = ConvexHull(pts); hv = hull.vertices
    pos_l = np.where(hv == np.argmin(x))[0][0]
    pos_r = np.where(hv == np.argmax(x))[0][0]
    if pos_l < pos_r:
        arc1 = hv[pos_l:pos_r + 1]
        arc2 = np.concatenate((hv[pos_r:], hv[:pos_l + 1]))
    else:
        arc1 = hv[pos_r:pos_l + 1]
        arc2 = np.concatenate((hv[pos_l:], hv[:pos_r + 1]))
    lower_arc = arc1 if np.mean(y[arc1]) < np.mean(y[arc2]) else arc2
    lower_arc = lower_arc[np.argsort(x[lower_arc])]
    xu, iu = np.unique(x[lower_arc], return_index=True)
    yu = y[lower_arc][iu]
    return np.interp(x, xu, yu)


def original(files, water_files, water_factor):
    spectra = []
    for filename in files:
        s = spc.Spectrum("", ROOT + filename, 'tsv')      # the OLD route
        spectra.append(s)
    total = spectra[0] - spectra[0]
    for s in spectra:
        total = total + s / len(spectra)

    waters = [spc.Spectrum("", ROOT + f, 'tsv') for f in water_files]
    wtotal = waters[0] - waters[0]
    for w in waters:
        wtotal = wtotal + w / len(waters)

    my = spc.Spectrum(total)
    region = (my.x >= 900) & (my.x <= 1800)
    my.x, my.y = my.x[region], my.y[region]
    wx, wy = wtotal.x[region], wtotal.y[region]

    my.y = my.y - water_factor * wy
    my.y = my.y - rubberband(my.x, my.y)
    norm_region = (my.x >= 1050) & (my.x <= 1080)
    my.y = my.y / my.y[norm_region].max()

    n2d = -savgol_filter(my.y, 10, 3, deriv=2)
    peaks, _ = find_peaks(n2d, height=0.00001, distance=10, prominence=0.0001)
    return my.x, my.y, my.x[peaks]


# --------------------------------------------------------------------------
# NEW API
# --------------------------------------------------------------------------

def rewritten(pattern, water_pattern, water_factor):
    water = SpectrumCollection.from_files(ROOT + water_pattern,
                                          technique='ATR-FTIR').mean()
    result = (SpectrumCollection.from_files(ROOT + pattern, technique='ATR-FTIR')
              .mean()
              .crop(900, 1800)
              .subtract_reference(water.crop(900, 1800), factor=water_factor)
              .baseline_correct('rubberband')
              .normalize('max', window=(1050, 1080)))
    peaks = result.find_peaks(height=0.00001, distance=10, prominence=0.0001)
    return result, peaks


CASES = [
    ("Glucose",  "13/Glucose.[1-6].dpt",  ["14/H2O.0.dpt", "14/H2O.1.dpt", "14/H2O.2.dpt"], 0.7),
    ("Chitin",   "13/Chitin.[1-2].dpt",   ["14/H2O.0.dpt", "14/H2O.1.dpt", "14/H2O.2.dpt"], 0.0),
    ("Alginate", "10/Alginate.1[2-4].dpt", ["14/H2O.0.dpt", "14/H2O.1.dpt", "14/H2O.2.dpt"], 0.1),
]

print("=" * 74)
print("Phase 1 validation -- FTIR_sugars pipeline, old code vs new API")
print("=" * 74)

import glob                                               # noqa: E402
all_ok = True
for name, pattern, water_files, factor in CASES:
    files = [p[len(ROOT):] for p in sorted(glob.glob(ROOT + pattern))]
    if not files:
        print(f"\n{name}: no files matched {pattern}, skipped")
        continue

    ox, oy, opeaks = original(files, water_files, factor)
    new, npeaks = rewritten(pattern, "14/H2O.[0-2].dpt", factor)

    same_x = np.allclose(ox, new.x)
    same_y = np.allclose(oy, new.y, atol=1e-9)
    same_p = np.array_equal(np.sort(opeaks), np.sort(npeaks.position))

    print(f"\n{name}  ({len(files)} replicates, water factor {factor})")
    print(f"  x axis identical      : {same_x}")
    print(f"  y values identical    : {same_y}"
          f"   (max |diff| = {np.max(np.abs(oy - new.y)):.3g})")
    print(f"  peaks identical       : {same_p}"
          f"   ({len(opeaks)} old vs {len(npeaks)} new)")
    if not same_p:
        print(f"      old: {np.round(np.sort(opeaks), 0)}")
        print(f"      new: {np.round(np.sort(npeaks.position), 0)}")
    all_ok &= same_x and same_y and same_p

    print("  recorded history:")
    for number, step in enumerate(new.history, 1):
        print(f"      {number}. {step}")

print("\n" + "=" * 74)
print("ALL CASES IDENTICAL" if all_ok else "DIFFERENCES FOUND -- see above")
print("=" * 74)
