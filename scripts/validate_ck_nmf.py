#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Phase 4 validation: reproduce Figure_CK_111225's NMF against the original
notebook code, on the real biofilm ATR-FTIR series.
"""
import os
import warnings

import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from sklearn.decomposition import NMF

warnings.simplefilter("ignore")

from spectroscopy import SpectrumCollection  # noqa: E402
from spectroscopy.processing import multivariate as mv  # noqa: E402

ROOT = os.path.expanduser("~/Documents/Research/Notebook/2025/")
FILES = ([ROOT + "10/27/" + n for n in
          ["PAK_av.txt", "MucA22_av.txt", "MucADP_av.txt", "MucADPplus_av.txt"]]
         + [ROOT + f"11/07/Sample{i}_av.txt" for i in range(1, 9)]
         + [ROOT + f"11/24/Sample{i}_av.txt" for i in range(1, 9)]
         + [ROOT + f"11/25/Sample{i}_av.txt" for i in (5, 6, 7)])
FILES = [f for f in FILES if os.path.exists(f)]
K = 6


def rubberband(x, y):
    pts = np.column_stack((x, y)); hull = ConvexHull(pts); hv = hull.vertices
    pl = np.where(hv == np.argmin(x))[0][0]; pr = np.where(hv == np.argmax(x))[0][0]
    if pl < pr:
        a1 = hv[pl:pr + 1]; a2 = np.concatenate((hv[pr:], hv[:pl + 1]))
    else:
        a1 = hv[pr:pl + 1]; a2 = np.concatenate((hv[pl:], hv[:pr + 1]))
    lo = a1 if np.mean(y[a1]) < np.mean(y[a2]) else a2
    lo = lo[np.argsort(x[lo])]; xu, iu = np.unique(x[lo], return_index=True)
    return np.interp(x, xu, y[lo][iu])


def original():
    aligned = []
    for path in FILES:
        frame = pd.read_csv(path, sep='\t', names=["Wavenumber", "Absorbance"])
        wn, ab = frame.iloc[:, 0].values, frame.iloc[:, 1].values
        region = (wn >= 900) & (wn <= 1800)
        wn, ab = wn[region], ab[region]
        ab = ab - rubberband(wn, ab)
        aligned.append(ab / np.max(ab))
    matrix = np.array(aligned)

    model = NMF(n_components=K, init='nndsvda', max_iter=2000, random_state=0)
    weights = model.fit_transform(matrix)
    components = model.components_
    spacing = np.mean(np.diff(wn))
    areas = components.sum(axis=1) * spacing
    scaled = weights * areas[np.newaxis, :]
    return (matrix, wn,
            1 - np.var(matrix - weights @ components) / np.var(matrix),
            scaled / scaled.sum(axis=1, keepdims=True))


def rewritten():
    collection = (SpectrumCollection.from_files(FILES, technique='ATR-FTIR')
                  .crop(900, 1800)
                  .baseline_correct('rubberband')
                  .normalize('max'))
    return collection, mv.decompose(collection, 'nmf', K)


if __name__ == '__main__':
    old_matrix, old_x, old_var, old_fractions = original()
    collection, result = rewritten()
    new_x, new_matrix = collection.to_matrix()

    print(f"{len(FILES)} spectra, k={K}")
    print(f"  matrix identical        : {np.allclose(old_matrix, new_matrix)} "
          f"(max |diff| {np.max(np.abs(old_matrix - new_matrix)):.2e})")
    print(f"  explained variance      : {old_var:.6f} / {result.explained_variance:.6f}")
    print(f"  contributions identical : "
          f"{np.allclose(old_fractions, result.contributions())} "
          f"(max |diff| {np.max(np.abs(old_fractions - result.contributions())):.2e})")

    print("\n  stability vs k on the real data:")
    for k in (2, 4, 6, 8):
        runs = mv.stability(collection, 'nmf', k, mode='runs', n_trials=10)
        boot = mv.stability(collection, 'nmf', k, mode='bootstrap', n_trials=10)
        print(f"    k={k}: runs {runs.overall:.3f}   bootstrap {boot.overall:.3f}")
