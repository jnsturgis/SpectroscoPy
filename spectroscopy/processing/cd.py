# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Secondary structure from a CD spectrum: several methods, and a way to compare
them.

Estimating structure from a far-UV CD spectrum means writing the measured
spectrum as a combination of reference spectra whose structures are known, and
reading the composition off the combination. Every published method does that;
they differ in *how they cope with there being more references than the data
can distinguish*, which is the whole difficulty.

A typical reference set has 70-130 proteins. A measurement from 190 to 240 nm
carries about 50 independent numbers. Fitting 130 unknowns to 50 equations has
infinitely many exact solutions, and picking one by least squares gives an
excellent residual and an arbitrary answer. The methods below are four
different answers to that.

The methods
-----------
``all-references``
    Plain non-negative least squares against every reference. Correct only
    when the data can actually determine them, and this refuses when it
    cannot. Included as the honest baseline rather than as a recommendation.

``nearest-shapes``
    Do not invert anything. Rank the references by how much the *shape*
    resembles the measurement -- cosine similarity on unit-normalised spectra
    -- and average the known structures of the closest few, weighted by
    similarity. There is no underdetermined system, so there is no arbitrary
    solution to pick.

``subset-average``
    The idea behind CDSSTR (Johnson 1999; see docs/references.md): repeatedly fit small
    random subsets of references, keep the solutions that pass a
    self-consistency test, and average them. Small subsets are determined by
    the data even when the whole set is not.

``ridge``
    Tikhonov-regularised non-negative least squares over every reference. The
    standard statistical treatment of an underdetermined system: penalise the
    size of the coefficients so that one solution is preferred among the many
    that fit.

Which needs matched units, and why it matters
---------------------------------------------
A reference set is in delta epsilon per residue. A measurement is in
millidegrees until somebody supplies a concentration, a path length and a
residue count. Methods that constrain the fitted weights to sum to one -- the
classical self-consistency test -- are comparing an amplitude, so they need
those three numbers to be right. Methods that normalise instead are blind to
amplitude and work on the shape alone.

============================  ==================  ==========================
method                        needs the amplitude  what it is sensitive to
============================  ==================  ==========================
``all-references``            no                   conditioning
``nearest-shapes``            no                   how well the set covers
                                                   the fold
``subset-average``            only with            the acceptance threshold
                              ``sum_to_one=True``
``ridge``                     no                   the penalty
============================  ==================  ==========================

That table is not a detail. Getting a path length wrong by a factor of ten is
easy and silent, and it moves any amplitude-dependent answer by the same
factor while leaving the spectrum looking entirely normal.

Choosing between them
---------------------
By held-out validation, not by which one agrees with a protein you already
know. :func:`benchmark` hides a fraction of the reference set, estimates those
proteins from the rest, and compares with their known structures. Run it on
your own set before trusting any of these on your own protein --
:doc:`the guide </guide/cd-methods>` reports what it gives on SP175 and
SMP180.
"""

from __future__ import annotations

import numpy as np

__all__ = ['METHODS', 'estimate', 'benchmark', 'design_matrix',
           'DEFAULT_REGION', 'DEFAULT_NEIGHBOURS', 'DEFAULT_SUBSET_SIZE']

#: Far-UV range used unless told otherwise. The lower limit is where most
#: instruments run out of light in an ordinary buffer; the reference sets
#: themselves reach 175 nm and are worth using to their limit when the
#: measurement does.
DEFAULT_REGION = (190.0, 240.0)

#: Neighbours averaged by ``nearest-shapes``. Held-out validation on SMP180
#: puts 3, 5 and 8 within 0.005 of each other on mean rmse and shows a slow
#: decline above that as more distant shapes are averaged in; five keeps a
#: usable spread estimate without reaching for them.
DEFAULT_NEIGHBOURS = 5

#: References per draw in ``subset-average``. Eight is CDSSTR's choice.
DEFAULT_SUBSET_SIZE = 8


def design_matrix(spectrum, basis, region=DEFAULT_REGION):
    """
    ``(x, measured, design)`` on a common axis at the **coarser** spacing.

    Coarser because a measurement cannot carry more information about a basis
    than the basis has. Resampling a 1 nm reference set onto a 0.1 nm spectrum
    makes a matrix of 439 rows whose rank is 46: ten times as many equations
    as there is information, which makes an underdetermined fit look well
    posed to anything that counts rows.
    """
    x = np.asarray(spectrum.x, dtype=float)
    spacings = [float(np.median(np.abs(np.diff(np.asarray(b.x, dtype=float)))))
                for b in basis]
    step = max(spacings + [float(np.median(np.abs(np.diff(x))))])

    low, high = sorted(region)
    low = max(low, float(x.min()), *[float(b.x.min()) for b in basis])
    high = min(high, float(x.max()), *[float(b.x.max()) for b in basis])
    if high - low < 5 * step:
        raise ValueError(
            f"the spectrum, the basis and the region {sorted(region)} overlap "
            f"over only {max(high - low, 0):g} nm, which is not enough to "
            f"compare shapes"
        )
    grid = np.arange(low, high + step / 2, step)
    measured = np.asarray(spectrum.resample(grid).y, dtype=float)
    design = np.column_stack([np.asarray(b.resample(grid).y, dtype=float)
                              for b in basis])
    return grid, measured, design


def _unit(values):
    norm = float(np.linalg.norm(values))
    return np.asarray(values, dtype=float) / (norm or 1.0)


def _nnls(design, measured):
    from scipy.optimize import nnls  # noqa: PLC0415

    scale = float(np.max(np.abs(measured))) or 1.0
    weights, _ = nnls(design / scale, measured / scale)
    return weights


def _relative_rms(measured, fitted):
    span = float(np.max(measured) - np.min(measured)) or 1.0
    return float(np.sqrt(np.mean((measured - fitted) ** 2)) / span)


# -- the methods ------------------------------------------------------------
#
# Each takes the measured spectrum, the design matrix, and the table of
# reference compositions (n_references x n_categories), and returns
# (fractions over categories, quality dict).

def _all_references(measured, design, table, **options):
    rank = int(np.linalg.matrix_rank(design))
    if rank < design.shape[1]:
        raise ValueError(
            f"{design.shape[1]} references cannot be determined from data of "
            f"rank {rank}. Least squares would still return an answer -- one "
            f"of infinitely many that fit equally well, with an excellent "
            f"residual. Use 'nearest-shapes', 'subset-average' or 'ridge', "
            f"or fewer references."
        )
    weights = _nnls(design, measured)
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("the fit put no weight on any reference")
    return (table.T @ weights) / total, {
        'rmsd_relative': _relative_rms(measured, design @ weights),
        'condition_number': float(np.linalg.cond(design)),
        'weight_sum': total,
    }


def _nearest_shapes(measured, design, table, *, neighbours=DEFAULT_NEIGHBOURS,
                    **options):
    target = _unit(measured)
    similarity = np.array([float(target @ _unit(design[:, column]))
                           for column in range(design.shape[1])])
    order = np.argsort(similarity)[::-1][:max(1, int(neighbours))]
    weights = np.clip(similarity[order], 0.0, None)
    weights = weights / (weights.sum() or 1.0)
    fractions = weights @ table[order]
    return fractions, {
        'n_neighbours': len(order),
        'best_similarity': float(similarity[order[0]]),
        'worst_similarity': float(similarity[order[-1]]),
        # Disagreement among the neighbours, which is the honest uncertainty:
        # if the closest shapes differ about sheet content, so does the answer.
        'spread': table[order].std(axis=0).tolist(),
    }


def _subset_average(measured, design, table, *, subset_size=DEFAULT_SUBSET_SIZE,
                    draws=2000, tolerance=0.03, sum_to_one=False, seed=0,
                    **options):
    generator = np.random.default_rng(seed)
    n_references = design.shape[1]
    size = min(int(subset_size), n_references)
    accepted, residuals = [], []
    for _ in range(int(draws)):
        pick = generator.choice(n_references, size=size, replace=False)
        weights = _nnls(design[:, pick], measured)
        total = float(weights.sum())
        if total <= 0:
            continue
        # The classical self-consistency test. It compares an amplitude, so it
        # is only meaningful when the measurement is in the reference set's
        # units -- hence off by default.
        if sum_to_one and not 0.9 <= total <= 1.1:
            continue
        residual = _relative_rms(measured, design[:, pick] @ weights)
        if residual > tolerance:
            continue
        accepted.append((table[pick].T @ weights) / total)
        residuals.append(residual)

    if not accepted:
        raise ValueError(
            f"no subset of {size} references fitted within {tolerance:.1%} of "
            f"range in {draws} draws"
            + (". With sum_to_one=True the weights must also sum to about 1, "
               "which requires the spectrum to be in the reference set's "
               "units -- delta epsilon per residue, not millidegrees."
               if sum_to_one else ".")
        )
    accepted = np.array(accepted)
    return accepted.mean(axis=0), {
        'n_accepted': len(accepted),
        'n_draws': int(draws),
        'subset_size': size,
        'median_rmsd_relative': float(np.median(residuals)),
        'spread': accepted.std(axis=0).tolist(),
    }


def _ridge(measured, design, table, *, penalty=0.05, **options):
    from scipy.optimize import nnls  # noqa: PLC0415

    scale = float(np.max(np.abs(measured))) or 1.0
    scaled = design / scale
    largest = float(np.linalg.svd(scaled, compute_uv=False)[0]) or 1.0
    strength = float(penalty) * largest
    augmented = np.vstack([scaled,
                           strength * np.eye(design.shape[1])])
    target = np.concatenate([measured / scale, np.zeros(design.shape[1])])
    weights, _ = nnls(augmented, target)
    total = float(weights.sum())
    if total <= 0:
        raise ValueError("the penalised fit put no weight on any reference")
    return (table.T @ weights) / total, {
        'penalty': float(penalty),
        'effective_strength': strength,
        'rmsd_relative': _relative_rms(measured, design @ weights),
        'n_used': int(np.sum(weights > 1e-8)),
    }


#: The available methods. Names are explicit because which one produced a
#: number is part of the number's meaning (ADR-0002 section 7).
METHODS = {
    'all-references': _all_references,
    'nearest-shapes': _nearest_shapes,
    'subset-average': _subset_average,
    'ridge': _ridge,
}

#: Which methods compare an amplitude rather than only a shape.
NEEDS_AMPLITUDE = {'subset-average': 'only with sum_to_one=True'}


def _table(compositions, categories):
    return np.array([[composition.fractions.get(category, 0.0) or 0.0
                      for category in categories]
                     for composition in compositions], dtype=float)


def estimate(spectrum, method, basis, compositions, *, region=DEFAULT_REGION,
             **options):
    """
    Estimate a composition by one named method. See :data:`METHODS`.

    Returns a :class:`~spectroscopy.processing.structure.Composition`; the
    method's own diagnostics are in its ``quality``.
    """
    from spectroscopy.processing.structure import Composition  # noqa: PLC0415

    if method not in METHODS:
        raise ValueError(
            f"unknown CD method {method!r}; available are "
            f"{sorted(METHODS)}"
        )
    categories = list(compositions[0].fractions)
    _, measured, design = design_matrix(spectrum, basis, region)
    fractions, quality = METHODS[method](measured, design,
                                         _table(compositions, categories),
                                         **options)
    quality = {**quality, 'method': method, 'n_references': len(basis)}
    return Composition(fractions=dict(zip(categories, map(float, fractions))),
                       method=method, technique=spectrum.technique or 'CD',
                       quality=quality, source=spectrum.name)


def benchmark(basis, compositions, methods=None, *, folds=10,
              region=DEFAULT_REGION, seed=0, options=None):
    """
    Hold out part of the reference set and estimate it from the rest.

    **The only honest way to choose between these methods.** Comparing them on
    one protein whose structure you happen to know selects whichever agrees,
    and any of them will agree with something. Here every reference protein is
    estimated in turn from a set that does not contain it, and the errors are
    against structures that were known before the fit.

    Parameters
    ----------
    basis, compositions : sequence
        The reference set, as :func:`~spectroscopy.library.load_dichroweb_basis`
        returns it.
    methods : sequence of str, optional
        Which to compare. All of :data:`METHODS` by default.
    folds : int
        How many parts to split the set into. Ten means each round hides 10%
        and estimates it from the other 90%.
    options : dict, optional
        ``{method: {keyword: value}}``, passed through per method.

    Returns
    -------
    dict
        ``{method: {category: {'rmse', 'bias', 'mae'}, ..., 'n': int,
        'failures': int}}``. **Bias is worth as much as rmse**: a method that
        pulls every protein towards the average composition of the reference
        set has a bias that grows with how unusual the protein is, and that is
        exactly when an estimate matters.
    """
    methods = list(methods or METHODS)
    options = options or {}
    categories = list(compositions[0].fractions)
    truth = _table(compositions, categories)

    generator = np.random.default_rng(seed)
    order = generator.permutation(len(basis))
    parts = np.array_split(order, int(folds))

    errors = {method: [] for method in methods}
    failures = {method: 0 for method in methods}
    for part in parts:
        keep = [index for index in range(len(basis)) if index not in set(part)]
        pool = [basis[index] for index in keep]
        pool_table = truth[keep]
        pool_compositions = [compositions[index] for index in keep]
        for held in part:
            _, measured, design = design_matrix(basis[held], pool, region)
            for method in methods:
                try:
                    fractions, _ = METHODS[method](
                        measured, design, pool_table,
                        **options.get(method, {}))
                except (ValueError, RuntimeError):
                    failures[method] += 1
                    continue
                errors[method].append(np.asarray(fractions) - truth[held])
        del pool_compositions

    summary = {}
    for method in methods:
        rows = np.array(errors[method]) if errors[method] else np.zeros((0, len(categories)))
        entry = {'n': len(rows), 'failures': failures[method]}
        for index, category in enumerate(categories):
            column = rows[:, index] if len(rows) else np.array([np.nan])
            entry[category.name] = {
                'rmse': float(np.sqrt(np.mean(column ** 2))),
                'bias': float(np.mean(column)),
                'mae': float(np.mean(np.abs(column))),
            }
        summary[method] = entry
    return summary
