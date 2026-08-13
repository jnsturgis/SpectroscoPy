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
           'uncertainty_from_replicates', 'NEEDS_AMPLITUDE',
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


def _references(references, argument='references'):
    """
    Check that what arrived is a reference set, and say so plainly if not.

    Every method here needs each spectrum's known structure. Taking that as a
    second, parallel list is what ADR-0004 removed: a list of compositions in
    the wrong order raises nothing and shifts the answer by 0.16 in helix.
    """
    from spectroscopy.library import ReferenceSet  # noqa: PLC0415

    if isinstance(references, ReferenceSet):
        return references
    raise TypeError(
        f"{argument} must be a library.ReferenceSet, not "
        f"{type(references).__name__}. A set of spectra and a separate list "
        f"of their compositions can fall out of step without raising "
        f"anything, so the pairing is made once, where it can be checked: "
        f"library.load_dichroweb_basis(directory), library.load_basis("
        f"manifest), or ReferenceSet.from_compositions(spectra, compositions)."
    )


def _truth(references):
    """``(categories, table)`` -- the known structure of a set, as a matrix."""
    compositions = references.compositions
    if not compositions:
        raise ValueError("the reference set is empty")
    categories = references.categories or list(compositions[0].fractions)
    return categories, _table(compositions, categories)


def uncertainty_from_replicates(collection):
    """
    Per-wavelength standard error from repeat scans.

    CD is nearly always measured as several accumulations, so the noise is
    **measured, not assumed** -- and it is strongly wavelength-dependent, since
    it grows wherever the photomultiplier is working hardest. Feeding that to
    :func:`estimate` weights each wavelength by how well it is known, and lets
    the uncertainty on the composition be propagated rather than guessed.

    Returns a :class:`~spectroscopy.spectra.Spectrum` of standard errors,
    which is just ``collection.sem()``; this exists to say so in one place.
    """
    return collection.sem()


def estimate(spectrum, method, references, *, region=DEFAULT_REGION,
             sigma=None, resamples=0, seed=0, **options):
    """
    Estimate a composition by one named method. See :data:`METHODS`.

    Parameters
    ----------
    references : ReferenceSet
        Spectra of proteins whose structure is known, carrying that structure
        with them -- from :func:`~spectroscopy.library.load_dichroweb_basis`,
        :func:`~spectroscopy.library.load_basis`, or
        :meth:`~spectroscopy.library.ReferenceSet.from_compositions`. A
        structural basis and a set of reference proteins are the same type
        here; the second is the first with a table that is not one-hot.
    sigma : Spectrum or array_like, optional
        Per-wavelength standard error, as
        :func:`uncertainty_from_replicates` returns from a set of repeat
        scans. Given, the fit is **weighted** by ``1/sigma``: a wavelength
        known to ten times the precision of another counts for ten times as
        much, which is the right thing to do and not what an unweighted fit
        does. Far-UV CD noise rises steeply towards the blue as the detector
        starves, so the weighting is far from uniform in practice.
    resamples : int
        With ``sigma``, refit this many times on the spectrum perturbed by its
        own measured noise, and report the spread as
        ``quality['uncertainty']``. **This is the honest error bar**: it says
        how much the answer moves for noise the size you actually measured,
        which no goodness-of-fit statistic can tell you.

    Returns a :class:`~spectroscopy.processing.structure.Composition`; the
    method's own diagnostics are in its ``quality``.
    """
    from spectroscopy.processing.structure import Composition  # noqa: PLC0415

    if method not in METHODS:
        raise ValueError(
            f"unknown CD method {method!r}; available are "
            f"{sorted(METHODS)}"
        )
    references = _references(references)
    categories, table = _truth(references)
    grid, measured, design = design_matrix(spectrum, references, region)

    errors = None
    if sigma is not None:
        errors = np.asarray(
            sigma.resample(grid).y if hasattr(sigma, 'resample')
            else np.interp(grid, grid, np.asarray(sigma, dtype=float)),
            dtype=float)
        if errors.shape != measured.shape:
            raise ValueError(
                f"sigma has {errors.shape} values for {measured.shape} "
                f"wavelengths"
            )
        floor = float(np.median(errors[errors > 0])) if np.any(errors > 0) else 1.0
        errors = np.where(errors > 0, errors, floor)

    def solve(values):
        if errors is None:
            return METHODS[method](values, design, table, **options)
        weight = 1.0 / errors
        return METHODS[method](values * weight, design * weight[:, None],
                               table, **options)

    fractions, quality = solve(measured)
    quality = {**quality, 'method': method, 'n_references': len(references),
               'weighted': errors is not None}
    # Which reference set an answer came from is part of the answer: the same
    # spectrum against SP175 and against SMP180 is two results, not one.
    provenance = references.info.get('source') or references.name
    if provenance:
        quality['references'] = provenance

    if errors is not None and int(resamples) > 0:
        generator = np.random.default_rng(seed)
        spread = []
        for _ in range(int(resamples)):
            perturbed = measured + errors * generator.normal(size=measured.size)
            try:
                candidate, _ = solve(perturbed)
            except (ValueError, RuntimeError):
                continue
            spread.append(candidate)
        if spread:
            spread = np.array(spread)
            quality['uncertainty'] = dict(zip(
                [c.name for c in categories], spread.std(axis=0).tolist()))
            quality['n_resamples'] = len(spread)

    return Composition(fractions=dict(zip(categories, map(float, fractions))),
                       method=method, technique=spectrum.technique or 'CD',
                       quality=quality, source=spectrum.name)


def benchmark(references, methods=None, *, folds=10,
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
    references : ReferenceSet
        The reference set, as :func:`~spectroscopy.library.load_dichroweb_basis`
        returns it. Each protein's structure travels with its spectrum, so
        hiding a fold cannot separate the two.
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
    references = _references(references)
    categories, truth = _truth(references)

    generator = np.random.default_rng(seed)
    order = generator.permutation(len(references))
    parts = np.array_split(order, int(folds))

    errors = {method: [] for method in methods}
    failures = {method: 0 for method in methods}
    for part in parts:
        held_out = set(part)
        keep = [index for index in range(len(references))
                if index not in held_out]
        pool = [references[index] for index in keep]
        pool_table = truth[keep]
        for held in part:
            _, measured, design = design_matrix(references[held], pool, region)
            for method in methods:
                try:
                    fractions, _ = METHODS[method](
                        measured, design, pool_table,
                        **options.get(method, {}))
                except (ValueError, RuntimeError):
                    failures[method] += 1
                    continue
                errors[method].append(np.asarray(fractions) - truth[held])

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


# -- SELCON ------------------------------------------------------------------
#
# The self-consistent method of Sreerama & Woody (1993), with the variable
# selection and solution rules of SELCON3 (Sreerama & Woody 2000). See
# docs/references.md.
#
# What is taken from the literature, and what is not, is set out in the
# docstring of _selcon below -- the primary papers are paywalled and the open
# re-implementations are NonCommercial, so parts of this are the published
# description rather than the published code, and the parts that are choices
# are named as choices.

#: SELCON3's rule on individual fractions: a solution may go slightly negative
#: but not far. Quoted in the literature as -0.025.
SELCON_MIN_FRACTION = -0.025

#: SELCON3's sum rule: fractions must sum to within 5% of one. **This compares
#: an amplitude**, so it only means anything when the spectrum is in the
#: reference set's units.
SELCON_SUM_BOUNDS = (0.95, 1.05)

#: Fit residual a solution must beat, relative to the spectrum's range.
SELCON_MAX_RMSD = 0.25


def _truncated_solutions(design, table, measured, limit):
    """
    ``f = F A+ a`` at every truncation from 1 to ``limit``, in one pass.

    SELCON collects solutions across truncations rather than choosing one, so
    the same matrix is inverted once per retained singular value. Decomposing
    it once and varying only the truncation is the difference between this
    being usable and not: on SMP180's 128 references the naive form -- one
    ``np.linalg.svd`` per truncation, of a matrix that had not changed -- took
    **281 seconds** for a single estimate against 5.

    **``limit`` is capped at the number of singular values.** Asking for more
    truncations than the matrix has does not produce more solutions; it
    produces the full-rank solution again. Collecting it repeatedly weights it
    more heavily in every average downstream, which is what the earlier
    per-truncation version did: on SMP180's 128 references over 51 wavelengths
    it counted the full-rank answer about eighty times per subset, so roughly
    half of SELCON's "accepted solutions" were one solution. Nothing in the
    published description asks for that -- it was the bound on a loop.

    Yields ``(keep, fractions, weights)``.
    """
    left, values, right = np.linalg.svd(design, full_matrices=False)
    projected = left.T @ measured
    threshold = values[0] * 1e-12 if len(values) else 0.0

    inverse = np.zeros_like(values)
    for keep in range(1, min(int(limit), len(values)) + 1):
        index = keep - 1
        if values[index] > threshold:
            inverse[index] = 1.0 / values[index]
        weights = right.T @ (inverse * projected)
        yield keep, table.T @ weights, weights


def _pseudo_inverse_solution(design, table, measured, keep):
    """One truncation of :func:`_truncated_solutions`, for a single answer."""
    for _, fractions, weights in _truncated_solutions(design, table, measured,
                                                      max(1, int(keep))):
        result = (fractions, weights)
    return result


def _selcon(measured, design, table, *, min_references=5, max_references=None,
            iterations=10, convergence=1e-4, sum_bounds=SELCON_SUM_BOUNDS,
            min_fraction=SELCON_MIN_FRACTION, max_rmsd=SELCON_MAX_RMSD,
            ignore_sum_rule=False, **options):
    """
    SELCON: the self-consistent method with variable selection.

    From the published description of Sreerama & Woody:

    1. **Self-consistency.** The unknown's own spectrum is added to the basis
       set carrying a *guess* at its structure. The augmented system is solved
       by singular value decomposition, the guess is replaced by the solution,
       and the process repeats until it stops moving. Including the unknown is
       what makes the method self-consistent, and it is the step that
       distinguishes SELCON from a plain fit.
    2. **Variable selection.** References are ordered by how close they are to
       the query, and solutions are sought using increasing numbers of the
       closest ones, rather than the whole set at once.
    3. **Selection rules.** A candidate is kept only if its fractions sum to
       within 5% of one and none is below -0.025. Surviving solutions are
       averaged.

    Deliberately named choices, because the primary papers are paywalled and
    the open re-implementations are NonCommercial, so these could not be read
    off either:

    * how many singular values to retain -- solutions are collected across
      every truncation from one to the subset size, which is the behaviour the
      description implies rather than a number anyone stated;
    * the residual threshold, ``max_rmsd``;
    * the convergence test on the self-consistent loop.

    .. warning::

       **SELCON needs the spectrum in the reference set's units, and there is
       no honest way round it.** The self-consistent step puts the query into
       the basis *as a column beside the references*, so its magnitude
       relative to them is part of the model. Measured on a synthetic set
       whose references peak near 30: the same shape scaled by 0.1, 1, 10 and
       800 gave helix 0.350, 0.562, 0.558 and 0.550, and the unscaled run
       refused outright at 0.1 and 10 while succeeding at 1 and 800.

       ``ignore_sum_rule=True`` drops the sum rule and renormalises, which
       lets a spectrum in millidegrees produce an answer. It does **not** make
       the method scale-free -- the drift above was measured with the flag on.
       Use it to get a number when you have no concentration, understanding
       that the number moves with an amplitude you have not pinned down. For a
       genuinely amplitude-blind estimate use ``ridge`` or ``nearest-shapes``,
       and convert with
       :meth:`~spectroscopy.spectra.Spectrum.to_mean_residue_ellipticity`
       when you can.
    """
    import warnings  # noqa: PLC0415

    if ignore_sum_rule:
        warnings.warn(
            "SELCON with ignore_sum_rule=True drops one of its three "
            "selection rules and still is not scale-free: the query sits in "
            "the basis beside the references, so the answer moves with an "
            "amplitude you have not supplied. Prefer ridge or nearest-shapes "
            "when the concentration is unknown.",
            UserWarning, stacklevel=3)
    n_references = design.shape[1]
    upper = min(int(max_references or n_references), n_references)
    lower = max(2, min(int(min_references), upper))

    target = measured

    # Variable selection: order the references by closeness to the query, on
    # shape, so that the ordering does not itself depend on the amplitude.
    query_direction = _unit(target)
    closeness = np.array([float(query_direction @ _unit(design[:, j]))
                          for j in range(n_references)])
    order = np.argsort(closeness)[::-1]

    accepted, diagnostics = [], []
    for size in range(lower, upper + 1):
        pick = order[:size]
        sub_design, sub_table = design[:, pick], table[pick]

        # The self-consistent loop. The unknown joins the basis carrying a
        # guess; the guess is whatever the closest references average to.
        guess = sub_table.mean(axis=0)
        for _ in range(int(iterations)):
            augmented = np.column_stack([sub_design, target])
            augmented_table = np.vstack([sub_table, guess])
            limit = min(size + 1, augmented.shape[1])
            solutions = [fractions for _, fractions, _
                         in _truncated_solutions(augmented, augmented_table,
                                                 target, limit)]
            updated = np.mean(solutions, axis=0)
            if np.max(np.abs(updated - guess)) < convergence:
                guess = updated
                break
            guess = updated

        # Apply the rules to every truncation, not only to the converged mean:
        # SELCON collects solutions, it does not return one.
        augmented = np.column_stack([sub_design, target])
        augmented_table = np.vstack([sub_table, guess])
        limit = min(size + 1, augmented.shape[1])
        for keep, fractions, weights in _truncated_solutions(
                augmented, augmented_table, target, limit):
            total = float(np.sum(fractions))
            if ignore_sum_rule:
                if abs(total) < 1e-9:
                    continue
                fractions = fractions / total
            elif not sum_bounds[0] <= total <= sum_bounds[1]:
                continue
            if float(np.min(fractions)) < min_fraction:
                continue
            residual = _relative_rms(target, augmented @ weights)
            if residual > max_rmsd:
                continue
            accepted.append(fractions)
            diagnostics.append((size, keep, residual))

    if not accepted:
        raise ValueError(
            "no SELCON solution passed the selection rules: fractions summing "
            f"to {sum_bounds[0]}-{sum_bounds[1]} with none below "
            f"{min_fraction}. The sum rule compares an **amplitude**, so this "
            "is what happens when the spectrum is in millidegrees and the "
            "basis in delta epsilon. Convert the spectrum with "
            "Spectrum.to_mean_residue_ellipticity. Passing "
            "ignore_sum_rule=True will produce a number, but it will be one "
            "that depends on the amplitude you have not supplied."
        )

    accepted = np.array(accepted)
    fractions = accepted.mean(axis=0)
    # SELCON reports the fractions renormalised, having required them to sum
    # to about one already.
    total = float(fractions.sum()) or 1.0
    return fractions / total, {
        'n_solutions': len(accepted),
        'subset_sizes': sorted({size for size, _, _ in diagnostics}),
        'spread': accepted.std(axis=0).tolist(),
        'median_rmsd_relative': float(np.median([r for _, _, r in diagnostics])),
        'sum_before_renormalising': total,
        'ignore_sum_rule': bool(ignore_sum_rule),
    }


METHODS['selcon'] = _selcon
NEEDS_AMPLITUDE['selcon'] = 'always -- the query joins the basis'
