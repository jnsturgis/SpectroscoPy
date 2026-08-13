# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Protein secondary structure, from spectra.

    >>> from spectroscopy.processing import structure

One output type, several inputs. :func:`from_ftir` estimates from an amide I
band, :func:`from_cd` from a far-UV CD spectrum against reference proteins,
and :func:`helix_from_theta222` from a single CD ordinate. All return the same
:class:`Composition` against the same vocabulary, which is what makes two
estimates of one sample comparable -- :meth:`Composition.compare` does that
comparison on the categories both methods can express.

:func:`nearest_references` answers a different and often better question: not
"what is this made of" but "what does it most look like, and what are those".
:func:`cd_shape_descriptors` reports band positions and ratios, which are
amplitude-free and need no reference set at all.

**DSSP is the baseline vocabulary** (ADR-0002). Every category is declared as
the set of DSSP states it claims, so comparing two methods is an operation on
set partitions rather than a table of judgement calls. See
:cite:`kabsch1983dssp` -- or, until this project has a citation extension,
``docs/references.md``.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

__all__ = [
    'Category', 'Composition', 'DSSP_STATES', 'AMIDE_I_BANDS',
    'from_ftir', 'FTIR_METHODS',
    'from_cd', 'CD_METHODS', 'helix_from_theta222',
    'cd_shape_descriptors', 'nearest_references',
]

#: The eight DSSP states, and what they mean. The vocabulary everything else
#: is defined against.
DSSP_STATES = {
    'H': 'alpha-helix',
    'G': '3-10 helix',
    'I': 'pi-helix',
    'E': 'extended strand in a beta-sheet',
    'B': 'isolated beta-bridge',
    'T': 'hydrogen-bonded turn',
    'S': 'bend',
    '-': 'none of the above',
}


@dataclass(frozen=True)
class Category:
    """
    One structural category, declared as the DSSP states it covers.

    Parameters
    ----------
    name : str
        What to call it -- ``'helix'``, ``'regular-helix'``, ``'aggregated'``.
    states : frozenset of str
        The DSSP states this category claims. **Empty means the category has
        no DSSP equivalent at all**, which is the honest answer for
        intermolecular aggregation: DSSP runs on one chain, and aggregation is
        a quaternary feature. Such categories are excluded from comparison
        rather than folded into a structural one.
    note : str, optional
        Why, when the mapping is not exact. The CD reference sets split helix
        into regular and distorted by position within the element, which is a
        subdivision DSSP does not encode; both carry ``states={'H'}`` and a
        note saying so.
    """

    name: str
    states: frozenset = frozenset()
    note: str | None = None

    def __str__(self) -> str:
        return self.name

    @property
    def is_structural(self) -> bool:
        """False for categories with no DSSP equivalent, e.g. aggregation."""
        return bool(self.states)


#: Amide I band assignments, as bands rather than a fixed mapping, because the
#: boundaries differ between authors and the choice belongs to the caller
#: (ADR-0002 section 7). Ranges are contiguous so that a fitted component
#: cannot fall into a gap and be silently lost.
#:
#: The aggregation band is deliberately first and deliberately not structural:
#: intermolecular beta-sheet absorbs near 1615-1625 cm^-1 with a weak partner
#: near 1690, and reporting it as sheet would turn a failed preparation into a
#: structural result.
AMIDE_I_BANDS = (
    (Category('aggregated', frozenset(),
              note='intermolecular beta-sheet; no DSSP state, DSSP runs on '
                   'one chain'), (1600.0, 1620.0)),
    (Category('sheet', frozenset({'E', 'B'})), (1620.0, 1640.0)),
    (Category('other', frozenset({'S', '-'})), (1640.0, 1648.0)),
    (Category('helix', frozenset({'H', 'G', 'I'})), (1648.0, 1660.0)),
    (Category('turn', frozenset({'T'})), (1660.0, 1680.0)),
    (Category('sheet', frozenset({'E', 'B'})), (1680.0, 1700.0)),
)

#: Estimators available to :func:`from_ftir`. Named at the call site, never
#: defaulted: which one is used is a scientific choice, not an implementation
#: detail (ADR-0002 section 7).
FTIR_METHODS = ('amide-i-curve-fit',)


@dataclass
class Composition:
    """
    A secondary structure estimate.

    The common output of every estimator, whatever the technique.

    Attributes
    ----------
    fractions : dict
        ``{Category: fraction}``. A fraction of ``None`` means *this method
        cannot estimate this category* -- never zero, which would be a claim
        the method never made.
    method, technique : str
        How, and from what. Both appear in :meth:`compare` output, because
        the provenance of a disagreement is most of its meaning.
    quality : dict
        How much to trust it. Method-specific, but never empty: an estimate
        that cannot say how good it is invites being read as fact.
    source : str, optional
        The sample or spectrum it came from.
    """

    fractions: dict
    method: str
    technique: str
    quality: dict = field(default_factory=dict)
    source: str | None = None

    def __len__(self) -> int:
        return len(self.fractions)

    def __repr__(self) -> str:
        where = f" of {self.source}" if self.source else ""
        return f"<Composition{where}: {self.method} on {self.technique}>"

    def __str__(self) -> str:
        lines = [f"{self.technique} / {self.method}"
                 + (f" -- {self.source}" if self.source else "")]
        for category, fraction in self.fractions.items():
            if fraction is None:
                shown = "not estimated"
            else:
                shown = f"{100 * fraction:5.1f} %"
            marker = "" if category.is_structural else "   (not a DSSP state)"
            lines.append(f"  {category.name:<14} {shown}{marker}")
        if self.quality:
            summary = ", ".join(
                f"{key}={value:.4g}" if isinstance(value, (int, float))
                else f"{key}={value}"
                for key, value in self.quality.items()
                if not isinstance(value, (list, tuple, np.ndarray))
            )
            if summary:
                lines.append(f"  quality: {summary}")
        return "\n".join(lines)

    # -- access -------------------------------------------------------------

    def get(self, name, default=None):
        """The fraction of the category called ``name``."""
        for category, fraction in self.fractions.items():
            if category.name == name:
                return fraction
        return default

    @property
    def estimated(self) -> dict:
        """Only the categories this method actually estimated."""
        return {category: fraction
                for category, fraction in self.fractions.items()
                if fraction is not None}

    # -- comparison ---------------------------------------------------------

    def compare(self, other) -> Comparison:
        """
        Compare with another estimate, on whatever both can express.

        Two methods rarely use the same categories. Rather than lining up the
        names and hoping -- which is what everyone does, and which is wrong by
        the 3-10 and pi content between amide I and the CD reference sets --
        this merges categories on both sides until the two partitions of the
        DSSP alphabet agree, then compares the merged groups.

        Categories with no DSSP states (aggregation) and fractions of ``None``
        take no part, and are reported separately rather than treated as zero.

        Returns
        -------
        Comparison
        """
        return _compare(self, other)


@dataclass
class Comparison:
    """
    What two estimates agree and disagree about.

    Attributes
    ----------
    groups : dict
        ``{group name: (fraction_a, fraction_b)}`` over the coarsest set of
        categories both estimates can express. **This is the rigorous
        comparison**, and it is sometimes disappointingly coarse -- see
        ``caveats``.
    nominal : dict
        ``{name: (fraction_a, fraction_b)}`` matching categories by name, with
        positional variants (``regular-``/``distorted-``) summed into their
        parent. **This is the comparison everybody actually makes**, and it is
        approximate exactly to the extent ``caveats`` describes. Offered
        because refusing to show it does not stop anyone making it.
    caveats : list of str
        Why the rigorous grouping is coarser than the nominal one: which DSSP
        states one method assigns somewhere the other does not. This is the
        actionable part -- it names the systematic error in the comparison a
        reader was about to make anyway.
    rmsd : float
        Root mean square difference over ``groups``.
    excluded : dict
        What took no part, and why -- non-structural categories, and
        categories a method did not estimate.
    a, b : Composition
    """

    groups: dict
    rmsd: float
    excluded: dict
    a: Composition
    b: Composition
    nominal: dict = field(default_factory=dict)
    caveats: list = field(default_factory=list)

    def __str__(self) -> str:
        lines = [f"{self.a.technique}/{self.a.method}  vs  "
                 f"{self.b.technique}/{self.b.method}"]

        if self.nominal:
            lines.append("by name (approximate -- see caveats):")
            lines.append(f"{'':<20}{'first':>9}{'second':>9}{'diff':>8}")
            for name, (first, second) in sorted(self.nominal.items()):
                lines.append(f"  {name:<18}{100 * first:8.1f}%"
                             f"{100 * second:8.1f}%{100 * (second - first):+8.1f}")

        lines.append("strictly comparable:")
        for name, (first, second) in self.groups.items():
            shown = name if len(name) <= 40 else name[:37] + "..."
            lines.append(f"  {shown:<18}{100 * first:8.1f}%"
                         f"{100 * second:8.1f}%{100 * (second - first):+8.1f}")
        lines.append(f"  RMSD {100 * self.rmsd:.1f} percentage points")

        for caveat in self.caveats:
            lines.append(f"  ! {caveat}")
        for name, reason in self.excluded.items():
            lines.append(f"  excluded: {name} -- {reason}")
        return "\n".join(lines)

    @property
    def largest_disagreement(self):
        """``(group, difference)`` for the group the two differ on most."""
        if not self.groups:
            return None, 0.0
        name = max(self.groups,
                   key=lambda key: abs(self.groups[key][1] - self.groups[key][0]))
        first, second = self.groups[name]
        return name, second - first


def _merged_groups(states_a, states_b):
    """Coarsest partition of the DSSP alphabet that is coarser than both.

    Two blocks join when they share a state; transitively, that is the
    connected components of the overlap graph.
    """
    blocks = [set(states) for states in list(states_a) + list(states_b)]
    merged = []
    for block in blocks:
        overlapping = [existing for existing in merged if existing & block]
        combined = set(block)
        for existing in overlapping:
            combined |= existing
            merged.remove(existing)
        merged.append(combined)
    return merged


def _compare(first, second) -> Comparison:
    excluded = {}

    def usable(composition, label):
        keep = {}
        for category, fraction in composition.fractions.items():
            if fraction is None:
                excluded[f"{label}:{category.name}"] = "not estimated"
            elif not category.is_structural:
                excluded[f"{label}:{category.name}"] = (
                    category.note or "no DSSP equivalent")
            else:
                keep[category] = keep.get(category, 0.0) + fraction
        return keep

    kept_a = usable(first, 'first')
    kept_b = usable(second, 'second')

    groups = _merged_groups([c.states for c in kept_a],
                            [c.states for c in kept_b])

    result = {}
    for group in groups:
        total_a = sum(f for c, f in kept_a.items() if c.states & group)
        total_b = sum(f for c, f in kept_b.items() if c.states & group)
        names = sorted({c.name for c in list(kept_a) + list(kept_b)
                        if c.states & group})
        result["+".join(names)] = (total_a, total_b)

    if result:
        differences = np.array([b - a for a, b in result.values()])
        rmsd = float(np.sqrt(np.mean(differences ** 2)))
    else:
        rmsd = 0.0

    return Comparison(groups=result, rmsd=rmsd, excluded=excluded,
                      a=first, b=second,
                      nominal=_by_name(kept_a, kept_b),
                      caveats=_caveats(kept_a, kept_b))


#: Prefixes marking a subdivision that is positional rather than structural,
#: summed into the parent for the nominal comparison.
_VARIANT_PREFIXES = ('regular-', 'distorted-')


def _parent_name(category):
    name = category.name
    for prefix in _VARIANT_PREFIXES:
        if name.startswith(prefix):
            return name[len(prefix):]
    return name.rstrip('s') if name.endswith('s') else name


def _by_name(kept_a, kept_b):
    """Match categories by name, summing positional variants into the parent.

    Approximate by construction: it is comparing labels, not state sets. It
    exists because this comparison gets made whether or not the library offers
    it, and offering it alongside the caveats is more use than withholding it.
    """
    def collapse(kept):
        totals = {}
        for category, fraction in kept.items():
            key = _parent_name(category)
            totals[key] = totals.get(key, 0.0) + fraction
        return totals

    first, second = collapse(kept_a), collapse(kept_b)
    return {name: (first[name], second[name])
            for name in sorted(set(first) & set(second))}


def _caveats(kept_a, kept_b):
    """Name the states that make the two vocabularies disagree.

    The actionable half of a comparison: not "these are incomparable" but
    "the second method files 3-10 and pi helix under unordered, so its helix
    number is lower than the first's for that reason alone".
    """
    notes = []
    for category_a in kept_a:
        for category_b in kept_b:
            if _parent_name(category_a) != _parent_name(category_b):
                continue
            only_a = category_a.states - category_b.states
            only_b = category_b.states - category_a.states
            if only_a:
                notes.append(
                    f"'{category_a.name}' claims "
                    f"{', '.join(sorted(only_a))} which '{category_b.name}' "
                    "does not; the first is higher for that reason alone"
                )
            if only_b:
                notes.append(
                    f"'{category_b.name}' claims "
                    f"{', '.join(sorted(only_b))} which '{category_a.name}' "
                    "does not; the second is higher for that reason alone"
                )
    return sorted(set(notes))


# ---------------------------------------------------------------------------
# FTIR
# ---------------------------------------------------------------------------

def from_ftir(spectrum, method=None, *, bands=AMIDE_I_BANDS,
              positions=None, region=(1480.0, 1720.0),
              model='voigt', derivative_weight=2.0,
              position_tolerance=4.0, **fit_kwargs) -> Composition:
    """
    Estimate secondary structure from an amide I band.

    Parameters
    ----------
    spectrum : Spectrum
        **Already baseline-corrected and water-subtracted.** This function does
        not attempt either: both need judgement -- which reference, what scale
        factor -- and doing them silently inside an estimator would hide the
        step that most affects the answer. What it does do is weight the fit
        towards the second derivative, which greatly reduces the damage done by
        whatever background is left (see ``derivative_weight``).
    method : str
        Required. Currently ``'amide-i-curve-fit'``. Named rather than
        defaulted because which estimator was used is part of the result.
    bands : sequence of (Category, (low, high))
        The assignment table. Defaults to :data:`AMIDE_I_BANDS`; the boundaries
        differ between authors, so they are a parameter rather than a constant.
    positions : array_like, optional
        Starting positions for the components. Found from the second
        derivative when omitted, which is the standard approach: overlapping
        amide I bands have no maxima of their own.
    region : tuple, default (1480, 1720)
        What to **fit**, which is deliberately wider than what is
        **interpreted**. The default spans amide I and amide II together.

        Fitting both bands constrains the baseline far better than fitting
        amide I alone: a sloping or curved residual cannot be absorbed into the
        amide I components when it also has to be consistent with amide II
        eighty wavenumbers away. Only the components falling inside the
        assignment table (``bands``) contribute to the composition; the amide
        II components are there to hold the baseline honest, and are reported
        in ``quality['outside_assignment']`` rather than silently dropped.
    model, derivative_weight, position_tolerance, **fit_kwargs
        Passed to :meth:`spectroscopy.spectra.Spectrum.fit_peaks`.

    Returns
    -------
    Composition

    Raises
    ------
    ValueError
        If ``method`` is missing or unknown, or if no bands could be found.

    Notes
    -----
    ``quality`` carries the fit's R-squared and RMSE, the number of components,
    and the per-component position uncertainties -- read them. A weak shoulder
    between two strong bands is barely determined by the data, and its area is
    a structure fraction.

    Examples
    --------
    >>> composition = from_ftir(spectrum,                    # doctest: +SKIP
    ...                         method='amide-i-curve-fit')
    >>> composition.get('helix')                             # doctest: +SKIP
    0.38
    """
    if method is None:
        raise ValueError(
            "from_ftir needs an explicit method: "
            f"{', '.join(FTIR_METHODS)}. Which estimator was used is part of "
            "the result, so it is not defaulted."
        )
    if method not in FTIR_METHODS:
        raise ValueError(
            f"unknown method {method!r}; available: {', '.join(FTIR_METHODS)}"
        )

    band = spectrum.crop(*region)
    fit = band.fit_peaks(positions, model=model,
                         derivative_weight=derivative_weight,
                         position_tolerance=position_tolerance, **fit_kwargs)

    # Interpretation happens only over the span the assignment table covers.
    # Everything else -- amide II, and anything else inside the fitted window --
    # supported the fit without contributing to the composition.
    assigned_low = min(low for _, (low, _) in bands)
    assigned_high = max(high for _, (_, high) in bands)
    inside = ((fit.position >= assigned_low) & (fit.position < assigned_high))

    # Sum component areas into the assignment ranges. Done here rather than via
    # FitResult.assign so that a category appearing twice in the table -- sheet
    # has a low and a high range -- accumulates instead of overwriting.
    totals = {}
    order = []
    for category, _ in bands:
        if category not in totals:
            totals[category] = 0.0
            order.append(category)

    # Renormalise over the interpreted region: a composition is a share of the
    # amide I band, not of everything that happened to be fitted.
    interpreted_area = float(fit.area[inside].sum())
    outside_area = float(fit.area[~inside].sum())
    if interpreted_area == 0:
        raise ValueError(
            f"no fitted component fell inside the assignment range "
            f"({assigned_low:g} to {assigned_high:g}). Check the region, the "
            "starting positions, and that the x axis is in the units the band "
            "table assumes."
        )

    unassigned = 0.0
    for position, area in zip(fit.position[inside], fit.area[inside]):
        fraction = float(area) / interpreted_area
        for category, (low, high) in bands:
            if low <= position < high:
                totals[category] += fraction
                break
        else:
            unassigned += fraction

    result = {category: totals[category] for category in order}
    if unassigned:
        result[Category('unassigned', frozenset(),
                        note='fitted outside every band in the table')] = unassigned

    quality = {
        'r_squared': fit.r_squared,
        'rmse': fit.rmse,
        'components': len(fit),
        'position_stderr': (None if fit.stderr is None
                            else fit.stderr['position']),
        'model': fit.model,
        'derivative_weight': derivative_weight,
        'fitted_region': tuple(region),
        'interpreted_region': (assigned_low, assigned_high),
        'outside_assignment': int((~inside).sum()),
        'outside_assignment_area': (outside_area / (interpreted_area + outside_area)
                                    if interpreted_area + outside_area else 0.0),
    }
    return Composition(fractions=result, method=method,
                       technique=spectrum.technique or 'FTIR',
                       quality=quality, source=spectrum.name)


# ---------------------------------------------------------------------------
# Circular dichroism -- structure from the shape of the spectrum
# ---------------------------------------------------------------------------

#: The two kinds of CD standard, which have to be told apart because they
#: answer with different things (ADR-0002 section 7.2).
#:
#: ``'basis-spectra'``
#:     The basis is one spectrum per **structural category** -- pure helix,
#:     pure sheet, coil. The fitted coefficients *are* the composition. Simple,
#:     and only ever as good as the basis.
#: ``'reference-proteins'``
#:     The basis is whole **proteins of known structure**, each carrying its
#:     own composition. The unknown is fitted as a combination of proteins and
#:     the same combination is taken of their compositions. This is where the
#:     field is -- SELCON, CDSSTR, CONTIN -- and those methods differ mainly in
#:     how they select and weight the reference set.
CD_METHODS = ('basis-spectra', 'reference-proteins')

def _reference_set(references):
    """
    Check that the standards arrived as a set that knows its own structures.

    Taking spectra and their compositions as two arguments is what ADR-0004
    removed: the lists can fall out of step, and the case that does not raise
    -- same length, wrong order -- moved AqpZ's helix estimate from 0.464 to
    0.307 without complaining.
    """
    from spectroscopy.library import ReferenceSet  # noqa: PLC0415

    if isinstance(references, ReferenceSet):
        return references
    raise TypeError(
        f"references must be a library.ReferenceSet, not "
        f"{type(references).__name__}. Build one with "
        f"ReferenceSet.from_compositions(spectra, compositions), or load one "
        f"with library.load_basis(manifest) or "
        f"library.load_dichroweb_basis(directory) -- each reference then "
        f"carries its own known structure, and there is no second list to "
        f"keep in step."
    )


def _cd_design(spectrum, basis, region):
    """
    Put sample and basis on a common axis, at the **coarser** of the two
    spacings, and crop to ``region``.

    Coarser, not the sample's, because a spectrum cannot carry more
    information about a basis than the basis itself has. Resampling a 1 nm
    reference set onto a 0.1 nm measurement produced a design matrix of 439
    rows whose numerical rank was 46 -- ten times as many equations as there
    was information, which made a hopelessly underdetermined fit look
    well-posed to any check that counts rows.
    """
    x = np.asarray(spectrum.x, dtype=float)
    spacings = [float(np.median(np.abs(np.diff(np.asarray(b.x, dtype=float)))))
                for b in basis]
    coarsest = max(spacings + [float(np.median(np.abs(np.diff(x))))])
    if coarsest > 0:
        low_end, high_end = float(np.min(x)), float(np.max(x))
        x = np.arange(low_end, high_end + coarsest / 2, coarsest)
    if region is not None:
        low, high = sorted(region)
        inside = (x >= low) & (x <= high)
        if inside.sum() < len(basis) + 1:
            raise ValueError(
                f"the region {low:g}-{high:g} nm holds {int(inside.sum())} "
                f"points, which cannot determine {len(basis)} components. "
                f"Widen it, or use fewer references."
            )
        x = x[inside]
    measured = np.asarray(spectrum.resample(x).y, dtype=float)
    design = np.column_stack([np.asarray(reference.resample(x).y, dtype=float)
                              for reference in basis])
    return x, measured, design


def _fit_fractions(design, measured):
    """
    Non-negative coefficients. The fractions are these, normalised.

    Non-negativity is physics: a negative fraction of helix is not a small
    number, it is a wrong model.

    Summing to one is **not** imposed as a constraint, and that is deliberate.
    The measured spectrum has an unknown overall amplitude -- it may be in
    millidegrees at any concentration -- so a constraint forcing the
    coefficients to sum to one fights that amplitude instead of describing the
    shape, and makes the answer depend on how loud the spectrum happens to be.
    Fitting freely and normalising afterwards imposes the same constraint for
    nothing and leaves the fit scale-free, which is what lets an unknown
    concentration still give a composition.

    The recovered sum is kept as ``weight_sum``: with a basis in absolute
    units it is the amplitude, and worth a look.
    """
    from scipy.optimize import nnls  # noqa: PLC0415

    scale = float(np.max(np.abs(measured))) or 1.0
    coefficients, _ = nnls(design / scale, measured / scale)
    return coefficients


def from_cd(spectrum, method=None, *, references=None,
            region=(190.0, 250.0)) -> Composition:
    """
    Estimate secondary structure from the **shape** of a CD spectrum.

    Shape, not size: the fit is scale-free, so the spectrum does not have to
    be in mean residue ellipticity and an unknown concentration does not
    prevent an answer. That is deliberate -- getting to MRE needs a
    concentration, a path length and a residue count, and requiring all three
    before any structure could be estimated would block the common case for a
    reason that does not apply to it.

    Parameters
    ----------
    spectrum : Spectrum
        The measured CD spectrum. Any ellipticity unit.
    method : str
        Required, one of :data:`CD_METHODS`. Named rather than defaulted
        because which kind of standard was used is most of what the answer
        means, and the two are not interchangeable.
    references : ReferenceSet
        The standards, carrying their own known structure -- from
        :func:`~spectroscopy.library.load_basis`,
        :func:`~spectroscopy.library.load_dichroweb_basis` or
        :meth:`~spectroscopy.library.ReferenceSet.from_compositions`.

        Both kinds of standard are one type. A **structural basis** gives each
        spectrum one category, so its known truth is one-hot and the fitted
        coefficients *are* the composition. A **reference-protein set** gives
        each spectrum a whole composition, and the answer is the same mixture
        of those. ``method`` must match which kind was supplied, and is checked
        against it: the two are not interchangeable and the result has to say
        which was used (ADR-0002 section 7.2).
    region : tuple, default (190, 250)
        Wavelength range to fit, nm. The default is the far-UV amide region.
        Below about 190 nm most instruments run out of light and the noise
        rises steeply, which is why the default stops there rather than at
        whatever the file happens to contain.

    Returns
    -------
    Composition

    Notes
    -----
    **No reference data ships with this package**, by decision (ADR-0002
    section 9): the published sets have redistribution terms that have not been
    checked, and inventing a basis would make a decorative answer look like a
    measurement. You supply the basis; the package supplies the arithmetic.

    The arithmetic is tested against synthetic mixtures of a known basis, which
    proves it recovers what it is given. It does not prove that any particular
    basis describes your protein -- that depends entirely on the standards, and
    the ``rmsd`` in :attr:`Composition.quality` is what tells you whether the
    fit was able to reproduce your spectrum at all.
    """
    if method not in CD_METHODS:
        raise ValueError(
            f"method must be one of {CD_METHODS}, got {method!r}. It is "
            f"required because a basis of pure structures and a set of "
            f"reference proteins give different answers from the same "
            f"spectrum, and the result has to say which was used."
        )
    if not references:
        raise ValueError(
            "from_cd needs references: no reference spectra ship with this "
            "package (ADR-0002 section 9). Pass references=... -- a "
            "library.ReferenceSet of pure-structure spectra for "
            "'basis-spectra', or of proteins of known structure for "
            "'reference-proteins'."
        )
    references = _reference_set(references)
    if not references.has_truth:
        references.compositions          # raises, naming what is missing
    if references.is_structural != (method == 'basis-spectra'):
        supplied = ('a structural basis' if references.is_structural
                    else 'reference proteins of known structure')
        wanted = ('basis-spectra' if references.is_structural
                  else 'reference-proteins')
        raise ValueError(
            f"method={method!r} does not match the standards supplied, which "
            f"are {supplied}. Use method={wanted!r}, or pass the other kind of "
            f"set. The two give different answers from the same spectrum, "
            f"which is why the method is named rather than inferred."
        )

    x, measured, design = _cd_design(spectrum, references, region)

    # Rows are not information. With more references than the data can
    # distinguish, least squares still returns *an* answer -- one of infinitely
    # many that fit equally well -- and its residual looks excellent. On a real
    # AqpZ spectrum against all 71 of SP175 this gave a condition number of
    # 2.6e17 and a composition 34 points from the crystal structure, at an rmsd
    # of 0.6%.
    rank = int(np.linalg.matrix_rank(design))
    if rank < design.shape[1]:
        raise ValueError(
            f"{design.shape[1]} references cannot be determined from data of "
            f"rank {rank}: the fit would be one of infinitely many that "
            f"describe the spectrum equally well, and its residual would look "
            f"excellent. Use fewer references -- see "
            f"structure.nearest_references for choosing them by shape -- or "
            f"measure further into the far UV."
        )
    coefficients = _fit_fractions(design, measured)

    total = float(np.sum(coefficients))
    if total <= 0:
        raise ValueError(
            "the fit put no weight on any reference, which means the basis "
            "cannot describe this spectrum at all -- check that both are in "
            "the same wavelength range and the same sign convention."
        )

    # One path for both kinds of standard. A structural basis has a one-hot
    # table, so this reduces to "the coefficients are the composition"; a
    # reference-protein set mixes whole compositions in the same proportions.
    # The third case really was the second with a degenerate table, which is
    # the argument for it not having had its own code path (ADR-0004 section 3).
    fractions = {}
    for composition, weight in zip(references.compositions, coefficients):
        for category, value in composition.fractions.items():
            if value is None:
                continue
            fractions[category] = (fractions.get(category, 0.0)
                                   + (weight / total) * value)

    residual = measured - design @ coefficients
    span = float(np.max(measured) - np.min(measured)) or 1.0
    quality = {
        'rmsd': float(np.sqrt(np.mean(residual ** 2))),
        'rmsd_relative': float(np.sqrt(np.mean(residual ** 2)) / span),
        'n_references': len(references),
        'fitted_points': int(len(x)),
        'region': (float(x.min()), float(x.max())),
        # What the coefficients summed to before being normalised. With a
        # basis in absolute units this is the spectrum's amplitude; with an
        # arbitrary one it is arbitrary too. Either way the fractions below
        # are a projection onto whatever the basis spans, so 'rmsd_relative'
        # is the number to read before believing any of them.
        'weight_sum': total,
        'condition_number': float(np.linalg.cond(design)),
    }
    return Composition(fractions=fractions, method=method,
                       technique=spectrum.technique or 'CD',
                       quality=quality, source=spectrum.name)


#: Chen, Yang & Chau (1974) reference ellipticity for a fully helical chain at
#: 222 nm, deg cm^2 dmol^-1. The chain-length term matters: a helix has two
#: ends that contribute no hydrogen bonds, so a short chain gives a weaker
#: signal per residue than a long one, and ignoring it overestimates helix in
#: small proteins and peptides.
HELIX_REFERENCE_222 = -39500.0
HELIX_CHAIN_CORRECTION = 2.57


def helix_from_theta222(spectrum, residues=None) -> Composition:
    """
    Helix fraction from the mean residue ellipticity at 222 nm.

    **This is not a decomposition and must not be read as one.** It estimates
    one number from one wavelength. The returned :class:`Composition` fills
    ``helix`` and leaves every other category ``None`` -- which means *this
    method cannot estimate this*, and is deliberately not zero, because zero
    would be a claim about sheet content that a single wavelength at 222 nm is
    in no position to make (ADR-0002 section 7.2).

    Use it when the far-UV data does not reach low enough for a shape fit --
    below about 200 nm a detergent-containing or high-salt buffer often
    saturates the detector, and 222 nm survives that when 195 nm does not. It
    is a real answer from a compromised spectrum, not a second-best version of
    :func:`from_cd`.

    Parameters
    ----------
    spectrum : Spectrum
        **In mean residue ellipticity.** Millidegrees are refused: the
        conversion needs a concentration, a path length and a residue count,
        and guessing any of them scales the answer silently. See
        :meth:`~spectroscopy.spectra.Spectrum.to_mean_residue_ellipticity`.
    residues : int, optional
        Residues per chain, for the chain-length correction. Taken from
        ``metadata['n_residues']`` when omitted; without either, the
        correction is skipped and that is recorded in ``quality``.

    Returns
    -------
    Composition
    """
    if spectrum.y_unit != 'deg cm^2 dmol^-1':
        raise ValueError(
            f"helix_from_theta222 needs mean residue ellipticity, and this "
            f"spectrum is in {spectrum.y_unit!r}. The conversion needs the "
            f"sample's concentration, path length and residue count -- see "
            f"Spectrum.to_mean_residue_ellipticity(). It is not applied "
            f"automatically because a guessed path length rescales the helix "
            f"content without changing how the spectrum looks."
        )

    x = np.asarray(spectrum.x, dtype=float)
    index = int(np.argmin(np.abs(x - 222.0)))
    if abs(x[index] - 222.0) > 2.0:
        raise ValueError(
            f"no point within 2 nm of 222 nm; the nearest is {x[index]:g} nm"
        )
    measured = float(np.asarray(spectrum.y, dtype=float)[index])

    residues = residues or spectrum.metadata.get('n_residues')
    if residues:
        reference = HELIX_REFERENCE_222 * (
            1.0 - HELIX_CHAIN_CORRECTION / float(residues))
    else:
        reference = HELIX_REFERENCE_222

    fraction = measured / reference
    helix = Category('helix', frozenset({'G', 'H', 'I'}))
    quality = {
        'theta222': measured,
        'reference_222': reference,
        'chain_length_corrected': bool(residues),
        'n_residues': residues,
        'single_wavelength': True,
    }
    if not 0.0 <= fraction <= 1.0:
        quality['out_of_range'] = True
        warnings.warn(
            f"theta222 gives a helix fraction of {fraction:.2f}, which is "
            f"outside 0-1 and so cannot be right. The usual cause is a wrong "
            f"path length, concentration or residue count in the conversion "
            f"to mean residue ellipticity -- each of them scales this "
            f"linearly.",
            UserWarning, stacklevel=2)

    return Composition(
        fractions={helix: fraction},
        method='theta-222', technique=spectrum.technique or 'CD',
        quality=quality, source=spectrum.name)


def cd_shape_descriptors(spectrum, region=None, edge_tolerance=1.0):
    """
    Amplitude-free descriptors of a far-UV CD spectrum.

    Every quantity here is a **wavelength or a ratio**, so none of them
    changes if the spectrum is multiplied by a constant. That makes them
    independent of concentration, path length, residue count and any pipetting
    error -- the four things that decide an absolute ellipticity and are the
    usual reason two labs disagree about the same protein.

    Over-reliance on amplitude is a standing weakness of UV-CD analysis:
    :func:`helix_from_theta222` is entirely an amplitude measurement, and it
    inherits every error in the three numbers its conversion needs. These
    descriptors, and :func:`from_cd`, which is scale-free by construction, say
    what the *shape* supports on its own.

    Parameters
    ----------
    spectrum : Spectrum
        Far-UV CD, any ellipticity unit. **Crop it to the range the
        photomultiplier could actually see** before calling this: on a JASCO
        the HT channel above about 600 V means the detector is starved and the
        signal there is not a measurement. Nothing here can detect that for
        you -- saturated noise has a shape too.
    region : tuple, optional
        Restrict the analysis, in nm.
    edge_tolerance : float
        How close to the end of the range a minimum may sit before it is
        called an edge artefact rather than a band.

    Returns
    -------
    dict
        ``zero_crossing`` -- where the spectrum crosses zero on the blue side,
        which moves from about 200 nm for helix towards 210 nm and beyond as
        helix is lost. A position, so nothing about amplitude enters it.

        ``minimum`` -- position of the most negative point, and
        ``minimum_at_edge``, which is ``True`` when that point is the end of
        the range rather than a turning point. Then it is not a band position
        at all, and any interpretation of it is an interpretation of where the
        data stopped.

        ``ratio_222_208`` -- around 0.8-0.9 for isolated helices, at or above
        1.0 for interacting ones in a bundle or coiled coil. Raised also by
        absorption flattening in a scattering or membrane sample, so it is
        evidence about helix packing, not proof.

        ``ratio_222_minimum``, and the ``region`` actually used.
    """
    x = np.asarray(spectrum.x, dtype=float)
    y = np.asarray(spectrum.y, dtype=float)
    order = np.argsort(x)
    x, y = x[order], y[order]
    if region is not None:
        low, high = sorted(region)
        inside = (x >= low) & (x <= high)
        x, y = x[inside], y[inside]
    if len(x) < 3:
        raise ValueError("not enough points in the region to describe a shape")

    def at(wavelength):
        index = int(np.argmin(np.abs(x - wavelength)))
        return (y[index] if abs(x[index] - wavelength) <= 2.0 else np.nan)

    changes = np.flatnonzero(np.diff(np.sign(y)) != 0)
    crossing = np.nan
    if len(changes):
        first = changes[0]
        crossing = float(np.interp(0.0, [y[first], y[first + 1]],
                                   [x[first], x[first + 1]])
                         if y[first] < y[first + 1] else
                         np.interp(0.0, [y[first + 1], y[first]],
                                   [x[first + 1], x[first]]))

    lowest = int(np.argmin(y))
    at_edge = bool(x[lowest] <= x[0] + edge_tolerance
                   or x[lowest] >= x[-1] - edge_tolerance)
    if at_edge:
        warnings.warn(
            f"the most negative point is at {x[lowest]:g} nm, the end of the "
            f"range {x[0]:g}-{x[-1]:g} nm, so it is where the data stops "
            f"rather than a band. The real minimum is outside what was "
            f"measured -- or outside what the detector could see.",
            UserWarning, stacklevel=2)

    theta222, theta208 = at(222.0), at(208.0)
    return {
        'zero_crossing': crossing,
        'minimum': float(x[lowest]),
        'minimum_at_edge': at_edge,
        'ratio_222_208': float(theta222 / theta208)
                         if np.isfinite(theta222 * theta208) and theta208
                         else np.nan,
        'ratio_222_minimum': float(theta222 / y[lowest])
                             if np.isfinite(theta222) and y[lowest] else np.nan,
        'region': (float(x[0]), float(x[-1])),
    }


def nearest_references(spectrum, references, count=8, region=(190.0, 240.0)):
    """
    Which reference proteins does this spectrum most *look* like, and what
    are their structures?

    Amplitude-free: both spectra are normalised to unit length before the
    comparison, so this asks only about shape. Nothing is fitted and nothing
    is inverted, which is the point -- there is no underdetermined system to
    have infinitely many answers to.

    **Prefer this over a whole-set fit when the two disagree.** Fitting a
    spectrum as a combination of many references is underdetermined whenever
    the reference set is larger than the data's rank, and averaging over the
    solutions that fit **regresses towards the mean composition of the
    reference set**. Measured on an AqpZ spectrum against SMP180, whose set
    mean is 33% helix: the averaged fit returned 40% helix and 21% sheet, and
    tightening the acceptance from 3% to the best 0.1% of subsets moved it
    only to 45%, where it plateaued. The eight nearest shapes gave a mean of
    **67% helix and no sheet**, against 70% and none from the crystal
    structure. The pull towards the set mean grows with how far the protein
    is from it, which is exactly when an answer matters.

    Parameters
    ----------
    spectrum : Spectrum
        The measured CD spectrum, any ellipticity unit.
    references : ReferenceSet
        Reference proteins. Where they carry their known structures -- as any
        set built by :mod:`spectroscopy.library` does -- the result also
        carries a composition averaged over the neighbours, weighted by
        similarity. A set without them still ranks by shape.
    count : int
        How many neighbours to report.
    region : tuple
        Wavelength range to compare over, nm.

    Returns
    -------
    dict
        ``neighbours`` -- ``(similarity, name, composition)``, most similar
        first. ``composition`` -- the similarity-weighted mean over them, if
        the references carry known structures; and ``spread``, the standard
        deviation across the neighbours, which is the honest uncertainty.
        ``region``.

    Notes
    -----
    A high similarity to a protein of known structure is evidence; it is not
    a measurement, and two different folds can share a far-UV shape. Read the
    neighbour list, not just the average -- if the eight nearest disagree
    wildly about sheet content, so does the answer.
    """
    x = np.asarray(spectrum.x, dtype=float)
    low, high = sorted(region)
    grid = np.arange(max(low, x.min()), min(high, x.max()) + 0.5, 1.0)
    if len(grid) < 10:
        raise ValueError(
            f"only {len(grid)} nm of overlap between the spectrum and "
            f"{low:g}-{high:g} nm; there is no shape to compare"
        )

    def unit(values):
        norm = float(np.linalg.norm(values))
        return np.asarray(values, dtype=float) / (norm or 1.0)

    references = _reference_set(references)
    compositions = references.compositions if references.has_truth else None

    measured = unit(spectrum.resample(grid).y)
    scored = []
    for index, reference in enumerate(references):
        similarity = float(measured @ unit(reference.resample(grid).y))
        scored.append((similarity, reference.name,
                       compositions[index] if compositions is not None else None))
    scored.sort(key=lambda row: -row[0])
    neighbours = scored[:count]

    result = {'neighbours': neighbours, 'region': (float(grid[0]),
                                                   float(grid[-1]))}
    if compositions is not None:
        weights = np.array([max(similarity, 0.0)
                            for similarity, _, _ in neighbours])
        weights = weights / (weights.sum() or 1.0)
        categories = list(neighbours[0][2].fractions)
        table = np.array([[composition.fractions.get(category, 0.0)
                           for category in categories]
                          for _, _, composition in neighbours], dtype=float)
        mean = weights @ table
        result['composition'] = Composition(
            fractions=dict(zip(categories, mean)),
            method='nearest-reference-shapes',
            technique=spectrum.technique or 'CD',
            quality={
                'n_neighbours': len(neighbours),
                'best_similarity': float(neighbours[0][0]),
                'worst_similarity': float(neighbours[-1][0]),
                'spread': dict(zip([c.name for c in categories],
                                   table.std(axis=0).tolist())),
                'single_wavelength': False,
            },
            source=spectrum.name)
        result['spread'] = dict(zip([c.name for c in categories],
                                    table.std(axis=0).tolist()))
    return result
