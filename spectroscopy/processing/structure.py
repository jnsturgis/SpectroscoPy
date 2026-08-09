# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Protein secondary structure, from spectra.

One output type, several inputs. :func:`from_ftir` is the first estimator;
circular dichroism follows on its own branch, and a structure from a PDB file
is a third. They return the same :class:`Composition` against the same
vocabulary, which is what makes two estimates of the same sample comparable.

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
    >>> composition = from_ftir(spectrum, method='amide-i-curve-fit')
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

def _cd_design(spectrum, basis, region):
    """Resample the basis onto the sample's axis and crop both to ``region``."""
    x = np.asarray(spectrum.x, dtype=float)
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


def from_cd(spectrum, method=None, *, basis=None, compositions=None,
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
    basis : sequence of Spectrum
        The standards. For ``'basis-spectra'`` each must carry the
        :class:`Category` it represents in ``metadata['category']``. For
        ``'reference-proteins'`` they are whole proteins and ``compositions``
        supplies their structures.
    compositions : sequence of Composition
        One per reference protein, required for ``'reference-proteins'``.
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
    if not basis:
        raise ValueError(
            "from_cd needs a basis: no reference spectra ship with this "
            "package (ADR-0002 section 9). Pass basis=[...] -- pure-structure "
            "spectra for 'basis-spectra', reference proteins plus their "
            "compositions for 'reference-proteins'."
        )

    basis = list(basis)
    x, measured, design = _cd_design(spectrum, basis, region)
    coefficients = _fit_fractions(design, measured)

    total = float(np.sum(coefficients))
    if total <= 0:
        raise ValueError(
            "the fit put no weight on any reference, which means the basis "
            "cannot describe this spectrum at all -- check that both are in "
            "the same wavelength range and the same sign convention."
        )

    if method == 'basis-spectra':
        fractions = {}
        for reference, weight in zip(basis, coefficients):
            category = reference.metadata.get('category')
            if category is None:
                raise ValueError(
                    f"basis spectrum {reference.name!r} has no "
                    f"metadata['category']; with method='basis-spectra' each "
                    f"basis spectrum must say which structural category it "
                    f"represents, as a structure.Category."
                )
            fractions[category] = fractions.get(category, 0.0) + weight / total
    else:
        if compositions is None or len(compositions) != len(basis):
            raise ValueError(
                f"method='reference-proteins' needs one Composition per "
                f"reference protein; got {len(basis)} references and "
                f"{0 if compositions is None else len(compositions)} "
                f"compositions."
            )
        fractions = {}
        for composition, weight in zip(compositions, coefficients):
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
        'n_references': len(basis),
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
