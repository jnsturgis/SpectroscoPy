# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Reference spectra of known things, to measure an unknown against.

Two questions use this module. "How much of each of these is in my sample?"
needs reference spectra of the pure components, and is answered by the unmix
module. "What structure does this protein have?" needs reference spectra of
proteins whose structures are already known, and is answered by the circular
dichroism module. Both are the same idea at different scales, so both use the
same container.

    >>> import spectroscopy as spc
    >>> references = spc.datasets.reference_set('sp175')
    >>> len(references)
    71

A ReferenceSet is a collection of spectra that each carry what is known about
them. Because the known property travels with its own spectrum there is no
second list to keep in step -- which matters more than it sounds, since a list
of properties in the wrong order raises nothing at all and simply gives a
wrong answer.

You can build one from spectra and their known compositions, load one from a
folder of files described by a small spreadsheet, or use one of the published
sets that ship with the package.

For ultraviolet work, a Library is the same idea keyed by species name, holding
each reference with its units and its measured uncertainty. If you have
measured a dilution series of something, from_series turns it into a reference
by running Beer-Lambert backwards: absorbance against concentration at every
wavelength, with the slope as the extinction coefficient.

What ships, and what does not
-----------------------------
The published circular dichroism reference sets -- SP175 and SMP180 -- ship
with the package, because their terms allow it. So do the scalar extinction
coefficients everyone doing ultraviolet work on protein or nucleic acid
already uses, with their sources recorded.

Full extinction spectra for ultraviolet components do not ship. Inventing a
plausible absorbance curve for DNA would be making up reference data, and made
up reference data is worse than none: it makes a decorative answer look like a
measurement. Measure your own with from_series, which is also the honest way to
get a library of the things your lab actually works with.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from spectroscopy.collection import SpectrumCollection

__all__ = [
    'ReferenceSet',
    'Reference', 'Library', 'COEFFICIENTS', 'coefficient',
    'concentration_from_absorbance', 'from_series',
    'load_basis', 'MANIFEST_COLUMNS',
    'load_dichroweb_basis', 'DICHROWEB_CATEGORIES',
]


@dataclass
class Coefficient:
    """
    A published extinction coefficient at one wavelength.

    .. warning::

       :attr:`value` is an **extinction coefficient**, so for the nucleic acid
       rules it is the *reciprocal* of the number people quote. "An A260 of
       1.0 is 50 ug/mL" is stored as ``0.02 (ug/mL)^-1 cm^-1``. Concentration
       is ``A / (value * l)`` -- use :meth:`concentration` or
       :func:`concentration_from_absorbance` and let them do it. Multiplying
       instead of dividing gives an answer wrong by the square of the
       coefficient, and for dsDNA that is a factor of 2500.

       :attr:`quoted_as` gives the familiar form, for reading and for
       checking against a textbook.
    """

    #: Wavelength, nm.
    wavelength: float
    #: Value in :attr:`unit`.
    value: float
    #: ``'M^-1 cm^-1'``, or a mass form such as ``'(mg/mL)^-1 cm^-1'``.
    unit: str
    #: Where the number comes from. Never leave this empty.
    source: str
    #: What it applies to, in words -- assumptions matter more than the number.
    note: str = ''

    @property
    def quoted_as(self):
        """
        ``1 / value`` -- the concentration giving an absorbance of 1.0 in
        1 cm, which is how these are stated in the literature and on the
        side of a kit box. 50.0 for dsDNA, where :attr:`value` is 0.02.
        """
        return 1.0 / self.value

    def concentration(self, absorbance, path_length=1.0):
        """Beer-Lambert, solved for concentration: ``c = A / (eps * l)``."""
        return np.asarray(absorbance, dtype=float) / (self.value * path_length)

    def __str__(self):
        concentration_unit = self.unit.split('^-1')[0].strip('()')
        return (f"{self.value:g} {self.unit} at {self.wavelength:g} nm "
                f"(an absorbance of 1.0 in 1 cm is "
                f"{self.quoted_as:g} {concentration_unit})")


#: Published scalar coefficients, keyed by ``(species, wavelength)``.
#:
#: These are the conventional numbers, not measurements of anyone's particular
#: sample, and the ``note`` on each says what it assumes. The nucleic acid
#: figures are the standard "an A260 of 1.0 corresponds to N ug/mL" rules,
#: expressed here as ``(ug/mL)^-1 cm^-1`` so Beer-Lambert applies directly.
COEFFICIENTS = {
    ('dsDNA', 260): Coefficient(
        260, 1.0 / 50.0, '(ug/mL)^-1 cm^-1',
        'Conventional value; an A260 of 1.0 is 50 ug/mL double-stranded DNA.',
        'Assumes average base composition. Sequence dependence is real: '
        'GC-rich duplexes absorb less per mass through hypochromicity.'),
    ('ssDNA', 260): Coefficient(
        260, 1.0 / 33.0, '(ug/mL)^-1 cm^-1',
        'Conventional value; an A260 of 1.0 is 33 ug/mL single-stranded DNA.',
        'Single strands absorb more per base than duplexes -- the difference '
        'between this and dsDNA is the hyperchromic effect, not an error.'),
    ('RNA', 260): Coefficient(
        260, 1.0 / 40.0, '(ug/mL)^-1 cm^-1',
        'Conventional value; an A260 of 1.0 is 40 ug/mL RNA.', ''),
    ('protein', 280): Coefficient(
        280, 1.0 / 1000.0, '(ug/mL)^-1 cm^-1',
        'Rule of thumb; an A280 of 1.0 is approximately 1 mg/mL protein.',
        'Only a rule of thumb. A280 comes from tryptophan, tyrosine and '
        'cystines, so it varies several-fold between proteins. Use a '
        'sequence-derived coefficient whenever the sequence is known.'),
}

#: Per-residue contributions to epsilon at 280 nm, for computing a protein's
#: coefficient from its sequence. This is the right way to get A280 for a
#: known protein, and it is far better than the 1 mg/mL rule of thumb above.
#:
#: Source: Pace, Vajdos, Fee, Grimsley & Gray (1995), "How to measure and
#: predict the molar absorption coefficient of a protein", Protein Science
#: 4:2411-2423, following Gill & von Hippel (1989).
RESIDUE_EPSILON_280 = {'W': 5500.0, 'Y': 1490.0, 'cystine': 125.0}


def coefficient(species, wavelength=None):
    """
    Look up a published coefficient.

    Raises rather than guessing when the species is unknown, and lists what is
    available -- a wrong extinction coefficient is a wrong concentration, and
    it fails silently.
    """
    if wavelength is None:
        matches = [key for key in COEFFICIENTS if key[0] == species]
        if len(matches) == 1:
            return COEFFICIENTS[matches[0]]
        if not matches:
            known = sorted({key[0] for key in COEFFICIENTS})
            raise KeyError(
                f"no coefficient for {species!r}; known species are "
                f"{', '.join(known)}. Sequence-derived coefficients are "
                f"better than any of them where the sequence is known -- see "
                f"protein_epsilon_280()."
            )
        raise KeyError(
            f"{species!r} has coefficients at several wavelengths "
            f"({', '.join(str(key[1]) for key in sorted(matches))}); say which"
        )
    try:
        return COEFFICIENTS[(species, wavelength)]
    except KeyError:
        available = [key for key in COEFFICIENTS if key[0] == species]
        if available:
            where = ', '.join(str(key[1]) for key in sorted(available))
            detail = f"{species!r} is only tabulated at {where} nm"
        else:
            known = ', '.join(sorted({key[0] for key in COEFFICIENTS}))
            detail = f"no coefficients for {species!r}; known species are {known}"
        raise KeyError(
            f"no coefficient for {species!r} at {wavelength} nm -- {detail}"
        ) from None


def protein_epsilon_280(tryptophan, tyrosine, cystine=0):
    """
    Molar epsilon at 280 nm from residue counts, in M^-1 cm^-1.

    ``eps = 5500*nW + 1490*nY + 125*n_cystine`` -- Pace et al. (1995). Count
    **cystines** (disulfide bridges), not cysteines: free thiols contribute
    essentially nothing at 280 nm.

    Worth preferring over the 1 mg/mL rule of thumb by a wide margin. A protein
    with no tryptophan absorbs several-fold less at 280 than the rule assumes,
    and the resulting concentration is wrong by that factor with nothing to
    show for it.
    """
    return (RESIDUE_EPSILON_280['W'] * tryptophan
            + RESIDUE_EPSILON_280['Y'] * tyrosine
            + RESIDUE_EPSILON_280['cystine'] * cystine)


def concentration_from_absorbance(absorbance, species, wavelength=None,
                                  path_length=1.0):
    """
    Concentration from a single absorbance reading, via a published
    coefficient.

    Convenience over :meth:`Coefficient.concentration`; the units are those of
    the coefficient, which :func:`coefficient` will tell you.
    """
    return coefficient(species, wavelength).concentration(absorbance,
                                                          path_length)


class ReferenceSet(SpectrumCollection):
    """
    Spectra of samples whose structure is already known by other means.

    A reference set **is** a :class:`~spectroscopy.collection.SpectrumCollection`
    -- not a wrapper around one -- so it arrives with ``to_matrix``, ``select``,
    ``resample``, ``crop``, indexing and iteration already working, and stays a
    ``ReferenceSet`` through all of them (ADR-0004).

    What it adds is *gathered views* of each reference's known truth:
    :attr:`compositions` and :attr:`categories`. The truth is stored in each
    spectrum's own ``metadata``, and these read it back. That is the whole
    point of the class. The arrangement it replaces was two parallel lists,
    ``(spectra, compositions)``, passed around together:

        compositions one short : ValueError from numpy about core dimensions
        compositions reversed  : no error. helix 0.464 -> 0.307

    The second is the dangerous one -- same length, wrong pairing, a plausible
    answer, no complaint. With one list there is nothing to reverse.

    Set-level facts -- ``source``, ``citation``, ``licence``, ``accession``,
    ``unit`` -- live in :attr:`~spectroscopy.collection.SpectrumCollection.info`
    and survive slicing, so a five-protein subset of SP175 still carries
    SP175's citation condition.

    Two kinds of set, one type
    --------------------------
    A **reference-protein set** (SP175, SMP180) gives each spectrum a whole
    composition. A **structural basis** gives each spectrum one category, which
    is the same thing with a degenerate table -- all of one structure. Both
    read back through :attr:`compositions`, which is why they no longer need
    two code paths.

    A UV-Vis unmixing set has no known truth at all: a spectrum in epsilon
    units is its own property and the fitted coefficient is the answer. Such a
    set is still a ``ReferenceSet``; the table is simply empty.
    """

    # -- construction ------------------------------------------------------

    @classmethod
    def from_compositions(cls, spectra, compositions, *, name=None, info=None,
                          known_from=None):
        """
        Build one from spectra and their known compositions, matched by
        position.

        **This is the only place the positional matching happens**, and it is
        deliberately a single, named, checkable step rather than an invariant
        every call site has to maintain. After this, there is one list.

        The spectra are copied, not modified: the truth is written onto the
        copies, so a caller's own spectra do not silently acquire metadata.
        """
        spectra = list(spectra)
        compositions = list(compositions)
        if len(spectra) != len(compositions):
            raise ValueError(
                f"got {len(spectra)} spectra and {len(compositions)} "
                f"compositions. They are matched by position, so the counts "
                f"must agree -- and note that agreeing counts do not make the "
                f"pairing right, which is why this is the only place it is done."
            )

        declared = _declared_categories(compositions)
        labelled = []
        for spectrum, composition in zip(spectra, compositions):
            copied = spectrum._derive()      # pylint: disable=protected-access
            copied.metadata['composition'] = {
                category.name: (None if fraction is None else float(fraction))
                for category, fraction in composition.fractions.items()
            }
            source = known_from or composition.method
            if source:
                copied.metadata['known_from'] = source
            labelled.append(copied)

        info = dict(info or {})
        info.setdefault('categories', declared)
        return cls(labelled, name=name, info=info)

    # -- gathered views ----------------------------------------------------

    @property
    def categories(self):
        """
        The structural categories this set's known truth is declared against.

        Held once, on the set, because it is a fact about the set: the stored
        form of a composition is ``{name: fraction}`` and of a category is its
        name, since ``.spy`` serialises metadata as JSON and turns a
        ``Category`` into a bare string with its DSSP states gone. These are
        what rebuild the objects (ADR-0004 section 2.5).
        """
        declared = self.info.get('categories')
        return list(declared) if declared else []

    @property
    def has_truth(self) -> bool:
        """True when every spectrum carries a composition or a category."""
        return bool(self) and all(
            'composition' in s.metadata or 'category' in s.metadata
            for s in self)

    @property
    def is_structural(self) -> bool:
        """
        True for a basis of pure structures -- each spectrum *is* one category.

        The distinction matters to the reader of a result and not to the
        arithmetic: :attr:`compositions` returns a one-hot composition either
        way, so nothing downstream has to branch on it.
        """
        return bool(self) and all('category' in s.metadata for s in self)

    @property
    def compositions(self):
        """
        Each reference's known composition, in order, as
        :class:`~spectroscopy.processing.structure.Composition` objects.

        Computed on every access rather than cached, which keeps it impossible
        for the view to be stale. The sets in use are at most a few hundred
        spectra; if that changed, caching with invalidation would be worth what
        it costs, and not before.
        """
        from spectroscopy.processing.structure import (  # noqa: PLC0415
            Category,
            Composition,
        )

        declared = {category.name: category for category in self.categories}
        missing = [s.name for s in self
                   if 'composition' not in s.metadata
                   and 'category' not in s.metadata]
        if missing:
            raise ValueError(
                f"{len(missing)} of {len(self)} references have no known "
                f"structure, so this set cannot say what its spectra are of: "
                f"{', '.join(str(n) for n in missing[:5])}"
                f"{' ...' if len(missing) > 5 else ''}. Build the set with "
                f"ReferenceSet.from_compositions(), or load it with "
                f"library.load_basis()."
            )

        result = []
        for spectrum in self:
            stored = spectrum.metadata.get('composition')
            if stored is None:
                label = str(spectrum.metadata['category'])
                stored = {name: (1.0 if name == label else 0.0)
                          for name in declared} or {label: 1.0}
            unknown = sorted(set(stored) - set(declared))
            if unknown:
                raise ValueError(
                    f"{spectrum.name} declares structure in categories "
                    f"{unknown}, which this set does not declare. Its "
                    f"categories are {sorted(declared) or 'undeclared'} -- "
                    f"set collection.info['categories'] to the Category "
                    f"objects the fractions mean, since a bare name does not "
                    f"say which DSSP states it covers."
                )
            result.append(Composition(
                fractions={declared.get(name, Category(name)): fraction
                           for name, fraction in stored.items()},
                method=spectrum.metadata.get(
                    'known_from', 'supplied with the reference set'),
                technique='reference',
                source=spectrum.name,
            ))
        return result

    def __repr__(self) -> str:
        label = f" {self.name!r}" if self.name else ""
        kind = ('structural basis' if self.is_structural
                else 'reference proteins' if self.has_truth
                else 'no known structure')
        where = self.info.get('source')
        return (f"<ReferenceSet{label}: {len(self)} references, {kind}"
                + (f", from {where}" if where else "") + ">")


def _declared_categories(compositions):
    """
    The category objects a set of compositions agree on.

    Two ``Category`` objects with the same name and different DSSP states are
    two different claims about what "helix" means, and averaging over them
    would be silently answering a question nobody asked.
    """
    declared: dict = {}
    for composition in compositions:
        for category in composition.fractions:
            seen = declared.setdefault(category.name, category)
            if seen.states != category.states:
                raise ValueError(
                    f"two references disagree about what {category.name!r} "
                    f"covers: DSSP states {sorted(seen.states)} and "
                    f"{sorted(category.states)}. A reference set has one "
                    f"vocabulary or it is two reference sets."
                )
    return list(declared.values())


@dataclass
class Reference:
    """
    One entry in a :class:`Library`: a spectrum that stands for a species.

    Attributes
    ----------
    name : str
        How the component is reported in an unmixing result.
    spectrum : Spectrum
        The reference itself. Ideally an **extinction spectrum**, epsilon
        against wavelength, in which case unmixing returns concentrations
        directly. A reference in arbitrary units still works and still
        separates the components -- the coefficients are then relative, and
        :attr:`is_absolute` says so.
    unit : str
        Units of the y axis, e.g. ``'M^-1 cm^-1'``. Empty when relative.
    source : str
        Where it came from: a citation, or the measurement that produced it.
    uncertainty : ndarray, optional
        Per-wavelength standard error, as :func:`from_series` produces.
    """

    name: str
    spectrum: object
    unit: str = ''
    source: str = ''
    uncertainty: np.ndarray | None = None
    metadata: dict = field(default_factory=dict)

    @property
    def is_absolute(self) -> bool:
        """True when the y axis is a real extinction coefficient."""
        return bool(self.unit)

    def __repr__(self) -> str:
        kind = self.unit if self.unit else 'relative'
        return f"<Reference {self.name!r} ({kind})>"


class Library:
    """
    A named set of reference spectra, keyed by species name, for unmixing.

    The UV-Vis specialisation of :class:`ReferenceSet`, and the one case with
    no known-truth table at all: a spectrum in epsilon units *is* its own
    property, and the fitted coefficient is the answer. That is why it holds
    :class:`Reference` entries -- name, unit, uncertainty -- rather than bare
    spectra, and why it is keyed rather than ordered.

    .. note::

       ADR-0004 records this as a subclass of ``ReferenceSet``. It is not one
       yet, deliberately: ``Library`` iterates over ``Reference`` objects and a
       collection iterates over ``Spectrum``, so inheriting would change what
       ``for reference in library`` yields -- and that is the loop inside
       :func:`~spectroscopy.processing.unmix.unmix`, whose signature freezes at
       1.0. The unification is worth doing on its own, with the UV-Vis tests
       watching, and not as a side effect of the CD work.

    Deliberately thin. The useful operations are selecting a subset and
    handing it to :func:`spectroscopy.processing.unmix.unmix`; anything
    cleverer belongs where the science is.

        >>> library = Library([water, dna], name='house')      # doctest: +SKIP
        >>> library.select('dna')                              # doctest: +SKIP
    """

    def __init__(self, references=(), name=''):
        self.name = name
        self._entries = {}
        for reference in references:
            self.add(reference)

    def add(self, reference):
        """Add a reference. Replacing an existing name is refused."""
        if reference.name in self._entries:
            raise KeyError(
                f"{reference.name!r} is already in this library; rename one of "
                f"them rather than shadowing, or drop it first"
            )
        self._entries[reference.name] = reference
        return self

    def __getitem__(self, name):
        try:
            return self._entries[name]
        except KeyError:
            raise KeyError(
                f"no reference named {name!r}; this library has "
                f"{', '.join(sorted(self._entries)) or 'nothing in it'}"
            ) from None

    def __len__(self):
        return len(self._entries)

    def __iter__(self):
        return iter(self._entries.values())

    def __contains__(self, name):
        return name in self._entries

    @property
    def names(self):
        return tuple(self._entries)

    def select(self, *names):
        """A new Library holding only the named references, in that order."""
        return Library([self[name] for name in names], name=self.name)

    def __repr__(self):
        where = f" {self.name!r}" if self.name else ""
        return f"<Library{where}: {len(self)} references: {', '.join(self.names)}>"

    def on(self, x):
        """
        Every reference resampled onto ``x``, as an ``(n_references, n_x)``
        matrix.

        Unmixing needs the references on the sample's own wavelength grid, and
        a reference measured on a different instrument never is. Resampling
        the references rather than the sample keeps the data being explained
        untouched.
        """
        x = np.asarray(x, dtype=float)
        return np.vstack([reference.spectrum.resample(x).y
                          for reference in self])


def from_series(collection, concentrations=None, name=None, *, path_length=1.0,
                unit='M^-1 cm^-1', source=''):
    """
    Build a :class:`Reference` from spectra of known concentration.

    This is Beer-Lambert used the other way round: instead of a concentration
    from an absorbance and a known epsilon, an **epsilon spectrum** from a set
    of absorbances at known concentrations. At each wavelength the slope of
    absorbance against ``concentration * path_length`` is the coefficient, and
    its standard error is the uncertainty on it.

    The fit is through the origin, which is the physics -- zero concentration
    absorbs nothing. A systematic offset therefore shows up as curvature in
    the residual rather than being absorbed into an intercept, which is what
    you want: it usually means scattering or a baseline that was not removed.

    Parameters
    ----------
    collection : SpectrumCollection
        Spectra of the same species at different known concentrations, on a
        common wavelength axis.
    concentrations : sequence of float, optional
        One per spectrum, in whatever unit ``unit`` is the inverse of.
        Defaults to the collection's own parameter, so a series loaded with
        ``from_files(parameter_from=...)`` or labelled with
        ``with_parameters(...)`` needs no second list -- and cannot fall out of
        step with one.
    name : str, optional
        What the reference is called. Defaults to the collection's name.
    path_length : float
        Cuvette path length, cm.

    Returns
    -------
    Reference
        With :attr:`Reference.uncertainty` set per wavelength.
    """
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

    if concentrations is None:
        concentrations = getattr(collection, 'parameters', None)
        if concentrations is None or np.isnan(concentrations).any():
            raise ValueError(
                "no concentrations given and the collection does not carry "
                "them. Either pass them, or attach them when loading with "
                "SpectrumCollection.from_files(parameter_from=...) or "
                "collection.with_parameters([...])."
            )
    if name is None:
        name = getattr(collection, 'name', None)
        if not name:
            raise ValueError(
                "no name given and the collection has none; a reference has "
                "to be called something to be selected from a Library later"
            )

    concentrations = np.asarray(concentrations, dtype=float)
    if len(concentrations) != len(collection):
        raise ValueError(
            f"{len(collection)} spectra but {len(concentrations)} "
            f"concentrations; they must correspond one to one"
        )
    if len(collection) < 2:
        raise ValueError(
            "need at least two concentrations to fit a slope; one spectrum "
            "gives a coefficient with no uncertainty and no way to notice "
            "that the response is not linear"
        )

    x = np.asarray(collection[0].x, dtype=float)
    absorbance = np.vstack([spectrum.resample(x).y for spectrum in collection])

    # Least squares through the origin, per wavelength: eps = sum(cA)/sum(c^2).
    load = concentrations * path_length
    denominator = float(np.sum(load ** 2))
    if denominator == 0:
        raise ValueError("all concentrations are zero; nothing to fit")
    epsilon = (load @ absorbance) / denominator

    residual = absorbance - np.outer(load, epsilon)
    if len(collection) > 1:
        variance = np.sum(residual ** 2, axis=0) / (len(collection) - 1)
        uncertainty = np.sqrt(variance / denominator)
    else:                                        # unreachable, kept explicit
        uncertainty = np.zeros_like(epsilon)

    spectrum = Spectrum(x, epsilon,
                        x_quantity=collection[0].x_quantity,
                        x_unit=collection[0].x_unit,
                        y_quantity='Extinction coefficient', y_unit=unit)
    spectrum.name = name
    return Reference(name=name, spectrum=spectrum, unit=unit,
                     source=source or f"fitted from {len(collection)} spectra",
                     uncertainty=uncertainty,
                     metadata={'path_length': path_length,
                               'concentrations': concentrations.tolist()})


# ---------------------------------------------------------------------------
# Loading a basis somebody else measured
# ---------------------------------------------------------------------------

#: Columns a basis manifest may carry. ``file`` and one of ``category`` or the
#: fraction columns are required; the rest are provenance and are optional but
#: strongly wanted -- a reference whose origin is not recorded cannot be cited,
#: and most published sets require citation as a condition of use.
MANIFEST_COLUMNS = ('file', 'name', 'category', 'known_from', 'source',
                    'citation', 'accession')


def load_basis(manifest, directory=None, file_type=None, **read_kwargs):
    """
    Build a CD basis from files you have, described by a small manifest.

    **No reference set ships with this package and none can**: the published
    sets are distributed through the PCDDB, whose terms grant access but never
    grant redistribution, and the SSCalcPy packaging of them is
    NonCommercial. So the arrangement is the one already used for the Galactic
    ``.spc`` samples -- you obtain the data, under whatever terms its provider
    sets, and this loads it. Nothing is redistributed by SpectroscoPy.

    That also makes this deliberately **not PCDDB-specific**. It reads whatever
    :func:`spectroscopy.read` reads, so a basis measured in your own lab, one
    exported from a supplier, and one downloaded from a public bank all load
    the same way.

    Parameters
    ----------
    manifest : str or Path
        A CSV with a header row. One line per reference spectrum:

        ``file``
            Path to the spectrum, relative to ``directory``.
        ``category``
            For a **structural basis**: which category this spectrum is of --
            ``helix``, ``sheet``, ``turn``, ``other``. Matched by name against
            :data:`~spectroscopy.processing.structure.DSSP_STATES`-backed
            categories.
        ``helix``, ``sheet``, ``turn``, ``other`` (any subset)
            For a **reference-protein set**: this protein's known composition,
            as fractions. Use these *instead of* ``category``.
        ``name``, ``source``, ``citation``, ``accession``
            Provenance. Carried into each reference and into the returned
            library, so a result can say where its basis came from.

    directory : str or Path, optional
        Where the files are. Defaults to the manifest's own directory, which
        is the usual arrangement.

        ``known_from``
            How the structure was determined -- ``'DSSP on 1RC2'``. Optional,
            and worth filling in: a composition taken from an earlier CD fit
            makes any set built on it circular, and nothing else records that.

    Returns
    -------
    ReferenceSet
        Either kind, as one type. A structural basis has
        :attr:`~ReferenceSet.is_structural` true and its
        :attr:`~ReferenceSet.compositions` read back as one-hot; a
        reference-protein set carries whole compositions. Both go into
        ``cd.estimate(spectrum, method, references)`` unchanged.

    Notes
    -----
    Nothing here checks that a basis is any *good*. That is what
    ``Composition.quality['rmsd_relative']`` is for, and what a protein of
    known structure is for.
    """
    import csv  # noqa: PLC0415
    from pathlib import Path  # noqa: PLC0415

    from spectroscopy import read as _read  # noqa: PLC0415
    from spectroscopy.processing.structure import (  # noqa: PLC0415
        DSSP_STATES,
        Category,
        Composition,
    )

    #: name -> the DSSP states it claims. Kept here rather than imported so a
    #: manifest may spell a category without importing anything.
    known = {
        'helix': frozenset({'G', 'H', 'I'}),
        'sheet': frozenset({'E', 'B'}),
        'turn': frozenset({'T'}),
        'other': frozenset({'S', '-'}),
    }
    assert set().union(*known.values()) <= set(DSSP_STATES)

    manifest = Path(manifest)
    directory = Path(directory) if directory is not None else manifest.parent
    with manifest.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{manifest} has a header but no reference rows")

    fields = {name.strip().lower() for name in rows[0]}
    if 'file' not in fields:
        raise ValueError(
            f"{manifest} needs a 'file' column; it has {sorted(fields)}"
        )
    fraction_columns = sorted(fields & set(known))
    structural = 'category' in fields
    if structural == bool(fraction_columns):
        raise ValueError(
            f"{manifest} must have either a 'category' column (a structural "
            f"basis: each spectrum is one kind of structure) or fraction "
            f"columns {sorted(known)} (reference proteins: each spectrum is a "
            f"whole protein of known composition) -- not both and not "
            f"neither. It has "
            f"{'both' if structural else 'neither'}."
        )

    spectra, compositions = [], []
    set_level: dict = {}
    for number, row in enumerate(rows, start=2):
        row = {key.strip().lower(): (value or '').strip()
               for key, value in row.items() if key}
        path = directory / row['file']
        if not path.is_file():
            raise FileNotFoundError(
                f"{manifest} line {number} refers to {path}, which is not "
                f"there. Reference data is not shipped with this package -- "
                f"see the docstring for why -- so the files have to be "
                f"fetched before the manifest can be used."
            )
        spectrum = _read(path, file_type, **read_kwargs)
        spectrum.name = row.get('name') or path.stem
        for key in ('source', 'citation', 'accession'):
            if row.get(key):
                spectrum.metadata[f'reference_{key}'] = row[key]
                # A value every row repeats is a fact about the set, and a
                # citation condition needs somewhere to live as one. Where the
                # rows differ it stays per-row and the set says nothing.
                set_level[key] = (row[key] if set_level.get(key, row[key])
                                  == row[key] else None)

        if structural:
            label = row['category']
            if label not in known:
                raise ValueError(
                    f"{manifest} line {number}: category {label!r} is not one "
                    f"of {sorted(known)}"
                )
            spectrum.metadata['category'] = label
            spectrum.metadata['known_from'] = (
                row.get('known_from') or f'{manifest.name}, line {number}')
        else:
            fractions = {}
            for label in fraction_columns:
                if row.get(label):
                    fractions[Category(label, known[label])] = float(row[label])
            total = sum(fractions.values())
            if not 0.9 <= total <= 1.1:
                raise ValueError(
                    f"{manifest} line {number}: the fractions sum to "
                    f"{total:.3f}, which is not a composition. They should be "
                    f"fractions of 1, not percentages."
                )
            compositions.append(Composition(
                fractions=fractions,
                method=(row.get('known_from')
                        or f'{manifest.name}, line {number}'),
                technique='reference', source=spectrum.name))
        spectra.append(spectrum)

    info = {key: value for key, value in set_level.items() if value}
    if structural:
        # Only the categories this basis actually contains. A basis with no
        # turn spectrum cannot estimate turn, and declaring the category would
        # make every fit report a confident 0.0 for it -- a claim the method
        # never made, which ADR-0002 section 7.2 distinguishes from None.
        info['categories'] = [
            Category(label, known[label])
            for label in dict.fromkeys(s.metadata['category'] for s in spectra)
        ]
        return ReferenceSet(spectra, name=manifest.stem, info=info)
    return ReferenceSet.from_compositions(spectra, compositions,
                                          name=manifest.stem, info=info)


#: How the DichroWebGit datasets spell their structural categories, mapped
#: onto the DSSP states this project defines everything against (ADR-0002).
#: "Disorder" is their name for what the FTIR side calls "other": DSSP's bend
#: and none-of-the-above.
DICHROWEB_CATEGORIES = {
    'helix': frozenset({'G', 'H', 'I'}),
    'sheet': frozenset({'E', 'B'}),
    'turn': frozenset({'T'}),
    'disorder': frozenset({'S', '-'}),
}

#: The DichroWebGit datasets all start at 240 nm and step down 1 nm per row.
DICHROWEB_FIRST_NM = 240.0
DICHROWEB_STEP_NM = 1.0


def load_dichroweb_basis(directory, first_nm=DICHROWEB_FIRST_NM,
                         step_nm=DICHROWEB_STEP_NM):
    """
    Load an SP175 / SMP180 / IDP175 reference set in DichroWebGit's layout.

    These are the published CD reference sets, and **they are redistributable**
    -- unlike the same data taken from the PCDDB website, whose terms grant
    access but never grant redistribution. The ``pcddb`` organisation
    publishes them on GitHub as part of DichroWebGit under the **MIT licence**,
    which does permit copying and redistribution provided the copyright notice
    travels with them. Keep ``LICENSE`` beside the data.

    Four files per dataset, as the DichroWebGit README describes them:

    ``A.txt``
        CD data in columns, one column per protein, starting at 240 nm and
        stepping down. Delta epsilon per residue.
    ``F.txt``
        Secondary structure fractions, one row per category, columns in the
        same protein order.
    ``lbl1.txt``
        The category names, one per line, in ``F.txt``'s row order.
    ``lbl2.txt``
        One tab-separated row of protein labels, in column order.

    Returns
    -------
    ReferenceSet
        Ready for ``cd.estimate(spectrum, method, references)``, carrying the
        MIT licence and the source in its ``info`` so a result can say where
        its basis came from.

    Notes
    -----
    The wavelength axis is **not in the files** -- it is implied by the first
    wavelength and the step, which is why both are parameters here rather than
    assumptions buried in the code. If a future dataset starts somewhere else,
    a silently wrong axis would shift every band and still fit something.
    """
    from pathlib import Path  # noqa: PLC0415

    from spectroscopy.processing.structure import (  # noqa: PLC0415
        Category,
        Composition,
    )
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

    directory = Path(directory)
    missing = [name for name in ('A.txt', 'F.txt', 'lbl1.txt', 'lbl2.txt')
               if not (directory / name).is_file()]
    if missing:
        raise FileNotFoundError(
            f"{directory} is not a DichroWebGit dataset: {', '.join(missing)} "
            f"missing. Expected the four-file layout A/F/lbl1/lbl2."
        )

    absorbance = np.loadtxt(directory / 'A.txt')
    fractions = np.atleast_2d(np.loadtxt(directory / 'F.txt'))
    categories = [line.strip() for line in
                  (directory / 'lbl1.txt').read_text().splitlines()
                  if line.strip()]
    names = [label.strip() for label in
             (directory / 'lbl2.txt').read_text().split('\t') if label.strip()]

    if absorbance.shape[1] != fractions.shape[1] or len(names) != absorbance.shape[1]:
        raise ValueError(
            f"{directory}: A.txt has {absorbance.shape[1]} proteins, F.txt "
            f"has {fractions.shape[1]} and lbl2.txt names {len(names)}. They "
            f"are matched by position, so all three must agree."
        )
    if len(categories) != fractions.shape[0]:
        raise ValueError(
            f"{directory}: lbl1.txt names {len(categories)} categories but "
            f"F.txt has {fractions.shape[0]} rows"
        )

    unknown = [name for name in categories
               if name.lower() not in DICHROWEB_CATEGORIES]
    if unknown:
        raise ValueError(
            f"{directory}: unrecognised structural categories {unknown}; "
            f"known are {sorted(DICHROWEB_CATEGORIES)}"
        )
    resolved = [Category(name.lower(), DICHROWEB_CATEGORIES[name.lower()])
                for name in categories]

    x = first_nm - step_nm * np.arange(absorbance.shape[0])
    order = np.argsort(x)
    spectra, compositions = [], []
    for column, name in enumerate(names):
        spectrum = Spectrum(x[order], absorbance[order, column],
                            technique='CD', name=name)
        spectrum.y_quantity = 'Delta epsilon'
        spectrum.y_unit = 'delta epsilon'
        spectra.append(spectrum)
        compositions.append(Composition(
            fractions={category: float(fractions[row, column])
                       for row, category in enumerate(resolved)},
            method='supplied with the reference set',
            technique='reference', source=name))

    return ReferenceSet.from_compositions(
        spectra, compositions, name=directory.name,
        known_from='supplied with the reference set',
        info={
            'source': f'DichroWebGit {directory.name}',
            'licence': 'MIT (c) 2023 Andy Miles, github.com/pcddb/DichroWebGit',
            'unit': 'delta epsilon',
            'categories': resolved,
        })
