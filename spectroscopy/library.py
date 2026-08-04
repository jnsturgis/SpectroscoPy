# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Reference spectra and extinction coefficients.

Roadmap section 13.3 ranks "library/reference lookup *and* decomposition"
first among the things the paper needs, on the grounds that the pairing is the
novelty: other packages match a spectrum against a library, and other packages
decompose, but identifying what a decomposition pulled out by matching it
against references is the step normally done by eye.

This module is the reference half. :mod:`spectroscopy.processing.unmix` is the
decomposition half.

What is and is not shipped
--------------------------
**Scalar extinction coefficients are shipped**, because they are published,
citable numbers that everyone doing UV-Vis on protein or nucleic acid already
uses. They live in :data:`COEFFICIENTS` with their source recorded.

**Full epsilon(lambda) reference spectra are not.** Inventing a plausible
absorbance curve for dsDNA would be fabricating reference data, and a
fabricated reference is worse than no reference: it makes an unmixing look
quantitative when it is decorative. Real ones come from measurement, which is
what :func:`from_series` is for -- a concentration series in, an extinction
spectrum with per-wavelength uncertainties out. That is also the honest route
to a house library of the things a particular lab actually works with.

The synthetic references used in the documentation are built in the page that
uses them and are labelled as synthetic there.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ['Reference', 'Library', 'COEFFICIENTS', 'coefficient',
           'concentration_from_absorbance']


@dataclass
class Coefficient:
    """A published extinction coefficient at one wavelength."""

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

    def concentration(self, absorbance, path_length=1.0):
        """Beer-Lambert, solved for concentration: ``c = A / (eps * l)``."""
        return np.asarray(absorbance, dtype=float) / (self.value * path_length)


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
    A named set of reference spectra.

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


def from_series(collection, concentrations, name, *, path_length=1.0,
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
    concentrations : sequence of float
        One per spectrum, in whatever unit ``unit`` is the inverse of.
    path_length : float
        Cuvette path length, cm.

    Returns
    -------
    Reference
        With :attr:`Reference.uncertainty` set per wavelength.
    """
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

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
