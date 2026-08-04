# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Separating a spectrum into known components.

The supervised counterpart of :mod:`spectroscopy.processing.multivariate`.
That module asks "how many things are varying here and what do they look
like"; this one asks "given that I know what these components look like, how
much of each is present". Both are useful and they answer different questions
-- a blind decomposition finds components nobody can name, and this names
components nobody has to find.

Roadmap section 13.3 puts the pairing first among the things the paper needs.

Why non-negative least squares
------------------------------
A concentration cannot be negative, and an ordinary least-squares fit will
happily return one. When it does, the number is not a small negative
concentration -- it is the fit compensating for a component that is missing
from the reference set, or for a background that was not removed. Constraining
the coefficients to be non-negative does not hide that; it pushes the problem
into the **residual**, where it is visible. See :attr:`UnmixResult.residual`.

The A260/A280 case
------------------
Splitting a sample into nucleic acid and protein is the two-component case of
exactly this, and :func:`nucleic_acid_and_protein` is a thin wrapper that says
so. The ratio itself is in :func:`absorbance_ratio`, with what it is and is
not good for written down beside it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import nnls

__all__ = ['unmix', 'UnmixResult', 'absorbance_ratio',
           'nucleic_acid_and_protein', 'best_wavelengths']


@dataclass
class UnmixResult:
    """
    What was found, and how well it explains the spectrum.

    Attributes
    ----------
    names : tuple of str
        Component names, in the order of every array here.
    amounts : ndarray
        How much of each component. A concentration when the references were
        extinction spectra, in the units implied by them; otherwise relative.
    stderr : ndarray
        Standard error on each amount, from the covariance of the fit. Assumes
        independent, equal-variance noise, which is optimistic -- the residual
        is the better guide to whether the model is right.
    residual : Spectrum
        Measured minus reconstructed. **The thing to look at.** A residual
        that looks like noise means the reference set was sufficient; one with
        structure in it means something is absorbing that was not offered to
        the fit, and the amounts are wrong by however much that something
        overlaps the components that were.
    r_squared : float
        Fraction of variance explained. High values are easy and mean little
        on a smooth spectrum -- two broad bands can reach 0.99 while being
        quantitatively wrong. Read the residual.
    """

    names: tuple
    amounts: np.ndarray
    stderr: np.ndarray
    residual: object
    reconstruction: object
    r_squared: float
    unit: str = ''
    condition_number: float = float('nan')
    #: Path length used to turn the fitted ``c*l`` into a concentration, cm.
    path_length: float = 1.0
    metadata: dict = field(default_factory=dict, repr=False)

    def __len__(self):
        return len(self.names)

    def __iter__(self):
        return iter(zip(self.names, self.amounts))

    def __getitem__(self, name):
        return self.amounts[self.names.index(name)]

    def fractions(self):
        """Amounts as fractions of their total. Meaningless unless the
        references share a unit -- two extinction spectra in M^-1 cm^-1 give
        mole fractions, one in mass units and one molar gives nonsense."""
        total = float(np.sum(self.amounts))
        return self.amounts / total if total else np.zeros_like(self.amounts)

    def __repr__(self):
        parts = ", ".join(f"{name}={value:.4g}"
                          for name, value in zip(self.names, self.amounts))
        return f"<UnmixResult {parts} (R^2={self.r_squared:.4f})>"

    def __str__(self):
        unit = f" {self.unit}" if self.unit else ""
        lines = [f"{'component':<20}{'amount':>12}{'+/-':>12}"]
        for name, amount, error in zip(self.names, self.amounts, self.stderr):
            lines.append(f"{name:<20}{amount:>12.5g}{error:>12.3g}")
        lines.append(f"\nR^2 = {self.r_squared:.5f}{unit}")
        if np.isfinite(self.condition_number):
            lines.append(f"condition number = {self.condition_number:.1f}")
        if self.path_length != 1.0:
            lines.append(f"path length = {self.path_length:g} cm")
        return "\n".join(lines)


def unmix(spectrum, library, *, path_length=None, non_negative=True,
          wavelengths=None):
    """
    Fit a spectrum as a sum of known references.

    Parameters
    ----------
    spectrum : Spectrum
        The measurement to explain. Baseline-correct and, for a scattering
        sample, scatter-correct it first: a sloping background is not one of
        the components, and the fit will distribute it across the ones it has.
    library : Library or sequence of Reference
        What it may be made of. References are resampled onto the spectrum's
        own wavelength grid.
    path_length : float, optional
        Cuvette path length in cm. Beer-Lambert is ``A = eps * c * l``, so
        fitting extinction spectra to an absorbance recovers ``c * l``, and
        the concentration needs the division. Defaults to
        ``spectrum.metadata['path_length']`` if that is set, and otherwise to
        the standard 1 cm.

        Worth stating whenever it is not 1 cm: a 1 mm cuvette left unstated
        makes every concentration here ten times too small, and nothing in the
        result looks wrong. It has no effect at all when the references are
        relative rather than absolute, since there is no concentration to get
        wrong.
    non_negative : bool
        Constrain the amounts to be at least zero, the default and the
        physical choice. Setting it False is occasionally useful for a
        difference spectrum, where a negative coefficient means a component
        was lost rather than gained.
    wavelengths : sequence of float, optional
        Fit at these wavelengths only, rather than the whole spectrum. See
        :func:`best_wavelengths` for choosing them, and note that using the
        whole spectrum is usually better -- restricting to a few points throws
        away exactly the information that would have revealed a missing
        component.

    Returns
    -------
    UnmixResult
    """
    from spectroscopy.library import Library  # noqa: PLC0415
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

    if not isinstance(library, Library):
        library = Library(list(library))
    if not len(library):
        raise ValueError("no references to fit against")

    x = np.asarray(spectrum.x, dtype=float)
    y = np.asarray(spectrum.y, dtype=float)
    if wavelengths is not None:
        indices = np.abs(x[None, :] - np.asarray(wavelengths, float)[:, None]
                         ).argmin(axis=1)
        fit_x, fit_y = x[indices], y[indices]
    else:
        fit_x, fit_y = x, y

    design = library.on(fit_x).T                    # (n_points, n_components)
    if design.shape[0] < design.shape[1]:
        raise ValueError(
            f"{design.shape[1]} components cannot be separated using "
            f"{design.shape[0]} wavelengths; the problem is underdetermined"
        )

    if non_negative:
        amounts, _ = nnls(design, fit_y)
    else:
        amounts, *_ = np.linalg.lstsq(design, fit_y, rcond=None)

    # The fit recovers c*l; Beer-Lambert wants c. Only meaningful when the
    # references are absolute -- with relative references there is no
    # concentration for a path length to be wrong about.
    if path_length is None:
        path_length = float(spectrum.metadata.get('path_length', 1.0))
    if path_length <= 0:
        raise ValueError(f"path length must be positive, got {path_length}")
    absolute = any(reference.is_absolute for reference in library)
    concentrations = amounts / path_length if absolute else amounts

    # Uncertainties from the covariance of the unconstrained problem. With the
    # non-negative constraint active on a component this understates it, which
    # is another reason the residual is the honest diagnostic.
    fitted = design @ amounts
    degrees = max(len(fit_y) - design.shape[1], 1)
    variance = float(np.sum((fit_y - fitted) ** 2)) / degrees
    try:
        covariance = np.linalg.inv(design.T @ design) * variance
        stderr = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    except np.linalg.LinAlgError:
        stderr = np.full(design.shape[1], np.nan)
    if absolute:
        stderr = stderr / path_length          # scales with the amounts

    # Reconstruct across the whole axis even when the fit used a few points,
    # so the residual shows what the model does everywhere and not only where
    # it was asked to agree.
    full = library.on(x).T @ amounts
    total = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - float(np.sum((y - full) ** 2)) / total if total else 0.0

    reconstruction = Spectrum(x, full, x_quantity=spectrum.x_quantity,
                              x_unit=spectrum.x_unit,
                              y_quantity=spectrum.y_quantity,
                              y_unit=spectrum.y_unit)
    reconstruction.name = f"{spectrum.name} (reconstruction)"
    residual = Spectrum(x, y - full, x_quantity=spectrum.x_quantity,
                        x_unit=spectrum.x_unit,
                        y_quantity=f"Residual {spectrum.y_quantity}",
                        y_unit=spectrum.y_unit)
    residual.name = f"{spectrum.name} (residual)"

    units = {reference.unit for reference in library if reference.unit}
    return UnmixResult(
        names=tuple(library.names), amounts=concentrations, stderr=stderr,
        residual=residual, reconstruction=reconstruction,
        r_squared=r_squared,
        unit=(units.pop() if len(units) == 1 else ''),
        condition_number=float(np.linalg.cond(design)),
        path_length=path_length,
        metadata={'non_negative': non_negative,
                  'n_points': len(fit_y),
                  'absolute': absolute,
                  'wavelengths': None if wavelengths is None
                  else list(np.asarray(fit_x))},
    )


def best_wavelengths(library, count=None, *, candidates=None):
    """
    Wavelengths at which the components are most separable.

    Choosing wavelengths by where each component happens to peak is the usual
    approach and it is the wrong criterion: two species can both peak at 280
    and be impossible to tell apart there. What matters is that the rows of
    the extinction matrix point in *different directions*, which is the
    condition number of that matrix -- and that is what this maximises, by
    greedy selection.

    Returns ``count`` wavelengths, defaulting to one per component, which is
    the minimum that can separate them at all.

    Use the whole spectrum when you can. This exists for the cases where you
    cannot -- a filter photometer, a plate reader with fixed wavelengths, or
    an instrument where reading every wavelength costs time you do not have.
    """
    from spectroscopy.library import Library  # noqa: PLC0415

    if not isinstance(library, Library):
        library = Library(list(library))
    n_components = len(library)
    if count is None:
        count = n_components
    if count < n_components:
        raise ValueError(
            f"{count} wavelengths cannot separate {n_components} components"
        )

    if candidates is None:
        candidates = np.asarray(next(iter(library)).spectrum.x, dtype=float)
    candidates = np.asarray(candidates, dtype=float)
    matrix = library.on(candidates).T               # (n_candidates, n_comp)

    chosen = []
    for _ in range(count):
        best_index, best_score = None, -np.inf
        for index in range(len(candidates)):
            if index in chosen:
                continue
            trial = matrix[chosen + [index], :]
            if trial.shape[0] < n_components:
                # Not yet square: keep the rows as mutually distinct as
                # possible, measured by the smallest singular value.
                score = float(np.linalg.svd(trial, compute_uv=False)[-1])
            else:
                singular = np.linalg.svd(trial, compute_uv=False)
                score = float(singular[-1] / singular[0])   # 1 / condition
            if score > best_score:
                best_index, best_score = index, score
        chosen.append(best_index)

    return np.sort(candidates[np.asarray(chosen)])


def absorbance_ratio(spectrum, numerator=260.0, denominator=280.0):
    """
    The ratio of absorbance at two wavelengths, e.g. A260/A280.

    Provided because it is asked for constantly, with the caveat attached
    because it is worth less than its popularity suggests.

    What it is good for: a quick, single-number flag that something is not
    what you think it is. Pure double-stranded DNA gives about 1.8, RNA about
    2.0, pure protein about 0.57.

    What it is not good for is quantifying the contamination, and that is the
    use it usually gets. The ratio is dominated by whichever species absorbs
    more, and nucleic acid absorbs far more per unit mass at 260 than protein
    does at 280 -- so a preparation that is a few percent nucleic acid by mass
    reads as substantially contaminated, while one that is heavily protein
    contaminated can still sit near 1.8. Two numbers from a spectrum of
    hundreds of points cannot say more than that.

    :func:`nucleic_acid_and_protein` uses the whole spectrum instead, and its
    residual will tell you when neither component fits -- which the ratio
    never can.
    """
    x = np.asarray(spectrum.x, dtype=float)
    y = np.asarray(spectrum.y, dtype=float)
    top = float(y[int(np.abs(x - numerator).argmin())])
    bottom = float(y[int(np.abs(x - denominator).argmin())])
    if bottom == 0:
        raise ValueError(
            f"absorbance at {denominator} nm is zero; the ratio is undefined"
        )
    return top / bottom


def nucleic_acid_and_protein(spectrum, nucleic_acid, protein, **kwargs):
    """
    The two-component case of :func:`unmix`, named for what it is usually for.

    ``nucleic_acid`` and ``protein`` are References -- extinction spectra if
    the amounts are to be concentrations. There is nothing here that
    :func:`unmix` does not do; it exists so the common case reads clearly and
    so the docstring has somewhere to say that this is strictly better than
    the ratio, because it uses every wavelength and shows you a residual.
    """
    from spectroscopy.library import Library  # noqa: PLC0415

    return unmix(spectrum, Library([nucleic_acid, protein]), **kwargs)
