# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Curve fitting: decomposing a band into overlapping components.

Peak *detection* answers "where are the maxima" (:mod:`spectroscopy.peaks`).
Fitting answers the harder question -- "this envelope is several overlapping
bands; how much of each?" -- which is what amide I secondary structure, band
area ratios and any quantitative deconvolution actually need.

The return type is a dataclass over numpy arrays rather than a DataFrame, for
the reason in review section 5.6: a pandas return type obliges every caller to
install pandas to do anything with the result.

ADR-0001 section 6.2 deferred this until there was a real caller to design it
against, on the grounds that a return type shaped by an imagined use fits real
ones badly. Secondary structure analysis is that caller.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from spectroscopy import lineshapes
from spectroscopy.processing import common

__all__ = ['FitResult', 'fit_components', 'MODELS']

#: Area of a unit-amplitude, unit-width component. Gaussian:
#: ``sqrt(pi / (4 ln 2))``. Lorentzian, for the parameterisation in
#: :mod:`spectroscopy.lineshapes`: ``pi / 2``.
_GAUSSIAN_AREA = math.sqrt(math.pi / (4.0 * math.log(2.0)))
_LORENTZIAN_AREA = math.pi / 2.0

#: Component shapes, by name. Each takes ``(x, position, fwhm, amplitude)``
#: plus, for pseudo-Voigt, a Gaussian fraction.
MODELS = ('gaussian', 'lorentzian', 'voigt')


def _evaluate(x, position, fwhm, amplitude, model, eta=None):
    """One component, evaluated on x."""
    if model == 'gaussian':
        return lineshapes.gauss(x, position, fwhm, amplitude)
    if model == 'lorentzian':
        return lineshapes.lorentz(x, position, fwhm, amplitude)
    if model == 'voigt':
        return lineshapes.spec_comp(x, position, fwhm, amplitude, eta)
    raise ValueError(f"unknown model {model!r}; known models are {MODELS}")


def _component_area(fwhm, amplitude, model, eta=None):
    """Analytic area, rather than a numerical integral of a sampled curve."""
    gaussian = _GAUSSIAN_AREA * fwhm * amplitude
    lorentzian = _LORENTZIAN_AREA * fwhm * amplitude
    if model == 'gaussian':
        return gaussian
    if model == 'lorentzian':
        return lorentzian
    return eta * gaussian + (1.0 - eta) * lorentzian


@dataclass
class FitResult:
    """
    A band decomposed into components.

    Attributes
    ----------
    position, fwhm, amplitude : ndarray
        The fitted parameters of each component, in the spectrum's units.
    area : ndarray
        Analytic area of each component -- the quantity secondary structure
        and band-ratio analyses actually use. Computed from the fitted shape,
        not by integrating a sampled curve, so it does not depend on the
        spacing of the x axis.
    eta : ndarray or None
        Gaussian fraction per component, for ``model='voigt'``.
    stderr : dict of ndarray or None
        One-sigma uncertainties on ``position``, ``fwhm`` and ``amplitude``,
        from the covariance matrix. ``None`` when the fit did not produce a
        usable one -- which is itself worth noticing.

        **Read these.** A weak component between two strong ones is barely
        determined by the data: with 0.5 % noise on a synthetic amide I
        envelope, a shoulder at 1645 cm^-1 drifts by 4 cm^-1 while its
        neighbours hold to 0.1, and its uncertainty grows by more than tenfold
        to say so. Fitting cannot recover information the spectrum does not
        contain; reporting a position without its uncertainty hides that.
    model : str
        ``'gaussian'``, ``'lorentzian'`` or ``'voigt'``.
    x, residual : ndarray
        The axis fitted over, and (data - fit) on it.
    rmse, r_squared : float
        Goodness of fit. ``r_squared`` is on the fitted window only.
    x_unit, y_unit, source : str
        Carried through so a result can be labelled without its spectrum.
    """

    position: np.ndarray
    fwhm: np.ndarray
    amplitude: np.ndarray
    area: np.ndarray
    model: str
    x: np.ndarray = field(repr=False)
    residual: np.ndarray = field(repr=False)
    rmse: float = 0.0
    r_squared: float = 0.0
    eta: np.ndarray | None = None
    stderr: dict | None = field(default=None, repr=False)
    x_unit: str = ''
    y_unit: str = ''
    source: str | None = None

    def __len__(self) -> int:
        return len(self.position)

    def __repr__(self) -> str:
        where = f" of {self.source}" if self.source else ""
        return (f"<FitResult: {len(self)} {self.model} components{where}, "
                f"R^2={self.r_squared:.4f}>")

    def __str__(self) -> str:
        unit = f" {self.x_unit}" if self.x_unit else ""
        lines = [f"{len(self)} {self.model} components, "
                 f"R^2 = {self.r_squared:.4f}, RMSE = {self.rmse:.4g}",
                 f"{'position':>12}{unit} {'fwhm':>10} {'amplitude':>12}"
                 f" {'area':>12} {'fraction':>9}"]
        for position, width, height, area, fraction in zip(
                self.position, self.fwhm, self.amplitude,
                self.area, self.fractions()):
            lines.append(f"{position:12.1f}{' ' * len(unit)} {width:10.2f} "
                         f"{height:12.5g} {area:12.5g} {100 * fraction:8.1f}%")
        return "\n".join(lines)

    # -- derived quantities -------------------------------------------------

    def fractions(self) -> np.ndarray:
        """Each component's area as a fraction of the total.

        This is the number secondary-structure analysis reports, and the one
        band-ratio work needs. Areas can be negative if a fit was allowed to
        go there, so this is a plain ratio and not forced onto [0, 1] -- a
        fraction outside it means the fit is wrong and should say so.
        """
        total = self.area.sum()
        if total == 0:
            return np.zeros_like(self.area)
        return self.area / total

    def components(self, x=None) -> np.ndarray:
        """Each component evaluated separately, shape (n_components, n_points)."""
        x = self.x if x is None else np.asarray(x, dtype=float)
        etas = self.eta if self.eta is not None else [None] * len(self)
        return np.array([
            _evaluate(x, position, width, height, self.model, eta)
            for position, width, height, eta
            in zip(self.position, self.fwhm, self.amplitude, etas)
        ])

    def curve(self, x=None) -> np.ndarray:
        """The sum of the components -- the fitted envelope."""
        return self.components(x).sum(axis=0)

    def assign(self, bands) -> dict:
        """
        Group component areas by named x ranges.

        Parameters
        ----------
        bands : dict
            ``{name: (low, high)}``. A component is assigned to the first
            range containing its position.

        Returns
        -------
        dict
            ``{name: summed fraction}``, plus ``'unassigned'`` if any
            component fell outside every range -- never silently dropped,
            because a band that matches nothing is the interesting one.

        Examples
        --------
        The amide I assignment, which is the reason this method exists::

            result.assign({'beta-sheet':   (1620, 1640),
                           'random':       (1640, 1650),
                           'alpha-helix':  (1650, 1658),
                           'turns':        (1660, 1680)})
        """
        fractions = self.fractions()
        totals = {name: 0.0 for name in bands}
        unassigned = 0.0
        for position, fraction in zip(self.position, fractions):
            for name, (low, high) in bands.items():
                if min(low, high) <= position <= max(low, high):
                    totals[name] += fraction
                    break
            else:
                unassigned += fraction
        if unassigned:
            totals['unassigned'] = unassigned
        return totals

    # -- interop ------------------------------------------------------------

    def to_dataframe(self):
        """As a pandas DataFrame. Imports pandas inside the method."""
        import pandas as pd  # pylint: disable=C0415

        columns = {'position': self.position, 'fwhm': self.fwhm,
                   'amplitude': self.amplitude, 'area': self.area,
                   'fraction': self.fractions()}
        if self.eta is not None:
            columns['eta'] = self.eta
        if self.stderr is not None:
            columns['position_stderr'] = self.stderr['position']
        return pd.DataFrame(columns)


def _uniform_spacing(x):
    """The step of a uniform axis, or a refusal explaining why it matters."""
    steps = np.diff(x)
    spacing = float(np.mean(steps))
    if not np.allclose(steps, spacing, rtol=1e-3):
        raise ValueError(
            "derivative weighting needs a uniformly spaced axis: the "
            "Savitzky-Golay filter assumes one, and on an uneven axis the "
            "second derivative is wrong by a factor that varies along the "
            "spectrum. Resample first, e.g. "
            "spectrum.resample(np.linspace(x.min(), x.max(), len(x)))."
        )
    return spacing


def fit_components(x, y, positions, *, model='voigt', fwhm=None,
                   position_tolerance=None, max_fwhm=None,
                   non_negative=True, derivative_weight=0.0,
                   derivative_window=11, derivative_polyorder=3,
                   maxfev=10000) -> FitResult:
    """
    Fit overlapping components to (x, y).

    Parameters
    ----------
    x, y : array_like
        The data to fit. Usually already cropped to the band of interest and
        baseline-corrected -- fitting a band that sits on a slope puts the
        slope into the component areas.
    positions : array_like
        Starting positions, one per component. Where they come from matters
        more than the fitting does: for amide I they are the minima of the
        second derivative, not the maxima of the spectrum, because the
        overlapping bands have no maxima of their own.
    model : {'voigt', 'gaussian', 'lorentzian'}, default 'voigt'
        The component shape. **The default is pseudo-Voigt deliberately, and
        picking a pure shape is usually a mistake.**

        The wrong shape does not announce itself. Fitting four bands generated
        with a 50/50 Voigt shape using pure Gaussians gives R^2 = 0.978 -- a fit
        that looks fine on a plot -- while the area fractions are wrong by 33
        percentage points, one 17 % component coming out at 50 %. Pure
        Lorentzian bands fitted as Gaussians reach R^2 = 0.981 with a 62-point
        error. **Goodness of fit does not validate the lineshape**, which makes
        a fixed shape a silent error rather than a visible one.

        Letting the mixing float costs nothing. With 0.5 % noise on genuinely
        Gaussian bands, fitting as Voigt is as accurate as fitting as Gaussian
        (2.3 against 2.7 points of composition error over 20 trials); on
        intermediate bands it is eight times better (3.5 against 28). The extra
        parameter per component pays for itself immediately.

        Use a pure shape only when there is a physical reason to fix it, or to
        show that it does not matter for a particular band.
    fwhm : float or array_like, optional
        Starting width(s). Defaults to the mean spacing between the given
        positions, which is a reasonable guess for a crowded band.
    position_tolerance : float, optional
        How far each component may move from its starting position. Amide I
        work normally constrains this hard (a few cm^-1); leaving it free
        lets two components converge onto one strong band and report a
        structure that is not there.
    max_fwhm : float, optional
        Upper bound on width. Defaults to the fitted range, which stops a
        component flattening into a pseudo-baseline.
    non_negative : bool, default True
        Constrain amplitudes to be >= 0. Turn it off only when fitting a
        difference spectrum, where negative bands are real.
    derivative_weight : float, default 0
        Fit the spectrum **and** its second derivative together, with this
        relative weight on the derivative. ``0`` fits the spectrum alone.

        This is what makes a crowded amide I fit trustworthy, and the reason
        is worth knowing, because it decides when to use it.

        **The second derivative annihilates a smooth background.** A residual
        slope or curvature left behind by imperfect water subtraction is
        invisible to it, while the absorbance envelope absorbs that residual
        straight into the component areas. Measured on two Gaussians 10 cm^-1
        apart -- close enough that the envelope has a single maximum -- with
        0.4 % noise, 30 trials, median error:

        =========================  ============  =============
        Background under the band  ``weight=0``  ``weight=2``
        =========================  ============  =============
        none                        0.12 cm^-1    0.16 cm^-1
        linear residual             3.37 cm^-1    0.72 cm^-1
        curved residual             4.20 cm^-1    0.63 cm^-1
        =========================  ============  =============

        The composition errors move with them: on the curved background the
        area fractions are wrong by 31 % unweighted and 5 % at ``weight=2``.

        So: **on a perfectly corrected band it costs a little** -- the model is
        already exactly right and Savitzky-Golay amplifies noise -- **and on a
        real one it is worth a factor of five**. Real amide I data has residual
        background, so the useful default is non-zero.

        Both blocks are scaled to unit magnitude before weighting, so ``1.0``
        means "the derivative matters as much as the spectrum" whatever the
        units. Improvement saturates around 4 (0.46 cm^-1 at 4, 0.45 at 8);
        **1 to 4 is the useful range**.
    derivative_window, derivative_polyorder : int
        Savitzky-Golay parameters for that derivative. **The same filter is
        applied to the model as to the data**, rather than differentiating the
        model analytically: SG smooths as it differentiates, so an analytic
        model derivative would be systematically sharper than the measured one
        and the fit would compensate by narrowing every band.
    maxfev : int
        Passed to :func:`scipy.optimize.curve_fit`.

    Returns
    -------
    FitResult

    Raises
    ------
    ValueError
        If no positions are given, or a position lies outside the data.
    RuntimeError
        If the fit does not converge -- not silently returned as a bad fit.
    """
    from scipy.optimize import curve_fit  # pylint: disable=C0415

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    positions = np.atleast_1d(np.asarray(positions, dtype=float))
    if positions.size == 0:
        raise ValueError("fitting needs at least one starting position")
    if model not in MODELS:
        raise ValueError(f"unknown model {model!r}; known models are {MODELS}")

    low, high = float(np.min(x)), float(np.max(x))
    outside = positions[(positions < low) | (positions > high)]
    if outside.size:
        raise ValueError(
            f"starting positions {outside.tolist()} lie outside the data "
            f"({low:g} to {high:g}). Crop the spectrum to the band you mean "
            "to fit, or check the units."
        )

    count = len(positions)
    span = high - low
    if fwhm is None:
        spacing = span / count if count > 1 else span / 4.0
        widths = np.full(count, spacing)
    else:
        widths = np.broadcast_to(np.asarray(fwhm, dtype=float), (count,)).copy()
    heights = np.interp(positions, x, y)
    # An all-zero or negative starting amplitude gives the optimiser nothing
    # to work with; fall back to the data range.
    fallback = max(float(np.max(y) - np.min(y)), 1e-12)
    heights = np.where(np.abs(heights) < 1e-12, fallback, heights)

    per_component = 4 if model == 'voigt' else 3
    guess = np.empty(count * per_component)
    lower = np.empty_like(guess)
    upper = np.empty_like(guess)

    width_ceiling = max_fwhm if max_fwhm is not None else span
    for index in range(count):
        base = index * per_component
        guess[base:base + 3] = (positions[index], widths[index], heights[index])
        if position_tolerance is None:
            lower[base], upper[base] = low, high
        else:
            lower[base] = positions[index] - position_tolerance
            upper[base] = positions[index] + position_tolerance
        lower[base + 1], upper[base + 1] = 1e-9, width_ceiling
        lower[base + 2] = 0.0 if non_negative else -np.inf
        upper[base + 2] = np.inf
        if model == 'voigt':
            guess[base + 3] = 0.5
            lower[base + 3], upper[base + 3] = 0.0, 1.0

    # Starting values must lie inside the bounds or curve_fit refuses to start;
    # a caller-supplied fwhm larger than the window is an easy way to trip it.
    guess = np.clip(guess, lower, upper)

    def envelope(x_values, *parameters):
        total = np.zeros_like(x_values)
        for component in range(count):
            base = component * per_component
            eta = parameters[base + 3] if model == 'voigt' else None
            total = total + _evaluate(x_values, parameters[base],
                                      parameters[base + 1],
                                      parameters[base + 2], model, eta)
        return total

    if derivative_weight:
        spacing = _uniform_spacing(x)

        def second_derivative(values):
            return common.derivative(values, order=2,
                                     window_length=derivative_window,
                                     polyorder=derivative_polyorder,
                                     delta=spacing)

        # Scale each block to unit magnitude so derivative_weight is a true
        # relative weight and not an accident of the units: d2A/dv2 on a
        # cm^-1 axis is smaller than A by roughly the square of a band width.
        target_derivative = second_derivative(y)
        amplitude_scale = max(float(np.max(np.abs(y))), 1e-30)
        derivative_scale = max(float(np.max(np.abs(target_derivative))), 1e-30)

        observed = np.concatenate([
            y / amplitude_scale,
            derivative_weight * target_derivative / derivative_scale,
        ])

        def stacked(x_values, *parameters):
            curve = envelope(x_values, *parameters)
            return np.concatenate([
                curve / amplitude_scale,
                derivative_weight * second_derivative(curve) / derivative_scale,
            ])

        objective, target = stacked, observed
    else:
        objective, target = envelope, y

    try:
        best, covariance = curve_fit(objective, x, target, p0=guess,
                                     bounds=(lower, upper), maxfev=maxfev)
    except RuntimeError as error:
        raise RuntimeError(
            f"the fit did not converge after {maxfev} evaluations. "
            "Common causes: starting positions far from the real bands, too "
            "many components for the data, or a baseline that was not "
            "removed first."
        ) from error

    best = np.asarray(best)
    fitted_positions = best[0::per_component]
    fitted_widths = best[1::per_component]
    fitted_heights = best[2::per_component]
    etas = best[3::per_component] if model == 'voigt' else None

    areas = np.array([
        _component_area(width, height, model,
                        None if etas is None else etas[index])
        for index, (width, height) in enumerate(zip(fitted_widths, fitted_heights))
    ])

    residual = y - envelope(x, *best)
    rmse = float(np.sqrt(np.mean(residual ** 2)))
    variance = float(np.sum((y - np.mean(y)) ** 2))
    r_squared = 1.0 - float(np.sum(residual ** 2)) / variance if variance else 0.0

    stderr = None
    if covariance is not None and np.all(np.isfinite(covariance)):
        sigma = np.sqrt(np.abs(np.diag(covariance)))
        stderr = {'position': sigma[0::per_component],
                  'fwhm': sigma[1::per_component],
                  'amplitude': sigma[2::per_component]}

    order = np.argsort(fitted_positions)
    return FitResult(
        position=fitted_positions[order],
        fwhm=fitted_widths[order],
        amplitude=fitted_heights[order],
        area=areas[order],
        eta=None if etas is None else etas[order],
        stderr=None if stderr is None else {k: v[order] for k, v in stderr.items()},
        model=model,
        x=x,
        residual=residual,
        rmse=rmse,
        r_squared=r_squared,
    )
