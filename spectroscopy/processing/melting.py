# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Thermal unfolding: melting temperature and van 't Hoff enthalpy.

A melt is a series of spectra measured along temperature, and the analysis is
the same whatever the technique: CD at 222 nm, tryptophan fluorescence,
absorbance. So this takes a signal against temperature and knows nothing about
which instrument produced it.

The model
---------
Two states, folded and unfolded, in equilibrium::

    K(T) = exp(-dH/R * (1/T - 1/Tm))        unfolding equilibrium constant
    f_U  = K / (1 + K)                      fraction unfolded
    y(T) = y_F(T) (1 - f_U) + y_U(T) f_U

with **sloping baselines** ``y_F = a_f + b_f T`` and ``y_U = a_u + b_u T``.
The slopes are not decoration: the signal of each state depends on temperature
in its own right, and a fit that holds the baselines flat pushes that
dependence into the transition, moving Tm and badly distorting dH. Six
parameters for what looks like a two-parameter problem, and all six are
needed.

Why the fit is not the evidence
-------------------------------
A two-state curve fits almost anything sigmoidal, including a three-state
transition with a short-lived intermediate. R^2 will be excellent either way,
which is the trap this project keeps meeting -- see the amide I estimator,
which reached R^2 = 0.999 while spreading the same protein over twenty
percentage points.

So two independent checks are provided, and they use information the fit does
not: :func:`isodichroic_point` and :func:`two_state_rank`. Both look across
the whole spectrum rather than at one wavelength, and both can fail while the
fit still looks perfect. **Report them beside any Tm.**
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np

__all__ = ['MeltResult', 'two_state', 'from_collection',
           'isodichroic_point', 'two_state_rank', 'GAS_CONSTANT']

#: J / (mol K). CODATA, exact since the 2019 SI redefinition.
GAS_CONSTANT = 8.314462618

#: Above this, a "Celsius" temperature is almost certainly Kelvin.
_IMPLAUSIBLE_CELSIUS = 150.0


@dataclass
class MeltResult:
    """
    A fitted two-state thermal transition.

    Attributes
    ----------
    tm : float
        Melting temperature, degrees C -- where half the protein is unfolded.
    enthalpy : float
        van 't Hoff enthalpy of unfolding, kJ/mol. This is the enthalpy the
        *transition shape* implies, which equals the calorimetric enthalpy
        only if the transition really is two-state. A van 't Hoff enthalpy
        well below the calorimetric one is the classic sign that it is not.
    tm_stderr, enthalpy_stderr : float
        One standard error, from the fit covariance. They describe the
        precision of the fit, not whether the model is right.
    folded, unfolded : tuple of float
        ``(intercept, slope)`` of each baseline, in signal units per degree.
    residual_rms : float
        Fit residual, in signal units.
    temperature, signal : ndarray
        The data, as fitted.
    """

    tm: float
    enthalpy: float
    tm_stderr: float
    enthalpy_stderr: float
    folded: tuple
    unfolded: tuple
    residual_rms: float
    temperature: np.ndarray = field(repr=False, default=None)
    signal: np.ndarray = field(repr=False, default=None)
    source: str | None = None

    def fraction_unfolded(self, temperature=None):
        """Fraction unfolded at ``temperature`` (default: the measured grid)."""
        t = self.temperature if temperature is None else np.asarray(
            temperature, dtype=float)
        kelvin = np.asarray(t, dtype=float) + 273.15
        tm_kelvin = self.tm + 273.15
        exponent = -(self.enthalpy * 1000.0 / GAS_CONSTANT) * (
            1.0 / kelvin - 1.0 / tm_kelvin)
        # exp overflows long before the fraction stops being 0 or 1.
        return 1.0 / (1.0 + np.exp(-np.clip(exponent, -700, 700)))

    def fitted(self, temperature=None):
        """The model curve, in signal units."""
        t = self.temperature if temperature is None else np.asarray(
            temperature, dtype=float)
        return _model(t, self.tm, self.enthalpy,
                      *self.folded, *self.unfolded)

    @property
    def entropy(self):
        """Unfolding entropy at Tm, J/(mol K). ``dS = dH / Tm``."""
        return self.enthalpy * 1000.0 / (self.tm + 273.15)

    def free_energy(self, temperature=25.0):
        """
        Unfolding free energy at ``temperature`` (default 25 C), kJ/mol.

        The Gibbs-Helmholtz form **without** a heat-capacity term, so it is
        only reliable near Tm. Extrapolating tens of degrees to quote a
        stability at 25 C from a melt alone is a well-known way to get a
        confident wrong number; dCp cannot be had from a single melting curve.
        """
        kelvin = np.asarray(temperature, dtype=float) + 273.15
        return self.enthalpy * (1.0 - kelvin / (self.tm + 273.15))

    def __str__(self) -> str:
        where = f" of {self.source}" if self.source else ""
        return (f"Two-state melt{where}:\n"
                f"  Tm        {self.tm:8.2f} +/- {self.tm_stderr:.2f} C\n"
                f"  dH        {self.enthalpy:8.1f} +/- "
                f"{self.enthalpy_stderr:.1f} kJ/mol\n"
                f"  dS(Tm)    {self.entropy:8.1f} J/(mol K)\n"
                f"  residual  {self.residual_rms:8.4g} (signal units)")


def _model(temperature, tm, enthalpy, folded_intercept, folded_slope,
           unfolded_intercept, unfolded_slope):
    """Two-state signal with linear pre- and post-transition baselines."""
    kelvin = np.asarray(temperature, dtype=float) + 273.15
    tm_kelvin = tm + 273.15
    exponent = -(enthalpy * 1000.0 / GAS_CONSTANT) * (
        1.0 / kelvin - 1.0 / tm_kelvin)
    fraction = 1.0 / (1.0 + np.exp(-np.clip(exponent, -700, 700)))
    lower = folded_intercept + folded_slope * temperature
    upper = unfolded_intercept + unfolded_slope * temperature
    return lower * (1.0 - fraction) + upper * fraction


def _initial_guess(temperature, signal):
    """Baselines from the ends, Tm from the midpoint, dH from the steepness."""
    order = np.argsort(temperature)
    t, y = temperature[order], signal[order]
    edge = max(2, len(t) // 5)

    folded = np.polyfit(t[:edge], y[:edge], 1)[::-1]        # (intercept, slope)
    unfolded = np.polyfit(t[-edge:], y[-edge:], 1)[::-1]

    midpoint = 0.5 * (y[:edge].mean() + y[-edge:].mean())
    tm = float(np.interp(0.5, np.linspace(0, 1, len(t)),
                         t) if y[0] == y[-1] else
               t[int(np.argmin(np.abs(y - midpoint)))])

    # At Tm the transition's slope is dH / (4 R Tm^2) in fraction per kelvin.
    span = y[-edge:].mean() - y[:edge].mean()
    with np.errstate(invalid='ignore', divide='ignore'):
        steepest = np.nanmax(np.abs(np.gradient(y, t)))
    tm_kelvin = tm + 273.15
    enthalpy = 200.0
    if span and np.isfinite(steepest) and steepest > 0:
        per_kelvin = steepest / abs(span)
        enthalpy = 4.0 * GAS_CONSTANT * tm_kelvin ** 2 * per_kelvin / 1000.0
        enthalpy = float(np.clip(enthalpy, 20.0, 2000.0))
    return [tm, enthalpy, folded[0], folded[1], unfolded[0], unfolded[1]]


def two_state(temperature, signal, *, source=None) -> MeltResult:
    """
    Fit a two-state transition to a signal measured against temperature.

    Parameters
    ----------
    temperature : array_like
        Degrees **Celsius**. Kelvin input is detected and refused rather than
        fitted: it produces a Tm near 330 with an unremarkable-looking
        residual.
    signal : array_like
        Any monotonic-in-state observable -- ellipticity at 222 nm, a
        fluorescence ratio, an absorbance. Direction does not matter; the
        baselines take care of it.
    source : str, optional
        What the data came from, carried into the result.

    Returns
    -------
    MeltResult
    """
    from scipy.optimize import curve_fit  # noqa: PLC0415

    temperature = np.asarray(temperature, dtype=float).ravel()
    signal = np.asarray(signal, dtype=float).ravel()
    if temperature.size != signal.size:
        raise ValueError(
            f"{temperature.size} temperatures for {signal.size} signal "
            f"values; they are matched by position"
        )
    if temperature.size < 10:
        raise ValueError(
            f"a two-state fit has six parameters and there are "
            f"{temperature.size} points. Below ten the baselines and the "
            f"transition cannot be told apart: on a real seven-point melt "
            f"this returned a Tm with a standard error of 1e17, which is the "
            f"fit saying it has no idea while still printing a number. "
            f"Fifteen or more is where the answer starts to mean something."
        )
    if np.nanmin(temperature) > _IMPLAUSIBLE_CELSIUS:
        raise ValueError(
            f"temperatures start at {np.nanmin(temperature):g}, which is "
            f"Kelvin rather than Celsius. Pass Celsius: the fit would "
            f"otherwise succeed and report a Tm around "
            f"{np.nanmedian(temperature):.0f} C."
        )

    guess = _initial_guess(temperature, signal)
    try:
        popt, pcov = curve_fit(_model, temperature, signal, p0=guess,
                               maxfev=20000)
    except RuntimeError as error:
        raise RuntimeError(
            f"the two-state fit did not converge ({error}). The usual causes "
            f"are a transition that is not complete at either end -- so a "
            f"baseline is unconstrained -- or a series too short to separate "
            f"the baselines from the transition."
        ) from error

    residual = signal - _model(temperature, *popt)
    errors = np.sqrt(np.diag(pcov)) if np.all(np.isfinite(pcov)) else \
        np.full(len(popt), np.nan)

    result = MeltResult(
        tm=float(popt[0]), enthalpy=float(popt[1]),
        tm_stderr=float(errors[0]), enthalpy_stderr=float(errors[1]),
        folded=(float(popt[2]), float(popt[3])),
        unfolded=(float(popt[4]), float(popt[5])),
        residual_rms=float(np.sqrt(np.mean(residual ** 2))),
        temperature=temperature, signal=signal, source=source)

    # Whether the transition was actually observed, which is a better
    # question than whether Tm happens to fall inside the range: with an
    # unobserved transition the baselines can fit the data exactly and put Tm
    # anywhere, inside the range included.
    traversed = float(result.fraction_unfolded(temperature.max())
                      - result.fraction_unfolded(temperature.min()))
    span_measured = float(temperature.max() - temperature.min())
    if not np.isfinite(result.tm_stderr) or result.tm_stderr > span_measured:
        warnings.warn(
            f"the fit is degenerate: Tm came back as {result.tm:.1f} +/- "
            f"{result.tm_stderr:.3g} C, an uncertainty larger than the "
            f"{span_measured:g} C that was measured. The six parameters are "
            f"not determined by this data, and the value above should not be "
            f"quoted.",
            UserWarning, stacklevel=2)
    elif abs(traversed) < 0.5:
        warnings.warn(
            f"only {100 * abs(traversed):.0f}% of the transition happened "
            f"inside the measured range {temperature.min():g}-"
            f"{temperature.max():g} C, so Tm = {result.tm:.1f} C is an "
            f"extrapolation from a curve that never turned over. Measure "
            f"further either side; neither Tm nor dH means much here.",
            UserWarning, stacklevel=2)
    elif result.enthalpy_stderr > 0.5 * abs(result.enthalpy):
        warnings.warn(
            f"dH = {result.enthalpy:.0f} +/- {result.enthalpy_stderr:.0f} "
            f"kJ/mol is barely determined, usually because one baseline is "
            f"too short to pin down. Extend the temperature range.",
            UserWarning, stacklevel=2)
    return result


def _signal_from(collection, wavelength):
    """One number per spectrum: at a wavelength, or the leading SVD score."""
    x, matrix, temperature = collection.to_matrix(with_parameter=True)
    if wavelength is not None:
        return temperature, matrix[:, int(np.argmin(np.abs(x - wavelength)))]
    # Whole-spectrum: the first left singular vector, which uses every
    # wavelength and so is far less noisy than any single one. Its sign and
    # scale are arbitrary, which does not matter -- the baselines absorb both.
    centred = matrix - matrix.mean(axis=0)
    scores, values, _ = np.linalg.svd(centred, full_matrices=False)
    return temperature, scores[:, 0] * values[0]


def from_collection(collection, wavelength=None, **kwargs) -> MeltResult:
    """
    Fit a melt from a :class:`~spectroscopy.collection.SpectrumCollection`
    whose parameter is temperature.

    ``wavelength=None`` uses the **whole spectrum** through its leading
    singular vector rather than one wavelength. That is usually the better
    choice: it uses every point, so it is much less noisy, and it does not
    require guessing which wavelength reports on unfolding. Give a wavelength
    when you want the classical single-wavelength number -- 222 nm for helix
    content -- or when only part of the spectrum is trustworthy.
    """
    temperature, signal = _signal_from(collection, wavelength)
    return two_state(temperature, signal,
                     source=collection.name or 'collection', **kwargs)


def isodichroic_point(collection, region=None):
    """
    Find the wavelength at which every spectrum in the series crosses.

    **The two-state test that the fit cannot do.** If a sample really passes
    between two states and nothing else, then every spectrum is a weighted
    average of the same two, so they all pass through one point -- the
    isodichroic point in CD, isosbestic in absorbance. An intermediate state
    smears that crossing out.

    This is evidence of a different kind from a good fit, because it uses the
    whole spectrum at every temperature, and a three-state transition can fit
    a two-state curve at one wavelength perfectly while having no isodichroic
    point at all.

    Returns
    -------
    dict
        ``wavelength`` of tightest crossing, the ``spread`` of signals there
        in signal units, that spread ``relative`` to the spectra's full range,
        and ``is_tight`` -- ``True`` when the crossing is under 5% of range,
        which is the rule of thumb, not a law.
    """
    x, matrix = collection.to_matrix()
    temperature = collection.parameters
    if region is not None:
        low, high = sorted(region)
        inside = (x >= low) & (x <= high)
        x, matrix = x[inside], matrix[:, inside]
    if matrix.shape[0] < 3:
        raise ValueError(
            f"an isodichroic point needs at least three spectra to be "
            f"meaningful; got {matrix.shape[0]}"
        )

    # A crossing is where the series *inverts*, not merely where the spectra
    # happen to be close together. Correlating signal against the parameter
    # wavelength by wavelength finds that: strongly one sign on one side,
    # strongly the other on the other, passing through zero at the crossing.
    #
    # Spread alone is not enough, and real data showed why. On an AqpZ melt it
    # picked 238.8 nm -- the red end, where every spectrum has decayed towards
    # zero and so they trivially agree. Spectra agreeing because there is no
    # signal is the absence of information, not evidence of two states.
    finite = np.isfinite(temperature)
    ordering = np.zeros(matrix.shape[1])
    if finite.sum() >= 3:
        centred_t = temperature[finite] - temperature[finite].mean()
        centred_y = matrix[finite] - matrix[finite].mean(axis=0)
        denominator = (np.sqrt(np.sum(centred_t ** 2))
                       * np.sqrt(np.sum(centred_y ** 2, axis=0)))
        with np.errstate(invalid='ignore', divide='ignore'):
            ordering = np.where(denominator > 0,
                                (centred_t @ centred_y) / denominator, 0.0)

    spread = matrix.max(axis=0) - matrix.min(axis=0)
    span = float(matrix.max() - matrix.min()) or 1.0

    # Only where there is signal to cross: the mean spectrum must carry a
    # decent fraction of its own maximum.
    mean_amplitude = np.abs(matrix.mean(axis=0))
    eligible = mean_amplitude > 0.25 * float(mean_amplitude.max())
    # ... and where the ordering actually inverts nearby, which is what makes
    # it a crossing rather than a convergence.
    inverts = np.zeros_like(eligible)
    signs = np.sign(ordering)
    changes = np.flatnonzero(np.diff(signs) != 0)
    for index in changes:
        inverts[max(index - 2, 0):index + 3] = True
    candidates = eligible & inverts
    if not candidates.any():
        candidates = eligible if eligible.any() else np.ones_like(eligible)

    index = int(np.flatnonzero(candidates)[np.argmin(spread[candidates])])
    return {
        'wavelength': float(x[index]),
        'spread': float(spread[index]),
        'relative': float(spread[index] / span),
        'is_tight': bool(spread[index] / span < 0.05 and candidates[index]
                         and inverts[index]),
        'inverts_here': bool(inverts[index]),
        'signal_here': float(mean_amplitude[index]
                             / (float(mean_amplitude.max()) or 1.0)),
    }


#: How far the third singular value must stand above the noise floor before a
#: third species is called. A rule of thumb, not a law.
THIRD_COMPONENT_THRESHOLD = 3.0


def two_state_rank(collection, region=None):
    """
    Count the species contributing to the series, by SVD.

    A mixture of two spectra in varying proportion is a rank-2 matrix, so a
    two-state melt has two singular values above the noise and the rest are
    noise. A third one standing clear of that floor means a third species: an
    intermediate, an aggregate, or a baseline that drifted with time rather
    than with temperature.

    The matrix is **not** mean-centred. Centring would remove one rank -- a
    two-state series centres to rank 1 -- so the count would no longer be the
    number of species, which is the thing worth reporting.

    The test is whether the third singular value stands above the noise, not
    whether it is small: it is small in absolute terms either way. The noise
    floor is estimated from the trailing values, which are noise by
    construction, so this needs a decent number of spectra -- at least five,
    and it sharpens with more.

    Returns
    -------
    dict
        ``third_over_noise`` -- the third singular value in units of the noise
        floor, which is the number to read. Around 1 means two states; a value
        of several means a third species. ``is_two_state`` applies
        :data:`THIRD_COMPONENT_THRESHOLD` to it. ``singular_values`` and
        ``noise_floor`` are given so the judgement can be made by eye.

    Examples
    --------
    On two synthetic series of 29 spectra with identical noise, one two-state
    and one passing through an intermediate, this gives ``third_over_noise``
    of about 1.6 and about 20 respectively.

    What makes that worth having is what the fit says meanwhile. On the
    three-state data a two-state curve at 222 nm fits with an unremarkable
    residual and reports **Tm = 61.7 C**, against true transitions at 42 and
    62 C: it finds the second, and the first leaves no trace in Tm at all.

    The fitted dH does carry a signature -- 39 kJ/mol against the 220 that
    built each step. **A van 't Hoff enthalpy far below what the transition
    should have is the classical sign that it is not two-state**, and it is
    worth checking for that reason. But it needs an expectation to compare
    against, which is exactly what an unknown protein does not come with.
    These two checks need no such expectation.
    """
    x, matrix = collection.to_matrix()
    if region is not None:
        low, high = sorted(region)
        inside = (x >= low) & (x <= high)
        matrix = matrix[:, inside]
    if matrix.shape[0] < 5:
        raise ValueError(
            f"estimating the noise floor needs the singular values beyond the "
            f"third, so at least five spectra; got {matrix.shape[0]}"
        )

    values = np.linalg.svd(matrix, compute_uv=False)
    noise_floor = float(np.median(values[3:])) or float(values[-1]) or 1.0
    third = float(values[2])
    return {
        'third_over_noise': third / noise_floor,
        'is_two_state': bool(third / noise_floor
                             < THIRD_COMPONENT_THRESHOLD),
        'noise_floor': noise_floor,
        'singular_values': values[:5].tolist(),
    }
