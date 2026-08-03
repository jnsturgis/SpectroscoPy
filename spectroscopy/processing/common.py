# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Technique-agnostic algorithms, as plain functions over numpy arrays.

Nothing here knows about :class:`~spectroscopy.spectra.Spectrum`: these take
``x`` and ``y`` and return arrays. The Spectrum methods are thin wrappers that
call these and record a :class:`~spectroscopy.history.ProcessingStep`. Keeping
the split means the numerics can be tested against synthetic ground truth
without constructing objects, and that ``processing`` never imports ``core``.

Priorities follow the notebook inventory (review section 1.3): rubberband
baseline and second-derivative peak detection first, because between them they
were copy-pasted into a dozen notebooks.
"""

from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks as _find_peaks
from scipy.signal import savgol_filter
from scipy.spatial import ConvexHull, QhullError

from spectroscopy import units

__all__ = [
    'als_baseline', 'poly_baseline', 'rubberband_baseline',
    'baseline', 'normalize', 'smooth', 'derivative',
    'detect_peaks', 'window_mask',
]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def window_mask(x, window=None):
    """
    Boolean mask selecting ``window`` = (low, high) on ``x``.

    The bounds are order-insensitive, which matters because FTIR axes are
    conventionally quoted high-to-low: (1800, 900) and (900, 1800) both select
    the amide region.
    """
    if window is None:
        return np.ones(len(x), dtype=bool)
    low, high = min(window), max(window)
    return (x >= low) & (x <= high)


def _as_odd_at_least(value, minimum):
    """Savitzky-Golay needs an odd window strictly greater than the order."""
    value = int(value)
    if value <= minimum:
        value = minimum + 1
    if value % 2 == 0:
        value += 1
    return value


# ---------------------------------------------------------------------------
# baselines
# ---------------------------------------------------------------------------

def rubberband_baseline(x, y, return_points=False, side='lower'):
    """
    Convex-hull ("rubberband") baseline: one hull arc, interpolated.

    This is the function that was pasted verbatim into twelve notebooks.

    Parameters
    ----------
    side : {'lower', 'upper'}
        Which envelope to take. **This depends on which way the bands point.**
        In absorbance the baseline is the *lower* hull, tending to zero; in
        transmittance it is the *upper* hull, tending to 100%. Taking the lower
        hull of a transmittance spectrum traces the tips of the absorption
        bands and is meaningless.

        :meth:`~spectroscopy.spectra.Spectrum.baseline` picks this from the
        spectrum's y unit, so it is usually not something to think about.
    return_points : bool
        Also return the (x, y) anchor points the baseline passes through, so
        they can be inspected, plotted, adjusted by eye, and fed back to
        :func:`poly_baseline` as explicit guide points.
    """
    if side not in ('lower', 'upper'):
        raise ValueError(f"side must be 'lower' or 'upper', not {side!r}")
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        edge = (y.min() if side == 'lower' else y.max()) if len(y) else 0.0
        base = np.full_like(y, edge)
        return (base, x.copy(), y.copy()) if return_points else base

    try:
        hull = ConvexHull(np.column_stack((x, y)))
    except QhullError:
        # Degenerate input: the points are collinear (a flat or perfectly
        # linear spectrum), so there is no 2-D hull. The rubberband of a
        # straight line is that line, which the endpoint interpolation below
        # reproduces exactly. Without this the user gets a wall of qhull
        # diagnostics instead of a baseline.
        ends = np.array([np.argmin(x), np.argmax(x)])
        anchor_x, anchor_y = x[ends], y[ends]
        base = np.interp(x, anchor_x, anchor_y)
        return (base, anchor_x, anchor_y) if return_points else base

    vertices = hull.vertices                       # counter-clockwise

    left = np.where(vertices == np.argmin(x))[0][0]
    right = np.where(vertices == np.argmax(x))[0][0]

    # The hull is a closed loop; the two arcs between leftmost and rightmost
    # point are the upper and lower envelopes.
    if left < right:
        arc_a = vertices[left:right + 1]
        arc_b = np.concatenate((vertices[right:], vertices[:left + 1]))
    else:
        arc_a = vertices[right:left + 1]
        arc_b = np.concatenate((vertices[left:], vertices[:right + 1]))

    if side == 'lower':
        chosen = arc_a if np.mean(y[arc_a]) < np.mean(y[arc_b]) else arc_b
    else:
        chosen = arc_a if np.mean(y[arc_a]) > np.mean(y[arc_b]) else arc_b

    chosen = chosen[np.argsort(x[chosen])]
    anchor_x, keep = np.unique(x[chosen], return_index=True)
    anchor_y = y[chosen][keep]

    base = np.interp(x, anchor_x, anchor_y)
    return (base, anchor_x, anchor_y) if return_points else base


def poly_baseline(x, y, degree=3, points=None, coefficients=None, halfwidth=0):
    """
    Polynomial baseline, from guide points or from known coefficients.

    Guide points are the intended everyday form (review section 5.4): give the
    x positions where the spectrum is known to be baseline, and a polynomial of
    ``degree`` is least-squares fitted through the spectrum's own y values
    there.

    Parameters
    ----------
    points : sequence of float, optional
        Guide-point x positions. The nearest sample to each is used; with
        ``halfwidth > 0`` the median over +/- that many samples is used
        instead, which is steadier on noisy data.
    coefficients : sequence of float, optional
        Use these polynomial coefficients directly instead of fitting, in
        numpy's highest-power-first order, as ``np.polyval`` expects.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if coefficients is not None:
        return np.polyval(np.asarray(coefficients, dtype=float), x)

    if points is None:
        raise ValueError(
            "poly_baseline needs either 'points' (guide points to fit through) "
            "or 'coefficients' (a known polynomial)"
        )

    points = np.atleast_1d(np.asarray(points, dtype=float))
    if len(points) <= degree:
        raise ValueError(
            f"need more than degree ({degree}) guide points to fit, got {len(points)}"
        )

    indices = np.abs(x[None, :] - points[:, None]).argmin(axis=1)
    if halfwidth > 0:
        guide_y = np.array([
            np.median(y[max(0, i - halfwidth):i + halfwidth + 1]) for i in indices
        ])
    else:
        guide_y = y[indices]

    return np.polyval(np.polyfit(x[indices], guide_y, degree), x)


def als_baseline(y, lam=1e6, p=0.01, niter=10):
    """
    Asymmetric least squares baseline (Eilers & Boelens).

    Points above the current baseline get weight ``p``, points below get
    ``1 - p``, so the fit is pulled under the peaks; ``lam`` sets smoothness.

    Note: the copy of this in the notebooks never ran -- it referenced
    ``sparse`` and ``spsolve`` without importing them, and was never called.
    This is the first working version.
    """
    from scipy import sparse  # pylint: disable=C0415
    from scipy.sparse.linalg import spsolve  # pylint: disable=C0415

    y = np.asarray(y, dtype=float)
    size = len(y)
    if size < 3:
        return y.copy()

    differences = sparse.diags([1, -2, 1], [0, 1, 2], shape=(size - 2, size))
    smoothness = lam * differences.T.dot(differences)

    weights = np.ones(size)
    base = y.copy()
    for _ in range(int(niter)):
        base = spsolve(sparse.diags(weights, 0) + smoothness, weights * y)
        weights = p * (y > base) + (1 - p) * (y <= base)
    return base


#: Units in which bands point *down* from a high baseline, so the rubberband
#: envelope is the upper hull. Everything else is treated as absorbance-like.
INVERTED_Y_UNITS = ('transmittance', '%t', 'reflectance')


def baseline_side_for(y_unit):
    """
    Which hull arc is the baseline for data in ``y_unit``.

    Note the entries are compared lower-cased, so '%T' has to be stored as
    '%t' -- an easy way to write a check that silently never matches.
    """
    return 'upper' if str(y_unit).lower() in INVERTED_Y_UNITS else 'lower'


#: Accepted spellings for each baseline method, so that both the legacy
#: uppercase codes and readable names work.
BASELINE_ALIASES = {
    'rb': 'rubberband', 'rubberband': 'rubberband', 'hull': 'rubberband',
    'poly': 'poly', 'polynomial': 'poly',
    'als': 'als', 'asls': 'als',
}


def baseline(x, y, method='rubberband', *, return_points=False, **kwargs):
    """Dispatch to a baseline algorithm by name. See the specific functions."""
    canonical = BASELINE_ALIASES.get(str(method).lower())
    if canonical is None:
        raise ValueError(
            f"Unknown baseline method {method!r}; "
            f"try one of {sorted(set(BASELINE_ALIASES.values()))}"
        )

    if canonical == 'rubberband':
        return rubberband_baseline(x, y, return_points=return_points,
                                   side=kwargs.pop('side', 'lower'))

    if canonical == 'poly':
        base = poly_baseline(x, y, **kwargs)
    else:
        base = als_baseline(y, **kwargs)

    if return_points:
        return base, np.asarray([]), np.asarray([])
    return base


# ---------------------------------------------------------------------------
# smoothing, derivatives, normalisation
# ---------------------------------------------------------------------------

def smooth(y, method='savgol', window_length=15, polyorder=3):
    """Smooth ``y``. Methods: ``savgol`` (default) or ``moving_average``."""
    y = np.asarray(y, dtype=float)
    name = str(method).lower()

    if name in ('savgol', 'sg', 'savitzky-golay', 'savitsky-golay'):
        return savgol_filter(y, _as_odd_at_least(window_length, polyorder), polyorder)
    if name in ('moving_average', 'boxcar', 'ma'):
        width = max(1, int(window_length))
        kernel = np.ones(width) / width
        return np.convolve(y, kernel, mode='same')
    raise ValueError(f"Unknown smoothing method {method!r}; try 'savgol'.")


def derivative(y, order=1, window_length=11, polyorder=3, delta=1.0):
    """
    Savitzky-Golay derivative of ``y``.

    The second derivative is the workhorse here: inverted, its peaks mark band
    positions including shoulders, which is how every peak-picking cell in the
    notebooks works.
    """
    y = np.asarray(y, dtype=float)
    window = _as_odd_at_least(window_length, polyorder)
    return savgol_filter(y, window, polyorder, deriv=int(order), delta=delta)


def normalize(x, y, method='max', window=None):
    """
    Scale ``y``. Methods:

    ``max``      divide by the maximum (optionally within ``window``) -- the
                 notebooks' usual "normalise to the 1050-1080 glycan band".
    ``area``     divide by the trapezoidal area within ``window``.
    ``minmax``   rescale to [0, 1].
    ``vector``   divide by the Euclidean norm.
    ``snv``      standard normal variate: subtract mean, divide by std.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = window_mask(x, window)
    if not mask.any():
        raise ValueError(f"normalisation window {window} selects no points")

    name = str(method).lower()
    if name == 'max':
        factor = np.nanmax(y[mask])
    elif name == 'area':
        order = np.argsort(x[mask])
        factor = np.trapezoid(y[mask][order], x[mask][order])
    elif name == 'minmax':
        low, high = np.nanmin(y[mask]), np.nanmax(y[mask])
        if high == low:
            raise ValueError("cannot min-max normalise a flat spectrum")
        return (y - low) / (high - low)
    elif name == 'vector':
        factor = np.linalg.norm(y[mask])
    elif name == 'snv':
        spread = np.nanstd(y[mask])
        if spread == 0:
            raise ValueError("cannot SNV normalise a flat spectrum")
        return (y - np.nanmean(y[mask])) / spread
    else:
        raise ValueError(f"Unknown normalisation method {method!r}.")

    if factor == 0 or not np.isfinite(factor):
        raise ValueError(f"normalisation factor is {factor}; cannot scale")
    return y / factor


# ---------------------------------------------------------------------------
# peak detection
# ---------------------------------------------------------------------------

#: find_peaks arguments that are a y-magnitude, and so need rescaling when
#: thresholds are given relative to the signal.
_MAGNITUDE_ARGUMENTS = ('height', 'prominence', 'threshold')


def detect_peaks(x, y, method='second_derivative', *, troughs=None,
                 window_length=11, polyorder=3, relative=False, y_unit=None,
                 **kwargs):
    """
    Locate peaks and return ``(indices, properties)``.

    ``method='second_derivative'`` (the default, and what every notebook does)
    finds peaks in the inverted second derivative, so shoulders on a broad band
    are found as well as maxima. ``method='direct'`` runs scipy's find_peaks on
    the signal itself.

    Extra keyword arguments go straight to :func:`scipy.signal.find_peaks`
    (``height``, ``distance``, ``prominence``, ``width`` ...).

    Parameters
    ----------
    troughs : bool, optional
        Look for minima instead of maxima. **Left unset, the direction is
        taken from** ``y_unit``: bands in transmittance, ``%T`` and reflectance
        are minima, so those are searched downward, and absorbance-like units
        upward. Pass it explicitly to override -- which is what a difference
        spectrum needs, where both directions are meaningful and the unit
        cannot tell you which one you meant.

        This defaulting exists because getting it wrong is silent. On a
        transmission spectrum, searching upward returns the two inflection
        points flanking each band and not the band, so a single band at
        1100 cm-1 is reported as a pair at 1086 and 1114 -- a plausible number
        of plausible-looking positions, none of them a band.
    y_unit : str, optional
        The y unit of ``y``, used only to choose the direction when ``troughs``
        is unset. Unknown or absent units keep the historical upward
        behaviour, which is right for the a.u. and counts cases that dominate
        Raman and fluorescence; the assumption is reported in the returned
        properties rather than guessed at silently.
    relative : bool
        Interpret ``height``, ``prominence`` and ``threshold`` as fractions of
        the detection signal's range rather than as absolute values.

        Worth knowing about: with the default method those thresholds apply to
        the *second derivative*, whose magnitude is typically four orders
        smaller than the spectrum. On a spectrum normalised to 1.0 the second
        derivative spans about +/-0.005, so ``prominence=0.02`` -- a reasonable
        guess in spectrum units -- silently finds nothing, while the value that
        works is ``0.0001``. That is where the unexplained small constants in
        the notebooks come from. With ``relative=True``, ``prominence=0.02``
        means 2% of the detection signal's range whichever method is in use.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    name = str(method).lower()
    if name in ('second_derivative', 'd2', '2nd'):
        signal = -derivative(y, order=2, window_length=window_length,
                             polyorder=polyorder)
    elif name == 'direct':
        signal = y
    else:
        raise ValueError(
            f"Unknown peak method {method!r}; try 'second_derivative' or 'direct'."
        )

    # Direction: an explicit troughs= always wins, because a difference
    # spectrum has bands both ways and only the caller knows which is meant.
    direction = units.band_direction(y_unit)
    if troughs is None:
        troughs = direction == 'down'
        chosen_by = 'y_unit' if direction != 'unknown' else 'assumed'
    else:
        troughs = bool(troughs)
        chosen_by = 'caller'

    if troughs:
        signal = -signal

    if relative:
        span = float(np.nanmax(signal) - np.nanmin(signal))
        kwargs = {key: (np.asarray(value) * span
                        if key in _MAGNITUDE_ARGUMENTS and value is not None
                        else value)
                  for key, value in kwargs.items()}

    indices, properties = _find_peaks(signal, **kwargs)
    # Recorded, not just acted on: which way the search ran and why is part of
    # reading the result, and it is the first thing to check when a peak list
    # looks wrong.
    properties = {**properties, 'troughs': troughs,
                  'direction_from': chosen_by, 'y_unit': y_unit}
    return indices, properties
