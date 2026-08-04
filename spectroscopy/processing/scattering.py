# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Removing a scattering background from an absorbance spectrum.

A turbid sample -- membranes, vesicles, inclusion bodies, anything that has
not fully dissolved -- attenuates light by scattering it as well as by
absorbing it, and a spectrophotometer cannot tell the difference. What it
reports as absorbance is the sum, and the scattering part rises steeply toward
the blue, so it sits underneath exactly the bands UV-Vis work depends on.

This belongs **before** unmixing or any concentration estimate, not after. A
scattering background is not one of the components, and
:func:`~spectroscopy.processing.unmix.unmix` will otherwise spread it across
whichever references happen to slope the same way -- silently, and with an
excellent R^2.

The model
---------
Scattering by particles small compared with the wavelength goes as
``lambda**-4`` (Rayleigh). Real samples are not that: as particles approach
the wavelength the exponent falls, through the Mie regime, toward 0 for
particles much larger than lambda. A single fixed exponent is therefore a
guess about particle size, and a bad one biases the whole correction.

So the default is a **basis of power laws** -- several exponents fitted
together with non-negative coefficients -- rather than one. That spans the
range of behaviour a real polydisperse sample shows, without anybody having to
declare a particle size, and it collapses to the single-exponent case when
that is genuinely what the data wants.

Measured backgrounds can be used instead, and are better when available: a
blank of the same particles without the chromophore is the honest reference.
See :func:`from_references`.

The fit region
--------------
Fitted where **nothing absorbs**, then extrapolated across the whole spectrum.
That is what makes it a correction rather than a curve subtraction: choosing a
region that contains a band lets the fit eat the band. For protein and nucleic
acid the usual window is 320-400 nm, which is why that is the default here,
but it is a default and not a law -- a sample with a haem, a carotenoid or any
other visible chromophore needs somewhere else.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import nnls

__all__ = ['scatter_baseline', 'correct_scattering', 'from_references',
           'DEFAULT_EXPONENTS', 'DEFAULT_REGION']

#: Exponents in the default basis. Spans Rayleigh (4) through the Mie regime to
#: nearly flat (0.5), which is what a polydisperse suspension actually shows.
DEFAULT_EXPONENTS = (4.0, 3.0, 2.0, 1.0, 0.5)

#: Where protein and nucleic acid stop absorbing, in nm. A default, not a law.
DEFAULT_REGION = (320.0, 400.0)


def _basis(x, exponents, reference_x):
    """Power-law columns, each normalised at ``reference_x`` so the fitted
    coefficients are comparable and the design matrix stays well conditioned."""
    x = np.asarray(x, dtype=float)
    return np.vstack([(x / reference_x) ** (-power) for power in exponents]).T


def scatter_baseline(x, y, *, region=DEFAULT_REGION, exponents=DEFAULT_EXPONENTS,
                     non_negative=True):
    """
    Fit a scattering background and return it across the whole axis.

    Parameters
    ----------
    x, y : array_like
        Wavelength in nm and absorbance.
    region : (float, float)
        Where to fit -- a window in which the sample is believed not to
        absorb. Order-insensitive.
    exponents : sequence of float
        Powers of ``1/lambda`` to fit. A single value gives the classical
        fixed-exponent correction; the default basis spans Rayleigh to nearly
        flat.
    non_negative : bool
        Keep the coefficients at or above zero. Scattering adds signal, it
        does not remove it, so a negative contribution is unphysical.

    Returns
    -------
    ndarray
        The background, on the same axis as ``x``.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    exponents = np.atleast_1d(np.asarray(exponents, dtype=float))

    low, high = sorted(region)
    inside = (x >= low) & (x <= high)
    if inside.sum() < max(len(exponents) + 1, 3):
        raise ValueError(
            f"the fit region {low:g}-{high:g} holds {int(inside.sum())} points, "
            f"which is not enough to fit {len(exponents)} components. Widen it, "
            f"or pass fewer exponents."
        )

    # Normalising at the red end keeps the columns O(1) over the fitted window
    # and makes the extrapolation to the blue the only place they grow.
    reference_x = float(np.max(x))
    design = _basis(x[inside], exponents, reference_x)

    if non_negative:
        coefficients, _ = nnls(design, y[inside])
    else:
        coefficients, *_ = np.linalg.lstsq(design, y[inside], rcond=None)

    return _basis(x, exponents, reference_x) @ coefficients


def correct_scattering(spectrum, *, region=DEFAULT_REGION,
                       exponents=DEFAULT_EXPONENTS, non_negative=True,
                       return_background=False):
    """
    Subtract a fitted scattering background from a spectrum.

    Returns the corrected :class:`~spectroscopy.spectra.Spectrum`, or
    ``(corrected, background)`` with ``return_background=True`` -- and looking
    at the background is worth the extra line, for the same reason
    :func:`~spectroscopy.viz.plot_baseline` exists.

    Warns when the correction leaves the fit region substantially negative,
    which means the background was fitted through something that absorbs.
    """
    from spectroscopy.history import ProcessingStep  # noqa: PLC0415

    x = np.asarray(spectrum.x, dtype=float)
    background = scatter_baseline(x, spectrum.y, region=region,
                                  exponents=exponents,
                                  non_negative=non_negative)

    corrected = spectrum._derive(
        y=np.asarray(spectrum.y, dtype=float) - background,
        step=ProcessingStep('correct_scattering', {
            'region': tuple(sorted(region)),
            'exponents': [float(e) for e in np.atleast_1d(exponents)],
            'non_negative': non_negative}))

    # In a region where nothing absorbs, subtracting the background should
    # leave approximately zero. Anything else means the region was not
    # absorption-free and the fit has eaten part of a band -- which is the one
    # way this correction goes badly wrong, and it does not announce itself.
    low, high = sorted(region)
    inside = (x >= low) & (x <= high)
    residual = corrected.y[inside]
    span = float(np.nanmax(corrected.y) - np.nanmin(corrected.y))
    if len(residual) and span > 0:
        leftover = float(np.nanmax(np.abs(residual))) / span
        if leftover > 0.1:
            warnings.warn(
                f"after correction the fit region {low:g}-{high:g} nm still "
                f"holds structure at {100 * leftover:.0f}% of the spectrum's "
                f"range, so it is probably not absorption-free and the "
                f"background has absorbed part of a band. Choose a region "
                f"where your sample does not absorb.",
                UserWarning, stacklevel=2,
            )

    if return_background:
        from spectroscopy.spectra import Spectrum  # noqa: PLC0415
        background_spectrum = Spectrum(
            x, background, x_quantity=spectrum.x_quantity,
            x_unit=spectrum.x_unit, y_quantity='Scattering',
            y_unit=spectrum.y_unit)
        background_spectrum.name = f"{spectrum.name} (scattering)"
        return corrected, background_spectrum
    return corrected


def from_references(spectrum, backgrounds, *, region=DEFAULT_REGION,
                    return_background=False):
    """
    Correct using **measured** scattering backgrounds rather than power laws.

    The better method where the blanks exist: a suspension of the same
    particles without the chromophore scatters the way the sample scatters,
    including whatever the power-law basis cannot express. Several blanks
    spanning a range of turbidity are interpolated between -- fitted, strictly
    -- so the amount of scattering does not have to be matched in advance.

    Parameters
    ----------
    backgrounds : sequence of Spectrum
        Measured scattering-only spectra. Resampled onto ``spectrum``.
    region : (float, float)
        Where to fit them, as for :func:`scatter_baseline`.
    """
    from spectroscopy.history import ProcessingStep  # noqa: PLC0415
    from spectroscopy.spectra import Spectrum  # noqa: PLC0415

    backgrounds = list(backgrounds)
    if not backgrounds:
        raise ValueError("no scattering backgrounds given")

    x = np.asarray(spectrum.x, dtype=float)
    matrix = np.vstack([b.resample(x).y for b in backgrounds]).T

    low, high = sorted(region)
    inside = (x >= low) & (x <= high)
    if inside.sum() < len(backgrounds) + 1:
        raise ValueError(
            f"the fit region {low:g}-{high:g} holds {int(inside.sum())} points "
            f"for {len(backgrounds)} backgrounds; widen it or use fewer"
        )

    coefficients, _ = nnls(matrix[inside], np.asarray(spectrum.y)[inside])
    background = matrix @ coefficients

    corrected = spectrum._derive(
        y=np.asarray(spectrum.y, dtype=float) - background,
        step=ProcessingStep('correct_scattering', {
            'method': 'measured backgrounds',
            'n_backgrounds': len(backgrounds),
            'weights': [round(float(c), 6) for c in coefficients],
            'region': (low, high)}))

    if return_background:
        background_spectrum = Spectrum(
            x, background, x_quantity=spectrum.x_quantity,
            x_unit=spectrum.x_unit, y_quantity='Scattering',
            y_unit=spectrum.y_unit)
        background_spectrum.name = f"{spectrum.name} (scattering)"
        return corrected, background_spectrum
    return corrected
