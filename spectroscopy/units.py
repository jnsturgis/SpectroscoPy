# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Axis unit conversion.

Roadmap section 2.2 asks for ``spectrum.to("nm")`` as a first-class, tested
operation, on the grounds that nm-vs-cm^-1 and absorbance-vs-%T confusion is
exactly the kind of bug that erodes trust in a shared tool.

Why not pint
------------
The roadmap suggests pint. Implemented natively instead, for two reasons that
only became clear once the conversions were written down:

1. **Half of what is needed is not a unit conversion at all.** Absorbance to
   transmittance is ``T = 10**-A``, a functional transform between two
   dimensionless quantities. pint does not do that; it would need custom code
   regardless.
2. **The x-axis conversions are reciprocal**, not scalar. nm to cm^-1 is
   ``1e7 / nm``. pint handles this only through its ``spectroscopy`` context,
   so the call site is no simpler than the formula.

The unit space here is small and closed, so a native table is exact,
dependency-free, and covers the y axis too. If the unit space grows -- path
lengths, concentrations, per-unit-time intensities -- revisit this: at that
point pint earns its place. See also review section 5.6 on dependency policy.
"""

from __future__ import annotations

import re

import numpy as np

__all__ = [
    'convert_x', 'convert_y', 'can_convert_x', 'can_convert_y',
    'band_direction', 'is_valley_pointing', 'is_extinction', 'is_bipolar',
    'X_UNITS', 'Y_UNITS', 'BAND_DIRECTION', 'EXTINCTION_UNITS',
    'BIPOLAR_UNITS',
]

#: nm per eV, i.e. h*c in eV*nm (CODATA 2018).
_NM_PER_EV = 1239.841984

#: x units, each with converters to and from nanometres as the common base.
X_UNITS = {
    'nm':    (lambda v: v,               lambda v: v),
    'um':    (lambda v: v * 1000.0,      lambda v: v / 1000.0),
    'cm^-1': (lambda v: 1.0e7 / v,       lambda v: 1.0e7 / v),
    'eV':    (lambda v: _NM_PER_EV / v,  lambda v: _NM_PER_EV / v),
}

#: y units, each with converters to and from transmittance (0-1) as the base.
#: 'counts' and 'a.u.' are deliberately absent: they are not convertible.
#:
#: Note that *absorbance* and *absorptance* are different quantities and the
#: distinction matters. Absorbance (optical density) is -log10(T). Absorptance
#: is the fraction of light actually absorbed, 1 - T. They agree only in the
#: weak-absorption limit, where 1 - T ~ 2.303 * A; by A = 0.3 they differ by
#: 15%. Absorptance is the one to compare against an excitation spectrum,
#: because the excitation signal is proportional to photons absorbed, not to
#: optical density.
Y_UNITS = {
    'transmittance': (lambda v: v,                    lambda v: v),
    '%T':            (lambda v: v / 100.0,            lambda v: v * 100.0),
    'absorbance':    (lambda v: np.power(10.0, -v),   lambda v: -np.log10(v)),
    'absorptance':   (lambda v: 1.0 - v,              lambda v: 1.0 - v),
    '%absorptance':  (lambda v: 1.0 - v / 100.0,      lambda v: 100.0 * (1.0 - v)),
}


#: Which way a band points in each y unit: ``'up'`` for a maximum, ``'down'``
#: for a minimum. A unit that is absent is not a statement that bands point
#: up -- it is an absence of knowledge, and :func:`band_direction` reports it
#: as ``'unknown'`` so a caller can decide what to do.
#:
#: This is not a detail of peak picking; it is a property of the quantity, and
#: it belongs with the units for the same reason ``absorbance`` vs ``%T``
#: does. Galactic's SPC format made the same call in 1993 -- ``SPC.H`` splits
#: its y-type codes at 128 with the comment "ALL HIGHER MUST HAVE VALLEYS",
#: because every program that draws or picks a spectrum needs to know.
BAND_DIRECTION = {
    'absorbance':    'up',
    'absorptance':   'up',
    '%absorptance':  'up',
    'transmittance': 'down',
    '%T':            'down',
    'reflectance':   'down',
    # Signed quantities: bands go both ways in the same spectrum, and both
    # signs are signal. See BIPOLAR_UNITS below.
    'mdeg':          'both',
    'deg':           'both',
    'anisotropy':    'both',
    'polarization':  'both',
    'dA':            'both',
    'delta absorbance': 'both',
}

#: The signed y quantities, listed together because they are a class rather
#: than a handful of special cases.
#:
#: A dichroism spectrum is a **difference** between two measurements -- left
#: minus right circular polarisation for CD, parallel minus perpendicular for
#: LD -- so its sign carries meaning that an absorbance never has. An
#: alpha-helix shows negative bands at 208 and 222 nm and a positive one at
#: 193 nm, in one spectrum, and none of the three is a baseline artefact.
#: Fluorescence anisotropy and polarisation are signed for the same reason.
#:
#: ``'unknown'`` would be the wrong answer for these. Unknown means the
#: quantity does not fix a direction; here it fixes it precisely, and the
#: answer is *both*.
#:
#: **A plain difference spectrum cannot be recognised this way**, because
#: subtracting two absorbance spectra leaves the unit saying ``absorbance``.
#: Label it ``'dA'`` to have it treated as signed, or say so at the call site
#: with ``troughs='both'`` -- which is why that spelling exists.
BIPOLAR_UNITS = ('mdeg', 'deg', 'anisotropy', 'polarization',
                 'dA', 'delta absorbance')


#: Extinction-coefficient units, matched exactly. Anything of the form
#: ``<concentration>^-1 cm^-1`` is also recognised by :func:`is_extinction`.
EXTINCTION_UNITS = (
    'M^-1 cm^-1', 'mM^-1 cm^-1', 'uM^-1 cm^-1', 'nM^-1 cm^-1',
    '(mg/mL)^-1 cm^-1', '(ug/mL)^-1 cm^-1', '(mg/mL)^-1', 'M^-1',
)

#: ``<concentration>^-1 cm^-1``. The ``cm^-1`` suffix is required, not
#: optional: without it the pattern also matches ``cm^-1`` itself, which is
#: wavenumber -- an x unit, and emphatically not an extinction coefficient.
#: Path-length-free spellings such as ``M^-1`` are listed explicitly instead.
_EXTINCTION_PATTERN = re.compile(r'^\(?[^()]+\)?\^-1\s+cm\^-1$')


def is_extinction(unit) -> bool:
    """
    True for an extinction-coefficient unit -- ``M^-1 cm^-1`` and friends.

    Deliberately **not** a member of :data:`Y_UNITS`. That table means "can be
    converted to transmittance", and epsilon cannot: going from epsilon to
    absorbance is ``A = eps * c * l``, which needs a concentration and a path
    length that are not properties of the spectrum. It is a different
    quantity, not another spelling of absorbance.

    It is absorbance-*shaped* though -- bands point up -- which is why it
    appears in :data:`BAND_DIRECTION`.
    """
    if unit in EXTINCTION_UNITS:
        return True
    return bool(unit) and bool(_EXTINCTION_PATTERN.match(str(unit).strip()))


def band_direction(unit) -> str:
    """
    ``'up'``, ``'down'``, ``'both'`` or ``'unknown'`` for a y unit.

    ``'both'`` is for the **signed** quantities of :data:`BIPOLAR_UNITS` --
    dichroism, anisotropy, an explicitly labelled difference spectrum -- where
    maxima and minima are equally real and a caller looking for one of them
    will miss half the spectrum.

    ``'unknown'`` covers the honestly ambiguous ones -- ``a.u.``, ``counts``,
    an emission intensity, an *unlabelled* difference spectrum -- where the
    quantity does not fix a direction. Callers should keep their existing
    behaviour there rather than guess, and say what they assumed.

    The distinction between the last two is worth keeping sharp: ``'unknown'``
    means nobody knows, ``'both'`` means the answer is known and is *both*.
    """
    if unit in BAND_DIRECTION:
        return BAND_DIRECTION[unit]
    # An extinction spectrum is absorbance divided by two constants, so its
    # bands point the same way absorbance's do. Recognised by shape rather than
    # by enumeration, because the concentration unit in front varies.
    return 'up' if is_extinction(unit) else 'unknown'


def is_valley_pointing(unit) -> bool:
    """
    True when bands in ``unit`` are minima.

    Worth more than the direction alone: these are also the units in which
    peak *heights*, *areas* and *ratios* are not linear in concentration.
    Beer-Lambert applies to absorbance, so a band ratio taken on ``%T`` is not
    a ratio of anything. Callers doing quantitative work should convert.

    False for the signed units of :data:`BIPOLAR_UNITS`, and correctly so on
    both counts: their bands are not *exclusively* minima, and a dichroism
    signal is linear in concentration in the way ``%T`` is not.
    """
    return band_direction(unit) == 'down'


def is_bipolar(unit) -> bool:
    """True when bands in ``unit`` go both ways and both signs are signal."""
    return band_direction(unit) == 'both'


def can_convert_x(unit) -> bool:
    """Whether ``unit`` is a known, convertible x unit."""
    return unit in X_UNITS


def can_convert_y(unit) -> bool:
    """Whether ``unit`` is a known, convertible y unit."""
    return unit in Y_UNITS


def convert_x(values, from_unit, to_unit):
    """
    Convert x values between units.

    Returns ``(new_values, reversed)``. ``reversed`` is True when the
    conversion inverted the ordering of the axis -- going from nm to cm^-1
    turns an ascending wavelength axis into a descending wavenumber one, and
    callers have to reorder y to match.
    """
    if from_unit == to_unit:
        return np.asarray(values, dtype=float), False
    for unit in (from_unit, to_unit):
        if unit not in X_UNITS:
            raise ValueError(
                f"cannot convert x unit {unit!r}; known units are "
                f"{', '.join(sorted(X_UNITS))}"
            )

    values = np.asarray(values, dtype=float)
    if np.any(values == 0):
        raise ValueError(
            f"cannot convert {from_unit!r} to {to_unit!r}: the axis contains "
            f"zero, and the conversion is reciprocal"
        )

    to_nm, _ = X_UNITS[from_unit]
    _, from_nm = X_UNITS[to_unit]
    converted = from_nm(to_nm(values))

    # An ordering flip is what a reciprocal conversion does; detect it from the
    # data rather than from a table of which pairs are reciprocal.
    was_ascending = values[0] < values[-1] if len(values) > 1 else True
    now_ascending = converted[0] < converted[-1] if len(converted) > 1 else True
    return converted, was_ascending != now_ascending


def convert_y(values, from_unit, to_unit):
    """
    Convert y values between units.

    Absorbance/transmittance is a functional transform, not a rescaling:
    ``A = -log10(T)``. Converting a spectrum with non-positive transmittance to
    absorbance is undefined and raises rather than silently producing inf.
    """
    if from_unit == to_unit:
        return np.asarray(values, dtype=float)
    for unit in (from_unit, to_unit):
        if unit not in Y_UNITS:
            raise ValueError(
                f"cannot convert y unit {unit!r}; convertible units are "
                f"{', '.join(sorted(Y_UNITS))}. 'counts' and 'a.u.' are "
                f"instrument-specific and have no defined conversion."
            )

    values = np.asarray(values, dtype=float)
    to_base, _ = Y_UNITS[from_unit]
    _, from_base = Y_UNITS[to_unit]

    base = to_base(values)
    if to_unit == 'absorbance' and np.any(base <= 0):
        raise ValueError(
            f"cannot convert {from_unit!r} to 'absorbance': the spectrum "
            f"contains non-positive transmittance, and log10 of that is "
            f"undefined. Baseline-correct or crop first."
        )
    return from_base(base)
