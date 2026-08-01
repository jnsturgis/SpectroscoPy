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

import numpy as np

__all__ = [
    'convert_x', 'convert_y', 'can_convert_x', 'can_convert_y',
    'X_UNITS', 'Y_UNITS',
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
