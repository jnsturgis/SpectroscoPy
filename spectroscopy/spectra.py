# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Spectra structure to hold a spectrum and associated metadata.

The spectrum structure
The metadata associated with the spectrum and default values:

name	 a name for the spectrum initially derived from the filename.
filename the name of the file the spectrum came from.
x_labels a list of 2 tuples giving (label, units) for each x dimension.
y_label  a 2 tuple describing the y data (name, units)

x        a numpy array of x values ordered along all axes
y        a numpy array of y values at positions given by the x values

         Renamed from x_data/y_data in 0.1 (roadmap section 2.5). They are
         plain attributes for now; Phase 1 makes them properties so that the
         backing store can change without breaking callers.

Acquisition metadata - for original untreated data it is interesting to have
information on the acquisition parameters.

Treatment metadata - is it interesting to have the history? Certainly not
for the moment.

Todo:
    Should define functions for jcamp required info and use it.
    Should make sanity checks and die gracefully

-------------------------------------------------------------

"""

# pylint: disable=W0511, W0107
import copy
import operator
import os
import warnings
from typing import TYPE_CHECKING

import numpy as np
from scipy.interpolate import CubicSpline

import spectroscopy.messages
from spectroscopy import units
from spectroscopy.history import ProcessingStep, describe_operand
from spectroscopy.io import registry
from spectroscopy.peaks import PeakTable
from spectroscopy.processing import common

# Annotation only: fit_peaks imports it lazily, which keeps the two modules
# acyclic (fitting imports lineshapes, and lineshapes imports nothing of ours).
if TYPE_CHECKING:
    from spectroscopy.fitting import FitResult

#: The one name this module exports. Everything below -- the lookup tables, the
#: label helper, the imports above -- is implementation. Without this list,
#: ``from .spectra import *`` in the package __init__ re-exported all of it,
#: including ``os``, ``copy``, ``np`` and ``CubicSpline``, which made
#: ``spectroscopy.np`` a public name the 1.0 freeze would have promised to keep.
#: See roadmap section 14.2.
__all__ = ['Spectrum']


#: Extension -> file type, and the set of known types. Both are now *derived*
#: from the format registry rather than maintained here. There used to be four
#: hand-kept tables (these two plus the match statements in reload() and
#: save()); they drifted apart, which is what defect D5 was. Registering a
#: reader is now the only step needed to teach Spectrum a new format.
def _file_exts():
    return registry.known_extensions()


def _known_file_types():
    return registry.known_types()


#: Kept as module attributes for backwards compatibility -- some notebooks
#: inspect them. They are snapshots; call the registry for the live view.
FILE_EXTS = registry.known_extensions()
KNOWNFILETYPES = registry.known_types()
#: Rendering of each unit for a matplotlib axis label. An empty string means
#: dimensionless, so the label is just the quantity with no parentheses.
UNIT_LABELS = {
    'cm^-1':         r'cm$^{-1}$',
    'nm':            'nm',
    'um':            r'$\mu$m',
    'eV':            'eV',
    'absorbance':    '',
    'transmittance': '',
    'absorptance':   '',
    '%absorptance':  '%',
    '%T':            '%',
    'counts':        'counts',
    'a.u.':          'a.u.',
    '':              '',
}

#: Per-technique axis defaults. Quantity and unit are stored separately so the
#: unit is machine readable -- 'cm^-1' vs 'nm' is what makes a future .to("nm")
#: possible, and it is what stops nm/cm^-1 confusion silently producing a
#: plausible-looking wrong plot (roadmap section 2.2).
KNOWNSPECTYPES = {
    'FTIR':         {'x_quantity': 'Wavenumber',   'x_unit': 'cm^-1',
                     'y_quantity': 'Absorbance',   'y_unit': 'absorbance'},
    'ATR-FTIR':     {'x_quantity': 'Wavenumber',   'x_unit': 'cm^-1',
                     'y_quantity': 'Absorbance',   'y_unit': 'absorbance'},
    'UV-Vis':       {'x_quantity': 'Wavelength',   'x_unit': 'nm',
                     'y_quantity': 'Absorbance',   'y_unit': 'absorbance'},
    'Raman':        {'x_quantity': 'Raman shift',  'x_unit': 'cm^-1',
                     'y_quantity': 'Signal',       'y_unit': 'a.u.'},
    'Fluorescence': {'x_quantity': 'Wavelength',   'x_unit': 'nm',
                     'y_quantity': 'Fluorescence', 'y_unit': 'a.u.'},
}

#: Techniques whose x axis is conventionally plotted high-to-low.
REVERSED_AXIS_TECHNIQUES = ('FTIR', 'ATR-FTIR', 'Raman')

#: The quantity a unit implies, used to relabel an axis after conversion so a
#: wavelength axis taken to cm^-1 does not stay called 'Wavelength'.
QUANTITY_FOR_X_UNIT = {
    'nm': 'Wavelength', 'um': 'Wavelength',
    'cm^-1': 'Wavenumber', 'eV': 'Photon energy',
}
QUANTITY_FOR_Y_UNIT = {
    'absorbance': 'Absorbance',
    'transmittance': 'Transmittance',
    '%T': 'Transmittance',
    'absorptance': 'Absorptance',
    '%absorptance': 'Absorptance',
}


def _looks_like_data(args):
    """True when ``args`` is the ``Spectrum(x, y, ...)`` form.

    Exactly two positional arguments, both sized, neither a string, a path or
    a Spectrum. A Spectrum has ``__len__`` too, which is why it is excluded
    explicitly rather than by luck of branch ordering.
    """
    if len(args) != 2:
        return False
    return all(
        not isinstance(a, (str, bytes, os.PathLike, Spectrum))
        and hasattr(a, '__len__')
        for a in args
    )


def compose_label(quantity, unit):
    """Build an axis label from a quantity and a unit, e.g. 'Wavenumber (cm^-1)'."""
    rendered = UNIT_LABELS.get(unit, unit)
    return f"{quantity} ({rendered})" if rendered else str(quantity)

def _infer_file_type( name ):
    """Work out if possible from the file extension the possible file types."""
    return registry.infer_file_type(name)

class Spectrum:
    """A class for spectra."""

    def __init__(self, *args, **kwargs):
        """Initialize a Spectrum object.

        Several forms, distinguished by what is passed:

        ``Spectrum()``
            An empty spectrum.

        ``Spectrum(x, y, **options)``
            **From data.** ``x`` and ``y`` are sequences of equal length --
            arrays, lists, anything :func:`numpy.asarray` accepts. See below
            for the options.

        ``Spectrum(path)``
            Read a file, inferring the format from its extension. Prefer
            :meth:`read`, which is explicit about doing I/O and takes the
            format as a named argument.

        ``Spectrum(path, name)``
            Read ``name`` from directory ``path``, inferring the format.

            .. warning::

               These two arguments are a **directory and a filename**, not a
               filename and a format. Until 0.1.1 this docstring claimed the
               latter, and ``Spectrum("ethanol.jdx", "jcamp")`` raised
               ``TypeError: Unknown filetype unknown`` as a result. Use
               ``Spectrum.read("ethanol.jdx", "jcamp")`` for that.

        ``Spectrum(path, name, file_type)``
            As above with the format given explicitly.

        ``Spectrum(other)``
            A deep copy of another spectrum.

        Parameters
        ----------
        x, y : array_like
            The axis and the intensities, of equal length. Data form only.
        technique : str, optional
            One of the known techniques (``"ATR-FTIR"``, ``"Raman"`` ...).
            Sets the conventional axis quantities and units for it.
        x_quantity, x_unit, y_quantity, y_unit : str, optional
            Set explicitly, and take precedence over ``technique``. Passing
            any of them marks the units authoritative, so a later
            :meth:`set_type` will not quietly relabel them.
        name : str, optional
            A label for the spectrum. Defaults to ``"unnamed"``.
        metadata : dict, optional
            Copied, not referenced.
        history : list of ProcessingStep, optional
            For rebuilding a spectrum whose provenance is already known -- a
            reader, or a round-trip through a file. Building from arrays does
            not itself record a step: creating data is not processing it.

        Raises
        ------
        ValueError
            If ``x`` and ``y`` differ in length, or either is not 1-D.
        TypeError
            If the arguments match none of the forms above.

        Examples
        --------
        >>> import numpy as np
        >>> x = np.linspace(900, 1800, 901)
        >>> s = Spectrum(x, np.exp(-((x - 1650) / 40) ** 2),
        ...              technique="ATR-FTIR", name="synthetic")
        >>> len(s)
        901
        """
        if kwargs and not _looks_like_data(args):
            raise TypeError(
                "keyword arguments are only accepted by the data form, "
                "Spectrum(x, y, ...); got "
                f"{', '.join(sorted(kwargs))}"
            )

        if len(args) == 0 :
            self._init_blank()

        elif _looks_like_data(args):
            self._init_blank()
            self._init_from_data(args[0], args[1], **kwargs)

        elif isinstance(args[0], Spectrum ) :
            # This is the copy method - initiate from another Spectrum()
            other = args[0]
            self.name        = str(other.name)
            self.fileinfo    = {'PATH':'','NAME':'','TYPE':'csv'}
            self._set_axes(other.x_quantity, other.x_unit,
                           other.y_quantity, other.y_unit)
            self._x_label_override = other._x_label_override
            self._y_label_override = other._y_label_override
            self.technique   = other.technique
            self.units_from_file = other.units_from_file
            self.history     = list(other.history)
            self.x      = np.copy(other.x)
            self.y      = np.copy(other.y)
            # deepcopy, not a bare assignment: this constructor is what every
            # arithmetic operator uses to build its result, so sharing the dict
            # meant `avg = (a+b)/2; avg.metadata['sample'] = "avg"` silently
            # renamed a as well (defect D2). Deep rather than shallow because
            # values can be containers -- e.g. metadata['file_header'] from the
            # .dpt reader -- and because this branch is documented above as
            # "a deepcopy of the original".
            self.metadata    = copy.deepcopy(other.metadata)

        elif isinstance(args[0], (str, os.PathLike)) :
            # This is the open method - initialize from a file
            # First parse the args to get filetype and full filename
            # Len 1: filename - infer type from 1
            # Len 2: path, filename - infer type  from 2
            # Len 3: path, filename, type
            args = tuple(os.fspath(a) if isinstance(a, os.PathLike) else a
                         for a in args)

            self._init_blank()
            self.name        = args[0]

            if len(args) == 3:
                self.fileinfo['PATH'] = args[0]
                self.fileinfo['NAME'] = args[1]
                self.fileinfo['TYPE'] = args[2]
                self.name             = args[1]
            elif len(args) == 2:
                self.fileinfo['PATH'] = args[0]
                self.fileinfo['NAME'] = args[1]
                self.fileinfo['TYPE'] = _infer_file_type(args[0]+args[1])
                self.name             = args[1]
            elif len(args) == 1:
                self.fileinfo['NAME'] = args[0]
                self.fileinfo['TYPE'] = _infer_file_type(args[0])

            # Now read the file according to type. Ask the registry, do not
            # consult KNOWNFILETYPES: that is a snapshot taken when this module
            # was imported, so a format registered afterwards -- which is the
            # whole point of register_reader being a decorator anyone can use --
            # would be inferred correctly and then rejected here.
            if self.fileinfo['TYPE'] in _known_file_types():
                self.reload()
            else:
                raise TypeError(
                    f"Unknown filetype {self.fileinfo['TYPE']!r} for "
                    f"{self.fileinfo['NAME']}; known types are "
                    f"{', '.join(_known_file_types())}"
                )
        else :
            raise TypeError("Unknown spectrum initializer")
        # TODO Should finish by checking that we have a well formed Spectrum()
        pass

##=============================================================================
#
#   Some Magic Dunder methods
#
#   Not comparison operations (does not make sense)
#   Not floordiv, mod or power operators
#
#   TODO add an iterator method for point in spectrum.
#
##=============================================================================

    def _init_blank(self) -> None:
        """The default state every constructor form starts from."""
        self.name        = 'unnamed'
        self.fileinfo    = {'PATH':'','NAME':'','TYPE':'csv'}
        self._set_axes('Wavelength', 'nm', 'Absorbance', 'absorbance')
        self.technique   = None
        self.history     = []
        self.x      = np.empty(1)
        self.y      = np.empty(1)
        self.metadata    = {}

    def _init_from_data(self, x, y, *, technique=None,
                        x_quantity=None, x_unit=None,
                        y_quantity=None, y_unit=None,
                        name=None, metadata=None, history=None) -> None:
        """Populate from arrays. See :meth:`__init__` for the arguments."""
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.ndim != 1 or y.ndim != 1:
            raise ValueError(
                f"x and y must be one-dimensional; got {x.ndim}-D and {y.ndim}-D. "
                "A collection of spectra is a SpectrumCollection."
            )
        if len(x) != len(y):
            raise ValueError(
                f"x and y must be the same length; got {len(x)} and {len(y)}"
            )
        self.x, self.y = x, y

        if technique is not None:
            self.set_type(technique)

        # Explicit axis arguments win over the technique's conventions, and
        # mark the units authoritative -- otherwise a later set_type() would
        # treat them as defaults and quietly relabel the axis.
        axes = (x_quantity, x_unit, y_quantity, y_unit)
        if any(a is not None for a in axes):
            self._set_axes(
                x_quantity if x_quantity is not None else self.x_quantity,
                x_unit if x_unit is not None else self.x_unit,
                y_quantity if y_quantity is not None else self.y_quantity,
                y_unit if y_unit is not None else self.y_unit,
            )
            self.units_from_file = True

        if name is not None:
            self.name = name
        if metadata is not None:
            self.metadata = copy.deepcopy(dict(metadata))
            # set_type() records the technique here too; keep them consistent.
            if technique is not None:
                self.metadata['spec_type'] = technique
        if history is not None:
            self.history = list(history)

    @classmethod
    def read(cls, path, file_type=None) -> "Spectrum":
        """Read a spectrum from a file.

        The explicit form of ``Spectrum(path)``: it says that it does I/O, and
        it takes the format as a named argument rather than as the third of
        three positional strings.

        Parameters
        ----------
        path : str or path-like
            The file to read.
        file_type : str, optional
            One of the registered formats (``"dpt"``, ``"jcamp"``, ``"csv"``,
            ``"spy"`` ...). Inferred from the extension when omitted.

        Returns
        -------
        Spectrum

        Examples
        --------
        >>> from spectroscopy import datasets
        >>> s = Spectrum.read(datasets.path('ethanol'))
        >>> len(s) > 0
        True

        See Also
        --------
        SpectrumCollection.from_files : many files at once.
        """
        path = os.fspath(path)
        if file_type is None:
            return cls(path)
        directory, filename = os.path.split(path)
        return cls(directory, filename, file_type)

    def __str__(self) -> str:       # Users string version of object
        return f"Spectrum :{self.name} {self.y_label} vs {self.x_label}"

    def __repr__(self) -> str:      # Developper string version of object
        span = f"{self.x[0]:g}-{self.x[-1]:g}" if len(self.x) else "empty"
        technique = f" {self.technique}" if self.technique else ""
        steps = f", {len(self.history)} steps" if self.history else ""
        return (f"<Spectrum {self.name!r}{technique}: {len(self)} points, "
                f"{span} {self.x_unit}{steps}>")

    def __len__(self) -> int:       # Length of object
        return len(self.x)

##=============================================================================
#
#   Axes: quantity + unit are the stored truth, labels are derived.
#
#   A label may still be set directly -- the csv reader does exactly that with
#   whatever text the file's header row contained -- in which case the override
#   wins until set_type() replaces it.
#
##=============================================================================

    def _set_axes(self, x_quantity, x_unit, y_quantity, y_unit) -> None:
        """Set both axes' quantity and unit, clearing any label overrides."""
        self.x_quantity = x_quantity
        self.x_unit = x_unit
        self.y_quantity = y_quantity
        self.y_unit = y_unit
        self._x_label_override = None
        self._y_label_override = None
        #: True once a reader has established units from the file's own
        #: metadata, which set_type() then leaves alone.
        self.units_from_file = False

    @property
    def x_label(self) -> str:
        """Display label for the x axis."""
        if self._x_label_override is not None:
            return self._x_label_override
        return compose_label(self.x_quantity, self.x_unit)

    @x_label.setter
    def x_label(self, value) -> None:
        self._x_label_override = value

    @property
    def y_label(self) -> str:
        """Display label for the y axis."""
        if self._y_label_override is not None:
            return self._y_label_override
        return compose_label(self.y_quantity, self.y_unit)

    @y_label.setter
    def y_label(self, value) -> None:
        self._y_label_override = value

    @property
    def path(self):
        """Where this spectrum was read from, or None if built in memory."""
        name = self.fileinfo.get('NAME')
        if not name:
            return None
        return os.path.join(self.fileinfo.get('PATH', ''), name)

    @property
    def reversed_x(self) -> bool:
        """Whether this technique is conventionally plotted high-to-low."""
        return self.technique in REVERSED_AXIS_TECHNIQUES

##=============================================================================
#
#   Provenance. Every operation that returns a new Spectrum records what it
#   did and with which parameters, so a chain built interactively can later be
#   turned back into a runnable Pipeline (roadmap section 2.3).
#
##=============================================================================

    def _derive(self, y=None, x=None, step=None, name=None) -> "Spectrum":
        """
        Build a new Spectrum from this one, optionally with new data, and
        append ``step`` to its history. The single funnel for every non-mutating
        operation, so nothing can quietly skip the history.
        """
        result = Spectrum(self)
        if x is not None:
            result.x = np.asarray(x)
        if y is not None:
            result.y = np.asarray(y)
        if name is not None:
            result.name = name
        if step is not None:
            result.history = list(result.history) + [step]
        return result

    def _arithmetic(self, other, operation, symbol, reflected=False):
        """Shared body for the arithmetic dunders, including history."""
        if isinstance(other, Spectrum):
            if len(other) != len(self):
                raise ValueError(
                    f"cannot combine spectra of different lengths "
                    f"({len(self)} and {len(other)}); resample one onto the "
                    f"other's axis first, e.g. b = b.resample(a.x)"
                )
            # Equal length is not the same as equal axis. Warn rather than
            # raise: spectra from one instrument routinely differ in the last
            # decimal, and refusing those would break real workflows. The
            # tolerance is loose enough for float noise, tight enough to catch
            # genuinely different axes (defect D3).
            if len(self) and not np.allclose(self.x, other.x, rtol=1e-6, atol=0.0):
                worst = np.max(np.abs(self.x - other.x))
                warnings.warn(
                    f"combining spectra sampled at different x positions "
                    f"(largest difference {worst:g} {self.x_unit}); the result "
                    f"keeps this spectrum's axis. Resample first to be sure: "
                    f"other.resample(self.x)",
                    RuntimeWarning, stacklevel=3,
                )
            operand_y = other.y
        elif isinstance(other, (int, float)):
            operand_y = other
        else:
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)

        y = operation(operand_y, self.y) if reflected else operation(self.y, operand_y)
        step = ProcessingStep(symbol, {"operand": describe_operand(other),
                                       "reflected": reflected})
        return self._derive(y=y, step=step)

    def __add__(self, other):
        return self._arithmetic(other, operator.add, "add")

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other ):
        return self._arithmetic(other, operator.sub, "subtract")

    def __rsub__(self, other ):
        return self._arithmetic(other, operator.sub, "subtract", reflected=True)

    def __neg__(self):
        return 0.0 - self

    def __pos__(self):
        return self + 0.0

    def __mul__(self, other):
        return self._arithmetic(other, operator.mul, "multiply")

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        return self._arithmetic(other, operator.truediv, "divide")

    def __rtruediv__(self, other):
        return self._arithmetic(other, operator.truediv, "divide", reflected=True)

##=============================================================================
#
#   Some samples that manipulate the data.
#
##=============================================================================
    def subtract_reference(self, other, factor=1.0) -> "Spectrum":
        """
        Subtract a scaled reference spectrum, returning a new Spectrum.

        The scale factor is the point of this method. Across the notebooks the
        reference subtraction is always written ``sample - 0.995 * buffer``,
        with the multiplier tuned by eye per sample and kept in a bare list
        positionally aligned to a file list. Passing it here means it lands in
        ``.history`` as a recorded parameter instead.

        Parameters
        ----------
        other : Spectrum
            The reference (buffer, water, water vapour, PEG ...).
        factor : float
            Scaling applied to ``other`` before subtraction.

        Returns
        -------
        Spectrum
            A new spectrum; the original is untouched.
        """
        if not isinstance(other, Spectrum):
            raise TypeError(spectroscopy.messages.SPEC_MATH_ERR)
        if len(other) != len(self):
            raise ValueError(
                f"reference has {len(other)} points but the spectrum has "
                f"{len(self)}; resample it first, e.g. ref.resample(spec.x)"
            )

        result = self._derive(
            y=self.y - factor * other.y,
            step=ProcessingStep("subtract_reference", {
                "reference": describe_operand(other),
                "factor": float(factor),
            }),
        )

        own = self.metadata.get('reference')
        theirs = other.metadata.get('reference')
        their_sample = other.metadata.get('sample')
        if own == theirs:
            # Same reference on both sides: it cancels, and what is left is
            # this sample relative to the other one.
            result.metadata['reference'] = their_sample
        elif own is None:
            result.metadata['reference'] = their_sample
        else:
            result.metadata['reference'] = f"{their_sample} + {own} - {theirs}"
        return result

##=============================================================================
#
#   Processing. Each of these returns a NEW Spectrum and records a step.
#
##=============================================================================

    def to(self, x_unit=None, y_unit=None) -> "Spectrum":
        """
        Convert to different axis units, returning a new Spectrum.

        ``spectrum.to("nm")`` and ``spectrum.to(y_unit="%T")`` are the two
        everyday forms; both axes can be converted at once.

        A reciprocal x conversion (nm to cm^-1, or either to eV) reverses the
        ordering of the axis, so the points are re-sorted ascending and y is
        reordered with them -- otherwise every later interpolation, hull and
        crop would be working on a descending axis.

        The quantity name is updated along with the unit where the pairing is
        unambiguous, so a wavelength axis converted to cm^-1 is relabelled
        'Wavenumber' rather than staying 'Wavelength (cm^-1)'.
        """
        if x_unit is None and y_unit is None:
            raise ValueError("to() needs an x_unit, a y_unit, or both")

        result = Spectrum(self)
        step_params = {}

        if x_unit is not None and x_unit != self.x_unit:
            converted, flipped = units.convert_x(self.x, self.x_unit, x_unit)
            result.x = converted
            result.y = np.copy(self.y)
            if flipped:
                order = np.argsort(converted)
                result.x = converted[order]
                result.y = result.y[order]
            result.x_unit = x_unit
            result.x_quantity = QUANTITY_FOR_X_UNIT.get(x_unit, self.x_quantity)
            result._x_label_override = None
            step_params.update({"x_unit": x_unit, "from_x_unit": self.x_unit})

        if y_unit is not None and y_unit != self.y_unit:
            result.y = units.convert_y(result.y, self.y_unit, y_unit)
            result.y_unit = y_unit
            result.y_quantity = QUANTITY_FOR_Y_UNIT.get(y_unit, self.y_quantity)
            result._y_label_override = None
            step_params.update({"y_unit": y_unit, "from_y_unit": self.y_unit})

        if not step_params:
            return result            # already in the requested units
        result.history = list(result.history) + [ProcessingStep("to", step_params)]
        return result

    def crop(self, x_min=None, x_max=None) -> "Spectrum":
        """
        Restrict the spectrum to an x range, returning a new Spectrum.

        Bounds are order-insensitive, so ``crop(1800, 900)`` and
        ``crop(900, 1800)`` agree -- FTIR ranges get quoted both ways.

        This replaces the four-line boolean-mask idiom that appears in 19
        notebooks, and unlike that idiom it cannot leave x and y out of step.
        """
        low = -np.inf if x_min is None else min(x_min, x_max if x_max is not None else x_min)
        high = np.inf if x_max is None else max(x_max, x_min if x_min is not None else x_max)
        if x_min is not None and x_max is not None:
            low, high = min(x_min, x_max), max(x_min, x_max)

        mask = (self.x >= low) & (self.x <= high)
        if not mask.any():
            raise ValueError(
                f"crop({x_min}, {x_max}) selects no points; the spectrum spans "
                f"{self.x.min():g} to {self.x.max():g} {self.x_unit}"
            )
        return self._derive(
            x=self.x[mask], y=self.y[mask],
            step=ProcessingStep("crop", {"x_min": x_min, "x_max": x_max}),
        )

    def clip(self, region) -> "Spectrum":
        """
        Deprecated alias for :meth:`crop` taking a boolean mask.

        The old version mutated in place and returned None, which made it
        impossible to record in history; it is now non-mutating like everything
        else. Prefer ``crop(x_min, x_max)``.
        """
        warnings.warn("Spectrum.clip() is deprecated and no longer mutates in "
                      "place; use crop(x_min, x_max).",
                      DeprecationWarning, stacklevel=2)
        mask = np.asarray(region, dtype=bool)
        return self._derive(
            x=self.x[mask], y=self.y[mask],
            step=ProcessingStep("clip", {"n_selected": int(mask.sum())}),
        )

    def smooth(self, method='savgol', parameters=None, *,
               window_length=15, polyorder=3) -> "Spectrum":
        """
        Return a smoothed copy.

        Parameters
        ----------
        method : str
            ``'savgol'`` (aliases ``'SG'``, ``'savitzky-golay'``) or
            ``'moving_average'``.
        parameters : sequence, optional
            Legacy positional form, ``[polyorder, window_length]``. Note the
            order is the reverse of scipy's; it is accepted for the existing
            notebooks and normalised into keywords before being recorded.
        """
        if parameters is not None:
            polyorder, window_length = parameters[0], parameters[1]

        y = common.smooth(self.y, method=method,
                          window_length=window_length, polyorder=polyorder)
        canonical = 'savgol' if str(method).lower() in (
            'savgol', 'sg', 'savitzky-golay', 'savitsky-golay') else str(method).lower()
        return self._derive(
            y=y, name=f"{self.name} smooth",
            step=ProcessingStep("smooth", {
                "method": canonical,
                "window_length": int(window_length),
                "polyorder": int(polyorder),
            }),
        )

    def derivative(self, order=2, *, window_length=11, polyorder=3) -> "Spectrum":
        """Savitzky-Golay derivative, as a new Spectrum."""
        y = common.derivative(self.y, order=order,
                              window_length=window_length, polyorder=polyorder)
        return self._derive(
            y=y, name=f"{self.name} d{order}",
            step=ProcessingStep("derivative", {
                "order": int(order),
                "window_length": int(window_length),
                "polyorder": int(polyorder),
            }),
        )

    def normalize(self, method='max', window=None) -> "Spectrum":
        """
        Rescale, returning a new Spectrum.

        ``method`` is one of ``max``, ``area``, ``minmax``, ``vector``, ``snv``.
        ``window`` is an optional ``(low, high)`` x range over which the scaling
        factor is computed -- e.g. ``normalize('max', (1050, 1080))`` for the
        glycan band.
        """
        y = common.normalize(self.x, self.y, method=method, window=window)
        return self._derive(
            y=y,
            step=ProcessingStep("normalize", {
                "method": str(method).lower(),
                "window": list(window) if window is not None else None,
            }),
        )

    def baseline(self, method='rubberband', parameters=None, *,
                 return_points=False, **kwargs) -> "Spectrum":
        """
        Construct a baseline as a new Spectrum -- it is *not* subtracted.

        Kept in this form deliberately (review section 5.3): it is useful to
        plot a baseline next to the data, and ``spectrum - spectrum.baseline()``
        reads well. Use :meth:`baseline_correct` when you want the corrected
        spectrum and a history entry describing the correction.

        Methods
        -------
        ``rubberband`` (aliases ``RB``, ``hull``)
            Convex-hull lower envelope. With ``return_points=True`` you also get
            the anchor points back, which can be adjusted and handed to the
            polynomial method as guide points.
        ``poly`` (alias ``POLY``)
            Give ``points=[...]`` to fit a polynomial of ``degree`` through the
            spectrum's values at those guide positions -- the everyday form --
            or ``coefficients=[...]`` to evaluate a polynomial you already know.
            The legacy positional ``baseline('POLY', [a, b, c])`` still means
            coefficients.
        ``als`` (alias ``asls``)
            Asymmetric least squares; ``lam`` sets smoothness, ``p`` asymmetry.
        """
        canonical = common.BASELINE_ALIASES.get(str(method).lower())
        if canonical is None:
            raise ValueError(
                f"Unknown baseline method {method!r}; "
                f"try 'rubberband', 'poly' or 'als'."
            )
        if parameters is not None and canonical == 'poly':
            kwargs.setdefault('coefficients', parameters)

        # Which hull arc is the baseline depends on which way the bands point:
        # absorbance rises from a baseline near zero, transmittance dips from
        # one near 100%. Taking the lower hull of a transmittance spectrum
        # traces the tips of the absorption bands, which is meaningless.
        if canonical == 'rubberband':
            kwargs.setdefault('side', common.baseline_side_for(self.y_unit))

        anchors = None
        if canonical == 'rubberband':
            if return_points:
                base, anchor_x, anchor_y = common.rubberband_baseline(
                    self.x, self.y, return_points=True, side=kwargs['side'])
                anchors = (anchor_x, anchor_y)
            else:
                base = common.rubberband_baseline(self.x, self.y,
                                                  side=kwargs['side'])
        else:
            base = common.baseline(self.x, self.y, method=canonical, **kwargs)

        recorded = {"method": canonical}
        recorded.update({k: (list(v) if isinstance(v, (list, tuple, np.ndarray)) else v)
                         for k, v in kwargs.items()})
        result = self._derive(
            y=base, name=f"{self.name} baseline",
            step=ProcessingStep("baseline", recorded),
        )
        return (result, *anchors) if anchors is not None else result

    def baseline_correct(self, method='rubberband', parameters=None,
                         mode=None, **kwargs) -> "Spectrum":
        """
        Remove a baseline, returning the corrected spectrum.

        Same arguments as :meth:`baseline`. This is the form that records the
        *correction* in history -- doing ``s - s.baseline()`` by hand records
        only an anonymous subtraction, which would not replay.

        Parameters
        ----------
        mode : {'subtract', 'divide'}, optional
            How to remove it. Absorbance is additive, so the baseline is
            subtracted; transmittance is multiplicative, so it is divided out
            and the result is a corrected transmittance still referenced to 1.
            Chosen from the y unit when not given.
        """
        base = self.baseline(method, parameters, **kwargs)
        canonical = common.BASELINE_ALIASES.get(str(method).lower(), method)

        if mode is None:
            mode = ('divide'
                    if common.baseline_side_for(self.y_unit) == 'upper'
                    else 'subtract')
        if mode == 'divide':
            with np.errstate(divide='ignore', invalid='ignore'):
                corrected = np.where(base.y != 0, self.y / base.y, np.nan)
        elif mode == 'subtract':
            corrected = self.y - base.y
        else:
            raise ValueError(
                f"mode must be 'subtract' or 'divide', not {mode!r}")

        # A subtractive baseline is supposed to sit *under* the data. If the
        # result is negative over much of the range, it did not -- almost
        # always because the guide points were not in baseline regions, or a
        # too-flexible polynomial bowed above the data between them. Saying so
        # beats letting someone carry a systematically depressed spectrum into
        # an integration.
        if mode == 'subtract' and len(corrected):
            span = float(np.nanmax(corrected) - np.nanmin(corrected))
            if span > 0:
                below = float(np.mean(corrected < -0.01 * span))
                if below > 0.25:
                    warnings.warn(
                        f"the baseline sits above the data over "
                        f"{100 * below:.0f}% of the range, so the corrected "
                        f"spectrum is mostly negative. Check that the guide "
                        f"points are in genuinely flat regions, or lower the "
                        f"polynomial degree.",
                        RuntimeWarning, stacklevel=2)

        recorded = {"method": canonical, "mode": mode}
        recorded.update({k: (list(v) if isinstance(v, (list, tuple, np.ndarray)) else v)
                         for k, v in kwargs.items()})
        if parameters is not None and canonical == 'poly':
            recorded.setdefault('coefficients', list(parameters))
        return self._derive(
            y=corrected,
            step=ProcessingStep("baseline_correct", recorded),
        )

    def find_peaks(self, method='second_derivative', *, troughs=None,
                   relative=False, **kwargs) -> PeakTable:
        """
        Detect peaks and return a :class:`~spectroscopy.peaks.PeakTable`.

        By default this works on the inverted second derivative, which is what
        every peak-picking cell in the notebooks does, so that shoulders on a
        broad band are found as well as maxima. Extra keywords (``height``,
        ``distance``, ``prominence``, ``width``) go to scipy's find_peaks.

        **The search direction follows the y unit.** Bands in transmittance,
        ``%T`` and reflectance are minima, so those spectra are searched
        downward without being asked; absorbance-like units are searched
        upward. Signed units -- dichroism, anisotropy -- are searched **both
        ways**, because for those a band pointing down is not a failure to
        point up. Pass ``troughs=`` to override, including ``troughs='both'``,
        which is what a difference spectrum wants: there bands go both ways
        and the unit still just says ``absorbance``, so only you know.
        ``result.properties['direction_from']`` says whether the direction came
        from the unit, from you, or was assumed.

        Getting this wrong used to be silent: searching a transmission
        spectrum upward finds the two inflection points that flank each band
        rather than the band, so one band at 1100 cm-1 came back as a pair at
        1086 and 1114.

        Unlike the old :meth:`peaks` stub, the result is a return value rather
        than something written into ``metadata``: analysis results and
        acquisition facts do not belong in the same dictionary.

        Note that ``height``/``prominence`` apply to the *detection signal*,
        which for the default method is the second derivative and so is orders
        of magnitude smaller than the spectrum. Pass ``relative=True`` to give
        them as fractions of that signal's range instead -- see
        :func:`spectroscopy.processing.common.detect_peaks`.
        """
        indices, properties = common.detect_peaks(
            self.x, self.y, method=method, troughs=troughs,
            relative=relative, y_unit=self.y_unit, **kwargs)
        return PeakTable(
            position=self.x[indices],
            height=self.y[indices],
            index=indices,
            prominence=properties.get('prominences'),
            width=properties.get('widths'),
            kind={'both': 'both', 'down': 'trough'}.get(
                properties['direction'], 'peak'),
            sign=properties['sign'],
            x_unit=self.x_unit, y_unit=self.y_unit,
            source=self.name,
            properties=properties,
        )

    def fit_peaks(self, positions=None, *, model='voigt', n_peaks=None,
                  **kwargs) -> "FitResult":
        """
        Decompose this spectrum into overlapping components.

        :meth:`find_peaks` answers "where are the maxima". This answers the
        question a crowded band actually poses -- "how much of each overlapping
        component is in here" -- and it is the quantity that band-ratio and
        secondary-structure work report.

        Parameters
        ----------
        positions : array_like, optional
            Starting positions, one per component. When omitted they are found
            with :meth:`find_peaks` using the second-derivative method, which
            is the right default: overlapping bands have no maxima of their
            own, so detecting on the spectrum itself finds too few.
        model : {'voigt', 'gaussian', 'lorentzian'}, default 'voigt'
            Pseudo-Voigt by default: a fixed pure shape can be wrong by
            tens of percentage points in composition while still fitting
            with R^2 > 0.97. See :func:`spectroscopy.fitting.fit_components`.
        n_peaks : int, optional
            With ``positions`` omitted, keep only the strongest this many.
        **kwargs
            Passed to :func:`spectroscopy.fitting.fit_components` --
            ``fwhm``, ``position_tolerance``, ``max_fwhm``, ``non_negative``.

        Returns
        -------
        FitResult

        Notes
        -----
        Fit a band that still sits on a background and the background ends up
        inside the component areas. Crop and baseline-correct first.

        Examples
        --------
        >>> import numpy as np
        >>> from spectroscopy.lineshapes import gauss
        >>> x = np.linspace(1600, 1700, 201)
        >>> y = gauss(x, 1652, 12, 1.0) + gauss(x, 1630, 14, 0.6)
        >>> fit = Spectrum(x, y).fit_peaks([1630, 1652])
        >>> len(fit)
        2
        """
        from spectroscopy.fitting import fit_components  # pylint: disable=C0415

        # Positions come out right either way now that detection follows the
        # unit, but areas do not. Beer-Lambert is linear in absorbance, so a
        # component area or a band ratio taken on transmittance is not a
        # quantity with a meaning -- and it will still fit beautifully, which
        # is the trap (see roadmap sections 19-20). Warn rather than convert:
        # silently changing what was fitted is its own kind of wrong.
        if units.is_valley_pointing(self.y_unit):
            warnings.warn(
                f"fitting a spectrum in {self.y_unit!r}: component areas and "
                f"band ratios are not linear in this unit, so the areas this "
                f"returns are not proportional to concentration. Convert "
                f"first with .to(y_unit='absorbance') if the amounts matter. "
                f"Positions are unaffected.",
                UserWarning, stacklevel=2,
            )

        if positions is None:
            found = self.find_peaks(method='second_derivative')
            if n_peaks is not None:
                found = found.strongest(n_peaks).sorted_by_position()
            if not len(found):
                raise ValueError(
                    "no starting positions found automatically; pass "
                    "positions= explicitly, or check that the spectrum has "
                    "been cropped to the band of interest"
                )
            positions = found.position

        result = fit_components(self.x, self.y, positions, model=model, **kwargs)
        result.x_unit = self.x_unit
        result.y_unit = self.y_unit
        result.source = self.name
        return result

##=============================================================================
#
#   Presentation
#
##=============================================================================

    def plot(self, ax=None, *args, **kwargs):
        """
        Plot on a matplotlib axes.

        A thin delegator to :func:`spectroscopy.viz.plot`, which sets the axis
        labels from this spectrum and reverses x for the techniques that are
        quoted high-to-low. Keeping the drawing in the viz layer is review
        crossing C2; ``ax`` stays first and optional so existing notebook calls
        (``spectrum.plot(ax)``) are unchanged.
        """
        from spectroscopy import viz  # pylint: disable=C0415
        return viz.plot(self, ax, *args, **kwargs)

##=============================================================================
#
#   Resampling
#
##=============================================================================

    def resample(self, x_values) -> "Spectrum":
        """
        Estimate the spectrum at a new set of x positions.

        Needed before arithmetic between spectra that were not sampled at the
        same points -- different spectrometer settings, or different machines.
        Arithmetic now raises a clear error pointing here rather than emitting a
        raw numpy broadcast failure.

        Uses a cubic spline, which extrapolates beyond the original range; be
        wary of the result outside it.
        """
        x_values = np.asarray(x_values, dtype=float)
        spline = CubicSpline(self.x, self.y)
        return self._derive(
            x=x_values, y=spline(x_values),
            step=ProcessingStep("resample", {
                "n_points": int(len(x_values)),
                "x_min": float(x_values.min()) if len(x_values) else None,
                "x_max": float(x_values.max()) if len(x_values) else None,
                "method": "cubic_spline",
            }),
        )

#
#   Some gui functions to implement
#
##=============================================================================

#   save_as( new_filename: str)
#   save()

#   edit_metadata()
#   fit_components()
#   adjust_components()
#   display()

    def set_sample( self, info ) -> None:
        """
        Set the information string about the sample in the metadata.
        """
        self.metadata['sample'] = info

    def set_reference( self, ref_name ) -> None:
        """
        Set the information string about the reference in the metadata.
        """
        self.metadata['reference'] = ref_name

    def set_parameter(self, value, name=None, unit=None) -> None:
        """
        Record the **continuous** condition this spectrum was measured at.

        ``sample`` is categorical -- two spectra either are the same sample or
        are not. A titration, a melt or a dilution series is not like that: the
        spectra differ along an axis with an order and a distance, and the
        analysis needs the number, not a label. -120 mV, 37 C, 5 uM.

        The value goes in ``metadata['parameter']`` and stays there, so that
        anything reading a collection knows where to look without being told;
        ``name`` and ``unit`` are stored alongside as labels for axes and
        reports, not as lookup keys.

        Parameters
        ----------
        value : float
            The condition. Must be a real number -- that is the whole point of
            this being separate from :meth:`set_sample`.
        name : str, optional
            What it is: ``'potential'``, ``'temperature'``, ``'concentration'``.
        unit : str, optional
            What it is in: ``'mV'``, ``'C'``, ``'uM'``.
        """
        try:
            value = float(value)
        except (TypeError, ValueError) as error:
            raise TypeError(
                f"parameter must be a number, got {value!r}. For a categorical "
                f"label use set_sample() instead -- the difference matters "
                f"because a parameter gets sorted and fitted along, and a "
                f"sample name does not."
            ) from error
        self.metadata['parameter'] = value
        if name is not None:
            self.metadata['parameter_name'] = name
        if unit is not None:
            self.metadata['parameter_unit'] = unit

    def set_type( self, spec_type, force_units=False ) -> None:
        """
        Set the technique, and with it the default axis quantities and units.

        Units the file itself declared are **kept**, not overwritten. A JCAMP
        file that says ``##YUNITS=TRANSMITTANCE`` holds transmittance; calling
        ``set_type('FTIR')`` should not relabel it as absorbance just because
        absorbance is the usual FTIR ordinate. Getting that wrong mislabels a
        figure and makes :meth:`to` silently do the wrong thing.

        Pass ``force_units=True`` to override anyway -- useful when a file's
        own metadata is known to be wrong.
        """
        if spec_type not in KNOWNSPECTYPES:
            raise TypeError(
                f"Unknown spectrum type {spec_type}; "
                f"known types are {', '.join(KNOWNSPECTYPES)}"
            )
        self.technique = spec_type
        self.metadata['spec_type'] = spec_type

        if self.units_from_file and not force_units:
            return

        data = KNOWNSPECTYPES[spec_type]
        self._set_axes(data['x_quantity'], data['x_unit'],
                       data['y_quantity'], data['y_unit'])

    def get_info( self ) -> str:
        """
        Get information about the spectrum as a (multiline) string
        """
        spec_type = self.technique or self.metadata.get('spec_type')
        sample    = self.metadata.get('sample')
        reference = self.metadata.get('reference')
        info = ""
        if spec_type:
            info = info + spec_type + ' spectrum'
        if sample:
            info = f'{info} of {sample} '
        if reference:
            info = f'{info} - {reference}'
        info = f'{info}\n{self.y_label} vs {self.x_label}'
        info = f'{info}\n{len(self)} points'
        if len(self):
            info = f'{info}, {self.x.min():g} to {self.x.max():g} {self.x_unit}'
        if self.history:
            info = f'{info}\nhistory:'
            for number, step in enumerate(self.history, start=1):
                info = f'{info}\n  {number}. {step}'
        info = f'{info}\n{self.metadata}'
        return info

    def describe_history( self ) -> str:
        """The processing history as a numbered, human-readable list."""
        if not self.history:
            return "No processing recorded."
        return "\n".join(f"{n}. {step}" for n, step in enumerate(self.history, 1))

##=============================================================================
#
#   Here are the routines to interface to different file formats and returned
#   data structures.
#
#   TODO: Populate the data structure properly from the dictionary
#
##=============================================================================
    def reload(self) -> None:
        """
        Load (or re-load) this spectrum from the file named in ``fileinfo``.

        Dispatch, encoding detection and the multi-spectrum question all live
        in :mod:`spectroscopy.io.registry` now; this is a thin adapter that
        keeps the existing in-place semantics.
        """
        filename = os.path.join(self.fileinfo['PATH'], self.fileinfo['NAME'])
        spectra = registry.read_spectra(filename, self.fileinfo['TYPE'])

        if len(spectra) != 1:
            raise ValueError(
                f"{self.fileinfo['NAME']} holds {len(spectra)} spectra; use "
                f"spectroscopy.io.read_spectra() to load them all"
            )

        # Copy everything the reader produced except the file bookkeeping,
        # which this object already owns. Deliberately not an explicit list of
        # attribute names: such a list has to be updated whenever the data
        # model grows a field, and one that silently falls out of step is what
        # defect D5 was.
        loaded = spectra[0]
        for attribute, value in vars(loaded).items():
            if attribute != 'fileinfo':
                setattr(self, attribute, value)

    def save_as(self, filename, file_type = 'spy', **kwargs) -> None:
        """
        Set the path name and type for the spectrum and then save the
        spectrum in the designated spot.

        Extra keywords go to the writer, so the options a format offers are
        reachable from here rather than only from
        :func:`spectroscopy.io.write_spectrum`::

            spectrum.save_as("pour_chloe.csv", "csv", decimal=',')
        """
        self.fileinfo['PATH'],self.fileinfo['NAME'] = os.path.split(filename)
        self.fileinfo['TYPE'] = file_type
        self.save(**kwargs)

    def save(self, **kwargs) -> None:
        """
        Write the spectrum to the file described by ``fileinfo``.

        The registry resolves the format *before* opening the file, so an
        unwritable type raises instead of truncating the target to nothing --
        the second half of defect D5.
        """
        filename = os.path.join(self.fileinfo['PATH'], self.fileinfo['NAME'])
        registry.write_spectrum(self, filename, self.fileinfo['TYPE'], **kwargs)
