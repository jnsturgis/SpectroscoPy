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

import numpy as np
from scipy.interpolate import CubicSpline

import spectroscopy.messages
from spectroscopy.history import ProcessingStep, describe_operand
from spectroscopy.io import csv, dpt, jcamp, spy
from spectroscopy.peaks import PeakTable
from spectroscopy.processing import common

#: Map file extension -> file type. Keys must be lower case; lookup lower-cases
#: the extension so that .DPT and .CSV work as well as .dpt and .csv.
FILE_EXTS = {
    '.csv':   'csv',
    '.tsv':   'tsv',
    '.txt':   'tsv',
    '.dpt':   'dpt',
    '.jcamp': 'jcamp',
    '.jdx':   'jcamp',
    '.dx':    'jcamp',
    '.spy':   'spy',
}

#: Every type here must be handled by BOTH Spectrum.reload() and
#: Spectrum.save(). The two match statements used to disagree, which meant an
#: unhandled type silently truncated the output file -- see review defect D5.
KNOWNFILETYPES = ('csv','tsv','dpt','jcamp','spy')
#: Rendering of each unit for a matplotlib axis label. An empty string means
#: dimensionless, so the label is just the quantity with no parentheses.
UNIT_LABELS = {
    'cm^-1':         r'cm$^{-1}$',
    'nm':            'nm',
    'um':            r'$\mu$m',
    'eV':            'eV',
    'absorbance':    '',
    'transmittance': '',
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


def compose_label(quantity, unit):
    """Build an axis label from a quantity and a unit, e.g. 'Wavenumber (cm^-1)'."""
    rendered = UNIT_LABELS.get(unit, unit)
    return f"{quantity} ({rendered})" if rendered else str(quantity)

def _infer_file_type( name ):
    """Work out if possible from the file extension the possible file types."""
    name, ext = os.path.splitext(name)
    return FILE_EXTS.get(ext.lower(), 'unknown')

class Spectrum:
    """A class for spectra."""

    def __init__(self, *args):
        """Initialize a Spectrum object.

        Simulate overloading the different possitilities for initializing
        a spectrum:

        Parameters
        ----------
        args a set of parameters that allow for different initiation routines
        depending on the parameters.


        Spectrum()                     : no arguments an empty spectrum.
        Spectrum( filename )           : a file
        Spectrum( filename, filetype ) : a file where filetype can be from
                 ('jcamp','csv')
        Spectrum( a_spectrum )         : a deepcopy of the original
        """

        if len(args) == 0 :
            #     This is the empty initializer
            self.name        = 'unnamed'
            self.fileinfo    = {'PATH':'','NAME':'','TYPE':'csv'}
            self._set_axes('Wavelength', 'nm', 'Absorbance', 'absorbance')
            self.technique   = None
            self.history     = []
            self.x      = np.empty(1)
            self.y      = np.empty(1)
            self.metadata    = {}

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

        elif isinstance(args[0], str ) :
            # This is the open method - initialize from a file
            # First parse the args to get filetype and full filename
            # Len 1: filename - infer type from 1
            # Len 2: path, filename - infer type  from 2
            # Len 3: path, filename, type

            self.name        = args[0]
            self.fileinfo    = {'PATH':'','NAME':'','TYPE':'csv'}
            self._set_axes('Wavelength', 'nm', 'Absorbance', 'absorbance')
            self.technique   = None
            self.history     = []
            self.x      = np.empty(1)
            self.y      = np.empty(1)
            self.metadata    = {}

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

            # Now read the file according to type
            if self.fileinfo['TYPE'] in KNOWNFILETYPES:
                self.reload()
            else:
                raise TypeError(f"Unknown filetype {self.fileinfo['TYPE']}")
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

        anchors = None
        if canonical == 'rubberband':
            if return_points:
                base, anchor_x, anchor_y = common.rubberband_baseline(
                    self.x, self.y, return_points=True)
                anchors = (anchor_x, anchor_y)
            else:
                base = common.rubberband_baseline(self.x, self.y)
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
                         **kwargs) -> "Spectrum":
        """
        Subtract a baseline, returning the corrected spectrum.

        Same arguments as :meth:`baseline`. This is the form that records the
        *correction* in history -- doing ``s - s.baseline()`` by hand records
        only an anonymous subtraction, which would not replay.
        """
        base = self.baseline(method, parameters, **kwargs)
        canonical = common.BASELINE_ALIASES.get(str(method).lower(), method)
        recorded = {"method": canonical}
        recorded.update({k: (list(v) if isinstance(v, (list, tuple, np.ndarray)) else v)
                         for k, v in kwargs.items()})
        if parameters is not None and canonical == 'poly':
            recorded.setdefault('coefficients', list(parameters))
        return self._derive(
            y=self.y - base.y,
            step=ProcessingStep("baseline_correct", recorded),
        )

    def find_peaks(self, method='second_derivative', *, troughs=False,
                   **kwargs) -> PeakTable:
        """
        Detect peaks and return a :class:`~spectroscopy.peaks.PeakTable`.

        By default this works on the inverted second derivative, which is what
        every peak-picking cell in the notebooks does, so that shoulders on a
        broad band are found as well as maxima. Extra keywords (``height``,
        ``distance``, ``prominence``, ``width``) go to scipy's find_peaks.

        Unlike the old :meth:`peaks` stub, the result is a return value rather
        than something written into ``metadata``: analysis results and
        acquisition facts do not belong in the same dictionary.
        """
        indices, properties = common.detect_peaks(
            self.x, self.y, method=method, troughs=troughs, **kwargs)
        return PeakTable(
            position=self.x[indices],
            height=self.y[indices],
            index=indices,
            prominence=properties.get('prominences'),
            width=properties.get('widths'),
            kind='trough' if troughs else 'peak',
            x_unit=self.x_unit, y_unit=self.y_unit,
            source=self.name,
            properties=properties,
        )

##=============================================================================
#
#   Presentation
#
##=============================================================================

    def plot(self, ax, *args, label=None, apply_labels=True, **kwargs):
        """
        Plot on a matplotlib axes.

        Also sets the axis labels from this spectrum's own quantity and unit,
        and reverses the x axis for techniques conventionally drawn high-to-low
        (FTIR, Raman) -- the two lines typed by hand in every figure cell in
        the notebooks. Pass ``apply_labels=False`` to suppress that.
        """
        lines = ax.plot(self.x, self.y, *args,
                        label=self.name if label is None else label, **kwargs)
        if apply_labels:
            ax.set_xlabel(self.x_label)
            ax.set_ylabel(self.y_label)
            if self.reversed_x and not ax.xaxis_inverted():
                ax.invert_xaxis()
        return lines

##=============================================================================
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

    def set_type( self, spec_type ) -> None:
        """
        Set the technique, and with it the axis quantities and units.

        Clears any axis-label overrides, so a label picked up from a file
        header is replaced by the technique's proper label.
        """
        if spec_type not in KNOWNSPECTYPES:
            raise TypeError(
                f"Unknown spectrum type {spec_type}; "
                f"known types are {', '.join(KNOWNSPECTYPES)}"
            )
        data = KNOWNSPECTYPES[spec_type]
        self._set_axes(data['x_quantity'], data['x_unit'],
                       data['y_quantity'], data['y_unit'])
        self.technique = spec_type
        self.metadata['spec_type'] = spec_type

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
        Reload the spectrum from the file, or load a first time after setting fileinfo
        """
        filename = os.path.join(self.fileinfo['PATH'],self.fileinfo['NAME'])
        with open( filename, encoding="utf-8") as f:
            match self.fileinfo['TYPE']:
                case 'jcamp':
                    jcamp.read(f, self)
                case 'csv':
                    csv.read(f, self )
                case 'tsv':
                    csv.read(f, self, delimiter='\t')
                case 'dpt':
                    # Separator is sniffed per file: 140 of the 825 .dpt files
                    # in the notebooks are comma separated, not tab. See
                    # spectroscopy/io/dpt.py.
                    dpt.read(f, self)
                case 'spy':
                    spy.read(f, self, format='0.0')
                case _:
                    raise ValueError(
                        f"Cannot read '{self.fileinfo['TYPE']}' files; "
                        f"known types are {', '.join(KNOWNFILETYPES)}"
                    )

    def save_as(self, filename, file_type = 'spy') -> None:
        """
        Set the path name and type for the spectrum and then save the
        spectrum in the designated spot.
        """
        self.fileinfo['PATH'],self.fileinfo['NAME'] = os.path.split(filename)
        self.fileinfo['TYPE'] = file_type
        self.save()

    def save(self) -> None:
        """
        Write the spectrum to the file described by file info.
        """
        file_type = self.fileinfo['TYPE']

        # Validate BEFORE opening. open(..., 'w') truncates immediately, so a
        # check inside the `with` would already have destroyed the target -- an
        # unhandled type used to leave a 0-byte file behind with no error.
        if file_type not in KNOWNFILETYPES:
            raise ValueError(
                f"Cannot write '{file_type}' files; "
                f"known types are {', '.join(KNOWNFILETYPES)}"
            )

        filename = os.path.join(self.fileinfo['PATH'],self.fileinfo['NAME'])
        with open(filename, 'w', encoding="utf-8") as f:
            match file_type:
                case 'jcamp':
                    jcamp.write(f, self)
                case 'csv':
                    csv.write(f, self )
                case 'tsv':
                    csv.write(f, self, delimiter='\t')
                case 'dpt':
                    dpt.write(f, self)
                case 'spy':
                    spy.write(f, self )
