# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Working with many spectra at once.

Most experiments give you a folder of files rather than one spectrum: three
replicates of each sample, a titration measured at a dozen potentials, a melt
recorded every second degree. A collection holds them together so you can
treat them as the one thing they are.

It behaves like a list. You can index it, slice it, loop over it and ask how
long it is. What it adds is that anything you do to it happens to every
spectrum in it, and you get a new collection back -- the originals are never
altered.

    >>> import spectroscopy as spc
    >>> folder = spc.datasets.replicate_directory()
    >>> spectra = spc.SpectrumCollection.from_files(folder + "/*.dpt")
    >>> len(spectra)
    9
    >>> averages = {name: group.mean()
    ...             for name, group in spectra.group_by('sample').items()}
    >>> sorted(averages)
    ['CelluloseX', 'Glucose', 'H2O']

Averaging replicates is the thing this exists for. Grouping by sample name
finds the replicates by what they are rather than by where they sit in a list,
so adding a file later cannot quietly spoil an average -- which is exactly what
happens when you write spectra[0] + spectra[1] + spectra[2] and then measure a
fourth. Use median() instead of mean() when you suspect one replicate went
wrong.

Cropping, baseline correction, smoothing, normalising and resampling all work
on a collection just as they do on a single spectrum. For anything else, map()
applies a function of your own to each one.

When you want the numbers rather than the objects, to_matrix() gives you the
common x axis and one row per spectrum as plain numpy -- which is the form PCA,
NMF and anything else from scikit-learn expect.

Two things worth knowing before you rely on them.

A series usually varies in some measured quantity: temperature, potential,
concentration. If you tell the collection what that quantity is, it will gather
the values for you and put the spectra in order along it. Do put them in order
explicitly, because loading sorts filenames as *text*, and text order puts
-120 mV before -20 mV -- a titration loaded that way plots as a scribble and
fits as nonsense.

Facts about the set as a whole -- what the varying quantity is called and what
it is measured in, or for a set of reference spectra where it came from and on
what terms -- belong on the collection, in ``info``. Facts about one spectrum
stay on that spectrum. Keeping them apart is what stops the two disagreeing:
there is one copy of anything that describes the whole set, so it cannot say
degrees in one place and kelvin in another.
"""

from __future__ import annotations

import glob as _glob
import os
import re
import warnings
from collections.abc import Sequence

import numpy as np

from spectroscopy.history import ProcessingStep
from spectroscopy.spectra import Spectrum

__all__ = ['SpectrumCollection']

#: Sentinel for "argument not given", where ``None`` is a legitimate value.
_UNSET = object()


def _same(first, second) -> bool:
    """
    Equality that survives numpy arrays in ``info``.

    ``a == b`` on two arrays is an array, and ``bool()`` of it raises. Merging
    the provenance of two sets is not the place to discover that.
    """
    try:
        return bool(np.all(np.asarray(first) == np.asarray(second)))
    except (ValueError, TypeError):
        return first is second


def _parameter_reader(spec):
    """
    Work out how to read the measured value out of each filename.

    A callable is used as given. A string is a regular expression with one
    capture group, searched against the whole path -- so a titration named
    ``run3/sample_-120mV.dpt`` gives up its potential to ``r'(-?\\d+)mV'``, and
    a temperature held in the directory is reachable too.
    """
    if callable(spec):
        return spec
    if not isinstance(spec, str):
        raise TypeError(
            f"parameter_from must be a callable or a regular expression "
            f"string, got {type(spec).__name__}"
        )

    pattern = re.compile(spec)
    if pattern.groups != 1:
        raise ValueError(
            f"the parameter_from pattern {spec!r} has {pattern.groups} capture "
            f"groups; it needs exactly one, around the number itself. For "
            f"'sample_-120mV.dpt' that is r'(-?\\d+)mV'."
        )

    def read(path):
        match = pattern.search(os.fspath(path))
        if match is None:
            raise ValueError(
                f"the parameter_from pattern {spec!r} does not match "
                f"{path!r}. Every file needs a parameter, otherwise the "
                f"spectra cannot be put in order along it."
            )
        return float(match.group(1))

    return read


class SpectrumCollection(Sequence):
    """
    An ordered collection of ``Spectrum``.

    Parameters
    ----------
    spectra : iterable of Spectrum
    name : str, optional
    info : dict, optional
        Facts about the **set** rather than about any spectrum in it: what the
        series parameter is called, where a reference set came from, the terms
        it arrived under. See ``spectroscopy.metadata.SET_LEVEL``.

        Per-item data belongs in that item's own ``metadata``, and the
        collection offers *gathered views* of it -- ``parameters``,
        ``samples``. That way there is only ever one list and it cannot
        fall out of step with the spectra. Set-level data has the opposite
        problem: stored per item it can disagree with itself, and the
        disagreement is invisible. Hence two homes.

        ``info`` survives ``crop``, ``select``, ``map`` and
        slicing, exactly as ``name`` does -- a subset of SP175 is still SP175
        data, and still carries SP175's citation condition.
    """

    def __init__(self, spectra=(), name=None, info=None):
        self._spectra = list(spectra)
        for item in self._spectra:
            if not isinstance(item, Spectrum):
                raise TypeError(
                    f"{type(self).__name__} takes Spectrum objects, "
                    f"got {type(item).__name__}"
                )
        self.name = name
        self.info = dict(info) if info else {}

    def _like(self, spectra, name=_UNSET):
        """
        A collection of the same class as this one, keeping ``name`` and
        ``info``.

        Every operation that returns a collection goes through here, so a
        subclass -- ``ReferenceSet`` -- stays itself
        through ``select``, ``crop`` and slicing rather than degrading to a
        plain collection and losing its provenance on the way.
        """
        return type(self)(spectra,
                          name=self.name if name is _UNSET else name,
                          info=dict(self.info))

    # -- Sequence protocol -------------------------------------------------

    def __len__(self) -> int:
        return len(self._spectra)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self._like(self._spectra[index])
        return self._spectra[index]

    def __repr__(self) -> str:
        label = f" {self.name!r}" if self.name else ""
        samples = sorted({s.metadata.get('sample') for s in self._spectra}
                         - {None})
        detail = f", {len(samples)} samples" if samples else ""
        values = self.parameters
        if len(values) and not np.isnan(values).all():
            what = self.parameter_name or 'parameter'
            unit = f" {self.parameter_unit}" if self.parameter_unit else ""
            detail += (f", {what} {np.nanmin(values):g} to "
                       f"{np.nanmax(values):g}{unit}")
        return f"<SpectrumCollection{label}: {len(self)} spectra{detail}>"

    def __add__(self, other) -> SpectrumCollection:
        if isinstance(other, SpectrumCollection):
            merged = {key: value for key, value in self.info.items()
                      if key in other.info and _same(other.info[key], value)}
            return SpectrumCollection(self._spectra + list(other), info=merged)
        if isinstance(other, Spectrum):
            return self._like(self._spectra + [other])
        return NotImplemented

    # -- construction ------------------------------------------------------

    @classmethod
    def from_files(cls, patterns, file_type=None, *, technique=None,
                   sample_from=None, parameter_from=None, parameter_name=None,
                   parameter_unit=None, sort=True, **kwargs):
        """
        Load many files at once.

        Parameters
        ----------
        patterns : str or iterable of str
            Glob patterns and/or plain paths. ``"data/*.dpt"`` and
            ``["a.dpt", "b.dpt"]`` both work.
        file_type : str, optional
            Force a type; otherwise inferred per file from the extension.
        technique : str, optional
            Passed to ``set_type`` on each spectrum, e.g. ``'ATR-FTIR'``.
        sample_from : callable, optional
            ``f(path) -> sample name``. Defaults to the basename up to the
            first dot, which turns ``08/PG_coli2.5.dpt`` into ``PG_coli2`` --
            the convention already used throughout the notebooks.
        parameter_from : str or callable, optional
            Where to find the **continuous** condition each spectrum was
            measured at -- the potential of a redox titration, the temperature
            of a melt, the concentration of a dilution series. Either a
            callable ``f(path) -> float`` or a regular expression with one
            capture group, searched against the path::

                SpectrumCollection.from_files(
                    "titration/*.dpt", parameter_from=r'(-?\\d+)mV',
                    parameter_name='potential', parameter_unit='mV')

            See ``Spectrum.set_parameter`` for why this is not just another
            ``sample``. Use ``sorted_by_parameter`` afterwards, because
            ``sort=True`` orders the
            files as *text*, which puts ``-120mV`` before ``-20mV``.
        parameter_name, parameter_unit : str, optional
            Labels for the parameter, for axes and reports.
        """
        if isinstance(patterns, (str, os.PathLike)):
            patterns = [patterns]

        paths = []
        for pattern in patterns:
            pattern = os.fspath(pattern)
            matches = _glob.glob(pattern, recursive=True)
            paths.extend(sorted(matches) if sort else matches)
            if not matches and os.path.exists(pattern):
                paths.append(pattern)
        if not paths:
            raise FileNotFoundError(f"no files matched {patterns!r}")

        if sample_from is None:
            def sample_from(path):
                return os.path.basename(path).split('.')[0]

        read_parameter = (None if parameter_from is None
                          else _parameter_reader(parameter_from))

        spectra = []
        for path in paths:
            directory, filename = os.path.split(path)
            arguments = [directory + os.sep if directory else "", filename]
            if file_type is not None:
                arguments.append(file_type)
            spectrum = Spectrum(*arguments, **kwargs)
            if technique is not None:
                spectrum.set_type(technique)
            spectrum.set_sample(sample_from(path))
            if read_parameter is not None:
                spectrum.set_parameter(read_parameter(path),
                                       name=parameter_name,
                                       unit=parameter_unit)
            spectra.append(spectrum)

        info = {key: value for key, value in
                (('parameter_name', parameter_name),
                 ('parameter_unit', parameter_unit)) if value is not None}
        return cls(spectra, info=info)

    def save_as(self, filename, file_type='spy', **kwargs) -> None:
        """
        Write the whole collection to one file. The partner of
        ``spectroscopy.io.read_spectra``, which reads it back.

        ``.spy`` holds either one spectrum or a set, and only it has room for
        what is true of the *set* -- its name, where it came from, the terms it
        arrived under, the units it is in. Asking for a format that stores one
        spectrum raises rather than quietly writing the first of forty.

        To write the spectra as separate files instead, loop::

            for spectrum in collection:
                spectrum.save_as(f"{spectrum.name}.spy")

        which is fine for the spectra and saves nothing about the set.
        """
        from spectroscopy.io import registry  # pylint: disable=C0415
        registry.write_collection(self, filename, file_type, **kwargs)

    # -- grouping and reduction -------------------------------------------

    def group_by(self, key='sample'):
        """
        Group into ``{value: SpectrumCollection}``, preserving order.

        ``key`` is a metadata field name, or a callable taking a Spectrum.
        """
        getter = key if callable(key) else (lambda s: s.metadata.get(key))
        groups: dict = {}
        for spectrum in self._spectra:
            groups.setdefault(getter(spectrum), []).append(spectrum)
        return {value: self._like(members, name=str(value))
                for value, members in groups.items()}

    def _stack(self):
        """Common x plus the y values as a (n_spectra, n_points) array."""
        if not self._spectra:
            raise ValueError("empty collection")
        lengths = {len(s) for s in self._spectra}
        if len(lengths) != 1:
            raise ValueError(
                f"spectra have different lengths {sorted(lengths)}; resample "
                f"them onto a common axis first, e.g. "
                f"collection.resample(collection[0].x)"
            )
        return self._spectra[0].x, np.vstack([s.y for s in self._spectra])

    def _reduce(self, function, step_name):
        x, matrix = self._stack()
        template = self._spectra[0]
        result = template._derive(              # pylint: disable=protected-access
            x=x, y=function(matrix, axis=0),
            step=ProcessingStep(step_name, {"n_spectra": len(self)}),
        )
        return result

    def mean(self) -> Spectrum:
        """Point-wise mean -- the replicate average."""
        result = self._reduce(np.mean, "mean")
        result.name = f"{self.name or 'collection'} mean"
        return result

    def std(self) -> Spectrum:
        """Point-wise standard deviation across the collection."""
        result = self._reduce(lambda m, axis: np.std(m, axis=axis, ddof=1), "std")
        result.name = f"{self.name or 'collection'} std"
        return result

    def sem(self) -> Spectrum:
        """Point-wise standard error of the mean."""
        result = self._reduce(
            lambda m, axis: np.std(m, axis=axis, ddof=1) / np.sqrt(m.shape[axis]),
            "sem")
        result.name = f"{self.name or 'collection'} sem"
        return result

    def median(self) -> Spectrum:
        """Point-wise median -- steadier than the mean against a bad replicate."""
        result = self._reduce(np.median, "median")
        result.name = f"{self.name or 'collection'} median"
        return result

    # -- batch operations --------------------------------------------------

    def map(self, function) -> SpectrumCollection:
        """Apply ``function`` to each spectrum, returning a new collection."""
        return self._like([function(s) for s in self._spectra])

    def _batch(self, method_name, *args, **kwargs) -> SpectrumCollection:
        return self.map(lambda s: getattr(s, method_name)(*args, **kwargs))

    def crop(self, x_min=None, x_max=None):
        """Crop every spectrum."""
        return self._batch('crop', x_min, x_max)

    def baseline_correct(self, method='rubberband', parameters=None, **kwargs):
        """Baseline-correct every spectrum."""
        return self._batch('baseline_correct', method, parameters, **kwargs)

    def normalize(self, method='max', window=None):
        """Normalise every spectrum."""
        return self._batch('normalize', method, window)

    def smooth(self, method='savgol', parameters=None, **kwargs):
        """Smooth every spectrum."""
        return self._batch('smooth', method, parameters, **kwargs)

    def resample(self, x_values):
        """Put every spectrum onto a common axis."""
        return self._batch('resample', x_values)

    def subtract_reference(self, reference, factor=1.0):
        """Subtract one reference from every spectrum, with a common factor."""
        return self._batch('subtract_reference', reference, factor)

    def select(self, predicate) -> SpectrumCollection:
        """The spectra for which ``predicate(spectrum)`` is true."""
        return self._like([s for s in self._spectra if predicate(s)])

    # -- aliases -----------------------------------------------------------
    #
    # The spellings a pandas user reaches for without thinking. Guessability
    # audit, 2026-08-05.

    def groupby(self, key='sample'):
        """Alias for ``group_by``, which is how pandas spells it."""
        return self.group_by(key)

    def filter(self, predicate) -> SpectrumCollection:
        """Alias for ``select``."""
        return self.select(predicate)

    def normalise(self, method='max', window=None):
        """Alias for ``normalize``, for British spelling."""
        return self.normalize(method, window)

    # -- interop -----------------------------------------------------------

    def to_matrix(self, with_parameter=False):
        """
        ``(x, X)`` with ``X`` of shape (n_spectra, n_points).

        This is the handover to PCA/NMF/ICA, and the reason those notebooks no
        longer need to assemble ``np.array(aligned)`` by hand.

        ``with_parameter=True`` returns ``(x, X, parameter)`` instead, which is
        what a titration, a melt or a calibration wants: those analyses fit
        along the parameter, so all three arrays travel together. It raises if
        any spectrum lacks one, rather than handing a fit a silent ``nan``.
        """
        x, matrix = self._stack()
        if not with_parameter:
            return x, matrix

        parameters = self.parameters
        missing = np.isnan(parameters)
        if missing.any():
            names = [self._spectra[i].name for i in np.flatnonzero(missing)][:5]
            raise ValueError(
                f"{int(missing.sum())} of {len(self)} spectra have no "
                f"parameter, so there is nothing to fit against: "
                f"{', '.join(str(n) for n in names)}"
                f"{' ...' if int(missing.sum()) > 5 else ''}. Set them with "
                f"from_files(parameter_from=...) or with_parameters()."
            )
        return x, matrix, parameters

    @property
    def samples(self):
        """Sample name of each spectrum, in order."""
        return [s.metadata.get('sample') for s in self._spectra]

    @property
    def parameters(self):
        """
        The continuous parameter of each spectrum, in order, as an array.

        ``nan`` where a spectrum has none -- unlike ``to_matrix``, which
        refuses. Reading the values is a reasonable thing to do on a
        part-labelled collection; fitting them is not.
        """
        return np.array([s.metadata.get('parameter', np.nan)
                         for s in self._spectra], dtype=float)

    @property
    def parameter_name(self):
        """What the parameter is -- 'potential', 'temperature'."""
        return self._parameter_label('parameter_name')

    @property
    def parameter_unit(self):
        """What the parameter is in -- 'mV', 'C'."""
        return self._parameter_label('parameter_unit')

    def _parameter_label(self, key):
        """
        From ``info``, falling back to the per-spectrum copies.

        The fallback is how this used to work, kept going for collections
        assembled by hand. It is the reason this is not simply
        ``self.info.get(key)``: a spectrum labelled by ``set_parameter`` still
        carries its own copy, and dropping the fallback would lose the label
        of every collection built that way.

        Where the copies disagree there is no answer, and the old behaviour --
        return ``None``, indistinguishable from *never set* -- hid a real
        conflict. It now says so.
        """
        if key in self.info:
            return self.info[key]
        labels = {s.metadata.get(key) for s in self._spectra} - {None}
        if len(labels) > 1:
            warnings.warn(
                f"the spectra in this collection disagree about {key}: "
                f"{sorted(labels)}. It describes the series rather than any "
                f"one spectrum, so set it once on the collection -- "
                f"collection.info[{key!r}] = ... -- rather than per spectrum.",
                stacklevel=3,
            )
            return None
        return labels.pop() if labels else None

    def with_parameters(self, values, name=None, unit=None):
        """
        Attach parameters from a sequence aligned with the spectra.

        The partner to ``from_files(parameter_from=...)``, for when the numbers
        live in a lab notebook rather than in the filenames -- which is the
        usual case, since the potentiostat does not write into the
        spectrophotometer's files.

        **Matched by position**, and nothing can check that for you: a list in
        the order the experiment was run is silently wrong if the files loaded
        in a different one, and they usually did, because ``from_files`` sorts
        paths as text. Look at the collection first. Where the number is in the
        filename, ``parameter_from`` avoids the question entirely by reading
        each one off the file it belongs to.

        Returns a **new** collection; the originals are untouched.
        """
        values = np.asarray(values, dtype=float).ravel()
        if len(values) != len(self):
            raise ValueError(
                f"got {len(values)} parameters for {len(self)} spectra. They "
                f"are matched by position, so the counts must agree -- check "
                f"whether a file was added or a blank crept into the glob."
            )
        spectra = []
        for spectrum, value in zip(self._spectra, values):
            copied = spectrum._derive()  # pylint: disable=protected-access
            copied.set_parameter(value, name=name, unit=unit)
            spectra.append(copied)
        labelled = self._like(spectra)
        for key, value in (('parameter_name', name), ('parameter_unit', unit)):
            if value is not None:
                labelled.info[key] = value
        return labelled

    def sorted_by_parameter(self, reverse=False):
        """
        In increasing order of the parameter.

        Worth doing explicitly, because ``from_files`` sorts paths as *text*:
        ``-120mV`` sorts before ``-20mV``, and a titration loaded in that order
        plots as a scribble and fits as nonsense. Spectra without a parameter
        are an error here rather than being quietly dropped.
        """
        parameters = self.parameters
        if np.isnan(parameters).any():
            raise ValueError(
                f"{int(np.isnan(parameters).sum())} of {len(self)} spectra "
                f"have no parameter, so the collection cannot be put in order "
                f"along one."
            )
        order = np.argsort(parameters)
        if reverse:
            order = order[::-1]
        return self._like([self._spectra[i] for i in order])

    def to_dataframe(self, orientation='wide'):
        """
        As a pandas DataFrame. Requires pandas, deliberately not a dependency
       .

        ``orientation='wide'`` gives one column per spectrum indexed by x;
        ``'long'`` gives tidy (sample, x, y) rows.
        """
        import pandas as pd  # pylint: disable=C0415

        x, matrix = self._stack()
        labels = [s.name for s in self._spectra]
        if orientation == 'wide':
            return pd.DataFrame(matrix.T, index=pd.Index(x, name='x'),
                                columns=labels)
        if orientation == 'long':
            return pd.DataFrame({
                'sample': np.repeat(self.samples, matrix.shape[1]),
                'name': np.repeat(labels, matrix.shape[1]),
                'x': np.tile(x, matrix.shape[0]),
                'y': matrix.ravel(),
            })
        raise ValueError("orientation must be 'wide' or 'long'")

    def plot(self, ax=None, **kwargs):
        """
        Overlay every spectrum on one axes, with a legend.

        See also ``spectroscopy.viz.stack`` for offset traces and
        ``spectroscopy.viz.grid`` for one panel per sample, which read
        better than an overlay past a handful of spectra.
        """
        from spectroscopy import viz  # pylint: disable=C0415
        return viz.plot_collection(self, ax, **kwargs)
