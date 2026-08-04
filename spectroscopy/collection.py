# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
SpectrumCollection -- an ordered set of spectra with batch operations.

This is the answer to friction item #1 in the review: replicate averaging done
by hard-coded index. Across the notebooks that is written as

    PG_coli   = (spectra[0] + spectra[1] + spectra[2]) / 3.0
    Glucose   = (spectra[32] + spectra[33] + ... + spectra[37]) / 6.0
    Cellulose = (spectra[23] + spectra[24] + ... + spectra[27]) / 5.0

which breaks silently the moment a file is added, and which the notebooks
themselves flag: *"This should be simpler and neater call"* and *"Should be
based on sample name and list of numbers"*. With a collection it is

    averages = collection.group_by('sample').mean()

``to_matrix()`` is the bridge to the multivariate work: it returns plain numpy,
which is what scikit-learn wants, without pandas entering the core (review
section 5.6).
"""

from __future__ import annotations

import glob as _glob
import os
import re
from collections.abc import Sequence

import numpy as np

from spectroscopy.history import ProcessingStep
from spectroscopy.spectra import Spectrum

__all__ = ['SpectrumCollection']


def _parameter_reader(spec):
    """
    Turn the ``parameter_from`` argument into ``f(path) -> float``.

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
    """An ordered collection of :class:`~spectroscopy.spectra.Spectrum`."""

    def __init__(self, spectra=(), name=None):
        self._spectra = list(spectra)
        for item in self._spectra:
            if not isinstance(item, Spectrum):
                raise TypeError(
                    f"SpectrumCollection takes Spectrum objects, got {type(item).__name__}"
                )
        self.name = name

    # -- Sequence protocol -------------------------------------------------

    def __len__(self) -> int:
        return len(self._spectra)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return SpectrumCollection(self._spectra[index], name=self.name)
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
            return SpectrumCollection(self._spectra + list(other))
        if isinstance(other, Spectrum):
            return SpectrumCollection(self._spectra + [other])
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

            See :meth:`~spectroscopy.spectra.Spectrum.set_parameter` for why
            this is not just another ``sample``. Use
            :meth:`sorted_by_parameter` afterwards: ``sort=True`` orders the
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
        return cls(spectra)

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
        return {value: SpectrumCollection(members, name=str(value))
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
        return SpectrumCollection([function(s) for s in self._spectra],
                                  name=self.name)

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
        return SpectrumCollection([s for s in self._spectra if predicate(s)],
                                  name=self.name)

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
                f"from_files(parameter_from=...) or set_parameters()."
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

        ``nan`` where a spectrum has none -- unlike :meth:`to_matrix`, which
        refuses. Reading the values is a reasonable thing to do on a
        part-labelled collection; fitting them is not.
        """
        return np.array([s.metadata.get('parameter', np.nan)
                         for s in self._spectra], dtype=float)

    @property
    def parameter_name(self):
        """What the parameter is, if the spectra agree on it."""
        return self._parameter_label('parameter_name')

    @property
    def parameter_unit(self):
        """What the parameter is in, if the spectra agree on it."""
        return self._parameter_label('parameter_unit')

    def _parameter_label(self, key):
        labels = {s.metadata.get(key) for s in self._spectra} - {None}
        return labels.pop() if len(labels) == 1 else None

    def set_parameters(self, values, name=None, unit=None):
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
        return SpectrumCollection(spectra, name=self.name)

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
        return SpectrumCollection([self._spectra[i] for i in order],
                                  name=self.name)

    def to_dataframe(self, orientation='wide'):
        """
        As a pandas DataFrame. Requires pandas, deliberately not a dependency
        (review section 5.6).

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

        See also :func:`spectroscopy.viz.stack` for offset traces and
        :func:`spectroscopy.viz.grid` for one panel per sample, which read
        better than an overlay past a handful of spectra.
        """
        from spectroscopy import viz  # pylint: disable=C0415
        return viz.plot_collection(self, ax, **kwargs)
