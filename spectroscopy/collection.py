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
from collections.abc import Sequence

import numpy as np

from spectroscopy.history import ProcessingStep
from spectroscopy.spectra import Spectrum

__all__ = ['SpectrumCollection']


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
                   sample_from=None, sort=True, **kwargs):
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

    def to_matrix(self):
        """
        ``(x, X)`` with ``X`` of shape (n_spectra, n_points).

        This is the handover to PCA/NMF/ICA, and the reason those notebooks no
        longer need to assemble ``np.array(aligned)`` by hand.
        """
        return self._stack()

    @property
    def samples(self):
        """Sample name of each spectrum, in order."""
        return [s.metadata.get('sample') for s in self._spectra]

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

    def plot(self, ax, *args, **kwargs):
        """Plot every spectrum on one axes."""
        for spectrum in self._spectra:
            spectrum.plot(ax, *args, **kwargs)
        return ax
