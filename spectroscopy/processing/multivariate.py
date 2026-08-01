# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Decomposition of a set of spectra: PCA, NMF, ICA -- with stability testing.

This is the analysis in ``Biofilm_CK_241125`` / ``Figure_CK_111225`` /
``Figures_CK``, which was the newest and most developed work in the notebooks
and lived entirely outside the library, operating on numpy matrices assembled
by hand. The stability testing in particular -- repeated fits, bootstrap
resampling, and Hungarian matching of components between runs -- is a genuinely
reusable piece of method work that was trapped in one notebook.

Two things it adds beyond wrapping scikit-learn:

* **Components come back as spectra.** An NMF component *is* a spectrum, so it
  arrives as a :class:`~spectroscopy.spectra.Spectrum` with the collection's x
  axis and units, and can be peak-picked, plotted and saved like any other.
* **Stability is a first-class question.** NMF has no unique solution; a
  k-component fit that moves under reseeding or resampling is not a finding.
  :func:`stability` answers that in one call.

scikit-learn is an optional dependency -- ``pip install spectroscopy[multivariate]``
(review section 5.6). It is imported lazily so the rest of the library does not
require it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = [
    'DecompositionResult', 'StabilityResult',
    'decompose', 'stability', 'match_components', 'scan_components',
]

METHODS = ('pca', 'nmf', 'ica')


def _sklearn():
    """Import scikit-learn, with an error that says how to get it."""
    try:
        from sklearn.decomposition import NMF, PCA, FastICA  # pylint: disable=C0415
    except ImportError as error:                             # pragma: no cover
        raise ImportError(
            "Multivariate analysis needs scikit-learn, which is an optional "
            "dependency. Install it with:\n"
            "    pip install 'spectroscopy[multivariate]'\n"
            "or\n"
            "    pip install scikit-learn"
        ) from error
    return {'pca': PCA, 'nmf': NMF, 'ica': FastICA}


def _matrix(source):
    """Accept a SpectrumCollection or a plain (x, X) pair."""
    if hasattr(source, 'to_matrix'):
        x, matrix = source.to_matrix()
        names = [s.name for s in source]
        samples = list(source.samples)
        template = source[0]
        return np.asarray(x, float), np.asarray(matrix, float), names, samples, template
    x, matrix = source
    matrix = np.asarray(matrix, float)
    names = [f"sample {i}" for i in range(matrix.shape[0])]
    return np.asarray(x, float), matrix, names, list(names), None


def _as_spectra(x, components, template, method):
    """Wrap component rows as Spectrum objects carrying the collection's axis."""
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415
    from spectroscopy.history import ProcessingStep  # pylint: disable=C0415
    from spectroscopy.spectra import Spectrum  # pylint: disable=C0415

    spectra = []
    for index, row in enumerate(components):
        spectrum = Spectrum(template) if template is not None else Spectrum()
        spectrum.x, spectrum.y = np.copy(x), np.copy(row)
        spectrum.name = f"{method.upper()} component {index + 1}"
        spectrum.metadata = {'component': index + 1, 'method': method}
        spectrum.history = [ProcessingStep(
            "decompose", {"method": method, "component": index + 1})]
        spectra.append(spectrum)
    return SpectrumCollection(spectra, name=f"{method.upper()} components")


@dataclass
class DecompositionResult:
    """
    The outcome of a decomposition.

    Attributes
    ----------
    components : SpectrumCollection
        The component spectra -- ``H`` for NMF, the loadings for PCA.
    weights : ndarray, shape (n_samples, n_components)
        How much of each component each sample carries -- ``W`` for NMF,
        the scores for PCA.
    explained_variance : float
        Fraction of the total variance the reconstruction accounts for.
    explained_variance_ratio : ndarray or None
        Per-component, where the method defines it (PCA).
    r2_per_sample : ndarray
        Reconstruction quality per sample; useful for spotting one bad spectrum.
    """

    method: str
    n_components: int
    components: object
    weights: np.ndarray
    x: np.ndarray
    reconstruction: np.ndarray
    residuals: np.ndarray
    explained_variance: float
    r2_per_sample: np.ndarray
    sample_names: list
    samples: list
    explained_variance_ratio: np.ndarray | None = None
    model: object = field(default=None, repr=False)
    _data: np.ndarray = field(default=None, repr=False)

    def __len__(self) -> int:
        return self.n_components

    def __repr__(self) -> str:
        return (f"<DecompositionResult {self.method.upper()} k={self.n_components}, "
                f"{len(self.sample_names)} samples, "
                f"{100 * self.explained_variance:.1f}% explained>")

    def contributions(self, normalize=True):
        """
        Fraction of each sample's spectral area coming from each component.

        This is the table the CK notebooks build and then run statistics on.
        Each component spectrum is integrated to give its area, the sample's
        weight is scaled by it, and the row is normalised to sum to 1.

        Returns
        -------
        ndarray, shape (n_samples, n_components)
        """
        spacing = float(np.mean(np.diff(self.x))) if len(self.x) > 1 else 1.0
        areas = np.array([np.abs(spectrum.y).sum() for spectrum in self.components])
        areas = areas * abs(spacing)
        scaled = self.weights * areas[np.newaxis, :]
        if not normalize:
            return scaled
        totals = scaled.sum(axis=1, keepdims=True)
        with np.errstate(invalid='ignore', divide='ignore'):
            return np.where(totals != 0, scaled / totals, np.nan)

    def reconstructed(self, index):
        """Sample ``index`` as rebuilt from the components, as a Spectrum."""
        return self._spectrum(self.reconstruction[index],
                              f"{self.sample_names[index]} reconstructed")

    def residual(self, index):
        """Observed minus reconstructed for sample ``index``, as a Spectrum."""
        return self._spectrum(self.residuals[index],
                              f"{self.sample_names[index]} residual")

    def _spectrum(self, values, name):
        from spectroscopy.spectra import Spectrum  # pylint: disable=C0415
        template = self.components[0] if len(self.components) else None
        spectrum = Spectrum(template) if template is not None else Spectrum()
        spectrum.x, spectrum.y = np.copy(self.x), np.copy(values)
        spectrum.name = name
        spectrum.metadata = {}
        spectrum.history = []
        return spectrum

    def to_dataframe(self):
        """
        Contributions as a pandas DataFrame, indexed by sample name.

        Requires pandas, which is deliberately not a dependency -- see review
        section 5.6. The arrays above are the supported interface.
        """
        import pandas as pd  # pylint: disable=C0415
        return pd.DataFrame(
            self.contributions(),
            index=pd.Index(self.sample_names, name='name'),
            columns=[f"Comp_{i + 1}" for i in range(self.n_components)],
        )


@dataclass
class StabilityResult:
    """How reproducible a decomposition is under reseeding or resampling."""

    method: str
    n_components: int
    mode: str
    n_trials: int
    correlations: np.ndarray            # (n_trials, n_components)
    mean_per_component: np.ndarray
    std_per_component: np.ndarray

    @property
    def overall(self) -> float:
        """Mean matched correlation across all components."""
        return float(np.mean(self.mean_per_component))

    def __repr__(self) -> str:
        return (f"<StabilityResult {self.method.upper()} k={self.n_components} "
                f"{self.mode}: overall {self.overall:.3f}>")

    def __str__(self) -> str:
        lines = [f"{self.method.upper()} k={self.n_components}, {self.mode}, "
                 f"{self.n_trials} trials"]
        for index, (mean, spread) in enumerate(
                zip(self.mean_per_component, self.std_per_component), start=1):
            lines.append(f"  component {index}: {mean:.3f} +/- {spread:.3f}")
        lines.append(f"  overall: {self.overall:.3f}")
        return "\n".join(lines)


def _check_non_negative(matrix, method):
    """
    NMF needs non-negative input; sklearn's own error does not say what to do.

    Baseline-corrected spectra routinely dip slightly below zero, so this is
    the first thing that goes wrong when someone tries NMF on real data.
    """
    if method != 'nmf':
        return matrix
    smallest = float(matrix.min())
    if smallest >= 0:
        return matrix
    raise ValueError(
        f"NMF needs non-negative data but the most negative value is "
        f"{smallest:.4g}. Baseline-correct first, or clip with "
        f"numpy.clip(X, 0, None) if the negatives are just noise about zero."
    )


def _fit(method, matrix, n_components, random_state=None, **kwargs):
    """Fit one model and return (components, weights, reconstruction, model)."""
    estimators = _sklearn()
    method = method.lower()
    if method not in METHODS:
        raise ValueError(f"Unknown method {method!r}; try one of {METHODS}")

    if method == 'pca':
        model = estimators['pca'](n_components=n_components,
                                  random_state=random_state, **kwargs)
        centre = matrix.mean(axis=0, keepdims=True)
        weights = model.fit_transform(matrix - centre)
        components = model.components_
        reconstruction = weights @ components + centre
    elif method == 'nmf':
        options = {'init': 'nndsvda', 'max_iter': 2000}
        options.update(kwargs)
        model = estimators['nmf'](n_components=n_components,
                                  random_state=random_state, **options)
        weights = model.fit_transform(matrix)
        components = model.components_
        reconstruction = weights @ components
    else:
        options = {'max_iter': 2500}
        options.update(kwargs)
        model = estimators['ica'](n_components=n_components,
                                  random_state=random_state, **options)
        weights = model.fit_transform(matrix)
        components = model.mixing_.T
        reconstruction = weights @ components + model.mean_

    return components, weights, reconstruction, model


def decompose(source, method='nmf', n_components=3, random_state=0, **kwargs):
    """
    Decompose a collection of spectra into ``n_components`` component spectra.

    Parameters
    ----------
    source : SpectrumCollection or (x, X)
        The spectra. A collection is the usual case; the raw pair is accepted
        so the function can be used on a matrix you already have.
    method : {'nmf', 'pca', 'ica'}
    n_components : int
    random_state : int or None
        Fixed by default, because NMF and ICA are seed-dependent and a result
        that changes between runs is not one you can put in a paper. Use
        :func:`stability` to find out how much it depends on the seed.

    Returns
    -------
    DecompositionResult
    """
    x, matrix, names, samples, template = _matrix(source)
    method = method.lower()
    matrix = _check_non_negative(matrix, method)

    components, weights, reconstruction, model = _fit(
        method, matrix, n_components, random_state, **kwargs)

    residuals = matrix - reconstruction
    total = np.sum((matrix - matrix.mean()) ** 2)
    explained = 1.0 - np.sum(residuals ** 2) / total if total else float('nan')

    per_sample = []
    for observed, rebuilt in zip(matrix, reconstruction):
        spread = np.sum((observed - observed.mean()) ** 2)
        per_sample.append(1.0 - np.sum((observed - rebuilt) ** 2) / spread
                          if spread else float('nan'))

    return DecompositionResult(
        method=method,
        n_components=n_components,
        components=_as_spectra(x, components, template, method),
        weights=weights,
        x=x,
        reconstruction=reconstruction,
        residuals=residuals,
        explained_variance=float(explained),
        explained_variance_ratio=getattr(model, 'explained_variance_ratio_', None),
        r2_per_sample=np.array(per_sample),
        sample_names=names,
        samples=samples,
        model=model,
        _data=matrix,
    )


def match_components(reference, target):
    """
    Pair up two sets of component spectra by correlation.

    A decomposition returns its components in arbitrary order and sign, so
    comparing two fits means solving an assignment problem first. Uses the
    Hungarian algorithm on the absolute correlation, keeping the sign in the
    reported value.

    Parameters
    ----------
    reference, target : ndarray, shape (n_components, n_points)

    Returns
    -------
    (rows, columns, correlations) : tuple of ndarray
    """
    from scipy.optimize import linear_sum_assignment  # pylint: disable=C0415
    from scipy.stats import pearsonr  # pylint: disable=C0415

    reference = np.atleast_2d(np.asarray(reference, float))
    target = np.atleast_2d(np.asarray(target, float))

    matrix = np.zeros((reference.shape[0], target.shape[0]))
    for row in range(reference.shape[0]):
        for column in range(target.shape[0]):
            matrix[row, column] = pearsonr(reference[row], target[column])[0]

    rows, columns = linear_sum_assignment(-np.abs(matrix))
    return rows, columns, matrix[rows, columns]


def _component_matrix(result):
    return np.vstack([spectrum.y for spectrum in result.components])


def stability(source, method='nmf', n_components=3, mode='bootstrap',
              n_trials=20, random_seed=0, **kwargs):
    """
    Measure how reproducible a decomposition is.

    NMF and ICA have no unique solution, so a k-component fit that moves when
    you drop a sample is not a finding. This refits repeatedly and reports the
    matched correlation of each component against a reference fit.

    Parameters
    ----------
    mode : {'bootstrap', 'runs'}
        ``bootstrap`` (the default) resamples the spectra with replacement --
        *does the answer depend on which samples I happen to have?* That is
        the scientific question, and on real data it is the probe that
        discriminates.

        ``runs`` refits the same data with different seeds -- *is the
        algorithm itself reproducible?* See the note below before relying
        on it.
    n_trials : int

    Returns
    -------
    StabilityResult

    Notes
    -----
    Mean matched correlations above ~0.95 indicate components you can talk
    about; below ~0.8 the number of components is probably wrong.

    **``mode='runs'`` is often uninformative for NMF.** The default
    initialisation is deterministic and coordinate descent converges to the
    same optimum, so reseeding returns an identical answer and the score is
    exactly 1.0 whatever ``k`` is. This function switches to the perturbed
    ``nndsvdar`` initialisation to give the seed something to change, but on
    well-conditioned data even that lands in the same place. Measured on the
    23 real ATR-FTIR spectra of the biofilm series:

    ======  ===============  ==================
    k       runs             bootstrap
    ======  ===============  ==================
    2       1.000            0.994
    4       1.000            0.939
    6       1.000            0.858
    8       1.000            0.826
    ======  ===============  ==================

    Only the bootstrap column separates a defensible ``k`` from an
    over-fitted one, which is why it is the default.
    """
    if mode not in ('runs', 'bootstrap'):
        raise ValueError(f"mode must be 'runs' or 'bootstrap', not {mode!r}")

    x, matrix, _, _, template = _matrix(source)
    method = method.lower()
    matrix = _check_non_negative(matrix, method)
    generator = np.random.RandomState(random_seed)

    # Reseeding a deterministic initialiser measures nothing. NMF's default
    # here (and in the notebooks) is 'nndsvda', which is fully determined by
    # the data -- so every "different seed" run returns the identical answer
    # and stability comes out as exactly 1.0 whatever k is. For the reseeding
    # probe, switch to 'nndsvdar', the same initialisation perturbed by noise,
    # so the seed genuinely matters. decompose() keeps the deterministic
    # default, which is what you want for a reproducible published fit.
    if mode == 'runs' and method == 'nmf' and 'init' not in kwargs:
        kwargs = {**kwargs, 'init': 'nndsvdar'}

    reference = decompose((x, matrix), method, n_components,
                          random_state=random_seed, **kwargs)
    reference_components = _component_matrix(reference)

    collected = []
    for _ in range(n_trials):
        seed = int(generator.randint(0, 2 ** 31 - 1))
        if mode == 'runs':
            trial_matrix = matrix
        else:
            rows = generator.randint(0, matrix.shape[0], size=matrix.shape[0])
            trial_matrix = matrix[rows, :]
        trial = decompose((x, trial_matrix), method, n_components,
                          random_state=seed, **kwargs)
        _, _, correlations = match_components(reference_components,
                                              _component_matrix(trial))
        collected.append(np.abs(correlations))

    collected = np.array(collected)
    _ = template
    return StabilityResult(
        method=method.lower(),
        n_components=n_components,
        mode=mode,
        n_trials=n_trials,
        correlations=collected,
        mean_per_component=collected.mean(axis=0),
        std_per_component=collected.std(axis=0),
    )


def scan_components(source, method='nmf', candidates=range(2, 8),
                    stability_trials=0, **kwargs):
    """
    Fit several component counts and report how each does.

    Choosing k is the awkward decision in this analysis; the notebooks do it by
    printing explained variance for a range and eyeballing it. This returns the
    same numbers as arrays, optionally with a stability figure -- which is the
    more informative criterion, since explained variance rises monotonically
    with k and so never tells you to stop.

    Returns
    -------
    dict with keys ``n_components``, ``explained_variance``,
    ``mean_r2``, and (when ``stability_trials`` > 0) ``stability``.
    """
    x, matrix, _, _, _ = _matrix(source)
    counts, explained, mean_r2, stabilities = [], [], [], []

    for count in candidates:
        result = decompose((x, matrix), method, count, **kwargs)
        counts.append(count)
        explained.append(result.explained_variance)
        mean_r2.append(float(np.nanmean(result.r2_per_sample)))
        if stability_trials:
            stabilities.append(stability((x, matrix), method, count,
                                         n_trials=stability_trials,
                                         **kwargs).overall)

    table = {
        'n_components': np.array(counts),
        'explained_variance': np.array(explained),
        'mean_r2': np.array(mean_r2),
    }
    if stability_trials:
        table['stability'] = np.array(stabilities)
    return table
