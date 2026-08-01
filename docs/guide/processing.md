# Processing

Every operation returns a new spectrum and records itself. All of them exist as
plain array functions in `spectroscopy.processing.common` too, if you would
rather work with numpy directly.

## Cropping

```python
spectrum.crop(900, 1800)
```

Order-insensitive, so `crop(1800, 900)` is the same — infrared ranges get
quoted both ways.

## Baselines

`baseline()` returns the baseline itself; `baseline_correct()` returns the
corrected spectrum and records the correction. Both take the same arguments.

```python
spectrum.baseline_correct('rubberband')
spectrum.baseline_correct('poly', degree=2, points=[600, 620, 900, 950])
spectrum.baseline_correct('als', lam=1e5, p=0.01)
```

**rubberband** — the convex hull. No parameters, which is its virtue, but it
can only ever be convex, so it cannot follow a background curving the other
way. `return_points=True` gives back the anchors, which can be adjusted and fed
to the polynomial method.

**poly** — a polynomial through guide points: name the x positions where you
know nothing absorbs. `halfwidth=3` takes a median around each rather than
trusting one reading. `coefficients=` uses a polynomial you already know.

**als** — asymmetric least squares. `lam` is stiffness, `p` asymmetry
(0.001–0.02). The right choice when the background is not convex.

:::{admonition} Which side is the baseline?
:class: note

In absorbance it is the *lower* hull, tending to zero. In transmittance it is
the *upper* hull, tending to 100%, and the correction is a division rather than
a subtraction. The library reads that off the y unit; `side=` and `mode=`
override.

If a subtractive correction leaves the result mostly negative you get a
warning: that means the baseline sat above the data, almost always because a
guide point was on a band shoulder rather than in a flat region.
:::

**Always look at the baseline before trusting it.** `viz.plot_baseline` exists
for this. Numeric checks pass on baselines that are visibly wrong.

## Smoothing and derivatives

```python
spectrum.smooth('savgol', window_length=15, polyorder=3)
spectrum.derivative(order=2)
```

Even window lengths are rounded up to odd. This matters: an even-length
Savitzky-Golay window is asymmetric about its centre and biases every recovered
peak position by half a sample.

Use as little smoothing as you can. It broadens bands and moves shoulders.

## Normalising

```python
spectrum.normalize('max', window=(1050, 1080))
```

`max`, `area`, `minmax`, `vector`, `snv`. The window restricts where the factor
is computed — normally a band the samples have in common.

## Peaks

```python
peaks = spectrum.find_peaks(prominence=0.05, relative=True, distance=15)
```

Works on the inverted second derivative by default, so shoulders on a broad
band are found as well as maxima. `method='direct'` uses the signal itself.

:::{admonition} Use relative=True
:class: tip

Without it, `height` and `prominence` apply to the *detection signal* — the
second derivative — which is orders of magnitude smaller than the spectrum. A
sensible-looking `prominence=0.05` then finds nothing and you have to guess
down to `0.00001`. With `relative=True` the number means a fraction of the
signal range.
:::

You get a `PeakTable`: `.position`, `.height`, `.index`, `.prominence`,
`.width`, plus `.strongest(n)`, `.within(low, high)`, `.sorted_by_position()`,
`.annotate(ax)` and `.to_dataframe()`.

## Reference subtraction

```python
corrected = sample.subtract_reference(buffer, factor=0.995)
```

The factor is an argument so it lands in the history. It is usually the one
number in an analysis chosen by eye, and the one most worth recording.

## Multivariate

Needs `pip install "spectroscopy[multivariate]"`.

```python
from spectroscopy.processing import multivariate as mv

result = mv.decompose(collection, 'nmf', n_components=4)
result.components          # a SpectrumCollection -- peak-pick them, plot them
result.contributions()     # (n_samples, n_components) area fractions
```

Components come back as spectra carrying the collection's axis and units, so
everything else in the library works on them.

:::{admonition} Check stability before believing a component count
:class: warning

```python
mv.stability(collection, 'nmf', n_components=4)      # bootstrap by default
```

Explained variance rises monotonically with k, so it can never tell you when to
stop. Resampling the samples can. Above ~0.95 you have components worth
discussing; below ~0.8 the count is probably wrong.

Reseeding (`mode='runs'`) is usually uninformative for NMF — the initialisation
is deterministic, so it returns 1.0 whatever k is.
:::
