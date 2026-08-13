---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Processing

Every operation returns a new spectrum and records itself. All of them exist as
plain array functions in `spectroscopy.processing.common` too, if you would
rather work with numpy directly.

Everything on this page runs when the documentation is built.

```{code-cell}
:tags: [hide-cell]

import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz
from spectroscopy.lineshapes import gauss

rng = np.random.default_rng(1)
x = np.linspace(1400, 1800, 801)

def make(amount=1.0, seed_noise=0.005):
    """An amide I/II pair on a sloping background, with noise."""
    y = (gauss(x, 1652, 32, amount) + gauss(x, 1548, 30, 0.6 * amount)
         + 0.10 + 1.2e-4 * (x - 1400)
         + seed_noise * rng.normal(size=x.size))
    return spc.Spectrum(x, y, technique='ATR-FTIR', name='sample')

spectrum = make()
```

## Cropping

```{code-cell}
spectrum.crop(1500, 1750)
```

Order-insensitive, so `crop(1750, 1500)` is the same — infrared ranges get
quoted both ways.

```{code-cell}
a = spectrum.crop(1500, 1750)
b = spectrum.crop(1750, 1500)
print("same result:", np.array_equal(a.x, b.x) and np.array_equal(a.y, b.y))
```

## Baselines

`baseline()` returns the baseline itself; `baseline_correct()` returns the
corrected spectrum and records the correction. Both take the same arguments.

```{code-cell}
fig, ax = plt.subplots()
viz.plot_baseline(spectrum, spectrum.baseline('rubberband'), ax);
```

```{code-cell}
corrected = spectrum.baseline_correct('rubberband')
print(corrected.history[-1])
```

**rubberband** — the convex hull. No parameters, which is its virtue, but it
can only ever be convex, so it cannot follow a background curving the other
way. `return_points=True` gives back the anchors, which can be adjusted and fed
to the polynomial method.

**poly** — a polynomial through guide points: name the x positions where you
know nothing absorbs. `halfwidth=3` takes a median around each rather than
trusting one reading. `coefficients=` uses a polynomial you already know.

Guide points can also be **(x, y) pairs**, giving the baseline value at each
position rather than reading it off the spectrum. That is the form to reach
for when the spectrum never actually reaches its own baseline — with bands
overlapping across the whole range, every scalar guide point sits on a
shoulder and drags the fit upward:

```{code-cell}
crowded_y = (gauss(x, 1650, 120, 1.0) + gauss(x, 1400, 150, 0.8)
             + gauss(x, 1100, 140, 0.7) + 0.15 + 2.0e-4 * (x - 1400))
crowded = spc.Spectrum(x, crowded_y, technique='ATR-FTIR', name='crowded')
truth = 0.15 + 2.0e-4 * (x - 1400)

from_positions = crowded.baseline('poly', degree=1, points=[1400, 1450, 1750, 1800])
from_pairs = crowded.baseline('poly', degree=1,
                              points=[(1400, 0.15), (1600, 0.19), (1800, 0.23)])

for name, estimate in (('x positions', from_positions), ('(x, y) pairs', from_pairs)):
    print(f"{name:<14} mean error {np.mean(np.abs(estimate.y - truth)):.4f}")
```

```{code-cell}
fig, ax = plt.subplots()
crowded.plot(ax, label='spectrum')
ax.plot(x, truth, 'k:', label='true background')
ax.plot(x, from_positions.y, label='from x positions')
ax.plot(x, from_pairs.y, label='from (x, y) pairs')
ax.legend(fontsize='small');
```

It also lets a background measured somewhere else — a buffer, a blank window,
a neighbouring sample — be imposed rather than guessed. Either way the anchors
land in the history, so the figure stays traceable to the numbers that made
it.

**als** — asymmetric least squares. `lam` is stiffness, `p` asymmetry
(0.001–0.02). The right choice when the background is not convex.

```{code-cell}
fig, ax = plt.subplots()
spectrum.plot(ax, label='raw')
for method, kwargs in (('rubberband', {}),
                       ('poly', dict(degree=2, points=[1400, 1450, 1750, 1800])),
                       ('als', dict(lam=1e5, p=0.01))):
    ax.plot(x, spectrum.baseline(method, **kwargs).y, lw=1, label=method)
ax.legend();
```

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

```{code-cell}
noisy = make(seed_noise=0.03).baseline_correct('rubberband')

fig, ax = plt.subplots()
noisy.plot(ax, label='raw', lw=0.8)
noisy.smooth('savgol', window_length=15, polyorder=3).plot(ax, label='smoothed')
ax.legend();
```

Even window lengths are rounded up to odd. This matters: an even-length
Savitzky–Golay window is asymmetric about its centre and biases every recovered
peak position by half a sample.

```{code-cell}
even = noisy.smooth('savgol', window_length=14, polyorder=3)
odd = noisy.smooth('savgol', window_length=15, polyorder=3)
print("14 was treated as 15:", np.allclose(even.y, odd.y))
```

Use as little smoothing as you can. It broadens bands and moves shoulders.

### Putting a spectrum on different x values

`resample` works out what a spectrum would read at positions you did not
measure. You need it before arithmetic between spectra recorded at different
settings, and every fit against a reference set does it for you.

Two methods. `spline` is the default: a cubic through every point, so it
reproduces the values you have exactly — and therefore reproduces the noise
exactly, since a spline is obliged to pass through it. `savgol` fits a
polynomial through a window of neighbours instead, so it averages noise down
as it interpolates, and you supply the window and the order.

The measured difference, against a Gaussian band of known shape, in root-mean-
square error over the interpolated points:

| points across a band | signal-to-noise | `savgol` window=5, order=3 | best window | that window's gain |
|---|---|---|---|---|
| 4 | 1000 | 1.16× better than spline | 5 (1.2 band widths) | 1.2× |
| 10 | 300 | 1.20× | 11 (1.1 band widths) | 2.3× |
| 20 | 100 | 1.17× | 31 (1.6 band widths) | 3.3× |
| 10 | 20 | 1.20× | 21 (2.1 band widths) | 3.1× |
| 20 | 10 | 1.16× | 41 (2.0 band widths) | 4.1× |

Three things to take from it.

**With real noise, savgol always wins.** The spline is better only on data with
no noise in it, which no instrument produces. This is why the guidance is not
"use the default and stop thinking".

**A five-point window is worth a flat 20%, whatever the noise.** Five points
fitting a cubic leaves one degree of freedom to average with, and one degree of
freedom is worth about a fifth and no more — at signal-to-noise of 1000 or of
10 alike. It is a safe choice precisely because it barely does anything.

**The real gain needs a window matched to the band**, one to four times the
number of points across it, wider when noisier — and that is worth up to four
times. Which is why the window is asked for and never guessed: it depends on
your sampling interval against your band width, and neither is recoverable from
the array of numbers.

So the default stays the spline, on the grounds that re-indexing a spectrum
onto a common axis is not the moment to be quietly denoising it by a fifth. If
you want the noise down, say so, with a window you chose — that is `smooth`,
above, and the same window guidance applies.

:::{admonition} If you really wanted it automatic
:class: note

The band width is in principle recoverable: take the power spectrum of the
spectrum and fit a Gaussian to it, and the narrowest real band falls out of the
high-frequency edge. It is a good technique and it is not implemented here,
because for ordinary re-indexing the 20% at stake does not pay for a
transform, a fit, and a new way to be wrong. It would earn its place in a
routine that has to process spectra it has never seen, unattended.
:::

```{code-cell}
fig, ax = plt.subplots()
corrected.derivative(order=2).plot(ax);
```

## Normalising

```{code-cell}
scaled = corrected.normalize('max', window=(1600, 1700))
print("max in window:", round(float(scaled.crop(1600, 1700).y.max()), 3))
```

`max`, `area`, `minmax`, `vector`, `snv`. The window restricts where the factor
is computed — normally a band the samples have in common.

## Peaks

```{code-cell}
peaks = corrected.find_peaks(prominence=0.10, relative=True, distance=40)
print(peaks)
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

```{code-cell}
print("absolute prominence=0.05 :", len(corrected.find_peaks(prominence=0.05)))
print("relative prominence=0.05 :",
      len(corrected.find_peaks(prominence=0.05, relative=True)))
```

You get a `PeakTable`: `.position`, `.height`, `.index`, `.prominence`,
`.width`, plus `.strongest(n)`, `.within(low, high)`, `.sorted_by_position()`,
`.annotate(ax)` and `.to_dataframe()`.

### Which way the bands point

Bands in transmittance, `%T` and reflectance are **minima**. `find_peaks` takes
the direction from the y unit, so those spectra are searched downward without
being asked, and `.kind` comes back as `'trough'`.

```{code-cell}
band = np.exp(-0.5 * ((x - 1650) / 20.0) ** 2)          # one band at 1650
absorbance = spc.Spectrum(x, 0.8 * band, x_unit='cm^-1',
                          y_quantity='Absorbance', y_unit='absorbance')
transmittance = spc.Spectrum(x, 1.0 - 0.8 * band, x_unit='cm^-1',
                             y_quantity='Transmittance', y_unit='transmittance')

for spec in (absorbance, transmittance):
    found = spec.find_peaks(prominence=0.05, relative=True)
    print(f"{spec.y_unit:<14} {np.round(found.position, 1)}  "
          f"kind={found.kind}  from={found.properties['direction_from']}")
```

Getting this wrong is silent rather than loud. Searching a transmission
spectrum *upward* finds the two inflection points that flank each band instead
of the band itself:

```{code-cell}
wrong = transmittance.find_peaks(prominence=0.05, relative=True, troughs=False)
print("forced upward:", np.round(wrong.position, 1), "— neither is the band")
```

A plausible number of plausible-looking positions, none of which is a band.
`.properties['direction_from']` says where the direction came from: `'y_unit'`,
`'caller'`, or `'assumed'` for units like `a.u.` and `counts` that do not fix a
direction. Those keep the upward default, which is right for the Raman and
fluorescence cases where they are usual.

Pass `troughs=` to decide for yourself. **A difference spectrum needs this** —
bands go both ways and the unit cannot know which you meant:

```{code-cell}
gained_band = np.exp(-0.5 * ((x - 1660) / 18.0) ** 2)
lost_band = np.exp(-0.5 * ((x - 1630) / 18.0) ** 2)
difference = spc.Spectrum(x, 0.4 * gained_band - 0.3 * lost_band,
                          x_unit='cm^-1', y_unit='absorbance',
                          y_quantity='Δ Absorbance')

gained = difference.find_peaks(prominence=0.2, relative=True, troughs=False)
lost = difference.find_peaks(prominence=0.2, relative=True, troughs=True)
print("appeared   :", np.round(gained.position, 1))
print("disappeared:", np.round(lost.position, 1))
```

:::{admonition} Amounts need absorbance, not transmittance
:class: warning

Positions are fine in either unit, but heights, areas and band **ratios** are
not. Beer–Lambert is linear in absorbance, so an area measured on `%T` is not
proportional to concentration — a ratio taken there is not a ratio of anything.
`fit_peaks` warns when you fit a valley-pointing unit:
:::

```{code-cell}
import warnings

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    transmittance.fit_peaks(n_peaks=1)
    print(caught[0].message)
```

Convert first, and the warning goes away:

```{code-cell}
fit = transmittance.to(y_unit='absorbance').fit_peaks(n_peaks=1)
print("fitted at", np.round(fit.position, 1), " R^2 =", round(fit.r_squared, 4))
```

The warning is deliberate rather than a silent conversion: quietly changing
what a fit was performed on is its own kind of wrong.

## Reference subtraction

```{code-cell}
buffer = spc.Spectrum(x, 0.30 * gauss(x, 1640, 60, 1.0) + 0.02,
                      technique='ATR-FTIR', name='buffer')
sample = spc.Spectrum(x, corrected.y + 0.995 * buffer.y,
                      technique='ATR-FTIR', name='sample in buffer')

cleaned = sample.subtract_reference(buffer, factor=0.995)
print(cleaned.history[-1])
```

The factor is an argument so it lands in the history. It is usually the one
number in an analysis chosen by eye, and the one most worth recording.

## Multivariate

Needs `pip install "spectroscopy[multivariate]"`.

```{code-cell}
:tags: [hide-cell]

components = [gauss(x, 1652, 32, 1.0), gauss(x, 1630, 26, 1.0),
              gauss(x, 1548, 30, 1.0)]
mixtures = [(1.00, 0.20, 0.60), (0.80, 0.45, 0.50), (0.50, 0.80, 0.35),
            (0.95, 0.25, 0.58), (0.75, 0.50, 0.48), (0.45, 0.85, 0.33)]

series = []
for index, weights in enumerate(mixtures):
    y = sum(w * c for w, c in zip(weights, components))
    y = y + 0.10 + 1.2e-4 * (x - 1400) + 0.006 * rng.normal(size=x.size)
    member = spc.Spectrum(x, y, technique='ATR-FTIR', name=f"sample {index + 1}")
    member.set_sample("ABC"[index % 3])
    series.append(member)

collection = (spc.SpectrumCollection(series, name='mixtures')
              .baseline_correct('rubberband').crop(1500, 1750))
```

```{code-cell}
from spectroscopy.processing import multivariate as mv

result = mv.decompose(collection, 'nmf', n_components=3)
print(result.components)
print(np.round(result.contributions(), 3))
```

Components come back as spectra carrying the collection's axis and units, so
everything else in the library works on them — including peak picking:

```{code-cell}
for component in result.components:
    found = component.find_peaks(prominence=0.2, relative=True)
    print(component.name, "->", np.round(found.position, 0))
```

:::{admonition} Check stability before believing a component count
:class: warning

Explained variance rises monotonically with k, so it can never tell you when to
stop. Resampling the samples can. Above ~0.95 you have components worth
discussing; below ~0.8 the count is probably wrong.
:::

```{code-cell}
print(mv.stability(collection, 'nmf', n_components=3, n_trials=5))
```

Reseeding (`mode='runs'`) is usually uninformative for NMF — the initialisation
is deterministic, so it returns 1.0 whatever k is.
