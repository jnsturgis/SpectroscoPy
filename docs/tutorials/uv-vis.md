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

# UV-Vis: comparing membrane fractions

Two fractions from a bacterial membrane preparation, measured 600–950 nm. Both
contain a light-harvesting complex with its characteristic pair of
near-infrared bands, and the question is how the two fractions differ.

The example spectra ship with the package, so every cell runs.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz

fractions = spc.datasets.load_pair()
print(fractions)

fig, ax = plt.subplots(figsize=(7, 3))
fractions.plot(ax)
```

Two bands, near 800 and 850 nm. Their ratio is the thing worth measuring.

## A word about units

These files are headerless CSV — two columns of numbers and nothing else. That
means the file cannot say what its y axis is, so the library falls back to the
technique's default and calls it absorbance:

```{code-cell} ipython3
first = fractions[0]
print(f"y unit:  {first.y_unit}")
print(f"y range: {first.y.min():.1f} to {first.y.max():.1f}")
```

A peak absorbance of 143 is not absorbance — these are almost certainly
milli-absorbance units, so about 0.14 AU, which is a sensible optical density
for a membrane fraction.

:::{admonition} Why the library makes a fuss about units
:class: note

Nothing here is wrong, but nothing here is *known* either. A JCAMP or `.spy`
file states its units and the library keeps what it is told, in preference to
any technique default. A bare two-column CSV cannot, so you are trusting a
guess. It is worth setting it straight yourself:

```python
spectrum.y_unit = 'absorbance'
spectrum.y = spectrum.y / 1000        # mAU -> AU
```

The analysis below normalises, so the absolute scale drops out — but that is
luck, not safety.
:::

## Correcting the baseline

Turbid samples scatter, which lifts the whole spectrum and tilts it towards the
blue. These fractions are fairly clear, so the offset is small — but taking it
off still matters, because a band *ratio* is sensitive to any pedestal
underneath both bands.

Guide points are the natural tool: name the wavelengths where you know there is
nothing absorbing, and a polynomial is fitted through the spectrum's own values
there.

Choosing them is the whole job, and it is easy to get wrong. Look for regions
that are genuinely flat before naming any:

```{code-cell} ipython3
first = fractions[0]
for low in range(600, 950, 50):
    window = first.crop(low, low + 50)
    print(f"  {low}-{low+50} nm: median {np.median(window.y):7.2f}  "
          f"spread {window.y.std():5.2f}")
```

Only 625–700 nm and 925–950 nm are flat. Everything from 700 nm upward is
already on the rising edge of the bands, and the extreme blue end is noisy.

```{code-cell} ipython3
guides = [630, 650, 670, 690, 930, 940, 950]

baseline = first.baseline('poly', degree=1, points=guides, halfwidth=3)
corrected = fractions.map(
    lambda s: s.baseline_correct('poly', degree=1, points=guides, halfwidth=3))

fig, ax = plt.subplots(figsize=(7, 3))
viz.plot_baseline(first, baseline, ax)
```

:::{admonition} A guide point in the wrong place ruins the result
:class: warning

Writing this page I first used `[600, 620, 640, 660, 900, 930, 950]` — round
numbers that looked "away from both bands". But 900 nm is still on the shoulder
of the 850 band (median 5.4) and 600 nm is a noisy rising edge. Those two
points dragged a degree-2 polynomial up to +2…+5 while the true baseline sits
at −2, so the correction pushed the whole spectrum *further* negative.

The library now warns when a subtractive baseline leaves the result mostly
negative, which is the signature of exactly this mistake. `halfwidth=3` also
helps: it takes the median over a few points either side of each guide instead
of trusting a single one.
:::

:::{admonition} If your sample really is turbid
:class: tip

Scattering falls off roughly as λ⁻⁴, so it rises steeply towards the blue. A
degree-3 or -4 polynomial through guide points either side of the bands handles
it; check the result by eye, because a polynomial that is too flexible will
happily eat the edge of a real band.
:::

## Smoothing

```{code-cell} ipython3
smooth = corrected.map(lambda s: s.smooth('savgol', window_length=9,
                                          polyorder=3))

fig, ax = plt.subplots(figsize=(7, 3))
viz.plot(corrected[0].crop(780, 880), ax, label='raw')
viz.plot(smooth[0].crop(780, 880), ax, label='smoothed')
ax.legend(frameon=False)
```

Use as little as you can get away with. Smoothing broadens bands and moves
shoulders, and a band ratio taken from an over-smoothed spectrum is a
measurement of the filter as much as of the sample.

## Normalising and comparing

Normalising to the stronger band puts both fractions on the same footing, so
what is left is the difference in composition rather than in concentration.

```{code-cell} ipython3
normalised = smooth.map(lambda s: s.normalize('max', window=(840, 860)))

fig, ax = plt.subplots(figsize=(7.5, 3.6))
viz.plot_collection(normalised, ax, labels=['fraction 1', 'fraction 2'],
                    frame=True)
viz.annotate_bands(ax, {800: "B800", 850: "B850"}, y=1.05)
ax.set_ylim(-0.05, 1.25)
```

The second fraction carries relatively more of the 800 nm band.

## Putting a number on it

```{code-cell} ipython3
def band_ratio(spectrum, first=(790, 810), second=(840, 860)):
    """Peak height in one window divided by peak height in another."""
    return (spectrum.crop(*first).y.max() /
            spectrum.crop(*second).y.max())

for name, spectrum in zip(['fraction 1', 'fraction 2'], normalised):
    print(f"  {name}: B800/B850 = {band_ratio(spectrum):.3f}")
```

A ratio like this is the usual quantitative readout, and it is worth writing as
a two-line function rather than reading numbers off a plot: it is reproducible,
and it does not change when you redraw the figure.

## Changing the x axis

Absorption bands are symmetric in energy, not in wavelength, so band-shape work
is often done in wavenumbers. One call, and the y values follow the reordering:

```{code-cell} ipython3
in_wavenumbers = normalised[0].to('cm^-1')

fig, ax = plt.subplots(figsize=(7, 3))
in_wavenumbers.plot(ax)
print(f"800 nm and 850 nm are {1e7/800:.0f} and {1e7/850:.0f} cm-1 "
      f"-- {1e7/800 - 1e7/850:.0f} cm-1 apart")
```

## The record

```{code-cell} ipython3
print(normalised[0].describe_history())
```

## What to change for your own data

| line | change to |
|---|---|
| `spc.datasets.load_pair()` | `spc.SpectrumCollection.from_files("fractions/*.csv", technique="UV-Vis")` |
| `guides` | wavelengths where your sample does not absorb |
| `window=(840, 860)` | your reference band |
| `band_ratio(...)` windows | the two bands you are comparing |
