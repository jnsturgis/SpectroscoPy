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

# Plotting

`spectroscopy.viz` holds the drawing. `Spectrum.plot`, `SpectrumCollection.plot`
and `PeakTable.annotate` delegate to it, so the short forms still work.

Every figure on this page is drawn when the documentation is built, so what you
see is what the code produces.

```{code-cell}
:tags: [hide-cell]

import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz
from spectroscopy.lineshapes import gauss

# Six ATR-FTIR-like spectra: three samples, two replicates each, on a
# sloping background. Amide I near 1650, amide II near 1548.
rng = np.random.default_rng(1)
x = np.linspace(1500, 1750, 501)

spectra = []
for index, (name, amount) in enumerate(
        [('control', 1.00), ('treated', 0.75), ('washed', 0.50)] * 2):
    y = (gauss(x, 1652, 32, amount) + gauss(x, 1548, 30, 0.6 * amount)
         + 0.10 + 1.2e-4 * (x - 1500) + 0.005 * rng.normal(size=x.size))
    spectrum = spc.Spectrum(x, y, technique='ATR-FTIR',
                            name=f"{name} {index // 3 + 1}")
    spectrum.set_sample(name)
    spectra.append(spectrum)

spectra = spc.SpectrumCollection(spectra, name='demo')
corrected = spectra.baseline_correct('rubberband').crop(1520, 1720)
```

## One spectrum

```{code-cell}
fig, ax = plt.subplots()
corrected[0].plot(ax);
```

The axis labels come from the spectrum, and the x axis is reversed for the
techniques that are conventionally drawn high-to-low — FTIR and Raman. Those
are the two lines typed by hand in every figure cell of a notebook.

`frame=True` draws a full box with inward ticks, as most physical-chemistry
journals expect. To make it the default for a session:

```{code-cell}
fig, ax = plt.subplots()
corrected[0].plot(ax, frame=True);
```

```python
from spectroscopy import viz
viz.FRAME_AXES = True          # for every plot from here on
```

## Several

```{code-cell}
fig, ax = plt.subplots()
viz.plot_collection(corrected, ax);
```

Overlay works to about half a dozen. Past that, `stack` or `grid` reads better
than any palette can.

```{code-cell}
fig, ax = plt.subplots(figsize=(7, 4.5))
viz.stack(corrected, ax);
```

`grid` is worth running on every new dataset before doing anything else: a
replicate that went wrong is obvious in a panel-per-sample overview and very
hard to spot after averaging.

```{code-cell}
fig, axes = viz.grid(corrected, key='sample')
```

## Annotation

```{code-cell}
peaks = corrected[0].find_peaks(prominence=0.10, relative=True, distance=40)
print(peaks)
```

Detection found five features; three of them are the bands worth discussing.
That gap is the normal case, not a failure — which is why the annotation step
takes `strongest(n)` rather than everything:

```{code-cell}
fig, ax = plt.subplots()
corrected[0].plot(ax)
peaks.strongest(3).annotate(ax, offset=0.03);
```

Label only the peaks you intend to discuss. Labelling every peak on a crowded
spectrum produces an overlapping smear that is worse than no labels — use
`strongest(n)` and `distance=` in the detection.

When you already know what the bands are, name them instead of numbering them:

```{code-cell}
fig, ax = plt.subplots()
corrected[0].plot(ax)
viz.annotate_bands(ax, {1652: "Amide I", 1548: "Amide II"});
```

A baseline should be looked at before it is trusted — numeric checks pass on
baselines that are visibly wrong:

```{code-cell}
raw = spectra[0]
fig, ax = plt.subplots()
viz.plot_baseline(raw, raw.baseline('rubberband'), ax);
```

## Colour

The default cycle is the six chromatic Okabe–Ito colours, ordered so adjacent
pairs stay separated under deuteranopia. The order was checked with a palette
validator rather than by eye; the obvious sequence puts pink next to green and
fails. Okabe–Ito's yellow and black are excluded — the yellow is too light
against white, the black has no chroma.

```{code-cell}
:tags: [hide-input]

fig, ax = plt.subplots(figsize=(7, 0.9))
for index, colour in enumerate(viz.PALETTE):
    ax.add_patch(plt.Rectangle((index, 0), 0.92, 1, color=colour))
    ax.text(index + 0.46, -0.28, colour, ha='center', va='top', fontsize=7)
ax.set_xlim(-0.1, len(viz.PALETTE)); ax.set_ylim(-0.6, 1)
ax.axis('off');
```

Colours are assigned in fixed order from the start of the list, never shuffled,
so a series keeps its colour when the set it is drawn with changes. Past six
the colours repeat but the line style advances, so the pair stays unique.

```{code-cell}
for index in (0, 1, 5, 6, 7):
    print(index, viz.series_style(index))
```

Three of the six sit below 3:1 contrast against white. That is why a legend is
always drawn for two or more series and why `stack` labels its traces directly:
identity never rests on colour alone.

For a **map** — an excitation-emission surface, say — use a sequential ramp
(one hue, light to dark, like `viridis`). A categorical palette encodes
identity; a map encodes magnitude. Never a rainbow: it invents boundaries where
the data is smooth and is unreadable in greyscale.

## Decomposition figures

```{code-cell}
from spectroscopy.processing import multivariate as mv

result = mv.decompose(corrected, 'nmf', n_components=3)
fig, axes = viz.plot_decomposition(result)
```

```{code-cell}
fig, ax = plt.subplots()
viz.plot_scores(result, ax);
```

## Saving

Ordinary matplotlib:

```python
fig.savefig("figure.pdf")      # vector, for a paper
fig.savefig("figure.png", dpi=300)
```

Nothing in this module prevents you doing anything matplotlib can do — every
function takes and returns axes, so drop down whenever you need to.
