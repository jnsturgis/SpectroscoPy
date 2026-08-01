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

# Fluorescence: an excitation-emission map

An emission spectrum recorded at one excitation wavelength tells you the shape
of the emission. Recorded at *many* excitation wavelengths it tells you which
absorbing species is responsible — and that is what an excitation-emission map
(EEM) is for.

This uses a real series recorded on a candidate flavoprotein: eighteen emission
spectra, excitation stepped from 290 to 455 nm.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz

series = spc.datasets.emission_series()
print(series)
print([f"{s.metadata['excitation_nm']:.0f}" for s in series])
```

## One file, many spectra

Spectrofluorimeters export a series as one wide table — here 42 columns of
paired wavelength/intensity, with the series names on a sparse first header row
and the column labels on a second. The generic table reader handles that shape:

```python
series = spc.io.read_spectra("J_Peri.csv", 'table', paired=True)
```

You get a `SpectrumCollection` back, not a single spectrum, because the file
genuinely holds many.

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7, 3.4))
viz.plot_collection(series[::4], ax, frame=True)
```

Emission peaks near 520 nm regardless of excitation, which is what you expect
if one emitting species is responsible for all of it.

## Building the map

```{code-cell} ipython3
excitation = np.array([s.metadata['excitation_nm'] for s in series])
emission, intensity = series.to_matrix()

print(f"{intensity.shape[0]} excitations x {intensity.shape[1]} emission points")
```

`to_matrix()` is the same call the multivariate analysis uses — a collection of
spectra on a common axis is a matrix, whatever you then do with it.

## Masking the scatter

An EEM is dominated by light that was never absorbed at all. Two features have
to be removed before the map means anything:

- **Rayleigh scatter**, where emission wavelength equals excitation wavelength.
- **Second-order scatter**, at twice the excitation wavelength, let through by
  the grating.

Neither is fluorescence. Both are far more intense than the signal, and left in
they set the colour scale so that everything real looks flat.

```{code-cell} ipython3
def mask_scatter(excitation, emission, intensity,
                 rayleigh=15.0, second_order=30.0):
    """Blank the scatter ridges. Returns a copy with NaN where they were."""
    masked = np.array(intensity, dtype=float)
    for row, ex in enumerate(excitation):
        masked[row, np.abs(emission - ex) < rayleigh] = np.nan
        masked[row, np.abs(2 * ex - emission) < second_order] = np.nan
    return masked

clean = mask_scatter(excitation, emission, intensity)
print(f"{100 * np.isnan(clean).mean():.0f}% of the map blanked as scatter")
```

:::{admonition} Why NaN and not zero
:class: note

Blanking to zero would put a hard black stripe through the map and drag any
average taken over it. NaN is excluded from contouring and from statistics,
which is what "there is no measurement here" should mean.
:::

## The map

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(6.5, 4))

mesh = ax.contourf(emission, excitation, clean, levels=30, cmap='viridis')
fig.colorbar(mesh, ax=ax, label='Fluorescence (a.u.)')

ax.set_xlabel('Emission wavelength (nm)')
ax.set_ylabel('Excitation wavelength (nm)')
viz.set_frame(ax, True)
```

The emission stays at one wavelength while the excitation varies, and the
intensity picks out which excitation bands feed it — the signature of a single
emitting species with more than one absorbing transition.

:::{admonition} Colour on a map is a different job
:class: tip

The categorical palette used elsewhere in this library is for telling series
apart. A map encodes *magnitude*, which needs a sequential ramp — one hue,
light to dark, like `viridis` here. Never a rainbow: it invents boundaries
where the data is smooth and it is unreadable in greyscale or with colour-vision
deficiency.
:::

## A slice through it

Reading one excitation back out is just indexing the collection, and the result
is an ordinary spectrum with everything else in the library available to it:

```{code-cell} ipython3
best = int(np.nanargmax(np.nanmax(clean, axis=1)))
slice_ = series[best]
print(f"strongest at excitation {slice_.metadata['excitation_nm']:.0f} nm")

peaks = slice_.find_peaks(prominence=0.1, relative=True, distance=10)

fig, ax = plt.subplots(figsize=(7, 3.2))
slice_.plot(ax, frame=True)
peaks.strongest(2).annotate(ax, offset=20)
```

## Still to come

The corrections that make an EEM quantitative are not written yet:

- **Inner-filter correction.** At absorbances above about 0.1 the excitation
  light is attenuated on the way in and the emission on the way out, so
  intensities are understated — increasingly so at the absorption maximum,
  which distorts exactly the comparison an EEM is for.
- **Transfer efficiency** from absorptance divided by an excitation spectrum.
  The `absorptance` unit exists for it; the analysis does not yet.

Until then, keep samples dilute (peak absorbance below ~0.1) and treat the map
as qualitative.

## What to change for your own data

| line | change to |
|---|---|
| `spc.datasets.emission_series()` | `spc.io.read_spectra("yours.csv", 'table', paired=True)` |
| `rayleigh=15, second_order=30` | wider if your slits are wide, narrower if not |
| `levels=30` | fewer for a cleaner map, more for fine structure |

If your instrument writes one shared wavelength column instead of paired
columns, use `read_spectra(..., 'table', x_col=0)`.
