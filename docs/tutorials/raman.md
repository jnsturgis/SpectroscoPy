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

# Raman: getting a spectrum out from under its background

Raman spectra arrive sitting on a background that is usually much larger than
the signal — fluorescence from the sample or its impurities, rising smoothly
and swamping the sharp Raman bands. Removing it without eating the bands is
most of the work.

This uses the tannic acid spectrum that ships with the package.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz

raman = spc.datasets.load('tannic_acid')
print(raman)

fig, ax = plt.subplots(figsize=(7, 3))
raman.plot(ax)
```

Sharp bands, and underneath them a background that climbs steeply towards the
low-wavenumber end — towards the laser line, where fluorescence and the wing of
the Rayleigh scattering are strongest. It is several times the size of the
signal. The x axis is drawn high-to-low because that is the Raman convention;
you did not have to ask for it.

## Two ways to take the background off

The library has three baseline methods, and the choice matters more here than
anywhere else.

**Rubberband** stretches a taut band under the spectrum — the lower convex
hull. It has no parameters, which is its great virtue, but it can only ever be
*convex*, so it cannot follow a background that curves the other way.

**Asymmetric least squares (ALS)** fits a smooth curve that is pulled below the
peaks by weighting points above it very lightly. Two parameters: `lam` for
stiffness and `p` for how hard it is pulled under.

```{code-cell} ipython3
hull = raman.baseline('rubberband')
als = raman.baseline('als', lam=1e3, p=0.005)

fig, ax = plt.subplots(figsize=(7.5, 3.6))
raman.plot(ax, label='measured', frame=True)
ax.plot(hull.x, hull.y, label='rubberband', lw=1.4, ls='--',
        color=viz.PALETTE[1])
ax.plot(als.x, als.y, label='ALS', lw=1.4, color=viz.PALETTE[3])
ax.legend(frameon=False)
```

:::{admonition} Choosing `lam` and `p`
:class: tip

`lam` sets stiffness: too small and the baseline climbs into the peaks, too
large and it becomes a gentle arc that cannot follow the background at all. `p`
sets asymmetry — 0.001 to 0.02 is the usual range; make it smaller if the
baseline is being dragged up by the bands, larger if it dips below them.

**Look at the picture.** Writing this page I first tried `lam=1e8`, checked
numerically that the baseline stayed under the data — it did, everywhere — and
only saw on plotting it that it was a flat arc sitting at 45 while the real
background rose past 270. Every numeric check it was given, it passed. The
figure took one second to disprove it.
:::

## The corrected spectrum

```{code-cell} ipython3
corrected = raman.baseline_correct('als', lam=1e3, p=0.005)

fig, ax = plt.subplots(figsize=(7, 3))
corrected.plot(ax)
```

Two checks worth making every time. The corrected baseline should sit at zero
between bands, and it should not go significantly negative — a negative trough
means the baseline was pulled through a real band.

```{code-cell} ipython3
quiet = corrected.crop(1900, 2600)          # a region with no bands
print(f"between bands : {quiet.y.mean():+.1f} +/- {quiet.y.std():.1f}")
print(f"most negative : {corrected.y.min():+.1f}")
print(f"tallest band  : {corrected.y.max():+.1f}")
```

The mean between bands should be near zero and the most negative value should
be comparable to the noise, not to the band heights. If it is not, the baseline
has been pulled through something real.

## Peaks

```{code-cell} ipython3
normalised = corrected.normalize('max')
peaks = normalised.find_peaks(prominence=0.04, relative=True, distance=12)

print(peaks.strongest(10).sorted_by_position())
```

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7.5, 3.8))
normalised.crop(200, 1800).plot(ax, frame=True)
ax.set_ylim(-0.05, 1.45)

peaks.within(200, 1800).strongest(7).annotate(ax, offset=0.03)
viz.annotate_bands(ax, {
    1610: "aromatic C=C",
    1320: "C-O / phenol",
    750:  "ring breathing",
}, y=1.18)
```

## Still to come

Two Raman-specific corrections are planned but not written:

- **Cosmic-ray removal.** A cosmic ray hitting the detector leaves a spike one
  or two pixels wide, far narrower than any real band, and it will wreck a peak
  fit. The usual detection compares repeated acquisitions of the same spot, or
  looks for features narrower than the instrument response.
- **Fluorescence background modelling** beyond the general-purpose baselines
  used here.

Until then, if you see a suspiciously narrow single-pixel spike, crop it out
and say so in your methods.

## What to change for your own data

| line | change to |
|---|---|
| `spc.datasets.load('tannic_acid')` | `spc.Spectrum("my_spectrum.txt")` |
| `lam=1e3, p=0.005` | tune against your own background, looking at the picture |
| `crop(1900, 2600)` | a genuinely empty region of *your* spectrum |
| the band dictionary | your own assignments |
