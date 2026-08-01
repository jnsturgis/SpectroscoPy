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

# ATR-FTIR: a complete analysis

This is a real analysis, start to finish: nine spectra recorded on a diamond
ATR — three replicates each of glucose, cellulose and water — worked up into a
figure comparing the two sugars with their bands labelled.

It is the workflow this library was written for, and it is the one that used to
take forty lines of hand-indexed arithmetic in a notebook.

The data ships with the package, so you can run every cell.

:::{admonition} What you will end up with
:class: tip

A labelled comparison of two polysaccharides, with the water contribution
removed, baselines corrected, and a record of every step that produced it.
:::

```{code-cell} ipython3
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz
```

## 1. Load the folder

Nine files, one call. Sample names come from the file name up to the first dot,
so `Glucose.1.dpt`, `Glucose.2.dpt` and `Glucose.3.dpt` are recognised as three
replicates of the same thing.

```{code-cell} ipython3
spectra = spc.datasets.ftir_replicates()

print(spectra)
for name, group in spectra.group_by('sample').items():
    print(f"  {name:<12} {len(group)} replicates")
```

With your own data that first line is:

```python
spectra = spc.SpectrumCollection.from_files("2025/10/13/*.dpt",
                                            technique="ATR-FTIR")
```

## 2. Look before you leap

Always. A replicate that went wrong — a bubble, a dry patch, a smear of the
previous sample — is obvious in a panel-per-sample overview and very hard to
spot after averaging.

```{code-cell} ipython3
fig, axes = viz.grid(spectra, key='sample', ncols=1, figsize=(7, 6))
```

The three replicates of each sample lie on top of one another, so nothing needs
throwing away. If one did, `spectra.select(...)` drops it.

## 3. Average the replicates

```{code-cell} ipython3
averages = {name: group.mean()
            for name, group in spectra.group_by('sample').items()}

water = averages['H2O']
print(water.describe_history())
```

That is the whole of it. The old way — `(spectra[32] + spectra[33] + ... +
spectra[37]) / 6.0` — silently gave the wrong answer the moment a file was
added to the folder.

## 4. Subtract the water

Sugars are measured wet, so every spectrum carries a water contribution that
has to come off. How much is a judgement made by eye: too little and the broad
O–H bending band near 1640 cm⁻¹ survives, too much and it turns into a trough.

The factor is an argument, so it ends up recorded rather than living as a bare
number in a cell you will not recognise later.

```{code-cell} ipython3
region = (900, 1800)

glucose = (averages['Glucose']
           .crop(*region)
           .subtract_reference(water.crop(*region), factor=0.70))

cellulose = (averages['CelluloseX']
             .crop(*region)
             .subtract_reference(water.crop(*region), factor=0.002))

fig, ax = plt.subplots(figsize=(7, 3))
viz.plot_collection(spc.SpectrumCollection([glucose, cellulose],), ax,
                    labels=['glucose', 'cellulose'])
```

:::{admonition} Choosing the factor
:class: note

Try a few and look. Increase it until the 1640 cm⁻¹ water band flattens, and
stop before it goes negative — a trough there means you have over-subtracted.
Cellulose here is nearly dry, hence the very small factor; glucose was measured
in solution.
:::

## 5. Baseline and normalise

ATR spectra sit on a sloping background from the changing penetration depth.
The rubberband baseline handles it without any parameters to choose.

Normalising to a band the samples have in common — here the C–O stretch region
around 1050–1080 cm⁻¹, which every polysaccharide has — makes the shapes
comparable regardless of how much material was on the crystal.

```{code-cell} ipython3
def workup(spectrum):
    return (spectrum
            .baseline_correct('rubberband')
            .normalize('max', window=(1050, 1080)))

glucose, cellulose = workup(glucose), workup(cellulose)

fig, ax = plt.subplots(figsize=(7, 3))
viz.plot_collection(spc.SpectrumCollection([glucose, cellulose]), ax,
                    labels=['glucose', 'cellulose'])
```

## 6. Find the peaks

Peak detection works on the second derivative, so shoulders on the broad
glycan envelope are found as well as the obvious maxima — which is the whole
reason for doing it that way on this kind of spectrum.

```{code-cell} ipython3
peaks = glucose.find_peaks(prominence=0.06, relative=True, distance=15)
print(peaks.within(950, 1200).sorted_by_position())
```

## 7. The figure

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7.5, 4))

viz.stack(spc.SpectrumCollection([cellulose, glucose]), ax,
          labels=['cellulose', 'glucose'], gap=1.15, frame=True)

viz.annotate_bands(ax, {
    1640: "H-O-H bend",
    1370: "C-H bend",
    1150: "C-O-C asym",
    1030: "C-O stretch",
}, y=2.35)
ax.set_ylim(-0.1, 2.9)
```

## 8. Keep the provenance

```{code-cell} ipython3
print(glucose.describe_history())
```

Saved as `.spy`, that history travels with the file:

```{code-cell} ipython3
import os, tempfile

target = os.path.join(tempfile.mkdtemp(), "glucose_worked_up.spy")
glucose.save_as(target, "spy")

print(spc.Spectrum(target).describe_history())
```

Six months from now, that answers "what did I do to get this?" without you
having to remember.

## What to change for your own data

| line | change to |
|---|---|
| `spc.datasets.ftir_replicates()` | `spc.SpectrumCollection.from_files("your/folder/*.dpt", technique="ATR-FTIR")` |
| `factor=0.70` | whatever flattens your water band |
| `window=(1050, 1080)` | a band your samples have in common |
| `region = (900, 1800)` | the range you care about |

Everything else stays as it is.
