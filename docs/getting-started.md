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

# Getting started

This page takes you from a fresh install to a finished, labelled figure. It
assumes you can open a terminal and copy text into it. It does **not** assume
you write Python.

Every block below is real, runnable code, and every picture underneath a block
was produced by running it — so you can copy any of them, change the file name,
and get the same thing with your own data. Nothing here needs data of your own:
the package ships a handful of example spectra.

:::{admonition} Which kind of user are you?
:class: tip

**"I just want a figure from my spectra."** Read this page and stop. The five
recipes here cover most of what the library is for.

**"I want to script a whole analysis."** Read this page, then the User guide.

**"I don't want to write any code at all."** Then this page is not really for
you yet, and that is a fair thing to want. A graphical application is planned;
until it exists, the closest thing is to copy a recipe from here and change the
file name. Ask for help — it is meant to be one line.
:::

## Install

```
pip install spectroscopy
```

If you also want the multivariate analysis (PCA, NMF), ask for it explicitly —
it pulls in a large extra package, so it is not installed by default:

```
pip install "spectroscopy[multivariate]"
```

:::{admonition} If pip refuses with "externally-managed-environment"
:class: note

Recent Linux distributions protect the system Python. The safe answer is a
virtual environment:

```
python3 -m venv ~/spectro-env
~/spectro-env/bin/pip install spectroscopy
```

then use `~/spectro-env/bin/python` (or select that environment in Jupyter).
:::

## Load a spectrum and look at it

The library ships some example spectra so that everything on this page runs
before you have loaded any of your own.

```{code-cell} ipython3
import spectroscopy as spc

print(spc.datasets.describe())
```

Loading one and plotting it takes two lines. The axis labels, the units, and
the right-to-left wavenumber axis that infrared spectra are conventionally
drawn with all come from the file — you do not set them.

```{code-cell} ipython3
import matplotlib.pyplot as plt

spectrum = spc.datasets.load('ethanol')

fig, ax = plt.subplots(figsize=(7, 3))
spectrum.plot(ax)
```

To use **your own** file instead, replace that one line. The format is worked
out from the file itself, so `.dpt`, `.jdx`, `.dx`, `.csv`, `.txt` and `.spy`
all just work:

```python
spectrum = spc.Spectrum("my_sample.dpt")
```

```{code-cell} ipython3
print(spectrum)
print(f"{len(spectrum)} points from {spectrum.x.min():.0f} to "
      f"{spectrum.x.max():.0f} {spectrum.x_unit}")
```

## Convert the units if you need to

This particular file is a **transmittance** spectrum — the library knows,
because the file says so, and the y axis above is labelled accordingly. Peaks
point downwards. Most analysis wants absorbance, where they point up:

```{code-cell} ipython3
absorbance = spectrum.to(y_unit='absorbance')

fig, ax = plt.subplots(figsize=(7, 3))
absorbance.plot(ax)
```

The same call converts the x axis — `spectrum.to('nm')` on an infrared
spectrum, for instance. This is why the library insists on knowing the units
rather than treating the axis as anonymous numbers: it is what stops a
wavelength spectrum being quietly plotted as if it were wavenumbers.

## Clean it up

The usual sequence is: cut out the region you care about, take the baseline
off, and scale it so spectra can be compared. Each step returns a **new**
spectrum, so the original is never damaged and you can always go back.

```{code-cell} ipython3
clean = (absorbance
         .crop(950, 1500)
         .baseline_correct('rubberband')
         .normalize('max'))

fig, ax = plt.subplots(figsize=(7, 3))
clean.plot(ax)
```

It is worth looking at what the baseline actually did before trusting it — a
baseline that has cut through a real band is obvious here and invisible in the
corrected spectrum on its own.

Which envelope counts as "the baseline" depends on which way the bands point.
In absorbance it is the *lower* hull, tending to zero; in transmittance it is
the *upper* hull, tending to 100%, and the correction is a division rather than
a subtraction. The library reads that off the y unit, so you get the right one
without asking — but it is worth knowing when you look at the picture.

```{code-cell} ipython3
from spectroscopy import viz

region = absorbance.crop(950, 1500)

fig, ax = plt.subplots(figsize=(7, 3))
viz.plot_baseline(region, region.baseline('rubberband'), ax)
```

## Find and label the peaks

`find_peaks` works on the second derivative by default, which finds shoulders
on a broad band as well as obvious maxima.

```{code-cell} ipython3
peaks = clean.find_peaks(prominence=0.05, relative=True, distance=20)

print(peaks.sorted_by_position())
```

:::{admonition} `relative=True` is worth keeping
:class: note

Without it, `prominence` is measured on the second derivative, which is
typically ten thousand times smaller than the spectrum — so a sensible-looking
`prominence=0.05` finds nothing at all and you have to guess your way down to
`0.00001`. With `relative=True` the number means "5% of the signal range" and
behaves the way you expect.
:::

Putting the labels on the figure is one call:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7.5, 3.4))
clean.plot(ax)
ax.set_ylim(-0.05, 1.6)

peaks.strongest(5).annotate(ax, offset=0.03)
viz.annotate_bands(ax, {1050: "C-O stretch"}, y=1.35)
```

:::{admonition} Label only what you can read
:class: tip

Two settings are doing real work. `distance=20` refuses peaks closer than 20
points together, and `strongest(5)` keeps only the five most prominent.
Labelling every peak on a crowded spectrum produces an overlapping smear that
is worse than no labels at all — pick the few you want to talk about.
:::

### Framed axes

By default the top and right edges are left off. Many physical-chemistry
journals want a full box, with the ticks turned inwards — pass `frame=True`:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7, 3))
clean.plot(ax, frame=True)
```

If you always want that, set it once and forget it:

```python
from spectroscopy import viz
viz.FRAME_AXES = True
```

## Several spectra at once

Most real work is a folder of files: replicates of a few samples. A
`SpectrumCollection` loads them together, groups them by sample and averages
each group — without you counting rows or indexing anything.

```python
spectra = spc.SpectrumCollection.from_files("data/*.dpt", technique="ATR-FTIR")
averages = {name: group.mean()
            for name, group in spectra.group_by('sample').items()}
```

The sample name is taken from the file name up to the first dot, so
`PG_coli.0.dpt`, `PG_coli.1.dpt` and `PG_coli.2.dpt` are recognised as three
replicates of `PG_coli`.

Here is the same idea with the two example spectra:

```{code-cell} ipython3
pair = spc.datasets.load_pair()
print(pair)

fig, ax = plt.subplots(figsize=(7, 3))
pair.plot(ax)
```

Stacking reads better than overlaying once there are more than a handful, and
labels each trace directly so you never have to match a colour to a key:

```{code-cell} ipython3
fig, ax = plt.subplots(figsize=(7, 3))
viz.stack(pair.crop(400, 900), ax)
```

## Subtract a reference

Subtracting a buffer or water spectrum, scaled by eye until the unwanted bands
cancel, is the most common single operation in this kind of work. The scale
factor is an argument, so it is **recorded** rather than left as a bare number
in a notebook cell:

```{code-cell} ipython3
sample, reference = pair[0], pair[1]

difference = sample.subtract_reference(reference, factor=0.98)

fig, ax = plt.subplots(figsize=(7, 3))
difference.crop(400, 900).plot(ax)
```

## What did I actually do?

Every operation records itself. This is the answer to "how did I get from the
raw file to this plot?" — six months later, or for a reviewer.

```{code-cell} ipython3
print(clean.describe_history())
```

Saving in the library's own `.spy` format keeps that history in the file, along
with the units and the sample name, so it survives being closed and reopened:

```{code-cell} ipython3
import tempfile, os

target = os.path.join(tempfile.mkdtemp(), "cleaned.spy")
clean.save_as(target, "spy")

reloaded = spc.Spectrum(target)
print(reloaded.describe_history())
```

Any other format works the same way — `clean.save_as("cleaned.csv", "csv")` —
but only `.spy` keeps the history.

## Where to go next

- **User guide** — one page per idea: the data model, reading files,
  processing, plotting.
- **Tutorials** — a complete worked analysis per technique, taken from real
  experiments, downloadable as notebooks you can run and edit.
- **API reference** — every function and what it takes. You should not need
  this to get work done; it is there when you do.

:::{admonition} Something did not work?
:class: warning

The library tries to say what to do rather than just what went wrong. If a
message is unhelpful, that is a bug worth reporting — unclear errors are the
main thing that stops people using a tool like this.
:::
