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

# The data model

Two objects carry everything: `Spectrum` and `SpectrumCollection`.

## Spectrum

A spectrum is two arrays plus what is needed to interpret them.

| attribute | what it is |
|---|---|
| `x`, `y` | numpy arrays, same length |
| `x_unit`, `y_unit` | machine-readable: `'cm^-1'`, `'nm'`, `'absorbance'`, `'transmittance'` … |
| `x_quantity`, `y_quantity` | what is being measured: `'Wavenumber'`, `'Absorbance'` … |
| `x_label`, `y_label` | display strings, composed from the two above unless overridden |
| `technique` | `'ATR-FTIR'`, `'UV-Vis'`, `'Raman'`, `'Fluorescence'` … |
| `metadata` | free-form dict: sample, reference, acquisition details |
| `history` | what has been done to it (see below) |
| `fileinfo`, `path` | where it came from |

The units are stored separately from the label because they have to be
*usable*, not just printable. That is what makes `spectrum.to('nm')` possible,
and what stops a wavelength axis being silently plotted as if it were
wavenumbers.

```python
spectrum = spc.Spectrum("sample.dpt")
spectrum.set_type("ATR-FTIR")      # sets quantities, units and conventions
spectrum.set_sample("PG_coli")
```

:::{admonition} set_type does not overwrite what a file declared
:class: note

If a file states its units — JCAMP does, `.spy` does — those are kept.
`set_type('FTIR')` will not relabel a transmittance spectrum as absorbance just
because absorbance is the usual FTIR ordinate. Pass `force_units=True` when a
file's own metadata is wrong.
:::

## Building one from data

Not every spectrum comes from a file. A simulated band, a model, a difference
worked out elsewhere — pass the arrays directly:

```{code-cell} python
import numpy as np
import spectroscopy as spc

x = np.linspace(1500, 1750, 251)
amide_i = 0.8 * np.exp(-((x - 1652) / 12) ** 2)      # α-helix
amide_ii = 0.5 * np.exp(-((x - 1548) / 15) ** 2)

model = spc.Spectrum(x, amide_i + amide_ii,
                     technique="ATR-FTIR",
                     name="two-band model",
                     metadata={"sample": "simulated"})
model
```

It behaves like any other spectrum from that point on — the same processing, the
same history, the same file formats:

```{code-cell} python
peaks = model.find_peaks(prominence=0.1, relative=True)
print(peaks)
```

`x` and `y` must be the same length and one-dimensional; anything
`numpy.asarray` accepts will do, so lists are fine.

**Units come from the technique unless you say otherwise.** Passing `x_unit`,
`y_unit`, `x_quantity` or `y_quantity` overrides the technique's conventions and
marks them authoritative, so a later `set_type()` will not relabel them — the
same protection a file's own declared units get:

```{code-cell} python
excitation = spc.Spectrum(np.linspace(400, 700, 301), np.zeros(301),
                          technique="Fluorescence",
                          x_quantity="Wavelength", x_unit="nm",
                          y_quantity="Intensity", y_unit="counts")
print(excitation.x_label, "/", excitation.y_label)
```

Building a spectrum records **no** history entry. Creating data is not
processing it, and a history that claimed otherwise would be lying about where
the numbers came from.

## Reading one from a file

```{code-cell} python
from spectroscopy import datasets

spectrum = spc.Spectrum.read(datasets.path("ethanol"))
print(spectrum.technique, len(spectrum), "points")
```

`Spectrum.read(path)` and `Spectrum(path)` do the same thing; `read` is
preferred in new code because it says that it touches the disk, and it takes the
format as a named argument when the extension is unhelpful:

```python
spectrum = spc.Spectrum.read("scan.txt", "dpt")
```

:::{admonition} Two positional strings mean a directory and a filename
:class: warning

`Spectrum("path/to/", "oddly_named_file")` reads a *file called
`oddly_named_file` in directory `path/to/`*. It does **not** mean
`(filename, format)` — that is `Spectrum.read(filename, format)`. The
constructor's docstring claimed the second reading until 0.1.1, and the call it
described raised `TypeError: Unknown filetype unknown`.
:::

## Everything returns a new spectrum

No operation modifies the spectrum it is called on.

```python
cropped = spectrum.crop(900, 1800)     # spectrum is unchanged
```

This costs a copy and buys the ability to branch, compare and go back — and it
is what makes the history meaningful, since each result carries its own.

## History

Every operation records its name and its arguments.

```python
result = (spectrum.crop(900, 1800)
                  .baseline_correct('rubberband')
                  .normalize('max', window=(1050, 1080)))

print(result.describe_history())
```

```
1. crop(x_max=1800, x_min=900)
2. baseline_correct(method='rubberband', mode='subtract')
3. normalize(method='max', window=[1050, 1080])
```

Steps store structured `(name, params)` pairs rather than sentences, so they
are machine-readable, JSON-serialisable, and survive a round trip through
`.spy`. Arithmetic is recorded too, including the scale factor on a reference
subtraction — which is usually the one number in an analysis chosen by
judgement rather than by rule.

The one gap worth knowing: a spectrum used as an *operand* is recorded by name
and sample, not by value. Replaying a reference subtraction needs the reference
supplied from outside.

## Arithmetic

`+ - * /` work between spectra and with numbers.

```python
difference = sample - 0.995 * buffer
average    = (a + b + c) / 3
```

Spectra of different lengths raise rather than broadcasting; same length but
different x positions warns. Use `resample` to put one onto another's axis:

```python
buffer = buffer.resample(sample.x)
```

## SpectrumCollection

An ordered set of spectra, which is what a folder of files actually is.

```python
spectra = spc.SpectrumCollection.from_files("data/*.dpt", technique="ATR-FTIR")

averages = {name: group.mean()
            for name, group in spectra.group_by('sample').items()}
```

The sample name comes from the file name up to the first dot, so
`PG_coli.0.dpt` … `PG_coli.2.dpt` group together. Pass `sample_from=` for a
different convention.

Reductions: `mean`, `median`, `std`, `sem`. Batch operations —
`crop`, `baseline_correct`, `normalize`, `smooth`, `resample`,
`subtract_reference` — apply to every member and return a new collection.
`select(predicate)` filters, `map(function)` does anything else.

`to_matrix()` returns `(x, X)` with `X` of shape (n_spectra, n_points), which
is the handover to any multivariate analysis.

## Units

```python
spectrum.to('cm^-1')              # x axis
spectrum.to(y_unit='%T')          # y axis
spectrum.to('nm', 'absorbance')   # both
```

x: `nm`, `um`, `cm^-1`, `eV`. y: `absorbance`, `transmittance`, `%T`,
`absorptance`, `%absorptance`.

A reciprocal conversion reverses the axis, so the points are re-sorted and y is
reordered with them — otherwise every later interpolation and crop would be
working on a descending axis.

:::{admonition} absorbance is not absorptance
:class: warning

Absorbance (optical density) is `-log10(T)`. Absorptance is the fraction of
light actually absorbed, `1 - T`. They agree only in the weak-absorption limit;
by A = 0.3 they differ by nearly 30%. Absorptance is the one to compare with an
excitation spectrum, because the excitation signal follows photons absorbed,
not optical density.
:::

`counts` and `a.u.` have no defined conversion and saying so is better than
inventing a factor.
