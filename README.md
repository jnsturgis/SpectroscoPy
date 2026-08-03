# SpectroscoPy

One way of handling spectra, whatever instrument they came from.

[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](LICENSE)

> ⚠️ **Pre-1.0.** The API is still settling and may change between 0.x
> releases. It is used for real work, but pin a version if you depend on it.

## Why

Every spectrometer writes its own file format and ships its own idiosyncratic
analysis program. The result is that the same operation — average the
replicates, take the baseline off, subtract the buffer, label the peaks — gets
re-implemented, slightly differently, in every notebook.

SpectroscoPy reads them all into one object, gives you those operations once,
and keeps a record of everything it did, so a figure can be traced back to the
raw files months later.

```python
import spectroscopy as spc

spectra = spc.SpectrumCollection.from_files("data/*.dpt", technique="ATR-FTIR")

result = (spectra.group_by('sample')['PG_coli']
          .mean()                                   # average the replicates
          .crop(900, 1800)
          .subtract_reference(water, factor=0.7)    # the factor is recorded
          .baseline_correct('rubberband')
          .normalize('max', window=(1050, 1080)))

peaks = result.find_peaks(prominence=0.05, relative=True)
print(result.describe_history())
```

```
1. mean(n_spectra=6)
2. crop(x_max=1800, x_min=900)
3. subtract_reference(factor=0.7, reference={'kind': 'spectrum', 'sample': 'H2O'})
4. baseline_correct(method='rubberband', mode='subtract')
5. normalize(method='max', window=[1050, 1080])
```

## Install

```bash
pip install spectroscopy                     # core
pip install "spectroscopy[multivariate]"     # + PCA / NMF / ICA
```

Python 3.10 or newer. If `pip` refuses with *externally-managed-environment*,
use a virtual environment — see Getting started.

## What it does

**Reads** `.dpt` (Bruker OPUS), JCAMP-DX (`.jdx`, `.dx`), delimited text
(`.csv`, `.tsv`, `.txt`), wide and paired multi-column exports, and its own
`.spy`. Separators, header rows and text encodings are detected rather than
assumed. Adding a format is one decorator.

**Processes** — cropping, baselines (rubberband, guide-point polynomial,
asymmetric least squares), smoothing and derivatives, five normalisations,
second-derivative peak detection, replicate averaging, scaled reference
subtraction, resampling, unit conversion.

**Analyses** — PCA, NMF and ICA over a set of spectra, with bootstrap stability
testing to say whether the number of components is defensible.

**Plots** — with the axis conventions each technique expects (including
right-to-left wavenumbers and optional framed axes), a colour-vision-safe
palette, peak and band annotation, stacked and faceted layouts.

**Remembers** — every operation records what it did and with which parameters,
and `.spy` files carry that history with them.

## Documentation

**<https://jnsturgis.github.io/SpectroscoPy/>** — published from `main` on every
push, so it always matches the code above it.

Start with [Getting Started](https://jnsturgis.github.io/SpectroscoPy/getting-started.html),
or, if you know what you want to find out rather than which function does it,
[What do you want to know?](https://jnsturgis.github.io/SpectroscoPy/what-do-you-want-to-know.html)

To build it yourself:

```bash
pip install -e ".[docs]"
cd docs && make html          # or `make rebuild` for a clean build
```

Pages execute their own code at build time, so an example that has gone stale
breaks the build rather than quietly misleading a reader. Use `make rebuild`
before trusting the page-to-page navigation: Sphinx caches the table of
contents, and an incremental build can leave a newly-added page unlinked.

Getting Started goes from `pip install` to a
labelled figure and assumes only that you can copy text into a terminal. Every
example runs against sample spectra that ship with the package, so it works
before you have any data of your own.

## Project layout

```
spectroscopy/
  spectra.py      Spectrum -- the core data model
  collection.py   SpectrumCollection -- many spectra, grouped and averaged
  io/             format readers, self-registering
  processing/     algorithms over plain arrays
  viz.py          plotting
  history.py      what was done, with which parameters
  units.py        axis unit conversion
  datasets.py     the example spectra
```

The layering runs one way: `io` → `core` → `processing` → `viz`. Nothing
reaches backwards.

## Status

Working: the data model, provenance, the format registry, the processing and
multivariate layers, plotting, and a getting-started guide.

Planned: more vendor formats, technique-specific corrections (ATR, inner-filter
and scatter removal for EEM, Beer-Lambert), hyperspectral maps, and a graphical
application for people who would rather not write code at all.

`SpectroscoPy_Development_Roadmap.md` has the plan;
`SpectroscoPy_Codebase_Review.md` records what has been done and why, including
several things found in the old code that had been quietly affecting results.

## Contributing

Bug reports are welcome, particularly *unclear error messages* — those stop
people using a tool more effectively than missing features do. **You do not
need a GitHub account**: every documentation page has a feedback link that
opens an email.

A file that fails to read is one of the most useful things you can send. Both
native binary readers exist because somebody supplied files.

```bash
pip install -e ".[dev]"
pytest
ruff check .
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the longer version, including how
to add a file format, and [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).

## Citing

[CITATION.cff](CITATION.cff) — or use GitHub's "Cite this repository" button.
Cite the repository URL rather than a `pip` name for now: the PyPI name is not
yet settled.

## Licence

[Mozilla Public License 2.0](LICENSE). Use it in anything; changes to these
files should be shared back.
