# ADR-0001 — The core data model

**Status:** Accepted, with four questions still open (§7) that must be closed
before the 1.0 API freeze.
**Date:** 2026-08-02
**Context:** Roadmap §2.4 asks for this before more code is written against the
model. It arrives late — after Phases 0–5 — so it records decisions that were
made under pressure from real code, rather than decisions taken in advance. That
makes it more useful as a record and less useful as a guide, and §7 is where the
difference shows.

This is the document to read before proposing a breaking change to `spectroscopy.spectra`
or `spectroscopy.collection`. If a proposal contradicts something here, the ADR
is what has to be argued with.

---

## 1. What the model has to do

Four techniques (ATR-FTIR, FTIR, UV-Vis, Raman, fluorescence), a dozen file
formats, and one requirement that shapes everything else: **a figure must be
traceable back to the raw file six months later**, by the person who made it or
by a reviewer who did not.

That requirement is why provenance is in the core object rather than in a
wrapper around it, and why serialisation had to carry it.

---

## 2. Decisions

### 2.1 A purpose-built class over numpy arrays, not xarray

**Decided: `Spectrum` holds `x` and `y` as plain `numpy` arrays.**

Roadmap §2.1 suggested backing the object with `xarray.DataArray` to get
interpolation, arithmetic, slicing and NetCDF serialisation for free. That is
not what was built, and the reason is worth stating because the roadmap's
argument was a good one.

- The free operations turned out not to be the expensive part. Interpolation is
  one `CubicSpline` call; arithmetic needed custom axis-compatibility checking
  anyway (defect D3); slicing on a spectral axis means `crop`, which is four
  lines of `numpy` indexing.
- The expensive part was **units and provenance**, and xarray gives neither.
  `.attrs` is an untyped dict that xarray does not interpret, so `x_unit` and
  `history` would have been our code either way.
- xarray brings pandas transitively, which conflicts with §2.6.

**Cost, recorded honestly:** we own the serialisation (`.spy`, §2.4), and we do
not get NetCDF/Zarr, dask chunking, or n-dimensional indexing. The last of those
is the one that matters — see §6.3 on hyperspectral maps.

**Revisit when:** spectral *maps* arrive (a spatial grid of spectra). At that
point n-dimensional indexing stops being a luxury and the arithmetic of
re-implementing it is worse than the dependency.

### 2.2 Units are a native table, not `pint`

**Decided: `spectroscopy.units` implements conversion directly.**

Roadmap §2.2 recommended `pint` from day one. Two things emerged once the
conversions were actually written:

1. **Half of what is needed is not a unit conversion.** Absorbance to
   transmittance is `T = 10**-A` — a functional transform between two
   dimensionless quantities. `pint` does not express that, so custom code was
   required regardless.
2. **The x conversions are reciprocal, not scalar.** nm to cm⁻¹ is `1e7 / nm`,
   which `pint` handles only through its `spectroscopy` context; the call site
   ends up no simpler than the formula.

The unit space here is small and closed, so a table is exact, dependency-free,
and covers the y axis too.

One consequence that is easy to get wrong and is therefore tested: **a
reciprocal conversion reverses the axis**, so `to()` re-sorts ascending and
reorders `y` with it. Without that, every later interpolation, hull and crop
silently operates on a descending axis.

**Revisit when:** the unit space stops being closed — path lengths,
concentrations, per-unit-time intensities. At that point `pint` earns its place.

### 2.3 Provenance is structured, not prose

**Decided: `ProcessingStep(name, params, timestamp)`, where `name` matches the
method that produced it and `params` are its normalised keyword arguments.**

This is the decision the whole "chaining now, `Pipeline` later" plan rests on
(roadmap §2.3). A step records `{"method": "savgol", "polyorder": 2,
"window_length": 15}` rather than the string "smoothed with Savitzky-Golay",
which is what makes a history losslessly convertible into a runnable pipeline —
whether or not that pipeline is ever built (§6.1).

Positional and abbreviated call forms are normalised *before* storage, so
`smooth('SG', [2, 15])` and `smooth(method='savgol', polyorder=2,
window_length=15)` record identically. Steps are frozen; `replace()` returns a
copy.

### 2.4 `.spy` 1.0: provenance that survives leaving memory

**Decided: a JSON header followed by tab-separated data.**

The numbers stay greppable and plottable with ordinary tools — `awk` and
`gnuplot` work on a `.spy` file — while the header carries name, technique,
quantities, units, label overrides, metadata and the full history with
timestamps. A processed spectrum round-trips to the hand-tuned water-subtraction
factor, not merely to the numbers it produced.

The format version is read **from the file**, not assumed. Legacy 0.0 files
still load.

### 2.5 Mutability: processing returns new objects, metadata setters do not

**Decided, and stated precisely because the honest answer is not "immutable":**

| Operation | Behaviour |
|---|---|
| `crop`, `smooth`, `normalize`, `baseline_correct`, `to`, `resample`, `subtract_reference`, arithmetic | Return a **new** `Spectrum` and append a `ProcessingStep`. Built through one private `_derive()` so this cannot drift per method |
| `set_sample`, `set_reference`, `set_type` | **Mutate in place.** They describe what the spectrum *is*, not something done to it, and they record no step |
| `reload`, `save` | Mutate / write in place, by definition |
| `.x`, `.y` | Public, writable attributes |

The split is defensible — labelling a sample is not a processing operation — but
it is a split, not a policy, and anyone extending the class needs to know which
side a new method falls on. The rule: **if it changes the numbers, it returns a
new object and records a step.**

`.x` and `.y` being writable is the loose end; see §7.3.

### 2.6 numpy-native core, pandas at arm's length

**Decided: `numpy` + `scipy` + `matplotlib` are required; `scikit-learn` is an
optional extra; `pandas` is not a dependency at all.**

A DataFrame return type is contagious: every caller then needs pandas to do
anything with the result, and its API surface leaks into docstrings, tests and
tutorials. So:

- Core returns numpy arrays and dataclasses. `PeakTable` is a dataclass over
  arrays. `SpectrumCollection.to_matrix()` returns `(x, X)` as plain arrays —
  which is what scikit-learn wants anyway.
- `to_dataframe()` survives as an explicit escape hatch that imports pandas
  *inside the method*. Never call it, never need pandas.
- **Nothing in `core`, `io` or `processing` may import pandas at module scope.**

### 2.7 `metadata` is about the sample; `fileinfo` is about the file

**Decided.** Provenance about where bytes came from lives in `fileinfo` and is
exposed as `Spectrum.path`. `metadata` describes the specimen and the
measurement.

This looks like housekeeping and is not: the registry briefly wrote
`source_path` into `metadata`, which broke `.spy` round-trip equality, because
reloading added a key the original did not have.

### 2.8 Baselines: one concept, three ways of choosing anchor points

**Decided.** `poly`, `rubberband` and `als` are not three algorithms but one
idea — fit a background through anchor points — differing only in how the
anchors are chosen: given explicitly, found by convex hull, or all points with
asymmetric weights.

`baseline()` returns *the baseline* (for plotting, and for the `s - s.baseline()`
idiom); `baseline_correct()` returns the corrected spectrum and records the
step. `return_points=True` hands back the anchors, so an automatic selection can
be adjusted by eye and fed back explicitly.

### 2.9 Layering, and a lazily imported `viz`

**Decided: `io` → `core` → `processing` → `viz`, no backward reach.**
`Spectrum.plot` delegates to `viz` rather than owning matplotlib.

`viz` is imported lazily through a module `__getattr__` (matplotlib costs
~350 ms, which someone reading a `.dpt` file should not pay).

### 2.10 The public surface is an explicit list

**Decided: ten names.** `Spectrum`, `SpectrumCollection`, `PeakTable`,
`ProcessingStep`, and the submodules `datasets`, `io`, `processing`, `units`,
`lineshapes`, `viz`.

Before this was written down the package re-exported whatever `spectra.py`
happened to import — `spectroscopy.os`, `spectroscopy.np`, forty names, seventeen
of them modules. A freeze would have promised to keep `spc.np` forever. The list
is pinned by `tests/test_public_api.py`, so adding to the public API takes a
deliberate line in a diff.

---

## 3. Rejected alternatives

| Alternative | Why not |
|---|---|
| `xarray.DataArray` backend | §2.1 — the free operations were not the expensive ones, and it brings pandas |
| Subclassing `numpy.ndarray` | Arithmetic on a spectrum needs axis compatibility checks and history; inheriting ndarray fights both |
| `pint` for units | §2.2 — half the conversions are not unit conversions |
| DataFrame return types | §2.6 — contagious |
| Free functions over an object | Provenance has to live somewhere; a free function has nowhere to append a step |
| Free-text history strings | §2.3 — cannot be replayed, and the format is frozen by `.spy` |

---

## 4. Consequences

**Good.** Provenance is automatic rather than remembered. A `.spy` file is
self-describing and still readable with `awk`. Unit confusion — nm vs cm⁻¹,
absorbance vs %T — is caught at the API rather than discovered in a figure. The
core installs with three dependencies.

**Bad, and accepted.** We maintain our own serialisation. We have no
n-dimensional indexing, which is why hyperspectral maps are deferred rather than
awkward. The mutability split (§2.5) has to be explained to every contributor.
Owning unit conversion means owning its edge cases — the axis-reversal
re-sort is one that had to be found by writing a test for it.

---

## 5. Layers this ADR does not cover

`io` (the format registry), `processing.multivariate` (decomposition and
stability) and `viz` (the palette) have their own rationales, recorded in the
codebase review §10, §11 and §12. A future ADR-0002 is expected for joint
decomposition across techniques.

---

## 6. Explicitly deferred — not oversights

Each of these is **additive**: it can arrive in 1.1 without breaking 1.0.

### 6.1 `Pipeline` and `Pipeline.from_history()`

Roadmap §2.3 calls the chaining/pipeline duality settled, and it is — the
expensive half is built. `ProcessingStep` already stores exactly what a pipeline
would need to replay.

**Deferred because** nothing has yet needed to *reuse* a recipe: the notebooks
chain, and the GUI that would batch-apply does not exist. Building it now would
be designing against an imagined caller.

**The constraint it places on 1.0:** the `.history` structure and its `.spy`
serialisation freeze now, and a later `Pipeline` must be reconstructible from
them. That is the reason §2.3 mattered enough to get right early.

### 6.2 `fit_peaks()` and `FitResult`

Roadmap §2.5 sketched both and left the return type open. `find_peaks()` and
`PeakTable` shipped; fitting did not. `lineshapes` provides the components, so
fitting is currently done by hand with `scipy.optimize`.

**Deferred because** the notebooks' fitting is heterogeneous — different
constraint schemes per problem — and a return type designed against one of them
would fit the others badly. When it arrives, §2.6 fixes its shape: a dataclass
over arrays, not a DataFrame.

### 6.3 `SpectralMap` / hyperspectral

Roadmap §2.1 lists it. Nothing in the four target techniques as practised here
produces one. Deferring it is what makes §2.1's numpy decision safe; it is also
the decision most likely to reopen §2.1.

### 6.4 `to_netcdf()`

Follows from §2.1. `.spy` plus `.csv` covers what is needed.

---

## 7. Open — must be settled before the 1.0 freeze

These are not deferrals. They are decisions that have not been made, and
freezing without making them would freeze an accident.

### 7.1 ✅ Resolved 2026-08-02 — a spectrum can now be constructed from arrays

`Spectrum(x, y, technique=..., x_unit=..., name=..., metadata=..., history=...)`
exists, alongside a `Spectrum.read(path, file_type=None)` classmethod. The
string forms still work, so nothing broke.

Three decisions taken while adding it, recorded because they are now API:

- **Explicit units beat the technique's conventions, and mark the units
  authoritative.** Passing `x_unit` or `y_unit` sets the same flag a file's own
  declared units set, so a later `set_type()` will not relabel them. A
  fluorescence excitation scan in nm stays in nm.
- **Building from arrays records no history entry.** Creating data is not
  processing it, and a history claiming otherwise would misdescribe where the
  numbers came from.
- **A 2-D input is refused by name**: the error says a collection of spectra is
  a `SpectrumCollection`, rather than failing later on a shape mismatch.

The original problem, kept because it is why the constructor exists:

### 7.1.1 The state before the fix

`Spectrum.__init__(self, *args)` accepts: nothing, a `Spectrum` (copy), or one
to three strings (file). **There is no way to build a spectrum from `x` and `y`
arrays.** The documented pattern — in the test suite, which is where synthetic
fixtures live — is:

```python
spec = Spectrum()
spec.x, spec.y = x, y          # empty-then-assign
spec.set_type("ATR-FTIR")
```

Roadmap §2.5's primary constructor, `Spectrum(x, y, x_unit=..., y_unit=...)`,
does not exist. Anyone who computes a spectrum — a simulation, a model, a
difference obtained outside the library — has to construct an empty object and
reach into its attributes.

**Recommendation:** add a keyword constructor and a `Spectrum.read(path)`
classmethod before the freeze, keeping the string forms working. Adding is
additive; what is *not* additive is freezing `*args` as the only entry point,
because the empty-then-assign pattern then becomes API that people write against
and that a later constructor cannot take back.

### 7.2 ✅ Resolved 2026-08-02 — the docstring was the wrong half

The behaviour stayed, the docstring was corrected, and `Spectrum.read(path,
file_type)` now provides what the docstring had wrongly promised. Two positional
strings still mean `(directory, filename)`, which is what the docs site and the
readers already assumed; changing that would have broken working calls to fix a
comment. The guide now carries the warning admonition as well, because the wrong
reading is the intuitive one.

The original finding:

### 7.2.1 The contradiction as found

The docstring says `Spectrum(filename, filetype)` with filetype "from
('jcamp','csv')". The code treats two arguments as `(path, name)` and infers the
type from their concatenation. So the documented call fails:

```python
>>> spc.Spectrum("ethanol.jdx", "jcamp")
TypeError: Unknown filetype unknown
```

Verified 2026-08-02. Either the docstring or the behaviour is wrong; the docs
site uses the three-argument form, which suggests the docstring is the stale
one. **This must be resolved before the signature is frozen**, and it is an
argument for §7.1's explicit constructor, which makes the whole overload
question go away.

### 7.3 `.x` and `.y` are writable

**Updated 2026-08-02:** the justification has gone. They had to be writable
while §7.1 left no other way to build a spectrum; now there is one, so this is a
real choice rather than a necessity. The question is live: does a frozen 1.0
promise that assigning `spectrum.y = something` works?

It bypasses history — a spectrum whose y values were replaced by hand records
nothing, which contradicts §1. **Recommendation:** keep them writable for 1.0
(read-only properties would break the fixtures and anything downstream doing the
same), but document them as the escape hatch they are, and revisit once a real
constructor removes the need.

### 7.4 The deprecated shims promise a version that will not exist

`calc`, `formats` and `tools_spc` emit a `DeprecationWarning` saying they will
be "removed in 0.2". The next release is 1.0. Either they go before the freeze,
or the promise is rewritten — they cannot silently survive into a version whose
whole point is that it keeps its promises.

---

## 8. Reopening this ADR

Expected triggers, in order of likelihood:

1. **Hyperspectral maps** → reopens §2.1 (numpy vs xarray).
2. **The unit space stops being closed** → reopens §2.2 (native vs `pint`).
3. **A GUI or batch runner needs to replay recipes** → builds §6.1, within the
   constraint that `.history` is already frozen.
4. **Tester feedback before the freeze** → may reopen anything in §7, which is
   why §7 exists and why the alpha went out before 1.0.
