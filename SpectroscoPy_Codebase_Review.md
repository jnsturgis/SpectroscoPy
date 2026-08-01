# SpectroscoPy — Codebase Review & Refactoring Plan

*Response to the brief in Roadmap §12. Order of operations as requested: inventory first (§1),
layer-boundary map (§2), data-model gap analysis (§3), then a sequenced plan (§4). No code has
been changed. GUI and technique-specific processing are deliberately out of scope.*

Reviewed on 2026-08-01 against commit `af883c8`.

---

## 0. What was read

**Code:** `spectroscopy/` (spectra.py, messages.py, main.py, simple_example.py), `formats/`
(jcamp.py, csv.py, spy.py), `calc.py`, `tools_spc/ftir_sidechains.py`, `README.md`, `TODO.txt`.

**Notebooks:** 64 found under `~/Documents/Research/` (the path in your brief is
`Notebook/`, singular, not `Notebooks/`), plus `Students/` and `Collaborations/`. 25 are
spectroscopy; the rest are sequence/structure/systems-biology work and are excluded. Note there
are near-duplicate pairs (`Notebook/Figure_PspA_structure.ipynb` vs `Notebook/2026/…`,
`Letitia_phage/Elina_FTIR.ipynb` vs `Untitled.ipynb`, `FTIR_sugars` in two places) — the
inventory treats each pair as one entry.

---

## 1. Notebook Inventory (Roadmap §11)

### 1.1 The inventory table

Ordered newest first. "Lib?" = does it use `import spectroscopy`.

| Notebook | Date | Technique | Raw input | Operations, in order | Lib? |
|---|---|---|---|---|---|
| `Research/FTIR_040326` | 2026-03 | ATR-FTIR | `.dpt` tab-sep, no header | load ×36 → group by (sample, wet/dry) → average replicates → crop 900–1800 → **rubberband** baseline → buffer subtract (hand factor 0.98–1.01) → area-normalise → 2nd-deriv peak pick → **amide I Gaussian fit weighted jointly on A and d²A** → plot | ✗ |
| `Notebook/2026/Figure_PspA_structure` | 2026-02 | ATR-FTIR | `.dpt`, `.txt` | as above + PCA/NMF + multi-Gaussian `curve_fit` | ✓ |
| `Notebook/2025/12/Figure_CK_111225` | 2025-12 | ATR-FTIR | 23 × `_av.txt` tsv | load → crop 900–1800 → rubberband → normalise to max → **PCA** → **NMF (k=6)** → per-component area fractions → **ANOVA + Tukey HSD + Cohen's d** → publication figure (Fig. S5) | ✗ |
| `Notebook/2025/10/Biofilm_CK_241125` | 2025-11 | ATR-FTIR | `.dpt` + `_av.txt` | load ×24 → `Average()` helper (water factor + index range + savetxt) → PCA → NMF → NMF **stability testing** (repeated runs + bootstrap, Hungarian matching) → ANOVA/Tukey | ✓ |
| `Notebook/2025/10/Figures_CK` | 2025-10 | ATR-FTIR | `.dpt` | same family as above | ✓ |
| `Notebook/2025/10/FTIR_sugars` | 2025-10 | ATR-FTIR | 43 + 11 + 44 + 21 `.dpt` | load → set sample/type/reference → grid-plot by sample → average replicates **by hard-coded index** → crop → water subtract (per-sample factor table) → rubberband → normalise to max in 1050–1080 → 2nd-deriv peaks → `savetxt` → reload the `.txt` later | ✓ |
| `Notebook/2025/10/PeaksDisplay` | 2025-10 | ATR-FTIR | `.txt` tsv | load → 2nd-deriv peaks → `find_peaks_cwt` cross-check → **manually inject known peaks** → seed multi-Gaussian `curve_fit` → area table → annotated plot | ✓ |
| `Notebook/2025/Membrane_Analysis_Roseobacter` | 2025-08 | ATR-FTIR **+ UV-Vis** | `.dpt`, `.csv` | FTIR: load ×37 → average by explicit sum/N → buffer subtract (0.991–0.998) → overlay. UV-Vis: `np.loadtxt` → build `Spectrum` by hand → `.baseline('POLY', coeffs)` scatter correction → `.smooth('SG',[2,15])` | ✓ |
| `Notebook/2025/Chloe` | 2025-10 | **Fluorescence EEM** + UV-Vis + FTIR | vendor `.csv` (col-indexed), `.dpt` | UV-Vis: 3 columns → baseline subtract. **EEM**: 35 excitations × 136 emissions → mask Rayleigh (\|λx−λm\|<15) and 2nd-order (\|2λx−λm\|<30) → contour map. FTIR: crop → rubberband → savgol smooth → normalise → 2nd-deriv peaks → band annotation → difference spectrum | ✓ |
| `Notebook/2025/AqpZ_Lipid` | 2025 | ATR-FTIR | `.dpt` via `np.loadtxt` | load → savgol → `curve_fit` ×8 → plots | ✓ |
| `Students/Candice Gomez/FTIR_200525` + `Untitled` | 2025-05/06 | ATR-FTIR | `.dpt` via `genfromtxt` | load ~25 detergent/buffer spectra → water subtract → **Beer-Lambert check: A at 8 fixed wavenumbers vs concentration, log-log** → pH titration series → `save_as('toto.spy')` round-trip test | ✓ |
| `Notebook/2025/FTIR_Tests` | 2025 | ATR-FTIR | `.dpt` | smoke test of `spc.Spectrum` | ✓ |
| `Notebook/2025/Repurification` | 2025 | UV-Vis | `.csv` | load → `curve_fit` → plot | ✓ |
| `Collaborations/Letitia_phage/FTIR_analysis2`, `Protein_FTIR` | 2025-05 | ATR-FTIR | **JCAMP `.dx`** | `import formats.jcamp` → per-sample PEG + water-vapour subtraction with a hand-tuned factor table → savgol 2nd derivative → overlay | ✓ (formats only) |
| `Collaborations/Letitia_phage/Elina_FTIR` / `Untitled` | 2025-05 | ATR-FTIR | `.dpt` | `genfromtxt` ×13 → subtract → `curve_fit` | ✗ |
| `Collaborations/Letitia_phage/FTIR_analysis` | 2024-04 | ATR-FTIR | JCAMP `.dx` | **whole `jcamp.py` pasted into the notebook** → reference spectra (water vapour, PEG, buffer) → 9-sample correction table → savgol 2nd deriv → 3 figures | ✗ |
| `Notebook/2024/SMALP_20240523/SMALP_purification` | 2024-05 | UV-Vis + FTIR + FPLC | ÄKTA `.csv` (utf-16!), UV `.csv`, JCAMP `.dx` | FPLC trace → fraction table → **A260/A280 two-component deconvolution to [SMA] and [protein]** → pooled-sample spectra → FTIR with **guide-point polynomial baseline** → region zooms. Contains an explicit *"How I would like it to work"* API sketch | ✓ |
| `Collaborations/Luca_Sapienza/GFP_binding` | 2023-11 | **Fluorescence** | vendor `.csv` wide format | emission spectra by column name → peak window 507.5–512.5 → average → subtract matched blank → **log-log standard curve** → error boxplots | ✗ |
| `Purple Bacteria/Roseobacter/chromatography/LH2_240115` | 2024-01 | UV-Vis / FPLC | ÄKTA `.csv` | pandas → column select → plot | ✗ |
| `Articles/…/Sevde/FigsNotebook`, `Students/Elissa Chams/Figures` | 2025 | mixed | `.csv` | figure assembly only | ✗ |

### 1.2 The canonical workflow

19 of the 25 notebooks are the same pipeline with different data. This *is* the library's spec:

```
load N replicate files          (.dpt tab-sep no header | JCAMP .dx | vendor .csv)
  → tag with sample name / technique / reference
  → group by sample, average replicates
  → subtract a reference (water, buffer, PEG, water vapour) × a hand-tuned scalar
  → crop to a region of interest      (900–1800, 1200–1800, 2600–3200 …)
  → baseline correct                  (rubberband ≫ polynomial ≫ ALS)
  → normalise                         (max in window | area in window | global max)
  → 2nd-derivative peak detection     (-savgol(deriv=2) → find_peaks)
  → { annotate a plot | Gaussian decomposition | PCA/NMF across samples }
  → export .txt / .pdf
```

Two things stand out. First, **the hand-tuned reference factor is the scientific crux** — the
0.99/0.98/0.995 numbers are chosen by eye, per sample, and re-typed constantly; that is the
single most valuable thing to make interactive and recordable. Second, **the multivariate
branch (PCA/NMF + stability testing + ANOVA) is your most recent and most developed analysis
and is 100 % outside the library** — it operates on bare numpy matrices assembled by hand.

### 1.3 Operation frequency (the §4 processing backlog, ranked by actual use)

| Rank | Operation | Notebooks | Current status |
|---|---|---|---|
| 1 | Load tab-sep no-header `.dpt` | 15 | works only via the `'tsv'` path — **and it is broken, see §3.2 D1** |
| 2 | Crop to x-region (boolean mask) | 19 | `Spectrum.clip(region)` exists but mutates; nobody uses it |
| 3 | Rubberband (ConvexHull) baseline | 12 | in `Spectrum.baseline('RB')` **and** pasted into 12 notebooks |
| 4 | 2nd-derivative peak detection | 14 | `Spectrum.peaks()` is an empty stub (`pass`) |
| 5 | Replicate averaging | 14 | done by `(s[0]+s[1]+s[2])/3` with hard-coded indices |
| 6 | Scaled reference subtraction | 13 | `subtract_reference()` exists, mutates, no scale factor |
| 7 | Normalise (max/area in a window) | 11 | absent |
| 8 | Annotate peaks on a plot | 10 | absent (6-line loop, always identical) |
| 9 | Export x/y as tsv | 8 | `csv.write` exists; notebooks use `np.savetxt` instead |
| 10 | Savitzky-Golay smoothing | 8 | `Spectrum.smooth('SG', [order, n])` exists, used twice |
| 11 | Multi-Gaussian peak fitting | 7 | absent; hand-rolled 3 different ways |
| 12 | PCA / NMF / ICA across a set | 5 | absent; no `SpectrumCollection` to hang it on |
| 13 | Polynomial baseline | 4 | `Spectrum.baseline('POLY', coeffs)` — takes *coefficients*, not guide points |
| 14 | JCAMP `.dx` load | 4 | `formats.jcamp` works; not reachable from a `.dx` extension (§3.2 D5) |
| 15 | ANOVA / Tukey / Cohen's d on components | 3 | out of scope for `core`, but needs a place to live |
| 16 | EEM Rayleigh/2nd-order masking | 1 | absent (Roadmap §4b, correctly deferred) |
| 17 | A260/A280 → concentration | 1 | absent (Roadmap §4b) |

### 1.4 Copy-paste hotspots — strongest signal for what belongs in the library

- **`rubberband()` — 12 verbatim copies.** Identical 40-line function, including the comments.
  It is *also* already implemented inside `Spectrum.baseline('RB')`, but no notebook uses that
  version. That gap is itself informative: the method form is harder to reach for than a free
  function, because it returns a `Spectrum` you then have to subtract, and because
  `spectra.py`'s copy has no docstring saying `'RB'` exists (the error message only mentions
  `'POLY'`).
- **`asls_baseline()` — 8 copies, 0 calls, and it cannot run.** It references `sparse` and
  `spsolve`, which none of the 8 notebooks import. You clearly *want* ALS (it is in Roadmap §4)
  but have never had a working one. A tested `method="als"` would be new capability, not a port.
- **`poly_baseline()` — 8 copies**, nearly unused, because it returns a `Polynomial` object
  rather than y-values, so it does not compose with anything.
- **The peak-annotation loop — ~10 copies**, always
  `ax.plot(x[peaks], y[peaks], "ro", ms=4)` + a `for` loop writing `f"{px:.0f}"` rotated 270°.
- **`find_nearest(array, value)` — 3 copies** (as `find_wn` in Candice's).
- **The 2nd-derivative peak call — 14 copies** of
  `-savgol_filter(y, 10, 3, deriv=2)` then `find_peaks(…, height=…, distance=10, prominence=…)`.
  Note `window_length=10` with `polyorder=3` is an even window; it works but is almost certainly
  a slip for 11, and it has propagated everywhere by copying.

### 1.5 Friction log — what the API should make *easier*

Ranked by how often the pain shows up:

1. **Replicate averaging by hard-coded index.** `(spectra[32]+…+spectra[37])/6` in FTIR_sugars,
   `range(29,32)` in the alginate cell, `Average(-0.0, range(0,3), …)` in Biofilm_CK. Every one
   of these breaks silently if a file is added. You already wrote the wish in a comment:
   *"This should be simpler and neater call"* and *"Should be based on sample name and list of
   numbers"*. **A `SpectrumCollection` with `.group_by('sample').mean()` removes this entirely.**
2. **Crop is done by mutating both arrays through a mask**, repeated 4 lines at a time, and it
   silently desynchronises a spectrum from its reference if you crop one and not the other.
   Biofilm_CK even has manual `print(len(H2O.x_data)); print(len(spectra[0].x_data))` guards.
3. **No axis compatibility check.** Arithmetic on differently-sampled spectra raises a raw numpy
   broadcast error (verified, §3.2 D3), and arithmetic on *same-length but differently-positioned*
   axes silently gives nonsense. `resample()` exists but is never called because you have to
   invoke it manually.
4. **`.set_sample()` / `.set_reference()` are called, then the result is thrown away** —
   nearly every average loses its metadata and gets `metadata['sample']` re-assigned by hand.
   §3.2 D2 shows why that re-assignment also corrupts the source spectrum.
5. **Round-tripping through `np.savetxt` + reload as `'tsv'`** is used as the interchange
   format between notebooks (FTIR_sugars → PeaksDisplay → Figure_CK). Every hop loses units,
   sample name, and the entire processing history. `.spy` exists to fix this but is used once,
   in a test cell — and does not round-trip (§3.2 D4).
6. **Reversed x-limits typed by hand** (`set_xlim((1800,900))`) in every single FTIR figure,
   along with the `r"Wavenumber (cm$^{-1}$)"` label — even though `KNOWNSPECTYPES` already
   knows the label and the technique implies the direction.
7. **The hand-tuned subtraction factor has no home.** It lives in a bare list
   (`peg_correction = [0.6, 0.0, …]`, `waters = [0.0, 0.0, 0.002, …]`) positionally aligned to a
   file list. This is exactly the kind of parameter `ProcessingStep(params=…)` exists to capture.
8. **`np.NaN` in SMALP_purification** — removed in NumPy 2.0; that notebook no longer runs.

---

## 2. Layer-Boundary Map (Roadmap §1 / §12.2)

Target: `io → core → processing → viz`, no backward reaches.

| File | Layers it currently occupies | Verdict |
|---|---|---|
| `spectroscopy/spectra.py` | **core + io + processing + viz, all four** | The central problem |
| `formats/jcamp.py`, `csv.py`, `spy.py` | io | Right layer, wrong direction of dependency |
| `calc.py` | processing (lineshapes) | Right layer, wrong package (top-level, not `spectroscopy.`) |
| `tools_spc/ftir_sidechains.py` | processing (FTIR-specific) + cli + viz | Mostly fine; `main()` mixes in plotting |
| `spectroscopy/messages.py` | cross-cutting | Holds error strings for `ftir.sidechains`, a module that isn't in the package |

### Specific crossings

**C1 — `core` reaches down into `io` (the important one).**
`spectra.py:33` does `from formats import jcamp, csv, spy`, and `Spectrum.__init__` dispatches on
file type. Meanwhile `formats/csv.py:11` and `formats/spy.py:17` do `import spectroscopy as spc`.
That is a circular import between the two layers; it currently resolves only because of module
ordering. Roadmap §3 wants a **reader registry** — that inverts this edge and breaks the cycle.

Concretely, the format dispatch is spread across four places that must be kept in sync by hand:
`FILE_EXTS` (line 40), `KNOWNFILETYPES` (line 45), the `match` in `reload()` (line 524), and the
`match` in `save()` (line 551). They already disagree — see §3.2 D5.

**C2 — `viz` inside `core`.** `Spectrum.plot(ax, *args)` (line 265) puts matplotlib in the data
model. It is used constantly in the notebooks so it should stay reachable, but per Roadmap §1 it
belongs in `spectroscopy.viz` and can be attached to `Spectrum` as a thin delegator.

**C3 — `processing` writes into `metadata`, not into `history`.** `Spectrum.peaks()`'s docstring
says *"The resulting peak list is stored in the metadata as 'peaks' or 'troughs'"*. Using
`metadata` as the output channel for analysis results is the pattern Roadmap §2.3 rules out — it
mixes acquisition facts with derived results, and it is the direct blocker to `PeakTable`.

**C4 — `io` reaches up into `core`.** `formats/*.py` `read()` takes a `my_spectrum` and mutates
it in place, rather than returning a `Spectrum`. That means a reader cannot be written or tested
without the core class, and `read(f, spectrum)` cannot report "this file has 3 spectra in it"
(relevant: JCAMP compound files, and the wide vendor CSVs in Chloe/GFP_binding where **one file
holds 35–70 spectra**).

**C5 — `formats/__init__.py` prints on import.** `print("Formats __init.py__ run.")` fires in
every notebook. A library should not write to stdout at import.

**Not a crossing, but worth flagging:** every notebook starts with
`sys.path.append('/home/james/src/Spectroscopy1.0')`. There is no `pyproject.toml` or `setup.py`
in the repo at all, despite Roadmap §8 assuming one exists. Phase 0 is genuinely from zero.

---

## 3. Data Model vs. the §2.5 Sketch (Roadmap §12.3)

### 3.1 Field-by-field

| §2.5 requires | Current `Spectrum` | Gap |
|---|---|---|
| `x`, `y` as properties | `x_data`, `y_data`, public attributes | Renaming later is breaking; notebooks touch `.x_data` directly ~200 times |
| `x_unit`, `y_unit` **explicit** | `x_label`, `y_label` — single strings mixing quantity, unit *and LaTeX*: `'Wavenumber (cm$^{-1})'` | No machine-readable unit. `.to("nm")` is impossible. Note the label in `KNOWNSPECTYPES` is missing its closing `$` |
| `technique` enum | `metadata['spec_type']`, free string, set via `set_type()` | Close; needs promotion to a field |
| `metadata` dict | ✓ `metadata` | Present, but doubles as the results channel (C3) and is aliased (D2) |
| `uncertainty` | absent | Deferred is fine — but note you *have* replicates, so a mean+SEM is naturally available |
| **`history`** | **absent entirely** | The core gap |
| Immutable, returns new | Mixed: `smooth`/`baseline`/`resample` return new; `clip`/`subtract_reference`/`set_*` mutate | Inconsistent |
| `to()`, `crop()`, `normalize()`, `find_peaks()`, `fit_peaks()` | absent (`peaks()` is `pass`) | |
| `Spectrum.read()` classmethod | constructor overloading on `*args` | See D5 |
| `SpectrumCollection` | absent | Highest-value missing type (§1.5 friction #1) |

### 3.2 Confirmed defects

Each of these was reproduced against the current tree; they matter because they affect data
already in your notebooks, not just future code.

**D1 — Loading a `.dpt` as `'tsv'` drops the first data point and corrupts both axis labels.**
This is the single most-used code path in the whole project (83+ call sites).
`csv.read` defaults to `skiprows=1`, so line 1 of a headerless `.dpt` is consumed as a header:

```
file has 4 points          →  loaded: 3
x_data                     →  [1001. 1002. 1003.]     (1000.0 lost)
x_label                    →  '1000.0'
y_label                    →  '1.0000'
```

Consequence: every `.dpt` spectrum in your notebooks is missing its first point, and any figure
that took its axis label from the spectrum would have shown a number. The plots look right
because the notebooks always set `set_xlabel()` by hand — which is *why* friction item #6 exists.

**D2 — Arithmetic results share the `metadata` dict with the left operand.**
`Spectrum.__init__` copy branch does `self.metadata = other.metadata` (line 102) — no copy, while
`x_data`/`y_data` *are* copied. So:

```python
avg = (a + b) / 2.0
avg.metadata['sample'] = "average"
a.metadata            # → {'sample': 'average'}   ← a was silently renamed
```

This fires in FTIR_sugars, Membrane_Analysis and Biofilm_CK, all of which do exactly
`X = (spectra[i]+spectra[j])/2; X.metadata['sample'] = …`. The sample labels on the *source*
spectra in those notebooks are not what they appear to be.

**D3 — No axis-compatibility handling.** `a + c` with different lengths raises a bare
`ValueError: operands could not be broadcast together`. With equal lengths but different x
positions it silently produces a wrong answer. The docstrings say `# Check that self and other
are compatible` and `# Resample data if necessary` — both are comments, not code. `resample()`
exists and would solve it.

**D4 — `.spy` does not round-trip.** `spy.write` writes name and labels; `spy.read` never parses
them (it skips to the `#` markers). Verified:

```
name    : 'MySpectrum'              → '/tmp/rt.spy'    (falls back to the filename)
x_label : 'Wavenumber (cm$^{-1})'   → 'Wavelength (nm)' (the empty-constructor default)
```

x/y and metadata survive. Since `.spy` exists precisely to be the lossless format (TODO.txt item
1), this needs fixing before it is used for anything real.

**D5 — Format dispatch tables disagree with each other and with reality.**
- `FILE_EXTS` maps `'.DX0'` (almost certainly a typo for `.dx`) but **not** `.dx`, `.jdx`,
  `.dpt`, or `.txt` — so `Spectrum("FdExperiment.04.dx")` raises `TypeError: Unknown filetype
  unknown`, even though `formats/jcamp.py` reads it. That is why SMALP_purification and the
  Letitia notebooks bypass `Spectrum` and call `formats.jcamp` directly.
- Extension matching is case-sensitive: `.DPT`, `.JDX`, `.CSV` all fail.
- `reload()` has a `case 'dpt':` branch, but `'dpt'` is not in `KNOWNFILETYPES`, so
  `Spectrum(path, name, 'dpt')` is rejected by the constructor before it can ever be reached.
  Dead branch — and the reason everyone types `'tsv'` and hits D1.
- **`reload()` and `save()` handle different sets of types** — `reload()` knows five
  (`jcamp`, `csv`, `tsv`, `dpt`, `spy`), `save()` only four (no `dpt`). An unhandled type falls
  straight through the `match` with the file already opened for writing, so it **silently
  truncates the target to zero bytes**, verified:

  ```
  save_as(..., 'tsv') ->  69 bytes
  save_as(..., 'csv') ->  69 bytes
  save_as(..., 'dpt') ->   0 bytes   <-- silently empty, no error
  ```

  Today this is latent, because `'dpt'` is unreachable through the constructor. It becomes a
  live data-loss bug the moment Phase 0.5 adds a working `'dpt'` type. **Fix both `match`
  statements together, and give `save()` a `case _:` that raises.**

  *(Correction: an earlier draft of this review claimed `save()` had no `'tsv'` case and that
  `'tsv'` spectra could not be re-saved. That was wrong — `'tsv'` writes correctly. The real
  defect is the reload/save asymmetry above.)*

**D6 — `subtract_reference` builds a tuple, not a string.**
Lines 254–256 use implicit concatenation of three f-strings inside `( … )` separated by commas,
so `metadata['reference']` becomes a 3-tuple. It also raises `KeyError` if `'reference'` was
never set, and mutates in place. (It happens not to fire in the notebooks because they use
`a - b * factor` instead — which is itself the signal that this method has the wrong shape: it
takes no scale factor, and the scale factor is the whole point.)

**D7 — Minor:** `resample()` names its result `… ' baseline'` (copy-paste from `baseline()`);
`__repr__` returns `None`; `Spectrum.peaks()` is entirely commented out and returns `None`;
`simple_example.py` calls three functions that do not exist (`spectra.write_csv_file`,
`spectra.display`, and `Spectrum(…, "xml")`), so the documented "hello world" fails.

### 3.3 Blockers to structured `ProcessingStep` (Roadmap §2.3 — settled requirement)

You asked specifically whether existing processing could log structured steps. Assessment:

**Good news — the shape is mostly right.** `smooth()`, `baseline()` and `resample()` already
(a) take a method name plus a parameter bundle, and (b) return a new `Spectrum`. That is exactly
the `(step_name, params)` pair `Pipeline.from_history()` needs. `smooth('SG', [2,15])` maps to
`ProcessingStep(name="smooth", params={"method":"SG","polyorder":2,"window":15})` almost
directly.

**Four things block it:**

- **B1 — Positional parameter lists.** `parameters` is an unnamed list whose meaning is
  method-dependent and, in `smooth`, **reversed relative to scipy** (`parameters[1]` is the
  window, `parameters[0]` the order). Serialising `[2, 15]` into history is not self-describing
  and cannot be replayed safely. Params must become a `dict` before any history is written.
- **B2 — `baseline()` returns a baseline instead of a corrected spectrum.** The user must then
  write `s - s.baseline(...)`. That subtraction is a plain `__sub__`, so the *correction* — the
  thing worth recording — is invisible to history. Either add `baseline_correct()` that returns
  the corrected spectrum and logs, or have the arithmetic propagate provenance.
- **B3 — Mutating methods have no return value to attach history to.** `clip()`,
  `subtract_reference()` and the `set_*()` family return `None`. Any of these in a chain breaks
  the record. `clip()` in particular is the crop operation, which appears in 19 notebooks.
- **B4 — Arithmetic is the most-used "processing" operation and is untracked.** The scaled
  reference subtraction `sample - 0.995 * buffer` *is* a processing step scientifically, and it
  carries the hand-tuned factor that most needs recording. `__sub__` currently copies `x_data`
  and `y_data` and nothing else. If history is added, `__add__`/`__sub__`/`__mul__`/`__truediv__`
  must append to it too — otherwise `Pipeline.from_history()` will reproduce a pipeline that
  skips the most important step.

**None of these is expensive to fix now.** All four get much more expensive once real `.history`
logs exist on disk, exactly as Roadmap §2.3 warns.

---

## 4. Proposed Refactoring Plan

Sequenced per Roadmap §9. Each phase is checkable against a named notebook before the next
starts. Nothing here touches the GUI or technique-specific corrections.

### Phase 0 — Foundations *(no behaviour change)* — **first, per §5.2**

0.1 Add `pyproject.toml` (hatchling + hatch-vcs), so `pip install -e .` replaces
`sys.path.append` in 17 notebooks. **This is a prerequisite for everything else** and does not
exist yet.
0.2 Move `calc.py` → `spectroscopy/lineshapes.py`; move `formats/` → `spectroscopy/io/`;
`tools_spc/ftir_sidechains.py` → `spectroscopy/processing/ftir.py` + a thin CLI entry point.
Keep the old import paths working with a deprecation shim for one release, so existing notebooks
do not break mid-analysis.
0.3 Remove the `print()` from `formats/__init__.py` (C5).
0.4 ruff + pre-commit + a GitHub Actions skeleton with two real tests.
0.5 Licence: **MPL-2.0** (§5.5). `LICENSE`, `license = "MPL-2.0"` in `pyproject.toml`, MPL header
in each source file, and correct the GPL v3 claim + the "2023-2924" typo in `messages.py`.

*Checkpoint:* every existing notebook still runs after deleting its `sys.path.append` line.
⚠️ D1 is still live during this phase — see §5.2.

### Phase 0.5 — Fix the confirmed defects *(before any redesign)*

Do D1–D6 as small independent commits with a regression test each. Rationale: D1 and D2 affect
data in analyses you have already published figures from, and each fix is a few lines. It would
be wrong to carry them into a new architecture, and wrong to discover them later and not know
which figures were affected.

D1 in particular deserves a `'dpt'` file type that actually works (`skiprows=0`, tab delimiter,
no header) rather than continuing to route through `'tsv'`.

*Checkpoint:* re-run `Notebook/2025/10/FTIR_sugars.ipynb` and confirm the first data point
reappears and no figure changes materially.

### Phase 1 — Core spike: ATR-FTIR only

The technique with by far the most notebook material (19 of 25).

1.1 `Spectrum` with `x`/`y` — **hard rename, no aliases** (§5.1) — `x_unit`/`y_unit` separated
from display labels, `technique` as a field. Codemod the repo (90 hits) in the same commit.
Because the label/unit split costs nothing in notebooks, do it here rather than deferring to
Phase 2.
1.2 `ProcessingStep(name, params: dict, timestamp)` and `Spectrum.history` — **including on the
arithmetic dunders** (B4). Everything immutable-returning-new; `clip` becomes `crop`, returning
new (B3).
1.3 The ~8 operations at the top of the §1.3 frequency table, with dict params (B1):
`crop`, `normalize(method="max"|"area", window=…)`, `smooth`, `derivative`,
`find_peaks(method="second_derivative")` → a `PeakTable` dataclass. Baselines per §5.3/§5.4:
`baseline(method="poly"|"rubberband"|"als", points=… | coefficients=…, return_points=…)`
returning the baseline, plus `baseline_correct(...)` with the same signature returning the
corrected spectrum and logging the step.
1.4 `SpectrumCollection` with `from_files(glob)`, `.group_by('sample')`, `.mean()`,
`.to_matrix()` — the last one is what PCA/NMF needs.
1.5 A working ALS baseline, since the copied one has never run.

*Validation (this is the point, not a formality):* re-implement two notebooks against the new
core — `Notebook/2025/10/FTIR_sugars.ipynb` (the full pipeline, replicate averaging, the
hand-tuned water factors) and `Notebook/2025/12/Figure_CK_111225.ipynb` (crop → rubberband →
normalise → matrix → PCA/NMF, and a real published figure). If FTIR_sugars' averaging block does
not collapse from ~40 lines of indexed arithmetic to a few lines of `group_by('sample').mean()`,
the collection API is wrong and should be redesigned before Phase 2.

### Phase 2 — Hardening + `0.1.0` alpha

Extend to UV-Vis, fluorescence and Raman based on Phase 1 friction. `pint`-backed conversion and
`.to("nm")`, building on the unit field split already done in Phase 1.1. Fix `.spy` round-trip
properly (D4) and make it version-aware, since it is the format that must carry `history`. Then
hand `0.1.0` to testers.

Candidate first testers from the inventory: Candice (FTIR, already using the library),
Chloé (fluorescence EEM + UV-Vis), and the Letitia phage group (JCAMP FTIR).

### Phase 3 — I/O registry

Replace the four-way dispatch (C1) with `@register_reader(extensions=[...], name=...)` returning
a `Spectrum` (fixing C4). Readers, in priority order derived from §1.3:
`.dpt` → JCAMP `.dx`/`.jdx` (already written, just unreachable — D5) → a **generic delimited-text
reader with column mapping**, which is what Chloe's 35-column EEM file, GFP_binding's wide
fluorescence CSV and the utf-16 ÄKTA exports all need → `.spy` → `.csv`/`.tsv`.

The wide-CSV case means `read_spectra()` must be able to return a `SpectrumCollection`, not just
one `Spectrum`. Worth designing in from the start.

### Phase 4 — Processing + multivariate

Everything left in §1.3, plus a `spectroscopy.processing.multivariate` module for PCA/NMF/ICA
over a `SpectrumCollection`, including the stability testing (repeated runs + bootstrap +
Hungarian matching) already written in Biofilm_CK_241125 — that is a genuinely reusable piece of
work currently trapped in a notebook.

Also here: `viz` module absorbing `plot()` (C2) with technique-aware defaults — reversed x-axis
and correct label for FTIR — plus `annotate_peaks(ax, peak_table)` (friction #6, #8).

### Phase 5+ — as per Roadmap §9. Docs, API freeze, then technique-specific, then GUI.

---

## 5. Settled Decisions

Answered 2026-08-01. Recorded here because each one changes the plan above; §4 has been updated
to match.

### 5.1 Rename `x_data`/`y_data` → `x`/`y` now, no alias period

**Settled: break it now.** The blast radius, measured rather than estimated:

| | notebooks | repo `.py` |
|---|---|---|
| `.x_data` / `.y_data` | **359** across 9 notebooks | 90 (69 in `spectra.py` alone) |
| `.x_label` / `.y_label` | **0** | 40 |
| `.fileinfo` | **0** | 8 |
| `.metadata` | 47 across 10 notebooks | 21 |

Two things this changes:

- **The rename is concentrated, not diffuse.** Four notebooks carry 313 of the 359 hits —
  FTIR_sugars (147), PeaksDisplay (73), Figures_CK (56), Biofilm_CK (37). The other five have
  ≤16 each. A regex codemod is safe here (`.x_data` → `.x` has no ambiguous matches), so this is
  an afternoon, not a project.
- **The genuinely invasive change is free.** Splitting `x_label` into `x_unit` + display label —
  the §3.1 gap that makes `.to("nm")` possible — touches **zero** notebook lines, because the
  notebooks always set axis labels by hand (friction #6). Same for `fileinfo`. So the label/unit
  redesign carries no migration cost at all and should be done in the same pass.
- `metadata` keeps its name.

*Open sub-question:* the 9 notebooks live outside this repo. I will codemod the repo; say the
word if you want the notebooks done too (recommend the four heavy ones at minimum, since they
are the Phase 1 validation targets).

### 5.2 Packaging (Phase 0) before defect fixes (Phase 0.5)

**Settled: packaging first**, as recommended — fixes land in an installable, testable state.

⚠️ **Consequence to be aware of:** D1 stays live through Phase 0. Until 0.5 lands, every `.dpt`
loaded via `'tsv'` is still missing its first data point. Phase 0 is small, but avoid starting a
new `.dpt` analysis in that window, or load with `np.loadtxt` directly as the older notebooks do.

### 5.3 `baseline()` keeps returning a baseline; add `baseline_correct()`

**Settled as recommended.** `baseline()` stays as-is for plotting and for the
`s - s.baseline(...)` idiom you already use; `baseline_correct()` returns the corrected spectrum
and appends the `ProcessingStep`. B2 was only ever a problem on the provenance path.

### 5.4 Guide-point baselines are the intended design

**Settled** — and this is more than a fix, because guide points turn out to unify all three
baseline methods. The distinction between them is just *how the anchor points are chosen*:

```python
s.baseline(method="poly", degree=3, points=[2600, 2625, ..., 3200])  # points given (SMALP)
s.baseline(method="rubberband")                                       # points found: convex hull
s.baseline(method="als", lam=1e6, p=0.01)                             # all points, asym. weights
```

The rubberband implementation *is* an automatic guide-point selector — it finds the lower hull
arc and interpolates through it. Worth exposing that: `return_points=True` should hand back the
anchors, so you can take the automatic ones, adjust them by eye, and feed them back as explicit
points. The copied notebook version already has exactly this flag (`return_hull=True`), unused —
so the notebooks anticipated the API before the library did.

The existing coefficient form is **retained**, not replaced: Membrane_Analysis uses
`baseline('POLY', [4e-7, -8.5e-4, 0.51])` for UV-Vis scatter correction with coefficients you
worked out deliberately. So `method="poly"` accepts *either* `points=` or `coefficients=`, with
`points=` as the documented default.

### 5.5 Licence: MPL-2.0

**Settled: leaning MPL-2.0.** Recording the implications so the choice is deliberate:

- **File-level copyleft.** Someone who modifies a SpectroscoPy source file must publish that
  file's changes; someone who merely imports SpectroscoPy from proprietary code does not. This
  is the middle ground between GPLv3 (currently claimed in `messages.py`) and MIT/BSD-3
  (Roadmap §8.6's recommendation), and it directly addresses §8.6's actual concern — adoption
  friction for labs running it alongside vendor tooling — while keeping improvements to the
  library itself open. For a tool you want labs to adopt *and* contribute back to, this is a
  defensible fit.
- **One caveat, not a blocker:** MPL-2.0 is less common in scientific Python than BSD-3
  (numpy, scipy, xarray are all BSD-3). It is OSI-approved and imposes nothing on downstream
  users, so it will not block `pip install`-and-use adoption. It is only a consideration if you
  later want code merged *into* a BSD-3 project.

Phase 0.5 actions: add `LICENSE` (MPL-2.0 full text), `license = "MPL-2.0"` in `pyproject.toml`,
the standard MPL header comment in each source file, and **fix the `ABOUT` string in
`messages.py`**, which currently asserts GPL v3 and also reads "Copyright (C) 2023-2924".

### 5.6 Dependency policy: numpy-native core, pandas at arm's length

**Settled: scikit-learn is an optional extra; pandas is not a dependency at all.**

```toml
dependencies  = ["numpy", "scipy", "matplotlib"]
[multivariate] = ["scikit-learn"]      # Phase 4 PCA/NMF/ICA only
```

The pandas half of this is a design decision, not just a packaging one, and it **resolves an
open question from Roadmap §2.5**. That section left `PeakTable`/`FitResult` as "lightweight
dataclasses *or* DataFrames" and sketched `to_dataframe()` on both `Spectrum` and
`SpectrumCollection`. The concern raised — *"if we use it it will get everywhere"* — is correct
and worth taking seriously: a DataFrame return type is contagious, because every caller then
needs pandas to do anything with the result, and pandas' API surface starts leaking into
docstrings, tests and tutorials.

So the rule for Phase 1 onward:

- **Core returns numpy arrays and dataclasses.** `PeakTable` and `FitResult` are dataclasses
  over numpy arrays, not DataFrames. `SpectrumCollection.to_matrix()` returns
  `(wavenumbers, X)` as plain arrays — which is exactly what scikit-learn wants anyway.
- **`to_dataframe()` survives as an explicit escape hatch**, on `Spectrum`, `SpectrumCollection`
  and `PeakTable`, importing pandas *inside the method*. If you never call it, you never need
  pandas installed.
- Nothing in `core`, `io` or `processing` may import pandas at module scope. Worth a lint rule
  or an AST test like the io/core layering check already in `tests/test_packaging.py`.

This costs one thing worth naming: the CK notebooks build their component-contribution tables in
pandas and run ANOVA/Tukey off them, so the Phase 4 statistics helpers will hand back arrays plus
a `to_dataframe()`, and those notebooks keep their own `import pandas` — which is the right place
for it, since that is analysis-specific rather than library behaviour.

---

## 6. Phase 0 — Status (done, branch `phase0-foundations`)

Delivered: `pyproject.toml` (hatchling + hatch-vcs, version resolves to
`0.0.1.dev11+gaf883c89a` from git), module moves with deprecation shims, MPL-2.0, ruff +
pre-commit + GitHub Actions, and 26 tests (22 pass, 4 `xfail(strict)` standing in for the
Phase 0.5 defect fixes). Wheel builds; `pip install -e .` done against the user site-packages
where the rest of your stack lives — reversible with `pip uninstall spectroscopy`.

**Checkpoint met.** Every import path the notebooks use resolves from outside the repo with no
`sys.path.append`: `import spectroscopy as spc`, `import spectroscopy.spectra as spec`,
`import formats.jcamp`, `from formats import jcamp`, `import calc`. The old names warn but work
and alias the *same* module objects. `spc-ftir-sidechains` is on `$PATH`.

Layer boundaries: C1 fixed (`csv.py` and `spy.py` never used their `spectroscopy` import at all;
`jcamp.py`'s single use is now function-scoped), C5 fixed, and the FTIR tool is split into a pure
`processing.ftir` plus a `cli.ftir_sidechains` that owns the argparse and matplotlib — so C2 no
longer applies to the tool, only to `Spectrum.plot`.

### 6.1 Found during Phase 0

- **`spectroscopy/gui/splash.py` exists** and was missed by the §0 inventory — my initial file
  listing was truncated by the `data/` tree. It is a June-2024 ChatGPT scratch snippet: no
  `__init__.py`, calls `open_file` before defining it, `sys` undefined, and it opens a Tk window
  at import. Excluded from the wheel and from lint; left in the tree untouched, since GUI work
  is late-phase and deleting it is your call.
- **`requires-python` is `>=3.12`, not `>=3.10`.** `subtract_reference` (defect D6) uses
  f-strings that reuse the outer quote character — PEP 701, 3.12+ only. So the package genuinely
  cannot run on 3.10/3.11 today. Declaring `>=3.10` would have been a lie that broke a tester.
  **Lowering the floor to 3.10 and widening the CI matrix to 3.10–3.13 (roadmap §6) is now a
  deliverable of the Phase 0.5 D6 fix**, and is noted as such in `pyproject.toml` and `ci.yml`.
- **`baseline('RB')` requires a `parameters` argument it then ignores** — you must write
  `s.baseline('RB', None)`. Add to the D7 minor list; `parameters` should default to `None`.
- **scikit-learn is not installed** in your user site-packages, so the PCA/NMF notebooks
  (`Figure_CK_111225`, `Biofilm_CK_241125`, `Figures_CK`) cannot currently run on this machine.
  That matters because `Figure_CK_111225` is one of the two Phase 1 validation targets — it will
  need `pip install --user scikit-learn` before it can be re-implemented against the new core.
- **D5 was partly wrong in the first draft** and has been corrected in place: `save()` handles
  `'tsv'` correctly. The real defect is that `reload()` knows five file types and `save()` only
  four, and an unhandled type silently truncates the output file to zero bytes.

### 6.2 The `x`/`y` rename — applied (Phase 1.1, first item)

`spectrum.x_data` / `spectrum.y_data` are now `spectrum.x` / `spectrum.y`, with no alias period
as settled in §5.1. Applied to both trees in one pass via `scripts/migrate_rename_xy.py`, which
handles `.py` and `.ipynb` (code-cell source only — outputs and markdown untouched):

```
repo       102 references across 7 files   (incl. tests)
notebooks  359 references across 9 files   (exactly the 9 predicted in §5.1)
```

Verified afterwards: 22 tests pass / 4 xfail; no `.x_data` or `.y_data` remains anywhere in the
package or the notebooks; all 64 notebooks in the tree are still valid JSON; all 129 code cells
across the 9 migrated notebooks still parse; and a live pipeline (load → average → mask-crop →
smooth → rubberband baseline → 2nd-derivative peak pick) runs on the renamed attributes.

`.bak` copies of all 9 notebooks sit alongside the originals — delete them once you are happy.

Two notes:

- **`scripts/` must be excluded from its own migration.** The script's regex literals
  (`r"\.x_data\b"`) contain the text `.x_data` and would be rewritten by their own pattern. The
  paths passed deliberately name `spectroscopy tools_spc calc.py tests` rather than the repo
  root.
- `x` and `y` are plain attributes for now, not properties. Making them properties — so the
  backing store can change without breaking callers — is the rest of Phase 1.1, along with the
  `x_unit`/`y_unit` split, which §5.1 established costs nothing in the notebooks.

---

## 7. Phase 0.5 — D1 fixed (`.dpt` gets a real reader)

`.dpt` is now a first-class file type with its own reader, `spectroscopy/io/dpt.py`, instead of
being routed through `csv.read` as `'tsv'`.

### 7.1 The format is not what it looked like

Surveying all the real `.dpt` files before writing the reader changed the design. They are not
uniformly tab separated:

| variant | count | consequence |
|---|---|---|
| tab separated | 681 | fine |
| **comma separated** | **140** | OPUS follows the machine locale — these **failed outright** as `'tsv'` |
| `#` metadata header | 2 | your own reference spectra, with Sample/Reference/Operator/Date/Machine |
| binary despite the extension | 1 | `eau2.0.dpt`, present in two directories |
| CRLF line endings | ~all | written on Windows |

So the separator is **sniffed per file** rather than assumed. Had I taken the `case 'dpt'` branch
that already existed in `reload()` at face value — it hardcodes `delimiter='\t'` — 140 files
would still have been unreadable.

That comma-separated set is Candice's entire `FTIR_spectra/` directory, which explains why those
notebooks use `np.genfromtxt` rather than the library: the library genuinely could not open them.

### 7.2 Verified against all 893 real files

```
new 'dpt' reader   891 ok / 2 fail
old 'tsv' route    753 ok / 140 fail
data points recovered: 753          (exactly one per file the old route could read)
```

The 2 remaining failures are the single binary file, twice, and now fail with a clear
`UnicodeDecodeError` rather than producing garbage. The `#` header is preserved verbatim in
`metadata['file_header']` rather than being interpreted — turning that block into real metadata
fields is an obvious follow-up, but guessing at it was not worth the risk here.

### 7.3 D5 came along for the ride, because D1 made it dangerous

Making `'dpt'` reachable turned a latent bug into a live one, so both had to land together:

- `FILE_EXTS` rewritten: `.dpt`, `.dx`, `.jdx`, `.txt`, `.spy` added, the `'.DX0'` typo removed,
  and lookup is now **case-insensitive**. `.dx`/`.jdx` files load through `Spectrum()` for the
  first time — the phage notebooks no longer need to call `formats.jcamp` directly.
- `reload()` and `save()` now cover the same five types, and both have a `case _:` that raises.
- **`save()` validates the file type *before* `open(..., 'w')`.** This is the important detail:
  the truncation happens at `open`, so a guard inside the `with` block would still have
  destroyed the target file before raising. There is a regression test that writes real content
  to a file, attempts an invalid save, and asserts the content survived.

### 7.4 D2 fixed — and no saved output was corrupted

The copy constructor now does `copy.deepcopy(other.metadata)` instead of a bare assignment. Deep
rather than shallow because metadata values can be containers — `metadata['file_header']` from
the new `.dpt` reader is a list — and because the constructor's own docstring already promised
*"a deepcopy of the original"*. The code simply never did it.

Because every arithmetic operator builds its result through this constructor, one line fixes all
of them; there is a parametrised test over `+ - * /`, the scalar forms, `__rsub__` and `__neg__`.

**Did this corrupt real results?** Worth answering rather than assuming, so I replayed the exact
metadata sequences from the affected notebooks under both semantics and diffed them:

| notebook | saved outputs | in-memory source labels |
|---|---|---|
| `FTIR_sugars` | **identical** — every `savetxt` filename correct | **corrupted**: `spectra[0]` `PG_coli2`→`PG_coli`, `spectra[17]` `Alginate`→`NaAlginate` |
| `Biofilm_CK` | **identical** — filenames passed explicitly | all 8 averages shared *one* dict (the water spectrum's) |
| `Membrane_Analysis` | **identical** — sources were inline temporaries | n/a |

So the conclusion is reassuring but narrowly so: **no figure or exported file is wrong.** In
FTIR_sugars the two corrupted labels were overwritten only *after* the last time they were read —
the alginate-filtering cell (`if spectrum.metadata['sample'] == "Alginate"`) runs earlier in the
notebook and saw the correct values, verified by counting matches both ways. That is luck of
cell ordering, not safety by design.

Biofilm_CK is the clearer warning. Its `Average()` helper starts each result from
`H2O_average * WFactor`, so under the old semantics **all eight sample averages shared a single
metadata dict**. Nothing broke only because that function never writes to metadata — it takes the
output filename as an argument. One `Sample1_av.metadata['sample'] = ...` would have silently
relabelled all eight.

### 7.5 Status of the review's defect list

| | | |
|---|---|---|
| D1 | `.dpt` loses first point / can't read comma files | **fixed** |
| D2 | metadata aliasing on arithmetic | **fixed** |
| D5 | dispatch tables disagree; silent truncation | **fixed** |
| D4 | `.spy` doesn't round-trip name/labels | open — Phase 2 |
| D3 | no axis compatibility checking | open — Phase 1 |
| D6 | `subtract_reference` tuple bug + PEP 701 f-strings | open — also gates `requires-python` back to 3.10 |
| D7 | assorted minor (incl. `baseline('RB')` needing a dummy arg) | open — Phase 1 |

Suite is 59 passed / 1 xfailed. **Phase 0.5 is complete** except for D6, which is bundled with
the Phase 1 rewrite of `subtract_reference` (it needs a scale-factor argument anyway, per
friction item #7).

---

## 8. Environment note — where scikit-learn actually is

It is **not** missing, it is in a conda environment:

```
/home/james/bin/miniconda3/envs/pylipid/    python 3.13, scikit-learn present
```

Meanwhile the notebooks run against `~/.local/lib/python3.12/site-packages` (numpy, scipy,
matplotlib, and now SpectroscoPy), and there is one unrelated venv at
`/hdd/james/src/3D_Simulation/.venv`. So the PCA/NMF notebooks currently cannot run in the same
interpreter as everything else — three Python environments, and the one with sklearn is not the
one with SpectroscoPy.

Since you mention wanting to rationalise this: the smallest step that unblocks Phase 1 is
`pip install --user --break-system-packages scikit-learn`, putting it alongside the rest of the
notebook stack. The tidier long-run answer is a single project venv that SpectroscoPy and the
notebooks both use, with `pip install -e .` into it — worth doing before testers arrive in
Phase 2, since "which Python?" is exactly the question you do not want them asking. Not done
here; it changes how you launch Jupyter and that should be your call.

---

## Appendix — files consulted

Notebook source was extracted to text for analysis; the extraction script and dumps are in this
session's scratchpad and are not part of the repo. Nothing outside `/hdd/james/src/Spectroscopy1.0`
was modified, and no notebook was executed.
