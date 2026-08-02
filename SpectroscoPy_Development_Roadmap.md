# SpectroscoPy — Development Roadmap

**Scope for v1:** Raman, FTIR/IR, UV-Vis, Fluorescence
**Target:** Public, open-source PyPI package
**Starting point:** existing rough package (setup.py/pyproject present)

---

## 1. Guiding Principles

These are worth pinning to the top of your README/CONTRIBUTING, because every architecture decision below traces back to them:

1. **One data model, many instruments.** The pain you're solving is "every spectrometer has its own tool." The fix isn't more format-specific tools — it's a single, well-designed in-memory representation that every format converts *into*, and every analysis tool operates *on*. Get this right once; everything else is a plugin around it.
2. **Separate concerns hard:** `io` (reading vendor files) → `core` (data structures) → `processing` (algorithms) → `viz`/`gui` (presentation). No layer should reach backward into a layer above it.
3. **Design for the format you don't support yet.** Even though v1 targets 4 techniques, the core object should not have "Raman-shaped" assumptions baked in (e.g. don't hardcode "wavenumber" as the only x-axis unit).
4. **Prefer boring, well-supported dependencies** (numpy, pandas/xarray, scipy) over clever custom code, especially in the core — this is what makes the library maintainable by others later.

---

## 2. Core Data Model (Priority #1)

This is where to spend the most design time, since it's the hardest thing to change post-1.0 without breaking every downstream user.

### 2.1 Proposed object hierarchy

```
Spectrum                # single spectrum: 1D intensity array + 1 axis + metadata
SpectrumCollection       # ordered/unordered set of Spectrum objects (batches, time series, replicates)
SpectralMap / HyperCube  # 2D/3D — spatial (x,y) + spectral axis, for mapping/imaging (Raman/FTIR mapping)
```

Consider building `Spectrum` as a thin, purpose-built wrapper rather than subclassing pandas/xarray directly — but back it internally with `xarray.DataArray` (axis = coordinate, intensity = data, metadata = `.attrs`). This gets you free interpolation, arithmetic, slicing, and serialization (NetCDF/Zarr) without reinventing them, while still giving you a clean, domain-specific API on top (`spectrum.baseline_correct()`, `spectrum.wavenumber`, etc.).

**Minimum required fields:**

| Field | Notes |
|---|---|
| `x` | axis values (wavenumber, wavelength, Raman shift, etc.) |
| `y` | intensity/absorbance/counts |
| `x_unit` | e.g. `cm^-1`, `nm` — **explicit, never assumed** |
| `y_unit` | e.g. `absorbance`, `counts`, `%T` |
| `technique` | enum: `raman`, `ftir`, `uvvis`, `fluorescence`, ... |
| `metadata` | free-form dict: instrument, integration time, laser wavelength, slit width, sample ID, timestamp, operator, original filename |
| `uncertainty` (optional) | per-point error/noise estimate, if available |
| `history` | append-only log of processing steps applied (critical for reproducibility — see §2.3) |

### 2.2 Unit handling

Strongly consider using `pint` for `x_unit`/`y_unit` from day one. Unit-confusion bugs (nm vs cm⁻¹, absorbance vs %T) are exactly the kind of thing that erodes trust in a shared tool. `spectrum.to("nm")` should be a first-class, tested operation.

### 2.3 Processing provenance

Every processing step (baseline correction, smoothing, normalization...) should return a **new** `Spectrum` (immutable-by-default, like pandas' non-inplace default) and append a record to `.history`: what was done, with what parameters, when. This is what makes the tool trustworthy for real analysis — a user (or reviewer) can always answer "how did I get from raw data to this plot?"

**Decision (settled): chaining for exploration, `Pipeline` for reuse/validation/reproducibility.** Interactive/notebook work stays method-chained (`spectrum.baseline_correct().smooth().normalize()`); the moment a recipe is worth sharing with a tester, batch-applying, or wiring into the GUI, it graduates into a `Pipeline` object via `Pipeline.from_history(result.history)` — no rewrite needed to get there. This means `.history` must carry more than a human-readable description: each `ProcessingStep` needs to record the **exact step type and parameters** (not just a log string), so it can be losslessly turned back into a runnable `Pipeline` step. Concretely: `ProcessingStep(name="baseline_correct", params={"method": "als", "lam": 1e5, "p": 0.01}, timestamp=...)` rather than `ProcessingStep(description="baseline corrected via ALS")`. Get this structured-vs-string distinction right early — it's the one thing that makes the chaining/pipeline duality actually work rather than being two disconnected APIs.

### 2.4 Architecture Decision Record

Before writing more code against this model, write a short ADR (`docs/adr/0001-core-data-model.md`) capturing: chosen backend (xarray vs custom), unit strategy, mutability policy, and what you explicitly deferred (e.g. multi-channel/hyperspectral support). This becomes the reference every contributor (including future-you) checks before proposing a breaking change to `core`.

### 2.5 Draft core API sketch

A concrete (deliberately incomplete) sketch, useful as a target to review the existing skeleton against — not a spec to implement verbatim:

```python
class Spectrum:
    def __init__(self, x, y, x_unit, y_unit, technique, metadata=None,
                 uncertainty=None, history=None):
        ...

    # -- properties --
    @property
    def x(self) -> np.ndarray: ...
    @property
    def y(self) -> np.ndarray: ...
    @property
    def technique(self) -> str: ...

    # -- units --
    def to(self, x_unit=None, y_unit=None) -> "Spectrum":
        """Return a new Spectrum converted to the given unit(s)."""

    # -- processing (each returns a NEW Spectrum, appends to .history) --
    def baseline_correct(self, method="als", **kwargs) -> "Spectrum": ...
    def smooth(self, method="savgol", **kwargs) -> "Spectrum": ...
    def normalize(self, method="area") -> "Spectrum": ...
    def crop(self, x_min=None, x_max=None) -> "Spectrum": ...
    def resample(self, new_x) -> "Spectrum": ...

    # -- analysis --
    def find_peaks(self, **kwargs) -> "PeakTable": ...
    def fit_peaks(self, model="gaussian", n_peaks=None, **kwargs) -> "FitResult": ...

    # -- provenance --
    @property
    def history(self) -> list[ProcessingStep]: ...

    # -- I/O --
    @classmethod
    def read(cls, path, **kwargs) -> "Spectrum": ...
    def to_csv(self, path) -> None: ...
    def to_netcdf(self, path) -> None: ...

    # -- interop --
    def plot(self, ax=None, **kwargs): ...
    def to_dataframe(self) -> pd.DataFrame: ...


class SpectrumCollection(Sequence[Spectrum]):
    """Ordered collection with batch operations, e.g.:
    coll.baseline_correct(...)  # applies to all, returns new SpectrumCollection
    coll.mean() -> Spectrum
    coll.to_dataframe() -> pd.DataFrame  # long or wide format
    """
```

One design question resolved, one still open:

- **Chaining vs. `Pipeline` — resolved: both.** Method chaining (`spectrum.baseline_correct().smooth().normalize()`) is the ergonomic surface for exploratory/notebook work; `Pipeline([BaselineCorrect(), Smooth(), Normalize()])` is the serialization/reuse layer for sharing recipes with testers, batch-applying, and driving the GUI. `Pipeline.from_history(spectrum.history)` converts one into the other, which is why `.history` must store structured `(step_name, params)` pairs, not free-text descriptions (see §2.3). A minimal `ProcessingStep` implementation should be part of the Phase 1 core spike, not deferred — it's cheap now and expensive to retrofit once real `.history` logs exist.
- **`PeakTable`/`FitResult` return types (still open):** worth designing now as lightweight dataclasses/DataFrames rather than ad hoc dicts, since these are exactly the outputs your notebooks probably already compute by hand.

---

## 3. I/O Layer

Design as a **registry of format readers**, not a chain of if/elif on file extension.

```python
@register_reader(extensions=[".spa"], vendor="thermo_ftir")
def read_thermo_spa(path) -> Spectrum: ...

read_spectrum("sample.spa")  # dispatches automatically
```

- Each reader's only job: parse vendor bytes → populate a `Spectrum` with correct units/metadata. No processing logic here.
- Start with the formats you personally have files for (fastest feedback loop), one per technique:
  - Raman: Renishaw `.wdf`, or a generic `.txt`/`.csv` fallback
  - FTIR: Thermo `.spa`/`.spc`, Bruker `.0` (OPUS), or `.dpt`/`.csv`
  - UV-Vis: vendor CSV/TXT (these vary a lot — worth a flexible delimited-text reader with a config/sniffing layer)
  - Fluorescence: vendor CSV/TXT, possibly with excitation-emission matrix (EEM) support later
- A **generic delimited-text reader** (configurable column/header mapping) as a fallback/escape hatch is high value early, since it unblocks "weird format my colleague sent me" cases without a dedicated parser.
- Writers (`Spectrum.to_csv()`, a native format like NetCDF via the xarray backend) matter too — for reproducible intermediate saves.

---

## 4. Processing / Analysis Layer

Split into:

**a) Technique-agnostic (build first — reusable across all 4 techniques):**
- Baseline correction (ALS, polynomial, rubberband)
- Smoothing/denoising (Savitzky-Golay, moving average)
- Normalization (min-max, SNV, area, vector norm)
- Peak detection & fitting (scipy.signal + lmfit/scipy.optimize for Gaussian/Lorentzian/Voigt fits)
- Interpolation/resampling, axis alignment/registration
- Basic arithmetic (subtract reference/blank, average replicates)

**b) Technique-specific (build after core algorithms are solid):**
- Raman: cosmic ray removal, fluorescence background subtraction
- FTIR: ATR correction, CO2/H2O artifact handling
- UV-Vis: Beer-Lambert concentration calculations, derivative spectroscopy
- Fluorescence: EEM correction (inner filter effect, Raman/Rayleigh scatter removal)

This split is also a natural module boundary: `spectroscopy.processing.common` vs `spectroscopy.processing.raman`, etc.

---

## 5. GUI Application

Deliberately sequenced **after** the library core is stable — a GUI built on a shaky data model means rewriting the GUI twice.

**Recommendation:** PySide6 (LGPL, avoids GPL licensing friction for open source) over PyQt. Consider a thin GUI that is *just* a client of the library's public API — if the GUI can only do what a user could also do by importing `spectroscopy` in a notebook, you'll avoid logic duplication and the GUI stays maintainable by contributors who aren't Qt experts.

Alternative worth weighing: a browser-based GUI (Dash, Panel, or Streamlit) — lower dev overhead, easier remote/lab-server use, weaker for a "native app" feel and offline packaging. Worth a quick spike (1-2 days) comparing PySide6 vs Panel before committing, given how much downstream work depends on this choice.

MVP feature set: load file(s) → view/overlay spectra → apply a processing pipeline interactively → export processed data + figure. Batch/GUI-driven scripting (i.e., GUI actions expressible as replayable code) is a nice differentiator vs. the vendor tools you're replacing, and falls naturally out of the `.history` provenance design in §2.3.

---

## 6. Testing Strategy

- **Framework:** pytest, with `pytest-cov` for coverage tracking.
- **Core data structure tests:** property-based testing (via `hypothesis`) for `Spectrum` arithmetic, unit conversion round-trips, and history immutability — these are the functions where subtle bugs are most costly.
- **I/O tests:** small, real (anonymized) sample files per format checked into `tests/data/`, one per vendor format, with hand-verified expected values (a handful of known x/y points, not full-file diffs).
- **Processing tests:** synthetic spectra with known ground truth (e.g. a synthetic Gaussian peak on a known baseline) to verify baseline correction/peak fitting recover the right answer within tolerance.
- **GUI tests:** `pytest-qt` for widget-level tests once GUI work starts; keep GUI logic thin enough that most logic is tested via the library API instead.
- **CI matrix:** GitHub Actions, test across Python 3.10–3.13 and at minimum Linux + one of macOS/Windows (spectroscopists are disproportionately on Windows — don't skip it).
- Target coverage as a guardrail, not a goal in itself (e.g. 80%+ on `core`/`processing`, more lenient on `gui`).

---

## 7. Documentation Strategy

- **Tooling:** Sphinx + `numpydoc`-style docstrings (standard in the scientific Python ecosystem — familiar to your future users) or MkDocs + `mkdocstrings` if you prefer Markdown-first. Either is fine; Sphinx has better scientific-Python-ecosystem precedent (numpy, scipy, xarray all use it), which lowers the learning curve for contributors.
- **Structure:**
  - *Getting started* (install, load a file, make a plot in 5 lines)
  - *User guide* (one page per concept: data model, I/O, processing, GUI)
  - *Tutorials* — one worked example per technique (Raman/FTIR/UV-Vis/Fluorescence), ideally as executable notebooks via `sphinx-gallery` or `myst-nb`
  - *API reference* — auto-generated from docstrings
  - *Architecture Decision Records* (from §2.4) — kept even after decisions are "settled," as a record for contributors
- **Hosting:** Read the Docs (free, standard, auto-builds from GitHub) or GitHub Pages.
- Write the "getting started" page **before** the GUI exists — forces you to validate the library API reads naturally standalone.

---

## 8. Packaging & Deployment

Since a rough `pyproject.toml`/`setup.py` already exists:

1. **Consolidate on `pyproject.toml`** only (drop `setup.py` if both currently exist) with a modern backend: `hatchling` or `setuptools>=64` with `pyproject.toml`-only config.
2. **Versioning:** `hatch-vcs` or `setuptools_scm` for git-tag-based versioning (avoids manually bumping version strings).
3. **Linting/formatting:** `ruff` (lint + format in one tool, fast) + `mypy` for the core/processing modules at minimum. Add as `pre-commit` hooks.
4. **CI/CD (GitHub Actions):**
   - On PR: lint, type-check, test matrix
   - On tag/release: build sdist+wheel, publish to PyPI via **Trusted Publishing** (OIDC — no API tokens to manage)
5. **Changelog:** `CHANGELOG.md` in Keep-a-Changelog format, updated per PR or generated from conventional commits.
6. **Community files:** `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue/PR templates, and a clear open-source license (MIT or BSD-3 are the norm for scientific Python tools and lower adoption friction vs GPL, especially since you're hoping labs will adopt this alongside their existing tooling). Feedback channel: mailto now, relay-service later.
7. **Optional GUI packaging:** once the GUI is stable, `pyinstaller`/`briefcase` for a standalone executable — many spectroscopists won't want to manage a Python environment.

---

## 9. Phased Roadmap (spike → validate → harden)

Given that the core needs real-world testing before being locked down, this is not a waterfall — it's a sequence of small, real spikes, each validated against actual use (your notebooks first, then testers) before the next layer builds on top. Treat everything pre-Phase 6 as **0.x**: semver where breaking API changes are expected and allowed, and should be communicated as such (a `⚠️ API unstable pre-1.0` note in the README is enough).

| Phase | Focus | Validated by | Rough deliverable |
|---|---|---|---|
| **0. Foundations** | Clean packaging config, ruff/pre-commit/CI skeleton, license, notebook inventory (§11) | — | Repo builds/lints/tests in CI; a written table of what your notebooks actually do |
| **1. Core spike** | `Spectrum` for **one** technique only (whichever has the most notebook material), §2.5 sketch | Re-implement 1–2 of your own notebook analyses using it | A core that already proved itself on a real analysis, not just unit tests |
| **2. Core hardening + first alpha** | Extend `Spectrum`/`SpectrumCollection` to all 4 techniques based on Phase 1 friction, `.history`, units | Hand **0.1.0** to 1–2 testers per technique | Tagged alpha release, feedback issues open |
| **3. I/O layer** | Reader registry, one real parser per technique + generic CSV fallback | Testers loading their own files | Real files → `Spectrum` objects, for all 4 techniques |
| **4. Common processing** | Baseline, smoothing, normalization, peak fitting — prioritized by frequency in your notebook inventory | Testers running their usual workflows through it | Reusable pipeline, tested against synthetic ground truth *and* tester workflows |
| **5. Docs & tutorials v1** | Getting started + notebook-derived tutorial per technique | — | Public-readable docs site, even pre-1.0 |
| **6. API freeze candidate** | Revisit ADR-0001 in light of tester feedback; lock `core`/`io` public API | Explicit tester sign-off / burn-in period | 1.0-candidate API surface |
| **7. Technique-specific processing** | Raman/FTIR/UV-Vis/Fluorescence-specific corrections | Testers | Feature parity with current vendor tools |
| **8. GUI MVP** | Load → view → process → export, built on the now-stable library API | Testers using it hands-on | Usable desktop app for daily lab work |
| **9. Packaging & public release** | PyPI publish, trusted publishing CI, changelog, community files | Wider public | `pip install spectroscopy` works for a stranger |
| **10. Post-1.0** | More formats, hyperspectral maps, plugin system | Community | Community contributions become possible |

The key structural change from a waterfall plan: **Phase 2 puts something in testers' hands before I/O, processing, or GUI exist**, purely to validate whether the core `Spectrum` object itself is pleasant and correct to use — the cheapest possible point to catch a bad design decision.

---

## 10. Immediate Next Steps

1. Do the notebook inventory (§11) — this is now the actual first step, before any refactoring.
2. Write ADR-0001 for the core data model, informed by that inventory.
3. Set up CI skeleton (even with 2 dummy tests) so every subsequent PR is tested from day one.
4. Implement `Spectrum` for one technique and use it to redo 1–2 existing notebook analyses (Phase 1) — this *is* the validation step, not a formality before one.
5. Once that feels right, brief Claude Code using §12 below.

---

## 11. Mining the Existing Notebooks

Before any refactoring, the notebooks are the highest-value asset in the codebase — they're evidence of real usefulness, not a guess. Worth doing this as an explicit, written exercise (a markdown table or spreadsheet), not just "keeping it in mind":

For each notebook, record:
- **Technique** it covers
- **Raw input**: what format/shape the data arrived in
- **Operations performed**, in order (baseline correction, smoothing, a specific normalization, a manual peak fit, a specific plot type, a unit conversion, etc.)
- **Anything done more than once by copy-pasting** — a strong signal it belongs in the library, not a script
- **Anything that was awkward/fiddly** — a strong signal for what the new API should make *easier*, specifically

This inventory does three jobs at once:
- Prioritizes §4's processing backlog by actual frequency of use rather than guesswork
- Becomes the seed content for tutorials (§7 docs), since a tutorial that mirrors "the thing I actually needed to do last month" is more useful than an invented example
- Gives Claude Code a concrete usage baseline to check the refactored core against — "does the new `Spectrum` API make this notebook's workflow shorter/clearer?" is a much sharper review question than "is this code good?"

---

## 12. Briefing Claude Code for the Initial Review

Rather than asking for open-ended "refactoring suggestions," a more useful review has a fixed order of operations. Suggested brief to hand Claude Code, alongside this roadmap document:

1. **Inventory first, don't refactor first.** Read through the existing skeleton and all notebooks. Produce the notebook-usage table from §11 before proposing any code changes.
2. **Map the skeleton against the layer boundaries** in §1 (`io` / `core` / `processing` / `gui`/viz). Flag anything that crosses layers (e.g. a plotting call inside a processing function, a hardcoded file format assumption inside what should be a generic `Spectrum` method).
3. **Compare the skeleton's data structures against the §2.5 API sketch.** Report gaps and mismatches — don't rewrite yet. In particular, check whether existing processing functions (if any) are structured in a way that could log structured `ProcessingStep` records (§2.3) — this is a settled requirement, not an open design question, so flag any existing pattern that would block it (e.g. processing done via free functions with no return-value/logging convention at all).
4. **Propose a refactoring plan as a reviewable document/PR description first**, sequenced according to §9's phases (core spike → hardening → I/O → processing), not as a single large rewrite. This keeps each step checkable against a real notebook re-implementation before moving to the next.
5. Explicitly hold off on GUI and technique-specific processing work until told to — those are late-phase by design.

Handing over this roadmap plus that ordered brief means Claude Code's first output is a *plan you can approve or redirect*, rather than a finished refactor you have to unpick.

---

*Share the existing codebase (zip or key files, notebooks included) whenever you'd like a second opinion here first — I can also do the §11 notebook-inventory pass and layer-boundary check directly, before you hand things to Claude Code, if that'd be useful as a sanity check on the brief itself.*

---

## 13. Validation Publication (added 2026-08, James)

> **Superseded in part.** `SpectroscoPy_Dual_Publication_Workflows.md` splits
> this into two papers — a JOSS software paper and a domain application paper —
> and adds a shared prerequisite this section does not mention (joint NMF
> decomposition across two techniques). The paper described below is that
> document's Workflow B. The feature list in §13.3 stands unchanged; the framing
> of "one paper" does not. See that document's §6.4 for what still needs
> resolving between the two.

A target that reorders §9 rather than extending it: **a methods/application paper,
written once a draft GUI exists**, based on multivariate analysis of fluorescence
and FTIR spectra of *Pseudomonas aeruginosa* biofilms.

### 13.1 Why this target changes the plan

§7 of the codebase review notes that the one credibility signal SpectroscoPy
lacks and RamanSPy has is a peer-reviewed citation. This supplies it. But it
also does something more useful than credibility: **a paper is a deadline with a
fixed feature list**, and the list is not the same as the one the phases would
have produced on their own. Everything in §13.3 is there because the paper needs
it, not because it was next.

The scientific angle is the sell: **NMF is little used in biology**, and a
biofilm system — where the spectra are genuinely mixtures of a few components
in varying proportion — is close to the ideal case for it. The library already
has decomposition with bootstrap stability testing (`stability(...)`), which is
the part most published NMF applications get wrong or skip, so the methods
section has something to say beyond "we ran scikit-learn".

### 13.2 What the paper must sell, and the one thing it cannot

- **Versatility** — one object model across fluorescence and FTIR in the same
  analysis, which is exactly what the biofilm study does anyway. Demonstrable.
- **Fluidity for scripting** — chaining, provenance, `.spy` round-trip. A figure
  produced from a `describe_history()` listing makes this visible on the page.
- **GUI ease of use** — **not demonstrable in an article.** A screenshot proves
  a GUI exists, not that it is pleasant. Worth stating plainly rather than
  padding the paper with interface figures: the honest move is to make the GUI
  the thing a reader can go and try, and spend the paper's space on the analysis
  that a reader can check.

### 13.3 Features the paper needs (priority follows from this list)

Ordered by how much the paper depends on them, not by difficulty.

| | Feature | Why the paper needs it |
|---|---|---|
| 1 | **Library/reference lookup *and* decomposition** | Competitors offer spectrum-vs-library matching; none couple it to decomposition. Identify the components NMF pulls out by matching them against references, rather than assigning them by eye — that pairing is the methodological novelty, not the lookup alone |
| 2 | **`concentration()` — Beer–Lambert** | Cross-validated against **Bradford protein assays** on the same samples. This turns a convenience helper into a validation figure: spectroscopic concentration vs a wet-chemistry standard |
| 3 | **OPUS native binary (`.0`, `.1`, `.2`)** | The biofilm FTIR data arrives this way. `.dpt` export is a step the paper's readers should not have to be told to perform, and review §3 already names this as the honest I/O gap |
| 4 | **Galactic `.spc`** | The other format reviewers will ask about; `spc-spectra` exists but is unmaintained since 2018, so this is a read-it-ourselves job |
| 5 | **A260/A280 decomposition** | Openly a hoop to jump through — the number is close to worthless as a purity measure and the two-component decomposition is three lines. It is asked for often enough that its absence reads as a gap, so it costs less to have than to keep explaining |
| 6 | **Fluorescence line narrowing / hole burning** | For the fun of it, and there is old data. Compared against Raman and FTIR of the same system, it is a genuinely unusual thing for a general-purpose package to handle, and it exercises the "one data model, many instruments" claim harder than four routine techniques do |

Items 1–3 are load-bearing for the paper. Items 4–6 strengthen it; 6 in
particular is the kind of thing a reviewer remembers.

### 13.4 Dependencies and ordering

- The GUI (§9 Phase 8) moves **before** the paper rather than after public
  release — a draft is enough, it does not need to be finished.
- The **name must be settled before submission**. A paper naming a package that
  later gets renamed on PyPI is a citation that stops resolving. See review §16.
- Bradford assays on the biofilm samples are a **wet-lab dependency**, not a
  code one, and they gate item 2's validation figure rather than the feature.
- Nothing here changes the §14.1 rule: each of these still needs a real example
  on real data before it is rolled out. The difference is that the paper
  supplies the examples, because they are the paper's own figures.

---

## 14. The road to 1.0.0 and JOSS (settled 2026-08-02, James)

Three decisions, which together replace §13's framing:

1. **The two papers are decoupled.** The JOSS paper does not wait on joint
   decomposition. See the dual-publication document §6.2.
2. **The domain paper is the biofilm work, enhanced** — a rewrite of the earlier
   FTIR-only publication with new data and an extra level of complexity: paired
   fluorescence and FTIR spectra of the same samples, decomposed jointly by NMF.
   This resolves the ambiguity flagged in that document's §6.4; there is one
   biological study, not two.
3. **JOSS is submitted against shipped capabilities at 1.0.0, in about three
   months** — target **early November 2026**.

### 14.1 What decoupling actually buys

It is worth being explicit about the consequence, because it is larger than it
looks: **the road to 1.0.0 now contains almost no new features.**

Everything on §13.3's list — library lookup, `concentration()`, OPUS binary,
`.spc`, A260/A280, line narrowing — is *additive*. None of it changes an
existing signature, so all of it can ship in 1.1 without breaking anybody. It
belongs to the domain-paper track, which has its own timetable.

What 1.0.0 actually is: **a promise not to break things**. So the next three
months are a scoping, documentation and freezing exercise, not a build. Three
months is realistic for that. It is not realistic for that *plus* features, and
the main way this deadline gets missed is by letting feature work leak into it
because features are more interesting than deciding what `spectroscopy.__all__`
should contain.

### 14.2 Blockers to an API freeze, found by checking

Not a wish list — these were verified against the tree on 2026-08-02.

**1. The public namespace leaks seventeen modules.** `spectroscopy/__init__.py`
does `from .spectra import *`, `spectra.py` has no `__all__`, and star-import
therefore re-exports everything the module imported. Today `import spectroscopy`
gives you `spectroscopy.os`, `spectroscopy.copy`, `spectroscopy.np`,
`spectroscopy.warnings`, `spectroscopy.operator` — 40 public names of which 17
are modules. **Freezing this at 1.0 makes `spc.np` part of the API**, and
someone will use it. The package needs an explicit top-level `__all__` and
`spectra.py` needs one, before anything is frozen. This is also §14.4's fourth
doc trap (a test that every exported name is documented), which cannot be
written until the export list exists.

**2. ADR-0001 — written 2026-08-02**, at `docs/adr/0001-core-data-model.md` and
published with the docs. It consolidates the decisions that were scattered
through review §5 and §9.1: numpy rather than xarray, native units rather than
`pint`, structured provenance, the `.spy` format, the mutability split, the
dependency policy, and the ten-name public surface. It also records `Pipeline`,
`fit_peaks`/`FitResult`, hyperspectral maps and `to_netcdf` as **explicit
deferrals** rather than leaving them to look like oversights — all four are
additive and can arrive in 1.1.

Writing it against the code rather than against the plan turned up **two freeze
blockers nobody had listed**, now items 6 and 7 below.

**6. ✅ Fixed 2026-08-02 — a `Spectrum` can be built from arrays.**
`Spectrum(x, y, technique=..., x_unit=..., name=..., metadata=...)` and
`Spectrum.read(path, file_type=None)` both exist; the string forms still work,
and 22 tests pin the behaviour. The guide page documenting it now executes at
build time, which is the first page of the user guide to do so (review §14.4
item 1). The problem it fixed: `__init__(self, *args)`
takes nothing, a `Spectrum`, or one to three strings. Roadmap §2.5's primary
constructor — `Spectrum(x, y, x_unit=...)` — does not exist, so anyone with
computed data must construct an empty object and assign `.x` and `.y`, which is
what the test fixtures do. Freezing `*args` as the only entry point makes
empty-then-assign the API. Fix before November: add a keyword constructor and a
`Spectrum.read()` classmethod, keeping the string forms.

**7. ✅ Fixed 2026-08-02 — the docstring was corrected, not the behaviour.**
`Spectrum.read(filename, file_type)` now does what the docstring wrongly claimed
the two-argument constructor did; two positional strings still mean
`(directory, filename)`, since the docs and readers already assumed that and
changing it would break working calls to fix a comment. The finding: The docstring says
`Spectrum(filename, filetype)`; the code reads two arguments as `(path, name)`.
The documented call raises `TypeError: Unknown filetype unknown` — verified.
One of the two is wrong and the signature cannot be frozen until it is decided
which.

**3. `Pipeline` is settled but unbuilt.** §2.3 and §2.5 call it decided, and
`ProcessingStep` exists and carries structured `(name, params)` — the expensive
part. `Pipeline.from_history()` does not exist. Since it is additive it *can*
wait, but the decision must be conscious, because the `.history` format and
`.spy` serialisation freeze at 1.0 and they are what a later `Pipeline` has to
be reconstructible from.

**4. `fit_peaks()` / `FitResult` do not exist.** §2.5 left the return type open;
`find_peaks()` and `PeakTable` shipped, peak *fitting* did not. Additive, so
defer — but say so explicitly in ADR-0001 rather than leaving it looking like an
oversight.

**5. The deprecated shims say "removed in 0.2".** `calc`, `formats` and
`tools_spc` emit a DeprecationWarning promising removal in 0.2. There is now no
0.2 — the next release is 1.0. Either they go before the freeze, or that promise
is rewritten. They cannot silently survive into a version that promises
stability.

### 14.3 Sequence, backwards from early November

| When | What | Why then |
|---|---|---|
| **August, immediately** | Hand 0.1.0 to the testers | The whole critical path. An API frozen without tester friction is a guess, and §9's entire structure exists to avoid that. Nothing else in this table can start it |
| August | Define the public surface: `__all__` at top level and in `spectra.py`; audit the 40 names | Blocker 1. Independent of feedback, so it can run in parallel |
| August–September | Collect and triage friction; relay reports into issues by proxy (review §14.5) | Feedback has to be *in hand* by mid-September to inform the freeze |
| September | Make every breaking change there is going to be: shims removed, names settled, signatures fixed | Last window. After this, breaking anything costs a major version |
| September | Write ADR-0001, including the explicit deferrals (`Pipeline`, `fit_peaks`) | Blocker 2. Needs the feedback to be honest about what was validated |
| October | JOSS-facing gaps: `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, `CITATION.cff`, repository description | Cheap, mechanical, and reviewers check them by name |
| October | The §14.4 documentation traps and the §14.2 comparison pages of the review | Documentation quality is an explicit JOSS review criterion, not a nicety |
| October | Draft `paper.md` (250–1000 words). The README "Why" is the statement of need; the comparison draft §1–§5 is the state of the field | Both already written, which is why this is a small job |
| **Early November** | Tag **1.0.0**, archive to Zenodo for a DOI, submit to JOSS | |

### 14.4 The risk, named

**The schedule has one dependency it does not control: the testers.** If
feedback has not arrived by mid-September, the choice is to slip 1.0 or to
freeze an API that no outsider has exercised — and the second is the option that
looks fine now and hurts for years, because §9's whole premise is that the core
must be validated by real use before it is locked.

Mitigation is unglamorous: hand the alpha over *this week*, ask each tester for
one concrete workflow rather than general impressions, and chase in early
September rather than waiting politely. The three named in review §9.5 cover the
four techniques between them, which is the minimum for a freeze that claims to
serve all four.

---

## 16. Redox titration analysis (added 2026-08-02, James — dataset to come)

A third target application, alongside secondary structure (§15, ADR-0002) and
the biofilm paper (§13). Dataset not yet in hand.

### 16.1 It is not a special case — it is the general one

A redox titration is a **series of spectra parameterised by a continuous
variable, fitted to a physical model**. Written that way, it is the same shape
as several things this library either already half-does or will be asked for:

| Series | Parameter | Model | Yields |
|---|---|---|---|
| Redox titration | potential (mV) | Nernst | midpoint potential *E*ₘ, electron count *n* |
| pH titration | pH | Henderson–Hasselbalch | p*K*ₐ |
| Thermal melt | temperature | van 't Hoff | *T*ₘ, ΔH |
| Ligand binding | concentration | binding isotherm | *K*d, stoichiometry |
| Kinetics | time | exponential | rate constants |

Two pieces of evidence that this generalisation is the right one rather than a
tidy-looking guess. `lineshapes.base(pka, ph)` **is** Henderson–Hasselbalch, and
has been in the library since before the rewrite. And `processing.ftir`'s
side-chain work already uses pKa values per residue. The pH case is half-built
already; nobody has called it a titration.

So: build `processing.titration` with the model named at the call site, exactly
as ADR-0002 requires for deconvolution methods, and make redox its first caller.
A bespoke `redox.py` would be the third place in this library to write "fit a
parameterised series", and the first two are already here.

### 16.2 The analysis, end to end

The pipeline is one the existing pieces nearly compose already:

1. Load the series; each spectrum carries its potential.
2. Reference and baseline — usually difference spectra against the fully
   oxidised or fully reduced end point.
3. **Decompose** to find how many spectral species are actually changing —
   `processing.multivariate` with `stability()` for the component count. This is
   the step that distinguishes a titration with two cytochromes from one with
   three, and it is the step normally done by eye.
4. **Fit** each component's contribution against potential with the Nernst
   equation, one or more transitions.
5. Report *E*ₘ and *n* per component, with uncertainties.

Step 4 is new; steps 1–3 exist. Step 3 is also where this library has something
to say that the standard tools do not — the usual practice is to pick a
wavelength, plot absorbance against potential, and fit, which silently assumes
the species are spectrally resolved at that wavelength.

### 16.3 ⚠️ A collection cannot yet carry a continuous parameter

`SpectrumCollection.from_files(..., sample_from=...)` attaches a **categorical
sample name** to each spectrum, and `group_by` groups on it. There is nowhere to
put "this spectrum was measured at −120 mV" except `metadata`, by hand, in a loop
after loading.

Every entry in §16.1's table needs the same thing, so this is not a redox
problem. Wanted, and additive: a `parameter_from=` argument that reads a number
from the filename or the file's own metadata, and an accessor on the collection
that returns it as an array aligned with the spectra. Worth doing before the
freeze because it is the natural partner to `to_matrix()` — the analysis wants
`(x, X, parameter)` and can currently only get the first two.

### 16.4 What the dataset needs to contain

Recorded now because the answer costs nothing at acquisition time and is
sometimes impossible to reconstruct afterwards:

- **The potential of every spectrum, and the reference electrode it was measured
  against.** vs SHE and vs Ag/AgCl differ by about 200 mV; a table of numbers
  with no electrode named is not usable, and no analysis can recover it.
- **The temperature.** The Nernst slope is *RT/nF* — 59.2 mV per decade at 25 °C
  and not at 4 °C. Fitting *n* with the wrong temperature returns a wrong
  electron count that looks perfectly reasonable.
- **Both directions.** Oxidative and reductive titrations of the same sample
  should agree; where they do not, the sample was not at equilibrium. This is
  the single best check on the data, and it only exists if both were measured.
- **The mediators used**, and their concentrations. They set what equilibrates,
  and several absorb in the visible.
- **One system with a published *E*ₘ**, for the same reason §15.2 wants a protein
  of known structure: without it, the code can be shown to be self-consistent but
  not correct.

### 16.5 Status

Nothing is built. Working agreement §14.1 applies as it does to OPUS, `.spc` and
CD: the dataset comes first, and it becomes the tutorial. What can be done ahead
of it is the `parameter_from=` gap in §16.3, which is needed regardless of which
of §16.1's five applications arrives first.

---

## 17. Scheduling revision: combined CD + FTIR is post-JOSS (James, 2026-08-02)

**The combined CD/FTIR secondary structure analysis moves after the JOSS
version**, i.e. after 1.0 in November. ADR-0002 stands as the design; what
changes is when the second half of it gets built.

This is the right call, and it makes the November plan credible again. §14.1
argued that the road to 1.0 contains almost no new features and that the way
the date slips is by letting feature work leak into it. CD is a fifth technique,
a units problem (§15.4), and a reference-data licensing question — three
unbounded things three months before a freeze.

### 17.1 What that leaves before November

| | Status |
|---|---|
| FTIR amide I secondary structure | **In scope.** Needs no new technique, no reference data and no new units; the fitter it depends on is built, with pseudo-Voigt and A/d²A weighting both settled by measurement |
| OPUS and `.spc` readers | **In scope, blocked on files** (§15.2) |
| `parameter_from=` on collections | **In scope**, unblocked, and needed by every titration application (§16.3) |
| The freeze work | ADR-0001 §7.3/§7.4, shim removal, `CONTRIBUTING`/`CODE_OF_CONDUCT`/`CITATION.cff`, `paper.md` |
| CD as a technique, deconvolution, the combined estimate | **After 1.0** |
| Redox titration | After the dataset (§16.5) |

### 17.2 One consequence worth deciding explicitly

Earlier, CD was scoped as a full CDSSTR/SELCON equivalent before the tester
emails. If the *combined* analysis is post-JOSS, most of the reason to rush the
CD half before November goes with it — because of what Hoffmann, Jones and
Rodger actually found (ADR-0002 §7.2): combining the techniques improves the
numbers by only about 2 %, and its real value is catching the cases where one
technique alone is badly wrong. That value only exists once both halves are
there.

CD deconvolution on its own is still worth having — people estimate structure
from CD alone every day. But building the expensive half before the payoff, in
the window reserved for not breaking things, is the shape of decision that makes
a November date slip.

**Confirmed (James, 2026-08-02): FTIR only for now, and CD develops on a
separate branch, with timing judged as it goes.**

The branch is a real mechanism rather than a postponement dressed up as one. CD
work can start whenever there is appetite or a dataset, without any of it
touching the release line or the freeze — and if it goes quickly, merging before
November becomes a decision taken on evidence rather than guessed at now. What
it must not acquire is the power to delay 1.0: `main` stays the thing that
ships, and the branch merges when it is ready *and* the freeze work is done, in
that order.

One constraint travels with it. ADR-0002 says both techniques return the same
`Composition` against the same DSSP vocabulary, so the branch cannot invent a
second answer type — and the FTIR estimator on `main` is what defines the
interface it has to satisfy. Building FTIR first is therefore not only
sequencing; it fixes the contract.

### 17.3 Unchanged

The two techniques still produce one `Composition`, and DSSP is still the
baseline vocabulary. Nothing in ADR-0002 is withdrawn — the FTIR estimator built
before November must return the same type, against the same vocabulary, that the
CD one will later fill in. That is the whole point of having settled the design
first: the half built now cannot foreclose the half built later.

---

## 18. Datasets located (2026-08-02)

Recorded so they are not lost between sessions. Working agreement §14.1 asks for
a real example before a feature is rolled out; these are the examples.

### 18.1 Cytochrome P450 reductase — used

`~/Documents/Research/Notebook/2026/03/04`, with the analysis begun in
`~/Documents/Research/FTIR_040326.ipynb`. `CPR_WT.{0..5}.dpt`,
`CPR_Q157.{0..5}.dpt`, `Buffer.{0..3}.dpt`; `.0-.2` wet, `.3-.5` dry.

Already used to correct two defaults (roadmap §17 commit): it is what showed
that fitting amide I alone reports more sheet than helix, and that an
11-point Savitzky-Golay window spans 21 cm⁻¹ on real sampling and fabricates an
aggregation band. No published composition, so it validates the *pipeline*, not
the numbers.

### 18.2 Candice Gomez's protein spectra — the validation set

`~/Documents/Research/Students/Candice Gomez/FTIR_spectra`, 179 files.

| Protein | Files | Why it matters |
|---|---|---|
| BSA | 1, 2, 5, 10 (+ `-evap` series) | Predominantly helical — the case Hoffmann et al. report infrared handling *worst* |
| Lysozyme | 1, 2, 5, 10 (+ `-evap`) | Mixed α/β, and about as well characterised as a protein gets |
| VDAC1 | 6 | A β-barrel: the opposite extreme, and the case CD handles worst |
| AqpZ | `-373`, `-390`, 3 each | Helical membrane protein, in detergent |

Plus the references any of it needs: `PBS`, `Tris` at several pH, `HEPES`,
`H2O`, and the detergents `DDM`, `LDAO`, `SDS`, `bOG`.

**Two independent validations are possible here, and the second needs no
literature at all.**

*Against published compositions.* BSA, lysozyme and VDAC1 all have reference
secondary-structure content in the literature and in the PDB. Those numbers must
be looked up and cited rather than remembered (§14.6), and the comparison is
the figure that turns "the pipeline runs" into "the pipeline is right".

*Against itself.* **The same protein at 1, 2, 5 and 10 must give the same
composition.** Concentration changes the amide I amplitude and nothing else, so
any variation across the series is method error, measured directly, with no
reference values required. That is a stronger internal check than a single
comparison against one published number — and the `-evap` series adds a second
axis, where change is expected and its size is the question.

Together these span helix-rich, mixed and sheet-rich, which is exactly the range
over which a single technique is known to fail.

### 18.3 AqpZ CD temperature scan — for the CD branch

A temperature scan of CD spectra on AqpZ, in the notebook for **2025-07-24**.
James has done the analysis; the notebook has not been found yet.

Doubly relevant. It is CD data, so it is the first real material for the branch
of §17.2. And it is a **series parameterised by temperature** — a thermal melt,
which is §16.1's table with a van 't Hoff model instead of Nernst. So it also
exercises the `parameter_from=` gap of §16.3, which remains the piece of work
that is unblocked, useful to several applications at once, and still not done.
