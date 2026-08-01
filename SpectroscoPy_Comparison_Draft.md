# SpectroscoPy vs. existing Python spectroscopy packages

**Status of this document:** a draft comparison, built from SpectroscoPy's actual
documentation (level-0 pages) against README/PyPI-level descriptions of the
other packages. Claims about SpectroscoPy are marked **[verified]** — read
directly from its docs. Claims about the others are marked **[as documented]**
— not yet confirmed by installing and running them. Before this goes public,
the **[as documented]** items should be checked hands-on, particularly for
SpectroChemPy and RamanSPy, the two closest comparisons (see §5).

---

## 1. Scope comparison

| Package | Techniques | Scope decision |
|---|---|---|
| **SpectroscoPy** | ATR-FTIR, FTIR, UV-Vis, Raman, fluorescence **[verified]** | Deliberately bounded to optical/vibrational spectroscopy common to a single lab bench |
| SpectroChemPy | NMR, IR, Raman, UV-Vis **[as documented]** | Broader — includes NMR, which SpectroscoPy explicitly does not attempt |
| RamanSPy | Raman only **[as documented]** | Narrow and deep — one technique, published methodology |
| rampy | Raman, IR, XAS **[as documented]** | Narrow, function-level, not a unified object model |
| Orange-Spectroscopy | Spectral/hyperspectral, technique-agnostic within Orange **[as documented]** | Broad but delivered as visual-programming nodes, not a scriptable library |
| pyspectra | NIR/FTIR/Raman via `.spc`/`.dx` **[as documented]** | Pandas-wrapper scope, not a dedicated object model |

SpectroscoPy's scope is narrower than SpectroChemPy's (no NMR) but broader than RamanSPy's or rampy's, and — unlike pyspectra or the original pyuvvis — is built around a single object model rather than a pandas convenience wrapper.

---

## 2. Core data model

**SpectroscoPy [verified]:** `Spectrum` and `SpectrumCollection`, with:
- explicit `x_unit`/`y_unit` and a `.to()` conversion method that changes both axis and representation (e.g. transmittance → absorbance, wavenumber → nm)
- baseline-direction (upper vs. lower hull) chosen automatically from the y-unit, rather than requiring the user to know which convention applies
- structured processing history (`describe_history()`), persisted through the native `.spy` format
- collection-level grouping/averaging (`group_by`, `.mean()`) driven by filename convention, not manual indexing

**SpectroChemPy [as documented]:** `NDDataset` object with labeled axes and metadata, project management across multiple datasets simultaneously. Structurally the closest thing to SpectroscoPy's `Spectrum`/`SpectrumCollection` split among all the packages surveyed. What's not yet confirmed: whether `NDDataset` tracks processing provenance the way `.describe_history()` does, or whether unit conversion is as central to its API as SpectroscoPy's `.to()`.

**RamanSPy, rampy, pyspectra [as documented]:** function- or DataFrame-based rather than a dedicated stateful object with its own methods; provenance and unit-safety are the user's responsibility unless done manually.

**Assessment:** SpectroscoPy's provenance-by-default design (every operation both returns a new object *and* logs itself) is not something the other packages document as a first-class feature. This is a genuine, verifiable differentiator — but it's also the exact kind of claim worth re-checking once you've read SpectroChemPy's actual `NDDataset` API docs rather than its front page.

---

## 3. File I/O

**SpectroscoPy [verified]:**

| Format | Extensions |
|---|---|
| Bruker OPUS data point table | `.dpt` |
| JCAMP-DX (incl. compound files) | `.jdx`, `.dx`, `.jcamp` |
| Delimited text (header detected) | `.csv`, `.tsv`, `.txt` |
| Wide/paired tables (many spectra per file) | any |
| Native (keeps units + history) | `.spy` |

One honest gap worth naming directly: SpectroscoPy currently reads OPUS via the
**exported `.dpt` text table**, not the raw binary OPUS file directly. SpectroChemPy's
ecosystem includes a dedicated `brukeropusreader` for the binary format
**[as documented]** — if a user's workflow exports to `.dpt` anyway this is a
non-issue, but if someone only has raw `.0` OPUS files, SpectroscoPy doesn't
yet read them without an export step. Worth stating plainly rather than
implying full OPUS support.

**Others [as documented]:** pyspectra reads `.spc`/`.dx`; SpectroChemPy's plugin
ecosystem appears to cover OMNIC and OPUS binary formats specifically. None of
the surveyed packages document a JCAMP-DX-with-compound-files reader or a
generic wide/paired-table auto-loader the way SpectroscoPy does — those look
like a genuine edge for handling "one file, many spectra" instrument exports.

---

## 4. Processing, analysis, and reproducibility

| Capability | SpectroscoPy | SpectroChemPy | RamanSPy | rampy |
|---|---|---|---|---|
| Baseline correction | rubberband, ALS **[verified]** | documented, methods TBD **[as documented]** | yes **[as documented]** | yes **[as documented]** |
| Peak finding | second-derivative default, relative prominence **[verified]** | unconfirmed | yes **[as documented]** | via scipy, manual |
| Peak fitting/decomposition | `lineshapes` module exists; general Gaussian/Lorentzian fitting present, but domain wrappers (e.g. amide I secondary-structure decomposition) explicitly **planned, not built** **[verified — self-reported gap]** | unconfirmed | yes, ML-integrated **[as documented]** | via lmfit, manual |
| Multivariate (PCA/NMF) | yes, optional extra (`decompose(...).contributions()`, `stability(...)` for component-count validation) **[verified]** | yes (SVD, PCA, MCR_ALS, EFA, PLS) **[as documented]** | via sklearn integration **[as documented]** | via sklearn integration |
| Reference subtraction with recorded scale factor | yes **[verified]** | unconfirmed | unconfirmed | unconfirmed |
| Reproducibility / provenance | `.describe_history()`, persisted in `.spy` **[verified]** | unconfirmed | unconfirmed | unconfirmed |

SpectroscoPy's own docs list several **planned, not built** analysis features worth
naming honestly in any public comparison, since they're real gaps against
mature tools:
- No `concentration()` / Beer–Lambert helper (UV-Vis) — read the peak manually today
- No A260/A280 two-component deconvolution
- No amide I secondary-structure decomposition — flagged in SpectroscoPy's own
  docs as its **highest-value gap**
- No spectral-library/database matching for compound identification
- No inner-filter correction for fluorescence, no absorptance-based
  energy-transfer ratio calculation

This is a case where the source of the gap list matters: these aren't guessed
weaknesses, they're the project's own stated backlog, generated from what
actual users asked for and hadn't gotten yet. That's a stronger, more credible
basis for "here's what's missing" than a competitor guessing at your gaps —
worth explicitly saying so in the published version, since it also implicitly
demonstrates the feedback loop described in §6 already works.

---

## 5. GUI

| Package | GUI status |
|---|---|
| SpectroscoPy | **Planned, not built [verified].** Docs explicitly tell code-averse users this isn't for them yet |
| Orange-Spectroscopy | **Shipped** — but as Orange3 visual-programming nodes, not a coded-library-plus-app pairing **[as documented]** |
| SpectroChemPy | Separate optional GUI package exists **[as documented]**, maturity unconfirmed |
| RamanSPy, rampy, pyspectra | None |

This is the most significant honest gap in the comparison: Orange-Spectroscopy
already has a working GUI for non-coders, and it's the package SpectroscoPy's
own "which kind of user are you?" section in Getting Started implicitly
concedes to right now ("a graphical application is planned; until it exists...").
The differentiator you're aiming for — one library, scriptable *and* GUI-driven
from the same API — is a real gap in the market against Orange's node-based
paradigm, but it's a gap SpectroscoPy hasn't closed yet either. This should be
framed as roadmap, not as a present advantage.

---

## 6. Documentation approach and user feedback loop

**SpectroscoPy [verified]:** Three-tier structure — a task-oriented "what do you
want to know?" index (mapped to real questions people have actually asked, with
an honest works/partly/planned status column), a linear "getting started" that
explicitly does not assume the reader writes Python, and a conventional API
reference kept deliberately last ("you should not need this to get work done").
The feedback address is repeated on every page and explicitly framed as
accepting *both* missing-feature reports and bugs, with an explicit statement
that the user doesn't need to know how a feature would be implemented to
request it.

**Others [as documented]:** conventional docs sites (Sphinx/ReadTheDocs-style
API references, tutorial galleries). None of the blurbs/README content
surveyed describe an equivalent task-first entry point or an explicit
status-labeled gap list as a documentation artifact in its own right.

This is arguably SpectroscoPy's clearest, most defensible differentiator right
now — not a feature of the library, but a design decision about how the
project relates to non-expert users, which is exactly the pain point that
motivated the project in the first place (vendor tools with "hard to learn
interfaces").

---

## 7. Maturity and credibility signals

| Package | Signal |
|---|---|
| SpectroscoPy | Pre-1.0, API explicitly not yet settled; docs state it is "used for real work" **[verified, self-reported]** — no external validation yet (no citations, no independent adopters documented) |
| SpectroChemPy | Also **explicitly pre-1.0 / actively changing**, despite being longer-established — "current design may undergo major changes" **[as documented]** |
| RamanSPy | Published in a peer-reviewed journal (Analytical Chemistry, 2024) **[as documented]** — the strongest external credibility signal among the surveyed packages |
| rampy, pyspectra, Orange-Spectroscopy | Established, in active-enough use, no peer-reviewed citation found in this search |

Worth noting directly in any public version: SpectroscoPy is not the only
pre-1.0, actively-changing package in this space — SpectroChemPy carries the
same warning despite more maturity elsewhere — so "pre-1.0" alone isn't a
weakness unique to SpectroscoPy. What *is* a real gap is external validation:
RamanSPy has a peer-reviewed paper behind it; SpectroscoPy currently has none.

---

## 8. Summary: honest positioning

**Where SpectroscoPy is currently ahead (verified, not just asserted):**
- Provenance-by-default core design (history tracked and persisted automatically)
- Unit-safety baked into the core object rather than left to the user
- A documentation model organized around real user questions with an honest, status-labeled gap list
- A "wide/paired table" and compound-JCAMP loader that none of the surveyed competitors document

**Where SpectroscoPy is currently behind, by its own admission:**
- No GUI (Orange-Spectroscopy has shipped one)
- No spectral-library matching, no `concentration()` helper, no secondary-structure decomposition — all flagged as its own top gaps
- No raw binary OPUS reader (only the exported `.dpt` table)
- No peer-reviewed validation (RamanSPy has this)

**Recommended framing for a public-facing version:** lead with the provenance
and unit-safety design, and the documentation philosophy — these are real,
checkable differentiators today. Present the GUI and analysis gaps as an
explicit, dated roadmap rather than omitting them; readers who check GitHub or
try the tool themselves will find the gaps regardless, and stating them
first is what makes the "ahead" claims credible.

---

## 9. Before publishing this

1. Install SpectroChemPy and RamanSPy and run one of your own spectra through
   each — confirms or corrects every **[as documented]** line above, especially
   §2 and §4.
2. Read SpectroChemPy's actual `NDDataset` API docs (not just its README) to
   check the provenance/history claim in §2 head-on — if `NDDataset` already
   does something equivalent, that's the single fact most likely to change
   this document's conclusions.
3. Decide where this document lives: a `COMPARISON.md` or docs page, not the
   README — keeps the README focused on your own value proposition rather than
   reading as a takedown of adjacent projects.

---

## 10. Notes added 2026-08-02

### 10.1 Where this document lives — decided: repository root, for now

It stays at the root beside `SpectroscoPy_Development_Roadmap.md` and
`SpectroscoPy_Codebase_Review.md`, because it is the same kind of object as
those two: an internal working document that argues, records what is not yet
checked, and is written to be revised. That is the right home for it *because*
of §9's caveat — most of the competitor claims are still `[as documented]`, and
a public page making unverified claims about neighbouring projects is a
credibility risk that lands on the project making them.

What should eventually be public is a **different, shorter document derived from
this one**: a "How this compares" page under `docs/`, written after the §9
verification, stating only what has been checked hands-on and dropping the
positioning advice entirely. Two documents, not one moved document — the honest
internal version needs to keep saying things ("Orange has already shipped the
GUI we are promising") that do not belong on a project's own website.

There is a third use for it, and it is the one that pays: **§1–§5 are the
related-work section of the roadmap §13 paper**, already half-written. A
comparison table that a reviewer accepts is one where the competitor claims were
verified by running them, which is exactly what §9 asks for. So the verification
work is not overhead on this document — it is paper work done early.

### 10.2 Release data, checked rather than assumed

§7 rated maturity from README-level impressions. The actual PyPI record:

| Package | Latest | Last release |
|---|---|---|
| spectrochempy | 0.11.0 | **2026-07-09** |
| orange-spectroscopy | 0.9.3 | **2026-07-27** |
| rampy | 0.6.4 | 2026-03-30 |
| ramanspy | 0.2.10 | 2024-06-01 |
| brukeropusreader | 1.3.4 | **2019-09-16** |
| spc-spectra | 0.4.0 | **2018-05-03** |

Two things follow. First, SpectroChemPy and Orange-Spectroscopy are **actively
released right now**, within the last month — §7's framing of SpectroChemPy as
"longer-established but also pre-1.0" understates it, and any claim of movement
relative to them has a moving target. RamanSPy, by contrast, has not released
since shortly after its 2024 paper: the peer-reviewed citation is real, the
development is not obviously continuing. That is worth knowing before treating
publication as the finish line — roadmap §13's paper is a milestone in the
project, not a destination for it.

Second, and directly actionable: **both libraries that would supply the formats
roadmap §13.3 wants are abandoned** — `brukeropusreader` since 2019, and
`spc-spectra` since 2018. Neither is a dependency worth taking on for a format
the project needs to be reliable about. OPUS binary and Galactic `.spc` are
read-them-ourselves jobs, which is also consistent with the §5.6 dependency
policy. The `.dpt` reader is the precedent: the format turned out not to be what
it looked like (review §7.1), and owning the parser is what made that
discoverable.

### 10.3 One more differentiator, from the paper

Roadmap §13 adds a claim this document does not yet make: **library lookup
coupled to decomposition**. §3 and §8 note that spectral-library matching is a
gap against the competition. The plan is not to close it by matching that
feature, but to pair it with NMF — decompose a mixture, then identify the
recovered components against references, rather than matching a whole spectrum
against a library and hoping it is pure. If that works as intended, the row in
§4's capability table is not "library matching: now yes too" but a capability
none of the surveyed packages document.
