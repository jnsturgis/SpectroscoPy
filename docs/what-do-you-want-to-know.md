# What do you want to know?

The tutorials answer *"how do I turn this data into information?"*. This page
answers the question you probably arrived with: **"I want to know X — what
should I do?"**

Find your question, check what it needs, follow the link.

:::{admonition} Status column
:class: note

**works** — supported now, with a worked example.
**partly** — the pieces exist; you assemble them yourself.
**planned** — not built. Listed because the list is a statement of what this
library is for, and the gaps are the roadmap.
:::

## Identifying and comparing materials

| I want to know… | I need | How | Status |
|---|---|---|---|
| What is this substance? | ATR-FTIR or Raman | Baseline, normalise, pick peaks, compare with a reference spectrum or an assignment table | **works** — [ATR-FTIR](tutorials/atr-ftir.md), [Raman](tutorials/raman.md) |
| Are these two samples different, and where? | any technique, both samples | Put on a common axis, normalise to a shared band, subtract | **works** — [ATR-FTIR](tutorials/atr-ftir.md) |
| Which bands differ between many samples? | a series | NMF or PCA across the set; read the component that separates them | **works** — see `processing.multivariate` in the [guide](guide/processing.md) |
| Is this the same batch as last time? | archived reference spectra | Load both, difference, quantify the residual | **partly** — no spectral-library matching |
| What is this contaminant? | the spectrum, a library | Peak list, search against a database | **planned** — needs a reference library format |

## Composition and quantity

| I want to know… | I need | How | Status |
|---|---|---|---|
| How much of X is in solution? | UV-Vis, an extinction coefficient | Baseline, read the peak, Beer–Lambert | **partly** — read the peak yourself; no `concentration()` yet |
| What is the protein : nucleic acid ratio? | UV-Vis 250–300 nm | A260/A280 two-component deconvolution | **planned** — the maths is three lines and belongs in the library |
| How much of each component in a mixture? | a series spanning the composition | NMF, then component area fractions | **works** — `decompose(...).contributions()` |
| Is the number of components I chose defensible? | the same series | Bootstrap stability, not explained variance | **works** — `stability(...)`; see the [guide](guide/processing.md) |
| What is the lipid : protein ratio in this membrane? | ATR-FTIR | Band areas: C–H stretch vs amide I | **partly** — crop and integrate yourself |
| How much did this band change? | before and after | Normalise to an invariant band, subtract, integrate | **partly** |

## Structure and state

| I want to know… | I need | How | Status |
|---|---|---|---|
| What is this protein's secondary structure? | ATR-FTIR, amide I (1600–1700 cm⁻¹) | Water subtraction, second-derivative band positions, Gaussian decomposition | **planned** — the highest-value gap; the notebooks do it by hand |
| Did my protein change conformation? | before and after | Difference spectrum in amide I, second derivative | **partly** |
| Is my sample aggregated or scattering? | UV-Vis | Rising λ⁻⁴ background towards the blue | **partly** — visible; not quantified |
| Is this cofactor oxidised or reduced? | UV-Vis, both states | Difference spectrum | **partly** |
| Which side chains are ionised? | ATR-FTIR, a pH series | Compare with calculated side-chain spectra | **works** — `processing.ftir.ftir_sidechain` |

## Photophysics

| I want to know… | I need | How | Status |
|---|---|---|---|
| What does this fluorophore emit? | fluorimeter | Emission scan | **works** — [fluorescence](tutorials/fluorescence.md) |
| Which absorbing species feeds the emission? | excitation-emission series | EEM map, scatter ridges masked | **works** — [fluorescence](tutorials/fluorescence.md) |
| How efficient is energy transfer? | absorption **and** excitation of the same sample | Convert absorption to absorptance (1−T), normalise both at a band where transfer is complete, divide | **planned** — `absorptance` exists; the ratio does not. See review §13.1 |
| Is my fluorescence being distorted by absorption? | absorbance of the same sample | Inner-filter correction | **planned** — keep A below ~0.1 meanwhile |

## Housekeeping

| I want to know… | How | Status |
|---|---|---|
| What did I do to get this figure? | `spectrum.describe_history()` | **works** |
| Can I reopen this in six months and still know? | Save as `.spy` | **works** |
| Which of my replicates went wrong? | `viz.grid(spectra, key='sample')` before averaging | **works** |
| Can my collaborator read my file? | Export `.csv`, or send the `.spy` | **works** |

## If your question is not here

That is useful information. The gaps above are the plan; a question that is not
listed at all is one nobody has told us about.

So if your question is missing, {{ feedback_question_link }} — what you want to
know, and what instrument you have. A question nobody reports stays missing,
and the entries on this page came from people who said what they were trying to
find out. You do not need to know how it would be implemented, or whether it is
reasonable to ask.

The same address takes bugs and feature requests; the links are at the foot of
every page.
