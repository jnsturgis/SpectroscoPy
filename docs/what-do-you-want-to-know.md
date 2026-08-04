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

**works** means the code is tested and documented. It does not mean the method
has been checked against measured samples whose answer was already known —
those are different claims, and a method can pass the first and fail the
second. Where the distinction matters, the row says so.
:::

## Identifying and comparing materials

| I want to know… | I need | How | Status |
|---|---|---|---|
| What is this substance? | ATR-FTIR or Raman | Baseline, normalise, pick peaks, compare with a reference spectrum or an assignment table | **works** — [ATR-FTIR](tutorials/atr-ftir.md), [Raman](tutorials/raman.md) |
| Are these two samples different, and where? | any technique, both samples | Put on a common axis, normalise to a shared band, subtract | **works** — [ATR-FTIR](tutorials/atr-ftir.md) |
| Which bands differ between many samples? | a series | NMF or PCA across the set; read the component that separates them | **works** — see `processing.multivariate` in the [guide](guide/processing.md) |
| Is this the same batch as last time? | archived reference spectra | Load both, difference, quantify the residual | **partly** — a reference library exists (`library.Library`); automatic matching does not |
| What is this contaminant? | the spectrum, a library | Peak list, search against a database | **partly** — you can build a library and fit against it ([UV-Vis components](guide/uv-vis-components.md)); there is no peak-list database search, and no reference spectra ship with the package |

## Composition and quantity

| I want to know… | I need | How | Status |
|---|---|---|---|
| How much of X is in solution? | UV-Vis, an extinction coefficient | Beer–Lambert, with the path length stated | **works** — `library.concentration_from_absorbance`; published coefficients for dsDNA, ssDNA, RNA and protein, or `protein_epsilon_280` from the sequence |
| What is the protein : nucleic acid ratio? | UV-Vis 240–340 nm | Fit both reference spectra at once, not one ratio | **works** — `unmix.nucleic_acid_and_protein`; the [guide](guide/uv-vis-components.md) shows why A260/A280 gives the same number for a clean and a contaminated sample |
| How much of each component, when I know what the components are? | reference spectra for each | Non-negative least squares against the references | **works** — `unmix`; the residual, not R², is the check that nothing is missing |
| How much of each component, when I do not? | a series spanning the composition | NMF, then component area fractions | **works** — `decompose(...).contributions()` |
| Is the number of components I chose defensible? | the same series | Bootstrap stability, not explained variance | **works** — `stability(...)`; see the [guide](guide/processing.md) |
| What is the extinction coefficient of my compound? | spectra at several known concentrations | Least squares through the origin, per wavelength | **works** — `library.from_series` returns an ε spectrum *and* its standard error at every wavelength |
| Which wavelengths best separate my components? | the reference spectra | Maximise the conditioning of the reference matrix | **works** — `unmix.best_wavelengths` |
| What is the lipid : protein ratio in this membrane? | ATR-FTIR | Band areas: C–H stretch vs amide I | **partly** — crop and integrate yourself |
| How much did this band change? | before and after | Normalise to an invariant band, subtract, integrate | **partly** |

:::{admonition} Two things to get right before believing any concentration
:class: warning

**State the path length** if it is not 1 cm. Nothing in the data records it, so
a 1 mm cuvette left unmentioned makes every concentration on this page ten
times too small, and nothing about the result looks wrong.

**Correct scattering first** if the sample is at all turbid — see *Structure
and state*, below. Scattering is not one of your
components, so a fit will spread it across whichever references happen to slope
the same way. In the guide's worked example that inflates the nucleic acid
eightfold, with a perfectly good-looking fit throughout.

The unmixing and extinction-coefficient routines are currently checked against
**synthetic** mixtures — built from invented references so the recovered
amounts can be compared with a known truth. That tests the arithmetic. It does
not yet test whether the methods survive a real instrument, and measurements to
close that gap are planned.
:::

## Structure and state

| I want to know… | I need | How | Status |
|---|---|---|---|
| What is this protein's secondary structure? | ATR-FTIR, amide I (1600–1700 cm⁻¹) | Water subtraction, second-derivative band positions, Gaussian decomposition | **partly** — it runs, and there is a [guide](guide/secondary-structure.md), but the same protein at four concentrations comes back spread over ±20 percentage points. Not yet a measurement |
| Did my protein change conformation? | before and after | Difference spectrum in amide I, second derivative | **partly** — the difference is trustworthy in a way the absolute composition above is not |
| Is my sample aggregated or scattering? | UV-Vis, a window where nothing absorbs | Fit a basis of power laws where nothing absorbs, extrapolate, subtract | **works** — `processing.scattering`; or `from_references` to fit measured blanks, which is better where you have them |
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
| Can this even read the file my instrument wrote? | `spc.Spectrum.read(path)` — the format is detected from the file, not the extension. Bruker OPUS and Galactic `.spc` binaries, JCAMP-DX, `.dpt`, and text with either decimal separator. See [reading files](guide/reading-files.md) | **works** |
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
