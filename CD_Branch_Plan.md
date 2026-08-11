# Circular dichroism — branch plan

Branch `circular-dichroism`, opened 2026-08-05. This is the branch roadmap
§17.2 called for: CD work that can start whenever there is appetite or a
dataset, **without any of it touching the release line**.

## The rules this branch inherits

Not re-argued here; these are settled and the plan is built on them.

- **`main` is what ships.** This branch merges when it is ready *and* the
  freeze work is done, in that order. It never acquires the power to delay
  1.0 (§17.2).
- **One `Composition`, one vocabulary.** Both techniques return the same type
  against the same DSSP categories, and the FTIR estimator on `main` defines
  the interface this must satisfy. The branch cannot invent a second answer
  type (ADR-0002, §17.3).
- **No reference data ships until its redistribution terms are checked**
  (ADR-0002 §9).
- **Reuse or improve, do not reimplement blind.** Hoffmann, Jones and Rodger
  (*QRB Discovery* **6**, 2025, [10.1017/qrd.2025.4](https://doi.org/10.1017/qrd.2025.4))
  published a Python SELCON3. Comparing against theirs beats writing a third
  one and assuming it agrees.
- **A real file before a reader.** The `.dpt` precedent stands: that reader was
  written against real files and the format turned out not to be what it looked
  like.

---

## ⚠️ Two decisions that must be taken on `main` before November

**This is the most urgent thing in this document, and the only part with a
deadline.** Both are cases where CD, though post-1.0, needs something from an
interface that freezes in November — exactly what §17.3 means by *the half
built now cannot foreclose the half built later*. Neither is expensive; both
become breaking changes if missed.

### D1. `band_direction` needs a fourth value, and CD is the reason

`units.band_direction()` documents exactly three returns — `'up'`, `'down'`,
`'unknown'`. **A CD spectrum is bipolar**: an α-helix gives negative bands at
208 and 222 nm and a positive one at 193 nm, in the same spectrum, and both
signs are real signal rather than one being a baseline artefact.

`'unknown'` is the wrong answer for it. Unknown means *the quantity does not
fix a direction*; for CD the quantity fixes it precisely, and the answer is
*both*. Peak finding on CD must return both maxima and minima, which is a case
the current model has no way to express — `detect_peaks` chooses one via
`troughs=`.

Adding a `'both'` return after 1.0 is a behaviour change for anyone who wrote
`if direction == 'up' ... else ...`. Adding it before costs a line.

**Recommendation:** add `'both'` to the documented contract on `main` before
the freeze, even with no CD unit yet using it, and have `detect_peaks` handle
it by returning both signs. ADR-0003 needs the same amendment.

### D2. The metadata keys for concentration, path length and residue count

`.spy` serialises `metadata`, so **metadata key names are part of the frozen
format** whether or not anyone has said so. CD's central conversion needs three
sample quantities that no spectrum currently carries in an agreed place:
concentration, path length, and residue count (or mean residue weight).

Part of this already exists by accident: `unmix()` and `correct_scattering()
`already read `metadata['path_length']`, established in August without being
written down as a schema. CD adds two more, and a redox titration will want
temperature and electrode.

**Recommendation:** agree the small set of sample-metadata keys on `main`
before November — names, units and what "absent" means — and document them in
ADR-0001. Deciding it now costs an afternoon; deciding it after 1.0 means
either a second set of names or a format migration.

---

## Work packages

### WP1 — CD as a technique ✅ done 2026-08-08

- `KNOWNSPECTYPES['CD']`: x Wavelength/nm, y Ellipticity/mdeg. **Not** in
  `REVERSED_AXIS_TECHNIQUES`.
- New y units, all of which are *quantities in their own right* rather than
  spellings of each other: `mdeg`, `deg`, mean residue ellipticity
  (deg·cm²·dmol⁻¹), and Δε (M⁻¹cm⁻¹ per residue).
- Register their band direction — see **D1**; this is where that bites first.
- The spectrum must know it is CD before any of the analysis makes sense, so
  this package gates everything else. It is also genuinely small.

### WP2 — the units problem (§15.4), the one non-additive piece

`units` converts with a table where every conversion is a function of `y`
alone. **CD's is not.** Getting from millidegrees to mean residue ellipticity
or to Δε per residue needs concentration, path length and residue count, and
`to()` has nowhere to get them.

Settled recommendation (§15.4): **a separate explicit method now, `pint` after
1.0** — it adds a name rather than changing the meaning of an existing one, and
a wrong answer here would have been frozen in November.

```
spectrum.to_mean_residue_ellipticity(concentration=..., path_length=...,
                                     residues=...)
```

**Before writing the code, write the conversion down and check it against a
published worked example.** Both conventions in circulation — MRE via mean
residue weight, and Δε per residue — and the factor-of-ten traps between mdeg,
deg, molar and per-residue are exactly the kind of thing that produces a
plausible wrong number. This is the `SPC_Format_Notes.md` pattern, which worked:
write the specification first so the code has something to be checked against.

Defaults are a trap here. If concentration or residue count is missing, this
must **raise**, not assume 1. A silently assumed concentration gives an
ellipticity wrong by whatever factor, with nothing on the face of it.

### WP3 — readers

CD instruments: JASCO (`.jws` binary; also text export), Aviv, Applied
Photophysics Chirascan, OLIS.

Order of work, following the working agreement: **text exports first**, binary
only with real files *and* documentation in hand. The proprietary-format rule
stands — no writer for a format that cannot be verified.

Blocked on knowing which instrument produced the AqpZ scan (§18.3). Finding
that notebook is the cheapest unblocking action available.

### WP4 — deconvolution 🔨 machinery done 2026-08-08

Per ADR-0002 §7.2, two kinds of standard, and the method name must say which:

1. **Basis spectra of pure structures** — the composition *is* the coefficient
   vector. Simple, and only as good as the basis.
2. **Reference proteins of known structure** — SELCON3, CDSSTR, CONTIN-LL. Fit
   as a combination of reference proteins, then take the same combination of
   their DSSP compositions.

Much of the machinery already exists on `main` and should be reused rather than
rewritten: `processing.unmix` is non-negative least squares against a reference
set with a residual diagnostic and condition-number-based wavelength selection,
which is structurally the same operation. The CD-specific parts are the
reference sets, the DSSP mapping, and the selection rules that distinguish
SELCON from CDSSTR.

**Check the licence of the published SELCON3 before reading its code.** This
project has already had one near miss — a GPL-3.0 reference implementation that
could not be adapted into an MPL-2.0 codebase, resolved by writing the
specification down independently. Same discipline applies.

Reference sets (SP175, SMP180 and similar) are **data with terms**, not
free-floating numbers. Nothing ships until checked; a helper documents where to
get them.

`theta-222` — helicity from ellipticity at 222 nm — is worth having and is
**not a decomposition**. It fills one category and leaves the rest `None`, and
must never be presented as an estimate of the whole composition.

### WP5 — validation, on §19's terms

**This is the package that decides whether any of it ships.** §19 measured the
FTIR estimator against a concentration series and found ±20 percentage points
at R² = 0.999. The lesson generalises: *goodness of fit does not validate a
composition.* CD gets the same treatment before any number is presented as
usable.

Three checks, in increasing strength:

1. **Invariance.** The same protein at different concentrations and path
   lengths must give the same composition. Computed from data alone, nothing to
   tune towards — this is the §19.3 objective function.
2. **Against DSSP.** CD reference proteins have solved structures; that is why
   they are reference proteins. A composition can be compared with the real
   answer, which FTIR never had.
3. **Against the published implementation.** Same input, same reference set,
   compare with Hoffmann *et al.*'s SELCON3. Disagreement is informative in
   both directions.

**Publish the spread beside any composition, always.** Same rule as `main`.

### WP6 — the AqpZ temperature melt (§18.3) 🔨 the fit is built

A CD temperature scan on AqpZ, notebook 2025-07-24, analysis already done by
James; the notebook has not been found.

This is the branch's best early target because it is **real data that exercises
new code immediately**: it is a series parameterised by temperature, which
`from_files(parameter_from=...)` and `sorted_by_parameter()` now handle
directly, as of 2026-08-04. Loading the melt, ordering it, and plotting it is
achievable in an afternoon and needs none of WP2 or WP4.

Then a two-state van 't Hoff fit — §16.1's table with van 't Hoff instead of
Nernst, and the same missing `processing.titration`.

### WP7 — the combined estimate (needs both halves)

`Composition.compare()`. Per ADR-0002 §7.2, combining IR and CD improves the
numbers by only ~2 %; **the real value is catching the cases where one
technique alone is badly wrong** — IR on high-helix proteins, CD on high-sheet.
So `compare` should report *where* two estimates diverge in terms a user can
act on, and the guide should say plainly that agreement is weak evidence both
are right, while disagreement is strong evidence one is wrong.

---

## Sequence

| Order | What | Depends on | Blocked? |
|---|---|---|---|
| 1 | **D1 and D2 raised on `main`** | — | No. **Deadline: November** |
| 2 | WP1, CD as a technique | D1 | No |
| 3 | WP6a, load and plot the AqpZ melt | WP1, WP3 text reader | On finding the notebook |
| 4 | WP2, the units conversion | WP1, D2 | No |
| 5 | WP3, readers proper | Real files | On instrument identification |
| 6 | WP4, deconvolution | WP2, reference sets | On licence checks |
| 7 | WP5, validation | WP4 | On reference proteins |
| 8 | WP6b, van 't Hoff | `processing.titration` | Unbuilt |
| 9 | WP7, combined | FTIR half working (§20) | On §19's diagnosis |

Items 1 and 2 are unblocked today. Item 3 is unblocked by finding one notebook.

## Merge criteria

Not "it works" — this branch merges when:

1. The freeze work on `main` is done and 1.0 is out.
2. WP5 has *measured* a spread, and it is small enough that a composition means
   something. If it is not, the code can still merge with the spread published
   and the numbers labelled as not yet usable — the FTIR precedent — but that
   has to be a deliberate statement, not silence.
3. No reference data is shipped whose terms have not been checked.
4. `Composition` is unchanged, or changed on `main` first with the FTIR
   estimator updated to match.

## What this branch must not do

- Touch the release line, or add a dependency to `main`.
- Change `Composition`, the DSSP vocabulary, or any frozen signature.
- Ship reference spectra of unchecked provenance.
- Become a reason 1.0 slips. If it starts to look like one, it stops.


---

## Progress, 2026-08-08

**WP1 done.** `KNOWNSPECTYPES['CD']` — Wavelength/nm, Ellipticity/mdeg, not
reversed. Mean residue ellipticity and delta-epsilon are registered as units
and as bipolar, though nothing converts to them yet (that is WP2). D1 does its
job unprompted: a helix CD spectrum peak-picks both ways and finds +193
alongside −208/−222.

**WP4, the machinery.** `structure.from_cd(spectrum, method, basis=...)`,
returning the same `Composition` as `from_ftir` per ADR-0002. Both kinds of
standard are supported and the method must be named, because they answer
differently: `'basis-spectra'` (coefficients *are* the composition) and
`'reference-proteins'` (fit proteins, then mix their known structures).

**The fit is scale-free**, and that turned out to matter more than expected.
The first version imposed "fractions sum to one" as a weighted constraint row;
it made the answer depend on how loud the spectrum was, because the constraint
fought the unknown amplitude instead of describing the shape. Fitting freely
and normalising afterwards imposes the same constraint for nothing and leaves
the result independent of concentration — which is what lets a protein of
unknown concentration still give a composition, and removes WP2 from the
critical path for the commonest question.

Still no reference data ships, by ADR-0002 §9. The arithmetic is tested against
synthetic mixtures of a known basis and recovers 0.550/0.200/0.250 from
0.55/0.20/0.25. That proves the arithmetic and nothing about any real basis;
`quality['rmsd_relative']` is what says whether a basis could describe your
spectrum at all.

**WP6, the fit.** `processing.melting` — deliberately technique-agnostic, since
a melt is the same analysis whether it is CD at 222 nm, tryptophan
fluorescence, or absorbance. Two-state van 't Hoff with **sloping baselines on
both states**, six parameters. Recovers Tm 51.99 ± 0.01 against a true 52.0 and
dH 250.3 ± 0.4 against 250.

### The part worth keeping: the fit is not the evidence

§19's lesson applied before it could bite. A two-state curve fits three-state
data at one wavelength with a residual of 1.4% of range and reports a single
Tm of 61.7 °C — against true transitions at **42 and 62 °C**. It finds the
second; the first leaves no trace in Tm at all.

So two checks that use information the fit does not, both across the whole
spectrum:

| | two-state | three-state |
|---|---|---|
| fit residual / range | 0.002 | **0.014** — unremarkable |
| isodichroic tightness | 0.022 ✓ | 0.050 ✗ |
| third singular value / noise | 1.6 ✓ | **20** ✗ |

`isodichroic_point` — if a sample really passes between two states, every
spectrum is a mixture of the same two and they all cross at one wavelength.
`two_state_rank` — a two-state series is a rank-2 matrix, so the third singular
value should sit at the noise floor.

Getting the second right needed a correction: the first version mean-centred
the matrix, which removes one rank, so it was testing for a *fourth* species
and reported "two-state" on the three-state series. Uncentred, and compared
against the noise floor estimated from the trailing singular values, it
separates 1.6 from 20.

The fitted dH also carries a signature — 39 kJ/mol against the 220 that built
each step, the classical sign of a non-two-state transition. It is worth
checking, but it needs an expectation to compare against, which an unknown
protein does not come with. The two structural checks need none.

### Next

- **WP2**, the units conversion, is now needed only for absolute quantities:
  `theta-222` helicity and anything comparing amplitudes between samples. Write
  the conversion down and check it against a published worked example before
  coding it — both MRE conventions and the mdeg/deg and molar/per-residue
  factors are exactly where a plausible wrong number comes from.
- **WP3**, a reader, still blocked on knowing which instrument wrote the AqpZ
  scan.
- **WP5** unchanged: none of this is validated on a real protein.


---

## The AqpZ data, found and looked at (2026-08-09)

The notebook §18.3 was waiting for is
`Documents/Research/Notebook/2025/AqpZ_Lipid.ipynb`, and the data is in
`2025/07/24/`: **JASCO J-815**, 280→180 nm, 1001 points, CD in mdeg with the
**HT voltage recorded as a second channel**. Two temperature series, 30–90 °C:
`3D_scan1.csv` at 0.5 mg/mL (20 µM) and `3D_scan2.csv` at 0.1 mg/mL (4 µM).
Sample is **AqpZ-W14A**, not wild type. WP3 is answered: the instrument is a
JASCO J-815 and the export is text, so no binary reader is needed to proceed.

### The HT channel decides what analysis is possible

The J-815 records the photomultiplier voltage alongside the signal, and above
about 600 V the detector is starved and the CD is not measurement any more.
Applying that:

| | usable to (HT < 600 V) | HT at 222 nm |
|---|---|---|
| 0.5 mg/mL | **209.7 nm** | 404 V |
| 0.1 mg/mL | **196.2 nm** | 342 V |

**The concentrated sample is the unusable one.** Below 210 nm it is absorbing
its own measurement; by 195 nm the HT is pinned at 1023.7 V, the instrument
ceiling, and the apparent +19 mdeg "band" at 193 nm is the detector flailing,
not protein.

That settles which spectrum a shape analysis can use, and it is not the
obvious one. A basis-spectra fit needs roughly 190–240 nm, because what
separates helix from sheet from coil lives below 210 — which is why the
published reference sets go to 175 nm. At 210 nm there is one negative band
and no shape left to fit. **The dilute scan, at 196 nm, is the only candidate.**

### What the structure says, so there is something to check against

1RC2 chain A: 231 residues, all resolved, **178 in helices** by the
crystallographer's HELIX records (no sheet at all — AqpZ is an all-helical
channel) = **77.1 %**.

The construct is longer than the crystal. A 12-residue N-terminal extension is
12 more residues the spectrophotometer sees and the crystal does not, so the
expected CD helix fraction is diluted:

| extension | chain | expected helix |
|---|---|---|
| none | 231 | 0.771 |
| 12 aa | 243 | 0.733 |
| **23 aa (confirmed)** | **254** | **0.701** |

**Confirmed 2026-08-09: 23 residues** — the full `MGHHHHHHHHHHSSGHIEGRHEF`
tag, uncut. (The 12 was a different construct. The tag is removed with TEV,
not Factor Xa; the IEGR in this sequence is not the site they use.) So the
expectation the CD must be checked against is **70.1 % helix, no sheet**.

### theta-222, and a path length recovered rather than recorded

The path length was not written down, which normally sinks any absolute
quantity. Here there are only a few plausible cuvettes, and the answer is
extremely sensitive to which:

| path | 20 µM | 4 µM |
|---|---|---|
| 0.2 mm | 172 % helix | 188 % |
| **0.5 mm** | **68.7 %** | **75.3 %** |
| 1.0 mm | 34.3 % | 37.7 % |

0.2 mm is impossible and 1.0 mm gives half the helix the crystal has.
**At 0.5 mm the two independent dilutions give 68.7 % and 75.3 %, bracketing
the crystal's 73.3 %.** That is a consistency argument rather than a
measurement — it assumes the structure to infer the cuvette — but the
agreement across a five-fold dilution is not something a wrong path length
would produce. **Worth confirming against the lab notebook.**

### Built for this

- `Spectrum.to_mean_residue_ellipticity(concentration, path_length, residues)`
  — WP2's conversion, and the plan's condition was met: the formula was
  written down and checked against a worked example *before* being coded. None
  of the three inputs is defaulted, because each scales the answer linearly
  and leaves a spectrum that still looks like a protein.
- `structure.helix_from_theta222` — fills `helix` and leaves every other
  category `None`, per ADR-0002 §7.2. Not a second-best decomposition: it is
  the right tool for a spectrum that dies at 210 nm.

### Still blocked, and now for one reason instead of two

The shape deconvolution needs a **reference basis**, and none ships
(ADR-0002 §9). The data limitation is resolved — the 4 µM scan reaches
196 nm — so the only remaining blocker on WP4 is choosing a basis and checking
its terms. AqpZ is a good first test case precisely because the answer is
known: 73 % helix, no sheet.


## Shape, not amplitude (James, 2026-08-09)

**No concentration supplied, deliberately.** The instruction is an estimate
from shape alone, on the grounds that over-reliance on amplitude is a standing
weakness of UV-CD analysis — and it is right, which retires the path-length
argument above. That argument assumed the crystal structure in order to infer
a cuvette, then reported agreement with the crystal structure. As evidence
about the protein it is circular; it is worth keeping only as a way of
recovering a lost instrument setting, and it is now labelled as that.

Everything below is a wavelength or a ratio, so multiplying the spectrum by
any constant leaves it unchanged.

| | 20 µM | 4 µM |
|---|---|---|
| usable range (HT < 600 V) | 209.7–250 nm | 196.2–250 nm |
| zero crossing | not reached | **203.0 nm** |
| minimum | 209.7 nm — **the crop itself** | 209.4 nm |
| theta222/theta208 | 0.947 (208 nm at the edge) | **1.007** |

Two dilutions, five-fold apart, agree on normalised shape to within about 5 %
across 210–250 nm. That is an internal reliability check with no amplitude in
it at all, and it is worth more than either spectrum alone.

**What the shape supports.** A zero crossing at 203 nm and a minimum near
209 nm are what a strongly helical protein looks like; there is no hint of the
~215 nm single minimum a sheet-rich protein would show, consistent with 1RC2's
zero sheet. theta222/theta208 slightly above 1 indicates helices packed
against each other rather than isolated — which is what a transmembrane
bundle is. That reading is weakened by the sample being a detergent-solubilised
membrane protein, where absorption flattening raises the same ratio.

**What the shape does not support: a number.** Turning these into percentages
needs reference shapes, and the package ships none. That is the WP4 blocker
and it has not moved; what has changed is that the data reaches far enough
(196 nm) for a fit to be possible once a basis exists.

`structure.cd_shape_descriptors` computes the above, and refuses to let the
edge artefact pass silently: on the 20 µM spectrum the most negative point is
the crop, and reading it as a band position would be reading the detector's
limit as a property of the protein.


---

## Reference sets: what is actually redistributable (checked 2026-08-09)

**Nothing, today.** Not "probably fine" — checked, and the answer is that no CD
reference set was found that SpectroscoPy may ship.

### Everything traces to one source

SP175 and SMP180 are distributed through the **PCDDB** (accessions
CD0000001000–CD0000071000 and CD00000099000–CD00000128000). BeStSel's basis was
optimised against SP175 plus further β-rich spectra *also* deposited in PCDDB.
SSCalcPy bundles SP175 and SMP180 "obtained from the PCDDB". So PCDDB's terms
govern nearly the whole field.

### PCDDB grants access, not reuse

Its terms and conditions are five clauses. They cover purpose, an as-is
disclaimer, the right to change content, browser compatibility, and this:

> *(4) It is intended that retrieval of information from the PCDDB will be
> freely accessible to all, without subscription.*

That is a statement about **access**, not about what may be done with the data
afterwards. There is no licence grant, no permission to copy, modify or
redistribute. Copyright is asserted over "the design and implementation of this
site" — which does not address the data either way. re3data records both the
database licence and the data licence as **"other"**.

There is also a citation condition, which is easy and right to honour:

> *A condition of use of the data in this website is that any publication or
> presentation using any data from the PCDDB must cite both the original
> reference for the data and the PCDDB.*

**Silence is not permission.** MPL-2.0 lets anyone redistribute this package,
including commercially; shipping data whose terms never grant that would be
assuming a licence nobody wrote.

*Caveat: the live site was unreachable and this is read from an archived copy
of the terms page. Confirm against the current text before relying on it.*

### SSCalcPy is doubly unusable

| | licence | why it fails |
|---|---|---|
| `AU-SRCD/SSCalcPy` | CC-BY-**NC**-SA 4.0 | NonCommercial bars a class of user MPL-2.0 must allow; ShareAlike would force its own terms on the package |
| `AU-SRCD/SSCalcPy-mAb` | CC-BY-NC-**ND** 4.0 | NoDerivatives forbids adapting it at all |

This is the same shape as the GPL-3.0 `.spc` reference implementation: useful to
**read and compare against**, impossible to incorporate. The comparison against
their SELCON3 in WP5 stands; vendoring their code or data does not.

### The way through, and it is one this project has already used

**Fetch, do not ship.** `scripts/fetch_spc_fixtures.py` downloads the Galactic
sample files into a gitignored directory, and `tests/conftest.py` skips when
they are absent. The same arrangement works here and is strictly better than
bundling:

- the user obtains the data from PCDDB themselves and accepts PCDDB's terms
  directly, so SpectroscoPy redistributes nothing;
- the data stays current, and citation metadata comes with it;
- **no API change is needed** — `from_cd(basis=...)` already takes a
  caller-supplied basis, which is why it was built that way.

**And ask.** A direct request to Wallace's group at Birkbeck for permission to
redistribute SP175 under a stated licence is likely to succeed for an academic
open-source tool, and costs an email. It worked for the PyPI name: ask rather
than assume, and take the answer. Until an answer arrives, the fetch route is
what ships.

### Consequence for the plan

WP4's blocker does not lift, but it changes shape. It is no longer "find a
redistributable set" — there isn't one — it is **build the fetch-and-load path,
and separately ask for permission**. The first is unblocked work; the second is
an email and a wait.


---

## Fetch-load built, and the first real-data analysis (2026-08-10)

### Built

`library.load_basis(manifest, directory)` — deliberately **not**
PCDDB-specific. It reads whatever `spc.read` reads, so a basis measured in
your own lab, one from a supplier and one downloaded from a public bank all
load identically. A small CSV declares what each file is: a `category` column
for a structural basis, or `helix`/`sheet`/`turn`/`other` fraction columns for
reference proteins, plus `source`/`citation`/`accession` carried into the
result so a composition can say where its basis came from. It catches
percentages given where fractions were meant, and a manifest that declares
both kinds or neither.

`scripts/fetch_pcddb_basis.py` downloads a named set into gitignored
`data/cd_reference/`, verifying each entry looks like a spectrum before
keeping it, and writes a manifest with accessions so the citation condition
can be met.

**The fetcher has never run against a live PCDDB.** The site has been
unreachable for two days — DNS resolves, both ports time out — so the URL
pattern and file layout are inferred from the site's own metadata. Under the
working agreement that makes it provisional, and it is marked as such in its
own docstring: the `.dpt` precedent is a reader written against a
specification where the format turned out not to be what it looked like. It
fails loudly per entry rather than writing rubbish quietly.

So **the deconvolution still cannot run**, and now for a purely practical
reason rather than a licensing one.

### What the real data does support, and two defects it found

Running the melt tools on the actual AqpZ series, the first time any of them
has seen real measurements:

| | 20 µM, 31 spectra | 4 µM, 7 spectra |
|---|---|---|
| usable at every temperature | 213–280 nm | 199–280 nm |
| Tm | **82.4 ± 0.2 °C** | fit refused |
| dH | 842 ± 58 kJ/mol | — |
| third singular value / noise | **8.2 — not two-state** | 2.0 |

**Two defects, both found by real data and both fixed.**

`isodichroic_point` reported a tight crossing at 238.8 nm. It is not one. That
is the red end, where every spectrum has decayed towards zero and so they
trivially agree — spectra agreeing because there is no signal is the absence
of information, not evidence of two states. A crossing is where the series
*inverts*, so the criterion is now the sign change of the correlation between
signal and temperature, restricted to wavelengths carrying real amplitude. On
this data it now returns *not tight*, correctly.

`two_state` accepted the 4 µM series — seven temperatures for six parameters
— and returned a Tm with a standard error of **10¹⁷ °C**: the fit saying it
has no idea while still printing a number. The minimum is now ten points with
an explanation, and a degenerate covariance is reported rather than dressed
up.

### On the melt disagreement

**Tm = 82.4 ± 0.2 °C on the 20 µM series.** The notebook fits at 222 nm and
on the 208/222 ratio with a two-component model, masked to T ≥ 60 °C and
seeded at 85. Those are not obviously far apart, and my earlier numbers were
against synthetic curves where I chose the answer — so there may be no
disagreement to resolve. What is worth comparing is the fitted values, on this
data.

**The rank test says this unfolding is not two-state** (third component at 8.2
times the noise floor), which agrees with the notebook already using a
two-component model rather than a single transition. Whatever Tm either of us
quotes describes one part of a process with more than two states in it.


---

## WP4 unblocked, and the analysis (2026-08-11)

**James found the datasets: <https://github.com/pcddb/DichroWebGit>.** The
`pcddb` organisation publishes SP175, SMP180 and IDP175 there under the **MIT
licence** — © 2023 Andy Miles, a co-author of the PCDDB and SP175 papers. MIT
permits redistribution provided the copyright notice travels with the data, so
the licensing blocker is gone. The repository README also **documents the file
format**, which satisfies the working agreement that a reader be written
against something real.

`library.load_dichroweb_basis` reads the four-file layout: `A.txt` (spectra in
columns, 240 nm downward), `F.txt` (fractions per category), `lbl1.txt`
(category names), `lbl2.txt` (protein names). Verified: SP175 is 71 proteins
over 240–175 nm, SMP180 is 128 over 240–180 nm, compositions summing to
0.99–1.01, in Δε per residue.

### The analysis, and what it took to get right

Fitting AqpZ against the whole of SP175 gave **helix 0.339, sheet 0.257** at an
rmsd of 0.6% — against a crystal truth of 0.701 and *no sheet at all*. Three
things were wrong, in increasing order of interest.

**1. Rows are not information.** Resampling a 1 nm reference set onto the 0.1
nm measurement made a design matrix of 439 rows whose numerical **rank was
46**, for 71 unknowns. The fit was one of infinitely many, and its residual
looked excellent. `_cd_design` now resamples onto the *coarser* of the two
grids, and `from_cd` refuses when the rank cannot determine the references.

**2. The classical self-consistency test needs the amplitude to be right.**
A CDSSTR-style subset search accepts solutions whose weights sum to ~1 — which
only holds if sample and basis are in the same units. Here the sample is mdeg
and the basis Δε, so the weights summed to 2.9 and **not one subset in 20 000
was accepted**. James's objection to amplitude-dependence lands on the standard
method itself, not only on θ222.

**3. Subset averaging regresses to the reference-set mean.** With acceptance on
shape alone, 28 000 subsets passed and returned helix 0.405 — against SMP180's
own set mean of **0.326**. Tightening from rmsd < 3% to the best 0.1% moved it
only to 0.452, where it plateaued. The pull is towards the average protein, and
it is strongest exactly when a protein is unusual, which is when the answer
matters.

### What works: ask what it looks like, not what it decomposes into

`structure.nearest_references` normalises both spectra to unit length and ranks
the reference set by shape. No inversion, so no underdetermined system, and no
amplitude anywhere.

| similarity | protein | helix | sheet |
|---|---|---|---|
| 0.9957 | INS | 0.67 | 0.00 |
| 0.9952 | FERR | 0.75 | 0.00 |
| 0.9947 | **LACY** (lactose permease) | 0.69 | 0.00 |
| 0.9940 | DHQS | 0.53 | 0.18 |
| 0.9931 | **MSC** (mechanosensitive channel) | 0.53 | 0.03 |
| 0.9927 | **NPSRII** (sensory rhodopsin II) | 0.77 | 0.02 |
| 0.9921 | EC1 | 0.65 | 0.00 |
| 0.9906 | **LEUT** (leucine transporter) | 0.77 | 0.01 |

**Similarity-weighted: helix 0.669 ± 0.092, sheet 0.030 ± 0.058, turn 0.120,
disorder 0.182.**

**Crystal 1RC2 with the 23-residue tag: helix 0.701, sheet 0.000.**

Agreement to three percentage points on helix, and the sheet content is
correctly near zero. Four of the eight nearest neighbours are α-helical
membrane transporters and channels — the method has recognised what kind of
protein this is, from shape alone, with no concentration supplied.

### What this does not settle

The neighbour list is evidence, not a decomposition: two folds can share a
far-UV shape, and the ±0.09 spread is the honest uncertainty. A proper
CDSSTR or SELCON — with the selection and self-consistency rules that make
subset fitting work — is still unbuilt, and the finding above says plainly why
a naive version should not be shipped in its place.
