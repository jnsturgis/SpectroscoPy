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

### WP1 — CD as a technique (small, additive, do first)

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

### WP4 — deconvolution

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

### WP6 — the AqpZ temperature melt (§18.3)

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
