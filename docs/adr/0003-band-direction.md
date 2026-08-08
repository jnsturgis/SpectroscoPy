# ADR-0003 — Band direction is a property of the y unit

**Status:** Accepted, 2026-08-03 (James: "yes for fix 1 and 2 and a warning for
3", and keep `troughs=` for difference spectra).
**Changes a default** in `Spectrum.find_peaks`, which is why this record
exists.

---

## 1. The problem

`find_peaks` detected maxima. Bands in transmittance are minima. Nothing
anywhere consulted the y unit, so peak-picking a transmission spectrum
returned the wrong answer — and returned it quietly.

One Gaussian band at 1100 cm⁻¹:

| spectrum | `find_peaks()` before | after |
|---|---|---|
| absorbance | **1100.0** ✓ | 1100.0 ✓ |
| transmittance | **1086.0, 1114.0** ✗ | 1100.0 ✓ |

The two spurious positions are the inflection points flanking the band, where
`−d²y/dx²` has its maxima. On the shipped ethanol spectrum this produced 68
plausible-looking positions at plausible wavenumbers — roughly twice the true
band count, which reads as sensitive detection rather than as a failure.

`fit_peaks` inherited it and was worse, because it seeds component positions
from `find_peaks`: the same band at 1100 was fitted with components at 1030.7
and 1169.3. A composition built on that is meaningless and still reports a
respectable R², which is the failure mode roadmap §19–20 is about.

This was found by converting `docs/guide/reading-files.md` to an executing
page and then asking what the peak documentation actually promised.

---

## 2. Decision

**Band direction is a property of the y unit, and belongs in
`spectroscopy.units` rather than in the peak finder.** `BAND_DIRECTION` maps
each unit to `'up'`, `'down'` or — by absence — `'unknown'`, exposed as
`band_direction()` and `is_valley_pointing()`.

Three consequences:

1. **`find_peaks` takes its direction from the unit.** `troughs` defaults to
   `None`, meaning "ask the unit". Transmittance, `%T` and reflectance are
   searched downward; absorbance-like units upward.
2. **An unknown unit keeps the historical upward behaviour and records that it
   was an assumption.** `a.u.` and `counts` dominate Raman and fluorescence,
   where upward is right, so warning on them would be noise. The result
   carries `properties['direction_from'] ∈ {'y_unit', 'caller', 'assumed'}`.
3. **`fit_peaks` warns on valley-pointing units** that areas and ratios are
   not linear there, and points at `.to(y_unit='absorbance')`.

---

## 3. What was rejected

**Warning instead of fixing the default.** A warning on every transmission
spectrum, telling the user to pass a flag to get the right answer, keeps a
wrong default and adds noise. The unit already carries the information.

**Auto-converting in `fit_peaks`.** Tempting, because absorbance is what the
areas need. Rejected: silently changing what a fit was performed on is its own
kind of wrong, and the returned `FitResult` would describe a spectrum the
caller never handed over. A warning naming the fix leaves the decision with
the person who knows what the numbers are for.

**Dropping `troughs=`.** Explicitly kept at James's request, and he is right:
**a difference spectrum has bands in both directions** and the unit cannot say
which was meant. `troughs=False` and `troughs=True` on the same difference
spectrum answer "what appeared" and "what disappeared" — both are wanted, and
an explicit value always overrides the unit.

**Inferring direction from the data** (say, by comparing the mean to the
median). It would work often and fail silently on a spectrum dominated by one
enormous band, on a difference spectrum, or on a baseline-uncorrected one.
The unit is a statement of fact; the shape of the data is evidence.

---

## 4. Prior art

Galactic's SPC format made the same call in 1993. `SPC.H` splits its y-type
codes at 128:

```c
#define YTRANS 128  /* Transmission (ALL HIGHER MUST HAVE VALLEYS!) */
```

Every code at or above 128 is valley-pointing, and the comment is shouted
because every program that draws or picks a spectrum needs it. That a file
format from thirty years ago treats "does this quantity have valleys" as
first-class information about the unit is good evidence this is the right
place for it. See `SPC_Format_Notes.md` §3.7.

---

## 5. Compatibility

`find_peaks(troughs=False)` and `find_peaks(troughs=True)` keep their exact
previous meanings. What changes is the *unspecified* case, and only for units
where the previous answer was wrong. Anybody who had discovered the problem
and was passing `troughs=True` on transmittance is unaffected.

`kind` is still `'peak'` or `'trough'`, describing what was located in the
spectrum. On a transmission spectrum the bands genuinely are troughs, so no
new vocabulary was invented.

---

## 6. What would justify revisiting this

- A unit that is genuinely direction-ambiguous in practice and common enough
  that assuming upward is wrong more often than right. Emission and counts
  were considered and are fine.
- Wanting the direction attached to a *spectrum* rather than a unit — a
  difference spectrum is the obvious candidate, since it is the one case where
  the unit is honestly insufficient. A `Spectrum.band_direction` attribute
  that defaults to the unit's would cover it, and would let a difference
  spectrum say so once instead of at every call.

---

## 7. Amendment, 2026-08-05 — `'both'`, for signed quantities

**Status:** Accepted. Extends §2 rather than replacing it; nothing above is
withdrawn.

§6 named the trigger for revisiting this, and it has fired — from the
direction it predicted, though for a wider class of spectra than it named.

### 7.1 What changed

`band_direction()` now returns **four** values: `'up'`, `'down'`, `'both'`,
`'unknown'`.

`'both'` is for the **signed** quantities — circular and linear dichroism,
fluorescence anisotropy and polarisation, an explicitly labelled difference
spectrum. These are listed in `units.BIPOLAR_UNITS` and answered by
`units.is_bipolar()`.

A dichroism spectrum is a *difference between two measurements* — left minus
right circular polarisation, parallel minus perpendicular — so its sign
carries meaning an absorbance never has. An α-helix shows negative bands at
208 and 222 nm and a positive one at 193 nm, in one spectrum, and none of the
three is a baseline artefact. Searching such a spectrum one way finds half of
it, and the half it finds looks like a complete answer.

### 7.2 Why not `'unknown'`

This was the tempting cheap option and it is wrong. **`'unknown'` means nobody
knows; `'both'` means the answer is known and is *both*.** Collapsing them
would throw away a fact the unit does determine, and would leave CD sharing a
category with `counts` — after which the only available behaviour is the
assumption of upward, which for CD is a silent half-answer.

### 7.3 Consequences

1. `find_peaks` on a signed unit searches both ways without being asked, and
   returns one table containing both.
2. `troughs=` gains a third accepted value, `'both'`. This is the case §6
   anticipated: a difference spectrum in plain `absorbance` **cannot** be
   recognised from its unit, because subtracting two absorbance spectra leaves
   the unit saying `absorbance`. It has to say so at the call site.
3. `PeakTable.kind` gains `'both'`, and the table gains a per-peak `sign`
   array with `maxima()` and `minima()` selectors. §5 previously said no new
   vocabulary was invented; that is no longer true, and the reason is that
   this is genuinely a third case rather than a relabelling of the first two.
4. `sign` is populated for *every* table, not only the mixed ones, so no
   caller has to special-case `kind`.

`sign` is necessary and not merely convenient: **a positive band sitting on a
negative offset still has a negative height**, so height cannot be used to
recover which is which.

### 7.4 Why now, when CD is post-1.0

`band_direction()` freezes at 1.0. Adding a fourth return value afterwards
breaks every caller that wrote `if direction == 'up' ... else ...` against a
documented three. Doing it before costs a line, and the change is useful on
its own account: difference spectra are not a CD feature, and neither is
fluorescence anisotropy.

This is roadmap **D1**, raised by the CD branch plan and settled on `main`
before any CD code exists — which is what §17.3 means by the half built now
not foreclosing the half built later.

### 7.5 Still open

The §6 suggestion of a `Spectrum.band_direction` **attribute**, letting a
difference spectrum declare itself once instead of at every call, is not
adopted here. `troughs='both'` covers the case; an attribute would be the
better answer if labelling difference spectra turns out to be common, and it
remains additive.

### 7.6 Two kinds of difference spectrum (James, 2026-08-05)

`'dA'` is not for every subtraction, and the distinction is about **intent**,
not arithmetic. Both come out of the same operator, which is exactly why
neither can be inferred.

**Corrections** subtract in order to arrive at a *better absolute spectrum* —
a baseline, a buffer blank, a scattering background. The result is still an
absorbance: its bands point up, and a negative excursion is an artefact worth
investigating rather than a measurement. These keep `'absorbance'`, and
should, because peak finding on them should behave exactly as it does on the
uncorrected spectrum.

**Comparisons** put the information *in the sign* — reduced minus oxidised,
ligand-bound minus free, illuminated minus dark. A band pointing down means a
species was lost; one pointing up means a species was gained; both are the
answer, and a tool that reports only one has described half the experiment.
These are `'dA'`.

This is why `subtract_reference` leaves the y unit alone rather than
promoting the result to `'dA'`. It cannot know which of the two you were
doing, and guessing wrong in the correction direction would make every
baseline-corrected spectrum report its noise as bands.
