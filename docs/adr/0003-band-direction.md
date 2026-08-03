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
