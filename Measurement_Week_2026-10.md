# Measurement week, 26–30 October 2026

> **Superseded by `Measurement_Week_2026-11.md`.** This week was not used, and the scope was cut on 23/10. Kept for its reasoning.

The operational half of `Reference_Spectra_Wanted.md`: what to run, in what
order, on which day. That document says *why* each sample earns its place and
stays the reference for the reasoning; this one is meant to be at the bench.

Dated on purpose. A generic worklist rots; a worklist for one week is either
used or thrown away.

Scheduling for this work lives in OpenProject under the codename **Farandole**,
which is where dates, order and dependencies are owned. **This repository cites
the codename and never a work-package number** — the numbering belongs to that
side and has already drifted twice between the two records. Cite nothing rather
than something that drifts.

Two jobs share the week:

- **A — UV-Vis references**, per `Reference_Spectra_Wanted.md`. Since
  2026-08-31 this is a **1.0.0 release criterion**, not a nice-to-have: `unmix`
  and `scattering` are currently validated against invented data only. Roadmap
  §23.5.
- **B — the amide I diagnosis** (roadmap §20). Answering why the FTIR estimator
  gives a ±20-point spread. **Part of this is a measurement, not an analysis** —
  see B1 below, which is the point of doing it in a lab week rather than at a
  desk.

---

## Before the week — the things that ruin it if left

| | When | |
|---|---|---|
| ✅ | 2026-09-01 | **Instrument free for the week** (James). The bench-access risk of §14.4 is closed |
| ✅ | 2026-09-01 | **Dichromate available; holmium not found** (James). Enough to proceed — see 1.1 below, which now splits into two checks rather than one |
| ✅ | 2026-09-08 | **Neodymium filter found** (James), in the instrument's own kit. *That* item needs no buying and no lead time. See 1.1 for what it does and does not cover |
| ☐ | **15–18/09** | **Order the consumables — the lead time is here, not on the filter.** Beads in three sizes, BSA, lysozyme, calf thymus DNA, yeast RNA, ATP or NAD⁺, latex or silica blanks. **Order the beads first:** three defined sizes, and one size out of stock costs 3.1 and 3.3, two of the five starred items. **Do not buy a wavelength standard** — the kit's neodymium filter is the standard, see 1.1 |
| | | *Moved from 01–02/10 (James, 2026-09-08).* Nothing on the list perishes and the order has no predecessor, so the only thing a later date bought was risk. Six weeks buys a named thing that three weeks and a bit did not: one full re-order cycle if a bead size comes back out of stock |
| ☐ | Before | Say the word on the §20.2 diagnostic script — one page per series, raw and reference, subtraction with the 2130 cm⁻¹ band marked, second derivative with detected positions, fit with components and residual, composition spread. It exists to make Wednesday about looking at spectra rather than writing plumbing |
| ☐ | Before | Fresh buffer, in quantity. Every blank in the week should come from the same bottle |

## Ground rules — the difference between usable and not

These are not fussiness; each one has cost this project something already.

- **Record the path length.** Every concentration downstream is wrong by the
  ratio if a 1 mm cuvette goes unrecorded. Put it in the filename if there is
  nowhere else.
- **Record concentrations as prepared, with the dilution arithmetic** — not the
  round number you meant to hit.
- **Save the blank as its own file.** Do not only subtract it in the instrument
  software. §21 found OPUS files where a subtraction had already happened and
  nothing in the export said so.
- **Native format plus a text export** of the same measurement wherever the
  instrument offers both. That pairing is the only reason the OPUS reader can be
  trusted.
- Same instrument, same day, same slit and scan speed **within a section**.
- **A file whose concentration was not written down is not a reference. It is
  just a spectrum.** This is the one thing worth being pedantic about.

---

## Monday 26/10 — instrument, then the foundation

Nothing measured later is believable until 1.1 passes.

**1.1 Instrument validation.** The original list said "dichromate *or* holmium
oxide" as if they were alternatives. They are not — they validate different
axes, and only one of them is in hand:

| Axis | Standard | Status |
|---|---|---|
| **Photometric accuracy and linearity** | Potassium dichromate, 5 concentrations in dilute acid | ✅ available |
| **Wavelength accuracy** | Neodymium glass filter, from the instrument's own validation kit | ✅ found 2026-09-08 |

**Dichromate alone is enough to run the week**, because the photometric axis is
the one this project's claims rest on: `from_series` inverts Beer–Lambert and
`unmix` fits amplitudes, so an absorbance-scale error propagates into every
number the library produces. A small wavelength error does not.

Both axes now have a standard. Run and **record** both on the Monday; the point
of 1.1 is a record that the instrument was right *on the day*, not a belief that
it generally is.

**Using the neodymium filter.** Two things decide whether the numbers mean
anything:

- **Take the band positions from its certificate or datasheet, not from a
  textbook.** Rare-earth glass filters differ between batches, and the
  certified values are quoted *at a stated spectral bandwidth*.
- **Record the slit width, and use the one the certificate assumes.** The
  apparent maximum of a rare-earth band shifts with bandwidth, so a filter
  measured at the wrong slit disagrees with its own certificate and tells you
  nothing about the monochromator. This is the same class of error as an
  unrecorded path length, and it costs the check rather than the sample.

**What it does not cover, and what to do about it.** Neodymium's sharp bands sit
in the visible. Everything this project actually claims is in the ultraviolet —
A260, A280, the DNA edge, `best_wavelengths()` choosing 224.5 nm. So the filter
verifies the monochromator where the work is not, which is worth stating in the
`data/uvvis/` README rather than glossing.

It is still worth having, for a reason that should be given rather than assumed:
a grating monochromator's wavelength error is a smooth function of drive
position, not a local defect, so agreement at several visible points is real
evidence about the UV — it is an extrapolation, but a well-founded one. Two
things narrow the gap further and both are free:

- **The deuterium lamp lines** at 656.1 and 486.0 nm, if the instrument's
  built-in test uses them.
- **The dichromate spectrum you are measuring anyway.** Its maxima near 257 and
  350 nm and minima near 235 and 313 nm are in the ultraviolet, and while the
  bands are too broad for a calibration they are a perfectly good coarse check
  that the UV end has not moved by several nanometres. It costs nothing: note
  the four positions when the photometric series is run.

If the dichromate is bench reagent rather than a certified material (NIST SRM
935a or equivalent), say so in the `data/uvvis/` README. An uncertified dilution
series still gives the **linearity** check, which is the more useful half here —
it is what says the absorbance scale is trustworthy across the range the
mixtures will actually use. What it cannot give is traceable absolute accuracy,
and no claim in §23 currently needs that.

**If either check disagrees with its expected values, stop and fix the
instrument** — everything downstream inherits the error.

Where wavelength accuracy bites hardest is **2.1**: A260/A280 is defined at two
specific wavelengths, and `best_wavelengths()` picks points on a steep part of
the DNA edge, where a nanometre of drive error moves the value. The *recovered
mixture ratios*, which are what 2.1 is actually testing, come from a whole-
spectrum fit and are far less sensitive. If the UV cross-checks above disagree
with each other, that ordering is what to fall back on: report the recovered
ratios, and mark the two-wavelength numbers provisional.

Then the extinction-coefficient standards. Each is a **dilution series of 5
concentrations spanning at least a factor of ten**, plus its buffer blank as a
separate file. Five is the minimum that makes a fitted uncertainty mean
anything.

| | Sample | Range | Note |
|---|---|---|---|
| 1.2 ★ | **BSA** | 0.05–1.0 mg/mL | ε₂₈₀ computable from sequence, so the fit has something to be checked against |
| 1.3 | **Lysozyme** | 5 concentrations | Very different Trp/Tyr ratio — tests that sequence-derived ε₂₈₀ generalises |
| 1.4 ★ | **Calf thymus DNA** | 5–50 µg/mL | A260 = 1 should come out at 50 µg/mL. That is a free check on the whole day |
| 1.5 | **Yeast RNA** | 5 concentrations | Makes RNA-vs-DNA separation testable |
| 1.6 | **ATP or NAD⁺** | 5 concentrations | Sharp, well-known ε — tests calibration on something that is not a broad polymer band |

## Tuesday 27/10 — the actual tests

**Section 2 — known mixtures.** Unmixing can only be validated against
mixtures whose composition you set. This is the section the exercise exists for.

| | Sample |
|---|---|
| 2.1 ★ | **BSA + DNA**, ~5 ratios from pure protein to pure DNA. The headline test: recover the ratios |
| 2.2 | One mixture at **three path lengths** (1 cm, 2 mm, 1 mm). Fastest way to catch a units slip |
| 2.3 | **BSA + DNA + a third absorber** (ATP, or haemoglobin). Unmix against only two references and confirm the residual shows the third at its own wavelength — currently demonstrated only synthetically |
| 2.4 | One mixture **on a second instrument**. Tests that references transfer between machines, which is the point of having a library |

**Section 3 — scattering.** The correction is fitted where nothing absorbs and
extrapolated, and real scattering is not a single power law.

| | Sample |
|---|---|
| 3.1 ★ | **Polystyrene beads**, 3 sizes (~100 nm, 500 nm, 2 µm) × 3 turbidities. Scattering with no absorption at all. The exponent should fall as size rises — if it does not, that is a finding |
| 3.2 | **Latex or silica blanks** matched to a real sample's turbidity, for `scattering.from_references` |
| 3.3 ★ | **BSA + beads**, known protein concentration, 3 turbidities. The real test: does correcting scattering recover a concentration you already know? |
| 3.4 | **Bacterial membrane fraction**, dilution series. The honest hard case — no ground truth, but what the method is actually for |
| 3.5 | Any of the above **before and after clarification** (spin or filter). Gives a "true" spectrum to compare the correction against |

## Wednesday 28/10 — FTIR, and the question that is a measurement

**B1. Remeasure a lysozyme concentration series with drying controlled.**

This is the part worth crossing town for. §20.1 asks, as the *first* question,
whether the four concentrations in §19 differ physically — the films held
varying water, subtraction factors ran 0.05 to 0.47 across the series. If they
were measured at different stages of drying they are not four measurements of
one thing, and §19's invariance test was testing an assumption that does not
hold.

**That cannot be answered from the existing files, because nobody recorded the
drying state.** It is answerable only by measuring it again on purpose:

| | |
|---|---|
| ☐ | One lysozyme stock, four concentrations in the ratio 1 : 2 : 5 : 10, matching Candice Gomez's series |
| ☐ | Each film measured at the **same recorded drying time**, timed rather than judged. This is the controlled version of §19 — if composition is now invariant across the series, §19 measured the samples and not the method |
| ☐ | Then one concentration as a **deliberate drying time-course** — several spectra as it dries, times recorded. This measures how much composition moves with water content, which is the size of the error the uncontrolled series was subject to |
| ☐ | Buffer and H₂O blanks from the same session, saved separately |

Cheap and worth having on the same films: your own **hand-tuned subtraction
factor judged by flatness at 1800–2000 cm⁻¹**, recorded alongside the value the
2130 cm⁻¹ anchor picks. §20.1 question 3 is whether the water subtraction is
recoverable, and that comparison answers it directly.

**B2. Then the desk half of §20.1**, in order, with the diagnostic script if it
is built: are the spectra good enough to decompose at all (S/N is 600:1, but
amide I after subtraction is ~0.02 A on a water background of 0.3, and residual
water structure is the size of the sub-bands being fitted); how many components
does the amide I actually support (automatic detection found two or three,
convention says five to seven); and only then, do the band assignments need
changing.

## Thursday–Friday 29–30/10 — fold in, and section 4

**The measurement is not the deliverable.** Until the files are in the test
suite with concentration, path length and blank recorded, the day has bought
nothing — and 1.0.0 is now waiting on it.

| | |
|---|---|
| ☐ | Files into `data/uvvis/`, committed. Unlike the `.spc` fixtures these are ours and small, so there is no licensing or size reason to keep them out |
| ☐ | A `README` in that directory giving, **per file**: sample, concentration, path length, instrument, date, and what the blank was |
| ☐ | Real-data tests for `unmix` and `scattering` alongside the synthetic ones — the synthetic tests stay, they prove the arithmetic |
| ☐ | Update roadmap §23 and §0 with what the real data showed, including anything it broke |

**Section 4 — collect if the week allows.** Measuring and analysing are
separable; do not skip a measurement because its analysis is unbuilt.

| | Sample | |
|---|---|---|
| 4.1 | **Protein melt, 20–90 °C** | Feeds §16's titration work, which has been "waiting on a dataset" since 2 August. `processing.titration` being unbuilt is a reason not to analyse it yet, not a reason to leave a warm instrument unused |
| 4.2 | Same sample at **1 mm and 10 mm, deliberately mislabelled** | A test fixture for catching the error, not for measuring correctly |
| 4.3 | A deliberately **saturated** spectrum, A > 3 | Beer–Lambert stops being linear; the library should say so rather than return a confident wrong number |
| 4.4 | **Empty cuvette, buffer, water** | Baseline behaviour, and useful for teaching |

---

## If the week is short

**1.2, 1.4, 2.1, 3.1, 3.3** — BSA and DNA references, a mixture series to test
them against, beads to test scattering with nothing absorbing, and beads plus
protein where the answer is known. About 30 spectra: half a day, and it is the
half that turns synthetic tests into real ones. Everything else strengthens the
case rather than making it.

Add **B1** to that list regardless of how short the week is. It is half a day
and it decides whether §19 measured the method or the samples, which is the
question the whole FTIR estimator is currently stuck behind.

## Rough size of the full list

About **73 UV-Vis spectra** plus blanks, excluding the melt — at four minutes a
sample including handling and washing, roughly five hours of pure bench time, so
two days with preparation is right rather than optimistic. The FTIR half is
about a dozen spectra and most of its cost is the drying time-course, which is
waiting rather than working.
