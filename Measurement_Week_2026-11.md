# Measurement week, 23–27 November 2026

This replaces `Measurement_Week_2026-10.md`, whose week was never used. It is
cut to what the PyPI release needs and to the scope James set on 23/10: **the
scattering standards are dropped** ("rather artificial anyway"). What is left is
section 1 (the instrument and the references) and section 2 (the known
mixtures, which are what the exercise is for). Beads are optional, and only if
some are already in the lab. Section 4 has been removed.

`Reference_Spectra_Wanted.md` explains why each sample is on the list. This
page is meant to be used at the bench. Scheduling is in OpenProject under
**Farandole**.

## Find or buy first

Check the shelves before ordering anything. The September order may never
have been placed, and most of these are routine lab stock.

| | Item | For | If missing |
|---|---|---|---|
| ✅ | Potassium dichromate | 1.1 | in hand |
| ✅ | Neodymium filter + its certificate (band positions **and** the slit width they assume) | 1.1 | in the instrument's kit |
| ☐ | **BSA** | 1.2, 2.1–2.3 | essential: buy it, or do 1.1 alone |
| ☐ | **Calf thymus DNA** | 1.4, 2.1–2.3 | essential: buy it, or do 1.1 alone |
| ☐ | Quartz cuvettes at **1 cm, 2 mm, 1 mm** | 2.2 | without the short paths, skip 2.2 |
| ☐ | Fresh buffer, **one bottle for the whole week** | every blank | make it on the first day |
| ☐ | Lysozyme | 1.3 | skip |
| ☐ | Yeast RNA | 1.5 | skip |
| ☐ | ATP or NAD⁺ | 1.6, 2.3 | skip 1.6, and do 2.3 with haemoglobin or not at all |
| ☐ | A second UV-Vis spectrophotometer you can use for an hour | 2.4 | skip |
| ☐ | Polystyrene beads, **any size already in the lab** | optional B | do not order |

**Only BSA and DNA are worth buying for this week.** Everything else is skipped
if it isn't in the lab.

Also check, before booking the instrument, whether Waltz's work needs the lab
that week.

## Ground rules

These are what make a file usable as test data. They carry over unchanged from
the October plan.

- Every sample has its **blank saved as a separate file**, from the same bottle of buffer.
- For each file, record **concentration, path length, slit width, instrument and date** at the time. A file with an unknown concentration can't be used as a test.
- Use 5 concentrations spanning **at least a factor of ten** for every dilution series. Five is the minimum for a fitted uncertainty to mean anything.
- Keep absorbance **below about 1.5**. Dilute rather than trust a high reading.

## In priority order: stop anywhere, and everything done so far is usable

### 1. The minimum (about half a day)

If the week shrinks to a single morning, do these four. They are enough for
real-data tests of `from_series` and `unmix`.

| | Sample | Detail |
|---|---|---|
| **1.1** | **Instrument check** | Dichromate, 5 concentrations in dilute acid, for linearity and photometric accuracy. Note the maxima near 257 and 350 nm and the minima near 235 and 313 nm: a free UV wavelength check. Then the neodymium filter, **at the slit width its certificate assumes**. If either disagrees with its expected values, stop and fix the instrument first |
| **1.2** | **BSA** | 0.05–1.0 mg/mL, 5 concentrations, plus blank. ε₂₈₀ can be computed from the sequence, so the fit can be checked |
| **1.4** | **Calf thymus DNA** | 5–50 µg/mL, 5 concentrations, plus blank. A260 = 1 should come out at 50 µg/mL, a check on the whole day |
| **2.1** | **BSA + DNA mixtures** | About 5 ratios from pure protein to pure DNA, total absorbance kept in range. **The headline test:** does `unmix` recover ratios you set? |

### 2. If there is time (another half day)

| | Sample | Detail |
|---|---|---|
| **2.2** | One 2.1 mixture at **1 cm, 2 mm and 1 mm** | The quickest way to catch a path-length or units mistake in the library |
| **2.3** | **BSA + DNA + ATP** (or haemoglobin) | Unmix against only BSA and DNA, and check that the residual shows the third absorber at its own wavelength. This has only been shown on synthetic data so far |
| **1.3** | Lysozyme, 5 concentrations | A very different Trp/Tyr ratio, so it tests whether ε₂₈₀ from sequence generalises |

### 3. Only if everything above is done

| | Sample | Detail |
|---|---|---|
| 1.5 | Yeast RNA, 5 concentrations | Makes it possible to test separating RNA from DNA |
| 1.6 | ATP or NAD⁺, 5 concentrations | A sharp, well-known ε, unlike the broad polymer bands |
| 2.4 | One 2.1 mixture **on a second instrument** | Tests whether references transfer between machines |
| B | BSA with **whatever beads are in the lab**, one size, 2–3 turbidities, plus the beads alone | A real-data check of `scattering`. Do it only if beads turn up; it is no longer worth ordering them |

**Not on the list.** The bead size series (old 3.1), latex/silica blanks (3.2),
the membrane fraction (3.4), clarification (3.5) and all of section 4 are
dropped. The FTIR lysozyme drying series (old Wednesday B1) is a separate
question: it matters for the experimental FTIR estimator, not for the release,
so it should get its own day when FTIR is picked up again.

## After the bench: the files are the deliverable

Until the data is in the repository with its metadata, the measuring has not
produced anything usable.

| | |
|---|---|
| ☐ | Files into `data/uvvis/`, committed |
| ☐ | `data/uvvis/README.md` giving, **per file**: sample, concentration, path length, slit width, instrument, date, and the blank. State whether the dichromate is certified material or bench reagent, and that the neodymium check covers the visible range, not the UV |
| ☐ | Real-data tests for `from_series` and `unmix` next to the synthetic ones (Claude can write these from the README) |
| ☐ | Roadmap §23 and §0 updated with what the real data showed, including anything it broke |
