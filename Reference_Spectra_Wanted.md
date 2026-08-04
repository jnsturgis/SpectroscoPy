# Spectra wanted — a day in the lab

Everything in `processing.unmix` and `processing.scattering` is currently
tested against **synthetic** data: mixtures built from invented references, so
that the amounts recovered can be checked against a known truth. That proves
the arithmetic. It does not prove the methods survive a real instrument, real
baselines, real turbidity, or a real operator.

This is the list that would fix that. It is written to be worked through in a
day, in the order below, because each section supplies what the next one
needs.

**Ground rules that make the difference between usable and not:**

- **Record the path length.** Every concentration downstream is wrong by the
  ratio if a 1 mm cuvette goes unrecorded. Put it in the filename if there is
  nowhere else.
- **Record concentrations as prepared, with the dilution arithmetic**, not as
  a round number you meant to hit.
- **Save the blank separately**, do not only subtract it in the instrument
  software. Section 21 found OPUS files where a subtraction had already
  happened and nothing in the export said so.
- **Native format plus a text export** of the same measurement, where the
  instrument offers both. That pairing is the only reason the OPUS reader can
  be trusted.
- Same instrument, same day, same slit and scan speed within a section.

---

## 1. Extinction coefficient standards — the foundation

For `library.from_series`: a dilution series gives an ε spectrum with a
standard error at every wavelength, and that becomes a reference everything
else is measured against.

| # | Sample | Series | Why |
|---|---|---|---|
| 1.1 | **Potassium dichromate** or **holmium oxide** | 5 concentrations | Certified reference. Validates the *instrument* — photometric accuracy and wavelength calibration — before anything else is believed |
| 1.2 | **BSA** | 5 concentrations, 0.05–1.0 mg/mL | The protein reference. ε₂₈₀ known and computable from sequence, so the fitted value has something to be checked against |
| 1.3 | **Lysozyme** | 5 concentrations | A second protein with a very different Trp/Tyr ratio. Tests that the sequence-derived ε₂₈₀ works generally rather than for BSA |
| 1.4 | **Calf thymus DNA** (or any clean dsDNA) | 5 concentrations, 5–50 µg/mL | The nucleic acid reference. A260 = 1 should come out at 50 µg/mL |
| 1.5 | **Yeast RNA** | 5 concentrations | Distinguishing RNA from DNA spectrally is a real question and the two references make it testable |
| 1.6 | **ATP or NAD⁺** | 5 concentrations | A small molecule with a sharp, well-known ε. Tests the calibration on something that is not a broad polymer band |

*Five concentrations is the minimum that makes the fitted uncertainty mean
anything; span at least a factor of ten. Include the buffer blank as a
separate file.*

## 2. Known mixtures — the actual test

This is the section the whole exercise exists for. Unmixing can only be
validated against **mixtures whose composition you set**.

| # | Sample | Why |
|---|---|---|
| 2.1 | BSA + DNA, ~5 ratios from pure protein to pure DNA | The A260/A280 case with a known answer. Recovering the ratios is the headline test |
| 2.2 | One mixture at 3 path lengths (1 cm, 2 mm, 1 mm) | Proves the path-length handling, and is the fastest way to catch a units slip |
| 2.3 | BSA + DNA + a third absorber (ATP, or haemoglobin) | The **missing component** test. Unmix against only two references and confirm the residual shows the third at its own wavelength — currently only demonstrated synthetically |
| 2.4 | A mixture measured on two instruments | Tests that references transfer between machines, which is the whole point of a library |

## 3. Scattering — for `processing.scattering`

The correction is fitted where nothing absorbs and extrapolated. It needs real
turbidity, because real scattering is not a single power law.

| # | Sample | Why |
|---|---|---|
| 3.1 | **Polystyrene bead suspensions**, 3 sizes (say 100 nm, 500 nm, 2 µm), 3 turbidities each | Scattering with *no* absorption at all, at known particle size. The cleanest possible test of the power-law basis, and it should show the exponent falling as size rises |
| 3.2 | **Latex or silica blanks** matched to a real sample's turbidity | The measured-background route, `scattering.from_references` |
| 3.3 | **BSA + beads**, known protein concentration, 3 turbidities | The real test: does correcting the scattering recover a concentration you already know? |
| 3.4 | **Bacterial membrane fraction**, dilution series | The honest hard case, and what the two shipped `uvvis_*` datasets already are. Not a ground truth, but it is what the method is actually for |
| 3.5 | Any of the above **before and after clarification** (spin or filter) | Gives a "true" absorption spectrum to compare the correction against |

## 4. Nice to have if the day goes well

| # | Sample | Why |
|---|---|---|
| 4.1 | A protein melt, 20–90 °C | Feeds `processing.titration` (§16), still unbuilt |
| 4.2 | The same sample on a **1 mm and 10 mm** cuvette, deliberately mislabelled | A test fixture for catching the error, not just for measuring correctly |
| 4.3 | A deliberately **saturated** spectrum (A > 3) | Beer–Lambert stops being linear; the library should be able to say so rather than returning a confident wrong number |
| 4.4 | Empty cuvette, buffer, and water | Baseline behaviour, and useful for teaching |

---

## What to do with the files

Real spectra used in tests would be committed if their size allows — these are
small, so unlike the `.spc` samples there is no licensing or size reason to
keep them out. Put them under `data/uvvis/` with a `README` giving, per file:
sample, concentration, path length, instrument, date, and what the blank was.

**A file whose concentration was not written down is not a reference**, it is
just a spectrum. That is the one thing worth being pedantic about on the day.

## Priority if the day is short

**1.2, 1.4, 2.1, 3.1, 3.3.** BSA and DNA references, a mixture series to test
them against, beads to test scattering with nothing absorbing, and beads plus
protein to test the correction where the answer is known. Everything else
strengthens the case; those five turn synthetic tests into real ones.
