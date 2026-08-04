---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Separating UV-Vis components

How much protein and how much nucleic acid is in this sample? The usual answer
is A260/A280, and it is a weaker answer than its popularity suggests. This page
does it properly, using the whole spectrum against known references — and the
two-component protein/nucleic-acid case turns out to be nothing more than the
smallest instance of that.

Everything here runs when the documentation is built.

```{code-cell}
:tags: [hide-cell]

import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy import viz
from spectroscopy.lineshapes import gauss
from spectroscopy.library import Library, Reference, from_series
from spectroscopy.processing.unmix import (
    unmix, absorbance_ratio, nucleic_acid_and_protein, best_wavelengths)

# ⚠️ SYNTHETIC references, invented for this page so that the numbers below can
# be checked against a known truth. They are the right shape and the wrong
# spectrum: do not use them for anything. Real ones are measured -- see
# "Making your own references" at the end.
x = np.linspace(220.0, 400.0, 361)
dna_eps = 20.0 * gauss(x, 260, 22, 1.0) + 8.0 * gauss(x, 230, 18, 1.0)
protein_eps = 1.0 * gauss(x, 280, 14, 1.0) + 12.0 * gauss(x, 225, 12, 1.0)
haem_eps = 30.0 * gauss(x, 300, 10, 1.0)

def uv(y, name):
    return spc.Spectrum(x, y, x_quantity='Wavelength', x_unit='nm',
                        y_quantity='Absorbance', y_unit='absorbance', name=name)

references = Library([
    Reference('dsDNA', uv(dna_eps, 'dsDNA'), unit='(ug/mL)^-1 cm^-1',
              source='SYNTHETIC — for documentation only'),
    Reference('protein', uv(protein_eps, 'protein'), unit='(ug/mL)^-1 cm^-1',
              source='SYNTHETIC — for documentation only'),
], name='documentation')

# The mixture used throughout: amounts we choose, so every number below can be
# checked against a known answer.
truth = {'dsDNA': 0.030, 'protein': 0.400}          # µg/mL
sample = uv(truth['dsDNA'] * dna_eps + truth['protein'] * protein_eps, 'sample')
```

## The references

A library is a named set of reference spectra. When those references are
**extinction spectra** — ε against wavelength — unmixing returns
concentrations directly rather than relative amounts.

```{code-cell}
print(references)

fig, ax = plt.subplots()
viz.plot_collection(spc.SpectrumCollection([r.spectrum for r in references]), ax)
ax.set_ylabel('ε  /  (µg/mL)⁻¹ cm⁻¹');
```

```{warning}
Those two curves are **invented for this page**. They are the right shape so
that the arithmetic below can be checked against a known answer, and they are
not real extinction spectra for anything. SpectroscoPy deliberately ships no
reference spectra: a fabricated reference makes an unmixing look quantitative
when it is decorative. It ships published *scalar* coefficients, and the means
to measure your own spectra.
```

## First: scattering

A turbid sample — membranes, vesicles, anything not fully dissolved —
attenuates light by scattering as well as absorbing, and the spectrophotometer
cannot tell the difference. **This has to be dealt with before unmixing**, and
here is why:

```{code-cell}
from spectroscopy.processing.scattering import correct_scattering

scatter = 0.6 * (x / 400.0) ** -3.2 + 0.25 * (x / 400.0) ** -1.5
turbid = uv(truth['dsDNA'] * dna_eps + truth['protein'] * protein_eps + scatter,
            'turbid')

spoiled = unmix(turbid, references)
print("true      ", truth)
print("uncorrected", {n: round(float(v), 4) for n, v in spoiled})
```

The nucleic acid comes out **eight times too high**. Scattering is not one of
the components, so the fit spreads it across the ones it has — and R² stays
respectable throughout.

The background is fitted where nothing absorbs (320–400 nm by default for
protein and nucleic acid) and extrapolated:

```{code-cell}
corrected, background = correct_scattering(turbid, return_background=True)

fig, ax = plt.subplots()
turbid.plot(ax, label='measured')
background.plot(ax, label='fitted scattering', linestyle='--')
corrected.plot(ax, label='corrected')
ax.axvspan(320, 400, color='k', alpha=0.07)
ax.legend(fontsize='small');
```

```{code-cell}
print("corrected ", {n: round(float(v), 4) for n, v in unmix(corrected, references)})
```

Rayleigh scattering goes as λ⁻⁴, but real samples are polydisperse and the
exponent falls toward zero as particles approach and exceed the wavelength. So
the default fits a **basis of power laws** rather than one, which spans that
range without anyone having to declare a particle size:

```{code-cell}
from spectroscopy.processing.scattering import scatter_baseline

basis = scatter_baseline(x, turbid.y)                       # the default
rayleigh = scatter_baseline(x, turbid.y, exponents=[4.0])   # fixed exponent
for label, estimate in (('power basis', basis), ('Rayleigh only', rayleigh)):
    print(f"{label:<14} max error {np.max(np.abs(estimate - scatter)):.4f}")
```

Where you have blanks — the same particles without the chromophore —
`scattering.from_references` fits those instead, and that is better still,
because a measured background scatters the way your sample does including
whatever a power law cannot express.

```{warning}
The fit region must be one where your sample really does not absorb. Fit
through a band and the background will eat it; the correction warns when the
region still holds structure afterwards, but it cannot know what your sample
contains. A haem or a carotenoid needs a different window.
```

## Unmixing a known mixture

The mixture was built from known amounts, so the test is whether they come
back:

```{code-cell}
print("true:", truth)
result = unmix(sample, references)
print(result)
```

The amounts are recovered exactly, because this mixture really is a sum of
those two references and nothing else. Real data is not like that, which is
what the rest of this page is about.

`nucleic_acid_and_protein(sample, dna, protein)` is the same call with a name
that says what it is for. There is no separate two-component algorithm.

## Why not A260/A280

```{code-cell}
print(f"A260/A280 = {absorbance_ratio(sample):.3f}")
```

As a flag that a preparation is not what you expected, that number is fine.
As a measurement it is weak, and the reason is worth seeing rather than being
told. Here is the same sample with a third species added — a haem, absorbing
at 300 nm, which neither reference knows about:

```{code-cell}
contaminated = uv(truth['dsDNA'] * dna_eps + truth['protein'] * protein_eps
                  + 0.02 * haem_eps, 'contaminated')

print(f"clean        A260/A280 = {absorbance_ratio(sample):.3f}")
print(f"contaminated A260/A280 = {absorbance_ratio(contaminated):.3f}")
```

Two numbers taken from 361 points cannot see it. The ratio barely moves.

## The residual is the diagnostic

Unmix the contaminated sample against the same two references:

```{code-cell}
bad = unmix(contaminated, references)
print(bad)
```

**R² is still 0.98.** On a smooth spectrum, two broad references can explain
almost anything to two decimal places, so a high R² is not evidence that the
model is right. What gives it away is the residual:

```{code-cell}
fig, axes = plt.subplots(2, 1, sharex=True, figsize=(7, 5),
                         gridspec_kw={'height_ratios': [2, 1]})
contaminated.plot(axes[0], label='measured')
bad.reconstruction.plot(axes[0], label='fitted', linestyle='--')
axes[0].legend(fontsize='small')
bad.residual.plot(axes[1], color='#D55E00')
axes[1].axhline(0, color='k', lw=0.6)
axes[1].set_title('residual', fontsize='small');
```

A band, sitting at 300 nm, exactly where the species that was left out
absorbs. That is what a missing component looks like, and no single ratio can
show it to you.

```{code-cell}
peak = bad.residual.x[int(np.argmax(np.abs(bad.residual.y)))]
print(f"residual peaks at {peak:.0f} nm")
```

## Amounts are never negative

A concentration below zero is not a small negative concentration; it is the
fit compensating for something. Constraining the amounts keeps that pressure
in the residual where it is visible:

```{code-cell}
odd = uv(0.40 * protein_eps - 0.05 * haem_eps, 'odd')
print("non-negative :", np.round(unmix(odd, references).amounts, 4))
print("unconstrained:", np.round(unmix(odd, references,
                                       non_negative=False).amounts, 4))
```

Lifting the constraint is right for a **difference spectrum**, where a
negative coefficient means a component was lost rather than gained.

## Choosing wavelengths

If you must work at a handful of wavelengths — a filter photometer, a plate
reader — the instinct is to pick where each component peaks. That is the wrong
criterion. Two species that both absorb at 280 nm cannot be told apart there
however strongly they absorb; what matters is that their extinction vectors
point in *different directions*.

```{code-cell}
chosen = best_wavelengths(references)
print("chosen        :", chosen)
print("naive (260/280):", np.array([260.0, 280.0]))

for label, pair in (('chosen', chosen), ('naive', [260.0, 280.0])):
    matrix = references.on(np.asarray(pair, dtype=float)).T
    print(f"  {label:<7} condition number {np.linalg.cond(matrix):8.1f}")
```

It picks 260 with something near 225 rather than 260 with 280, because the
protein reference has a strong band at 225 that the nucleic acid one lacks —
so that pair separates the two far better than the traditional one.

Use the whole spectrum whenever you can. Restricting to a few points throws
away precisely the information that revealed the contamination above.

## Making your own references

Measure a dilution series of known concentration and let Beer–Lambert run
backwards: instead of a concentration from an absorbance and a known ε, an **ε
spectrum from absorbances at known concentrations**, with a standard error at
every wavelength.

```{code-cell}
rng = np.random.default_rng(0)
concentrations = [2.0, 5.0, 10.0, 20.0, 40.0]          # µg/mL
series = spc.SpectrumCollection(
    [uv(c * dna_eps + 0.002 * rng.normal(size=x.size), f"{c} µg/mL")
     for c in concentrations])

measured = from_series(series, concentrations, 'dsDNA',
                       unit='(ug/mL)^-1 cm^-1', path_length=1.0)
print(measured)
print("median uncertainty on ε:", f"{np.median(measured.uncertainty):.2e}")
```

```{code-cell}
fig, ax = plt.subplots()
ax.plot(x, dna_eps, 'k:', label='true ε (synthetic)')
ax.plot(x, measured.spectrum.y, label='fitted from the series')
ax.fill_between(x, measured.spectrum.y - 3 * measured.uncertainty,
                measured.spectrum.y + 3 * measured.uncertainty,
                alpha=0.3, label='±3 s.e.')
ax.set_xlabel('Wavelength / nm'); ax.set_ylabel('ε  /  (µg/mL)⁻¹ cm⁻¹')
ax.legend(fontsize='small');
```

The fit goes through the origin, because zero concentration absorbs nothing.
That is deliberate: a systematic offset then shows up as structure in the
residual instead of being quietly absorbed into an intercept, and an offset in
a UV-Vis dilution series usually means scattering or a baseline that was never
removed.

A reference measured this way goes straight back into a library and is used
like any other — which closes the loop.

## Path length

Beer–Lambert is `A = ε·c·l`, so fitting extinction spectra recovers `c·l` and
the concentration needs the division. The default is the standard 1 cm cuvette;
anything else has to be said, because nothing about the result looks wrong when
it is not:

```{code-cell}
print("1 cm  :", np.round(unmix(sample, references).amounts, 4))
print("1 mm  :", np.round(unmix(sample, references, path_length=0.1).amounts, 4))
```

A spectrum that records its own path length is used without being asked:

```{code-cell}
in_a_short_cell = sample._derive()
in_a_short_cell.metadata['path_length'] = 0.1
print(unmix(in_a_short_cell, references).amounts.round(4),
      "from metadata, path length",
      unmix(in_a_short_cell, references).path_length, "cm")
```

## Is ε a unit?

Not a convertible one, and the distinction is deliberate. `spectroscopy.units`
knows about `M^-1 cm^-1` and its relatives, but they are **not** in the
convertible `Y_UNITS` table:

```{code-cell}
from spectroscopy import units

for unit in ('M^-1 cm^-1', '(ug/mL)^-1 cm^-1', 'absorbance', 'cm^-1'):
    print(f"{unit:<18} extinction={units.is_extinction(unit)!s:<6}"
          f" bands={units.band_direction(unit):<8}"
          f" convertible={units.can_convert_y(unit)}")
```

`Y_UNITS` means "can be converted to transmittance". ε cannot: going from ε to
absorbance is `A = ε·c·l`, which needs two numbers that are not properties of
the spectrum. It is a different quantity, not another spelling of absorbance.

It *is* absorbance-shaped, though — bands point up — so it is registered in
`BAND_DIRECTION` and peak-finding on an ε spectrum does the right thing rather
than falling back on an assumption. Note that `cm^-1` is excluded: that is
wavenumber, an x unit.

## Published coefficients

For the single-wavelength cases, the conventional numbers ship with their
sources attached:

```{code-cell}
from spectroscopy import library

print(library.concentration_from_absorbance(1.0, 'dsDNA', 260), "µg/mL at A260 = 1")
print(library.coefficient('dsDNA', 260).note)
```

For a protein whose sequence you know, compute ε rather than using the
1 mg/mL rule of thumb — the rule is wrong by a factor of several for a protein
with no tryptophan:

```{code-cell}
# Lysozyme: 6 Trp, 3 Tyr, 4 cystines
print(library.protein_epsilon_280(6, 3, 4), "M⁻¹ cm⁻¹")
```
