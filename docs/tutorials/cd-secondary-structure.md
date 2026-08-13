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

# CD: secondary structure from a far-UV spectrum

You have a far-UV CD spectrum of a protein and you want to know how much of it
is helix. This page does that, and then spends most of its length on the harder
question: **how much to believe the number.**

Everything runs on data that ships with the package — the published reference
sets and a real membrane-protein spectrum — so you can work through it before
loading anything of your own.

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt
import spectroscopy as spc
from spectroscopy.processing import cd

references = spc.datasets.reference_set('sp175')
print(references)
print(references.info['licence'])
print(references.info['citation'])
```

## What a reference set is

Seventy-one proteins whose structures were solved by crystallography, each with
a CD spectrum measured on the same instrument and a composition from DSSP on
its structure. The spectrum and the structure travel together, so there is no
second list to keep in step:

```{code-cell} ipython3
for spectrum, composition in list(zip(references, references.compositions))[:4]:
    print(f"{spectrum.name:<8} helix {composition.get('helix'):.2f}   "
          f"sheet {composition.get('sheet'):.2f}   "
          f"({composition.method})")
```

The shapes are the whole basis of the method — a helical protein and a sheet
protein do not look alike:

```{code-cell} ipython3
pairs = list(zip(references, references.compositions))
helical = max(pairs, key=lambda pair: pair[1].get('helix'))[0]
sheety = max(pairs, key=lambda pair: pair[1].get('sheet'))[0]

fig, ax = plt.subplots(figsize=(7, 3.5))
for spectrum in (helical, sheety):
    ax.plot(spectrum.x, spectrum.y, label=spectrum.name)
ax.axhline(0, color='0.7', lw=0.8)
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel(r'$\Delta\varepsilon$ (M$^{-1}$cm$^{-1}$ per residue)')
ax.legend()
```

A helix gives two negative bands near 208 and 222 nm and a strong positive one
near 193; a sheet gives a single broader minimum near 217. Estimating structure
means asking which mixture of the reference shapes reproduces yours.

## First, on a protein whose answer we know

The honest way to find out whether a method works is to give it a protein it
has not seen. Take one out of the set, estimate it from the other seventy, and
compare with the structure that was known before the fit:

```{code-cell} ipython3
index = 0
unknown = references[index]
truth = references.compositions[index]
rest = references.select(lambda s: s is not unknown)

print(f"{unknown.name}: holding it out leaves {len(rest)} references")
print(f"DSSP says   helix {truth.get('helix'):.3f}   sheet {truth.get('sheet'):.3f}")
```

```{code-cell} ipython3
rows = []
for method in ('nearest-shapes', 'ridge'):
    estimate = cd.estimate(unknown, method, rest)
    rows.append((method, estimate.get('helix'), estimate.get('sheet')))
    print(f"{method:<16} helix {estimate.get('helix'):.3f}   "
          f"sheet {estimate.get('sheet'):.3f}")
```

One of those is close and the other is out by about 0.14 — on a protein whose
answer we know, from a set of seventy references, with no membrane or detergent
anywhere near it. **That is the accuracy of the technique, not a bad day.**

And a single protein cannot tell you which method to trust, because the one
that won here may have been lucky. That is what
{func}`~spectroscopy.processing.cd.benchmark` is for: it repeats the hold-out
across the whole set, so the error is measured rather than sampled from one
draw.

```{code-cell} ipython3
scores = cd.benchmark(references, ['nearest-shapes', 'ridge'], folds=10)
for method, entry in scores.items():
    helix = entry['helix']
    print(f"{method:<16} n={entry['n']:<3} helix rmse {helix['rmse']:.3f}  "
          f"bias {helix['bias']:+.3f}")
```

**Read the bias as carefully as the error.** A method that pulls every protein
towards the average composition of the reference set has a bias that grows with
how unusual the protein is — which is exactly when you needed the answer.

:::{admonition} Why not just fit all 71 references at once?
:class: note

Because it cannot be done. A measurement from 190 to 240 nm carries roughly 50
independent numbers, and 71 unknowns fitted to 50 equations has infinitely many
exact solutions. Least squares returns one of them, with an excellent residual
and an arbitrary composition. `cd.estimate(..., 'all-references', ...)` refuses
rather than answering:

```{code-cell} ipython3
try:
    cd.estimate(unknown, 'all-references', rest)
except ValueError as error:
    print(error)
```
:::

## Now a real one

AqpZ-W14A, an aquaporin — an all-helical membrane channel, measured on a JASCO
J-815 at 4 µM in detergent. This is a spectrum with a real problem in it.

```{code-cell} ipython3
protein = spc.datasets.load('aqpz')
print(protein)
print(f"{len(protein)} points, {protein.x.min():.0f}-{protein.x.max():.0f} nm, "
      f"in {protein.y_unit}")

fig, ax = plt.subplots(figsize=(7, 3.5))
protein.plot(ax)
ax.axhline(0, color='0.7', lw=0.8)
```

### Where the spectrum stops being a measurement

It starts at 197 nm, and that is not where the instrument stopped scanning. The
J-815 records the photomultiplier voltage alongside the signal, and above about
600 V the detector is starved: the trace continues, but it is the detector
struggling rather than the protein absorbing. The scan ran to 180 nm; **only
the part above 197 nm is measurement**, and the rest has been cropped out
before shipping.

That matters more than it sounds, because what separates helix from sheet from
coil lives below 210 nm. A spectrum cropped at 210 has one negative band and no
shape left to fit at all.

### Use a reference set that contains proteins like yours

```{code-cell} ipython3
membrane = spc.datasets.reference_set('smp180')
print(membrane)
print(membrane.info['citation'])
```

SMP180 is SP175 extended with membrane proteins — 30 of its 128. That is not
bookkeeping. A membrane protein's far-UV spectrum resembles *other membrane
proteins* about twice as often as its secondary structure alone would explain:
measured on SMP180, the fraction of a protein's five nearest shapes that are
membrane proteins is 0.45 for membrane queries against 0.20 for soluble ones,
and it stays that way when you control for helix content. Detergent
scattering, absorption flattening and the long straight helices of a
transmembrane bundle all pull the same way. A reference set with no membrane
proteins in it has nothing close to yours.

### Three methods, three answers

```{code-cell} ipython3
import warnings

answers = {}
for method, options in (('selcon', {'ignore_sum_rule': True}),
                        ('nearest-shapes', {}),
                        ('ridge', {})):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        answers[method] = cd.estimate(protein, method, membrane, **options)

print(f"{'method':<16}{'helix':>8}{'sheet':>8}")
for method, estimate in answers.items():
    print(f"{method:<16}{estimate.get('helix'):>8.3f}{estimate.get('sheet'):>8.3f}")
```

The crystal structure of AqpZ is 1RC2: 178 of 231 residues helical, **no sheet
at all**. The construct measured here carries a 23-residue N-terminal tag that
the crystal does not, and those residues absorb too — so the fraction to expect
from a CD spectrum is diluted:

```{code-cell} ipython3
helical_residues, crystal_chain, tag = 178, 231, 23
expected = helical_residues / (crystal_chain + tag)
print(f"1RC2:                 {helical_residues}/{crystal_chain} helical "
      f"= {helical_residues / crystal_chain:.3f}")
print(f"with the {tag}-residue tag: {helical_residues}/{crystal_chain + tag} "
      f"= {expected:.3f} helix, 0.000 sheet")

for method, estimate in answers.items():
    print(f"  {method:<16} off by {estimate.get('helix') - expected:+.3f}")
```

### The disagreement is the result

The methods differ by more than the benchmark says any of them should, and
**the method with the best held-out score is the furthest from the truth.**

That is not a reason to pick the one that agrees. Choosing the method that
matches a structure you already know is selection on a single data point, and
it is how you end up recommending a method that only works on proteins whose
answer you had anyway. The reason to report the disagreement is that it is
telling you something true: **AqpZ is not well described by this reference
set.** A detergent-solubilised membrane channel has absorption flattening that
suppresses the 208 nm band and lifts the 222/208 ratio above one, and no
reference protein measured in buffer does that.

```{code-cell} ipython3
from spectroscopy.processing import structure

shape = structure.cd_shape_descriptors(protein, region=(197.0, 250.0))
for key in ('zero_crossing', 'minimum', 'ratio_222_208'):
    if shape.get(key) is not None:
        print(f"{key:<16} {shape[key]:.3f}")
```

A zero crossing near 203 nm and a minimum near 209 nm are what a strongly
helical protein looks like; there is no sign of the single ~215 nm minimum a
sheet-rich protein would give, which agrees with 1RC2 having no sheet. The
222/208 ratio slightly above one says the helices are packed against each other
rather than isolated — which is what a transmembrane bundle is, *and* what
absorption flattening imitates. Two explanations, one number, and the
measurement cannot separate them.

## What to report

Not one number. The three things worth writing down are:

1. **Which reference set**, because the same spectrum against SP175 and against
   SMP180 is two results, not one.
2. **The spread across methods**, when they disagree by more than the benchmark
   led you to expect.
3. **The wavelength range that was really measured**, and how you decided —
   here, HT below 600 V.

The first of those the result carries itself, so it cannot be lost on the way
into a figure legend:

```{code-cell} ipython3
print(answers['ridge'].quality['references'])
print(answers['ridge'].quality['method'], '/',
      answers['ridge'].quality['n_references'], 'references')
```

### The error bar worth having

Repeat scans give a better error bar than any goodness-of-fit statistic,
because the noise is then *measured*. Four AqpZ scans recorded below the
unfolding transition ship with the package:

```{code-cell} ipython3
repeats = spc.datasets.aqpz_near_native()
sigma = cd.uncertainty_from_replicates(repeats)

for wavelength in (215.0, 222.0, 240.0):
    index = int(np.argmin(np.abs(sigma.x - wavelength)))
    print(f"{wavelength:.0f} nm:  {sigma.y[index]:.3f} mdeg")
```

The error is far from uniform, and it is worst where the detector is working
hardest. An unweighted fit treats a wavelength known to three times the
precision of another as though they were equally good; `sigma=` weights each
by how well it is known.

```{code-cell} ipython3
weighted = cd.estimate(protein, 'ridge', membrane, sigma=sigma, resamples=200)
print(f"helix {weighted.get('helix'):.3f}")
print("moved by, under the measured noise:",
      {k: round(v, 3) for k, v in weighted.quality['uncertainty'].items()})
```

`resamples` refits on the spectrum perturbed by its own noise and reports how
far the answer moves — which is the question an error bar should answer, and
one no residual can.

:::{admonition} These are not true replicates, and the numbers say so
:class: warning

They are the 30, 40, 50 and 60 °C scans of a melt with a midpoint above 80 °C,
so the protein is folded in all four — but any real change with temperature
below the transition is being counted as noise. Look at the shape of the
result: the error at 222 nm comes out *larger* than at 215 nm, which the
detector alone would not do, because 222 nm is exactly where unfolding shows
first. So this overestimates the error, and overestimates it most at the
wavelength you care about. Genuine repeats of one sample at one temperature
would be better, and are what to record if you are planning the experiment.
:::

## Your own spectrum

Two lines change:

| this page | your data |
|---|---|
| `spc.datasets.load('aqpz')` | `spc.read("my_protein.dx")` |
| `spc.datasets.reference_set('smp180')` | the same, or your own via `library.load_basis` |

Crop to where your instrument was still measuring before you fit anything — for
a JASCO that means the HT channel, for other instruments the dynode voltage or
the recorded absorbance of the buffer. Then run more than one method, and if
they agree, say so; if they do not, that is the result.

The methods themselves, what each is sensitive to, and the full held-out
validation on both reference sets are in
[the CD methods guide](../guide/cd-methods.md).
