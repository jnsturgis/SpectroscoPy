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

# Secondary structure from CD: four methods, and how to choose

Estimating structure from a far-UV CD spectrum means writing the measurement
as a combination of reference spectra whose structures are known. Every
published method does that. They differ in how they cope with the one hard
fact of the problem:

**A reference set has 70–130 proteins. A spectrum from 190 to 240 nm carries
about 50 independent numbers.** Fitting 130 unknowns to 50 equations has
infinitely many exact solutions, and least squares will hand you one of them
with an excellent residual and no warning.

```{code-cell}
import numpy as np
import spectroscopy as spc
from spectroscopy.processing import cd

print(sorted(cd.METHODS))
```

## The algorithms, briefly

**`all-references`** — non-negative least squares against every reference at
once. Right when the data can determine them, and this **refuses** when it
cannot rather than returning an arbitrary member of the solution set. It is
the honest baseline, not a recommendation.

**`nearest-shapes`** — do not invert anything. Normalise both the measurement
and each reference to unit length, rank the references by cosine similarity,
and average the known structures of the closest few weighted by similarity.
No system to invert means no arbitrary solution to pick. It cannot describe a
protein as a mixture the set does not contain, which is both its weakness and
the reason it does not hallucinate one.

**`subset-average`** — the idea behind CDSSTR ([Johnson 1999](../references.md)):
fit many small random subsets, keep those that pass a self-consistency test,
average the survivors. Eight references *are* determined by fifty numbers even
when a hundred and thirty are not.

**`selcon`** — the self-consistent method of Sreerama & Woody, with SELCON3's
variable selection and rules. Three ideas together: the unknown's own spectrum
joins the basis carrying a *guess* at its structure, the system is solved by
SVD and the guess replaced by the solution until it stops moving; references
are ordered by closeness to the query and increasing numbers of the closest
are tried; and a candidate is kept only if its fractions sum to within 5 % of
one and none falls below −0.025. Surviving solutions are averaged.

**`ridge`** — Tikhonov-regularised non-negative least squares over the whole
set, the standard statistical treatment of an underdetermined system and the
same idea as CONTIN ([Provencher & Glöckner 1981](../references.md)). Penalising the size of the
coefficients picks one solution among the many that fit, and the penalty is
what decides which.

## The constraint that decides half of it: Δε or mdeg?

A reference set is in **Δε per residue**. Your spectrum is in **millidegrees**
until you supply a concentration, a path length and a residue count. That
difference is not cosmetic:

| method | needs the amplitude? | what it is sensitive to |
|---|---|---|
| `all-references` | no | conditioning |
| `nearest-shapes` | no | how well the set covers the fold |
| `subset-average` | **only with `sum_to_one=True`** | the acceptance threshold |
| `selcon` | **always** | the amplitude, unavoidably |
| `ridge` | no | the penalty |

**`selcon` cannot be made scale-free, and this is worth being clear about.**
Its self-consistent step puts the query into the basis *as a column beside the
references*, so the query's magnitude relative to them is part of the model.
Measured on a synthetic set whose references peak near 30, the same shape
scaled by 0.1, 1, 10 and 800 gave helix **0.350, 0.562, 0.558, 0.550** — and
the run without `ignore_sum_rule` refused outright at 0.1 and 10 while
answering at 1 and 800. Renormalising the fractions does not remove the
dependence; it only hides the refusal.

The classical self-consistency test — accept a subset if its weights sum to
about one — is **comparing an amplitude**. It is meaningful only when your
spectrum is already in the reference set's units.

This is not hypothetical. Running `subset-average` with `sum_to_one=True` on a
spectrum in millidegrees against a basis in Δε, the weights came out summing
to 2.9, and **not one subset in twenty thousand was accepted**. The ratio is
not an error in the fit; it is the amplitude scale nobody supplied.

```{code-cell}
:tags: [hide-input]
print("A path length entered as 1 cm when the cell was 1 mm moves any")
print("amplitude-dependent answer by ten, and the spectrum still looks fine.")
```

So: if you have a trustworthy concentration and path length, amplitude
information is real information and worth using. If you do not — or if you
suspect them — the shape-only methods still work, and that is the whole reason
`from_cd` and everything here is scale-free by construction.

## Choosing between them, honestly

Not by which one agrees with a protein you already know. **Any of them will
agree with something.** Choosing a method because it matched your one
well-characterised sample, and then explaining why it should have, is
selection on a single data point.

`benchmark` hides a fraction of the reference set, estimates those proteins
from the rest, and compares with structures that were known before the fit.

```{code-cell}
from spectroscopy.library import ReferenceSet
from spectroscopy.processing.structure import Category, Composition

# A small synthetic reference set, so this page runs without downloading
# anything. The real numbers are in the table below.
x = np.linspace(190.0, 240.0, 51)
def band(centre, width, height):
    return height * np.exp(-((x - centre) ** 2) / (2 * width ** 2))

HELIX = Category('helix', frozenset({'G', 'H', 'I'}))
SHEET = Category('sheet', frozenset({'E', 'B'}))
COIL = Category('other', frozenset({'S', '-'}))
pure = {HELIX: band(208, 7, -37) + band(222, 9, -37) + band(193, 7, 60),
        SHEET: band(217, 9, -18) + band(195, 8, 32),
        COIL: band(198, 7, -40) + band(220, 10, 3)}

rng = np.random.default_rng(0)
spectra, compositions = [], []
for index in range(24):
    weights = rng.dirichlet([1.4, 1.0, 1.0])
    y = sum(w * pure[c] for w, c in zip(weights, pure))
    spectra.append(spc.Spectrum(x, y + 0.3 * rng.normal(size=x.size),
                                technique='CD', name=f'ref{index}'))
    compositions.append(Composition(
        fractions=dict(zip(pure, weights)), method='known', technique='X-ray'))

# The spectra and their known structures are paired once, here, and travel
# together from then on -- see ADR-0004 for the failure that motivated it.
references = ReferenceSet.from_compositions(spectra, compositions)

scores = cd.benchmark(references, folds=6,
                      options={'subset-average': {'draws': 200}})
for method, entry in scores.items():
    if entry['n']:
        print(f"{method:<16} helix rmse {entry['helix']['rmse']:.3f}  "
              f"bias {entry['helix']['bias']:+.3f}   (n={entry['n']})")
```

**Bias is worth as much as rmse.** A method that pulls every protein towards
the average composition of its reference set has an error that grows with how
unusual the protein is — which is exactly when the estimate matters.

## Measured: SP175 and SMP180, ten-fold held out, 190–240 nm

Not from this page — these are runs against the real reference sets, which
you fetch yourself (see below). Errors are in fraction units, so 0.14 is
fourteen percentage points.

**SP175** (71 soluble proteins):

| method | helix rmse | helix bias | sheet rmse | sheet bias |
|---|---|---|---|---|
| `all-references` | — | — | — | — refused, all 71 |
| `nearest-shapes` | 0.174 | +0.050 | 0.136 | −0.035 |
| `subset-average` | 0.164 | −0.077 | 0.118 | +0.059 |
| **`ridge`** | **0.143** | **−0.019** | **0.117** | +0.019 |

**SMP180** (128 soluble and membrane proteins):

| method | n answered | helix rmse | helix bias | sheet rmse | sheet bias |
|---|---|---|---|---|---|
| `all-references` | 0 | — | — | — | — refused, all 128 |
| `nearest-shapes` | 128 | 0.160 | +0.048 | 0.137 | −0.036 |
| `subset-average` | 108 | 0.155 | −0.076 | 0.119 | +0.052 |
| `ridge` | 128 | 0.143 | −0.013 | 0.121 | +0.010 |
| **`selcon`** | **95** | **0.107** | **+0.001** | **0.103** | **−0.002** |

**`selcon` is the most accurate and the least biased** — and it declines to
answer for 26 % of the set, which is the trade. That refusal rate is a feature
rather than a shortfall: the alternative is a confident number from a fit that
failed its own consistency tests.

Its rmse is conditional on succeeding, so it needs the fair comparison. On the
**same 95 proteins** SELCON answered, `ridge` scores 0.137 against SELCON's
0.107 — so SELCON is genuinely better, not merely selective. The 33 it refused
are somewhat harder for `ridge` too (0.161), so there is a mild selection
effect on top.

Three things worth reading off these.

**`ridge` is the recommendation**, on both sets, on both rmse and bias. Not by
a wide margin, but consistently.

**`subset-average` has a real bias**: −0.077 on helix, +0.059 on sheet. It
pulls towards the reference set's own mean composition, which for SMP180 is
33 % helix. That is the regression-to-the-mean effect above, measured rather
than argued.

**`nearest-shapes` is biased the other way**, +0.05 on helix. It has no
mechanism for describing a protein as *less* helical than its nearest
neighbours.

**`all-references` never once succeeded.** It refused every fold, correctly:
50 numbers cannot determine 71 references, let alone 128.

### What the wavelength range costs

SMP180, ten-fold, helix:

| region | `ridge` rmse | `nearest-shapes` rmse |
|---|---|---|
| 180–240 nm | 0.126 | 0.127 |
| 190–240 nm | 0.143 | 0.147 |
| 197–240 nm | 0.136 | 0.133 |
| 205–240 nm | 0.145 | 0.141 |

Reaching to 180 nm is worth about 0.02 in rmse. Losing the range down to
205 nm costs surprisingly little — **the far-UV cut-off matters less than the
choice of method or the reference set**, which is not what we expected before
measuring it.

## Getting the reference sets

They are **not shipped**. SP175, SMP180 and IDP175 are published by the PCDDB
organisation as part of [DichroWebGit](https://github.com/pcddb/DichroWebGit)
under the **MIT licence**, which does permit redistribution — but the same data
taken from the PCDDB website carries no such grant, so the provenance matters
and is worth keeping with the files.

```
git clone https://github.com/pcddb/DichroWebGit
```

```python
from spectroscopy.library import load_dichroweb_basis

references = load_dichroweb_basis('DichroWebGit/Datasets/SMP180')
result = cd.estimate(my_spectrum, 'ridge', references, region=(190.0, 240.0))

print(references)                    # what it is and where it came from
print(references.info['licence'])    # and the terms it arrived under
print(result.quality['references'])  # which set this answer used
```

Cite [Lees *et al.* 2006](../references.md) for SP175, [Abdul-Gader *et al.* 2011](../references.md) for
SMP180, and [Miles *et al.* 2022](../references.md) for DichroWeb itself.

For a basis you measured or obtained elsewhere, `library.load_basis` takes a
small CSV manifest and any file format the library reads.

## A worked disagreement

AqpZ-W14A, an all-helical membrane channel, measured 196–280 nm. Its crystal
structure (1RC2) is **178 of 231 residues helical with no sheet at all**;
diluted by a 23-residue N-terminal tag the construct carries, that is
**70.1 % helix**.

| method (SMP180) | helix | sheet |
|---|---|---|
| `selcon` (`ignore_sum_rule=True`) | **0.657 ± 0.030** | 0.041 |
| `nearest-shapes` | 0.633 | 0.043 |
| `ridge` | 0.449 | 0.180 |
| `subset-average` | 0.406 | 0.206 |
| **crystal structure** | **0.701** | **0.000** |

`selcon` — the best method by held-out validation, chosen before AqpZ was
looked at — lands within 0.045 of the crystal structure, and its sheet content
is correctly near zero. But note `ignore_sum_rule=True`: no concentration was
supplied, so this is the amplitude-dependent method run without its amplitude,
and the number should be read with that in mind.

The methods disagree by 0.25 in helix — nearly twice the 0.14 rmse the
benchmark says to expect — and the best of them is 0.07 from the truth while
the benchmark's recommended one is 0.25 away.

**That disagreement is the result.** It says this protein is not well
described by the reference set, which is unsurprising: a detergent-solubilised
membrane channel, where absorption flattening suppresses the 208 nm band and
raises θ222/θ208 above 1. Reporting `ridge`'s 0.449 alone, with its excellent
residual, would have been a confident wrong number.

Run more than one method. When they agree, the answer is probably safe. When
they disagree by more than the benchmark's rmse, believe the disagreement.

## Membrane proteins look like each other

Asked which reference proteins AqpZ most resembles, the answer came back
dominated by membrane proteins — lactose permease, a mechanosensitive channel,
sensory rhodopsin II, a leucine transporter. That could be the method
recognising the fold, or it could be membrane spectra sharing a distortion.
It is testable, and it is the second.

In SMP180, 30 of 128 proteins are membrane proteins — 23 %. Taking each
protein in turn and asking what fraction of its five nearest shapes are
membrane proteins:

| query | membrane among its 5 nearest |
|---|---|
| membrane proteins | **0.45** |
| soluble proteins | 0.20 |
| chance | 0.23 |

**And it survives controlling for structure.** Within a single helix-content
band, where both classes are present:

| helix content | membrane query | soluble query |
|---|---|---|
| 0.0–0.2 | 0.43 | 0.19 |
| 0.4–0.6 | 0.44 | 0.20 |
| 0.6–1.0 | 0.56 | 0.48 |

At the same helix content, a membrane protein's spectrum resembles other
membrane proteins about twice as often as a soluble protein's does. So there
is a **class signature in membrane CD beyond secondary structure** —
absorption flattening, differential scattering, and the longer straighter
helices of a transmembrane bundle all point the same way.

Two consequences. Using SMP180 rather than SP175 for a membrane protein is
right, and not only because it contains more of them. And a shape-similarity
result on a membrane protein is partly recognising the *class*, so some of
what looks like structural agreement is agreement about being a membrane
protein — which is worth knowing before quoting it.

## Noise you measured, rather than noise you assumed

CD is almost always recorded as several accumulations, so the error at each
wavelength is a measurement. It is also strongly wavelength-dependent: on four
near-native AqpZ scans the standard error is **0.45 mdeg at 215 nm and 0.16 at
240**, because the photomultiplier is working far harder at the blue end.

```python
sigma = cd.uncertainty_from_replicates(my_repeat_scans)   # a SpectrumCollection
result = cd.estimate(spectrum, 'ridge', references,
                     sigma=sigma, resamples=200)
print(result.quality['uncertainty'])
```

`sigma` weights each wavelength by `1/sigma`, so a point known ten times as
well counts ten times as much — which an unweighted fit does not do.
`resamples` refits on the spectrum perturbed by its own measured noise and
reports the spread. **That is the honest error bar**: how far the answer moves
for noise the size you actually have, which no goodness-of-fit statistic can
tell you.
