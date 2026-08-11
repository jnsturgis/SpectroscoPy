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
| `ridge` | no | the penalty |

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
basis, compositions = [], []
for index in range(24):
    weights = rng.dirichlet([1.4, 1.0, 1.0])
    y = sum(w * pure[c] for w, c in zip(weights, pure))
    basis.append(spc.Spectrum(x, y + 0.3 * rng.normal(size=x.size),
                              technique='CD', name=f'ref{index}'))
    compositions.append(Composition(
        fractions=dict(zip(pure, weights)), method='known', technique='X-ray'))

scores = cd.benchmark(basis, compositions, folds=6,
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

| method | helix rmse | helix bias | sheet rmse | sheet bias |
|---|---|---|---|---|
| `all-references` | — | — | — | — refused, all 128 |
| `nearest-shapes` | 0.160 | +0.048 | 0.137 | −0.036 |
| `subset-average` | 0.155 | −0.076 | 0.119 | +0.052 |
| **`ridge`** | **0.143** | **−0.013** | 0.121 | +0.010 |

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

basis, compositions = load_dichroweb_basis('DichroWebGit/Datasets/SMP180')
result = cd.estimate(my_spectrum, 'ridge', basis, compositions,
                     region=(190.0, 240.0))
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
| `ridge` | 0.449 | 0.180 |
| `nearest-shapes` | 0.633 | 0.043 |
| `subset-average` | 0.406 | 0.206 |
| **crystal structure** | **0.701** | **0.000** |

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
