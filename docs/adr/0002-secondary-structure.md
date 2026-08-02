# ADR-0002 — Secondary structure: one output, several inputs

**Status:** Proposed, 2026-08-02. The shape is settled (James); the methods are
not all built.
**Supersedes nothing.** Depends on ADR-0001 §6.2 (`FitResult`), which was built
for this.

---

## 1. The requirement, as stated

> Different inputs FTIR/CD, same output Composition. Where there are different
> deconvolution methods in the literature the choice should be upfront.
> Potentially reference to a PDB structure.

Three things follow, and they are the whole design:

1. **`Composition` is the common currency.** FTIR, CD and a PDB structure are
   three ways of answering one question, so they return one type. A caller that
   has both an FTIR and a CD estimate should be able to compare them without
   caring how either was obtained.
2. **The method is a required, named argument.** Not a default that varies by
   technique, and not a heuristic that picks one. The literature disagrees, so
   the code should make the user say which disagreement they are on the side of.
3. **A PDB structure is a third estimator, not a special case.** It answers the
   same question with a different instrument — X-ray or NMR rather than a
   spectrometer — so it returns a `Composition` like the others, and it is what
   the spectroscopic estimates get validated against.

---

## 2. `Composition`

```python
@dataclass
class Composition:
    helix: float          # all helical: alpha, 3-10, pi
    sheet: float          # all extended: parallel and antiparallel
    turn: float
    other: float          # unordered, random coil, everything left

    detail: dict          # method-native categories, unmerged
    method: str           # 'amide-i-curve-fit', 'cdsstr', 'selcon3', 'dssp' ...
    technique: str        # 'ATR-FTIR', 'CD', 'PDB'
    source: str | None
    quality: dict         # method-specific fit diagnostics
```

**Four canonical categories, and a `detail` dict that keeps everything else.**
This is the load-bearing compromise, so it is worth being explicit about why.

The methods do not agree on categories. Amide I curve fitting yields bands
assigned to α-helix, β-sheet, turns and random, sometimes splitting antiparallel
β and aggregation. CD reference sets in the CDPro family report six — regular
and distorted helix, regular and distorted sheet, turns, unordered. DSSP reports
eight states. Forcing all of these into four loses real information; keeping all
of them means no two estimates can be compared, which defeats the requirement.

So: **four coarse categories that every method can produce honestly, plus the
method's own categories preserved in `detail`.** Comparison happens on the four;
nothing is thrown away.

The merge rules are fixed and documented, not per-method:

| Canonical | DSSP | CDPro-family | Amide I |
|---|---|---|---|
| `helix` | H, G, I | regular + distorted helix | 1650–1658 |
| `sheet` | E, B | regular + distorted sheet | 1620–1640, 1680–1695 |
| `turn` | T, S | turns | 1660–1680 |
| `other` | `-` | unordered | 1640–1650 |

**`quality` is not optional.** Every method reports something about how much to
trust the answer — an R² and residual for curve fitting, an NRMSD for CD
deconvolution, nothing at all for DSSP because it is not an estimate. A
`Composition` that cannot say how good it is invites being read as fact.

---

## 3. The estimators

```python
from spectroscopy.processing import structure

a = structure.from_ftir(spectrum, method='amide-i-curve-fit', ...)
b = structure.from_cd(spectrum, method='cdsstr', basis=...)
c = structure.from_pdb('1LYZ', method='dssp')

a.compare(c)        # -> per-category differences and an RMSD
```

**Method is required, with no default.** `from_cd(spectrum)` raises, listing the
methods it knows and what they assume. This is deliberate friction: choosing
CDSSTR over SELCON3 is a scientific decision, and a library that picks silently
has made it on the user's behalf and hidden it in a paper's methods section.

### 3.1 FTIR — `amide-i-curve-fit`

Water and vapour subtraction, crop to 1600–1700, second derivative for band
positions, constrained fit (`FitResult`), assign by frequency, integrate.

The uncertainty question from `FitResult` matters here: a weak shoulder between
two strong bands is barely determined, and its area — a structure fraction —
inherits that. `quality` carries the per-band position uncertainties, and a
fraction whose band is soft should be reported as soft.

### 3.2 CD — `cdsstr`, `selcon3`, `contin-ll`, plus `theta-222`

The CDPro-family methods are a reference set and a constrained least-squares
solve, differing in how they select and weight reference proteins. They are
buildable from their published descriptions; what needs checking before anything
ships is the **redistribution terms of the reference sets** (SP175, SMP180 and
relatives). Default position: **the package ships no reference data**, the user
supplies a basis, and a helper documents where the published sets come from.

`theta-222` — helicity from mean residue ellipticity at 222 nm — is a one-line
estimate of one number. It belongs here because it is what most people actually
want when they ask whether a protein is folded, but it fills only `helix` and
`other`, and must say so rather than reporting zeros for `sheet` and `turn`.
**A category a method cannot estimate is `None`, not zero.**

### 3.3 PDB — `dssp`

Two routes: parse the `HELIX`/`SHEET` records from the file header, which is
free and only as good as the depositor; or assign from geometry, which is
correct and more work. The header route is worth having first, clearly labelled
as `method='pdb-header'` rather than `'dssp'`, because the two are not the same
thing and a user comparing against "the crystal structure" deserves to know
which they got.

---

## 4. What this rules out

- **A single `secondary_structure(spectrum)` function that dispatches on
  technique.** It would have to choose a method, which §1.2 forbids.
- **Returning a bare dict.** A dict cannot carry `quality`, cannot say which
  method produced it, and cannot refuse to be compared with an incompatible one.
- **Merging CD's regular/distorted split silently.** It goes into `detail`,
  because the distinction is real and the four-category merge is a convenience
  for comparison, not a claim that the finer categories are noise.

---

## 5. Open

- **Reference-set licensing** gates §3.2 and needs answering before code that
  bundles data is written. Ships nothing until it is answered.
- **Does `stability()` generalise** from NMF to CD basis-set deconvolution? If
  it does, CD estimates get bootstrap uncertainty on the same machinery, which
  would be a genuine advantage over the standard tools — they report one number
  with no error bar.
- **`Composition.compare()` output shape** — per-category differences plus an
  RMSD is the obvious thing; whether it also needs a plot belongs with `viz`.
- **Where the module lives.** `processing.structure` follows the existing
  `processing.ftir` precedent, but the CD half is not FTIR-specific and the PDB
  half is not spectroscopy at all.
