# ADR-0002 — Secondary structure: DSSP as the baseline vocabulary

**Status:** Proposed, revised 2026-08-02 (James: DSSP as the baseline, and CD is
decomposition into standard spectra). The shape is settled; the methods are not
built.
**Depends on** ADR-0001 §6.2 (`FitResult`), which was built for this.

---

## 1. The requirement

> Different inputs FTIR/CD, same output Composition. Where there are different
> deconvolution methods in the literature the choice should be upfront.
> Potentially reference to a PDB structure.
> Stick with as baseline DSSP and then how this is folded into method dependent
> categories.

The last sentence is the one that makes the rest work, and it is what this
revision is about.

---

## 2. Why DSSP is the right baseline

An earlier draft proposed four canonical categories (helix, sheet, turn, other)
with a merge table. That was backwards: it invented a vocabulary and then asked
each method to approximate it. DSSP is better for three reasons.

**It is a definition, not an estimate.** DSSP assigns from hydrogen bonding in a
solved structure. Every spectroscopic method is trying to estimate something;
DSSP states what that something *is*. A vocabulary should be owned by the method
that can be checked, not by the ones being checked.

**The CD reference sets are already DSSP.** SP175, SMP180 and the CDPro family
derive their category fractions by running DSSP over the reference proteins'
crystal structures and grouping the states. So a CD estimate is already reported
in a coarsened DSSP vocabulary, whether or not anyone says so. Choosing DSSP as
the baseline is therefore not an imposition on CD — it is naming what CD already
does.

**It fixes the comparison problem at the root.** If every method's categories are
declared as *sets of DSSP states*, then comparing two methods is a well-defined
operation on set partitions rather than a table of judgement calls.

### 2.1 The vocabulary

| Code | State |
|---|---|
| `H` | α-helix |
| `G` | 3₁₀-helix |
| `I` | π-helix |
| `E` | extended strand, in a β-sheet |
| `B` | isolated β-bridge |
| `T` | hydrogen-bonded turn |
| `S` | bend |
| `-` | none of the above |

---

## 3. `Composition`

```python
@dataclass(frozen=True)
class Category:
    name: str                 # 'alpha-helix', 'regular-helix', 'aggregated'
    states: frozenset[str]    # DSSP states it covers; empty if none apply
    note: str | None = None   # why, when the mapping is not exact (see §5)

@dataclass
class Composition:
    fractions: dict[Category, float | None]
    method: str               # 'amide-i-curve-fit', 'cdsstr', 'dssp' ...
    technique: str            # 'ATR-FTIR', 'CD', 'PDB'
    quality: dict             # R^2, NRMSD, per-band uncertainty ...
    source: str | None = None
```

**A category is a name plus the DSSP states it claims.** That is the whole
mechanism. `helix` is `{H, G, I}`; DSSP's own output is eight categories of one
state each; a method that resolves nothing beyond helix-and-not-helix has two
categories, `{H,G,I}` and the rest.

**A fraction a method cannot estimate is `None`, never zero.** Helicity from
[θ]₂₂₂ estimates helix and nothing else; reporting `sheet: 0.0` would be a claim
that the protein has no sheet, which it never made.

**`quality` is required.** Every estimate says how much to trust it — R² and
per-band position uncertainty for curve fitting (ADR-0001's `FitResult` carries
both), NRMSD for CD decomposition, nothing for DSSP because it is not an
estimate. A `Composition` that cannot say how good it is invites being read as
fact.

---

## 4. Folding: how each method maps onto DSSP

| Method | Category | DSSP states |
|---|---|---|
| `dssp` | each state | one each — the identity mapping |
| `amide-i-curve-fit` | helix (1650–1658) | `{H, G, I}` |
| | sheet (1620–1640, 1680–1695) | `{E, B}` |
| | turn (1660–1680) | `{T}` |
| | other (1640–1650) | `{S, -}` |
| `cdsstr`, `selcon3`, `contin-ll` | regular helix | `{H}` — see §5 |
| | distorted helix | `{H}` — see §5 |
| | regular sheet | `{E}` — see §5 |
| | distorted sheet | `{E}` — see §5 |
| | turns | `{T}` |
| | unordered | `{S, -, B, I, G}` |
| `theta-222` | helix | `{H, G, I}` |
| | *(everything else)* | `None` |

Two things this makes visible that a merge table hides. Amide I's "other" and
CD's "unordered" are **not the same set** — CD's absorbs `B`, `I` and `G`, and
amide I puts `G` and `I` with the helix. Comparing an amide I helix fraction with
a CD helix fraction is therefore not quite comparing like with like, and the size
of the discrepancy is exactly the 3₁₀ and π content. That is a real effect,
usually small, and it has been invisible in every comparison either of us has
read.

---

## 5. Where the mapping breaks, and what to do about it

Two method categories are **not** unions of DSSP states, and pretending otherwise
would be the one dishonest move available here.

**Regular versus distorted (CD).** The CDPro-family sets split helix and sheet by
*position within the element* — the terminal residues of each helix are counted
as distorted. That is a subdivision of `H` on a criterion DSSP does not encode.
Both map to `{H}`, with `note` recording that the split is positional and not a
DSSP distinction. Consequence, and it is the right one: **`regular-helix` and
`distorted-helix` can only be compared with anything else after being summed.**

**Aggregation (FTIR).** The 1620 cm⁻¹ and 1690 cm⁻¹ pair from intermolecular
β-sheet has **no DSSP state at all** — it is a quaternary feature, and DSSP runs
on one chain. Its category carries `states=frozenset()` and is excluded from any
comparison, rather than being quietly folded into sheet.

For biofilm and membrane-protein work this is not an edge case; it is often the
band people care about most. A vocabulary that cannot say "this is not a
secondary structure" would force it into `sheet` and report aggregation as
structure.

---

## 6. Comparing two compositions

`a.compare(b)` projects both onto the **coarsest partition both can express**,
then reports per-category differences and an RMSD over them.

Concretely: an amide I estimate (`{H,G,I}` helix) against a CDSSTR estimate
(`{H}` regular + `{H}` distorted, `{S,-,B,I,G}` unordered) coarsens to
`{H,G,I,S,-,B}` versus... nothing usable, because CD's unordered has swallowed
`G` and `I`. The comparison that *is* defined is the union — helix+other against
helix+unordered — and `compare` should say so rather than silently pairing the
helix numbers.

**This is the payoff of the whole design.** The naive comparison — line up the
helix numbers — is what everyone does, and it is wrong by the 3₁₀/π content in
one direction and by the regular/distorted split in the other. Making categories
sets of DSSP states means the code can compute what is comparable instead of
assuming it.

Categories with no DSSP states (§5) and fractions of `None` are excluded from the
projection and reported separately, never treated as zero.

---

## 7. The estimators

```python
from spectroscopy.processing import structure

a = structure.from_ftir(spectrum, method='amide-i-curve-fit', ...)
b = structure.from_cd(spectrum, method='cdsstr', basis=...)
c = structure.from_pdb('1LYZ', method='dssp')

a.compare(c)
```

**Method is a required argument with no default.** `from_cd(spectrum)` raises and
lists what it knows. Choosing CDSSTR over SELCON3 is a scientific decision; a
library that picks silently has made it for the user and hidden it in a methods
section.

### 7.1 FTIR — `amide-i-curve-fit`

Water and vapour subtraction, crop to 1600–1700, second derivative for band
positions, then the fit. ADR-0001's `FitResult` does the work, and two of its
defaults exist because of this caller: components are pseudo-Voigt, because a
fixed lineshape can be wrong by tens of percentage points at R² > 0.97; and the
fit is weighted between A and d²A/dν², because the derivative annihilates the
residual background that water subtraction always leaves and that otherwise ends
up inside the band areas.

`quality` carries R² and the per-band position uncertainties, and a fraction
whose band is poorly determined must be reported as such — a weak shoulder
between two strong bands is barely constrained by the data, whatever the fit says.

### 7.2 CD — decomposition into standard spectra

The standard approach, and James's: express the measured spectrum as a
combination of standard spectra, then read the composition off the combination.
Two kinds of standard, which the method name must distinguish:

- **Basis spectra of pure structures.** The composition *is* the coefficient
  vector. Simple, and only as good as the basis.
- **Reference protein spectra with known structures** — SELCON3, CDSSTR,
  CONTIN-LL. Fit the unknown as a combination of reference proteins, then take
  the same combination of their DSSP-derived compositions. The methods differ in
  how they select and weight the reference set, which is precisely why §7's rule
  makes the choice explicit.

The second is where the field is, and it is also why §2's argument holds: those
reference compositions came from DSSP in the first place.

#### Prior art: this has been done, in 2025

Hoffmann, Jones and Rodger (*QRB Discovery* **6**, 2025,
[10.1017/qrd.2025.4](https://doi.org/10.1017/qrd.2025.4)) determine secondary
structure from IR and CD **independently and integrated**, using SELCON — which
is this ADR's requirement, published, with a year's head start. Three things it
settles that were open here:

1. **They publish a Python SELCON3.** So §7.2 is reuse-or-improve rather than
   reimplementation, and the honest thing is to compare against theirs rather
   than to write a third SELCON and assume it agrees.
2. **The scaling for combining the two techniques is specified**: CD as Δε per
   amino-acid residue, IR amide I normalised to a maximum absorbance of 15,
   which they found gives the best general performance. That is an empirical
   constant somebody had to determine, and it would have cost weeks to
   rediscover.
3. **Combining gains only about 2 %** in helix and sheet.

The third finding is the one that should change the design, and it is good news
rather than bad. If the point of using two techniques were a better number, 2 %
would not justify the work. What they actually report is that combining
**identifies anomalously large errors** — IR alone goes badly wrong on
high-helix proteins like haemoglobin, CD alone on high-sheet proteins.

So `Composition.compare()` (§6) is not a convenience for people who happen to
have run both. **It is the clinical use of having two techniques**, and the
disagreement is the signal, not an inconvenience. That argues for `compare`
reporting *where* two estimates diverge in terms a user can act on, rather than
returning a single distance — and for the guide saying plainly that agreement
between IR and CD is weak evidence that both are right, while disagreement is
strong evidence that one is wrong.

It also confirms §15.4 of the roadmap: Δε per residue is the unit their method
needs, and getting there from millidegrees requires concentration, path length
and residue count — which is exactly the conversion `units.py` cannot currently
express.

**No reference data ships until its redistribution terms are checked** (§9). The
user supplies a basis; a helper documents where the published sets come from.

`theta-222` — helicity from mean residue ellipticity at 222 nm — is not a
decomposition and must not be presented as one. It fills one category and leaves
the rest `None`.

### 7.3 PDB — `dssp` and `pdb-header`

Two routes, and they are **different methods with different names**, because a
user comparing against "the crystal structure" deserves to know which they got:

- `pdb-header` parses the depositor's `HELIX`/`SHEET` records. Free, and only as
  good as whoever deposited it.
- `dssp` assigns from geometry. Correct, more work, and the only one that
  produces the full eight states.

---

## 8. Rejected

| Alternative | Why not |
|---|---|
| Four fixed canonical categories | §2 — invents a vocabulary, then makes every method approximate it, and hides which comparisons are ill-defined |
| Letting each method keep its own categories | Then no two estimates can be compared, which defeats the requirement |
| Silently summing regular + distorted | §5 — it is a positional split, not a structural one; summing is right but must be visible |
| Folding aggregation into sheet | §5 — it is not secondary structure, and for aggregating samples it is the main result |
| Zero for a category a method cannot estimate | A claim the method never made |
| A single `secondary_structure(spectrum)` | Would have to choose a method, which §7 forbids |

---

## 9. Open

- **Reference-set licensing** gates §7.2 and must be answered before any code
  that bundles data is written.
- **Which DSSP.** A dependency (`mdtraj`, `biotite`, the `mkdssp` binary) or our
  own implementation. Given §5.6's dependency policy and that this is the
  baseline everything else is defined against, it deserves its own decision.
- **Does `stability()` generalise** from NMF to CD basis decomposition? If it
  does, CD estimates get bootstrap uncertainty on machinery that already exists —
  and the standard tools report one number with no error bar at all.
- **Where the module lives.** `processing.structure` follows the
  `processing.ftir` precedent, but the CD half is not FTIR and the PDB half is
  not spectroscopy.
