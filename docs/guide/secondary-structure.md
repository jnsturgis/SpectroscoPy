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

# Secondary structure

How much helix, how much sheet — estimated from a spectrum, reported in a
vocabulary that survives being compared with a different technique.

:::{admonition} Status
:class: warning

FTIR is implemented. Circular dichroism is designed and being built on a
branch; a structure read from a PDB file is planned. All three return the same
object, which is the point of
[ADR-0002](../adr/0002-secondary-structure.md).

The estimator has been checked against synthetic bands of known composition. It
has **not** yet been checked against a protein whose structure is published —
until it has, treat the numbers as arithmetic rather than as biology.
:::

## The recipe

```{code-cell} python
import numpy as np
import spectroscopy as spc
from spectroscopy.lineshapes import spec_comp
from spectroscopy.processing import structure

# A synthetic amide I envelope: four overlapping bands, known composition.
x = np.linspace(1580, 1720, 561)
bands = [(1630, 16, 0.50), (1645, 14, 0.25), (1653, 13, 0.85), (1670, 15, 0.30)]
y = sum(spec_comp(x, position, width, height, 0.6)
        for position, width, height in bands)

protein = spc.Spectrum(x, y, technique="ATR-FTIR", name="synthetic protein")

composition = structure.from_ftir(protein, method="amide-i-curve-fit")
print(composition)
```

Three things about that call.

**The method is named, and there is no default.** Which estimator produced a
number is part of the number. A library that picks one silently has made a
scientific choice on your behalf and hidden it in your methods section.

**The spectrum must already be baseline-corrected and water-subtracted.**
`from_ftir` does not attempt either, because both need judgement — which
reference, what scale factor — and doing them silently inside an estimator
hides the step that most affects the answer.

**The fit is weighted towards the second derivative.** Whatever background
survives your water subtraction is invisible to a second derivative but is
absorbed straight into band areas by the absorbance envelope. On a synthetic
band with a curved residual left under it, that weighting is the difference
between composition errors of 31 % and 5 %.

## Reading the quality, not just the answer

```{code-cell} python
q = composition.quality
print(f"R² = {q['r_squared']:.5f}   components = {q['components']}")
print("position uncertainties (cm⁻¹):", np.round(q['position_stderr'], 2))
```

A high R² does **not** mean the composition is right. Fitting bands of the
wrong shape gives R² above 0.97 with fractions wrong by tens of percentage
points, which is why the components are pseudo-Voigt by default. And a weak
component between two strong ones is barely determined by the data: its
position uncertainty will be several times its neighbours', and its area is a
structure fraction. Read the uncertainties.

## Aggregation is not a structure

Intermolecular β-sheet absorbs near 1615 cm⁻¹. It is not secondary structure —
it is a quaternary feature, and DSSP, which defines this vocabulary, runs on one
chain. So it gets a category with **no DSSP states**, is reported separately,
and never counts as sheet:

```{code-cell} python
aggregated = spc.Spectrum(
    x,
    spec_comp(x, 1615, 14, 0.6, 0.6) + spec_comp(x, 1653, 13, 0.6, 0.6),
    technique="ATR-FTIR", name="aggregated sample")

print(structure.from_ftir(aggregated, method="amide-i-curve-fit"))
```

For a sample that has aggregated, that band is usually the result you needed to
see. A vocabulary that could not say "this is not secondary structure" would
have reported a failed preparation as β-sheet.

## Comparing two techniques

This is the part that justifies the design. Categories are declared as **sets of
DSSP states**, so two estimates can be compared on what they can both actually
express rather than by lining up labels.

```{code-cell} python
from spectroscopy.processing.structure import Category, Composition

# What a CDSSTR-style estimate looks like: positional splits, and an
# "unordered" category that also claims 3-10 helix, π-helix and β-bridge.
cd = Composition(
    fractions={
        Category("regular-helix",   frozenset({"H"}), note="positional"): 0.30,
        Category("distorted-helix", frozenset({"H"}), note="positional"): 0.12,
        Category("regular-sheet",   frozenset({"E"})): 0.16,
        Category("distorted-sheet", frozenset({"E"})): 0.09,
        Category("turns",           frozenset({"T"})): 0.14,
        Category("unordered",       frozenset({"S", "-", "B", "I", "G"})): 0.19,
    },
    method="cdsstr", technique="CD", quality={"nrmsd": 0.03})

print(structure.from_ftir(protein, method="amide-i-curve-fit").compare(cd))
```

Two views, deliberately:

**By name** is the comparison everyone makes — helix against helix, with
`regular-` and `distorted-` summed into their parent. It is offered because
refusing to show it does not stop anyone making it.

**Strictly comparable** is the honest one. CD's *unordered* claims `G`, `I` and
`B`, which amide I files under helix and sheet — so those categories overlap,
and the only thing the two vocabularies separate cleanly is turns. That is a
disappointing answer and it is the true one.

The `!` lines are the useful part: they name the states responsible and the
direction of the bias, so "our FTIR helix runs high against CD" stops being a
mystery and becomes the 3₁₀ and π content.

:::{admonition} What agreement and disagreement mean
:class: tip

Combining infrared and CD improves the numbers themselves by only about 2 %
(Hoffmann, Jones & Rodger, 2025 — see [references](../references.md)). The
value of having both is not a better estimate; it is that **one technique alone
goes badly wrong on certain proteins** — infrared on high-helix proteins,
CD on high-sheet ones — and the disagreement is what tells you.

So: agreement is weak evidence that both are right. Disagreement is strong
evidence that one is wrong.
:::
