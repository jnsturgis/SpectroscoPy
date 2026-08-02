# References

Every method here implements something somebody published. This page says
which, so that a result can be attributed and a method checked against its
source rather than against this library's description of it.

Organised by what implements what, because that is the direction the question
usually runs: *"where does this number come from?"*

:::{admonition} Entries marked ⚠ are not yet verified
:class: warning

They were written from memory and their volume, page and DOI fields have not
been checked against the publisher. Do not cite them from here until the mark
is gone. The machine-readable source is `docs/references.bib`, which carries the
same marks and will seed the citation file for the software paper.
:::

## Secondary structure

The vocabulary itself, not merely a method — every category this library
reports is defined as a set of DSSP states. See
[ADR-0002](adr/0002-secondary-structure.md).

- **DSSP.** Kabsch, W. & Sander, C. (1983). Dictionary of protein secondary
  structure: pattern recognition of hydrogen-bonded and geometrical features.
  *Biopolymers* **22**, 2577–2637.
  [10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211)

- **CD deconvolution — CONTIN, SELCON, CDSSTR.** Sreerama, N. & Woody, R. W.
  (2000). Estimation of protein secondary structure from circular dichroism
  spectra: comparison of CONTIN, SELCON, and CDSSTR methods with an expanded
  reference set. *Analytical Biochemistry* **287**, 252–260.
  [PMID 11112271](https://pubmed.ncbi.nlm.nih.gov/11112271/)

- **FTIR and CD combined.** Hoffmann, S. V., Jones, N. C. & Rodger, A. (2025).
  Protein secondary structure determined from independent and integrated
  infra-red absorbance and circular dichroism data using the algorithm SELCON.
  *QRB Discovery* **6**.
  [10.1017/qrd.2025.4](https://doi.org/10.1017/qrd.2025.4)

  Prior art for exactly the design in [ADR-0002](adr/0002-secondary-structure.md)
  — one composition from two techniques, estimated both separately and jointly.
  They publish a Python SELCON3; they scale CD to Δε per residue and normalise
  the IR amide I band to a maximum absorbance of 15; and they report that
  combining the two gains only about 2 % in helix and sheet, its real value
  being to catch the cases where one technique alone is badly wrong.

- ⚠ **Amide I band assignment.** Byler, D. M. & Susi, H. (1986). Examination of
  the secondary structure of proteins by deconvolved FTIR spectra.
  *Biopolymers*.

- ⚠ **Protein infrared spectroscopy, review.** Barth, A. (2007). Infrared
  spectroscopy of proteins. *Biochimica et Biophysica Acta — Bioenergetics*.

## Processing

- ⚠ **Savitzky–Golay smoothing and derivatives** — `smooth`, `derivative`, and
  the derivative weighting in
  [`fit_components`](api.md). Savitzky, A. & Golay, M. J. E. (1964). Smoothing
  and differentiation of data by simplified least squares procedures.
  *Analytical Chemistry* **36**, 1627–1639.

- ⚠ **Asymmetric least squares baseline** — `baseline(method='als')`. Eilers,
  P. H. C. & Boelens, H. F. M. (2005). Baseline correction with asymmetric
  least squares smoothing. Leiden University Medical Centre.

## Multivariate analysis

- ⚠ **Non-negative matrix factorisation** — `decompose(method='nmf')`. Lee,
  D. D. & Seung, H. S. (1999). Learning the parts of objects by non-negative
  matrix factorization. *Nature* **401**, 788–791.

## Still to cite

Methods implemented here that have no reference on this page yet. The list is
public because it is a defect list, and a short one is easier to finish than a
forgotten one.

- Rubberband / convex-hull baseline
- Standard normal variate and vector normalisation
- Second-derivative peak detection, as a band-finding technique
- Principal component analysis and FastICA, as used in `processing.multivariate`
- Bootstrap resampling, behind `stability()`
- The Okabe–Ito colour-vision-deficiency-safe palette used by `viz`
- JCAMP-DX, the format specification the reader implements
- Nernst and Henderson–Hasselbalch, once `processing.titration` exists

## Applications

Papers that used this library. Empty, for now — the first entries will be the
software paper and the biofilm application described in the roadmap.

If you publish something that used SpectroscoPy, we would like to list it here:
the feedback links at the foot of any page reach a human.
