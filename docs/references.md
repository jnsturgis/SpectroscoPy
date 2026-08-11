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

  A general review of protein infrared spectroscopy — Barth's later one — is
  also wanted for the amide I assignments; the 2000 paper listed below covers
  side chains, which is a different question.

- **Amino-acid side-chain absorption** — the reference behind
  `processing.ftir.ftir_sidechain`. Barth, A. (2000). The infrared absorption
  of amino acid side chains. *Progress in Biophysics and Molecular Biology*
  **74**, 141–173.
  [10.1016/S0079-6107(00)00021-3](https://doi.org/10.1016/S0079-6107(00)00021-3)

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

- **Non-negative matrix factorisation** — `decompose(method='nmf')`. Lee,
  D. D. & Seung, H. S. (1999). Learning the parts of objects by non-negative
  matrix factorization. *Nature* **401**, 788–791.
  [10.1038/44565](https://doi.org/10.1038/44565)

  The non-negativity constraint is the reason NMF rather than PCA is the
  default for mixtures: additive-only combinations are what make the recovered
  components readable as species rather than as contrasts.

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


## Circular dichroism — the methods compared in `processing.cd`

Each entry says which method in
[the CD methods guide](guide/cd-methods.md) it stands behind. The comparison
between them in this library is by held-out validation on the reference sets,
not by these papers' own reported accuracies, which were measured on different
sets under different conditions.

- ⚠ **CDSSTR — the subset idea behind `subset-average`.** Johnson, W. C.
  (1999). Analyzing protein circular dichroism spectra for accurate secondary
  structures. *Proteins* **35**, 307–312.

The comparison of CONTIN, SELCON and CDSSTR by Sreerama & Woody (2000) is
listed under *Secondary structure* above and stands behind all three of the
classical methods here.

- ⚠ **CONTIN — regularised inversion, the ancestor of `ridge`.** Provencher,
  S. W. & Glöckner, J. (1981). Estimation of globular protein secondary
  structure from circular dichroism. *Biochemistry* **20**, 33–37.

- ⚠ **SP175, the soluble-protein reference set.** Lees, J. G., Miles, A. J.,
  Wien, F. & Wallace, B. A. (2006). A reference database for circular
  dichroism spectroscopy covering fold and secondary structure space.
  *Bioinformatics* **22**, 1955–1962.

- ⚠ **SMP180, adding membrane proteins.** Abdul-Gader, A., Miles, A. J. &
  Wallace, B. A. (2011). A reference dataset for the analyses of membrane
  protein secondary structures and transmembrane residues using circular
  dichroism spectroscopy. *Bioinformatics* **27**, 1630–1636.

- ⚠ **DichroWeb.** Miles, A. J., Ramalli, S. G. & Wallace, B. A. (2022).
  DichroWeb, a website for calculating protein secondary structure from
  circular dichroism spectroscopic data. *Protein Science* **31**, 37–46.

  The server whose reference sets are republished as
  [DichroWebGit](https://github.com/pcddb/DichroWebGit) under the MIT licence
  — which is what makes them usable here at all. The same data taken from the
  PCDDB website carries no redistribution grant.

- ⚠ **θ₂₂₂ helicity and its chain-length correction.** Chen, Y.-H., Yang, J. T.
  & Chau, K. H. (1974). Determination of the helix and β form of proteins in
  aqueous solution by circular dichroism. *Biochemistry* **13**, 3350–3359.

  Behind `structure.helix_from_theta222`, including the
  `(1 − 2.57/n)` term: a helix has two ends that make no hydrogen bonds, so a
  short chain signals less per residue than a long one.
