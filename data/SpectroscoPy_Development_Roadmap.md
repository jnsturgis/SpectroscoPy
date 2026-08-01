
---

## 15. Late-August capabilities, before the testers (added 2026-08-02, James)

**The emails go out late August** — earlier is pointless, nobody is at a bench in
the first three weeks. Two capabilities are wanted before they do:

1. **Secondary structure analysis**, from FTIR **and CD**.
2. **OPUS and `.spc` input.**

Both are defensible as pre-tester work rather than 1.1: a tester who cannot read
their own files tests nothing, and secondary structure is the question §14.3's
index already names as the largest gap.

### 15.1 The scope change nobody has said out loud yet

**CD is a fifth technique.** The header of this document says four. Adding it is
mostly additive to the data model — a `KNOWNSPECTYPES` entry, axis conventions,
readers — with one exception that is not additive at all, in §15.4.

Worth stating plainly: this lands three weeks before testers and three months
before an API freeze. Everything in §15 is *additive* except the unit question,
which is why that one is called out separately.

### 15.2 Neither format has a file to build against

Checked 2026-08-02: the repository holds 79 JCAMP, the `.dpt` replicates, some
CSV and text. **No OPUS binary, no `.spc`, no CD data anywhere on this machine.**

This collides with working agreement §14.1 — a feature earns its place with a
real example on real data, supplied before the feature is rolled out — and the
`.dpt` precedent says the rule is right. That reader was written against real
files and the format turned out not to be what it looked like (review §7.1); had
it been written against the specification alone it would have been wrong in a way
nothing would have caught.

**A reader written against a format specification and tested with files I
generated myself tests only my reading of the specification.** So: files first.

Wanted, and small numbers are enough:

| | What | Why that many |
|---|---|---|
| OPUS | 3–5 `.0`/`.1`/`.2`, ideally from more than one instrument | Block structure varies with acquisition settings; one file proves one path |
| `.spc` | 3–5, ideally one multi-subfile | The multi-file variant is where the format gets interesting, and it is the one that maps to `SpectrumCollection` |
| CD | 3–5 in whatever the instrument writes, plus one protein of known structure | The known one is the validation target |
| FTIR protein | One amide I spectrum whose answer is already known | Same reason |

### 15.3 What each piece actually needs

**OPUS binary.** Reverse-engineered rather than published: a header, a block
directory, then parameter and data blocks. `brukeropusreader` is abandoned since
2019 and `spc-spectra` since 2018 (comparison draft §10.2), so both are
read-them-ourselves jobs, consistent with the §5.6 dependency policy. Fits the
existing registry: a reader function, one `register_reader` call, no core change.

**`.spc`.** Genuinely documented, which makes it the easier of the two. Multiple
subfiles map onto the registry's existing "many spectra from one file" path
(review §10.2), so that machinery already exists.

**FTIR secondary structure.** The recipe is settled practice: water and vapour
subtraction, crop to 1600–1700, second derivative to find band positions,
constrained lineshape fit, assign bands to structure by frequency, integrate,
report fractions.

The part that matters for the freeze: **this needs `fit_peaks()` and
`FitResult`, which ADR-0001 §6.2 explicitly deferred.** Building secondary
structure pulls them forward. That is good news rather than bad — a return type
designed against a real analysis is exactly what §6.2 said was missing, and it is
far better to add it before 1.0 than after. It should be built as a general
fitter with secondary structure as its first caller, not as a protein-shaped
special case.

**CD secondary structure.** Two independent problems. The estimate itself is a
spectrum-times-basis-set least squares — the same shape as the NMF work already
in `processing.multivariate`, and `stability()` may generalise to validating it.
The reference basis sets (SP175, SMP180 and relatives) are published but their
redistribution terms need checking before any of them ship inside the package;
the safe default is that the user supplies a basis set and we ship none.

There is also a much cheaper thing worth having regardless: **helicity from
[θ]₂₂₂**, a one-line standard estimate. It is not a deconvolution and should not
pretend to be, but it answers "is my protein folded" for most people who ask.

### 15.4 ⚠️ CD breaks a stated assumption about units

`spectroscopy.units` converts with a table: every conversion is a function of
`y` alone. CD's units are not.

Millidegrees to mean residue ellipticity needs **concentration, path length and
the number of residues**. That is not a unit conversion in the sense the module
was built for — it is a calculation using sample metadata, and there is nowhere
in the current design for `to()` to get it.

ADR-0001 §2.2 named exactly this as the trigger for revisiting the native-table
decision: *"revisit when the unit space stops being closed — path lengths,
concentrations, per-unit-time intensities."* CD arrives with two of the three.

Three ways out, to be decided rather than drifted into:

1. **A separate method**, `to_mean_residue_ellipticity(concentration, path,
   residues)`, leaving `to()` alone. Honest, explicit, no design churn — and it
   makes the dependence on sample metadata visible at the call site.
2. **`to()` reads from metadata** when the conversion needs it, raising a
   pointed error when the metadata is absent. Convenient; hides a dependency;
   makes the same call succeed or fail depending on a dict.
3. **Adopt `pint`**, which is what §2.2 anticipated. Correct in the long run and
   the wrong size of job for August.

Recommendation: **(1) now, revisit (3) after 1.0.** It is the only one of the
three that cannot be wrong, because it adds a name rather than changing the
meaning of an existing one — and a wrong answer here is frozen in November.
