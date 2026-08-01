# Two Parallel Publication Workflows: JOSS + Domain Application

**Purpose of this document:** two workflows, run in parallel, producing two
non-overlapping outputs from one codebase — a JOSS software paper and a domain
application paper. The separation is structural, not just a matter of tone:
each paper is scoped to a different question, so there is nothing to
double-count.

- **JOSS paper** answers: *is this well-designed, well-engineered research
  software?* Content = architecture, design decisions, functionality.
  Explicitly must **not** contain new research results.
- **Domain paper** answers: *what did we learn about a biological system?*
  Content = a scientific question, methods, results, discussion. The software
  is a cited tool, not the subject.

If each paper stays inside its own question, there is no basis for a
double-publication concern — this is the standard, well-established pattern
(most software papers work this way), not a workaround.

---

## 0. Shared prerequisite: the joint-decomposition capability

Both workflows depend on one piece of new core functionality that doesn't
exist yet: **joint NMF decomposition across two spectral techniques with
different axes** (FTIR wavenumber + fluorescence wavelength/emission). This
needs its own small spike *before* either paper can build on it, following
the same spike → validate → harden approach as the rest of the roadmap.

**Design questions to resolve first (own ADR — e.g. `docs/adr/0002-joint-decomposition.md`):**
- How are two `SpectrumCollection` objects (one per technique) paired per
  sample — a new `PairedCollection`/`MultiTechniqueCollection` object, or a
  lighter join-by-sample-ID operation at call time?
- Concatenation strategy before NMF: simple horizontal stacking of
  (suitably normalized) intensity vectors per sample, or a weighted/block
  scheme so one technique's larger dynamic range doesn't dominate the
  factorization?
- Output shape: does `decompose(...).contributions()` extend unchanged, or
  does joint decomposition need its own return type (e.g. per-technique
  loading vectors alongside the shared per-sample contributions)?
- Validation: `stability(...)` already exists for component-count validation
  on single-technique series (per your docs) — confirm it generalizes to the
  joint case, or note explicitly if it needs adapting.

**This belongs in the core library, described generically** (the algorithm,
the design trade-off, why concatenation-then-NMF rather than some other joint
method) — that generic description is what goes in the JOSS paper. The
specific biological result of *applying* it is what's reserved for the domain
paper. This is the cleanest possible instance of the firewall in §3.

---

## 1. Workflow A — JOSS submission

**Governing constraint (recap from policy check):** at least 6 months of
*public* development history, tagged releases, and evidence of real external
engagement (issues/PRs from people other than you) before submission is
realistic. This is a calendar floor, not a writing-effort floor — start the
public clock now if the repo isn't already public.

| Step | Task | Notes |
|---|---|---|
| 1 | Confirm repo is public, license is OSI-approved, issue tracker is open without registration | Non-negotiable JOSS eligibility requirement |
| 2 | Tag releases regularly (e.g. each phase-completion from the core roadmap) | Gives JOSS reviewers, and future readers, a visible development timeline |
| 3 | Keep the tester feedback loop (§2 of the original roadmap) running in the open | External issues/PRs are exactly the "real engagement" signal JOSS now weighs |
| 4 | Build/finish the GUI to whatever MVP state you want represented in the paper | JOSS paper should describe what's shipped, not what's planned — keep the "planned" list out of paper claims |
| 5 | Build the joint-decomposition capability (§0) as one of the architectural features described | This is your strongest "design thinking" content — trade-offs weighed, why this approach |
| 6 | Draft `paper.md` (250–1000 words): summary, statement of need, comparison to existing tools, functionality overview | The comparison draft from earlier conversation is most of the "statement of need" section, already written |
| 7 | Reference the domain paper (preprint or published) as evidence of real-world reuse, if available by submission time | Satisfies JOSS's "research impact/significance" criterion directly — see §4 timing note |
| 8 | Submit once the 6-month public-history floor is met and the paper reflects only shipped functionality | Review itself is lightweight (GitHub checklist, weeks not months) once eligible |

**What must NOT appear in the JOSS paper:** any biological finding, dataset-specific
result, or interpretation from the domain application. Generic/synthetic
example data (or the library's own built-in tutorial datasets, e.g. the
`ethanol` example already in Getting Started) should illustrate functionality
instead — not the novel biological dataset reserved for the domain paper.

---

## 2. Workflow B — Domain application paper

**Novel contribution (the thing that makes this a new publication, not a
restatement of the earlier FTIR-only work):** a joint FTIR + fluorescence
analysis of the same/complementary sample set, using NMF-based joint
decomposition as the analytical method — something the prior, already-published
FTIR-only paper didn't do and couldn't have done without new fluorescence data.

| Step | Task | Notes |
|---|---|---|
| 1 | Define the biological question the *combination* answers that FTIR alone didn't | This is the paper's actual thesis — worth writing as one sentence before anything else, since it's what justifies the paper existing separately from the earlier FTIR work |
| 2 | Acquire/source the complementary fluorescence data on the same or matched samples | Real data generation step — flagged in your answer as still needed |
| 3 | Confirm joint-decomposition capability (§0) is built and validated on synthetic/toy data first | Don't debug the method on your only real biological dataset |
| 4 | Run the joint NMF analysis; validate component count via `stability(...)` | Reuses existing library validation tooling — a good demonstration the core design generalizes |
| 5 | Compare joint (FTIR+fluorescence) decomposition against FTIR-only decomposition on the same samples | This comparison *is* the evidence that combination adds value — likely your key result figure |
| 6 | Draft the paper: Introduction (biological question, why FTIR-only was insufficient), Methods (samples, acquisition, joint NMF — cite SpectroscoPy/JOSS paper as the tool, described at method-citation depth, not re-explained architecturally), Results (biological findings), Discussion | Software description here should be a citation + one sentence, not a re-description — that's the JOSS paper's job |
| 7 | Target journal: **Journal of Bacteriology** (ASM) — best topical fit for a cell-envelope/peptidoglycan structural finding, IF ~2.8, Q2 | Depends on a prior citation being established in Mol Microbiol or Biophys J first — owned by you and collaborators, not a Claude Code task |
| 8 | Consider a preprint (e.g. bioRxiv) once results are solid | Independently useful for the paper's own timeline, and doubles as citable "research impact" evidence for the JOSS submission if it lands first or concurrently |

---

## 3. The double-tap firewall — explicit checklist

Run this check before either paper is finalized:

- [ ] JOSS paper contains **zero** biological results, findings, or interpretation — only software description
- [ ] Domain paper's description of the software is a citation + brief method mention, not a re-explanation of its architecture or design trade-offs
- [ ] Figures are not shared between papers using the *same* dataset for a *different* purpose — JOSS uses generic/tutorial data, domain paper uses the real biological data
- [ ] Domain paper explicitly states what's new relative to the earlier FTIR-only publication (the one-sentence thesis from Workflow B, step 1)
- [ ] The joint-decomposition ADR (§0) is written from the software-design angle only — no biological interpretation embedded in it
- [ ] Author lists reflect actual contribution to each output — software-only contributors needn't appear on the domain paper, and vice versa, unless genuinely involved in both

---

## 4. Timing coordination note

Since you're running these in parallel rather than sequentially, one detail is
worth deciding deliberately rather than letting it fall out by accident: if
the domain paper (or its preprint) is available *before* the JOSS paper is
submitted, it becomes citable evidence of real research impact in the JOSS
submission — which is now an explicit criterion JOSS weighs. This isn't a
reason to rush the domain paper, but if there's flexibility in either
timeline, having the domain preprint land first (or simultaneously) is
strictly better for the JOSS submission than the reverse order.

---

## 5. What to hand Claude Code

Two separate briefs, corresponding to the two workflows — keep them as
separate tasks/sessions so the "no overlap" boundary is enforced structurally,
not just by instruction:

**Brief A (JOSS):** "Build and validate the joint-decomposition capability
(§0 ADR first). Prepare the repo for JOSS eligibility (license, public issue
tracker, tagged releases). Draft `paper.md` describing only shipped
functionality, using generic/tutorial datasets for any illustrative figures.
Do not reference specific biological findings."

**Brief B (domain paper):** "Once the joint-decomposition capability is
validated (do not begin on real data before then), run the joint FTIR+fluorescence
NMF analysis on [dataset], validate against `stability(...)`, and compare
against the FTIR-only decomposition. Draft results/figures for a domain
paper, targeted at Journal of Bacteriology. Cite SpectroscoPy via its JOSS
paper/repository; do not re-describe its internal architecture."

**Note (not a Claude Code task):** the Journal of Bacteriology submission
depends on a prior citation being established in Molecular Microbiology or
Biophysical Journal first. This is owned by you and your collaborators —
flagged here only so it's visible as a dependency on Workflow B's critical
path, not something to hand off.

---

## 6. Notes added 2026-08-02 — checked against the repository

### 6.1 JOSS eligibility: the calendar floor is already met

§1 assumes the public clock starts now. It does not — the repository has been
public since **2025-02-18**, seventeen months. Checked, not assumed:

| Requirement | Status |
|---|---|
| Public repository | ✅ since 2025-02-18 |
| OSI-approved licence | ✅ MPL-2.0, and GitHub's licence detection sees it |
| Open issue tracker, no registration wall | ✅ enabled |
| Tests, CI, documentation site | ✅ 264 tests, 8-job matrix, Sphinx site |
| Tagged releases | ✅ 0.1.0, more to come per §1 step 2 |
| Statement of need | ✅ README "Why" is exactly this, and is reusable in `paper.md` |
| Community guidelines | ⚠️ README has a short Contributing section; no `CONTRIBUTING.md`, no `CODE_OF_CONDUCT.md`. Reviewers check for these by name |
| `CITATION.cff` | ❌ absent. Not required, but it is what makes GitHub render a "Cite this repository" box, and it costs ten minutes |
| Repository description | ⚠️ still "Spectrum analysis and management programmes" — pre-dates the rewrite |
| **External engagement** | ❌ **0 stars, 0 forks, 0 issues, no outside PRs** |

The last row is the real floor, and it is the one a calendar cannot fix. JOSS
weighs evidence that people other than the author use the software; seventeen
months of public history with no external trace is weaker than six months with
three engaged testers. **This makes the alpha hand-off (review §9.5 — Candice,
Chloé, the Letitia phage group) a JOSS prerequisite, not just a design
validation step.**

**The mailto feedback of review §14.5 does not change** (James, 2026-08-02).
An earlier draft of this section suggested steering testers to the public issue
tracker instead. That was wrong, and for a reason worth keeping: Candice, Chloé
and the phage group **do not have GitHub accounts**, and the email links are how
real feedback actually arrives. Making the JOSS criterion drive the feedback
design would optimise the channel for a reviewer rather than for the scientists
it exists to serve — and would predictably yield less feedback, not more.

**Instead, the public trace is generated by proxy.** Email reports get relayed
into GitHub issues by James, crediting the reporter and saying plainly that the
report came in by mail. This satisfies the criterion honestly — a reviewer
seeing "reported by a tester by email" learns exactly what happened, which is
that outside scientists are using the software — while the people doing the
reporting never have to acquire an account to be heard.

It is also better for the project than the alternative, independently of JOSS:
relaying makes each report searchable, linkable from a commit, and visible to
the next person who hits the same thing. The cost is one paste per report, borne
by the person who already reads the mail.

### 6.2 The JOSS paper should not depend on §0

§0 makes joint decomposition a shared prerequisite for both workflows, and §1
step 5 puts it in the JOSS critical path. **This is a self-imposed risk worth
removing.** A JOSS paper describes shipped functionality: the library as it
stands — one data model across four techniques, provenance that survives
serialisation, a format registry, NMF with bootstrap stability testing — is
already a substantial scholarly effort by JOSS's own criterion. Joint
decomposition would be a nice paragraph in it; it is not what makes the paper
eligible.

Coupling them means a spike that turns out harder than expected (see §6.3 for
why it might) blocks *both* papers instead of one. Recommend: build §0 for the
domain paper, and describe it in the JOSS paper **if it exists by then**. The
firewall in §3 is unaffected either way.

### 6.3 Three things the §0 spike will hit

Answers to some of §0's open questions, and one trap it does not mention.

**The trap: NMF requires non-negative input, and baseline-corrected FTIR is
not.** ALS and rubberband correction both routinely leave small negative
excursions in the noise, and `sklearn`'s NMF raises on negative values rather
than quietly coping. So the spike has to make a decision — clip at zero, offset
the whole block, or exclude negative regions — and every one of those choices
biases the factorisation differently. Whatever is chosen must be a recorded
`ProcessingStep`, not something done silently inside `decompose`, or the
provenance claim quietly develops a hole exactly where the paper's method sits.

**Prior art the reviewers will raise.** Concatenating blocks and factorising is
**low-level data fusion** in the chemometrics literature, and the established
method for this problem is **multiset/multiblock MCR-ALS**, which SpectroChemPy
already implements (comparison draft §4). A joint-NMF paper that does not say
why NMF rather than MCR-ALS here will be asked, in both venues. There is a good
answer — NMF's parts-based decomposition and the bootstrap stability machinery
already in the library — but it needs to be written down deliberately rather
than discovered during review.

**Block scaling is the whole ball game.** §0 asks whether simple horizontal
stacking suffices. It does not: an FTIR block with ~1700 variables in absorbance
units and a fluorescence block with a few hundred in counts will let one
technique dominate the factorisation entirely through sheer variance. The
standard fix is to scale each block — by its Frobenius norm, or by
`sqrt(n_variables)` — so each contributes comparably before concatenation. This
is the single design decision most likely to change the biological conclusion,
which makes it the thing the ADR must justify most carefully, and a good
candidate for a sensitivity figure in the domain paper.

**On the other two design questions:** join by sample ID at call time rather
than introducing a `PairedCollection` type — the collection already groups by
filename convention, and a new object should be earned by friction the spike
actually finds. `stability(...)` should generalise, with one change: the
bootstrap must resample **paired** samples jointly across blocks, never
independently per block, or the pairing that the whole method rests on is
destroyed by the validation that is supposed to confirm it.

### 6.4 Reconciling with roadmap §13

Roadmap §13, written yesterday, describes a single application paper on
"multivariate analysis of fluorescence and FTIR of *Pseudomonas aeruginosa*
biofilms". This document describes Workflow B as a joint FTIR + fluorescence
analysis whose novelty is measured against "the earlier FTIR-only publication",
targeted at *Journal of Bacteriology* on a "cell-envelope/peptidoglycan
structural finding".

**Resolved (James, 2026-08-02): one study, the biofilm one.** Workflow B is an
enhanced rewrite of the earlier FTIR-only biofilm publication, with new data and
an added level of complexity — paired fluorescence and FTIR spectra of the same
samples, decomposed jointly. So §2 step 2 is one experiment: acquire the
fluorescence half of the pairs on the biofilm samples.

Decided at the same time: **the two workflows are decoupled** (per §6.2), and
**JOSS is submitted against shipped capabilities at 1.0.0, targeted early
November 2026**. Roadmap §14 is the schedule that follows from that, including
the five things found to be blocking an API freeze. The consequence worth
carrying back here: joint decomposition is no longer on the JOSS critical path
at all, so §1 step 5 lapses and §0 belongs entirely to Workflow B.

Everything else in §13 survives unchanged: the feature list there (library
lookup coupled to decomposition, Beer–Lambert against Bradford, OPUS binary,
`.spc`, A260/A280, line narrowing) is orthogonal to the joint-decomposition
work, and items 1–3 are still what the domain paper needs regardless of which
study it turns out to be.
