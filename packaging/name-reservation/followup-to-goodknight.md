# Email 2 — the follow-up, five weeks on

**To:** joey.goodknight@gmail.com
**From:** james.sturgis@univ-amu.fr
**Subject:** PyPI `spectroscopy` — how is the account recovery going?

**Drafted:** 2026-09-08. **Sent:** 2026-09-08, as it stands (James).

Sequel to `email-to-goodknight.md`, which went on 2026-08-03 and was answered
agreeably on 2026-08-04. He consented to the transfer and then found he was
locked out of his PyPI account; his recovery request has not moved since. See
review §16 for the full record.

---

Dear Joey,

I am writing to ask how you are getting on with recovering your PyPI account,
rather than to hurry you along. When we exchanged messages at the beginning of
August you very kindly agreed to transfer the `spectroscopy` project, and then
found you were locked out. I imagine the request has been sitting in a queue
ever since; August is not the month for that sort of thing.

The reason I ask is that I have a decision to make in the middle of November. I
registered `pyspectroscopy` as a fallback before I first wrote to you, so
nothing is blocked and the library will be released either way. The only open
question is which name it goes out under, and since the software paper has to
cite the name, whichever one is in it is the name for good. So I have given
myself until Friday 13 November to choose.

If you have heard anything from PyPI support, or have a feeling for whether this
is weeks or months, that would be enough for me to decide. And if you would
rather not spend more time on it, please just say so; I will use the fallback and
that is a perfectly good outcome. You have already done the part that mattered by
agreeing.

Either way what I offered stands: your 0.10 release stays where it is so anything
installing it keeps working, and the origin of the name is credited in the
documentation.

With best wishes,

James N. Sturgis
Laboratoire d'Ingénierie des Systèmes Macromoléculaires
Aix-Marseille Université / CNRS
james.sturgis@univ-amu.fr

---

## Why it is written this way

**Ask, not escalate.** The delay is not his doing: `pypi/support#11748` is a
volunteer queue and it was filed in August. Nothing in the message asks him to
hurry, and the one thing that could read as pressure — the November date — is
stated as James's own constraint and immediately paired with a release from it.

**The November date is in it on purpose.** Leaving it out would be safer and
would also remove his reason to answer with anything useful. What is wanted from
him is not action but information: has support responded at all, and is this
weeks or months. That is what the decision on 13/11 turns on, and either answer
settles it.

**Why the name cannot stay open.** `paper.md` has to name the package, and a
paper naming a package that is later renamed is a citation that stops resolving
(roadmap §13.4). That is the whole reason an open wait acquired a deadline
rather than staying a hope — see roadmap §14.3.

**The two commitments are repeated** because they were made in return for his
agreement and should not depend on either party's memory: version 0.10 stays on
PyPI so anything pinned to it keeps working, and he is credited for the name in
`README.md`, on the documentation site and in `CITATION.cff`. They stand whether
or not the transfer ever completes.

**First name, not "Dear Dr Goodknight".** He replied warmly and quickly to the
first letter.

**The support ticket number is not quoted back at him.** He filed it; citing it
reads as chasing.

## After sending

This is the sent text, kept as the record. Sent mail is the evidence PEP 541
would want if this ever had to be escalated, though on current form it will not:
the owner consents, and the only obstacle is PyPI restoring his access.

If he answers with a timescale, or does not answer at all, the decision on
**Friday 13 November** is unchanged either way — ship as `spectroscopy` if the
transfer has landed, otherwise as `pyspectroscopy`.
