# Claiming `spectroscopy` — status and the request, if one is still needed

## ✅ Status: the owner has agreed (2026-08-04)

Joey Goodknight replied to the 2026-08-03 email and **is willing to hand over
the name**. He cannot act on it yet: he is locked out of his PyPI account and
has filed an account recovery request.

So the position is no longer "dormant project, silent owner". It is "owner
consents, and the transfer is waiting on PyPI restoring his access".

⚠️ **Do not file the request below as it was originally written.** Its central
argument was that the owner had not responded. That is now false, and filing
it would be both wrong and needlessly hostile to somebody who has been
helpful.

### Keep the email

His reply is the single most valuable document in this whole business.
Preserve it — with full headers — because it is what turns a contested claim
into a formality. If the fallback below is ever needed, that email *is* the
request.

## What happens next

**Path 1 — his recovery succeeds (expected).** No PyPI admin involvement is
needed at all. Once he is back in, he can do it himself:

> `https://pypi.org/manage/project/spectroscopy/collaboration/` → invite
> `jnsturgis` as **Owner**. He can then remove himself, or not; either way the
> name is usable.

That is the clean route and it is worth saying so to him explicitly, because
it is not obvious that a project can be handed over without involving support.

**Path 2 — his recovery stalls.** Then a support request is still needed, but
it is a *completely different and much easier* one than the draft below: not
"take this name from an absent owner" but "the owner agrees and cannot reach
his account". PyPI's volunteers handle that routinely and sympathetically.

Both his recovery request and any transfer request live in the same tracker,
<https://github.com/pypi/support>, which makes cross-referencing easy — quote
his recovery issue number if he will share it.

**When to chase.** Account recovery is usually days to a few weeks. If nothing
has moved by **2026-09-15**, ask him how it is going rather than escalating;
he is on our side and the delay is not his doing.

---

## Fallback text — *owner consents but cannot access the account*

Use this, not the section after it, if Path 2 becomes necessary.

> **Project:** `spectroscopy` — https://pypi.org/project/spectroscopy/
> **Requesting user:** `jnsturgis` — https://pypi.org/user/jnsturgis/
>
> I am asking to take over the `spectroscopy` project name **with the current
> owner's agreement**. Joey Goodknight, who uploaded the only release in 2017,
> confirmed by email on 2026-08-04 that he is willing to transfer it. He is
> unable to do so himself because he cannot log in to his PyPI account and has
> filed an account recovery request [issue number, if known].
>
> I am happy to forward his message with full headers, or to have him confirm
> here once his access is restored — whichever you prefer.
>
> The replacement package is SpectroscoPy,
> https://github.com/jnsturgis/SpectroscoPy: an actively maintained MPL-2.0
> library for working with measured spectroscopic data (Raman, FTIR, UV-Vis,
> fluorescence), with continuous integration across Python 3.10–3.13, a
> published documentation site and a tagged release. It is developed at the
> Laboratoire d'Ingénierie des Systèmes Macromoléculaires (Aix-Marseille
> Université / CNRS) and used for laboratory work.
>
> I have no wish to break anyone's install: I am content for version 0.10 to
> remain available and would publish well above it.

That is the whole thing. With the owner's consent, none of the dormancy or
download-statistics argument below is needed — and leading with it would
invite a debate nobody is having.

---

## Original text — *contested claim, only if consent is ever withdrawn*

**Kept for the record.** This was written on the assumption of an absent
owner. It is now historical, and would only be relevant if the situation
reversed entirely.

**Where it would go:** a new issue at <https://github.com/pypi/support>,
choosing the **"PEP 541 Request"** template. Realistic expectation: months.
That is why `pyspectroscopy` was reserved — nothing is blocked on the outcome.

### Project to be claimed

`spectroscopy`: https://pypi.org/project/spectroscopy/

### Requesting user

`jnsturgis`: https://pypi.org/user/jnsturgis/

### Reasons for the request

The project is abandoned, and I would like to use the name for an actively
developed package in a different area of the same field.

**The existing project has been dormant for nine years.** Its only release,
version 0.10, was uploaded on 2017-04-04. There have been no releases, and no
activity in its linked repository (https://github.com/jgoodknight/spectroscopy/),
since.

**It is essentially unused.** PyPI download statistics currently show 60
downloads in the last month, 6 in the last week and 0 in the last day — a
profile consistent with mirroring and automated crawlers rather than human
installs. I am not aware of any package on PyPI that depends on it.

**The owner has not responded to a direct request.** — *No longer true; see
the top of this file. He replied on 2026-08-04 and agreed.*

**The name does not describe the existing project especially well, and does
describe mine.** The current package is a simulation tool for quantum
vibrational dynamics ("Spectroscopy of systems with explicit vibrational degrees
of Freedom"). The package I maintain is a general-purpose library for working
with measured spectroscopic data — Raman, FTIR, UV-Vis and fluorescence — which
is what a user typing `pip install spectroscopy` is overwhelmingly likely to be
looking for.

**The replacement is real and maintained.** SpectroscoPy is at
https://github.com/jnsturgis/SpectroscoPy under MPL-2.0, with continuous
integration across Python 3.10–3.13 on Linux and Windows, a published
documentation site, and a tagged 0.1.0 release. It is developed at the
Laboratoire d'Ingénierie des Systèmes Macromoléculaires (Aix-Marseille
Université / CNRS) and is in use for laboratory work. I intend to maintain it
for the foreseeable future and to publish it under this name.

**On preserving the existing releases:** I have no wish to break anyone's
install. If the transfer is granted I am content for version 0.10 to remain
available, and I would publish under a version number well above it.

---

## Before opening anything

- [x] ~~Fill in the PyPI username~~ — `jnsturgis`. **Check this is right**; it
      was taken from the account-recovery filename, not confirmed directly.
- [x] ~~Fill in the date the email was sent~~ — 2026-08-03, replied 2026-08-04.
- [ ] Ask Joey for his account-recovery issue number, so a transfer request
      can point at it.
- [ ] Tell him about the collaboration page — he may not know he can do it
      himself once he is back in.
- [ ] If quoting the test count or CI matrix, take the current numbers. A
      stale figure costs credibility in a request whose argument is that the
      replacement is real and maintained.
