# PEP 541 request — only if the email goes unanswered

**Where it goes:** a new issue at <https://github.com/pypi/support>, choosing the
**"PEP 541 Request"** issue template. Do not open this until roughly four weeks
after the email in `email-to-goodknight.md`, and fill in the date you sent it.

**Realistic expectation:** these take months, and PyPI's volunteers work through
a queue. That is the reason for reserving `pyspectroscopy` in the meantime — the
project must not be blocked waiting on this, and the paper (roadmap §13) needs a
name that is settled before submission, whichever one it turns out to be.

---

## Project to be claimed

`spectroscopy`: https://pypi.org/project/spectroscopy/

## Your PyPI username

`[YOUR PYPI USERNAME]`: https://pypi.org/user/[YOUR PYPI USERNAME]/

## Reasons for the request

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

**The owner has not responded to a direct request.** I emailed
joey.goodknight@gmail.com — the address published in the project's own metadata —
on [DATE SENT], asking whether he would be willing to transfer the project. I
have had no reply in the [N] weeks since. I am happy to forward that message to
the PyPI admins if it is useful.

**The name does not describe the existing project especially well, and does
describe mine.** The current package is a simulation tool for quantum
vibrational dynamics ("Spectroscopy of systems with explicit vibrational degrees
of Freedom"). The package I maintain is a general-purpose library for working
with measured spectroscopic data — Raman, FTIR, UV-Vis and fluorescence — which
is what a user typing `pip install spectroscopy` is overwhelmingly likely to be
looking for.

**The replacement is real and maintained.** SpectroscoPy is at
https://github.com/jnsturgis/SpectroscoPy under MPL-2.0, with continuous
integration across Python 3.10–3.13 on Linux and Windows, 264 tests, a published
documentation site, and a tagged 0.1.0 release. It is developed at the
Laboratoire d'Ingénierie des Systèmes Macromoléculaires (Aix-Marseille
Université / CNRS) and is in use for laboratory work. I intend to maintain it
for the foreseeable future and to publish it under this name.

**On preserving the existing releases:** I have no wish to break anyone's
install. If the transfer is granted I am content for version 0.10 to remain
available, and I would publish under a version number well above it.

---

## Before opening the issue

- [ ] Fill in `[YOUR PYPI USERNAME]` in both places.
- [ ] Fill in `[DATE SENT]` and `[N]` weeks from the email you actually sent.
- [ ] Check the download figures again and update them — they move, and a stale
      number is the kind of thing that costs you credibility in a request whose
      whole argument is "this project is inactive".
      `curl -s https://pypistats.org/api/packages/spectroscopy/recent`
- [ ] Confirm the GitHub repository is still dormant before asserting it.
- [ ] Have the 0.1.0 release and the docs site publicly visible first. The
      strongest part of this request is that the replacement obviously exists.
