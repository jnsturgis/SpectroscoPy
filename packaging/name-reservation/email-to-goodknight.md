# Email 1 — to the current owner of the PyPI name

**To:** joey.goodknight@gmail.com
**From:** james.sturgis@univ-amu.fr  *(see note below on which address to send from)*
**Subject:** PyPI name `spectroscopy` — would you consider transferring it?

---

Dear Dr Goodknight,

I am writing about the PyPI project `spectroscopy`, which you published in 2017
alongside your work on vibrational dynamics.

I lead a small open-source project, currently called SpectroscoPy, that provides
a single data model for Raman, FTIR, UV-Vis and fluorescence spectra — the aim
being that a biologist with four instruments and four incompatible vendor tools
can do their analysis in one place, with the processing history recorded so the
result can be traced back. It is at
https://github.com/jnsturgis/SpectroscoPy and is approaching its first public
release.

The name `spectroscopy` on PyPI would be the natural one for it, and I would
like to ask whether you would be willing to transfer the project to me. I can
see that your package has not had a release since April 2017, and I appreciate
that it may still matter to you or be cited somewhere — if that is the case,
please just say so and I will use a different name; I would rather ask than
assume.

If you are willing, the transfer is done from the project's page on PyPI: under
"Manage" → "Collaborators", adding me as Owner is enough, and you can then remove
yourself. My PyPI username is [YOUR PYPI USERNAME]. I am happy to keep your
releases in place so that any existing install of version 0.10 keeps working,
and to credit the origin of the name in the project's documentation if you would
like that.

Thank you for considering it, and either way, thank you for leaving the work
public.

With best wishes,

James N. Sturgis
Laboratoire d'Ingénierie des Systèmes Macromoléculaires
Aix-Marseille Université / CNRS
james.sturgis@univ-amu.fr

---

## Before sending

1. **Fill in `[YOUR PYPI USERNAME]`.** If you do not have a PyPI account yet,
   register one first — the request is much easier to act on when he can add
   you immediately rather than having to write back and wait.
2. **Send from the university address, not Gmail.** The whole request rests on
   you being a real research group with a real project; `univ-amu.fr` in the
   From line does that work before he has read a word. It also matches the
   feedback address in the documentation.
3. **Keep the message.** If he does not reply, PyPI will want evidence that you
   tried to contact him before they will consider a PEP 541 request (see
   `pep541-request.md`). Sent-mail is that evidence — note the date you sent it.
4. **Wait about four weeks** before escalating. That is not a formal
   requirement, but the PEP 541 process asks for a reasonable effort at contact,
   and a month of silence on a nine-year-dormant project is reasonable by any
   reading.

## If he replies and says no

Use the fallback name, and do it without argument. He is under no obligation,
and the project loses nothing that matters: the displayed name stays
SpectroscoPy either way, and only the `pip install` line changes.
