# What I need from testers

Two audiences, two different asks. Do not send the same message to both.

- **Scientist testers** (Candice, Chloé, the Letitia phage group) — people with
  real spectra and a real analysis to do. What is needed from them is *use*, and
  their reports arrive by email. Nothing below asks them to touch GitHub.
- **GitHub-registered testers** — colleagues with accounts who can be persuaded
  to do a bounded piece of work. What is needed from them is a **public,
  attributable trace** of someone other than the author using the software.

§1 explains why the second group matters. §2 is the brief to send them, ready to
copy. §3 is what James does with what comes back.

---

## 1. What the public trace is for, and what would make it worthless

JOSS review asks, in effect, whether anyone but the author uses this. The repo
has been public since February 2025 and shows **0 stars, 0 forks, 0 issues, no
outside pull requests** (dual-publication document §6.1). Seventeen months of
silence reads worse than three months with four engaged users.

**The trace has to be real.** Not because a reviewer would necessarily catch
otherwise, but because JOSS review is peer review, and manufacturing the
appearance of a user community is the same category of act as manufacturing a
figure. It is also self-defeating: the point of testers is to find the API
mistakes before 1.0 freezes them, and an issue written by someone who never ran
the code finds nothing.

So the rule for everything below: **ask people to actually use it, then ask them
to report what actually happened, in public.** If someone has nothing to report
because it all worked, the honest artefact is a one-line issue saying "installed
on Windows/3.12, ran the FTIR tutorial, worked" — which is genuinely useful
evidence and takes them a minute.

What a reviewer finds convincing, roughly in order:

| Artefact | Why it counts |
|---|---|
| An issue that gets closed by a commit | Shows the whole loop works: report → fix → release. The strongest single signal |
| A bug report with a real traceback and a real file | Obviously came from running the thing |
| A pull request, even a one-word typo fix | External *contribution*, which outranks external comment |
| "Works on my machine" install reports across OS/Python versions | Cheap, honest, and directly relevant to a package claiming to support Windows |
| Stars and forks | Weak, non-zero, free |
| Empty praise, or six issues opened in one hour by new accounts | Worse than nothing. Reads as manufactured, and invites a closer look at everything else |

**Volume needed: less than you think.** Four to six people, one or two
substantive items each, spread over the weeks between now and the November
submission. It is the spread and the substance that read as real, not the count.

### Second-degree testers are worth the most

Email A asks the GitHub group to pass it on to anyone they know who actually
measures spectra. This is the highest-value part of the whole exercise, for two
separate reasons.

**Scientifically:** someone at one remove has no investment in the project and no
idea what it is supposed to do. They hit the things a colleague has already
stopped being able to see — and they are the population the library claims to
serve, which nobody in the first ring is.

**For the paper:** a report from someone James has never met is qualitatively
different evidence from a report by a colleague doing him a favour. It is the
difference between "his friends tried it" and "it reached someone", and that
distinction is exactly what a JOSS reviewer is looking for when asking whether
the software has users.

The relay works the same way as in Email B, just in the other direction: the
non-nerdy user tells the nerdy one, who files it. Both directions of the same
pairing — the technically confident person is a translator, not a gatekeeper,
and the report should be in the user's own words with the relay stated plainly.

---

## 2. The brief — copy this to your GitHub-registered colleagues

> **The sendable version is `SpectroscoPy_Tester_Email.md`.** Use that one; it
> is the same four steps written as a personal email, with the personalisation
> points marked and a note on the one thing that must be done before sending
> (the documentation is not published anywhere yet). What follows is the
> earlier, terser draft, kept because it states the asks compactly.

> **Subject: 45 minutes of nerd favour, in exchange for a coffee**
>
> I'm submitting a Python spectroscopy library to a peer-reviewed software
> journal in November, and before I freeze the API I need people who are not me
> to actually try it and say what breaks. You have a GitHub account, which makes
> you rarer than you think.
>
> **The ask: about 45 minutes, and one GitHub issue at the end.** You do not
> need to know any spectroscopy — a good half of what I need to find out is
> what makes no sense to someone who doesn't.
>
> The repository: https://github.com/jnsturgis/SpectroscoPy
> The documentation: [link to the built docs]
>
> **Step 1 — install it (10 min).** In a fresh virtual environment:
>
> ```
> pip install git+https://github.com/jnsturgis/SpectroscoPy
> python -c "import spectroscopy as spc; print(spc.__version__)"
> ```
>
> Note your OS and Python version. If anything at all goes wrong here, stop and
> report that — installation problems are the ones that cost real users most and
> that I am least able to see from my own machine.
>
> **Step 2 — follow Getting Started (20 min).** Work through it with the data
> that ships with the library; you don't need any of your own. Do it *literally*
> — type what it says, don't fix my mistakes silently in your head. Keep a note
> of every point where you had to guess what was meant, or where what happened
> wasn't what you expected.
>
> **Step 3 — try to break one thing (10 min).** Pick whichever appeals: feed it
> a file of the wrong sort, ask for a unit conversion that makes no sense, pass a
> negative number where a positive one is expected. I want to know whether the
> error message tells you what to do next, or just fails.
>
> **Step 4 — open one issue (5 min).** One issue, on the single most annoying
> thing you hit. Not a list — one. If literally nothing annoyed you, then an
> issue saying "installed on <OS>, Python <version>, worked through Getting
> Started, no problems" is genuinely useful and I would like to have it.
>
> Template, for whatever you report:
>
> ```
> What I did:
> What I expected:
> What happened instead:
> OS and Python version:
> ```
>
> **What would help even more, if you have the appetite:** if you spot a typo or
> a sentence that misled you in the docs, send it as a pull request rather than
> an issue. Editing the file on GitHub's website makes one, and it takes about
> as long as writing the issue would.
>
> **Please don't:** open issues for things you haven't hit, pad the list to be
> helpful, or write anything complimentary you don't mean. A short honest report
> is worth more to me than a long generous one, and this goes in front of
> reviewers.

---

## 3. What James does with it

- **Answer every issue within a day or two.** Reviewers look at response times,
  and it is the cheapest credibility there is. It is also just good manners
  towards someone doing you a favour.
- **Fix at least a few, and close them with the commit that fixes them.** An
  issue closed by a linked commit is the artefact that shows the loop works.
- **Relay the email reports from the scientist testers into issues by proxy**,
  crediting the reporter and saying plainly that it came in by mail (review
  §14.5). "Reported by a tester by email" is honest and still counts.
- **Do not batch it.** Four issues in one evening reads differently from four
  issues over five weeks, and the second is what actually happens when people
  use software.
- **Chase in early September.** The API freeze needs feedback in hand by
  mid-September (roadmap §14.4); the failure mode is polite silence until
  October, by which point the choice is to slip 1.0 or freeze it blind.

## 4. Who to ask, from the inventory

The scientist testers cover the techniques — Candice on FTIR, Chloé on
fluorescence and UV-Vis, the Letitia phage group on JCAMP. The GitHub-registered
group does not need to know spectroscopy at all, and it is better if some of
them do not: the person who has never seen an ATR crystal is the one who finds
out that Getting Started assumes you know what a background spectrum is.
