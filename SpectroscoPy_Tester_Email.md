# Email to the GitHub-registered testers

Ready to send once the documentation is published (see the note at the bottom —
this is the one thing that must happen first). Personalise the bracketed bits;
the rest is written to be sent as it stands.

Send them **individually, not as a group mail with everyone in the To line**.
Four people who each got a personal note reply at a different rate from four
people who can see that three others were asked the same thing.

---

**Subject:** A 45-minute favour, and a coffee in it for you

Dear [NAME],

I have a favour to ask that involves you breaking something, which I know is
the kind of thing you enjoy.

I have spent this year turning the pile of spectroscopy scripts I have been
accumulating into an actual Python library — one way of handling spectra
whatever instrument they came from, instead of four vendor tools that disagree.
It is at the point where it works, it is documented, and I am submitting it to a
peer-reviewed software journal in November.

Before that I have to freeze the API — commit to not changing the names and
arguments people depend on. Doing that on my own judgement alone is how you end
up stuck for years with a decision that made sense to exactly one person. So I
need people who are not me to try it and tell me what makes no sense.

You have a GitHub account, which makes you rarer among my colleagues than you
might think, and you have no particular investment in my being right.

**What I am asking for: about 45 minutes, and one GitHub issue at the end.** You
do not need to know any spectroscopy. Honestly, it is better if you don't —
half of what I need to find out is which bits only make sense to someone who
already knows the answer.

Repository: https://github.com/jnsturgis/SpectroscoPy
Documentation: [DOCS URL]

**1. Install it (10 minutes).** In a fresh virtual environment:

```
pip install git+https://github.com/jnsturgis/SpectroscoPy
python -c "import spectroscopy as spc; print(spc.__version__)"
```

Please note which operating system and Python version you are on. If anything
goes wrong at this step, stop there and tell me — installation problems are the
ones that cost real users the most and that I am least able to see from my own
machine, where everything is already installed.

**2. Work through Getting Started (20 minutes).** It uses data that ships with
the library, so you need none of your own. Do it *literally* — type what it
says, and please do not silently fix my mistakes in your head, which is what
competent people do automatically. Note every point where you had to guess what
I meant, or where what happened was not what you expected.

**3. Try to break one thing (10 minutes).** Whatever appeals: feed it a file of
the wrong kind, ask for a unit conversion that makes no physical sense, pass a
negative number where a positive one obviously belongs. What I want to know is
whether the error message tells you what to do next, or just fails at you.

**4. Open one issue (5 minutes).** One — on the single most annoying thing you
hit, not a list. Roughly this shape:

```
What I did:
What I expected:
What happened instead:
OS and Python version:
```

If nothing annoyed you, then an issue saying "installed on [OS], Python [x.y],
worked through Getting Started, no problems" is genuinely useful to me and takes
a minute. That is a real result, not a non-answer.

**If you have the appetite for slightly more:** if you spot a typo, or a
sentence that sent you the wrong way, send it as a pull request instead of an
issue. Editing the file on GitHub's website makes one, and it takes about as
long as writing the issue would have.

**One thing I would ask you not to do:** please do not open issues for problems
you did not actually hit, pad the list to be helpful, or say anything
complimentary you do not mean. A short honest report is worth considerably more
to me than a long generous one — this material goes in front of reviewers, and
more to the point, an issue written by someone who did not run the code finds
nothing, which defeats the entire object.

If you can get to it before **[EARLY SEPTEMBER DATE]** that would fit the
timetable; after that I am into freezing things and it is too late for your
report to change anything, which would be a waste of your 45 minutes.

Thank you — and the coffee offer is serious.

James

---

## Before sending: the documentation must be published

**There is currently no [DOCS URL] to fill in.** CI builds the documentation and
uploads it as a workflow artifact, which is only reachable by someone logged in
to GitHub who knows to look in the Actions tab. There is no GitHub Pages site
and no Read the Docs.

Do not send this pointing at the markdown files in the repository instead. The
pages execute their code at build time, so the figures and outputs that make
Getting Started comprehensible only exist in the built site; the raw source
shows the code with nothing coming out of it. That would generate a batch of
confused issues about the documentation being broken, which is noise rather than
signal, and a poor first impression from people doing you a favour.

Publishing is a small job — a `gh-pages` deployment step in the existing CI
workflow, since the build already produces the site. It is also required for
JOSS regardless: reviewers must be able to read the documentation without
cloning anything.

## Who to send it to

Four to six people is plenty (tester brief §1). Aim for a mix: at least one who
has never touched spectroscopy, at least one on Windows, and at least one who
will enjoy step 3 far too much.
