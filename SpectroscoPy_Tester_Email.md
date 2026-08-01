# Emails to testers

Two emails for the two audiences of the tester brief. **Both are ready to
send** — the documentation is published, at
<https://jnsturgis.github.io/SpectroscoPy/>. Personalise the bracketed bits; the
rest is written to be sent as it stands.

- **Email A** — GitHub-registered colleagues. The ask is a bounded 45 minutes
  and a public issue.
- **Email B** — scientists with their own spectra. The ask is one real analysis,
  reported by mail, with a named nerdy colleague offered as a way round the
  parts they would rather not do.

Send both **individually, not as a group mail with everyone in the To line**.
Four people who each got a personal note reply at a different rate from four
people who can see that three others were asked the same thing.

---

# Email A — for the GitHub-registered

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
Documentation: https://jnsturgis.github.io/SpectroscoPy/

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

**And one thing worth more than your own 45 minutes:** if you know someone who
actually measures spectra — a colleague, a student, anyone who owns a
spectrometer and does not enjoy the software that came with it — please point
them at it. They do not need to know Python or want to learn it; the
documentation starts from "you can copy text into a terminal" and no further.

If they hit something and would rather not write it up, they can tell you and
you can pass it on — that route works fine, and I would rather have the finding
second-hand than not at all. Their report is worth more to me than either of
ours, because they will hit the things you and I have both stopped being able to
see.

If you can get to it before **[EARLY SEPTEMBER DATE]** that would fit the
timetable; after that I am into freezing things and it is too late for your
report to change anything, which would be a waste of your 45 minutes.

Thank you — and the coffee offer is serious.

James

---

# Email B — for the scientists

For Candice, Chloé, the Letitia phage group, and anyone else with their own
spectra. Different ask, different tone: they are not testing software as a
favour, they are being asked whether a tool is any good for their actual work,
which is a question they are better placed to answer than anyone.

Fill in **[NERD]** with whichever technically confident colleague *they* already
get on with — the point of naming a specific person is that "ask someone in
computing" is a wall, and "ask [NERD], who is expecting it" is a door. Check
with the nerd first.

**Subject:** Would you try this on your own spectra?

Dear [NAME],

You know the spectrum-handling code I have been building. It now reads [THEIR
FORMAT], handles [THEIR TECHNIQUE], and keeps a record of every processing step
so you can still tell in six months how you got a figure.

I would like to ask you to use it for **one real piece of analysis** — something
you would otherwise do in [THEIR CURRENT TOOL] — and then tell me where it
fought you.

The reason for asking now, rather than when it is finished: in November I have
to freeze the interface, which means committing to not changing what things are
called or how they are used. After that, changing them breaks other people's
work. So this is the last point at which "that is the wrong way round" from
someone who actually does the measurements can still change anything. Once it is
frozen, it is frozen for years.

**What I am asking for**, and it should be about an hour with your own data:

1. Get it installed — or say the word and I will do it on your machine.
2. Do one analysis you would normally do. Your data, your usual workflow.
3. Tell me where you got stuck, what you expected something to be called that it
   is not, and anything that took more steps than it obviously should.

**How to tell me:** every page of the documentation has *Report a bug*,
*Request a feature* and *Ask a question* links at the bottom. They open an email
to me with the page and version already filled in, so you just write what
happened in your own words. You do not need to know what caused it. "It did
something strange when I did X" is a perfectly good report — working out why is
my job, not yours.

**If the installing or the writing-up is the off-putting part:** [NERD] has
offered to help, and knows the project. They can get it running on your machine,
sit with you for the first analysis if that is useful, and pass anything you
find on to me in the form I need. Nothing is lost by going through them — I would
rather have your findings second-hand than not have them.

**What is most useful to me is what annoyed you.** Not the parts that worked.
If something felt clumsy, or you had to look up how to do a thing you do every
week, that is exactly the report I want. And if the honest answer turns out to be
that it is not worth your time compared with what you already use — that is the
single most useful thing you could tell me, and I would rather hear it now than
after I have published a paper saying otherwise.

If you can get to it before **[EARLY SEPTEMBER DATE]**, it can still change the
design.

[Optional, and true: I will acknowledge you in the paper — the people who shape
what the thing looks like deserve to be named in it.]

Thank you,

James

---

## The documentation is published — resolved 2026-08-02

<https://jnsturgis.github.io/SpectroscoPy/>, redeployed from `main` on every
push. This section previously said the emails could not go out because there was
no URL: CI built the site and uploaded it as a workflow artifact reachable only
from the Actions tab.

Kept as a note because the reasoning still applies to anyone tempted to point
readers at the markdown in the repository instead. The pages execute their code
at build time, so the figures and outputs that make Getting Started
comprehensible exist only in the built site; the source shows code with nothing
coming out of it.

## Who to send them to

**Email A:** four to six people is plenty (tester brief §1). Aim for a mix — at
least one who has never touched spectroscopy, at least one on Windows, and at
least one who will enjoy step 3 far too much.

**Email B:** the three from review §9.5, who cover the four techniques between
them. Pair each with a [NERD] they already know, and warn the nerd first — being
volunteered by email is a poor way to find out.

The pairing does two jobs at once. It removes the installation and write-up
barrier for people whose interest is their spectra rather than their terminal,
and it routes their findings towards someone who can file them publicly. That is
the proxy relay of review §14.5, with a name attached instead of a good
intention — and it is honest, because the report really did come from the
scientist and the issue says so.
