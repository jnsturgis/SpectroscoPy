# ADR-0004 — Reference sets are collections, not a parallel hierarchy

**Status:** **Accepted and implemented, 2026-08-13.** Proposed 2026-08-12
(James, reviewing the CD branch: *"is this well structured or just ad hoc?"*).
Ad hoc — this records what the structure should be instead. Section 5 was left
open and was resolved on 2026-08-13; section 6 records what was built and the
one thing that was deliberately deferred.
**Depends on** ADR-0001 (the core data model) and the metadata schema added
for roadmap D2.
**Affects** `library`, `processing.cd`, `processing.unmix`,
`processing.structure.from_cd`, and `SpectrumCollection`.

---

## 1. What went wrong

Three separate pieces of work each needed "a set of reference spectra with
known properties", and each invented its own container:

```
UV-Vis unmixing   Library([Reference(name, spectrum, unit, source, ...), ...])
CD                (list[Spectrum], list[Composition])       <- parallel lists
everything else   SpectrumCollection([Spectrum, ...])
```

None of this was designed. It accreted, and the CD form is the worst of the
three because **two parallel lists can fall out of step**. Measured on SP175:

```
compositions one short : ValueError from numpy about core dimensions
compositions reversed  : no error. helix 0.464 -> 0.307
```

The second is the dangerous one — same length, wrong pairing, a plausible
answer, no complaint. This is the same failure class as the positional
matching in `with_parameters`, which was found and documented a week earlier,
and then reintroduced somewhere else.

Two further defects came out of examining it, and both are in
`SpectrumCollection` rather than in the CD code:

**Set-level data is stored per-item and can disagree.** `parameter_name` and
`parameter_unit` read like properties of a series, but every spectrum keeps
its own copy and the collection gathers them. When they disagree the gather
silently returns `None`:

```
two items, units 'C' and 'K':  collection.parameter_unit = None
```

**Non-JSON metadata degrades silently through `.spy`.** A `Category` stored in
`metadata` comes back as a `str` — the frozenset of DSSP states is gone, and
nothing says so.

---

## 2. Decision

### 2.1 A reference set **is** a `SpectrumCollection`

Not a sibling of one, not a wrapper around one. `SpectrumCollection` already
has the whole shape: ordered spectra, `to_matrix`, `select`, `filter`, `map`,
`resample`, indexing, iteration, and per-spectrum metadata.

```
ReferenceSet(SpectrumCollection)
    inherits  everything above
    adds      gathered views of each reference's known properties
              set-level provenance (see 2.3)
```

### 2.2 Store on the item, view from the set

**Per-item data lives in that item's own metadata. The collection provides
gathered views of it.** `collection.parameters` already works exactly this
way, and it is the reason the parameter cannot desynchronise from its
spectrum: there is only one list.

A reference set therefore gains `.compositions` and `.categories` as gathered
views, in the same style. **No parallel list exists to reverse**, so the
0.464 → 0.307 failure becomes unrepresentable rather than merely tested
against.

Call sites lose an argument:

```python
cd.estimate(spectrum, 'selcon', references)        # not (..., basis, compositions)
cd.benchmark(references, folds=10)
unmix(spectrum, references)
```

### 2.3 Per-item data has kinds, and one of them is new

What is in `metadata` today is not one sort of thing. The distinction that
matters is *what the datum is about*:

| kind | examples | belongs to |
|---|---|---|
| measurement conditions | path length, concentration, temperature | the item |
| position in a series | `parameter` | the item |
| **series description** | `parameter_name`, `parameter_unit` | **the set** |
| provenance | `opus_*`, `file_*` | the item; `reference_source` is the set's |
| **known truth** | composition, category | **the item — a new kind** |

**"Known truth" is not a measurement condition** and should not be filed as
one. "This protein is 46 % helix" is not a fact about how the spectrum was
recorded; it is an answer obtained by a different technique, and it is exactly
what makes a spectrum a *reference* rather than just a spectrum. It gets its
own group in `spectroscopy/metadata.py`, alongside sample conditions,
identification, acquisition and provenance.

The practical consequence of the group is a rule: **known-truth values must
name their source**. A composition from DSSP on a crystal structure and one
from a previous CD fit are not the same kind of evidence, and a reference set
built from the second is circular.

### 2.4 Set-level data becomes real

`SpectrumCollection` gains a small `info` dict for facts about the set rather
than about any item — surviving `crop`, `select`, `map` and slicing exactly as
`name` already does. For a reference set that is where `source`, `citation`,
`licence`, `accession` and the ordinate `unit` live.

Two reasons this is not cosmetic:

- **A citation condition needs somewhere to live.** Several published
  reference sets require citation as a condition of use. Smearing that across
  128 spectra as `reference_source` is not a record of a licence obligation.
- **Facts that are true of the set cannot then disagree.** The `'C'` vs `'K'`
  case above stops being possible rather than being detected.

`parameter_name` and `parameter_unit` move to `info`, with the gathered
versions kept as deprecated readers until 1.0.

### 2.5 Stored metadata is JSON-native

Because `.spy` serialises `metadata` as JSON and silently degrades anything
else, the **stored** form of a composition is a plain `{name: fraction}` dict
and the stored form of a category is its name. `Category` objects are
reconstructed from the set-level category declaration when the set is loaded.

This is not a workaround. It is what forces a single canonical
representation, and it is why the gathered views are computed rather than
cached.

---

## 3. What this does not change

`Reference` and `Library` stay, as a thin specialisation of `ReferenceSet` for
the UV-Vis case, so existing unmixing code keeps working. The UV-Vis case
needs no known-truth table at all — a spectrum in ε units *is* its own
property, and the fitted coefficient is the answer — which is why it fitted
into a different container in the first place.

The three cases are then one type with a table that is sometimes empty:

| | per-item known truth | what a fit returns |
|---|---|---|
| UV-Vis unmixing | none | concentrations |
| CD reference proteins | a composition | a composition |
| CD structural basis | a category = a one-hot composition | a composition |

The third is the second with a degenerate table, which is the argument for it
not having had its own code path.

---

## 4. Rejected

**A per-set table of properties (columns × items), DataFrame-style.** It reads
well and reintroduces exactly the defect being fixed: a table beside a list of
spectra is a parallel structure, and slicing, filtering or reordering the
collection must remember to slice the table identically. `select` and `map`
would each become an opportunity to desynchronise.

**Leaving CD as parallel lists and validating lengths.** A length check
catches the case that already raises and misses the case that does not. The
reversed-compositions failure has the right length.

**A new top-level `Library` hierarchy independent of `SpectrumCollection`.**
It would duplicate `to_matrix`, `select`, `resample` and the batch operations,
and produce two ways to hold spectra that behave almost but not quite alike.

---

## 5. Resolved — `.spy` holds a set as well as a spectrum (James, 2026-08-13)

**The problem.** There was no collection writer; spectra were saved one at a
time. So set-level data had nowhere to persist, and a `ReferenceSet` loaded
from `load_dichroweb_basis` could not be written back out with its licence and
citation attached.

**The decision: one format, not two.** `.spy` gains an optional `# collection`
block followed by repeated spectrum blocks. Not a second format beside it —
a set of spectra and a spectrum are the same kind of document, and a caller
should not have to know which they have before they can open it.

```
# spy format 1.0                    # spy format 1.0
# header                            # collection
{spectrum json}                     {"name": ..., "kind": ..., "info": {...}}
# data                              # spectrum
...                                 {spectrum json}
                                    # data
                                    ...
                                    # spectrum
                                    ...
```

`# header` and `# spectrum` are the same marker, so every file ever written
still reads and a single spectrum is still written byte-for-byte as before.
One function writes a spectrum block and both paths call it.

**No version bump**, and the reason is stronger than "nobody has files yet":
a bump would not have protected anyone. `_detect_version` accepts any `1.x`
and dispatches to the same reader, so a `1.1` collection handed to today's
code would have parsed as one spectrum with every block's numbers run
together. Only a major bump would have tripped it. What protects a caller is
the marker, which is checked.

**Set-level `info` is JSON**, as `metadata` is. The one value JSON has no form
for is a list of `Category`, and it is the one that must survive: a category is
a name *plus the DSSP states it claims*, and `'helix'` alone does not say
whether it covers 3-10 and pi. Categories are written as name + states.

**`kind` records the class**, so a saved `ReferenceSet` comes back one rather
than a plain collection with `.compositions` gone and the truth sitting unread
in each spectrum's metadata. Restored from a small explicit table — a file
should not be able to name an arbitrary class and have it constructed — and an
unrecognised name loads as a plain collection *with a warning*, so a file from
a later version is still readable for its spectra and its provenance.

The catalogue has exactly two entries today, and the rule for admission is
worth stating because otherwise it grows by habit: **a kind earns an entry when
it reads data the base class stores but does not interpret.** `ReferenceSet`
qualifies. A titration does not, and is deliberately not a class: it is a
`SpectrumCollection` with a parameter per spectrum and its name and unit in
`info`, both of which the base class already reads. Needing somewhere to put
*that* is what produced `info` in the first place.

**API.** `collection.save_as(path)` writes; `io.read_spectra(path)` reads, and
already returned a `SpectrumCollection`, so nothing new was added to the frozen
top-level surface. `read_spectrum` keeps its contract unchanged — one spectrum
out, an error if the file holds several. A collection file holding exactly one
spectrum therefore just reads, with no special rule for it: a set of one is a
set, and asking it for its single spectrum is a fair question.

---

## 6. What was built, 2026-08-13

`library.ReferenceSet(SpectrumCollection)`, with `.compositions` and
`.categories` as gathered views and `from_compositions` as **the one place**
spectra and structures are matched by position. `load_basis` and
`load_dichroweb_basis` return one; `cd.estimate`, `cd.benchmark`,
`structure.from_cd` and `structure.nearest_references` each lost an argument
and now refuse a bare list with a message saying where to build the set.

`SpectrumCollection` gained `info`, which survives `crop`, `select`, `map`,
`group_by` and slicing — and so does the subclass, since every operation now
derives through one helper rather than hard-coding `SpectrumCollection(...)`.
A subset of SP175 is still a `ReferenceSet` and still carries SP175's licence.

`metadata.KNOWN_TRUTH` is the new group, JSON-native as section 2.5 requires,
and `metadata.SET_LEVEL` documents what belongs on the set instead. The
`'C'` versus `'K'` case now **warns** rather than returning `None`: silence was
the defect, and the fix is not only to provide a better home but to stop the
old home lying about it.

Three things fell out that the ADR did not anticipate:

- **`from_cd`'s two branches collapsed into one.** A structural basis reads
  back as a one-hot composition, so "the coefficients are the composition" is
  the same arithmetic as mixing whole compositions. Section 3 predicted this;
  what it did not say is that the method argument then has nothing to switch
  on — so it is now *checked against* the kind of set supplied, which catches
  running SMP180 as though it were a basis of pure structures.
- **A composition can be asked for in a vocabulary the set never declared.**
  A bare category name does not say which DSSP states it covers, so
  `compositions` refuses rather than inventing a `Category` with no states.
- **A structural basis declares only the categories it contains.** Declaring
  the full four-category vocabulary would have made every fit report a
  confident `0.0` for a category the basis has no spectrum for — a claim the
  method never made, which ADR-0002 section 7.2 distinguishes from `None`.

**Deferred, deliberately: `Library` is not yet a subclass.** Section 3 says it
should be. It iterates over `Reference` objects where a collection iterates
over `Spectrum`, and that loop is inside `unmix()`, whose signature freezes at
1.0. Changing it as a side effect of CD work would put a frozen UV-Vis
signature at risk for no CD benefit. What was done instead is the useful half:
`unmix()` **accepts** a `ReferenceSet`, so the three cases already meet at the
call site. The inheritance is worth finishing on its own, with the UV-Vis
tests watching.

## 7. Open — the abstractions may not all be right yet (James, 2026-08-13)

**To be resolved before the format freezes.** Raised while reviewing the above:
imagine a panel of double mutants, two sites varying — AA, AC, AS, CA, CC, CS,
SA, SC, SS. A perfectly good set, and one this model has no name for.

What it is not: it is not a `ReferenceSet` (no known truth per spectrum), and
under the admission rule in section 6 it does not earn a `kind` either, because
nothing would interpret data the base class does not. So the strain is not on
the new class. It is on **`parameter`**, which the schema defines as *one
continuous scalar*:

| | today |
|---|---|
| `collection.parameters` | all `nan` — there is no number to put there |
| `sorted_by_parameter()` | raises, and there is nothing to sort along anyway |
| `to_matrix(with_parameter=True)` | raises |
| `group_by('sample')` | works only by flattening `('A','C')` into `"AC"` |
| `group_by(lambda s: s.metadata['site1'])` | **works today**, with keys the schema does not know |

So the model is not wrong, it is **incomplete**: there is a named,
schema-blessed concept for one continuous scalar and nothing for *k*
categorical factors. A factorial design is usable now, with caller-chosen keys.

**Making `parameter` a list is the wrong fix**, on this project's own grounds.
`parameters` returns a float array, and `sorted_by_parameter`, the van 't Hoff
fit and `library.from_series` all do arithmetic on it. A `parameter` that is
sometimes a 2-tuple makes those fail in the interesting way rather than the
loud way — the same failure class this ADR exists to remove. A factorial design
wants a *different* concept beside `parameter`, not inside it: something like
`metadata['factors'] = {'site1': 'A', 'site2': 'C'}` with the factor names at
set level.

**Why it has a deadline.** Metadata keys freeze with the format at 1.0 — the
same argument as roadmap D2. If `factors` is the right name and shape, agreeing
it costs an afternoon now; agreeing it after 1.0 means either a second set of
names or a format migration. Deciding *not* to have one is also a decision, and
also has to be made before November rather than arrived at.

Deliberately not built: there is no dataset in hand to check a design against,
and this project's working agreement is that a reader comes after a real file.

## 8. What would justify revisiting this

- A reference set whose per-item known truth is not per-item — a set where the
  structures are only known jointly, for example as a covariance rather than a
  value each. Nothing like that is in use here.
- Per-item data large enough that gathering it on every access matters. The
  gathered views are O(n) over a list of at most a few hundred; if reference
  sets reached thousands, caching with invalidation would be worth the
  complexity it costs.
