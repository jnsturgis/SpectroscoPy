# ADR-0004 — Reference sets are collections, not a parallel hierarchy

**Status:** **Accepted and implemented, 2026-08-13.** Proposed 2026-08-12
(James, reviewing the CD branch: *"is this well structured or just ad hoc?"*).
Ad hoc — this records what the structure should be instead. Section 6 records
what was built and the one thing that was deliberately deferred; section 5
remains open and is on the 1.0 critical path.
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

## 5. Open — and it needs deciding before 1.0

**Collections cannot be saved.** There is no collection writer; spectra are
saved one at a time. So set-level data has nowhere to persist, and a
`ReferenceSet` loaded from `load_dichroweb_basis` cannot be written back out
with its licence and citation attached.

This is not a reference-set problem. It is a `.spy` question, and `.spy`
freezes at 1.0 (roadmap §14.2), so it is on the critical path whatever is
decided:

- a container format holding several spectra plus set-level data;
- or a sidecar manifest beside per-spectrum files, which is what
  `library.load_basis` already reads;
- or an explicit decision that collections are assembled at load time and
  never serialised, in which case set-level data is always reconstructed
  from its source and never round-trips.

The third is defensible and is the cheapest, but it should be *chosen* rather
than arrived at by not implementing the other two.

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

## 7. What would justify revisiting this

- A reference set whose per-item known truth is not per-item — a set where the
  structures are only known jointly, for example as a covariance rather than a
  value each. Nothing like that is in use here.
- Per-item data large enough that gathering it on every access matters. The
  gathered views are O(n) over a list of at most a few hundred; if reference
  sets reached thousands, caching with invalidation would be worth the
  complexity it costs.
