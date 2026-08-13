# ADR-0005 — A spectrum carries quantities; the set says what they mean

**Status:** Proposed, 2026-08-13 (James). **Supersedes parts of ADR-0004** —
section 2.3's split of per-item metadata into "conditions" and "known truth",
and the reasoning that made `ReferenceSet` a class rather than a convenience.
**Depends on** ADR-0001 (the core data model) and the metadata schema of
roadmap D2.
**Affects** `metadata`, `SpectrumCollection`, `Spectrum.set_parameter`, the
`.spy` format, and eventually every analysis that reads a per-member value.

**Decided here:** where per-member quantities live, and how a set declares what
each one is. **Not decided here:** any generic machinery to manipulate them.
Section 6 says why that separation is the whole point.

---

## 1. What prompted it

Reviewing ADR-0004, James put the case that the taxonomy was wrong:

> the redox potential is a known truth provided by an external technique as
> much as the secondary structure breakdown … the secondary structure breakdown
> is just a point in a multidimensional space like redox potential is a point in
> a one dimensional space. So … spectra can have different metadata associated
> with them that can take different forms (dimensionality, units etc) and within
> a collection some of these make sense for different analyses.

That is right, and the codebase already contains the proof.

### 1.1 The same operation, at two dimensionalities

`library.from_series` and `processing.cd.estimate` are relatives that do not
know it:

| | `from_series` → `unmix` | `cd.estimate` |
|---|---|---|
| training pairs | (spectrum, concentration) | (spectrum, composition) |
| the value's space | ℝ, in M | the simplex in ℝ⁴, dimensionless |
| what it does | build ε, then invert it for an unknown | build a design, then invert it for an unknown |

**A reference set is a calibration set.** ADR-0004 treated CD reference proteins
as a special kind of object; they are the four-dimensional case of something the
UV-Vis half already did in one dimension.

### 1.2 Five spellings of one idea

`metadata.py` already flagged part of this, and understated it: `excitation_nm`
and `z_value` are described as "the same idea as `parameter`", with the note
that *three spellings of one concept is how a format ends up with three*. It is
not three. What distinguishes one member of a set from the next is spelled:

| set | what distinguishes members | spelled |
|---|---|---|
| thermal melt | temperature | `parameter` (+ name, unit) |
| redox titration | potential | `parameter` |
| EEM | excitation wavelength | `excitation_nm` |
| multi-subfile SPC | subfile z | `z_value` |
| reference proteins | protein identity | `sample` |
| double-mutant panel | (site 1, site 2) | **nothing** |

The last row is not an oversight in the schema. It is the schema being the
wrong shape: five keys for one role, and none of them able to hold a pair.

### 1.3 The finding that decides it: role is per-analysis, not per-key

ADR-0004 filed per-member data by *what the datum is about* — measurement
conditions, position in a series, known truth. That partition does not exist.

**`concentration` is a fixed measurement condition in a CD scan and the thing
you regress against in a dilution series.** Same key, same value, different
role. `temperature` is a condition in an FTIR spectrum and the axis of a melt.
The role is assigned by the analysis being run, not carried by the datum.

So a namespace organised by role must mis-file something, and did.

---

## 2. Decision

**A spectrum carries named quantities. The set declares what each one is.**

```python
# per member -- the values
spectrum.metadata['quantities'] = {
    'temperature': 45.0,
    'composition': {'helix': 0.69, 'sheet': 0.00, 'turn': 0.11,
                    'disorder': 0.20},
}

# per set -- what they mean
collection.info['quantities'] = {
    'temperature': {'space': 'scalar', 'unit': 'C'},
    'composition': {'space': 'simplex', 'unit': None,
                    'over': [Category('helix', {'G', 'H', 'I'}), ...],
                    'known_from': 'DSSP on the deposited structures'},
}
```

Three things follow, and they are the reason for the shape:

**The declaration is set-level because it cannot then disagree with itself.**
This is ADR-0004 section 2.4's argument, which survives intact: two members
labelled `'C'` and `'K'` gathered to a unit of `None`, indistinguishable from
never having been set.

**The value is member-level because there is then one list.** ADR-0004 section
2.2's argument, which also survives: a parallel table beside the spectra can be
reordered out of step, and the case that does not raise is the dangerous one.

**Analyses ask by name and by space.** `cd.estimate` needs a quantity whose
space is a simplex over DSSP categories; `melting.two_state` needs a scalar in
degrees. Each asks, and refuses clearly when the set has not got one. That is
James's third point from the ADR-0004 discussion — *parts of this program know
how to handle this data* — stated as a requirement on the analysis rather than
as a subclass of the collection.

### 2.1 What `space` may be

Deliberately a short closed list, not a type system (see section 6):

| space | value | example |
|---|---|---|
| `scalar` | a number | temperature, potential, concentration |
| `categorical` | a name | protein identity, technique, mutant |
| `vector` | `{component: number}` | a composition before normalisation |
| `simplex` | `{component: number}` summing to ~1 | a secondary structure breakdown |
| `tuple` | an ordered list of the above | a double mutant, `('A', 'C')` |

`unit` is a string or `None`, and belongs to the declaration rather than to the
value for the reason `metadata.py` already gives: a value and a unit stored as
two keys can disagree, and the pair that disagrees is indistinguishable from the
pair that does not. `over` names the components of a `vector` or `simplex` —
this is ADR-0004's `info['categories']`, generalised, and it is what stops a
bare `'helix'` being read without knowing whether it covers 3-10 and pi.

### 2.2 What moves, and what does not

**Moves into `quantities`:** anything whose name, space or unit is the caller's
to choose — `parameter`, `excitation_nm`, `z_value`, `composition`, `category`.
These are exactly the five spellings of section 1.2 plus the known-truth group.

**Stays where it is:** the sample conditions whose unit the schema fixes —
`path_length` in cm, `pH`, `mean_residue_weight` in Da. They need no
declaration, because there is nothing for a caller to declare; and
`metadata['path_length']` is read by `unmix` and `scattering` today.

*This line is the part of the decision I am least sure of, and it is drawn for
migration cost rather than from principle: `path_length` is as much a quantity
with a space and a unit as temperature is. The defensible version of the line is
that a quantity needs a declaration only when its meaning is not already fixed
by the schema. If that turns out to be a distinction without a difference, the
conditions can join later — additively, since they would move into a namespace
that already exists.*

**`sample` stays**, as a top-level key and as the default of `group_by`. In a
mutant panel the sample and the distinguishing quantity are the same value; in
a melt they are not, because every scan is of one sample. Keeping it costs one
key and answers a question the quantities cannot: *what is this a spectrum of*.

### 2.3 Consequence: `ReferenceSet` becomes a convenience, not a kind

Everything the class does — `.compositions`, `.categories`, `from_compositions`
— is "read a declared quantity of space `simplex`". Under this ADR that is a
general operation, so the class stops being the thing that makes reference sets
possible and becomes a shorthand for the commonest case.

ADR-0004 section 6 gave a rule for what earns its own class: *it reads data the
base class stores but does not interpret*. That rule survives, but nothing meets
it any more, because the base class now interprets declared quantities. The
`kind` recorded in a `.spy` collection block stays useful for round-tripping a
convenience class; it stops being load-bearing.

---

## 3. What this does not decide

- **No generic quantity machinery.** No `collection.quantity()`, no automatic
  unit conversion, no analysis rewritten to take "any quantity of space X".
  Analyses keep today's signatures and look up the name they need.
- **`ReferenceSet` is not deleted.** Section 2.3 says what it becomes, not when.
- **`unmix` and `cd.estimate` are not merged**, though section 1.1 says they
  could be.

---

## 4. When

**Recommendation: agree now, implement in September.** The keys freeze in
November and this changes keys, so it cannot wait past the September
breaking-change window (roadmap section 14.3). It should not land before then
either: the testers get 0.1.0 on 20 August, and moving the metadata schema
under them mid-trial costs the feedback the freeze depends on.

The risk of that choice, named: September is already carrying every other
breaking change, and an agreed ADR that nobody has time to execute is worse
than a smaller change made now.

---

## 5. Rejected

**Keeping `parameter` numeric and adding `factors` for the categorical case.**
Two mechanisms for one role, which is how there came to be five. It also
re-files by role, which section 1.3 shows cannot be done.

**Making `parameter` itself polymorphic.** `sorted_by_parameter`, the van 't
Hoff fit and `from_series` all do arithmetic on `collection.parameters`. A value
that is sometimes a pair makes those fail in the interesting way rather than the
loud way. The general concept needs a general home; it must not be smuggled into
a key whose consumers assume a float. **The type check belongs at the operation,
not at the key** — which is the opposite of what this project argued a day
earlier, when it defended the key instead of the fit.

**Flat keys with declarations only** — leaving `temperature` and `composition`
as top-level keys and adding just the set-level declaration. Smallest change to
the frozen surface, and rejected because it leaves no way to tell a declared
quantity from any other metadata key: `unknown_keys()` could not distinguish
"a quantity this set declares" from "a typo", which is most of what it is for.

**A full quantity type system with declared spaces and unit algebra.** ADR-0001
rejected `pint` for native units; section 15.4 of the roadmap deferred it again
after 1.0. Building one in eleven weeks, against datasets not yet in hand,
contradicts both — and the working agreement that a reader comes after a real
file applies to a schema as much as to a format.

---

## 6. Why the split between deciding and building

The namespace and the declaration format **freeze in November**. Everything else
— the lookup helpers, the analyses learning to ask by space, the merging of
`unmix` and `cd.estimate` — is **additive and can ship in 1.1** without breaking
a file or a signature.

So the expensive half is the cheap half to build, and the cheap half is the one
with the deadline. Deciding the namespace now costs an afternoon; deciding it
after 1.0 means a second set of names or a format migration.

---

## 7. What would justify revisiting this

- **A quantity whose space is not in section 2.1's list** and is not expressible
  as one — a per-member covariance, a spectrum-valued property, a distribution
  rather than a point.
- **The `path_length` line in section 2.2 proving to be arbitrary in practice**,
  which would argue for the sample conditions joining `quantities` as well.
- **Tester feedback showing the nesting hurts** — `metadata['quantities']['t']`
  is two lookups where there was one, and if people find that consistently worse
  than a flat key, that is evidence about the interface rather than about the
  model.
