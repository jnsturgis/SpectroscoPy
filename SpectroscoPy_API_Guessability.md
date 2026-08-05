# API guessability audit (2026-08-05)

**One question, asked of every public name: what would somebody who knows
numpy, scipy and pandas — or a language model that half-remembers this
library — *guess* it was called?**

The reason to care is adoption. A user who guesses right writes working code
on the first attempt; a user who guesses wrong pays for it every time. The same
is true of an assistant writing code against this library, which is
increasingly how a small scientific package first gets used. But the test is
not an LLM trick: if a change only helps the model and not the person, it is
gaming, and it does not belong here. Every recommendation below improves the
human ergonomics identically.

The reason to do it **now** is §14.3: September is the last window in which a
public name can change without costing a major version. After the freeze these
become permanent.

## How findings are graded

| | Meaning | Why it matters |
|---|---|---|
| **A** | The wrong guess **works**, and gives a wrong result | The expensive kind. Nothing fails, so nothing is investigated |
| **B** | The wrong guess raises | Costs one attempt and a look at the docs. Cheap to fix with an alias |
| **C** | Cosmetic or already documented | Note it, decide once, move on |

---

## A1. `spectrum.technique = 'ATR-FTIR'` silently mislabels the axes

The worst finding, and a real hazard rather than a naming preference.
`technique` is a plain writable attribute, so the guessable line runs — and
does about half of what `set_type` does:

```
direct .technique =  technique='ATR-FTIR'  x_unit='nm'     x_quantity='Wavelength'   spec_type=None
set_type()           technique='ATR-FTIR'  x_unit='cm^-1'  x_quantity='Wavenumber'   spec_type='ATR-FTIR'
```

An ATR-FTIR spectrum labelled in **nm**. Every plot from it is wrong, `to()`
converts against the wrong axis, and `metadata['spec_type']` is never set, so
anything reading that key silently sees an untyped spectrum.

Nothing raises. Nothing warns. `spectrum.technique` afterwards reports exactly
what was asked for, so the one check a user would think to make passes.

**Recommended:** make `technique` a property whose setter calls `set_type`, so
the guessable line becomes the correct line. Keep `set_type(..., force_units=)`
for the case that needs the extra argument. This is strictly additive for
existing code and removes the trap entirely.

*Alternative if that is too much: make the attribute read-only, so the guess
raises with a message naming `set_type`. Worse for ergonomics, but honest.*

## A2. `set_parameters()` looks like it mutates, and does not

`Spectrum.set_sample()`, `set_reference()`, `set_type()` and `set_parameter()`
all mutate in place and return `None`. `SpectrumCollection.set_parameters()`,
added 2026-08-04, returns a **new** collection and leaves the original
untouched. Same prefix, opposite contract.

```python
collection.set_parameters([...])      # looks done; the collection is unchanged
```

Returning a new collection is the right behaviour — it matches the
non-mutating style of every other collection method. The name is what is
wrong, and it is my own from yesterday.

**Recommended:** rename to `with_parameters()`, which reads as returning a
copy, and reserve the `set_` prefix for the mutating `Spectrum` family. Cheap
now, permanent in November.

## A3. `baseline()` returns the baseline, not the corrected spectrum

`baseline()` gives the estimated background; `baseline_correct()` gives the
spectrum with it removed. Both return a `Spectrum`, so the wrong one plots
happily and looks plausible — a smooth curve through the data is not obviously
an error at a glance. **This one has already been made once in this project**,
in the user guide, where `spectrum.baseline(...)` was plotted as though it were
the corrected trace.

Graded A because the failure is quiet, though it is at the milder end: someone
looking closely at the figure would notice.

**Recommended:** keep both, and make the pairing self-describing —
`estimate_baseline()` and `subtract_baseline()`, with the current names kept as
aliases. If only one change is made, renaming `baseline()` is the one that
pays, because it is the name that reads like a verb but is not.

---

## B. Wrong guesses that fail loudly

Each costs a user one failed attempt. All are fixable with an alias, which is
additive and freeze-safe.

| Current | What will be guessed | Note |
|---|---|---|
| `collection.group_by('sample')` | `groupby` | pandas muscle memory, typed without thinking. The highest-frequency item in this table |
| `spectrum.normalize(...)` | `normalise` | The documentation prose says "normalise" throughout while the method is American. The docs themselves teach the wrong spelling |
| `collection.select(predicate)` | `filter` | `select` means something different in pandas; `filter` is the verb for "keep the ones matching" |
| `collection.to_matrix()` | `to_numpy`, `.values` | `to_numpy()` is the modern convention |
| `spectrum.save_as(path)` | `write`, `save`, `to_file` | Pairs with `Spectrum.read()`, so read/write or load/save would be symmetrical; `read`/`save_as` is neither |
| `spectrum.get_info()` | `describe()`, `info()` | Java-shaped, and `describe_history()` already exists, which makes `describe()` the natural guess |
| `library.coefficient(...)` | confused with `library.Coefficient` | A function and a class differing only in case. Both are public and both are plausible completions |
| `peaks.strongest(5)` | `nlargest`, `top`, `largest` | Minor; `strongest` is good spectroscopy English |

**Recommended:** add aliases for `groupby`, `normalise` and `filter`; decide
`to_numpy` and the save/write symmetry deliberately; rename `get_info` to
`describe` with an alias. Aliases cost nothing and never need removing.

---

## C. Noted, low priority

- **`processing.common`** is not a name anyone guesses. The array-level
  functions it holds are genuinely useful and effectively undiscoverable
  except through the guide.
- **`set_sample(info)`** — the parameter is called `info`; it is a name.
- **`clip()`** — already a documented deprecated alias for `crop`. Fine.
- **The top-level shims `calc`, `formats`, `tools_spc`** still emit
  "will be removed in 0.2", and there is no 0.2 — the next release is 1.0.
  This is §14.2 blocker 5, still open. They also occupy three guessable
  top-level names. Either they go before the freeze or that promise is
  rewritten; they cannot survive silently into a version that promises
  stability.

---

## What is already right, and should be protected in the freeze

Worth writing down so it does not get "tidied" later:

- **`find_peaks(...)`** matches `scipy.signal.find_peaks`, and passes unknown
  keywords through to it. Someone who knows scipy is already fluent here.
- **`mean()`, `std()`, `sem()`, `median()`** on a collection — exactly the
  numpy names, doing exactly the numpy thing.
- **`crop`, `resample`, `derivative`, `smooth`, `x_label`, `y_label`,
  `to_dataframe`** all land where a reader expects.
- **`read_spectrum` / `read_spectra`** — the plural distinction is guessable in
  both directions, and the singular form raising on a multi-spectrum file is
  the right call.
- **`sorted_by_parameter()` / `sorted_by_position()`** — consistent, explicit,
  and say that they return something rather than reorder in place.
- **The error messages.** Repeatedly, an error here says what to do instead.
  That is worth as much as any name: a wrong first attempt that is corrected by
  reading the exception costs nothing much, whereas `ValueError: invalid input`
  costs a documentation search. This is a genuine strength and should be an
  explicit review criterion for new code.

## The one that dwarfs the rest

**`pip install spectroscopy` and `import spectroscopy` must be the same word,
and must be this library.** Every guess above is a rounding error next to the
install line being wrong. That is the PyPI transfer, in progress with Joey
Knightbrook and blocked on his account recovery
([pypi/support#11748](https://github.com/pypi/support/issues/11748)); chase
15 September. Until it lands, the most natural thing anybody types — human or
assistant — gets the wrong package.

## Method, for repeating this

For each public name, in order: (1) write down the guess before looking at the
code; (2) run the guess; (3) grade by what happened — silently wrong, raised,
or worked. Step 2 is the one that cannot be skipped. A1 was found by running
the guess and reading the resulting units, not by thinking about the name.
