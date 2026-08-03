# Galactic / Thermo `.spc` — verified format notes

Working notes for the `.spc` reader that roadmap §15.3 has been waiting to
write. Recorded 2026-08-03, **before** any code, so that the code can be
checked against something written down independently of it.

Nothing here is copied from an implementation. See §4 for why that matters.

## 1. Provenance

Three sources, in descending order of authority. They agree.

**A. "Galactic Universal Data Format Specification", 9/4/97 — ✅ the
authoritative one, and it contains `SPC.H` itself.** Recovered from the
Internet Archive's copy of Galactic's own site:

```
https://web.archive.org/web/20040116104459id_/http://www.thermogalactic.com:80/instruments/spcdata/gspc_udf.pdf
```

50-odd pages, with the complete `SPC.H` header reproduced in Appendix A —
every `#define`, every structure, every enumerated code. This is the document
the brief guide keeps referring to, and it settles everything §3.7 was
previously guessing at. Local copy: `temp/spc-reference/galactic-udf-spec.pdf`
(with extracted text alongside it, which is the more useful form).

**B. Thermo Galactic, "A Brief Guide to SPC File Format and Using GSPCIO"
(2001).** The short version — block structure plus field tables. Copies at
<https://docs.c6h6.org/docs/assets/files/spc-3bb9ec9e4c158c5418bcfcc970be77f1.pdf>
and, via the docuri viewer, at `docuri.com/downloadFile/59c1d322f581710b28653306`
(the "download doesn't work" page: the real file is behind the pdf.js embed).
Both are byte-for-byte the same document, and it is the one whose printed
offsets are wrong — §2.

The Met Office `ensembles-eu` copy that search still turns up now 301s to an
unrelated page. Dead end.

**C. The `struct` format strings in a working reader** (Rohan Isaac's `spc`,
in `temp/spc2.py`). Used only as a *cross-check on the layout*, and now
superseded by A — which is just as well, because A shows it has a real bug
(§3.9). See §4 for why it could not have been used anyway.

**D. Twelve real `.spc` files written by Galactic's own software** — §5.
Every claim in §3 has been executed against them.

## 2. The finding that matters: the vendor's printed offsets are wrong

The PDF's "Byte #" column does not add up, and does not agree with itself.
It prints `8 double first X` then `15 double last X` — eight bytes apart would
be 16. From that point on the printed offsets run 2 low, with a further glitch
at the source-instrument field, and the final row (`324 char * 187`) ends at
511, not the 512 the same page declares.

The **field order, data types and sizes** are unambiguous and agree exactly
between sources. Laid out with no padding they sum to exactly 512, 32 and 64
bytes respectively. So:

> Take the field sequence from the vendor document; compute the offsets
> yourself. Do not transcribe the printed byte numbers.

Somebody working straight from the PDF would produce a reader that is wrong
from `flast` onward, and it would fail in a nasty way — plausible-looking
numbers, wrong axis.

## 3. Verified layout

### 3.1 Main header — 512 bytes, little-endian, no padding

`fversn` `0x4B` selects this. Offsets are computed, not transcribed.

| Off | Size | struct | Field | Meaning |
|----:|-----:|:------:|:------|:--------|
| 0 | 1 | `c` | `ftflgs` | flag bits, §3.4 |
| 1 | 1 | `c` | `fversn` | `0x4B` new LSB, `0x4C` new MSB, `0x4D` old |
| 2 | 1 | `c` | `fexper` | experiment type |
| 3 | 1 | `c` | `fexp` | Y exponent; `0x80` = already float, §3.5 |
| 4 | 4 | `i` | `fnpts` | points — **or a byte offset**, §3.6 |
| 8 | 8 | `d` | `ffirst` | first X |
| 16 | 8 | `d` | `flast` | last X |
| 24 | 4 | `i` | `fnsub` | number of subfiles |
| 28 | 1 | `c` | `fxtype` | X unit code, §3.7 |
| 29 | 1 | `c` | `fytype` | Y unit code, §3.7 |
| 30 | 1 | `c` | `fztype` | Z unit code |
| 31 | 1 | `c` | `fpost` | posting disposition |
| 32 | 4 | `i` | `fdate` | packed date, §3.8 |
| 36 | 9 | `9s` | `fres` | resolution text |
| 45 | 9 | `9s` | `fsource` | source instrument text |
| 54 | 2 | `h` | `fpeakpt` | interferogram peak point |
| 56 | 32 | `32s` | `fspare` | 8 spare floats |
| 88 | 130 | `130s` | `fcmnt` | memo |
| 218 | 30 | `30s` | `fcatxt` | custom X/Y/Z axis strings, NUL-separated |
| 248 | 4 | `i` | `flogoff` | byte offset to log block; 0 = none |
| 252 | 4 | `i` | `fmods` | file modification flags |
| 256 | 1 | `c` | `fprocs` | processing code |
| 257 | 1 | `c` | `flevel` | calibration level + 1 |
| 258 | 2 | `h` | `fsampin` | sub-method sample injection number |
| 260 | 4 | `f` | `ffactor` | concentration factor |
| 264 | 48 | `48s` | `fmethod` | method file |
| 312 | 4 | `f` | `fzinc` | Z increment for even-Z multifiles |
| 316 | 4 | `i` | `fwplanes` | number of W planes |
| 320 | 4 | `f` | `fwinc` | W increment |
| 324 | 1 | `c` | `fwtype` | W axis unit code |
| 325 | 187 | `187s` | `freserv` | reserved |

`<cccciddicccci9s9sh32s130s30siicchf48sfifc187s` — `struct.calcsize` = 512.

### 3.2 Subheader — 32 bytes, one per subfile, present even in single files

| Off | Size | struct | Field | Meaning |
|----:|-----:|:------:|:------|:--------|
| 0 | 1 | `c` | `subflgs` | 1 changed, 8 no peak table, 128 arithmetic |
| 1 | 1 | `c` | `subexp` | Y exponent for *this* subfile; `0x80` = float |
| 2 | 2 | `h` | `subindx` | subfile index |
| 4 | 4 | `f` | `subtime` | starting Z |
| 8 | 4 | `f` | `subnext` | ending Z |
| 12 | 4 | `f` | `subnois` | noise value |
| 16 | 4 | `i` | `subnpts` | points, **XYXY only**; 0 means use `fnpts` |
| 20 | 4 | `i` | `subscan` | co-added scans |
| 24 | 4 | `f` | `subwlevel` | W value |
| 28 | 4 | `4s` | `subresv` | reserved |

`<cchfffiif4s` — 32 bytes. Note the order: **three** floats, then two longs,
then one float. Guessing four-floats-then-two-longs also sums to 32 and is
wrong; it silently swaps noise/points/scans/W.

### 3.3 Log header — 64 bytes, at `flogoff`

`<iiiii44s`: `logsizd`, `logsizm`, `logtxto`, `logbins`, `logdsks`, reserved.
The ASCII text starts at `flogoff + logtxto` and runs `logsizd` bytes. In
practice it is `key=value` lines and is where instrument acquisition settings
live — the `.spc` counterpart of what the OPUS reader keeps as
`metadata['opus_history']`.

### 3.4 `ftflgs` bits

| Bit | Name | Meaning |
|----:|:-----|:--------|
| 1 | `TSPREC` | Y stored 16-bit, not 32-bit |
| 2 | `TCGRAM` | use experiment extension |
| 4 | `TMULTI` | multifile |
| 8 | `TRANDM` | multifile, Z randomly ordered |
| 16 | `TORDRD` | multifile, Z ordered but uneven |
| 32 | `TALABS` | use custom axis labels in `fcatxt` (obsolete) |
| 64 | `TXYXYS` | each subfile has its own X array |
| 128 | `TXVALS` | XY file — a global X array precedes the Y data |

The three storage layouts follow from the last two bits, and this is the whole
of the decision the reader has to make:

- neither → **`gx-y`**: X generated as `linspace(ffirst, flast, fnpts)`.
- `TXVALS` only → **`x-y`**: one float32 X array of `fnpts` after the header,
  shared by every subfile.
- `TXYXYS` → **`-xy`**: every subfile carries its own X and Y, and there is a
  directory block.

### 3.5 The Y exponent

Y is stored as fixed point unless `fexp` (or `subexp`) is `0x80`:

```
float_y = (2 ** exp) * int_y / (2 ** 32)     # 32-bit
float_y = (2 ** exp) * int_y / (2 ** 16)     # 16-bit, TSPREC set
```

Getting this wrong scales the whole spectrum by a power of two — which looks
like a units problem, not a bug, so it is worth a test.

### 3.6 `fnpts` is overloaded

When `TXYXYS` is set, `fnpts` is **not** a point count: it holds the byte
offset of the subfile directory. The directory is `fnsub` entries of
`<iif` — offset, size, Z time. Point counts then come from each `subnpts`.

This is the sharpest edge in the format and the reason a reader should branch
on the flag early rather than treat `fnpts` as a count everywhere.

⚠️ **And `fnpts = 0` means there is no directory at all** — verified against a
real file in §5.1. The subfiles then follow the header directly and each
`subnpts` carries its own count. Seeking unconditionally to a directory when
`TXYXYS` is set will read the main header as directory entries.

### 3.7 Axis unit codes — ✅ authoritative, from `SPC.H`

`fxtype`, `fztype` and `fwtype` share one table:

| | | | | | |
|--:|:--|--:|:--|--:|:--|
| 0 | Arbitrary | 11 | Days | 22 | Data points |
| 1 | **Wavenumber (cm⁻¹)** | 12 | Years | 23 | Milliseconds |
| 2 | Micrometers (µm) | 13 | **Raman shift (cm⁻¹)** | 24 | Microseconds |
| 3 | **Nanometers (nm)** | 14 | eV | 25 | Nanoseconds |
| 4 | Seconds | 15 | *(text labels, old `0x4D` only)* | 26 | Gigahertz |
| 5 | Minutes | 16 | Diode number | 27 | Centimeters |
| 6 | Hertz | 17 | Channel | 28 | Meters |
| 7 | Kilohertz | 18 | Degrees | 29 | Millimeters |
| 8 | Megahertz | 19 | Temperature (°F) | 30 | Hours |
| 9 | Mass (m/z) | 20 | Temperature (°C) | 255 | Double interferogram |
| 10 | Parts per million | 21 | Temperature (K) | | |

`fytype`:

| | | | | | |
|--:|:--|--:|:--|--:|:--|
| 0 | Arbitrary intensity | 9 | Millivolts | 21 | Temperature (K) |
| 1 | Interferogram | 10 | log(1/R) | 22 | Refractive index |
| 2 | **Absorbance** | 11 | Percent | 23 | Extinction coeff. |
| 3 | Kubelka-Munk | 12 | Intensity | 24 | Real |
| 4 | **Counts** | 13 | Relative intensity | 25 | Imaginary |
| 5 | Volts | 14 | Energy | 26 | Complex |
| 6 | Degrees | 16 | Decibel | 128 | **Transmission** |
| 7 | Milliamps | 19 | Temperature (°F) | 129 | Reflectance |
| 8 | Millimeters | 20 | Temperature (°C) | 130 | Arbitrary / single beam |
| | | | | 131 | Emission |

Note the gaps — 15, 17, 18 in the Y table are genuinely undefined, not an
omission here. `SPC.H` flags 128 as the threshold above which *"ALL HIGHER
MUST HAVE VALLEYS"*: codes ≥ 128 are quantities whose bands point down. That
is a real display/processing distinction and worth preserving in metadata even
though our own `units` vocabulary does not currently encode it.

These should map onto the project's `units` vocabulary rather than being
stored as raw integers.

### 3.8 Experiment type `fexper` — and why not to rely on it

| | | | |
|--:|:--|--:|:--|
| 0 | General SPC | 8 | X-ray diffraction |
| 1 | Gas chromatogram | 9 | Mass spectrum |
| 2 | General chromatogram | 10 | NMR spectrum or FID |
| 3 | HPLC chromatogram | 11 | Raman spectrum |
| 4 | FT-IR / FT-NIR / FT-Raman | 12 | Fluorescence spectrum |
| 5 | NIR spectrum | 13 | Atomic spectrum |
| 7 | UV-VIS spectrum | 14 | Chromatography diode array |

**There is no code 6.** The enumeration skips it. See §3.9.

Empirically (§5): **all twelve real files report `fexper = 0`**, including the
Raman, NMR, mass spectrum and X-ray ones. `SPC.H` explains why — in older
software `TCGRAM` had to be set for `fexper` to mean anything. So:

> Do not use `fexper` to decide what technique a file holds. Use `fxtype` and
> `fytype`, which are populated correctly in every real file tested.

### 3.9 A bug in the reference implementation

`SPC.H` skips code 6, but `temp/spc2.py` stores the experiment names in a
plain Python list and indexes it directly. From 6 onward every label is
therefore shifted by one: a file marked `9` (mass spectrum) is reported as
*NMR*, `10` (NMR) as *Raman*, and so on.

We would have inherited this silently by copying. It is a small vindication of
writing §3 from the specification first, and a reason to be wary of the other
enumerations in that file.

### 3.10 Packed date

`fdate` is a bitfield in one 32-bit word: year 12, month 4, day 5, hour 5,
minute 6, packed from the top. Zero means no date. Verified: the sample files
decode to 1987-05-28, 1989-07-28, 1992-03-04, 1997-10-01 and so on — all
plausible, none absurd, which a wrong bit layout would not manage.

### 3.11 `fmods` — what was done to the spectrum

A bitfield, one letter per operation, which is the `.spc` counterpart of the
OPUS history block: `A` averaging (2¹), `B` baseline (2²), `C` interferogram→
spectrum (2³), `D` derivative/integrate (2⁴), `E` resolution enhancement (2⁶),
`I` interpolation (2⁹), `N` smoothing (2¹⁴), `O` other arithmetic (2¹⁵),
`S` spectral subtraction (2¹⁹), `T` truncation (2²⁰), `W` collection time
modified (2²³), `X` X-units conversion (2²⁴), `Y` Y-units conversion (2²⁵),
`Z` zap (2²⁶).

Worth decoding into `metadata` rather than dropping. §21 is the precedent:
the OPUS work turned up files that had already been reference-subtracted
inside the instrument software, which nothing in the text export revealed.
`fmods` is where `.spc` records the same thing.

## 4. Licensing — the constraint on how this gets built

**Rohan Isaac's `spc` (`temp/spc2.py`) is GPL-3.0.** SpectroscoPy is MPL-2.0.
GPL-3.0 code cannot be taken into an MPL-2.0 project — the compatibility runs
the other way (MPL code may be combined into a GPL work, not the reverse). So
that file must not be copied from, adapted, or committed.

`temp/` is therefore in `.gitignore`.

What *was* taken from it is the field sequence and sizes — facts about a
third-party binary format, originating in Thermo's `SPC.H`, and used here only
to confirm that the vendor document had been read correctly. Facts about an
interface are not the copyrightable part. The reader will be written from §3
of this document, which is why §3 exists.

`temp/spc.py` is a different thing again: a BSD-3 `specio` plugin by Guillaume
Lemaitre. Permissively licensed, but it carries no format knowledge at all —
it just calls `spc.File(...)`. It is useful only as a worked example of what
an integration layer chooses to expose (it maps `gx-y`/`x-y` to one spectrum
with a shared axis, and `-xy` to a list — the same split our registry's
`multi=True` path already makes).

`spc-spectra`, the maintained fork, publishes no license metadata on PyPI
while being a derivative of GPL-3.0 code. Not usable either.

## 5. ✅ Real files — obtained, and the layout is confirmed against them

The same Internet Archive crawl that yielded the specification also preserved
**Galactic's own demonstration data set**, twelve files written by GRAMS
itself, in `temp/spc-reference/samples/`:

```
http://www.thermogalactic.com/products/dataviewer/spcfiles/*.spc
```

These are better than anything generated from the specification, because they
were written by the software that defined the format.

| File | `ftflgs` | npts | Range | X | Y |
|:--|--:|--:|:--|:--|:--|
| `ir-nh4.spc` | 0 | 8289 | 4398.10 → 401.23 | cm⁻¹ | Transmission |
| `raman.spc` | 0 | 1260 | 96.31 → 1728.00 | Raman cm⁻¹ | Counts |
| `uv-holm.spc` | 0 | 901 | 250 → 700 | nm | Absorbance |
| `nir-poly.spc` | 0 | 700 | 1100 → 2498 | nm | Absorbance |
| `vis-mirr.spc` | 0 | 692 | 890 → 199 | nm | Reflectance |
| `nmr-unk.spc` | 0 | 8192 | 14.83 → −5.24 | ppm | Arbitrary |
| `fluor.spc` | 0 | 491 | 255 → 500 | nm | Arbitrary |
| `gamma.spc` | 0 | 4076 | 0 → 4075 | Arbitrary | Arbitrary |
| `xraydiff.spc` | `0x20` | 18014 | 9.88 → 99.93 | Arbitrary | Counts |
| `ms-barb.spc` | `0xF0` | *(XYXY)* | 27.05 → 160.90 | Arbitrary | Arbitrary |
| `gc_gasln.cgm` | `0x02` | 54000 | 0 → 30 | min | Millivolts |
| `hplc.cgm` | `0x02` | 27200 | 0 → 34 | min | Arbitrary |

**The decisive check.** For every evenly-spaced file, the size predicted from
the parsed header —

```
512 + (4 × fnpts if TXVALS) + fnsub × (32 + 4 × fnpts)
```

— equals the actual file length **exactly**, to the byte, for all eleven. A
single wrong offset anywhere in §3.1 would put garbage in `fnpts` or `fnsub`
and this would not close. The axis values are independently sensible: an IR
spectrum descending 4398→401 cm⁻¹, NMR descending 14.83→−5.24 ppm, a holmium
oxide UV standard over 250–700 nm.

### 5.1 What `ms-barb.spc` settled

The one XYXY file resolved two things that the specification alone left
ambiguous, and both are traps:

**With `fnsub = 1` there is no directory block.** `fnpts` is `0` — meaning *no
directory offset*, not *no points* — and the subheader sits directly at 512
with the real count in `subnpts` (37). Size works out as
`512 + 32 + 8 × 37 = 840`, exactly. A reader that unconditionally seeks to a
directory when `TXYXYS` is set will read the file header as directory entries.

**The fixed-point exponent is not hypothetical.** `subexp = 15`, not `0x80`,
so the Y values are integers needing `2^15 × y / 2^32`. Read naively as
float32 they come out around 8×10¹⁵. Applied correctly they become 562 …
11798 — and, tellingly, *every one is an exact integer*, which is the
signature of correctly recovered fixed-point counts and makes a good assertion
in the test. The resulting spectrum is chemically sound for a barbiturate:
base peak m/z 119, with 91 (tropylium), 105, 57, and m/z 28/32 air background.

### 5.2 Still not covered

The set exercises `gx-y`, `x-y`, XYXY, custom axis labels and fixed-point Y.
It does **not** cover:

- `TSPREC` — 16-bit Y storage (`/2^16` rather than `/2^32`).
- A plain `TMULTI` multifile with several evenly-spaced subfiles — the
  ordinary many-spectra case the registry's `multi=True` path exists for.
- The old `0x4D` 256-byte-header format.
- `0x4C`, new-format MSB-first. Rare; detect and refuse clearly.
- A log block with real `key=value` acquisition text.

So James's archive is still worth finding — but it is no longer blocking.
Anything from it that is a **multifile** or an **old-format** file is what
would add coverage; another single-spectrum FTIR file would not.

## 6. Consequences already agreed

- The reader **sniffs, it does not trust the extension**: `.spc` is shared
  with Bruker EPR data (26 such files on this machine, in EasySpin's tests).
  Check `fversn` ∈ {`0x4B`, `0x4C`, `0x4D`} at offset 1, and say *"this looks
  like Bruker EPR data, which is a different format"* rather than "corrupt".
- Under working agreement §14.7, `.spc` is documented, so a **writer** is
  permissible here — the only format in the project where that is true.
- `0x4C` (new MSB-first) is big-endian and rare; detect it and raise a clear
  "not supported, please send the file" rather than mis-parsing.
- Shimadzu's `0xCF` variant is not SPC and is out of scope.

## 7. Open question for James: the sample files as test fixtures

The twelve files in `temp/spc-reference/samples/` are exactly what the test
suite wants, and they are ideal for it — vendor-written, varied, and covering
the two paths most likely to break.

They are also **Thermo Galactic's demonstration data from around 2003**, and
their licensing is unstated. They were freely published for developers working
with the format, which is the use here, and the company and its website are
long gone; but "no stated licence" is not the same as "ours to redistribute",
and committing them puts them in an MPL-2.0 repository and onto PyPI in every
sdist.

Three options, in the order I would rank them:

1. **Keep them out of the repository; check them into the test suite as an
   optional fixture** that skips when absent, with a download script pointing
   at the Internet Archive URLs. Costs a skipped test in CI unless the archive
   is reachable, keeps the distribution clean. This is what §14.7's caution
   about other people's material would suggest.
2. **Commit them** under a clearly labelled `tests/data/spc/` with provenance
   and the archive URLs recorded. Simplest, best CI, small residual risk.
3. **Generate our own fixtures** from the now-verified writer once it exists,
   and use the Thermo files only locally during development. Clean, but the
   tests then only prove self-consistency — the exact weakness §15.2 warned
   about.

My recommendation is 1, moving to 2 if the skipping proves annoying. Either
way the reader gets built from the specification, and this only decides what
ships.
