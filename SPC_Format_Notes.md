# Galactic / Thermo `.spc` — verified format notes

Working notes for the `.spc` reader that roadmap §15.3 has been waiting to
write. Recorded 2026-08-03, **before** any code, so that the code can be
checked against something written down independently of it.

Nothing here is copied from an implementation. See §4 for why that matters.

## 1. Provenance

Two independent sources were used, and they agree.

**A. Thermo Galactic, "A Brief Guide to SPC File Format and Using GSPCIO"
(2001).** Publicly downloadable; a copy was retrieved from
<https://docs.c6h6.org/docs/assets/files/spc-3bb9ec9e4c158c5418bcfcc970be77f1.pdf>.
This is the vendor's own document. It gives the block structure and a
field-by-field table of the 512-byte main header, the 32-byte subheader and
the 64-byte log header. It refers to `SPC.H` and to the fuller "Galactic
Universal Data Format Specification" for the enumerated code values; both ship
in the SPC File Developer's Kit, whose original download URL
(`thermogalactic.com/instruments/spcdata/spc_sdk.zip`) is long dead.

The other copy that search turns up, on the Met Office's ensembles-eu host,
now 301s to an unrelated page. If a second copy is ever needed, that is a dead
end.

**B. The `struct` format strings in a working reader** (Rohan Isaac's `spc`,
seen in `temp/spc2.py`). Used only as a *cross-check on the layout* — see §4.

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

### 3.7 Axis unit codes

The `fxtype`/`fytype`/`fztype` code tables are enumerations in `SPC.H`; the
brief guide does not reproduce them. Relevant to us: X `1` = wavenumber cm⁻¹,
`13` = Raman shift cm⁻¹, `3` = nm; Y `2` = absorbance, `1` = interferogram,
`4` = counts, `128` = transmission. These need checking against a real file or
`SPC.H` before being trusted — they are the one part of these notes with only
one source behind them, and they should map onto our existing `units`
vocabulary rather than being stored as raw integers.

### 3.8 Packed date

`fdate` is a bitfield in one 32-bit word: year 12, month 4, day 5, hour 5,
minute 6, packed from the top. Zero means no date.

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

## 5. Still outstanding

**Real `.spc` files.** The specification closes the gap §15.3 opened, but not
the one §15.2 cares about: a reader written from a specification and tested
against files I generated from the same specification only proves I read it
consistently. The OPUS reader is trustworthy because 43 real files came with
paired `.dpt` exports of the same measurements.

What would make this immediately buildable, in order of usefulness:

1. Three to five real `.spc` files with a text export of the same data.
2. One multi-subfile example — ideally XYXY, since §3.6 is where a reader
   breaks and where the registry's many-spectra path gets exercised.
3. One old-format (`0x4D`, 256-byte header) file, if any exist in the archive.

**`SPC.H` itself**, for the §3.7 code tables and the `fexper` enumeration.
The SDK zip is gone from Thermo's site; a copy may survive in a Linux distro
package or a mirror of the developer's kit.

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
