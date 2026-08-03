---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Reading and writing files

Every example on this page runs when the documentation is built. If a call
here stops working, the build fails rather than the page quietly becoming a
lie.

```{code-cell}
import numpy as np
import spectroscopy as spc
from spectroscopy import datasets
```

## Just open it

```{code-cell}
spectrum = spc.Spectrum.read(datasets.path("ethanol"))
spectrum
```

The format is worked out from the extension, and the details from the file
itself: separators, header rows and text encoding are all detected rather than
assumed.

```{code-cell}
print(spectrum.x_quantity, "in", spectrum.x_unit)
print(spectrum.y_quantity, "in", spectrum.y_unit)
print(len(spectrum.x), "points")
```

## What can be read

Generated from the registry, so it cannot drift from what the code does:

```{code-cell}
print(spc.io.describe_formats())
```

`r-` means read-only. The two binary instrument formats are read-only on
purpose: SpectroscoPy will not write a format it cannot verify it is writing
correctly. OPUS is undocumented, so nothing could establish that a file we
wrote is a file OPUS would accept.

## The binary formats sniff, they do not trust the extension

Two instrument formats have extensions that carry no reliable type
information, so both check the file's own magic value first.

* OPUS numbers its files by measurement — `sample.0`, `sample.1` — so the
  extension says nothing at all about what is inside.
* `.spc` is used by **Bruker for EPR data** as well as by Galactic. They share
  only the three letters.

Hand the reader an EPR file and it says so, by name, instead of calling it
corrupt:

```{code-cell}
import tempfile, pathlib

tmp = pathlib.Path(tempfile.mkdtemp())
fake_epr = tmp / "sample.spc"
fake_epr.write_bytes(bytes([0x00, 0xA0, 0x9B, 0xBE]) * 256)   # Bruker EPR

try:
    spc.Spectrum.read(fake_epr)
except ValueError as error:
    print(error)
```

## When one file holds several spectra

A native OPUS file often holds the same quantity twice: what was measured, and
what was left after something was done to it inside the instrument software.
Which one a `.dpt` export contains is not recoverable from the export, so the
reader returns them all and records what it knows.

```python
spectra = spc.io.read_spectra("sample.0")
for s in spectra:
    print(s.metadata["opus_block"], s.metadata.get("opus_history", ""))
```

```{note}
That block is not executed — it needs one of your own instrument files. The
multi-spectrum machinery below is exercised for real.
```

A multi-subfile `.spc` behaves the same way, and each subfile keeps its own Z
value — a time, a temperature, a depth — in `metadata['z_value']`.

## Many spectra from one file

Some files hold a whole series: a spectrofluorimeter export, a chromatography
run, a plate reader. Those go through the generic table reader.

```{code-cell}
wide = tmp / "series.csv"
x = np.linspace(400, 700, 6)
rows = ["wavelength,control,treated,washed"]
for i, value in enumerate(x):
    rows.append(f"{value},{i * 0.10:.3f},{i * 0.17:.3f},{i * 0.04:.3f}")
wide.write_text("\n".join(rows))

series = spc.io.read_spectra(wide, 'table', x_col=0)
print(len(series), "spectra:", [s.name for s in series])
```

The reader decided two things by looking rather than assuming: which header row
names the series — the one whose values actually distinguish them — and which
row supplies axis labels, one whose x and y entries differ, so a run name
repeated across every column is not mistaken for a label.

`read_spectra` always returns a `SpectrumCollection`. `read_spectrum` is the
convenience form for a single-spectrum file, and it **raises** if the file
turns out to hold several:

```{code-cell}
try:
    spc.io.read_spectrum(wide, 'table', x_col=0)
except ValueError as error:
    print(error)
```

Returning the first of three silently would be worse than failing.

## A `.csv` is often neither comma separated nor dot decimal

A spreadsheet exported under a French, German, Spanish or Italian locale
writes `0,1234`. The comma is then taken by the decimal, so Excel exports `;`
as the field separator instead. The file is still called `.csv`, and nothing
inside it says any of this.

Both are worked out from the numbers — and they have to be worked out
*together*, because neither can be settled alone:

```{code-cell}
locales = {
    'from Paris':   "Longueur d'onde;Absorbance\n400,5;0,1234\n401,0;0,2345\n",
    'from Berlin':  "Wellenlaenge\tExtinktion\n400,5\t0,1234\n401,0\t0,2345\n",
    'from Boston':  "Wavelength,Absorbance\n400.5,0.1234\n401.0,0.2345\n",
}

for origin, text in locales.items():
    path = tmp / f"{origin.replace(' ', '_')}.csv"
    path.write_text(text)
    spectrum = spc.Spectrum.read(path)
    print(f"{origin:<12} x={spectrum.x}  y={spectrum.y}")
```

Splitting `400,5;0,1234` on the comma finds `400` and `1234`, which look like
two perfectly good numbers until you notice the `5;0` stranded between them.
That is why the separator is scored on what it *leaves behind* as well as what
it produces.

Thousands separators are handled too — `1.400,5` and `1 400,5` are both
fourteen hundred and a half, including the non-breaking space Excel actually
emits:

```{code-cell}
from spectroscopy.io.table import parse_number, sniff_format

print(parse_number('1.400,5', ','), parse_number('1 400,5', ','))
print(sniff_format(["400,5;0,1234", "401,0;0,2345"]))
```

Pass `delimiter=` or `decimal=` to override the sniffing. A `.tsv` pins the
separator to a tab and still sniffs the decimal, because the locale problem is
independent of the separator.

Writing works the same way round, for when a file has to go back to somebody
whose spreadsheet expects it:

```{code-cell}
out = tmp / "pour_chloe.csv"
spectrum.save_as(str(out), "csv", decimal=',')
print(out.read_text()[:60])
```

Asking for a comma decimal moves the separator to `;` on its own, because a
comma cannot be both: `400,5,0,1234` is four fields or two and nothing can
tell which. Setting both to a comma explicitly raises rather than writing a
file no reader could get right — including this one.

```{code-cell}
try:
    spectrum.save_as(str(tmp / "bad.csv"), "csv", decimal=',', delimiter=',')
except ValueError as error:
    print(error)
```

## Encoding

Byte-order mark first, then UTF-8, then latin-1 as a backstop that maps every
byte and so cannot itself fail. Chromatography exports are commonly UTF-16 and
older instruments put latin-1 bytes in comment fields; neither should stop you
reading numbers that are perfectly readable.

```{code-cell}
utf16 = tmp / "akta-style.csv"
utf16.write_text("wavelength,absorbance\n400,0.10\n500,0.22\n600,0.31\n",
                 encoding="utf-16")

print("detected:", spc.io.detect_encoding(utf16))
print(len(spc.Spectrum.read(utf16).x), "points read")
```

## Writing

```{code-cell}
spectrum.save_as(str(tmp / "cleaned.spy"), "spy")   # keeps units and history
spectrum.save_as(str(tmp / "cleaned.csv"), "csv")   # just the numbers

restored = spc.Spectrum.read(tmp / "cleaned.spy")
print("units survived:", restored.x_unit, restored.y_unit)
```

Only `.spy` is lossless. Everything else keeps the data and loses the
provenance, which is fine for handing numbers to somebody else and not fine as
your own working format.

`save_as` passes anything else it is given to the writer, so a format's
options are reachable from here — `decimal=','` above being the one most
likely to be wanted.

The file type is validated *before* the file is opened, so asking for a type
that cannot be written raises instead of truncating the target to nothing:

```{code-cell}
try:
    spectrum.save_as(str(tmp / "nope.0"), "opus")
except ValueError as error:
    print(error)
```

## Adding a format

One decorator, and everything downstream follows — extension inference,
`Spectrum.read()`, the collection loader. Readers take an open text handle and
a spectrum to fill:

```{code-cell}
from spectroscopy.io import register_reader

@register_reader('demo', extensions=['.demo'],
                 description='two columns, semicolon separated')
def read_demo(handle, spectrum, **kwargs):
    rows = [line.split(';') for line in handle.read().splitlines()
            if line.strip()]
    spectrum.x = np.array([float(row[0]) for row in rows])
    spectrum.y = np.array([float(row[1]) for row in rows])
```

That is the whole registration. Nothing else needs editing — no dispatch
table, no list of extensions:

```{code-cell}
mine = tmp / "instrument.demo"
mine.write_text("400;0.10\n500;0.25\n600;0.40\n")

print("inferred as:", spc.io.infer_file_type(mine))
spc.Spectrum.read(mine)
```

For a format holding several spectra, pass `multi=True` and return a list
instead. For a binary format, pass `binary=True` so the registry does not try
to sniff a text encoding and corrupt it.

Please contribute one back if you write it — a format nobody else can read is
the problem this library exists to solve. See
[CONTRIBUTING.md](https://github.com/jnsturgis/SpectroscoPy/blob/main/CONTRIBUTING.md).
