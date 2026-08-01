# Reading and writing files

## Just open it

```python
spectrum = spc.Spectrum("sample.dpt")
```

The format is worked out from the extension, and the details from the file
itself: separators, header rows and text encoding are all detected rather than
assumed.

| format | extensions | notes |
|---|---|---|
| Bruker OPUS data point table | `.dpt` | separator detected per file |
| JCAMP-DX | `.jdx` `.dx` `.jcamp` | including compound files |
| Delimited text | `.csv` `.tsv` `.txt` | header row detected |
| Generic table | — | many spectra from one wide file |
| SpectroscoPy native | `.spy` | keeps units and processing history |

`spc.io.describe_formats()` prints the live list.

Force a format when the extension lies:

```python
spectrum = spc.Spectrum("path/", "oddly_named_file", "dpt")
```

## Many spectra from one file

Some files hold a whole series — a spectrofluorimeter export, a chromatography
run. Those go through the generic table reader.

```python
series = spc.io.read_spectra("J_Peri.csv", 'table', paired=True)   # (x,y) pairs
wide   = spc.io.read_spectra("samples.csv", 'table', x_col=0)      # shared x
```

`read_spectra` always returns a `SpectrumCollection`. `read_spectrum` is the
convenience form for a single-spectrum file, and it **raises** if the file turns
out to hold several — returning the first of eighty-six silently would be worse.

The reader decides two things by looking rather than assuming: which header row
names the series (the one whose values actually distinguish them), and which
row supplies axis labels (one whose x and y entries differ, so a run name
repeated across every column is not mistaken for a label).

## Encoding

Byte-order mark first, then UTF-8, then latin-1 as a backstop that maps every
byte and so cannot itself fail. Chromatography exports are commonly UTF-16, and
older instruments put latin-1 bytes in comment fields; neither should stop you
reading numbers that are perfectly readable.

## Writing

```python
spectrum.save_as("cleaned.spy", "spy")     # keeps units and history
spectrum.save_as("cleaned.csv", "csv")     # just the numbers
```

Only `.spy` is lossless. Everything else keeps the data and loses the
provenance, which is fine for handing numbers to someone else and not fine as
your own working format.

The file type is validated *before* the file is opened, so an unwritable type
raises instead of truncating the target to nothing.

## Adding a format

One decorator, and everything downstream follows — extension inference,
`Spectrum(...)`, the collection loader:

```python
from spectroscopy.io import register_reader

@register_reader('mine', extensions=['.mine'], description='my instrument')
def read(handle, spectrum, **kwargs):
    ...                       # fill spectrum.x, spectrum.y, spectrum.metadata
```

Readers take an open text handle and a spectrum to fill. For a format holding
several spectra, pass `multi=True` and return a list instead.

Please contribute one back if you write it — a format nobody else can read is
the problem this library exists to solve.
