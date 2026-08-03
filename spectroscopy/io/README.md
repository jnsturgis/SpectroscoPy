# Formats

Readers and writers for spectroscopy file formats. Nothing here dispatches on
file extension by hand: each module registers itself with the decorators in
`registry.py`, and `read_spectrum()` / `read_spectra()` / `write_spectrum()`
are derived from that one table.

A reader is either

* `read(handle, spectrum, **kwargs)` — fills a `Spectrum` handed to it, or
* `read(handle, **kwargs)` returning a list, registered with `multi=True`,
  for formats that hold several spectra in one file.

Binary formats register with `binary=True`, which stops the registry sniffing
a text encoding first — that attempt corrupts them.

**Layering:** this package must not import `spectroscopy.spectra` at module
scope. The dependency runs core → io, never back. Where a reader has to build
a `Spectrum` the import goes inside the function.

## Supported

| Module | Extensions | R/W | Notes |
|:--|:--|:--:|:--|
| `jcamp` | `.jdx` `.dx` `.jcamp` | rw | JCAMP-DX exchange format, including compound files |
| `csv` | `.csv` | rw | header row sniffed |
| `csv` | `.tsv` `.txt` | rw | tab separated |
| `table` | — | r- | generic delimited text, wide or paired columns |
| `dpt` | `.dpt` | rw | Bruker OPUS data point table, separator sniffed |
| `opus` | `.0` … `.9` | r- | Bruker OPUS native binary |
| `spc` | `.spc` `.cgm` | r- | Galactic/Thermo SPC binary |
| `spy` | `.spy` | rw | native; carries units and processing history |

## Two formats that need sniffing, not extensions

`opus` and `spc` both have extensions that carry no reliable type information,
so both check a magic value before trusting the name:

* OPUS numbers files by measurement (`sample.0`, `sample.1`), so the extension
  says nothing at all. It checks for `0xFEFE0A0A`.
* `.spc` is used by **Bruker for EPR data** as well as by Galactic — a
  different format that shares only the three letters. `.cgm` collides with
  Computer Graphics Metafile. The reader checks the version byte and, when it
  fails, says the file looks like EPR data rather than calling it corrupt.

## On writing proprietary formats

Working agreement §14.7: we do not write formats that are undocumented or that
we cannot verify against a large body of real files. That is why `opus` is
read-only. `spc` is documented — Galactic published the specification — so a
writer would be permissible there; it simply has not been needed yet.

`SPC_Format_Notes.md` at the repository root records the SPC layout, where the
vendor documentation is wrong, and how the reader was validated.
