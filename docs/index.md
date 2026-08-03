# SpectroscoPy

One way of handling spectra, whatever instrument they came from.

Spectroscopy data arrives in a different file format from every machine, and
each comes with its own idiosyncratic analysis program. SpectroscoPy reads them
all into one object, gives you the operations you actually use — baseline,
normalise, average replicates, find peaks, subtract a reference — and keeps a
record of everything it did, so a figure can always be traced back to the raw
files.

```python
import spectroscopy as spc

spectra = spc.SpectrumCollection.from_files("data/*.dpt", technique="ATR-FTIR")

result = (spectra.group_by('sample')['PG_coli']
          .mean()
          .crop(900, 1800)
          .baseline_correct('rubberband')
          .normalize('max', window=(1050, 1080)))

peaks = result.find_peaks(prominence=0.05, relative=True)
```

**New here? Start with [Getting started](getting-started.md).** If you know what
you want to find out but not how, try
[What do you want to know?](what-do-you-want-to-know.md). It goes from
`pip install` to a labelled figure, and assumes only that you can copy text
into a terminal.

:::{admonition} Pre-1.0
:class: warning

The API is still settling and may change between 0.x releases. It is usable
for real work — it is used for real work — but pin a version if you depend
on it.
:::

## Supported formats

| Format | Extensions | Notes |
|---|---|---|
| **Bruker OPUS native** | `.0` … `.9` | read-only; keeps what OPUS already did to the spectrum |
| **Galactic / Thermo SPC** | `.spc` `.cgm` | read-only; single and multi-subfile |
| Bruker OPUS data point table | `.dpt` | separator detected per file |
| JCAMP-DX | `.jdx` `.dx` `.jcamp` | including compound files |
| Delimited text | `.csv` `.tsv` `.txt` | separator and decimal comma detected |
| Wide / paired tables | any | many spectra from one file |
| SpectroscoPy native | `.spy` | keeps units and processing history |

The two binary formats are read-only on purpose: SpectroscoPy will not write a
format it cannot verify it is writing correctly.

Techniques with built-in axis conventions: ATR-FTIR, FTIR, UV-Vis, Raman,
fluorescence.

```{toctree}
:maxdepth: 2
:hidden:

getting-started
what-do-you-want-to-know
tutorials/index
guide/index
api
adr/index
references
contributing
```
