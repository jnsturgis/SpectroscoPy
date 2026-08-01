# Plotting

`spectroscopy.viz` holds the drawing. `Spectrum.plot`, `SpectrumCollection.plot`
and `PeakTable.annotate` delegate to it, so the short forms still work.

## One spectrum

```python
fig, ax = plt.subplots()
spectrum.plot(ax)
```

The axis labels come from the spectrum, and the x axis is reversed for the
techniques that are conventionally drawn high-to-low — FTIR and Raman. Those
are the two lines typed by hand in every figure cell of a notebook.

`frame=True` draws a full box with inward ticks, as most physical-chemistry
journals expect. To make it the default for a session:

```python
from spectroscopy import viz
viz.FRAME_AXES = True
```

## Several

```python
viz.plot_collection(spectra, ax)              # overlaid, with a legend
viz.stack(spectra, ax)                        # offset, labelled directly
fig, axes = viz.grid(spectra, key='sample')   # one panel per sample
```

Overlay works to about half a dozen. Past that, `stack` or `grid` reads better
than any palette can.

`grid` is worth running on every new dataset before doing anything else: a
replicate that went wrong is obvious in a panel-per-sample overview and very
hard to spot after averaging.

## Annotation

```python
peaks.strongest(6).annotate(ax, offset=0.03)
viz.annotate_bands(ax, {1650: "Amide I", 1550: "Amide II"})
viz.plot_baseline(spectrum, spectrum.baseline('rubberband'), ax)
```

Label only the peaks you intend to discuss. Labelling every peak on a crowded
spectrum produces an overlapping smear that is worse than no labels — use
`strongest(n)` and `distance=` in the detection.

## Colour

The default cycle is the six chromatic Okabe–Ito colours, ordered so adjacent
pairs stay separated under deuteranopia. The order was checked with a palette
validator rather than by eye; the obvious sequence puts pink next to green and
fails. Okabe–Ito's yellow and black are excluded — the yellow is too light
against white, the black has no chroma.

Colours are assigned in fixed order from the start of the list, never shuffled,
so a series keeps its colour when the set it is drawn with changes. Past six
the colours repeat but the line style advances, so the pair stays unique.

Three of the six sit below 3:1 contrast against white. That is why a legend is
always drawn for two or more series and why `stack` labels its traces directly:
identity never rests on colour alone.

For a **map** — an excitation-emission surface, say — use a sequential ramp
(one hue, light to dark, like `viridis`). A categorical palette encodes
identity; a map encodes magnitude. Never a rainbow: it invents boundaries where
the data is smooth and is unreadable in greyscale.

## Decomposition figures

```python
fig, axes = viz.plot_decomposition(result)   # components, one fit, residuals
viz.plot_scores(result, ax)                  # samples in component space
```

## Saving

Ordinary matplotlib:

```python
fig.savefig("figure.pdf")      # vector, for a paper
fig.savefig("figure.png", dpi=300)
```

Nothing in this module prevents you doing anything matplotlib can do — every
function takes and returns axes, so drop down whenever you need to.
