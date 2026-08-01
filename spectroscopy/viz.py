# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Plotting -- the figure cells of the notebooks, as functions.

This is the ``viz`` layer of the roadmap's io -> core -> processing -> viz
split. ``Spectrum.plot`` used to put matplotlib inside the data model (review
crossing C2); it now delegates here.

What is here is what the notebooks actually draw, repeatedly and by hand:

======================  ==================================================
:func:`plot`            one spectrum, with the axis labels and the reversed
                        x axis that FTIR and Raman want
:func:`plot_collection` several overlaid, with a legend
:func:`stack`           offset traces -- the Figure S5 panel A layout
:func:`grid`            one panel per sample
:func:`annotate_peaks`  the six-line label loop from ~10 notebooks
:func:`annotate_bands`  ``{1650: "Amide I", ...}`` as marked assignments
:func:`plot_baseline`   a spectrum with its baseline and the correction
:func:`plot_decomposition`  NMF/PCA components, a fit, and the residuals
:func:`plot_scores`     samples in component space
======================  ==================================================

Colour
------
The default cycle is the six chromatic Okabe-Ito colours -- the standard
colour-vision-deficiency-safe qualitative set in scientific publishing --
ordered so that adjacent pairs stay far apart under deuteranopia. The order was
checked with a palette validator rather than by eye:

    lightness band      PASS
    chroma floor        PASS
    CVD separation      PASS  worst adjacent pair dE 9.6 (deutan)
    normal-vision floor PASS  worst adjacent pair dE 20.0
    contrast vs white   WARN  three hues below 3:1

The contrast warning is why a legend is drawn whenever there is more than one
series, and why :func:`stack` direct-labels its traces: identity is never
carried by colour alone. Okabe-Ito's yellow and black are deliberately excluded
-- the yellow is too light against white and the black has no chroma.

Beyond six series the colours repeat but the line style advances, so the
(colour, style) pair stays unique to 24 traces. Past a handful of overlaid
spectra :func:`stack` or :func:`grid` reads better than any palette.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    'PALETTE', 'LINE_STYLES', 'FRAME_AXES', 'apply_axes', 'set_frame',
    'series_style',
    'plot', 'plot_collection', 'stack', 'grid',
    'annotate_peaks', 'annotate_bands', 'plot_baseline',
    'plot_decomposition', 'plot_scores',
]

#: Fixed order -- assign from the start of the list, never shuffle. See above.
PALETTE = ('#0072B2', '#E69F00', '#009E73', '#D55E00', '#56B4E9', '#CC79A7')

#: Advanced when the colours wrap, so a 7th series is not a repeat of the 1st.
LINE_STYLES = ('-', '--', '-.', ':')

DEFAULT_LINEWIDTH = 1.2

#: Draw a full box (all four spines) rather than just left and bottom.
#: Physical-chemistry journals overwhelmingly frame their axes; a lot of modern
#: plotting advice does not. Set ``viz.FRAME_AXES = True`` once and every figure
#: from then on is framed, or pass ``frame=True`` to a single call.
FRAME_AXES = False


def _pyplot():
    """Imported lazily: a library import should not choose a backend."""
    import matplotlib.pyplot as plt  # pylint: disable=C0415
    return plt


def _axes(ax):
    return ax if ax is not None else _pyplot().subplots()[1]


def series_style(index):
    """Colour and line style for series ``index``, counting from 0."""
    return {'color': PALETTE[index % len(PALETTE)],
            'linestyle': LINE_STYLES[(index // len(PALETTE)) % len(LINE_STYLES)]}


def apply_axes(ax, spectrum, reverse=None, frame=None):
    """
    Label the axes from the spectrum, and reverse x where convention wants it.

    FTIR and Raman are quoted high-to-low; every figure cell in the notebooks
    sets ``ax.set_xlim((1800, 900))`` and the wavenumber label by hand.

    Parameters
    ----------
    frame : bool, optional
        Draw a full box around the plot instead of only the left and bottom
        spines. Defaults to :data:`FRAME_AXES`, which you can set once at the
        top of a session rather than passing this everywhere. Ticks are drawn
        inwards on all four sides when framed, which is the usual convention
        that goes with it.
    """
    ax.set_xlabel(spectrum.x_label)
    ax.set_ylabel(spectrum.y_label)
    should_reverse = spectrum.reversed_x if reverse is None else reverse
    if should_reverse and not ax.xaxis_inverted():
        ax.invert_xaxis()
    set_frame(ax, FRAME_AXES if frame is None else frame)
    return ax


def set_frame(ax, framed=True):
    """
    Turn the box around a plot on or off.

    Useful on an axes this module did not draw -- an inset, or a panel you
    built yourself.
    """
    for side in ('top', 'right'):
        ax.spines[side].set_visible(bool(framed))
    if framed:
        ax.tick_params(direction='in', top=True, right=True)
    else:
        ax.tick_params(direction='out', top=False, right=False)
    return ax


def plot(spectrum, ax=None, *args, label=None, apply_labels=True, frame=None,
         **kwargs):
    """
    Plot one spectrum. Returns the matplotlib line list.

    ``frame=True`` draws a full box around the plot; see :data:`FRAME_AXES` to
    make that the default for a whole session.
    """
    ax = _axes(ax)
    kwargs.setdefault('linewidth', DEFAULT_LINEWIDTH)
    lines = ax.plot(spectrum.x, spectrum.y, *args,
                    label=spectrum.name if label is None else label, **kwargs)
    if apply_labels:
        apply_axes(ax, spectrum, frame=frame)
    return lines


def _labels_for(collection, labels):
    if labels is not None:
        return list(labels)
    return [spectrum.name for spectrum in collection]


def plot_collection(collection, ax=None, labels=None, colors=None,
                    legend=None, frame=None, **kwargs):
    """
    Overlay a collection on one axes.

    A legend appears whenever there is more than one series and never for a
    single one -- the title or the axis already names that.
    """
    ax = _axes(ax)
    names = _labels_for(collection, labels)
    kwargs.setdefault('linewidth', DEFAULT_LINEWIDTH)

    for index, (spectrum, name) in enumerate(zip(collection, names)):
        style = series_style(index)
        if colors is not None:
            style['color'] = colors[index % len(colors)]
        plot(spectrum, ax, label=name, apply_labels=(index == 0),
             frame=frame, **{**style, **kwargs})

    if legend is None:
        legend = len(collection) > 1
    if legend:
        ax.legend(frameon=False, fontsize='small')
    return ax


def stack(collection, ax=None, offsets=None, gap=None, labels=None,
          colors=None, direct_labels=True, frame=None, **kwargs):
    """
    Offset traces vertically -- the stacked-spectra figure.

    ``offsets`` gives an explicit shift per spectrum; otherwise they are spaced
    by ``gap`` (default: a little more than the tallest trace, so nothing
    collides). Traces are labelled at their right-hand end rather than in a
    legend box, so identity does not depend on matching colours to a key.
    """
    ax = _axes(ax)
    names = _labels_for(collection, labels)
    kwargs.setdefault('linewidth', DEFAULT_LINEWIDTH)

    # Set the axis direction first: the direct labels below ask the axes which
    # end is on the right, and inverting afterwards would put them all on the
    # wrong side.
    if len(collection):
        apply_axes(ax, collection[0], frame=frame)

    if offsets is None:
        if gap is None:
            spans = [float(np.nanmax(s.y) - np.nanmin(s.y)) for s in collection]
            gap = 1.05 * max(spans) if spans else 1.0
        offsets = [index * gap for index in range(len(collection))]

    for index, (spectrum, name, offset) in enumerate(
            zip(collection, names, offsets)):
        style = series_style(index)
        if colors is not None:
            style['color'] = colors[index % len(colors)]
        ax.plot(spectrum.x, spectrum.y + offset, label=name,
                **{**style, **kwargs})
        if direct_labels:
            # Anchor the label to whichever end of the trace is drawn on the
            # right, and to that end's own y value. Using x[-1] would assume
            # the data is stored in ascending order, which plenty of vendor
            # files are not -- UV-Vis exports commonly run high to low, and the
            # labels then land on the wrong side of the figure.
            end = (int(np.argmin(spectrum.x)) if ax.xaxis_inverted()
                   else int(np.argmax(spectrum.x)))
            ax.annotate(name, xy=(spectrum.x[end], spectrum.y[end] + offset),
                        xytext=(4, 0), textcoords='offset points',
                        va='center', fontsize='small', annotation_clip=False)

    ax.set_yticks([])
    if not (FRAME_AXES if frame is None else frame):
        # The y scale is arbitrary once traces are offset, so the left spine
        # is noise -- unless a full box was asked for, where it is the point.
        ax.spines['left'].set_visible(False)
    else:
        ax.tick_params(left=False, right=False)
    if not direct_labels and len(collection) > 1:
        ax.legend(frameon=False, fontsize='small')
    return ax


def grid(collection, key='sample', ncols=2, figsize=None, sharex=True,
         frame=None, **kwargs):
    """
    One panel per group -- the overview figure at the top of every notebook.

    ``key`` is a metadata field or a callable, as for
    :meth:`SpectrumCollection.group_by`. Returns ``(figure, axes)``.
    """
    plt = _pyplot()
    groups = collection.group_by(key)
    count = len(groups)
    nrows = int(np.ceil(count / ncols)) if count else 1

    figure, axes = plt.subplots(nrows, ncols, sharex=sharex,
                                figsize=figsize or (6 * ncols, 2.2 * nrows),
                                squeeze=False)
    flat = axes.ravel()
    for panel, (name, group) in zip(flat, groups.items()):
        for index, spectrum in enumerate(group):
            plot(spectrum, panel, apply_labels=(index == 0), label=None,
                 frame=frame, **{**series_style(index), **kwargs})
        panel.set_title(str(name), fontsize='small', loc='left')
    for panel in flat[count:]:
        panel.set_visible(False)
    figure.tight_layout()
    return figure, axes


def annotate_peaks(ax, peaks, fmt="{:.0f}", offset=0.0, rotation=270,
                   marker='o', markersize=4, color=None, fontsize=8,
                   label_peaks=True):
    """
    Mark and label peaks.

    Replaces the loop pasted into about ten notebooks::

        ax.plot(x[peaks], y[peaks], "ro", ms=4)
        for i in range(len(peaks)):
            ax.text(x[peaks[i]], y[peaks[i]] + dy, f"{x[peaks[i]]:.0f}", ...)
    """
    ax.plot(peaks.position, peaks.height, marker, ms=markersize,
            color=color or '#444444', linestyle='none')
    if label_peaks:
        for position, height in zip(peaks.position, peaks.height):
            ax.text(position, height + offset, fmt.format(position),
                    ha='center', va='bottom', fontsize=fontsize,
                    rotation=rotation)
    return ax


def annotate_bands(ax, bands, y=None, color='#999999', fontsize=8,
                   linestyle=':', linewidth=0.8, rotation=90):
    """
    Mark known band assignments: ``{1650: "Amide I", 1550: "Amide II"}``.

    From the peptidoglycan figure in the Chloe notebook, where the assignment
    table is the point of the figure.
    """
    if y is None:
        y = ax.get_ylim()[1]
    for position, name in bands.items():
        ax.axvline(x=position, color=color, linestyle=linestyle,
                   linewidth=linewidth, zorder=0)
        ax.text(position, y, str(name), rotation=rotation, va='bottom',
                ha='center', fontsize=fontsize)
    return ax


def plot_baseline(spectrum, baseline, ax=None, corrected=True, **kwargs):
    """
    A spectrum with its baseline, and optionally the corrected result.

    Worth looking at every time: a rubberband baseline that has cut a band is
    obvious here and invisible in the corrected spectrum alone.
    """
    ax = _axes(ax)
    kwargs.setdefault('linewidth', DEFAULT_LINEWIDTH)
    plot(spectrum, ax, label='measured', color=PALETTE[0], **kwargs)
    ax.plot(baseline.x, baseline.y, label='baseline', color=PALETTE[1],
            linestyle='--', **kwargs)
    if corrected:
        ax.plot(spectrum.x, spectrum.y - baseline.y, label='corrected',
                color=PALETTE[2], **kwargs)
    ax.legend(frameon=False, fontsize='small')
    return ax


def plot_decomposition(result, sample=0, figsize=(7.5, 5.5)):
    """
    The three-panel decomposition figure: components, one fit, all residuals.

    This is the layout built by hand in ``Biofilm_CK_241125`` and
    ``Figure_CK_111225``. Returns ``(figure, (components, fit, residuals))``.
    """
    plt = _pyplot()
    figure = plt.figure(figsize=figsize)
    spec = figure.add_gridspec(2, 2)
    top = figure.add_subplot(spec[0, :])
    left = figure.add_subplot(spec[1, 0])
    right = figure.add_subplot(spec[1, 1])

    plot_collection(result.components, top)
    top.set_title(f"{result.method.upper()} components (k={result.n_components})",
                  fontsize='small', loc='left')

    observed = result.components[0]
    left.plot(result.x, result._data[sample] if result._data is not None      # noqa: SLF001
              else result.reconstruction[sample],
              color='#333333', lw=DEFAULT_LINEWIDTH, label='observed')
    left.plot(result.x, result.reconstruction[sample], color=PALETTE[3],
              lw=DEFAULT_LINEWIDTH, label='reconstructed')
    left.set_title(f"{result.sample_names[sample]}", fontsize='small', loc='left')
    left.legend(frameon=False, fontsize='small')
    apply_axes(left, observed)

    for residual in result.residuals:
        right.plot(result.x, residual, lw=0.7, color='#888888', alpha=0.7)
    right.set_title("residuals", fontsize='small', loc='left')
    apply_axes(right, observed)

    figure.tight_layout()
    return figure, (top, left, right)


def plot_scores(result, ax=None, components=(0, 1), labels=None, annotate=True):
    """
    Samples in component space -- the PCA/NMF scores scatter.

    Points are labelled directly rather than coloured by group, so the figure
    does not need a key and does not depend on colour to identify a sample.
    """
    ax = _axes(ax)
    first, second = components
    scores = result.weights
    names = labels if labels is not None else result.sample_names

    ax.scatter(scores[:, first], scores[:, second], s=42,
               color=PALETTE[0], zorder=3)
    if annotate:
        for index, name in enumerate(names):
            ax.annotate(str(name), (scores[index, first], scores[index, second]),
                        xytext=(4, 3), textcoords='offset points', fontsize=8)

    ratios = result.explained_variance_ratio
    def _axis_label(index):
        if ratios is not None and index < len(ratios):
            return f"Component {index + 1} ({100 * ratios[index]:.1f}%)"
        return f"Component {index + 1}"

    ax.set_xlabel(_axis_label(first))
    ax.set_ylabel(_axis_label(second))
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    return ax
