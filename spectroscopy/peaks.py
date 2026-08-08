# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
PeakTable -- the result type for peak detection.

A dataclass over numpy arrays, not a DataFrame (review section 5.6): a pandas
return type would oblige every caller to have pandas installed to do anything
with the result. :meth:`PeakTable.to_dataframe` is the opt-in escape hatch and
imports pandas inside the method.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ['PeakTable']


@dataclass
class PeakTable:
    """
    Detected peaks in a spectrum.

    Attributes
    ----------
    position : ndarray
        x value of each peak, in the spectrum's x units.
    height : ndarray
        y value of the spectrum at each peak.
    index : ndarray
        Index into the parent spectrum's arrays.
    prominence, width : ndarray or None
        From scipy, when the detection call asked for them.
    kind : str
        ``'peak'``, ``'trough'``, or ``'both'`` for a signed spectrum --
        dichroism, anisotropy, a difference spectrum -- where maxima and
        minima are equally real and one table holds both.
    sign : ndarray or None
        ``+1`` for a maximum, ``-1`` for a minimum, per peak. Always populated
        for a table from :meth:`~spectroscopy.spectra.Spectrum.find_peaks`;
        the only way to tell the two apart when ``kind`` is ``'both'``, since
        a positive band sitting on a negative offset still has a negative
        height.
    x_unit, y_unit : str
        Carried through so a table can be labelled without the spectrum.
    """

    position: np.ndarray
    height: np.ndarray
    index: np.ndarray
    prominence: np.ndarray | None = None
    width: np.ndarray | None = None
    kind: str = 'peak'
    sign: np.ndarray | None = None
    x_unit: str = ''
    y_unit: str = ''
    source: str | None = None
    properties: dict = field(default_factory=dict, repr=False)

    def __len__(self) -> int:
        return len(self.position)

    def __iter__(self):
        """Iterate as (position, height) pairs."""
        return iter(zip(self.position, self.height))

    def __repr__(self) -> str:
        where = f" of {self.source}" if self.source else ""
        if self.kind == 'both' and self.sign is not None:
            up = int(np.sum(np.asarray(self.sign) > 0))
            return (f"<PeakTable: {up} peaks and {len(self) - up} troughs"
                    f"{where}>")
        return f"<PeakTable: {len(self)} {self.kind}s{where}>"

    def __str__(self) -> str:
        if not len(self):
            return f"No {self.kind}s found."
        signed = self.kind == 'both' and self.sign is not None
        header = f"{'position':>12} {'height':>12}"
        lines = [header + ("  kind" if signed else "")]
        for row, (p, h) in enumerate(zip(self.position, self.height)):
            line = f"{p:12.1f} {h:12.5g}"
            if signed:
                line += "  peak" if self.sign[row] > 0 else "  trough"
            lines.append(line)
        return "\n".join(lines)

    # -- selection ---------------------------------------------------------

    def _take(self, selection) -> PeakTable:
        def sub(array):
            return None if array is None else np.asarray(array)[selection]
        return PeakTable(
            position=self.position[selection],
            height=self.height[selection],
            index=self.index[selection],
            prominence=sub(self.prominence),
            width=sub(self.width),
            kind=self.kind, sign=sub(self.sign),
            x_unit=self.x_unit, y_unit=self.y_unit,
            source=self.source,
        )

    def maxima(self) -> PeakTable:
        """Only the bands pointing up. The whole table unless ``kind`` is both."""
        return self._signed(1, 'peak')

    def minima(self) -> PeakTable:
        """Only the bands pointing down."""
        return self._signed(-1, 'trough')

    def _signed(self, wanted, kind) -> PeakTable:
        if self.sign is None:
            return self if self.kind == kind else self._take(
                np.zeros(len(self), dtype=bool))
        selected = self._take(np.asarray(self.sign) == wanted)
        selected.kind = kind
        return selected

    def strongest(self, count) -> PeakTable:
        """
        The ``count`` strongest peaks, **strongest first**.

        Ranked by prominence where scipy supplied it, otherwise by height.
        Call :meth:`sorted_by_position` afterwards for a spectrum-order list;
        the tutorials do exactly that.

        Until 2026-08-05 this returned them in *position* order despite the
        name, so ``strongest(3)`` gave the right three peaks in the wrong
        order and anything reading ``position[0]`` as "the strongest" was
        quietly wrong.

        .. note::

           With the default ``method='second_derivative'``, prominence is
           measured on the second derivative, so it ranks by how *sharp* a
           band is rather than how intense. For intensity, detect with
           ``method='direct'``.
        """
        ranking = self.prominence if self.prominence is not None else self.height
        return self._take(np.argsort(ranking)[::-1][:count])

    def within(self, low, high) -> PeakTable:
        """Peaks whose position lies in [low, high]; order-insensitive."""
        lower, upper = min(low, high), max(low, high)
        return self._take((self.position >= lower) & (self.position <= upper))

    def sorted_by_position(self) -> PeakTable:
        """Peaks in ascending position order."""
        return self._take(np.argsort(self.position))

    # -- interop -----------------------------------------------------------

    def to_dataframe(self):
        """
        As a pandas DataFrame. Requires pandas, which is deliberately not a
        dependency -- see review section 5.6.
        """
        import pandas as pd  # pylint: disable=C0415

        columns = {'position': self.position, 'height': self.height,
                   'index': self.index}
        if self.prominence is not None:
            columns['prominence'] = self.prominence
        if self.width is not None:
            columns['width'] = self.width
        return pd.DataFrame(columns)

    def annotate(self, ax, **kwargs):
        """
        Mark and label these peaks on a matplotlib axes.

        Delegates to :func:`spectroscopy.viz.annotate_peaks`; replaces the
        six-line label loop copy-pasted into about ten notebooks.
        """
        from spectroscopy import viz  # pylint: disable=C0415
        return viz.annotate_peaks(ax, self, **kwargs)
