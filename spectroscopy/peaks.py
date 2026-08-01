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
        ``'peak'`` or ``'trough'``.
    x_unit, y_unit : str
        Carried through so a table can be labelled without the spectrum.
    """

    position: np.ndarray
    height: np.ndarray
    index: np.ndarray
    prominence: np.ndarray | None = None
    width: np.ndarray | None = None
    kind: str = 'peak'
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
        return f"<PeakTable: {len(self)} {self.kind}s{where}>"

    def __str__(self) -> str:
        if not len(self):
            return f"No {self.kind}s found."
        lines = [f"{'position':>12} {'height':>12}"]
        lines += [f"{p:12.1f} {h:12.5g}" for p, h in zip(self.position, self.height)]
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
            kind=self.kind, x_unit=self.x_unit, y_unit=self.y_unit,
            source=self.source,
        )

    def strongest(self, count) -> PeakTable:
        """The ``count`` most prominent peaks (falling back to height)."""
        ranking = self.prominence if self.prominence is not None else self.height
        order = np.argsort(ranking)[::-1][:count]
        return self._take(np.sort(order))

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
