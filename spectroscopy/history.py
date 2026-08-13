# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The record of what was done to a spectrum.

Every processing step a spectrum goes through appends an entry here, and the
entry keeps the numbers the step actually used -- not a sentence describing it.
A water subtraction records the factor you settled on; a baseline correction
records which algorithm and with what settings.

    ProcessingStep("baseline_correct", {"method": "als", "lam": 1e5, "p": 0.01})

rather than

    ProcessingStep("baseline corrected via ALS")

The difference matters for two reasons. A figure made a year ago can say what
produced it, down to the factor chosen by eye that no one wrote in a notebook.
And because the record is numbers rather than prose, an analysis built up
interactively can later be turned back into something that runs.

The values have to be things that survive being written to a file and read
back -- numbers, text, true and false, and lists or dictionaries of those. A
step that wants to record a whole spectrum records enough to identify it
instead.
"""

from __future__ import annotations

import datetime as _datetime
from dataclasses import dataclass, field, replace
from typing import Any

__all__ = ['ProcessingStep', 'describe_operand']


def describe_operand(operand):
    """
    Render an arithmetic operand for storage in ``params``.

    Scalars are stored as themselves; a Spectrum is stored as a dict naming it,
    since embedding a whole spectrum in the history would make it unbounded.
    """
    if isinstance(operand, (int, float)):
        return {"kind": "scalar", "value": float(operand)}
    name = getattr(operand, "name", None)
    sample = None
    metadata = getattr(operand, "metadata", None)
    if isinstance(metadata, dict):
        sample = metadata.get("sample")
    return {"kind": "spectrum", "name": name, "sample": sample}


@dataclass(frozen=True)
class ProcessingStep:
    """
    One recorded processing operation.

    Attributes
    ----------
    name : str
        The operation, matching the method that produced it -- ``"smooth"``,
        ``"baseline_correct"``, ``"crop"``, ``"subtract"`` ...
    params : dict
        Its arguments, by keyword. Positional/ambiguous forms are normalised
        before being stored, so ``smooth('SG', [2, 15])`` records
        ``{"method": "savgol", "polyorder": 2, "window_length": 15}``.
    timestamp : datetime.datetime
        When it was applied (timezone-aware, UTC).
    """

    name: str
    params: dict[str, Any] = field(default_factory=dict)
    timestamp: _datetime.datetime = field(
        default_factory=lambda: _datetime.datetime.now(_datetime.timezone.utc)
    )

    def __str__(self) -> str:
        arguments = ", ".join(f"{k}={v!r}" for k, v in sorted(self.params.items()))
        return f"{self.name}({arguments})"

    def to_dict(self) -> dict[str, Any]:
        """Serialise for .spy / JSON. Round-trips through ``from_dict``."""
        return {
            "name": self.name,
            "params": dict(self.params),
            "timestamp": self.timestamp.isoformat(),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ProcessingStep:
        """Rebuild a step from ``to_dict`` output."""
        return cls(
            name=data["name"],
            params=dict(data.get("params", {})),
            timestamp=_datetime.datetime.fromisoformat(data["timestamp"]),
        )

    def replace(self, **changes) -> ProcessingStep:
        """Return a copy with some fields changed (steps are frozen)."""
        return replace(self, **changes)
