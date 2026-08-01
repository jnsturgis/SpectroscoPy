# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Processing provenance: what was done to a spectrum, with which parameters.

Every operation that returns a new Spectrum appends a :class:`ProcessingStep`.
Steps record a **structured** (name, params) pair rather than a human-readable
sentence, because that is what lets an interactively-built chain be turned back
into a runnable Pipeline later (roadmap section 2.3):

    ProcessingStep("baseline_correct", {"method": "als", "lam": 1e5, "p": 0.01})

not

    ProcessingStep("baseline corrected via ALS")

Params must be JSON-ish -- numbers, strings, bools, None, and lists/dicts of
those -- so history survives a round trip through .spy or NetCDF. Where an
operation takes another Spectrum as an operand (a reference subtraction, say)
the params record a *reference* to it by name, not the data; replaying such a
step needs the operand supplied from outside.
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
        """Serialise for .spy / JSON. Round-trips through :meth:`from_dict`."""
        return {
            "name": self.name,
            "params": dict(self.params),
            "timestamp": self.timestamp.isoformat(),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ProcessingStep:
        """Rebuild a step from :meth:`to_dict` output."""
        return cls(
            name=data["name"],
            params=dict(data.get("params", {})),
            timestamp=_datetime.datetime.fromisoformat(data["timestamp"]),
        )

    def replace(self, **changes) -> ProcessingStep:
        """Return a copy with some fields changed (steps are frozen)."""
        return replace(self, **changes)
