# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The agreed metadata keys, and what they mean.

``Spectrum.metadata`` is a plain dictionary and anything may be put in it. But
``.spy`` serialises it verbatim, so **the names of the keys the library itself
reads are part of the file format**, and they freeze when the format does. A
key that is invented independently in three modules is three keys.

Some of this had already happened by accident. ``path_length`` was introduced
by :func:`~spectroscopy.processing.unmix.unmix` and read again by
:mod:`~spectroscopy.processing.scattering` without ever being written down;
``parameter`` arrived with continuous series. Circular dichroism needs
concentration and a residue count for its central conversion, and a redox
titration will want a temperature and a reference electrode. That is four
modules converging on the same handful of sample facts, which is the point at
which they should be agreed once rather than five times.

Four groups, and the difference matters:

**Sample conditions** -- what was in the cuvette and how it was measured.
These are the ones an analysis reads back to turn a signal into a quantity,
and getting them wrong is the expensive kind of error, because the arithmetic
succeeds either way.

**Identification** -- what this spectrum is of, and its place in a series.

**Acquisition** -- instrument settings general enough to deserve a shared
name rather than a per-reader one.

**Provenance** -- what the instrument said. Free-form by nature; recorded here
so that the reserved prefixes are known and readers do not collide.

Nothing here is enforced. A missing key is normal and every consumer decides
what to do about it -- but the rule for the sample conditions is that a
consumer **states its assumption or refuses**, never silently supplies a
default that changes the answer. A path length assumed to be 1 cm when it was
1 mm makes every concentration ten times too small, and nothing about the
result looks wrong.
"""

from __future__ import annotations

__all__ = ['SAMPLE_CONDITIONS', 'IDENTIFICATION', 'ACQUISITION',
           'PROVENANCE_PREFIXES',
           'KNOWN_KEYS', 'describe', 'unknown_keys']


#: Physical conditions of the measurement. ``(key, unit, meaning)``.
#:
#: The unit is fixed rather than carried alongside, deliberately: a value and
#: a unit stored as two keys can disagree, and the pair that disagrees is
#: indistinguishable from the pair that does not. Where a quantity genuinely
#: needs a caller-chosen unit -- a series parameter can be mV or degrees or
#: micromolar -- the unit is a separate labelling key and the value has no
#: fixed meaning without it. That is why ``parameter`` sits in IDENTIFICATION
#: and not here.
SAMPLE_CONDITIONS = (
    ('path_length',   'cm',
     "Optical path length of the cell. Read by unmix() and by any "
     "Beer-Lambert calculation. Absent means unknown, not 1."),
    ('concentration', 'M',
     "Molar concentration of the species the spectrum is of. Needed for an "
     "extinction coefficient, and for CD's mean residue ellipticity."),
    ('mass_concentration', 'mg/mL',
     "Where a molar concentration is not known -- proteins quoted by mass, "
     "nucleic acids in ug/mL (which is 0.001 of this)."),
    ('n_residues',    'count',
     "Residues per chain. Needed to get from ellipticity to *mean residue* "
     "ellipticity; a molar quantity divided by this is a per-residue one."),
    ('mean_residue_weight', 'Da',
     "Alternative to n_residues, and what CD literature usually quotes: "
     "molecular weight divided by residue count, typically near 110."),
    ('temperature',   'C',
     "Sample temperature. A melt varies it; a Nernst or van 't Hoff fit "
     "needs it, and using the wrong one returns a plausible wrong answer."),
    ('pH',            'pH',
     "Sample pH."),
    ('reference_electrode', 'name',
     "What a potential was measured against -- 'SHE', 'Ag/AgCl'. These "
     "differ by about 200 mV, and no analysis can recover which was meant."),
)

#: What the spectrum is of, and where it sits in a series.
IDENTIFICATION = (
    ('sample',    'name',  "What this is a spectrum of. Categorical; "
                           "group_by() groups on it."),
    ('reference', 'name',  "What it was measured against, or the reference "
                           "subtracted from it."),
    ('spec_type', 'name',  "The technique, as set by set_type()."),
    ('parameter', 'number', "The continuous condition this spectrum was "
                            "measured at -- see Spectrum.set_parameter()."),
    ('parameter_name', 'name', "What that number is: 'potential', "
                               "'temperature', 'concentration'."),
    ('parameter_unit', 'name', "What it is in: 'mV', 'C', 'uM'."),
)

#: Instrument settings that are generic rather than vendor-specific, so they
#: are worth a shared name instead of a prefix per reader.
#:
#: .. warning::
#:
#:    ``excitation_nm`` and ``z_value`` are **the same idea as** ``parameter``
#:    above: each is the continuous condition that distinguishes one spectrum
#:    in a series from the next -- excitation wavelength for an EEM, subfile z
#:    for a multi-subfile SPC, potential for a titration. They arrived
#:    separately, from a dataset helper and a binary reader, before
#:    ``parameter`` existed.
#:
#:    Three spellings of one concept is how a format ends up with three, so
#:    this needs deciding before the freeze rather than after. The specific
#:    names are worth keeping -- an excitation wavelength really is an
#:    excitation wavelength -- but a reader that sets one should arguably set
#:    ``parameter`` as well, so that ``sorted_by_parameter()`` and
#:    ``to_matrix(with_parameter=True)`` work on an EEM and on an SPC series
#:    without the caller knowing which reader produced them.
ACQUISITION = (
    ('excitation_nm', 'nm',
     "Excitation wavelength of an emission scan. The parameter of an EEM."),
    ('z_value',    'varies',
     "Third-axis value of one spectrum in a multi-spectrum file, in "
     "z_quantity's units. The parameter of an SPC subfile series."),
    ('z_quantity', 'name',
     "What z_value measures, as the file declared it."),
    ('scans',      'count',
     "Accumulations averaged into this spectrum. Bears on the noise, and so "
     "on whether two spectra are comparable."),
)

#: Reader-specific provenance. Keys beginning with these are the reader's own
#: and no analysis should depend on them; they exist so that what the
#: instrument said is not lost.
PROVENANCE_PREFIXES = ('opus_', 'spc_', 'jcamp_', 'file_')

#: Every key this library reads by name.
KNOWN_KEYS = frozenset(
    [key for key, _, _ in SAMPLE_CONDITIONS]
    + [key for key, _, _ in IDENTIFICATION]
    + [key for key, _, _ in ACQUISITION]
)


def describe() -> str:
    """The schema as a table, for a docstring or a terminal."""
    lines = []
    for title, group in (("Sample conditions", SAMPLE_CONDITIONS),
                         ("Identification", IDENTIFICATION),
                         ("Acquisition", ACQUISITION)):
        lines.append(f"{title}:")
        for key, unit, meaning in group:
            lines.append(f"  {key:<22} [{unit}]  {meaning}")
        lines.append("")
    lines.append("Provenance prefixes: " + ", ".join(PROVENANCE_PREFIXES))
    return "\n".join(lines)


def unknown_keys(spectrum_or_metadata):
    """
    Keys that are neither in the schema nor reader provenance.

    Not a validator -- extra keys are perfectly legitimate, and this returns
    them rather than complaining about them. It is for finding the near miss:
    ``pathlength`` where ``path_length`` was meant is silently ignored by
    every consumer, and looks identical to not having set it.
    """
    metadata = getattr(spectrum_or_metadata, 'metadata', spectrum_or_metadata)
    return sorted(
        key for key in metadata
        if key not in KNOWN_KEYS
        and not key.startswith(PROVENANCE_PREFIXES)
    )
