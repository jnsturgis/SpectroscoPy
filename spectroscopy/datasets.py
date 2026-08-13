# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Small example spectra that ship with the package.

Every code example in the documentation runs against these, so someone who has
just installed SpectroscoPy and has no data of their own can still follow the
whole of Getting Started -- copy, paste, run, see the figure. That matters more
than it sounds: "the tutorial doesn't work on my machine" is where most people
stop.

    >>> import spectroscopy as spc
    >>> spectrum = spc.datasets.load('ethanol')
    >>> len(spectrum)
    3570

Use ``available`` to see what there is. These are a deliberately tiny
selection; the fuller collection lives in ``data/`` in the source repository
and is not shipped in the wheel.

``reference_set`` is the exception to "these are examples": SP175 and
SMP180 are the real published CD reference sets, shipped because their terms
allow it, and are meant to be used for actual work rather than only for
following a page.
"""

from __future__ import annotations

import os

__all__ = ['available', 'describe', 'path', 'load', 'load_pair',
           'ftir_replicates', 'replicate_directory', 'emission_series',
           'reference_set', 'DATASETS', 'REFERENCE_SETS']

#: name -> (relative path, technique, one-line description)
DATASETS = {
    'ethanol': (
        'infrared_spectra/ethanol.jdx', 'FTIR',
        'Ethanol vapour-phase IR spectrum (JCAMP-DX, Coblentz Society)'),
    'toluene_uv': (
        'uvvis_spectra/toluene.jdx', 'UV-Vis',
        'Toluene UV absorption spectrum (JCAMP-DX)'),
    'tannic_acid': (
        'raman_spectra/tannic_acid.jdx', 'Raman',
        'Tannic acid Raman spectrum (JCAMP-DX)'),
    'uvvis_1': (
        'uvvis_spectra/Spectrum1.csv', 'UV-Vis',
        'Bacterial membrane fraction, UV-Vis (headerless CSV)'),
    'uvvis_2': (
        'uvvis_spectra/Spectrum2.csv', 'UV-Vis',
        'Second membrane fraction, UV-Vis (headerless CSV)'),
    'aqpz': (
        'cd_spectra/aqpz_w14a.spy', 'CD',
        'AqpZ-W14A far-UV CD, 4 uM, 30 C (JASCO J-815, HT < 600 V)'),
}

#: The published CD reference sets that ship, ``name -> (directory, citation)``.
#:
#: **These are somebody else's data and they ship because their terms allow
#: it.** The `pcddb organisation <https://github.com/pcddb/DichroWebGit>`_
#: publishes them under the MIT licence, © 2023 Andy Miles, which permits
#: redistribution provided the notice travels with them --
#: ``LICENSE.DichroWebGit`` sits beside the data and is installed with it.
#: Taken from the PCDDB website instead they would carry no such grant: those
#: terms give access and say nothing about reuse (ADR-0002 §9).
#:
#: Citation is a condition of use of the underlying data, so it travels in the
#: set's ``info`` rather than in a docstring nobody reads.
#:
#: .. warning::
#:
#:    The volume and page numbers in ``docs/references.md`` are marked
#:    **unverified** -- written from memory and not yet checked against the
#:    publishers. So what travels with the data here is author, year and
#:    journal, which are what identify the paper, and the reader is sent to the
#:    reference list rather than being handed page numbers this package cannot
#:    vouch for. Meeting a citation condition with a citation that might be
#:    wrong is worse than not printing one.
REFERENCE_SETS = {
    'sp175': (
        'cd_reference/sp175',
        'Lees, Miles, Wien & Wallace (2006), Bioinformatics -- SP175. '
        'Full reference: docs/references.md'),
    'smp180': (
        'cd_reference/smp180',
        'Abdul-Gader, Miles & Wallace (2011), Bioinformatics -- SMP180. '
        'Full reference: docs/references.md'),
}


def _root():
    """
    Where the sample files live.

    Installed, they sit inside the package; in a source checkout they are in
    ``data/`` at the top level. Both are supported so the documentation builds
    from a checkout and runs from an install.
    """
    here = os.path.dirname(os.path.abspath(__file__))
    packaged = os.path.join(here, 'data')
    if os.path.isdir(packaged):
        return packaged
    return os.path.join(os.path.dirname(here), 'data')


def available():
    """The names that can be passed to ``load``."""
    return tuple(sorted(DATASETS))


def describe():
    """A readable table of the example spectra -- handy in a notebook."""
    lines = [f"{'name':<14} {'technique':<10} description"]
    for name in available():
        _, technique, description = DATASETS[name]
        lines.append(f"{name:<14} {technique:<10} {description}")
    return "\n".join(lines)


def path(name):
    """Filesystem path of an example file, without reading it."""
    if name not in DATASETS:
        raise KeyError(
            f"No example dataset {name!r}; available: {', '.join(available())}"
        )
    relative, _, _ = DATASETS[name]
    full = os.path.join(_root(), relative)
    if not os.path.exists(full):
        raise FileNotFoundError(
            f"Example data for {name!r} is missing (looked in {full}). If you "
            f"are running from a source checkout, the files are in data/."
        )
    return full


def load(name):
    """
    Load an example spectrum, with its technique already set.

    Parameters
    ----------
    name : str
        One of ``available``.

    Returns
    -------
    Spectrum
    """
    from spectroscopy.io import read_spectrum  # pylint: disable=C0415

    _, technique, _ = DATASETS[name]
    spectrum = read_spectrum(path(name))
    if technique:
        spectrum.set_type(technique)
    # Only name it after the key when the file did not name itself. Most of
    # these are JCAMP or bare CSV and cannot; a .spy states its own name and
    # sample, and overwriting them here would contradict what the rest of the
    # library promises -- that a file's own statement beats any default.
    if spectrum.name in (None, '', 'unnamed') or spectrum.name.endswith(
            ('.jdx', '.csv', '.dx', '.dpt')):
        spectrum.name = name
    if not spectrum.metadata.get('sample'):
        spectrum.set_sample(name)
    return spectrum


def replicate_directory():
    """
    Folder holding the ATR-FTIR replicate files, for glob-based loading.

        >>> import spectroscopy as spc
        >>> folder = spc.datasets.replicate_directory()
        >>> spectra = spc.SpectrumCollection.from_files(folder + "/*.dpt")
    """
    folder = os.path.join(_root(), 'ftir_replicates')
    if not os.path.isdir(folder):
        raise FileNotFoundError(f"Example replicates are missing ({folder})")
    return folder


def ftir_replicates():
    """
    Nine real ATR-FTIR spectra: three replicates each of glucose, cellulose
    and water.

    Enough to demonstrate the whole workflow the library exists for -- group by
    sample, average the replicates, subtract the water contribution, baseline
    correct, normalise, pick peaks -- without needing a folder of your own.

    Returns
    -------
    SpectrumCollection
    """
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415

    return SpectrumCollection.from_files(
        os.path.join(replicate_directory(), '*.dpt'), technique='ATR-FTIR')


def emission_series():
    """
    A fluorescence excitation-emission series: 18 emission spectra, one per
    excitation wavelength from 290 to 455 nm.

    Recorded by Chloe (Sturgis group) on a candidate flavoprotein, and used
    with permission. The file is a wide export with paired (wavelength,
    intensity) columns, which is what makes it a good demonstration of the
    generic table reader.

    Returns
    -------
    SpectrumCollection
    """
    from spectroscopy.io import read_spectra  # pylint: disable=C0415

    path_ = os.path.join(_root(), 'fluorescence', 'J_Peri.csv')
    if not os.path.exists(path_):
        raise FileNotFoundError(f"Example emission series is missing ({path_})")

    series = read_spectra(path_, 'table', paired=True)
    wanted = [s for s in series if '_EX_' in s.name]
    for spectrum in wanted:
        spectrum.set_type('Fluorescence')
        spectrum.metadata['excitation_nm'] = float(
            spectrum.name.rsplit('_', 1)[-1])
        spectrum.name = f"ex {spectrum.metadata['excitation_nm']:.0f} nm"
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415
    return SpectrumCollection(wanted, name='J-peri emission series')


def reference_set(name='smp180'):
    """
    A published CD reference set, as a
    ``ReferenceSet``.

    ``'sp175'``
        71 soluble proteins, 175-240 nm. The standard set for a soluble
        protein, and the one most published comparisons use.
    ``'smp180'``
        128 proteins, 180-240 nm, **30 of them membrane proteins**. Use this
        one if your protein is in detergent or lipid: a membrane protein's
        far-UV spectrum resembles other membrane proteins about twice as often
        as its structure alone would explain, so a set with none in it has
        nothing close to yours.

    Both are in delta epsilon per residue, and both carry their licence and
    citation in ``info`` -- citation is a condition of use of the underlying
    data, so it travels with the set rather than living in a docstring.

    Deliberately **not** ``load()``: that returns one
    ``Spectrum``, and a function whose return type
    depends on the string you pass it is the kind of thing the guessability
    audit exists to remove.

        >>> import spectroscopy as spc
        >>> references = spc.datasets.reference_set('sp175')   # doctest: +SKIP
        >>> len(references)                                    # doctest: +SKIP
        71

    Returns
    -------
    ReferenceSet
    """
    from spectroscopy.library import load_dichroweb_basis  # pylint: disable=C0415

    if name not in REFERENCE_SETS:
        raise KeyError(
            f"No reference set {name!r}; available: "
            f"{', '.join(sorted(REFERENCE_SETS))}"
        )
    relative, citation = REFERENCE_SETS[name]
    folder = os.path.join(_root(), relative)
    if not os.path.isdir(folder):
        raise FileNotFoundError(
            f"Reference set {name!r} is missing (looked in {folder}). If you "
            f"are running from a source checkout, it is in data/cd_reference/."
        )
    references = load_dichroweb_basis(folder)
    references.info['citation'] = citation
    return references


def load_pair():
    """
    Two comparable UV-Vis spectra, as a ``SpectrumCollection``.

    Enough to demonstrate averaging, arithmetic and overlay plotting without
    needing a folder of your own replicates.
    """
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415

    return SpectrumCollection([load('uvvis_1'), load('uvvis_2')],
                              name='example pair')
