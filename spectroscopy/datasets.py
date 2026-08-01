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

Use :func:`available` to see what there is. These are a deliberately tiny
selection (about 100 kB); the fuller collection lives in ``data/`` in the
source repository and is not shipped in the wheel.
"""

from __future__ import annotations

import os

__all__ = ['available', 'describe', 'path', 'load', 'load_pair',
           'ftir_replicates', 'replicate_directory', 'emission_series',
           'DATASETS']

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
    """The names that can be passed to :func:`load`."""
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
        One of :func:`available`.

    Returns
    -------
    Spectrum
    """
    from spectroscopy.io import read_spectrum  # pylint: disable=C0415

    _, technique, _ = DATASETS[name]
    spectrum = read_spectrum(path(name))
    if technique:
        spectrum.set_type(technique)
    spectrum.name = name
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


def load_pair():
    """
    Two comparable UV-Vis spectra, as a :class:`SpectrumCollection`.

    Enough to demonstrate averaging, arithmetic and overlay plotting without
    needing a folder of your own replicates.
    """
    from spectroscopy.collection import SpectrumCollection  # pylint: disable=C0415

    return SpectrumCollection([load('uvvis_1'), load('uvvis_2')],
                              name='example pair')
