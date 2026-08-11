# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Processing and analysis algorithms.

``cd``
    Secondary structure from circular dichroism: four published or
    principled methods, and held-out validation to choose between them.
``common``
    Technique-agnostic: baseline, smoothing, normalisation, peak detection.
    Every one is also a :class:`~spectroscopy.spectra.Spectrum` method.
``ftir``
    Amide I side-chain contributions and residual spectra.
``melting``
    Two-state thermal unfolding: Tm, van 't Hoff enthalpy, and two checks
    that the transition really is two-state.
``multivariate``
    PCA / NMF / ICA across a collection, with bootstrap stability.
``scattering``
    Removing a scattering background before anything is quantified.
``structure``
    Secondary structure, as a :class:`~.structure.Composition`.
``unmix``
    Supervised separation against known reference spectra.
"""

__all__ = ['cd', 'common', 'ftir', 'melting', 'multivariate',
           'scattering', 'structure', 'unmix']
