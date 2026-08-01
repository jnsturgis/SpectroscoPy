# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Plotting.

These check structure and the decisions that are easy to get silently wrong --
axis direction, labelling, palette order, legend presence -- rather than pixels.
Matplotlib is driven through the Agg backend so nothing tries to open a window.
"""

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from spectroscopy import Spectrum, SpectrumCollection, viz  # noqa: E402
from spectroscopy.lineshapes import gauss  # noqa: E402


@pytest.fixture(autouse=True)
def close_figures():
    yield
    plt.close('all')


def _spectrum(technique="ATR-FTIR", name="sample", scale=1.0):
    spectrum = Spectrum()
    spectrum.x = np.linspace(900.0, 1800.0, 451)
    spectrum.y = scale * (gauss(spectrum.x, 1650.0, 50.0, 1.0) + 0.1)
    spectrum.set_type(technique)
    spectrum.set_sample(name)
    spectrum.name = name
    return spectrum


@pytest.fixture
def collection():
    return SpectrumCollection([
        _spectrum(name="a", scale=1.0), _spectrum(name="a", scale=1.1),
        _spectrum(name="b", scale=2.0), _spectrum(name="b", scale=2.2),
    ])


# ---------------------------------------------------------------------------
# axis conventions
# ---------------------------------------------------------------------------

def test_ftir_x_axis_is_reversed():
    """Every FTIR figure cell types set_xlim((1800, 900)) by hand."""
    _, ax = plt.subplots()
    viz.plot(_spectrum("ATR-FTIR"), ax)
    assert ax.xaxis_inverted()


def test_uvvis_x_axis_is_not_reversed():
    _, ax = plt.subplots()
    viz.plot(_spectrum("UV-Vis"), ax)
    assert not ax.xaxis_inverted()


def test_axis_labels_come_from_the_spectrum():
    _, ax = plt.subplots()
    spectrum = _spectrum("ATR-FTIR")
    viz.plot(spectrum, ax)
    assert ax.get_xlabel() == spectrum.x_label
    assert "cm$^{-1}$" in ax.get_xlabel()
    assert ax.get_ylabel() == spectrum.y_label


def test_plotting_twice_does_not_re_invert():
    """Two spectra on one axes must not flip it back."""
    _, ax = plt.subplots()
    viz.plot(_spectrum(), ax)
    viz.plot(_spectrum(), ax)
    assert ax.xaxis_inverted()


def test_labels_can_be_suppressed():
    _, ax = plt.subplots()
    viz.plot(_spectrum(), ax, apply_labels=False)
    assert ax.get_xlabel() == ""
    assert not ax.xaxis_inverted()


# ---------------------------------------------------------------------------
# palette
# ---------------------------------------------------------------------------

def test_palette_is_assigned_in_fixed_order():
    """Colour follows the entity's position, never a shuffle."""
    assert [viz.series_style(i)['color'] for i in range(3)] == list(viz.PALETTE[:3])


def test_a_seventh_series_differs_in_style_not_just_colour():
    """
    The colours wrap at six, so the line style advances -- a 7th trace is not
    visually identical to the 1st.
    """
    first, seventh = viz.series_style(0), viz.series_style(len(viz.PALETTE))
    assert first['color'] == seventh['color']
    assert first['linestyle'] != seventh['linestyle']


def test_palette_excludes_the_unusable_okabe_ito_members():
    """Yellow is too light against white; black has no chroma."""
    assert '#F0E442' not in viz.PALETTE
    assert '#000000' not in viz.PALETTE
    assert len(viz.PALETTE) == len(set(viz.PALETTE))


# ---------------------------------------------------------------------------
# collections
# ---------------------------------------------------------------------------

def test_overlay_draws_every_spectrum_with_a_legend(collection):
    _, ax = plt.subplots()
    viz.plot_collection(collection, ax)
    assert len(ax.lines) == len(collection)
    assert ax.get_legend() is not None


def test_a_single_series_gets_no_legend_box():
    """The title or axis already names it; a one-entry key is noise."""
    _, ax = plt.subplots()
    viz.plot_collection(SpectrumCollection([_spectrum()]), ax)
    assert ax.get_legend() is None


def test_stack_offsets_traces_and_labels_them_directly(collection):
    _, ax = plt.subplots()
    viz.stack(collection, ax)

    assert len(ax.lines) == len(collection)
    tops = [line.get_ydata().max() for line in ax.lines]
    assert tops == sorted(tops), "each trace should sit above the previous one"
    assert len(ax.texts) == len(collection)
    assert ax.get_yticks().size == 0


def test_stack_accepts_explicit_offsets(collection):
    _, ax = plt.subplots()
    viz.stack(collection, ax, offsets=[0, 5, 10, 15])
    assert ax.lines[1].get_ydata().min() > ax.lines[0].get_ydata().min()


def test_grid_makes_one_panel_per_group(collection):
    figure, axes = viz.grid(collection, key='sample', ncols=2)
    visible = [panel for panel in axes.ravel() if panel.get_visible()]
    assert len(visible) == 2                     # samples 'a' and 'b'
    assert {panel.get_title(loc='left') for panel in visible} == {'a', 'b'}
    assert figure is not None


def test_grid_hides_unused_panels():
    collection = SpectrumCollection([_spectrum(name=n) for n in "abc"])
    _, axes = viz.grid(collection, ncols=2)
    assert sum(not panel.get_visible() for panel in axes.ravel()) == 1


# ---------------------------------------------------------------------------
# annotation
# ---------------------------------------------------------------------------

def test_annotate_peaks_marks_and_labels():
    spectrum = _spectrum()
    peaks = spectrum.find_peaks(prominence=1e-4)
    _, ax = plt.subplots()
    viz.plot(spectrum, ax)
    before = len(ax.lines)
    viz.annotate_peaks(ax, peaks)
    assert len(ax.lines) == before + 1
    assert len(ax.texts) == len(peaks)


def test_peaks_can_be_marked_without_labels():
    spectrum = _spectrum()
    peaks = spectrum.find_peaks(prominence=1e-4)
    _, ax = plt.subplots()
    viz.annotate_peaks(ax, peaks, label_peaks=False)
    assert len(ax.texts) == 0


def test_peak_table_annotate_delegates_to_viz():
    spectrum = _spectrum()
    peaks = spectrum.find_peaks(prominence=1e-4)
    _, ax = plt.subplots()
    peaks.annotate(ax)
    assert len(ax.texts) == len(peaks)


def test_annotate_bands_draws_a_line_and_a_label_each():
    _, ax = plt.subplots()
    viz.plot(_spectrum(), ax)
    bands = {1650: "Amide I", 1550: "Amide II"}
    viz.annotate_bands(ax, bands)
    assert len(ax.texts) == 2
    positions = sorted(line.get_xdata()[0] for line in ax.lines[1:])
    assert positions == [1550, 1650]


def test_plot_baseline_shows_all_three_traces():
    spectrum = _spectrum()
    baseline = spectrum.baseline('rubberband')
    _, ax = plt.subplots()
    viz.plot_baseline(spectrum, baseline, ax)
    assert len(ax.lines) == 3
    assert ax.get_legend() is not None


# ---------------------------------------------------------------------------
# delegation -- core no longer knows about matplotlib
# ---------------------------------------------------------------------------

def test_spectrum_plot_still_takes_an_axes_positionally():
    """The notebooks all call spectrum.plot(ax); that must keep working."""
    _, ax = plt.subplots()
    _spectrum().plot(ax)
    assert len(ax.lines) == 1


def test_spectrum_plot_accepts_a_format_string():
    """plot_spectra(ax, spectra, ["g-", "r-"]) is the Membrane_Analysis idiom."""
    _, ax = plt.subplots()
    _spectrum().plot(ax, "g-")
    assert ax.lines[0].get_color() == 'g'


def test_collection_plot_delegates(collection):
    _, ax = plt.subplots()
    collection.plot(ax)
    assert len(ax.lines) == len(collection)


def test_core_modules_do_not_import_matplotlib_at_module_scope():
    """
    Review crossing C2: the data model should not carry the plotting library.
    Checked statically -- importing cannot tell us, since viz imports it.
    """
    import ast
    import pathlib

    import spectroscopy.spectra

    offenders = []
    root = pathlib.Path(spectroscopy.spectra.__file__).parent
    for name in ("spectra.py", "collection.py", "peaks.py", "history.py",
                 "units.py", "processing/common.py"):
        tree = ast.parse((root / name).read_text(), filename=name)
        for node in tree.body:
            modules = []
            if isinstance(node, ast.Import):
                modules = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                modules = [node.module or ""]
            offenders += [f"{name}: {m}" for m in modules
                          if m.startswith("matplotlib")]
    assert not offenders, offenders
