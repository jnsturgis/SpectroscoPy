# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The format registry, and the generic delimited-table reader.

Roadmap section 3: "a registry of format readers, not a chain of if/elif on
file extension". The four hand-kept tables it replaces had already drifted
apart -- that was defect D5.
"""

import numpy as np
import pytest

from spectroscopy import io
from spectroscopy.io import registry, table
from spectroscopy.spectra import Spectrum

# ---------------------------------------------------------------------------
# registry
# ---------------------------------------------------------------------------

def test_the_built_in_formats_are_registered():
    assert set(io.known_types()) >= {'csv', 'dpt', 'jcamp', 'spy', 'tsv', 'table'}


@pytest.mark.parametrize("name,expected", [
    ("a.dpt", "dpt"), ("a.DPT", "dpt"),
    ("a.jdx", "jcamp"), ("a.dx", "jcamp"), ("a.jcamp", "jcamp"),
    ("a.csv", "csv"), ("a.CSV", "csv"),
    ("a.tsv", "tsv"), ("a.txt", "tsv"),
    ("a.spy", "spy"),
    ("a.qqq", "unknown"),
])
def test_extension_inference_comes_from_the_registry(name, expected):
    assert io.infer_file_type(name) == expected


def test_spectrum_file_type_inference_agrees_with_the_registry():
    """
    Spectrum used to keep its own FILE_EXTS table. It now derives from the
    registry, so the two cannot drift apart again.
    """
    from spectroscopy.spectra import _infer_file_type
    for name in ("x.dpt", "x.jdx", "x.csv", "x.spy", "x.nope"):
        assert _infer_file_type(name) == io.infer_file_type(name)


def test_registering_a_new_format_teaches_spectrum_about_it():
    """The point of a registry: one step, and everything downstream follows."""
    from spectroscopy.spectra import _infer_file_type

    @registry.register_reader('demo_fmt', extensions=['.demo'],
                              description='test only')
    def _read(handle, spectrum, **kwargs):        # pragma: no cover - trivial
        values = [float(line) for line in handle if line.strip()]
        spectrum.x = np.arange(len(values), dtype=float)
        spectrum.y = np.array(values)

    try:
        assert 'demo_fmt' in io.known_types()
        assert io.infer_file_type("thing.demo") == 'demo_fmt'
        assert _infer_file_type("thing.demo") == 'demo_fmt'
    finally:
        registry.REGISTRY.pop('demo_fmt', None)


def test_a_read_only_format_cannot_be_written(tmp_path):
    spectrum = Spectrum()
    spectrum.x, spectrum.y = np.array([1.0]), np.array([2.0])
    with pytest.raises(ValueError, match="Cannot write"):
        registry.write_spectrum(spectrum, str(tmp_path / "x.out"), 'table')


def test_unknown_type_is_a_type_error(tmp_path):
    with pytest.raises(TypeError, match="Unknown filetype"):
        registry.read_spectra(str(tmp_path / "x.nope"))


def test_describe_formats_lists_capabilities():
    text = io.describe_formats()
    assert 'dpt' in text and 'jcamp' in text
    assert 'r-' in text          # 'table' is read-only


def test_read_spectrum_refuses_a_multi_spectrum_file(tmp_path):
    """Silently returning the first of many is exactly the quiet-wrong to avoid."""
    path = tmp_path / "wide.csv"
    path.write_text("x,a,b\n1,10,20\n2,11,21\n")
    with pytest.raises(ValueError, match="read_spectra"):
        registry.read_spectrum(str(path), 'table')

    assert len(registry.read_spectra(str(path), 'table')) == 2


# ---------------------------------------------------------------------------
# encoding detection
# ---------------------------------------------------------------------------

def test_utf16_with_a_bom_is_detected(tmp_path):
    """The AKTA chromatography exports are UTF-16-LE; UTF-8 gives nonsense."""
    path = tmp_path / "akta.csv"
    path.write_bytes("x\ty\n1.0\t2.0\n".encode("utf-16-le"))
    with open(path, "wb") as handle:
        handle.write(b"\xff\xfe" + "x\ty\n1.0\t2.0\n".encode("utf-16-le"))
    assert registry.detect_encoding(str(path)) == 'utf-16-le'
    assert len(registry.read_spectra(str(path), 'table')) == 1


def test_latin1_is_the_backstop(tmp_path):
    path = tmp_path / "old.csv"
    path.write_bytes(b"# caf\x97\n1.0,2.0\n3.0,4.0\n")
    assert registry.detect_encoding(str(path)) == 'latin-1'
    assert len(registry.read_spectra(str(path), 'csv')[0]) == 2


# ---------------------------------------------------------------------------
# the generic table reader
# ---------------------------------------------------------------------------

def test_shared_x_column_names_series_from_the_y_header(tmp_path):
    """The GFP layout: one wavelength column, many named sample columns."""
    path = tmp_path / "wide.csv"
    path.write_text("Wavelength (nm),GFP 0.5uM,GFP 100nM\n"
                    "500,436,63\n501,400,60\n")
    spectra = registry.read_spectra(str(path), 'table', x_col=0)
    assert [s.name for s in spectra] == ["GFP 0.5uM", "GFP 100nM"]
    assert np.allclose(spectra[0].x, [500, 501])
    assert np.allclose(spectra[0].y, [436, 400])


def test_paired_columns_name_series_from_the_x_header(tmp_path):
    """Chloe's fluorimeter layout: (x, y) pairs with sparse series names."""
    path = tmp_path / "paired.csv"
    path.write_text("Sample A,,Sample B,\n"
                    "Wavelength (nm),Intensity (a.u.),Wavelength (nm),Intensity (a.u.)\n"
                    "465,13,465,678\n466,14,466,679\n")
    spectra = registry.read_spectra(str(path), 'table', paired=True)
    assert [s.name for s in spectra] == ["Sample A", "Sample B"]
    assert spectra[0].x_label == "Wavelength (nm)"
    assert spectra[0].y_label == "Intensity (a.u.)"
    assert np.allclose(spectra[1].y, [678, 679])


def test_the_naming_header_row_is_chosen_by_distinctness(tmp_path):
    """
    Two real files disagree about which row names the series, so the reader
    picks the row that actually distinguishes them rather than a fixed one.
    Here row 1 is the discriminating one, the reverse of the case above.
    """
    path = tmp_path / "akta.csv"
    path.write_text("Chrom.1\t\tChrom.1\t\n"
                    "UV\tml\tConductivity\tml\n"
                    "0.0\t1.0\t0.0\t5.0\n0.1\t1.1\t0.1\t5.1\n")
    spectra = registry.read_spectra(str(path), 'table', paired=True)
    assert [s.name for s in spectra] == ["UV", "Conductivity"]


def test_a_byte_order_mark_does_not_leak_into_a_name(tmp_path):
    path = tmp_path / "bom.csv"
    path.write_bytes(b"\xef\xbb\xbf" + b"Trace A,\nx,y\n1,2\n3,4\n")
    spectra = registry.read_spectra(str(path), 'table', paired=True)
    assert spectra[0].name == "Trace A"


def test_delimiter_and_header_count_are_sniffed(tmp_path):
    path = tmp_path / "tabs.txt"
    path.write_text("a\tb\n1\t2\n3\t4\n")
    spectra = registry.read_spectra(str(path), 'table', x_col=0)
    assert len(spectra) == 1
    assert np.allclose(spectra[0].y, [2, 4])


def test_ragged_and_blank_cells_are_dropped_per_series(tmp_path):
    """Series in one file often have different lengths."""
    path = tmp_path / "ragged.csv"
    path.write_text("A,,B,\nx,y,x,y\n1,10,1,20\n2,11,,\n3,12,,\n")
    spectra = registry.read_spectra(str(path), 'table', paired=True)
    assert len(spectra[0]) == 3
    assert len(spectra[1]) == 1


def test_explicit_names_win(tmp_path):
    path = tmp_path / "wide.csv"
    path.write_text("x,a,b\n1,10,20\n2,11,21\n")
    spectra = registry.read_spectra(str(path), 'table', x_col=0,
                                    names=["first", "second"])
    assert [s.name for s in spectra] == ["first", "second"]


def test_an_empty_table_is_refused(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text("# only a comment\n")
    with pytest.raises(ValueError, match="no data"):
        registry.read_spectra(str(path), 'table')


def test_sniff_helpers():
    assert table.sniff_delimiter(["1\t2\t3"]) == '\t'
    assert table.sniff_delimiter(["1,2,3"]) == ','
    assert table.sniff_header_rows(["a,b", "1,2"], ',') == 1
    assert table.sniff_header_rows(["1,2", "3,4"], ',') == 0
