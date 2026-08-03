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


def test_a_format_registered_after_import_is_usable(tmp_path):
    """
    The registry must be consulted live, not snapshotted.

    ``spectra.KNOWNFILETYPES`` is filled in when the module is imported. A
    format registered after that -- which is the normal case, since
    ``register_reader`` is a decorator anyone can apply in their own code --
    was inferred correctly from the extension and then rejected by the
    Spectrum constructor as an unknown type. The documented promise that one
    decorator is enough was not true.
    """
    import numpy as np

    from spectroscopy.spectra import Spectrum

    @registry.register_reader('latecomer', extensions=['.late'],
                              description='registered after import')
    def _read(handle, spectrum, **kwargs):
        rows = [line.split(';') for line in handle.read().splitlines()
                if line.strip()]
        spectrum.x = np.array([float(row[0]) for row in rows])
        spectrum.y = np.array([float(row[1]) for row in rows])

    path = tmp_path / "sample.late"
    path.write_text("400;0.10\n500;0.25\n600;0.40\n")

    assert registry.infer_file_type(path) == 'latecomer'
    spectrum = Spectrum.read(path)          # used to raise TypeError
    assert len(spectrum.x) == 3
    assert spectrum.y[-1] == 0.40


def test_unknown_filetype_error_lists_what_is_known(tmp_path):
    from spectroscopy.spectra import Spectrum

    path = tmp_path / "sample.nonsense"
    path.write_text("1 2\n")
    with pytest.raises(TypeError) as caught:
        Spectrum.read(path)
    assert 'jcamp' in str(caught.value)


# ---------------------------------------------------------------------------
# Locale: comma decimals and the separators that come with them
# ---------------------------------------------------------------------------

LOCALE_CASES = {
    'fr_semicolon': "Longueur d'onde;Absorbance\n400,5;0,1234\n401,0;0,2345\n402,5;0,3456\n",
    'us_comma':     "Wavelength,Absorbance\n400.5,0.1234\n401.0,0.2345\n402.5,0.3456\n",
    'de_tab':       "Wellenlaenge\tExtinktion\n400,5\t0,1234\n401,0\t0,2345\n402,5\t0,3456\n",
    'fr_noheader':  "400,5;0,1234\n401,0;0,2345\n402,5;0,3456\n",
    'scientific':   "x,y\n400.5,1.2e-3\n401.0,2.3e-3\n402.5,3.4e-3\n",
}


@pytest.mark.parametrize('name', sorted(LOCALE_CASES))
def test_a_csv_is_read_whatever_the_locale_wrote(name, tmp_path):
    """
    A '.csv' from a European colleague is routinely neither comma separated
    nor dot decimal: the locale takes the comma for the decimal, so Excel
    exports ';' instead. Nothing in the file says so.
    """
    from spectroscopy.spectra import Spectrum

    path = tmp_path / f"{name}.csv"
    path.write_text(LOCALE_CASES[name])
    spectrum = Spectrum.read(path)

    assert np.allclose(spectrum.x, [400.5, 401.0, 402.5])
    expected = ([1.2e-3, 2.3e-3, 3.4e-3] if name == 'scientific'
                else [0.1234, 0.2345, 0.3456])
    assert np.allclose(spectrum.y, expected)


def test_grouping_separators_are_removed(tmp_path):
    """1.400,5 and 1 400,5 are both fourteen hundred and a half."""
    from spectroscopy.spectra import Spectrum

    for text in ("x;y\n1.400,5;0,1\n1.401,0;0,2\n1.402,5;0,3\n",
                 "x;y\n1 400,5;0,1\n1 401,0;0,2\n1 402,5;0,3\n"):
        path = tmp_path / "grouped.csv"
        path.write_text(text)
        assert np.allclose(Spectrum.read(path).x, [1400.5, 1401.0, 1402.5])


def test_an_explicit_delimiter_or_decimal_wins(tmp_path):
    from spectroscopy.spectra import Spectrum

    path = tmp_path / "fr.csv"
    path.write_text("400,5;0,1234\n401,0;0,2345\n")
    spectrum = Spectrum.read(path)
    assert np.allclose(spectrum.x, [400.5, 401.0])

    forced = registry.read_spectrum(path, 'csv', delimiter=';', decimal=',')
    assert np.allclose(forced.x, [400.5, 401.0])


def test_the_decimal_is_sniffed_even_when_the_separator_is_known(tmp_path):
    """A .tsv pins the tab but says nothing about the locale."""
    from spectroscopy.spectra import Spectrum

    path = tmp_path / "german.tsv"
    path.write_text("x\ty\n400,5\t0,1234\n401,0\t0,2345\n")
    spectrum = Spectrum.read(path)
    assert np.allclose(spectrum.x, [400.5, 401.0])
    assert np.allclose(spectrum.y, [0.1234, 0.2345])


def test_a_wide_table_can_be_comma_decimal(tmp_path):
    path = tmp_path / "series.csv"
    path.write_text("lambda;A;B\n400,0;0,10;0,20\n401,0;0,11;0,22\n402,0;0,12;0,24\n")
    collection = registry.read_spectra(path, 'table', x_col=0)
    assert len(collection) == 2
    assert np.allclose(collection[0].x, [400.0, 401.0, 402.0])
    assert np.allclose(collection[1].y, [0.20, 0.22, 0.24])


def test_number_parsing_helpers():
    assert table.parse_number('0,1234', ',') == pytest.approx(0.1234)
    assert table.parse_number('1.400,5', ',') == pytest.approx(1400.5)
    assert table.parse_number('400.5', '.') == pytest.approx(400.5)
    assert table.sniff_format(["400,5;0,1234", "401,0;0,2345"]) == (';', ',')
    assert table.sniff_format(["400.5,0.1234", "401.0,0.2345"]) == (',', '.')
    assert table.sniff_format(["400,5\t0,1234", "401,0\t0,2345"]) == ('\t', ',')


@pytest.mark.parametrize('file_type,kwargs,expected_sep', [
    ('csv', {}, ','),
    ('csv', {'decimal': ','}, ';'),
    ('csv', {'decimal': ',', 'delimiter': '\t'}, '\t'),
    ('tsv', {'decimal': ','}, '\t'),
])
def test_a_comma_decimal_can_be_written_and_read_back(file_type, kwargs,
                                                      expected_sep, tmp_path):
    """
    Round trip is the test that matters: whatever we write, we must read.
    """
    from spectroscopy.spectra import Spectrum

    original = Spectrum(np.array([400.5, 401.0, 402.5]),
                        np.array([0.1234, 0.2345, 0.3456]),
                        x_unit='nm', y_unit='absorbance')
    path = tmp_path / f"out.{file_type}"
    original.save_as(str(path), file_type, **kwargs)

    text = path.read_text()
    assert expected_sep in text
    if kwargs.get('decimal') == ',':
        assert '400,500' in text
        assert '400.500' not in text

    restored = Spectrum.read(path)
    assert np.allclose(restored.x, original.x)
    assert np.allclose(restored.y, original.y)


def test_writing_a_comma_for_both_separator_and_decimal_is_refused(tmp_path):
    """
    '400,5,0,1234' is four fields or two and nothing can tell which -- so this
    would produce a file no reader could get right, including ours.
    """
    from spectroscopy.spectra import Spectrum

    spectrum = Spectrum(np.array([1.0, 2.0]), np.array([3.0, 4.0]))
    with pytest.raises(ValueError, match="both the field separator"):
        spectrum.save_as(str(tmp_path / "bad.csv"), 'csv',
                         decimal=',', delimiter=',')


def test_an_unknown_decimal_is_refused(tmp_path):
    from spectroscopy.spectra import Spectrum

    spectrum = Spectrum(np.array([1.0, 2.0]), np.array([3.0, 4.0]))
    with pytest.raises(ValueError, match=r"decimal must be"):
        spectrum.save_as(str(tmp_path / "bad.csv"), 'csv', decimal='x')


def test_save_as_passes_options_through_to_the_writer(tmp_path):
    """Otherwise decimal= is only reachable from the low-level function."""
    from spectroscopy.spectra import Spectrum

    path = tmp_path / "out.csv"
    Spectrum(np.array([1.0]), np.array([2.0])).save_as(str(path), 'csv',
                                                       decimal=',')
    assert '1,000;2,00000' in path.read_text()
