"""End-to-end regression: the committed scan must reproduce the committed contour."""
import shutil
import sys

from bc5 import paths
from bc5.analysis import curves


def test_curves_reproduce_committed_files(tmp_path, monkeypatch):
    shutil.copy(paths.OUT_DATA / "signal_scan.json", tmp_path / "signal_scan.json")
    monkeypatch.setattr(curves, "DATA", str(tmp_path))
    monkeypatch.setattr(curves, "FIG", str(tmp_path))
    monkeypatch.setenv("BC5_SCAN_IN", "signal_scan.json")
    monkeypatch.setattr(sys, "argv", ["make_curves.py", "--s", "2.3"])
    curves.main()
    for name in ("exclusion_boundary.csv", "exclusion_s_curves.csv"):
        expected = (paths.OUT_DATA / name).read_text()
        got = (tmp_path / name).read_text()
        assert got == expected, name
