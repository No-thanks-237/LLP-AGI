"""The committed scan must name the exact inputs it was produced from."""
import hashlib
import json

from bc5 import paths


def _resolve(rel):
    if rel.startswith(("samples/", "hadron_lib/")):
        return paths.WORK / "magnus" / rel
    return paths.ROOT / rel


def test_fingerprints_match_files():
    scan = json.loads((paths.OUT_DATA / "signal_scan.json").read_text())
    inputs = scan["inputs"]
    assert len(inputs) == 16
    bad = []
    for rel, digest in inputs.items():
        path = _resolve(rel)
        if not path.is_file():
            bad.append((rel, "missing"))
        elif hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            bad.append((rel, "mismatch"))
    assert not bad, bad


def test_scan_shape():
    scan = json.loads((paths.OUT_DATA / "signal_scan.json").read_text())
    assert len(scan["masses"]) == 67
    assert len(scan["theta2"]) == 240
    assert scan["model_tag"].startswith("v4.0")
    assert min(scan["nsamp_by_mass"].values()) >= 1200
