#!/usr/bin/env python
"""Verify inputs/source.sha256 against the files in inputs/."""
import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "inputs" / "source.sha256"


def main():
    problems = []
    checked = 0
    for line in MANIFEST.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, rel = line.split(maxsplit=1)
        path = ROOT / rel
        checked += 1
        if not path.is_file():
            problems.append((rel, "missing"))
            continue
        got = hashlib.sha256(path.read_bytes()).hexdigest()
        if got != digest:
            problems.append((rel, "mismatch"))
    for rel, why in problems:
        print("%s: %s" % (why, rel))
    print("checked %d file(s), %d problem(s)" % (checked, len(problems)))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
