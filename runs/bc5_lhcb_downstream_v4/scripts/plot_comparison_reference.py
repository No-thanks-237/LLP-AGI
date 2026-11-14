#!/usr/bin/env python
"""Overlay our contour on the digitised BC5 reference curve."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.analysis.plots.reference import main

if __name__ == "__main__":
    main()
