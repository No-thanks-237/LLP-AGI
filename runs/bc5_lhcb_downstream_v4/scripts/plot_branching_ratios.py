#!/usr/bin/env python
"""Branching ratios of the BC5 scalar versus its mass."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.analysis.plots.branching import main

if __name__ == "__main__":
    main()
