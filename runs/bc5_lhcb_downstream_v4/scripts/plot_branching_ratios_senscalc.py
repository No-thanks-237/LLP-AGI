#!/usr/bin/env python
"""Branching ratios from the SensCalc reference tables."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.analysis.plots.senscalc import main

if __name__ == "__main__":
    main()
