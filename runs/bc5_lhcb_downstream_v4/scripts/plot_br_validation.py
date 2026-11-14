#!/usr/bin/env python
"""Validation figure: our widths versus the reference table."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.analysis.plots.br_validation import main

if __name__ == "__main__":
    main()
