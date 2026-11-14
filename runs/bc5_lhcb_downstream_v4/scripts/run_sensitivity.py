#!/usr/bin/env python
"""Run the BC5 downstream sensitivity end to end (scan + exclusion curves)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.analysis.sensitivity import main

if __name__ == "__main__":
    main()
