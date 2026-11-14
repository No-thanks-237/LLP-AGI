#!/usr/bin/env python
"""Parent-particle production parameter-space figures."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.analysis.plots.production import main

if __name__ == "__main__":
    main()
