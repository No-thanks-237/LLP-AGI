#!/usr/bin/env python
"""Audit every width channel against the local reference data."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.model.audit import main

if __name__ == "__main__":
    main()
