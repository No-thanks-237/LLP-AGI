#!/usr/bin/env python
"""Validate the width model against the digitised reference curves."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from bc5.model.validation import main

if __name__ == "__main__":
    main()
