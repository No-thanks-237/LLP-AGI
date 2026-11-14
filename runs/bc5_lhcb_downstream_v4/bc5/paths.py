"""Filesystem layout of one run.

Everything resolves through this module, so the directory can be copied anywhere and every
script keeps working.  Data locations can be overridden with environment variables:
BC5_SAMPLES, BC5_HADRON_LIB.
"""
from __future__ import annotations

import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

INPUTS = ROOT / "inputs"
WORK = ROOT / "work"
SAMPLES = Path(os.environ.get("BC5_SAMPLES", str(WORK / "magnus" / "samples")))
HADRON_LIB = Path(os.environ.get("BC5_HADRON_LIB", str(WORK / "magnus" / "hadron_lib")))
MODEL_DATA = WORK / "model"
REFERENCE_DATA = MODEL_DATA / "reference"
OUT_DATA = ROOT / "outputs" / "data"
OUT_FIG = ROOT / "outputs" / "figures"
LOGS = ROOT / "logs"
