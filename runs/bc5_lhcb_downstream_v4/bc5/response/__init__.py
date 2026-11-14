"""LHCb downstream detector response and lifetime integral."""
from . import kernel
from .lifetime import (Z_DEC_MIN, Z_DEC_MAX, ETA_S_MIN, ETA_S_MAX, RECO_EFF, ZGRID,
                        acceptance_profile, p_window, lifetime_probability)

__all__ = [
    "kernel",
    "Z_DEC_MIN", "Z_DEC_MAX", "ETA_S_MIN", "ETA_S_MAX", "RECO_EFF", "ZGRID",
    "acceptance_profile", "p_window", "lifetime_probability",
]
