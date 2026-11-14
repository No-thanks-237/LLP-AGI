"""S decay: phase-space sampling, channel composition, hadronisation library."""
from .kinematics import (PI, M_E, M_MU, M_TAU, M_PI, M_PI0, M_K, M_K0L, M_ETA,
                         boost, two_body, nbody, jet_fragment)
from .hadronlib import HadronLib, NODE_MTOL, load_hadron_lib
from .channels import CHANNELS, CTAU_KS, CTAU_KL, s_decay_products, channel_branching

__all__ = [
    "PI", "M_E", "M_MU", "M_TAU", "M_PI", "M_PI0", "M_K", "M_K0L", "M_ETA",
    "boost", "two_body", "nbody", "jet_fragment",
    "HadronLib", "NODE_MTOL", "load_hadron_lib",
    "CHANNELS", "CTAU_KS", "CTAU_KL", "s_decay_products", "channel_branching",
]
