"""The decay-in-window integral must be exact for an analytic A(z)."""
import math

import numpy as np

from bc5.response import ZGRID, lifetime_probability, p_window


def _analytic(alpha, z0, zmin=1.0, zmax=2.5):
    lo = max(zmin, z0)
    if z0 >= zmax:
        return 0.0
    return math.exp(-alpha * (lo - z0)) - math.exp(-alpha * (zmax - z0))


def test_window_integral_matches_analytic():
    A = np.ones(len(ZGRID))
    for alpha in (1e-3, 0.1, 1.0, 10.0):
        for z0 in (0.0, 0.5, 1.2, 2.4):
            # p_window takes (ms, pz, z0, A, ctau); choose ms/ctau so that alpha = ms/ctau/pz
            pz = 1.0
            ctau = 1.0 / alpha
            got = float(p_window(1.0, np.array([pz]), np.array([z0]), A, ctau)[0])
            assert abs(got - _analytic(alpha, z0)) < 1e-12, (alpha, z0, got, _analytic(alpha, z0))


def test_scalar_helper_matches_vector():
    A = np.linspace(0.2, 0.9, len(ZGRID))
    for ms, pz, z0, ctau in ((0.5, 2.0, 0.3, 0.05), (5.0, 30.0, 1.1, 0.2)):
        vec = float(p_window(ms, np.array([pz]), np.array([z0]), A, ctau)[0])
        scal = lifetime_probability(ms, (1.0, 0.0, 0.0, pz), z0, A, ctau)
        assert abs(vec - scal) < 1e-15


def test_daughter_after_window_never_decays_inside():
    A = np.ones(len(ZGRID))
    # born downstream of the window -> probability exactly zero
    assert float(p_window(1.0, np.array([1.0]), np.array([2.6]), A, 0.1)[0]) == 0.0
