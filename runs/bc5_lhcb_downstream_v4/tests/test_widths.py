"""Model-level checks: branching ratios, thresholds, lifetime."""
import math

import numpy as np

from bc5.model import widths as sw


def test_branching_ratios_sum_to_one():
    for m in (0.021, 0.1, 0.2802, 0.5, 1.0, 2.0, 2.1, 5.0, 30.0, 62.49):
        w = sw.widths(m)
        tot = sum(w.values())
        assert tot > 0
        assert abs(sum(w.values()) / tot - 1.0) < 1e-12
        assert all(v >= 0 for v in w.values())


def test_thresholds():
    # pi pi and K K are zero below their two-meson thresholds
    assert sw.widths(0.27)['pipi'] == 0.0
    assert sw.widths(0.98)['KK'] == 0.0
    assert sw.widths(0.29)['pipi'] > 0.0
    # c cbar opens at 2 m_D, b bbar at 2 m_B
    assert sw.widths(3.7)['cc'] == 0.0
    assert sw.widths(3.8)['cc'] > 0.0
    assert sw.widths(10.5)['bb'] == 0.0
    assert sw.widths(10.7)['bb'] > 0.0


def test_ctau_scales_with_theta2():
    m = 1.0
    a = sw.ctau_m(m, 1e-8)
    b = sw.ctau_m(m, 1e-9)
    assert abs(b / a - 10.0) < 1e-12


def test_leptonic_width_closed_form():
    m, ml = 1.0, sw.M_MU
    beta = math.sqrt(1.0 - 4.0 * ml ** 2 / m ** 2)
    exact = sw.G_F * m * ml ** 2 * beta ** 3 / (4.0 * math.sqrt(2.0) * math.pi)
    assert abs(sw.gamma_ll(m, ml) / exact - 1.0) < 1e-12


def test_hadronic_table_loaded():
    assert sw.HADRONIC_TABLE is not None
    m = np.asarray(sw.HADRONIC_TABLE[0])
    assert m[0] < 0.3 and m[-1] > 1.9
