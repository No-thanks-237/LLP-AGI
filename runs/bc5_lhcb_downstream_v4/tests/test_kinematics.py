"""Phase-space sampler invariants (fast, no experimental data needed)."""
import math
import random

import numpy as np

from bc5.decay import nbody, two_body, boost


def _p4_sum(parts):
    e = sum(p[1] for p in parts)
    px = sum(p[2] for p in parts)
    py = sum(p[3] for p in parts)
    pz = sum(p[4] for p in parts)
    return e, px, py, pz


def test_two_body_threshold_and_momentum():
    rng = random.Random(1)
    assert two_body(0.2, 0.15, 0.15, rng) is None
    parts = two_body(1.0, 0.1, 0.2, rng)
    p = math.sqrt((1.0 ** 2 - (0.1 + 0.2) ** 2) * (1.0 ** 2 - (0.1 - 0.2) ** 2)) / 2.0
    assert abs(math.hypot(parts[0][2], parts[0][3], parts[0][4]) - p) < 1e-12


def test_four_momentum_conservation():
    rng = random.Random(2)
    for n, masses in ((2, [0.1, 0.2]), (3, [0.14, 0.14, 0.14]), (4, [0.14] * 4)):
        parts = nbody(1.5, masses, rng)
        assert parts is not None and len(parts) == n
        e, px, py, pz = _p4_sum(parts)
        assert abs(e - 1.5) < 1e-9
        assert abs(px) < 1e-9 and abs(py) < 1e-9 and abs(pz) < 1e-9


def test_massless_four_body_energy_sharing():
    """Uniform 4-body phase space gives <E_i/M> = 1/4 for every daughter."""
    rng = random.Random(3)
    acc = np.zeros(4)
    n = 4000
    for _ in range(n):
        parts = nbody(1.0, [0.0] * 4, rng)
        acc += np.array([p[1] for p in parts])
    acc /= n
    assert np.allclose(acc, 0.25, atol=0.01), acc


def test_massive_three_body_s23():
    """<s23>/M^2 for 3 pi at M = 1 GeV is 0.3528 (exact uniform phase space)."""
    rng = random.Random(4)
    n = 4000
    m = 0.13957039
    vals = []
    for _ in range(n):
        parts = nbody(1.0, [m, m, m], rng)
        e2, e3 = parts[1][1], parts[2][1]
        p2, p3 = parts[1][2:], parts[2][2:]
        s23 = (e2 + e3) ** 2 - sum((a + b) ** 2 for a, b in zip(p2, p3))
        vals.append(s23)
    assert abs(float(np.mean(vals)) - 0.352813) < 0.005, float(np.mean(vals))


def test_boost_round_trip():
    rng = random.Random(5)
    parts = two_body(2.0, 0.3, 0.4, rng)
    lab = boost(parts, 3.0, (0.0, 0.0, 2.0))
    back = boost(lab, 3.0, (0.0, 0.0, -2.0))
    for a, b in zip(parts, back):
        assert abs(a[1] - b[1]) < 1e-12 and abs(a[4] - b[4]) < 1e-12
