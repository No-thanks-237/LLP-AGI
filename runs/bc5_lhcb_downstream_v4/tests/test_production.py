"""Production rates against the Boiarska tables."""
import pytest

from bc5.production import rates as P

PAPER_B = {'K': 0.430, 'pi': 0.013, 'K*892': 0.470, 'K*1410': 0.210, 'K*1680': 0.130,
           'K0*700': 0.760, 'K0*1430': 0.081, 'K1_1270': 0.910, 'K1_1400': 0.016,
           'K2*1430': 0.300}


def test_exclusive_branching_ratios():
    br = P.br_B_to_S(1e-6)
    for name, paper in PAPER_B.items():
        assert br[name] / paper == pytest.approx(1.0, abs=0.08), (name, br[name], paper)
    assert sum(br.values()) == pytest.approx(3.356, abs=0.25)


def test_b0_bplus_width_ratio():
    b_plus = sum(P.br_B_to_S(1e-6, 'B+').values())
    b_zero = sum(P.br_B_to_S(1e-6, 'B0').values())
    assert b_zero / b_plus == pytest.approx(0.9275, abs=0.002)


def test_kaon_rates():
    assert P.br_K_to_S(1e-6, 'K+') == pytest.approx(1.7e-3, rel=0.08)
    assert P.br_K_to_S(1e-6, 'K0L') == pytest.approx(7.0e-3, rel=0.08)
    assert P.br_K_to_S(1e-6, 'K0S') == 0.0


def test_quartic_channels_close():
    for ms in (0.1, 0.5, 1.0, 2.0):
        assert all(w >= 0 for _, _, _, w, _ in P.quartic_channels(ms))
    assert P.quartic_channels(10.0) == []


def test_higgs_portal_alpha():
    # BR(h -> SS) = 0.01 must be reproduced by the alpha returned for m_S -> 0
    ms = 1.0
    alpha = P.alpha_from_bh(ms)
    beta = (1.0 - 4.0 * ms ** 2 / P.M_H ** 2) ** 0.5
    gamma_ss = alpha ** 2 * beta / (32.0 * 3.141592653589793 * P.M_H)
    br = gamma_ss / (P.GAM_H_SM + gamma_ss)
    assert br == pytest.approx(P.B_H, rel=1e-6)
