# -*- coding: utf-8 -*-
"""Validate the meson -> S branching ratios against Boiarska et al. tab:BR."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bc5_production as P

EXPECT_B = {'K': 0.430, 'pi': 0.013, 'K*892': 0.470, 'K*1410': 0.210, 'K*1680': 0.130,
            'K0*700': 0.760, 'K0*1430': 0.081, 'K1_1270': 0.910, 'K1_1400': 0.016,
            'K2*1430': 0.300}
EXPECT_TOTAL = 3.3
EXPECT_K = {'K+': 1.7e-3, 'K0L': 7.0e-3}


def main():
    br = P.br_B_to_S(1e-6)
    rep = {'paper': 'Boiarska et al. arXiv:1904.10447 tab:BR (m_S -> 0, theta^2-normalised)',
           'B_to_XsS': {}, 'total_exclusive': {'ours': sum(br.values()), 'paper': EXPECT_TOTAL},
           'K_to_piS': {}}
    for k, v in br.items():
        rep['B_to_XsS'][k] = {'ours': v, 'paper': EXPECT_B[k],
                              'ratio': v / EXPECT_B[k] if EXPECT_B[k] else None}
    for which, tag in (('K+', 'K+'), ('K0L', 'K0L')):
        v = P.br_K_to_S(1e-6, which)
        rep['K_to_piS'][which] = {'ours': v, 'paper': EXPECT_K[tag], 'ratio': v / EXPECT_K[tag]}
    rep['alpha_from_bh_1GeV'] = P.alpha_from_bh(1.0)
    rep['sigma_h_total_pb'] = P.SIGMA_H_TOTAL
    out = os.path.join(HERE, 'production_validation.json')
    with open(out, 'w') as f:
        json.dump(rep, f, indent=2)
    print(json.dumps(rep, indent=2))


if __name__ == '__main__':
    main()
