# -*- coding: utf-8 -*-
"""Write work/production/production_reference.csv: production BRs / Higgs sigmas on the
run mass grid for a few tanb values (reference table for review; the scan calls the
module directly).  Runs in ~1 s."""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import twohdm_production as P

TANB = [5.0, 10.0, 30.0, 100.0, 300.0, 1e3, 1e4, 1e5]
COLS = ['mH', 'tanb', 'br_B_foresee', 'br_B_paper_mb', 'br_B_exclusive',
        'br_K_to_piH_Kp', 'br_K_to_piH_K0L', 'br_K_to_Henu_Kp', 'br_pi_to_Henu',
        'BR_hHH', 'Gamma_h_other', 'sigma_total', 'br_B_validity_flag']


def main():
    out = os.path.join(HERE, 'production_reference.csv')
    with open(out, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(COLS)
        for mH in P.mass_grid():
            for tb in TANB:
                e = P.evaluate(mH, tb)
                w.writerow(['%.6g' % mH, '%.6g' % tb,
                            '%.8e' % e['br_B_to_XsH'],
                            '%.8e' % e['br_B_to_XsH_paper_mb'],
                            '%.8e' % e['br_B_to_XsH_exclusive'],
                            '%.8e' % e['br_K_to_piH_K+'],
                            '%.8e' % e['br_K_to_piH_K0L'],
                            '%.8e' % e['br_K_to_Henu_K+'],
                            '%.8e' % e['br_pi_to_Henu'],
                            '%.8e' % e['BR_hHH'],
                            '%.8e' % e['Gamma_h_other'],
                            '%.8f' % e['sigma_total'],
                            e['br_B_to_XsH_validity']])
    print('wrote', out)


if __name__ == '__main__':
    main()
