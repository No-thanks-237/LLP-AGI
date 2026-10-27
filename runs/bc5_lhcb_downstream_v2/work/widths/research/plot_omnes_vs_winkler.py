# -*- coding: utf-8 -*-
"""
Comparison figure for the research/ Omnes work: the single-channel Omnes pi pi rate driven by
the CFD phase shift, against the digitised Winkler curve that the production path uses.

    python research/plot_omnes_vs_winkler.py

Output: outputs/figures/widths_pipi_omnes_vs_winkler.{png,pdf}
Requires numpy, scipy and matplotlib (scipy is only needed for this optional figure).
"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
WIDTHS = os.path.dirname(HERE)
RUN = os.path.abspath(os.path.join(WIDTHS, '..', '..'))
sys.path.insert(0, HERE)
import phases as P
import omnes as O

M_PI, G_F = 0.13957039, 1.1663788e-5
S_PIPI = (2 * M_PI) ** 2


def delta_rad(t, s_to_2pi):
    return P.delta_guided(t, s_to_2pi)


def gamma_pipi_omnes(ms, s_to_2pi=25.0, p=0.73):
    s = ms * ms
    beta = math.sqrt(max(1.0 - S_PIPI / s, 0.0))
    Om = O.omnes_single(s, delta_rad, S_PIPI, s_to_2pi)
    A = Om * ((7.0 / 9.0) * M_PI ** 2 + (2.0 / 9.0) * (2 * M_PI ** 2 + p * s))
    return 3.0 * G_F / (16.0 * math.sqrt(2.0) * math.pi * ms) * beta * abs(A) ** 2


def main():
    tab = np.genfromtxt(os.path.join(WIDTHS, 'hadronic_table.csv'), delimiter=',', names=True)
    ms = np.exp(np.linspace(math.log(0.29), math.log(0.985), 40))
    g_om = np.array([gamma_pipi_omnes(m) for m in ms])
    g_w = np.interp(ms, tab['mass_GeV'], tab['gamma_pipi_GeV'])

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(ms, g_w, 'k-', lw=2, label='digitised Winkler Fig. 2 (production input)')
    ax.plot(ms, g_om, 'r--', lw=2, label=r'single-channel Omnes, CFD $\delta_0^0$ (research)')
    ax.set_yscale('log')
    ax.set_xlabel(r'$m_S$ [GeV]')
    ax.set_ylabel(r'$\Gamma(S\to\pi\pi)/\theta^2$ [GeV]')
    ax.set_title('Omnes dispersion vs the digitised reference (elastic region)')
    ax.grid(alpha=0.3, which='both')
    ax.legend()
    fig.tight_layout()
    out = os.path.join(RUN, 'outputs', 'figures', 'widths_pipi_omnes_vs_winkler')
    fig.savefig(out + '.png', dpi=150)
    fig.savefig(out + '.pdf')
    r = g_om / g_w
    print('wrote %s.png / .pdf' % out)
    print('ratio Omnes/Winkler: min %.3f  max %.3f ; mean |ratio-1| = %.3f'
          % (r.min(), r.max(), np.mean(np.abs(r - 1.0))))


if __name__ == '__main__':
    main()
