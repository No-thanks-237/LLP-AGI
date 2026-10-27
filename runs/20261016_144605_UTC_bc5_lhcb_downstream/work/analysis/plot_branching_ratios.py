# -*- coding: utf-8 -*-
"""
Branching ratios of the BC5 scalar versus its mass — complete version.

Every channel is drawn over the whole scanned mass range 0.021 - 62.49 GeV.
  * solid  : the channel IS the production input in that mass region
  * dashed : the curve is shown but is NOT part of the model there
             - ss, gg below m_S = 2 GeV: the perturbative spectator model is only used above
               2 GeV; below that the hadronic width is dispersive (pi pi, K K, multi-meson).
               cc and bb are exactly zero below their hadronic thresholds 2 m_D = 3.74 GeV and
               2 m_B = 10.56 GeV, so they simply start there.
             - pi+ pi-, K+ K-, pi+pi-pi+pi- above m_S = 2 GeV: these channels do not exist
               there (the hadronic width is perturbative), so nothing is drawn.
             - s sbar, c cbar, b bbar, g g below m_S = 2 GeV: likewise, these are the perturbative
               description and are not drawn where the model uses the dispersive channels.
           Nothing is extrapolated: every curve covers exactly the mass range in which the model
           produces that channel.

Only exclusive channels are shown, so the curves do not sum to one:
  pi+ pi-      = 2/3 of the total pi pi width   (Winkler's Gamma_pipi is the isoscalar sum)
  K+ K-        = 1/2 of the total K K width     (the other half is K_L K_S)
  multi-meson  = the FULL toy channel C m_S^3 beta_2pi, which stands for 4pi + eta eta + rho rho
Not shown: pi0 pi0, K_L K_S, gamma gamma.

Each channel is drawn only over the mass range in which the model produces it: the dispersive
pi pi / K K / multi-meson channels below the 2 GeV matching scale, and s sbar, c cbar, b bbar,
g g only above it.  Nothing is extrapolated.

Outputs: outputs/figures/bc5_branching_ratios.{png,pdf}
         outputs/data/bc5_branching_ratios.csv
"""
from __future__ import annotations
import csv, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(RUN, 'work', 'widths'))
sys.path.insert(0, HERE)
import scalar_widths as sw

sw.load_hadronic_table(os.path.join(RUN, 'work', 'widths', 'hadronic_table.csv'))

DATA = os.path.join(RUN, 'outputs', 'data')
FIG = os.path.join(RUN, 'outputs', 'figures')
os.makedirs(FIG, exist_ok=True)

PERT = 2.0        # the Watson/Winkler matching scale

# key, label, colour, exclusive weight, "mode":
#   'lep'   -> open above its own threshold, solid wherever non-zero
#   'dis'   -> dispersive exclusive channel, solid below PERT, 1/m^2 dashed above
#   'pert'  -> perturbative channel, solid above PERT, dashed down to its threshold
CHANNELS = [
    ('e',     r'$e^+e^-$',                 'black',          1.0,       'lep'),
    ('mu',    r'$\mu^+\mu^-$',             'red',            1.0,       'lep'),
    ('pipi',  r'$\pi^+\pi^-$',             '#2ca02c',        2.0 / 3.0, 'dis'),
    ('KK',    r'$K^+K^-$',                 'blue',           0.5,       'dis'),
    ('ss',    r'$s\bar s$',                '#7B68EE',        1.0,       'pert'),
    ('multi', r'multi-meson ($4\pi+\eta\eta+\rho\rho$)', '#FFA500', 1.0, 'dis'),
    ('gg',    r'$gg$',                     'magenta',        1.0,       'pert'),
    ('tau',   r'$\tau^+\tau^-$',           'cyan',           1.0,       'lep'),
    ('cc',    r'$c\bar c$',                '#9ACD32',        1.0,       'pert'),
]
EXTRA_BB = ('bb', r'$b\bar b$', '#8B4513', 1.0, 'pert')

MASS_MIN, MASS_MAX = 0.021, 62.49
YMIN, YMAX = 1e-3, 1.1


def raw_widths(ms):
    """Every channel evaluated at this mass, including the ones the scheme does not use here.

    ss, cc, bb and gg are evaluated with the perturbative expressions at all masses; they are
    zero below their hadronic thresholds and undefined below alpha_s ~ 1 GeV (the running grid
    ends there), in which case zero is returned.
    """
    w = {'e': sw.gamma_ll(ms, sw.M_E), 'mu': sw.gamma_ll(ms, sw.M_MU),
         'tau': sw.gamma_ll(ms, sw.M_TAU), 'gamma': sw.gamma_gammagamma(ms)}
    gp, gk = sw.gamma_pipi_KK(ms)                 # digitised dispersive table, 0 above 2 GeV
    w['pipi'] = gp
    w['KK'] = gk
    w['multi'] = sw.gamma_multimeson(ms)
    calls = (('ss', lambda: sw.gamma_qq(ms, sw.MSB_S_2GEV, 2.0, sw.M_K)),
             ('cc', lambda: sw.gamma_qq(ms, sw.MSB_C_2GEV, 2.0, sw.M_D)),
             ('bb', lambda: sw.gamma_qq(ms, sw.MSB_B_MB, sw.MSB_B_MB, sw.M_B)),
             ('gg', lambda: sw.gamma_gg(ms)))
    for key, fn in calls:
        try:
            w[key] = fn()
        except Exception:
            w[key] = 0.0
    return w


def curve(key, ms, w, mode):
    """Width of one channel, restricted to the mass range in which the model produces it.

    Nothing is drawn outside that range: the dispersive pi pi / K K / 4pi channels exist only
    below the 2 GeV matching scale, and the perturbative s sbar / c cbar / b bbar / g g channels
    only above it.  Each channel therefore begins and ends at the boundary of its own validity
    region, which is what the model actually says.
    """
    if mode == 'pert':
        return (w.get(key, 0.0) if ms >= PERT else 0.0), ms >= PERT
    if mode == 'dis':
        return (w.get(key, 0.0) if ms <= PERT else 0.0), ms <= PERT
    return w.get(key, 0.0), True


def main():
    masses = np.exp(np.linspace(np.log(MASS_MIN), np.log(MASS_MAX), 700))
    allch = CHANNELS + [EXTRA_BB]
    keys = [c[0] for c in allch]
    val = {k: np.zeros(len(masses)) for k in keys}
    used = {k: np.zeros(len(masses), dtype=bool) for k in keys}
    tot = np.zeros(len(masses))

    for i, m in enumerate(masses):
        tot[i] = sw.gamma_total(m)      # the PRODUCTION total, not the sum of everything drawn

    for i, m in enumerate(masses):
        w = raw_widths(m)
        for k, _, _, fac, mode in allch:
            v, u = curve(k, m, w, mode)
            val[k][i] = fac * v / tot[i]
            used[k][i] = u

    with open(os.path.join(DATA, 'bc5_branching_ratios.csv'), 'w', newline='') as f:
        wtr = csv.writer(f)
        wtr.writerow(['mass_GeV'] + keys + [k + '_is_production_input' for k in keys])
        for i, m in enumerate(masses):
            wtr.writerow(['%.6g' % m] + ['%.8e' % val[k][i] for k in keys]
                         + ['%d' % used[k][i] for k in keys])

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--xmin', type=float, default=0.0)
    ap.add_argument('--xmax', type=float, default=4.5)
    ap.add_argument('--xscale', type=str, default='linear', choices=('linear', 'log'))
    ap.add_argument('--tag', type=str, default='')
    ap.add_argument('--show-note', action='store_true', help='draw the explanatory note box')
    args = ap.parse_args()

    def make(name, chans, note=None, xlim=None, xscale=None):
        xlim = xlim if xlim is not None else (args.xmin, args.xmax)
        xscale = xscale if xscale is not None else args.xscale
        fig, ax = plt.subplots(figsize=(8.0, 6.0))
        for k, lab, col, _, _ in chans:
            # every channel is drawn solid over the whole range where it is defined; the
            # _is_production_input column in the CSV records which points are model input
            y = np.where(val[k] > 0, val[k], np.nan)
            ax.plot(masses, y, color=col, lw=1.9, label=lab)
        ax.axvline(PERT, color='0.55', lw=0.9, ls=':')
        xoff = 1.06 if xscale == 'log' else 1.02
        ax.text(PERT * xoff, 1.02, r'$m_S=2$ GeV: scheme switch', fontsize=7, color='0.35')
        ax.set_xscale(xscale); ax.set_yscale('log')
        ax.set_xlim(xlim[0], xlim[1]); ax.set_ylim(YMIN, YMAX)
        ax.set_xlabel(r'$m_S$ [GeV]'); ax.set_ylabel(r'$\mathrm{BR}(S\to X)$')
        ax.set_title(r'BC5 dark scalar branching ratios ($\theta^2=1$ normalised)')
        ax.grid(alpha=0.3, which='both')
        ax.legend(loc='lower right', fontsize=8.5, ncol=1, framealpha=0.95)
        if note:
            ax.text(0.015, 0.02, note, transform=ax.transAxes, fontsize=6.8, va='bottom',
                    bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.92))
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, name + '.png'), dpi=150)
        fig.savefig(os.path.join(FIG, name + '.pdf'))
        plt.close(fig)

    note = ('pi+pi- = 2/3 Gamma_pipi,  K+K- = 1/2 Gamma_KK,  multi-meson = the full toy channel\n'
            'C m_S^3 beta_2pi standing for 4pi + eta eta + rho rho.  pi0pi0, K_L K_S and gamma gamma\n'
            'are omitted.  Each curve covers only the mass range in which the model produces it.')
    # requested view: the 1 - 4.5 GeV window around the scheme matching scale
    note_out = note if args.show_note else None
    make('bc5_branching_ratios', CHANNELS, note_out)
    # preserved log-scale full-range version
    make('bc5_branching_ratios_fullrange', CHANNELS, note,
         xlim=(MASS_MIN, MASS_MAX), xscale='log')
    make('bc5_branching_ratios_with_bb', list(CHANNELS) + [EXTRA_BB],
         'same as bc5_branching_ratios_fullrange plus b bbar, which opens at 2 m_B = 10.56 GeV\n'
         'and carries up to ~86% of the width at the high-mass end.')
    print('canonical figure: x = %.3g - %.3g GeV, %s scale' % (args.xmin, args.xmax, args.xscale))
    print('wrote figures and BC5 branching-ratio CSV')
    for m in (0.1, 0.3, 1.0, 2.0, 5.0, 10.0, 62.49):
        i = int(np.argmin(np.abs(masses - m)))
        print('  m=%7.3f  ' % masses[i] + '  '.join('%s=%.3g%s' % (k, val[k][i], '' if used[k][i] else '*') for k in keys))


if __name__ == '__main__':
    main()
