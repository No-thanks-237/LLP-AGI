# -*- coding: utf-8 -*-
"""Random access to the Pythia 8.316 hadronisation / decay library.

`work/magnus/hadron_lib/hadron_lib_<channel>.txt.gz` holds, for every channel and every mass on
the scan grid, 1500 full Pythia hadronisations of the isolated colour-singlet S decay, with the
stable charged particles given in the S rest frame.  When present it replaces the analytic jet
model and the leading-track tau model entirely.
"""
from __future__ import annotations

import os

import numpy as np

# Mass node tolerance of the library.  The library is generated on the scan grid, so a real
# mismatch means the grid changed and the library must be regenerated.
NODE_MTOL = 2e-3

class HadronLib:
    """Random access to a Pythia-generated S-decay library."""

    def __init__(self, paths):
        self.chan = {}
        for name, path in paths.items():
            d = np.loadtxt(path, comments='#')
            if d.ndim == 1:
                d = d.reshape(1, -1)
            mass, samp, pid, chg = d[:, 0], d[:, 1].astype(np.int64), d[:, 2].astype(int), d[:, 3].astype(int)
            e, px, py, pz = d[:, 4], d[:, 5], d[:, 6], d[:, 7]
            m = np.sqrt(np.maximum(e ** 2 - px ** 2 - py ** 2 - pz ** 2, 0.0))
            # group rows into (mass node, sample) blocks; the file is written mass- then sample-major
            key = mass * 1e6 + samp
            _, start, counts = np.unique(key, return_index=True, return_counts=True)
            order = np.argsort(start)
            self.chan[name] = {
                'mass': mass, 'key': key, 'start': start[order], 'counts': counts[order],
                'e': e, 'px': px, 'py': py, 'pz': pz, 'm': m, 'chg': chg,
                'mass_nodes': np.unique(mass),
                # number of generation ATTEMPTS per node.  Groups whose decay produced no
                # charged particle are absent from the file, so they must be sampled as
                # explicit empty outcomes instead of being normalised away.
                'n_attempt': int(samp.max()) + 1 if len(samp) else 0,
            }

    def sample(self, ms, channel, rng):
        """One random S -> channel decay as [(m, e, px, py, pz, charge), ...] in the S rest frame."""
        c = self.chan.get(channel)
        if c is None:
            return None
        i = int(np.argmin(np.abs(c['mass_nodes'] - ms)))
        node = float(c['mass_nodes'][i])
        rel = abs(node - ms) / max(ms, node, 1e-30)
        if rel > 0.20:
            return None            # channel closed at this mass: the library has no entry
        if rel > NODE_MTOL:
            # N06 (review 2026-10-25): the earlier fallback scaled (m,E,p) of every stable
            # daughter by m_S/m_node.  That matches the parent four-momentum but changes the
            # PHYSICAL masses of pions, kaons and protons, which is not a legitimate on-shell
            # mapping.  The library is now generated on the exact scan grid, so a mismatch means
            # the grid changed: refuse instead of silently producing wrong final states.
            raise RuntimeError(
                'hadron library has no node for %s at m_S = %.6f (nearest %.6f, %.2f%% away); '
                'regenerate the library on this mass grid (work/magnus/bundles_stage_hadronlib3)'
                % (channel, ms, node, 100.0 * rel))
        # Keys are mass*1e6 + sample_index with sample_index < NSAMPLE_MAX.  The upper bound
        # must be the NEXT MASS NODE, not node*1e6 + 1e6: in key units +1e6 is +1 GeV, which
        # swept in every node within 1 GeV above and mixed their spectra into this one.
        nodes = c['mass_nodes']
        nxt = float(nodes[i + 1]) if i + 1 < len(nodes) else node + 1.0
        # Draw one generation ATTEMPT uniformly over the node's attempts.  Attempts whose decay
        # produced no charged particle have no rows in the file, so they must come back as an
        # explicit empty outcome; drawing only among the non-empty groups would renormalise the
        # visible fraction upwards (review P2: 198 of 1500 attempts missing for the gg node).
        n_att = int(c.get('n_attempt', 0))
        if n_att <= 0:
            return None
        tgt = node * 1e6 + rng.randrange(n_att)
        row = int(np.searchsorted(c['key'], tgt - 0.5, side='left'))
        if row >= len(c['key']) or abs(float(c['key'][row]) - tgt) > 0.4:
            return []
        # row is an index into the ROW array; start/counts are indexed by GROUP, so convert
        g = int(np.searchsorted(c['start'], row, side='right')) - 1
        if g < 0 or g >= len(c['start']):
            return []
        s, n = int(c['start'][g]), int(c['counts'][g])
        return [(float(c['m'][k]), float(c['e'][k]), float(c['px'][k]),
                 float(c['py'][k]), float(c['pz'][k]), int(c['chg'][k]))
                for k in range(s, s + n)]


HADRON_LIB = None


def load_hadron_lib(directory):
    """Load hadron_lib_<channel>.txt(.gz) from a directory; returns True on success."""
    global HADRON_LIB
    paths = {}
    for ch in ('gg', 'ss', 'cc', 'bb', 'tautau'):
        for ext in ('.txt.gz', '.txt'):
            f = os.path.join(directory, 'hadron_lib_%s%s' % (ch, ext))
            if os.path.isfile(f):
                paths[ch] = f
                break
    if not paths:
        return False
    HADRON_LIB = HadronLib(paths)
    return True
