"""Prepare user-authorized decay probabilities without modifying raw tables.

No SensCalc code is imported or executed. Standard library only.
"""
import argparse
import bisect
import csv
import hashlib
import json
import math
from pathlib import Path

SOURCE = '1809.01876'
CHANNELS = ['e_plus_e_minus', 'mu_plus_mu_minus', 'tau_plus_tau_minus',
            'pi_plus_pi_minus', 'pi0_pi0', 'K_plus_K_minus', 'K_L_K_S',
            '2pi_plus_2pi_minus', 'pi_plus_pi_minus_2pi0',
            'Jets_cc', 'Jets_GG', 'Jets_bb', 'Jets_ss']
# Explicit preparation inputs, not inferred from display strings in the MX.
MASSES_GEV = {'electron': 0.00051099895, 'muon': 0.1056583755,
              'tau': 1.77686, 'pi_charged': 0.13957039, 'pi0': 0.1349768,
              'K_charged': 0.493677, 'K_neutral': 0.497611}
THRESHOLDS = {
    'e_plus_e_minus': 2 * MASSES_GEV['electron'],
    'mu_plus_mu_minus': 2 * MASSES_GEV['muon'],
    'tau_plus_tau_minus': 2 * MASSES_GEV['tau'],
    'pi_plus_pi_minus': 2 * MASSES_GEV['pi_charged'],
    'pi0_pi0': 2 * MASSES_GEV['pi0'],
    'K_plus_K_minus': 2 * MASSES_GEV['K_charged'],
    'K_L_K_S': 2 * MASSES_GEV['K_neutral'],
    '2pi_plus_2pi_minus': 4 * MASSES_GEV['pi_charged'],
    'pi_plus_pi_minus_2pi0': 2 * MASSES_GEV['pi_charged'] + 2 * MASSES_GEV['pi0'],
}
HBAR_C_GEV_M = 1.973269804593025e-16


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_source(path):
    with Path(path).open(encoding='utf-8', newline='') as f:
        rows = [r for r in csv.DictReader(f) if r['source'] == SOURCE]
    if not rows:
        raise ValueError('Selected source is absent')
    rows.sort(key=lambda r: float(r['mass_GeV']))
    masses = [float(r['mass_GeV']) for r in rows]
    if any(not math.isfinite(m) or m <= 0 for m in masses):
        raise ValueError('Invalid mass')
    if any(b <= a for a, b in zip(masses, masses[1:])):
        raise ValueError('Duplicate or unordered mass nodes')
    return rows, masses


def linear(masses, values, mass):
    if not math.isfinite(mass) or not masses[0] <= mass <= masses[-1]:
        raise ValueError('Extrapolation outside data domain is forbidden')
    i = bisect.bisect_left(masses, mass)
    if i < len(masses) and masses[i] == mass:
        return values[i]
    f = (mass - masses[i - 1]) / (masses[i] - masses[i - 1])
    return values[i - 1] * (1 - f) + values[i] * f


class DecayData:
    def __init__(self, branching_path, lifetime_path):
        rows, self.masses = read_source(branching_path)
        self.raw = {c: [float(r[c]) for r in rows] for c in CHANNELS}
        self.cleaned = {c: list(v) for c, v in self.raw.items()}
        self.repairs = []
        for channel in CHANNELS:
            original = self.raw[channel]
            if any(not math.isfinite(v) for v in original):
                raise ValueError('Nonfinite branching ratio')
            threshold = THRESHOLDS.get(channel, 0.0)
            for i, mass in enumerate(self.masses):
                old = original[i]
                if mass <= threshold:
                    new, reason = 0.0, 'closed_exclusive_channel'
                elif old < 0:
                    # Only the observed kaon interpolation defects are supported.
                    if channel not in ('K_plus_K_minus', 'K_L_K_S'):
                        raise ValueError('Unexpected negative channel; no generic clamp allowed')
                    valid = [(m, v) for m, v in zip(self.masses, original)
                             if m > threshold and v >= 0]
                    xs, ys = [threshold] + [x for x, _ in valid], [0.0] + [y for _, y in valid]
                    new = linear(xs, ys, mass)
                    reason = 'linear_bridge_between_nonnegative_kaon_nodes'
                else:
                    continue
                self.cleaned[channel][i] = new
                if new != old:
                    self.repairs.append({'mass_GeV': mass, 'channel': channel,
                                         'original': old, 'prepared': new, 'reason': reason})
        life_rows, self.life_masses = read_source(lifetime_path)
        self.life = [float(r['ctau_reference_m']) for r in life_rows]
        if any(not math.isfinite(v) or v <= 0 for v in self.life):
            raise ValueError('Nonpositive lifetime')
        self.log_life_masses = [math.log(m) for m in self.life_masses]
        self.log_life = [math.log(v) for v in self.life]

    def probabilities(self, mass):
        values = {}
        for c in CHANNELS:
            threshold = THRESHOLDS.get(c, 0.0)
            if mass <= threshold:
                values[c] = 0.0
                continue
            # Insert an exact threshold anchor to avoid feeding a newly open
            # exclusive channel from a below-threshold node.
            i = bisect.bisect_right(self.masses, threshold)
            if i and i < len(self.masses) and mass < self.masses[i]:
                values[c] = self.cleaned[c][i] * (mass - threshold) / (self.masses[i] - threshold)
            else:
                values[c] = linear(self.masses, self.cleaned[c], mass)
        if any(v < 0 or not math.isfinite(v) for v in values.values()):
            raise ValueError('Invalid sampling probability')
        total = math.fsum(values.values())
        if total <= 0:
            raise ValueError('No provided open decay channel at requested mass')
        return {c: v / total for c, v in values.items()}, total

    def ctau_reference(self, mass):
        return math.exp(linear(self.log_life_masses, self.log_life, math.log(mass)))


def prepare(branching, lifetime, destination, minimum=.021, maximum=62.49):
    data = DecayData(branching, lifetime)
    grid = sorted(set([minimum, maximum] + [m for m in data.masses if minimum <= m <= maximum]
                      + [m for m in THRESHOLDS.values() if minimum <= m <= maximum]))
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    try:
        with (destination / 'branching_ratios.csv').open('w', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerow(['source', 'mass_GeV', *CHANNELS, 'sum_before_normalization', 'normalization_factor'])
            max_error = 0.0
            for mass in grid:
                p, total = data.probabilities(mass)
                max_error = max(max_error, abs(math.fsum(p.values()) - 1))
                w.writerow([SOURCE, format(mass, '.17g'),
                            *[format(p[c], '.17g') for c in CHANNELS],
                            format(total, '.17g'), format(1 / total, '.17g')])
        with (destination / 'lifetime.csv').open('w', encoding='utf-8', newline='') as f:
            w = csv.writer(f)
            w.writerow(['source', 'mass_GeV', 'theta_squared_reference', 'ctau_reference_m', 'Gamma_reference_GeV'])
            for mass in grid:
                ctau = data.ctau_reference(mass)
                w.writerow([SOURCE, format(mass, '.17g'), 1, format(ctau, '.17g'), format(HBAR_C_GEV_M / ctau, '.17g')])
        with (destination / 'repairs.csv').open('w', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=['mass_GeV', 'channel', 'original', 'prepared', 'reason'])
            w.writeheader()
            w.writerows(r for r in data.repairs if minimum <= r['mass_GeV'] <= maximum)
        report = {'source': SOURCE, 'mass_min_GeV': minimum, 'mass_max_GeV': maximum,
                  'rows': len(grid), 'normalization_includes_invisible': True,
                  'reference_lifetime_changed': False, 'theta_squared_reference': 1,
                  'lifetime_interpolation': 'log-log linear, no extrapolation',
                  'branching_interpolation': 'linear per channel, exact exclusive thresholds, normalize at query mass',
                  'negative_rule': 'closed channels zero; open kaon defects bridged linearly between nonnegative nodes',
                  'masses_GeV': MASSES_GEV, 'exclusive_thresholds_GeV': THRESHOLDS,
                  'mass_input_note': 'Explicit particle mass values for threshold preparation; generation must use compatible masses.',
                  'max_abs_sum_minus_one': max_error,
                  'repair_count_in_scan': sum(minimum <= r['mass_GeV'] <= maximum for r in data.repairs),
                  'input_sha256': {Path(branching).name: sha(branching), Path(lifetime).name: sha(lifetime)},
                  'output_sha256': {p.name: sha(p) for p in destination.glob('*.csv')},
                  'assumption': 'Normalize the user-selected listed channels; does not establish completeness of physical decays.'}
        (destination / 'preparation_report.json').write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
        print(json.dumps({'prepared': str(destination), 'rows': len(grid), 'repairs': report['repair_count_in_scan'],
                          'max_abs_sum_minus_one': max_error}, ensure_ascii=False))
    except Exception:
        # Leave partial output for diagnosis; never delete or silently overwrite it.
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.input_dir / 'BrRatioScalar_wide.csv', args.input_dir / 'ctauS.csv', args.output_dir)
