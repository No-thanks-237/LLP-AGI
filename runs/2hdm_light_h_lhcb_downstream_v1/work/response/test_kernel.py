#!/usr/bin/env python
"""Self-contained validation suite for downstream_ctypes / downstream_kernel.dll.

Usage:
    D:/py/python.exe test_kernel.py

transport_track() returns (accepted: bool, status_name: str, front, end).
Prints one PASS/FAIL line per assertion and exits 0 only if every assertion
passes.  The independent reference in this file re-implements the exact uniform
field arc and integrates it with a fine-step composition, so the DLL is checked
against something other than itself.
"""

from __future__ import annotations

import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import downstream_ctypes as dk

K = dk.MAGNETIC_CONSTANT
FRONT_Z = dk.DEFAULT_CONFIG["frontZ"]
BACK_Z = dk.DEFAULT_CONFIG["backZ"]
FIELD_ZMIN = dk.DEFAULT_CONFIG["fieldZMin"]
FIELD_ZMAX = dk.DEFAULT_CONFIG["fieldZMax"]

_RESULTS = []


def check(name, ok, detail=""):
    ok = bool(ok)
    _RESULTS.append((name, ok))
    line = ("PASS: " if ok else "FAIL: ") + name
    if detail:
        line += "  [" + detail + "]"
    print(line, flush=True)


def _raises(exc, fn):
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


# ---------------------------------------------------------------------------
# Independent reference: exact uniform-field arc plus a fine-step composition.
# State tuple: (x, y, z, px, py, pz, energy, ct, path, charge).
# ---------------------------------------------------------------------------

def _mom(s):
    return math.sqrt(s[3] * s[3] + s[4] * s[4] + s[5] * s[5])


def _straight(s, z_target):
    x, y, z, px, py, pz, energy, ct, path, charge = s
    p = _mom(s)
    ds = (z_target - z) * p / pz
    return (x + ds * px / p, y + ds * py / p, z_target, px, py, pz,
            energy, ct + ds * energy / p, path + ds, charge)


def _arc(s, by, ds):
    x, y, z, px, py, pz, energy, ct, path, charge = s
    p = _mom(s)
    a = K * charge * by / p * ds
    if abs(a) < 1e-5:
        sinc = 1 - a * a / 6 + a ** 4 / 120
        cosc = a / 2 - a ** 3 / 24 + a ** 5 / 720
    else:
        sinc = math.sin(a) / a
        cosc = 2 * math.sin(a / 2) ** 2 / a
    return (x + ds * (px * sinc - pz * cosc) / p,
            y + ds * py / p,
            z + ds * (px * cosc + pz * sinc) / p,
            px * math.cos(a) - pz * math.sin(a),
            py,
            px * math.sin(a) + pz * math.cos(a),
            energy, ct + ds * energy / p, path + ds, charge)


def _fine_step_chain(s, by, ds_total, steps):
    """Compose many small exact arcs: an independent integrator of the ODE."""
    step = ds_total / steps
    current = s
    for _ in range(steps):
        current = _arc(current, by, step)
    return current


def _path_to_z(s, by, z_target, iterations=160):
    """Path length needed to reach z_target, by bisection on the exact arc."""
    p = _mom(s)
    w = K * s[9] * by / p
    if w == 0.0:
        return (z_target - s[2]) * p / s[5]
    high = (z_target - s[2]) * p / s[5]
    if not high > 0:
        return None
    for _ in range(400):
        if _arc(s, by, high)[2] >= z_target:
            break
        high *= 2.0
        if not math.isfinite(high):
            return None
    else:
        return None
    low = 0.0
    for _ in range(iterations):
        mid = 0.5 * (low + high)
        if _arc(s, by, mid)[2] < z_target:
            low = mid
        else:
            high = mid
    return 0.5 * (low + high)


def reference_front(x, y, z, px, py, pz, energy, charge, by,
                    front_z=FRONT_Z, field_zmin=FIELD_ZMIN, field_zmax=FIELD_ZMAX):
    state = (x, y, z, px, py, pz, energy, 0.0, 0.0, charge)
    if state[2] < field_zmin:
        state = _straight(state, field_zmin)
    before = state
    field_path = 0.0
    if state[2] < field_zmax:
        field_path = _path_to_z(state, by, field_zmax)
        if field_path is None:
            return None
        state = _arc(state, by, field_path)
        state = state[:2] + (field_zmax,) + state[3:]
    front = _straight(state, front_z)
    return {"front": front, "before": before, "field_path": field_path}


# ---------------------------------------------------------------------------
# (a) straight neutral-charge-like track
# ---------------------------------------------------------------------------
print("-- (a) straight neutral track --", flush=True)

pz_hi = 100.0
energy_hi = math.hypot(pz_hi, 0.105)
accepted, name, front, end = dk.transport_track(0.0, 0.0, 0.0, 0.0, 0.0, pz_hi, energy_hi, 0.0)
check("a1 charge-0 track along +z reaches the SciFi front plane",
      name in ("Accepted", "OutsideSciFi", "BeamHole") and abs(front.z - FRONT_Z) < 1e-12,
      "front.z=" + repr(front.z) + " status=" + name)
check("a2 charge-0 track at x=y=0 is accepted",
      accepted is True,
      "status=" + name + ": the axis lies inside the beam hole, so the kernel returns BeamHole")
accepted, name, front, end = dk.transport_track(0.15, 0.0, 0.0, 0.0, 0.0, pz_hi, energy_hi, 0.0)
check("a3 charge-0 track offset 0.15 m is accepted (clears the hole)",
      accepted is True and abs(front.x - 0.15) < 1e-12 and abs(front.z - FRONT_Z) < 1e-12,
      "status=" + name + " front.x=" + repr(front.x))
accepted, name, front, end = dk.transport_track(0.0, 0.0, 0.0, 0.0, 0.0, pz_hi, energy_hi, 0.0, holeRadius=0.0)
check("a4 the same axis track is accepted once the hole is disabled (holeRadius=0)",
      accepted is True and abs(front.x) < 1e-12 and abs(front.z - FRONT_Z) < 1e-12,
      "status=" + name + ": isolates the beam hole as the sole cause of the a2 rejection")

# ---------------------------------------------------------------------------
# (b) zero-field limit
# ---------------------------------------------------------------------------
print("-- (b) zero-field limit --", flush=True)

px, py, pz = 1.0, 0.5, 20.0
energy = math.hypot(math.hypot(px, py, pz), 0.105)
accepted, name, front, end = dk.transport_track(0.5, 0.2, 0.0, px, py, pz, energy, 1.0, by=0.0)
check("b1 charged track in zero field is accepted", accepted is True, "status=" + name)
expected_x = 0.5 + FRONT_Z * px / pz
expected_y = 0.2 + FRONT_Z * py / pz
check("b2 zero-field front equals the analytic straight line",
      abs(front.x - expected_x) < 1e-12 and abs(front.y - expected_y) < 1e-12,
      "x=" + repr(front.x) + " vs " + repr(expected_x) + ", y=" + repr(front.y) + " vs " + repr(expected_y))
check("b3 zero-field momentum is untouched",
      abs(front.px - px) < 1e-12 and abs(front.py - py) < 1e-12 and abs(front.pz - pz) < 1e-12,
      "p=" + repr((front.px, front.py, front.pz)))
check("b4 zero-field end equals the straight line at backZ",
      abs(end.x - (0.5 + BACK_Z * px / pz)) < 1e-12 and abs(end.z - BACK_Z) < 1e-12,
      "end=(" + repr(end.x) + ", " + repr(end.z) + ")")

# ---------------------------------------------------------------------------
# (c) born inside the SciFi
# ---------------------------------------------------------------------------
print("-- (c) born inside the SciFi --", flush=True)
accepted, name, front, end = dk.transport_track(0.5, 0.2, 8.0, 0.0, 0.0, 20.0, math.hypot(20.0, 0.105), 1.0)
check("c track born at z=8.0 m is BornAfterFront", name == "BornAfterFront", "status=" + name)

# ---------------------------------------------------------------------------
# (d) forward requirement
# ---------------------------------------------------------------------------
print("-- (d) forward requirement --", flush=True)
accepted, name, front, end = dk.transport_track(0.5, 0.2, 0.0, 1.0, 0.0, 0.0, 1.05, 1.0)
check("d1 pz=0 is NotForward", name == "NotForward", "status=" + name)
accepted, name, front, end = dk.transport_track(0.5, 0.2, 0.0, 1.0, 0.0, -5.0, 5.2, 1.0)
check("d2 pz<0 is NotForward", name == "NotForward", "status=" + name)

# ---------------------------------------------------------------------------
# (e) SciFi rectangle
# ---------------------------------------------------------------------------
print("-- (e) SciFi rectangle --", flush=True)
energy20 = math.hypot(20.0, 0.105)
accepted, name, front, end = dk.transport_track(3.30, 0.0, 0.0, 0.0, 0.0, 20.0, energy20, 0.0)
check("e1 x=3.30 m is OutsideSciFi", name == "OutsideSciFi", "status=" + name)
accepted, name, front, end = dk.transport_track(3.20, 0.0, 0.0, 0.0, 0.0, 20.0, energy20, 0.0)
check("e2 x=3.20 m is accepted", accepted is True, "status=" + name)

# ---------------------------------------------------------------------------
# (f) beam hole
# ---------------------------------------------------------------------------
print("-- (f) beam hole --", flush=True)
energy10 = math.hypot(10.0, 0.105)
accepted, name, front, end = dk.transport_track(0.05, 0.0, 0.0, 0.01, 0.0, 10.0, energy10, 0.0)
check("f1 x=0.05 m with a tiny angle is BeamHole", name == "BeamHole", "status=" + name)
accepted, name, front, end = dk.transport_track(0.15, 0.0, 0.0, 0.01, 0.0, 10.0, energy10, 0.0)
check("f2 x=0.15 m with a tiny angle is accepted", accepted is True,
      "status=" + name + " front.x=" + repr(front.x))

# ---------------------------------------------------------------------------
# (g) bend sign convention and independent arc integration
# ---------------------------------------------------------------------------
print("-- (g) bend sign and independent arc integration --", flush=True)

# Down-quark-like: charge = -1/3, by = -1, so K * charge * by > 0 and the
# header predicts px_out = px*cos(a) - pz*sin(a) < 0 for px = 0.
q_down = -1.0 / 3.0
by_down = -1.0
px0, py0, pz0 = 0.0, 0.0, 2.0
energy0 = math.hypot(pz0, 0.0047)
ref = reference_front(0.0, 0.0, 0.0, px0, py0, pz0, energy0, q_down, by_down)
p0 = math.hypot(px0, py0, pz0)
w = K * q_down * by_down / p0
a = w * ref["field_path"]
predicted_px = px0 * math.cos(a) - pz0 * math.sin(a)
accepted, name, front, end = dk.transport_track(0.0, 0.0, 0.0, px0, py0, pz0, energy0, q_down, by=by_down)
check("g1 down-quark-like track (K*charge*by > 0) bends towards -x",
      w > 0 and predicted_px < 0 and front.x < 0 and front.px < 0,
      "w=" + repr(w) + " a=" + repr(a) + " x=" + repr(front.x) + " px=" + repr(front.px))
check("g2 kernel front momentum matches px*cos(a)-pz*sin(a)",
      abs(front.px - predicted_px) < 1e-9,
      "px=" + repr(front.px) + " predicted=" + repr(predicted_px))
check("g3 kernel front matches the independent reference (down-quark-like)",
      abs(front.x - ref["front"][0]) < 1e-6 and abs(front.y - ref["front"][1]) < 1e-6
      and abs(front.z - ref["front"][2]) < 1e-12,
      "dx=" + repr(front.x - ref["front"][0]))

accepted_m, name_m, front_m, end_m = dk.transport_track(0.0, 0.0, 0.0, px0, py0, pz0, energy0, -q_down, by=by_down)
check("g4 flipping the charge mirrors the bend in x",
      front_m.x > 0 and abs(front_m.x + front.x) < 1e-12,
      "x(+)=" + repr(front_m.x) + " x(-)=" + repr(front.x))

accepted_p, name_p, front_p, end_p = dk.transport_track(
    0.0, 0.0, 0.0, 0.0, 0.0, 2.0, math.hypot(2.0, 0.105), 1.0, by=1.0)
check("g5 positive charge in positive field also bends towards -x",
      front_p.x < 0 and front_p.px < 0,
      "x=" + repr(front_p.x) + " px=" + repr(front_p.px))

rng = random.Random(271828)
worst_front = 0.0
worst_chain = 0.0
random_ok = True
for index in range(6):
    ref = None
    for _attempt in range(200):
        x0 = rng.uniform(-1.0, 1.0)
        y0 = rng.uniform(-1.0, 1.0)
        z0 = rng.uniform(0.0, 3.0)
        rpx = rng.uniform(-2.0, 2.0)
        rpy = rng.uniform(-2.0, 2.0)
        rpz = rng.uniform(5.0, 50.0)
        renergy = math.hypot(math.hypot(rpx, rpy), rpz) * 1.0001
        rcharge = rng.choice([-1.0, 1.0, -1.0 / 3.0, 2.0 / 3.0])
        rby = rng.choice([-1.0, 1.0])
        ref = reference_front(x0, y0, z0, rpx, rpy, rpz, renergy, rcharge, rby)
        if ref is not None:
            break
    if ref is None:
        check("g6 random track " + str(index) + " reference solved", False, "no crossing")
        random_ok = False
        continue

    exact_arc = _arc(ref["before"], rby, ref["field_path"])
    chain = _fine_step_chain(ref["before"], rby, ref["field_path"], 20000)
    chain_diff = max(abs(chain[i] - exact_arc[i]) for i in (0, 1, 2, 3, 5))
    worst_chain = max(worst_chain, chain_diff)

    accepted, name, front, end = dk.transport_track(x0, y0, z0, rpx, rpy, rpz, renergy, rcharge, by=rby)
    if name not in ("Accepted", "OutsideSciFi", "BeamHole"):
        random_ok = False
        check("g6 random track " + str(index) + " yielded a front state", False, "status=" + name)
        continue
    diff = max(abs(front.x - ref["front"][0]), abs(front.y - ref["front"][1]),
               abs(front.z - ref["front"][2]))
    worst_front = max(worst_front, diff)

check("g6 fine-step arc composition agrees with the exact arc",
      random_ok and worst_chain < 1e-9, "worst|chain-exact|=" + repr(worst_chain))
check("g7 random SciFi-front crossing agrees with the reference to < 1e-6 m",
      random_ok and worst_front < 1e-6, "worst|diff|=" + repr(worst_front))

# ---------------------------------------------------------------------------
# (h) strict birth thresholds
# ---------------------------------------------------------------------------
print("-- (h) strict birth thresholds --", flush=True)
check("h1 BIRTH_CUTS matches the production card",
      dk.BIRTH_CUTS == {"minE": 5.0, "minPT": 0.5}, repr(dk.BIRTH_CUTS))
check("h2 energy exactly 5.0 is rejected",
      dk.passes_birth_kinematics(1.0, 0.0, 4.0, 5.0, 5.0, 0.5) is False)
check("h3 energy 5.0000001 is accepted",
      dk.passes_birth_kinematics(1.0, 0.0, 4.0, 5.0000001, 5.0, 0.5) is True)
check("h4 pT exactly 0.5 is rejected",
      dk.passes_birth_kinematics(0.5, 0.0, 4.0, 6.0, 5.0, 0.5) is False)
check("h5 pT 0.5000001 is accepted",
      dk.passes_birth_kinematics(0.5000001, 0.0, 4.0, 6.0, 5.0, 0.5) is True)
check("h6 invalid thresholds raise instead of crossing the C boundary",
      _raises(ValueError, lambda: dk.passes_birth_kinematics(1.0, 0.0, 4.0, 6.0, -1.0, 0.5)))

# ---------------------------------------------------------------------------
# (i) determinism
# ---------------------------------------------------------------------------
print("-- (i) determinism --", flush=True)
raw_energy = math.hypot(math.hypot(1.3, 0.7), 8.0) * 1.0001
status1, name1, bytes1 = dk.transport_raw(0.4, 0.3, 0.0, 1.3, -0.7, 8.0, raw_energy, 1.0)
status2, name2, bytes2 = dk.transport_raw(0.4, 0.3, 0.0, 1.3, -0.7, 8.0, raw_energy, 1.0)
check("i1 repeated calls return byte-identical buffers",
      status1 == status2 and name1 == name2 and bytes1 == bytes2 and len(bytes1) == 12 * 8,
      "status=" + name1 + " len=" + str(len(bytes1)))
status3, name3, bytes3 = dk.transport_raw(0.4, 0.3, 0.0, 1.3, -0.7, 8.0,
                                          math.hypot(math.hypot(1.3, 0.7), 8.0) * 1.0001, 1.0)
check("i2 determinism survives argument reconstruction", bytes1 == bytes3, "status=" + name3)

# ---------------------------------------------------------------------------
# (j) error handling across the C boundary
# ---------------------------------------------------------------------------
print("-- (j) error handling --", flush=True)
status, name, raw = dk.transport_raw(0.0, 0.0, 0.0, 0.0, 0.0, 20.0, 20.1, 1.0, fieldZMax=8.0)
check("j1 invalid geometry returns a negative status without crashing",
      status == -2, "status=" + name)
accepted, name, front, end = dk.transport_track(float("nan"), 0.0, 0.0, 0.0, 0.0, 20.0, 20.1, 1.0)
check("j2 NaN state is InvalidState", name == "InvalidState", "status=" + name)
accepted, name, front, end = dk.transport_track(0.5, 0.2, 0.0, 0.0, 0.0, 20.0, 1.0, 1.0)
check("j3 energy below momentum is InvalidState", name == "InvalidState", "status=" + name)

# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------
passed = sum(1 for _n, ok in _RESULTS if ok)
failed = [n for n, ok in _RESULTS if not ok]
print("", flush=True)
print("=" * 72, flush=True)
print("SUMMARY: " + str(passed) + "/" + str(len(_RESULTS)) + " checks passed", flush=True)
for name in failed:
    print("  FAILED: " + name, flush=True)
if failed:
    print("", flush=True)
    print("NOTE: a2 is a specification conflict, not a kernel defect.  A straight", flush=True)
    print("      track at x=y=0 lies exactly on the beam axis, so the front->back", flush=True)
    print("      segment passes within holeRadius=0.09 m of the origin and", flush=True)
    print("      DownstreamTransport.h correctly reports BeamHole (requirement f).", flush=True)
print("EXIT: " + ("0" if not failed else "1"), flush=True)
sys.exit(0 if not failed else 1)
