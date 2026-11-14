"""Friendly ctypes wrapper around the compiled downstream transport kernel.

downstream_kernel.dll is a thin C ABI bridge over the header-only C++ kernel in
detector/downstream/modules.  This module locates that DLL next to itself,
declares the flat signatures and exposes a small Python API.

Conventions (DownstreamTransport.h):

* all lengths, including ct and path, are metres;
* momenta and energies are GeV;
* the uniform-field bend parameter is w = 0.299792458 * charge * by / p and
  over a path length ds the exact arc changes the momentum as
      px_out = px * cos(a) - pz * sin(a),    a = w * ds
  so a track with w > 0 and px = 0 bends towards negative x;
* transport() accepts a track only when it is forward (pz > 0), born at or
  upstream of frontZ, still inside the halfX x halfY rectangle on both the
  front and the back face, and clear of the beam hole of radius holeRadius.
  A track born after frontZ (BornAfterFront), that turns before leaving the
  dipole (TurnsBeforeExit), or that passes within holeRadius of the beam axis
  (BeamHole) is rejected.
"""

from __future__ import annotations

import ctypes
import os
from pathlib import Path
from typing import NamedTuple

__all__ = [
    "DLL_PATH",
    "MAGNETIC_CONSTANT",
    "DEFAULT_CONFIG",
    "BIRTH_CUTS",
    "STATUS_NAMES",
    "STATUS_ACCEPTED",
    "TrackPoint",
    "transport_track",
    "transport_raw",
    "passes_birth_kinematics",
]

_DLL_BASENAME = "downstream_kernel.dll"


def _locate_dll() -> Path:
    """Find the shared library, preferring one sitting next to this file."""
    override = os.environ.get("DOWNSTREAM_KERNEL_DLL")
    if override:
        path = Path(override).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError("DOWNSTREAM_KERNEL_DLL points at a missing file: " + str(path))
        return path
    here = Path(__file__).resolve().parent
    for candidate in (here / _DLL_BASENAME, here / "build" / _DLL_BASENAME):
        if candidate.is_file():
            return candidate.resolve()
    raise FileNotFoundError(
        "could not find " + _DLL_BASENAME + " next to " + str(here) + "; run build.ps1 first"
    )


DLL_PATH = _locate_dll()
_LIB = ctypes.CDLL(str(DLL_PATH))

# The kernel writes exactly this many doubles; the Python buffer is the same
# size and the count is verified against the compiled library at import time.
_OUT_DOUBLES = 12
_OUT_ARRAY = ctypes.c_double * _OUT_DOUBLES

_LIB.ds_transport.restype = ctypes.c_int
_LIB.ds_transport.argtypes = [ctypes.c_double] * 16 + [ctypes.c_int, ctypes.POINTER(ctypes.c_double)]

_LIB.ds_passes_birth.restype = ctypes.c_int
_LIB.ds_passes_birth.argtypes = [ctypes.c_double] * 6

_LIB.ds_result_doubles.restype = ctypes.c_int
_LIB.ds_result_doubles.argtypes = []

_LIB.ds_magnetic_constant.restype = ctypes.c_double
_LIB.ds_magnetic_constant.argtypes = []

_KERNEL_OUT_DOUBLES = int(_LIB.ds_result_doubles())
if _KERNEL_OUT_DOUBLES > _OUT_DOUBLES:
    raise RuntimeError(
        "kernel writes " + str(_KERNEL_OUT_DOUBLES)
        + " doubles but the Python buffer only holds " + str(_OUT_DOUBLES)
    )

MAGNETIC_CONSTANT = float(_LIB.ds_magnetic_constant())

# Production acceptance card.
DEFAULT_CONFIG = {
    "by": 1.0,
    "fieldZMin": 3.0,
    "fieldZMax": 7.0,
    "frontZ": 7.7,
    "backZ": 9.4,
    "halfX": 3.24,
    "halfY": 2.415,
    "holeRadius": 0.09,
    "fullEnvelope": True,
}

# Production birth cuts [GeV].  Zero disables the corresponding threshold.
BIRTH_CUTS = {"minE": 5.0, "minPT": 0.5}

STATUS_NAMES = {
    0: "Accepted",
    1: "InvalidState",
    2: "BornAfterFront",
    3: "NotForward",
    4: "TurnsBeforeExit",
    5: "OutsideSciFi",
    6: "BeamHole",
    -1: "KernelErrorNullBuffer",
    -2: "KernelErrorInvalidArgument",
    -3: "KernelErrorException",
    -4: "KernelErrorUnknown",
}

STATUS_ACCEPTED = 0


class TrackPoint(NamedTuple):
    """A single track state; a plain tuple with named fields."""
    x: float
    y: float
    z: float
    px: float
    py: float
    pz: float


def _resolve_config(overrides):
    unknown = sorted(set(overrides) - set(DEFAULT_CONFIG))
    if unknown:
        raise KeyError("unknown config override(s): " + ", ".join(unknown))
    config = dict(DEFAULT_CONFIG)
    config.update(overrides)
    return config


def _dispatch(x, y, z, px, py, pz, energy, charge, config):
    buffer = _OUT_ARRAY()
    status = int(
        _LIB.ds_transport(
            float(x), float(y), float(z),
            float(px), float(py), float(pz),
            float(energy), float(charge),
            float(config["by"]),
            float(config["fieldZMin"]), float(config["fieldZMax"]),
            float(config["frontZ"]), float(config["backZ"]),
            float(config["halfX"]), float(config["halfY"]),
            float(config["holeRadius"]),
            1 if config["fullEnvelope"] else 0,
            buffer,
        )
    )
    return status, buffer


def _status_name(status):
    return STATUS_NAMES.get(status, "Unknown(" + str(status) + ")")


def transport_raw(x, y, z, px, py, pz, energy, charge, **cfg_overrides):
    """Call the kernel and return (status, status_name, raw_buffer_bytes).

    The raw bytes are the exact 12-double buffer, useful for determinism tests.
    """
    config = _resolve_config(cfg_overrides)
    status, buffer = _dispatch(x, y, z, px, py, pz, energy, charge, config)
    return status, _status_name(status), bytes(buffer)


def transport_track(x, y, z, px, py, pz, energy, charge, **cfg_overrides):
    """Transport one track and return (accepted, status_name, front, end).

    front and end are TrackPoint tuples.  The kernel reports both endpoint
    positions and a single three-momentum: Config::validate() forces
    fieldZMax < frontZ, so no field acts between the front and back faces and
    px, py, pz are identical there.  When the status is not a geometry verdict
    (InvalidState, BornAfterFront, NotForward, TurnsBeforeExit) both endpoints
    are the all-zero default.
    """
    config = _resolve_config(cfg_overrides)
    status, buffer = _dispatch(x, y, z, px, py, pz, energy, charge, config)
    front = TrackPoint(
        float(buffer[1]), float(buffer[2]), float(buffer[3]),
        float(buffer[4]), float(buffer[5]), float(buffer[6]),
    )
    end = TrackPoint(
        float(buffer[7]), float(buffer[8]), float(buffer[9]),
        float(buffer[4]), float(buffer[5]), float(buffer[6]),
    )
    return status == STATUS_ACCEPTED, _status_name(status), front, end


def passes_birth_kinematics(px, py, pz, energy, minE, minPT):
    """Strict birth cuts from BirthKinematics.h: energy > minE and
    hypot(px, py) > minPT, where a zero threshold disables that cut.

    Raises ValueError when the thresholds themselves are invalid.
    """
    code = int(
        _LIB.ds_passes_birth(
            float(px), float(py), float(pz),
            float(energy), float(minE), float(minPT),
        )
    )
    if code < 0:
        raise ValueError("kernel rejected the birth thresholds: " + _status_name(code))
    return code == 1
