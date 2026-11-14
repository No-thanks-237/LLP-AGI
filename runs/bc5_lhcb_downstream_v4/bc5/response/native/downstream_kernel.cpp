// SPDX-License-Identifier: GPL-3.0-or-later
//
// Thin, allocation-free C ABI bridge over the header-only downstream detector
// kernel in detector/downstream/modules.  The headers are consumed in place
// through the compile-time include directory; nothing is copied into this tree.
//
// No C++ exception may cross the C boundary: every entry point catches and
// reports a negative status code instead.
//
// Build (see build.ps1):
//   g++ -std=c++17 -O2 -Wall -Wextra -shared -static -static-libgcc
//   -static-libstdc++ -I E:/AAAwork/LLP-AGI/detector/downstream/modules
//   downstream_kernel.cpp -o downstream_kernel.dll

#include "DownstreamTransport.h"

#include "BirthKinematics.h"

#include <cstddef>
#include <exception>
#include <stdexcept>

#if defined(_WIN32)
#define DS_API extern "C" __declspec(dllexport)
#else
#define DS_API extern "C" __attribute__((visibility("default")))
#endif

namespace {

// Flat result buffer written by ds_transport().  Every entry is a double:
//
//   [0]      status code (identical to the int return value)
//   [1..3]   front  x, y, z        [m]
//   [4..6]   front  px, py, pz     [GeV]
//   [7..9]   end    x, y, z        [m]
//   [10]     front  ct             [m]
//   [11]     end    ct             [m]
//
// Config::validate() enforces fieldZMax < frontZ, so no field acts between the
// front and back faces and the straight front->back leg cannot change px, py,
// pz; the end momentum is therefore identical to the front momentum and is
// reported once.  The transient State at the front is the state used for both
// the front and the (possibly advanced) end position.
constexpr std::size_t kOutDoubles = 12;

// Negative codes returned when the kernel cannot produce a geometry verdict.
constexpr int kErrNullBuffer = -1;
constexpr int kErrInvalidArgument = -2;
constexpr int kErrStdException = -3;
constexpr int kErrUnknownException = -4;

}  // namespace

// Number of doubles ds_transport() writes into its caller-provided out buffer.
DS_API int ds_result_doubles(void) { return static_cast<int>(kOutDoubles); }

// DownstreamTransport.h magneticConstant [GeV / (T m)], exposed so callers can
// reproduce the bend sign without re-declaring the constant.
DS_API double ds_magnetic_constant(void) { return downstream::magneticConstant; }

DS_API int ds_transport(double x, double y, double z, double px, double py, double pz, double energy, double charge,
                        double by, double fieldZMin, double fieldZMax, double frontZ, double backZ, double halfX,
                        double halfY, double holeRadius, int fullEnvelope, double *out) {
  if (out == nullptr) return kErrNullBuffer;
  for (std::size_t i = 0; i < kOutDoubles; ++i) out[i] = 0.0;
  try {
    downstream::State s;
    s.x = x;
    s.y = y;
    s.z = z;
    s.px = px;
    s.py = py;
    s.pz = pz;
    s.energy = energy;
    s.charge = charge;

    downstream::Config c;
    c.by = by;
    c.fieldZMin = fieldZMin;
    c.fieldZMax = fieldZMax;
    c.frontZ = frontZ;
    c.backZ = backZ;
    c.halfX = halfX;
    c.halfY = halfY;
    c.holeRadius = holeRadius;
    c.mode = fullEnvelope ? downstream::AcceptanceMode::FullEnvelope : downstream::AcceptanceMode::FrontFace;

    const downstream::Result r = downstream::transport(s, c);
    const int status = static_cast<int>(r.status);
    out[0] = static_cast<double>(status);
    out[1] = r.front.x;
    out[2] = r.front.y;
    out[3] = r.front.z;
    out[4] = r.front.px;
    out[5] = r.front.py;
    out[6] = r.front.pz;
    out[7] = r.end.x;
    out[8] = r.end.y;
    out[9] = r.end.z;
    out[10] = r.front.ct;
    out[11] = r.end.ct;
    return status;
  } catch (const std::invalid_argument &) {
    return kErrInvalidArgument;
  } catch (const std::exception &) {
    return kErrStdException;
  } catch (...) {
    return kErrUnknownException;
  }
}

// Birth cuts are applied to the ORIGINAL generator momentum, before any
// magnetic propagation (see BirthKinematics.h).  Returns 1/0, or a negative
// code when the thresholds themselves are invalid.
DS_API int ds_passes_birth(double px, double py, double pz, double energy, double minE, double minPT) {
  try {
    downstream::BirthCuts cuts;
    cuts.minEnergy = minE;
    cuts.minPT = minPT;
    return downstream::passesBirthKinematics(px, py, pz, energy, cuts) ? 1 : 0;
  } catch (const std::invalid_argument &) {
    return kErrInvalidArgument;
  } catch (const std::exception &) {
    return kErrStdException;
  } catch (...) {
    return kErrUnknownException;
  }
}
