// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LLP_AGI_BIRTH_KINEMATICS_H
#define LLP_AGI_BIRTH_KINEMATICS_H
#include <cmath>
#include <stdexcept>
namespace downstream {
struct BirthCuts {
  // GeV. Zero explicitly disables that threshold (geometry-only mode).
  double minEnergy=0, minPT=0;
  void validate() const {
    if(!std::isfinite(minEnergy) || !std::isfinite(minPT) || minEnergy<0 || minPT<0)
      throw std::invalid_argument("birth thresholds must be finite and nonnegative");
  }
};
// Apply to original generator momentum BEFORE any magnetic propagation.
// Positive configured thresholds use strict >, with no numerical margin.
inline bool passesBirthKinematics(double px,double py,double pz,double energy,const BirthCuts &cuts) {
  cuts.validate();
  if(!std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) || !std::isfinite(energy) || energy<=0)
    return false;
  const double pt=std::hypot(px,py);
  return std::isfinite(pt) && (cuts.minEnergy==0 || energy>cuts.minEnergy)
    && (cuts.minPT==0 || pt>cuts.minPT);
}
}
#endif
