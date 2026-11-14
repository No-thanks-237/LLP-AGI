// SPDX-License-Identifier: GPL-3.0-or-later
// Deterministic single-particle transport, not a sensitivity calculation.
#ifndef LLP_AGI_DOWNSTREAM_TRANSPORT_H
#define LLP_AGI_DOWNSTREAM_TRANSPORT_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace downstream {
constexpr double magneticConstant = 0.299792458; // GeV / (T m)
constexpr double pi = 3.14159265358979323846;

// All lengths (including ct and path) are m; momenta and energy are GeV.
struct State {
  double x=0, y=0, z=0, px=0, py=0, pz=0, energy=0;
  double ct=0, path=0, charge=0;
  double momentum() const { return std::hypot(std::hypot(px,py),pz); }
};
enum class AcceptanceMode { FullEnvelope, FrontFace };
struct Config {
  double by=1, fieldZMin=3, fieldZMax=7;
  double frontZ=7.7, backZ=9.4, halfX=3.24, halfY=2.415, holeRadius=0.09;
  AcceptanceMode mode=AcceptanceMode::FullEnvelope;
  void validate() const {
    const double values[]={by,fieldZMin,fieldZMax,frontZ,backZ,halfX,halfY,holeRadius};
    for(double v:values) if(!std::isfinite(v)) throw std::invalid_argument("nonfinite detector parameter");
    if(!(fieldZMin<fieldZMax && fieldZMax<frontZ && frontZ<backZ && halfX>0 && halfY>0 && holeRadius>=0 && holeRadius<std::min(halfX,halfY)))
      throw std::invalid_argument("invalid field/SciFi geometry");
    if(mode!=AcceptanceMode::FullEnvelope && mode!=AcceptanceMode::FrontFace)
      throw std::invalid_argument("invalid acceptance mode");
  }
};
enum class Status { Accepted, InvalidState, BornAfterFront, NotForward, TurnsBeforeExit, OutsideSciFi, BeamHole };
struct Result {
  Status status=Status::InvalidState;
  State front, end;
  bool accepted() const { return status==Status::Accepted; }
};
inline bool valid(const State &s) {
  const double values[]={s.x,s.y,s.z,s.px,s.py,s.pz,s.energy,s.ct,s.path,s.charge};
  for(double v:values) if(!std::isfinite(v)) return false;
  const double p=s.momentum();
  return std::isfinite(p) && p>0 && s.energy>0 && s.energy+1e-10*std::max(p,s.energy)>=p && s.path>=0;
}
inline bool straightToZ(State &s, double z) {
  if(s.pz<=0 || z<s.z) return false;
  const double p=s.momentum(), ds=(z-s.z)*p/s.pz;
  s.x+=ds*s.px/p; s.y+=ds*s.py/p; s.z=z;
  s.ct+=ds*s.energy/p; s.path+=ds;
  return valid(s);
}

// Exact uniform By solution at path length ds. Stable at arbitrarily small B.
inline State arc(const State &s, double by, double ds) {
  const double p=s.momentum(), w=magneticConstant*s.charge*by/p, a=w*ds;
  const double sinc=std::abs(a)<1e-5 ? 1-a*a/6+a*a*a*a/120 : std::sin(a)/a;
  const double cosc=std::abs(a)<1e-5 ? a/2-a*a*a/24+a*a*a*a*a/720 : 2*std::sin(a/2)*std::sin(a/2)/a;
  State out=s;
  out.x+=ds*(s.px*sinc-s.pz*cosc)/p;
  out.y+=ds*s.py/p;
  out.z+=ds*(s.px*cosc+s.pz*sinc)/p;
  out.px=s.px*std::cos(a)-s.pz*std::sin(a);
  out.pz=s.px*std::sin(a)+s.pz*std::cos(a);
  out.path+=ds; out.ct+=ds*s.energy/p;
  return out;
}

// Earliest forward crossing only. A track that turns before the downstream
// boundary is rejected: this module does not reconstruct looping/backward tracks.
inline bool fieldToZ(State &s, double by, double targetZ) {
  if(s.pz<=0 || targetZ<s.z) return false;
  if(targetZ==s.z) return true;
  const double p=s.momentum(), w=magneticConstant*s.charge*by/p;
  if(w==0) return straightToZ(s,targetZ);
  const double theta=std::atan2(s.px,s.pz);
  const double turn=(w>0 ? theta+pi/2 : pi/2-theta)/std::abs(w);
  // Bracket near the requested plane, rather than evaluating an astronomically
  // distant turn for tiny nonzero B (whose path length can overflow double).
  double low=0, high=std::min(turn,(targetZ-s.z)/(s.pz/p));
  State upper=arc(s,by,high);
  while(upper.z<targetZ) {
    if(high>=turn) return false;
    const double next=std::min(turn,2*high);
    if(!std::isfinite(next) || next<=high) return false;
    high=next;upper=arc(s,by,high);
  }
  if(!valid(upper)) return false;
  // Strict inequality excludes a tangent arrival with no forward exit velocity.
  if(high==turn && !(targetZ<upper.z)) return false;
  for(int i=0;i<90;++i) {
    const double mid=low+(high-low)/2;
    if(arc(s,by,mid).z<targetZ) low=mid; else high=mid;
  }
  s=arc(s,by,low+(high-low)/2);
  s.z=targetZ; // remove root-finding roundoff at the hard edge
  return valid(s) && s.pz>0;
}
inline bool insideRectangle(const State &s, const Config &c) {
  return std::abs(s.x)<=c.halfX && std::abs(s.y)<=c.halfY;
}
inline double minimumRadiusSquared(const State &a, const State &b) {
  const double dx=b.x-a.x,dy=b.y-a.y,den=dx*dx+dy*dy;
  const double t=den>0 ? std::max(0.0,std::min(1.0,-(a.x*dx+a.y*dy)/den)) : 0;
  const double x=a.x+t*dx,y=a.y+t*dy;
  return x*x+y*y;
}
inline Result transport(State s, const Config &c) {
  c.validate();
  Result out;
  if(!valid(s)) return out;
  if(s.z>c.frontZ) { out.status=Status::BornAfterFront; return out; }
  if(s.pz<=0) { out.status=Status::NotForward; return out; }
  if(s.z<c.fieldZMin && !straightToZ(s,c.fieldZMin)) return out;
  if(s.z<c.fieldZMax && !fieldToZ(s,c.by,c.fieldZMax)) {
    out.status=Status::TurnsBeforeExit; return out;
  }
  if(!straightToZ(s,c.frontZ)) return out;
  out.front=s; out.end=s;
  if(c.mode==AcceptanceMode::FullEnvelope && !straightToZ(out.end,c.backZ)) return out;
  if(!insideRectangle(out.front,c) || !insideRectangle(out.end,c)) {
    out.status=Status::OutsideSciFi; return out;
  }
  if(minimumRadiusSquared(out.front,out.end)<c.holeRadius*c.holeRadius) {
    out.status=Status::BeamHole; return out;
  }
  out.status=Status::Accepted;
  return out;
}
} // namespace downstream
#endif
