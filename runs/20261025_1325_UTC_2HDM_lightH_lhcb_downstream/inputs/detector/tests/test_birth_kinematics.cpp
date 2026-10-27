// Ordinary C++ software assertions, no generators or sensitivity calculation.
#include "../modules/BirthKinematics.h"
#include "../modules/DownstreamTransport.h"
#include <cstdlib>
#include <iostream>
#include <limits>
using namespace downstream;
int checks=0;
void check(bool ok,const char *name) {
  ++checks;if(!ok) { std::cerr<<"FAIL: "<<name<<'\n';std::exit(1); }
}
int main() {
  BirthCuts cuts;cuts.minEnergy=5;cuts.minPT=.5;
  const double abovePT=std::nextafter(.5,1.),belowPT=std::nextafter(.5,0.);
  const double aboveE=std::nextafter(5.,6.),belowE=std::nextafter(5.,4.);
  check(!passesBirthKinematics(belowPT,0,4,6,cuts),"PT below threshold");
  check(!passesBirthKinematics(.5,0,4,6,cuts),"PT equality rejected");
  check(passesBirthKinematics(abovePT,0,4,6,cuts),"PT nextafter above accepted");
  check(!passesBirthKinematics(1,0,4,belowE,cuts),"E below threshold");
  check(!passesBirthKinematics(1,0,4,5,cuts),"E equality rejected");
  check(passesBirthKinematics(1,0,4,aboveE,cuts),"E nextafter above accepted");
  check(passesBirthKinematics(-abovePT,0,4,6,cuts),"negative px accepted by magnitude");
  check(passesBirthKinematics(0,-abovePT,4,6,cuts),"negative py accepted by magnitude");
  check(!passesBirthKinematics(.3,.4,4,6,cuts),"combined PT exact equality rejected");
  // Massive on-shell example: |p|<5, E>5 must pass the ENERGY requirement.
  check(passesBirthKinematics(1,0,4.5,5.1,cuts),"E is not total momentum");
  const double nan=std::numeric_limits<double>::quiet_NaN(),inf=std::numeric_limits<double>::infinity();
  check(!passesBirthKinematics(nan,0,4,6,cuts),"NaN px rejected");
  check(!passesBirthKinematics(1,inf,4,6,cuts),"infinite py rejected");
  check(!passesBirthKinematics(1,0,nan,6,cuts),"NaN pz rejected");
  check(!passesBirthKinematics(1,0,4,inf,cuts),"infinite energy rejected");
  BirthCuts geometry;check(passesBirthKinematics(0,0,20,20.1,geometry),"zero thresholds disable cuts, including PT zero");
  check(!passesBirthKinematics(0,0,20,0,geometry),"nonpositive energy invalid even in geometry mode");
  for(double value:{-1.,inf,nan}) {
    BirthCuts invalid;invalid.minPT=value;bool threw=false;
    try { invalid.validate(); }catch(const std::invalid_argument&) { threw=true; }
    check(threw,"invalid configuration throws");
  }
  State initial;initial.x=.5;initial.y=.2;initial.z=2;initial.pz=20;
  initial.energy=20.1;initial.charge=1;
  const Result after=transport(initial,Config{});
  check(after.accepted() && std::hypot(after.end.px,after.end.py)>1.19,"geometry accepted and magnetic PT above 0.5");
  check(!passesBirthKinematics(initial.px,initial.py,initial.pz,initial.energy,cuts),"birth PT zero still rejected");
  check(passesBirthKinematics(after.end.px,after.end.py,after.end.pz,after.end.energy,cuts),"exit-cut counterexample would wrongly pass");
  // Conversely, birth PT above 0.5 may bend to exit PT below 0.5.
  initial.px=1;initial.energy=20.1;
  const Result cancellation=transport(initial,Config{});
  check(cancellation.accepted() && std::hypot(cancellation.end.px,cancellation.end.py)<.5,"field cancellation geometry");
  check(passesBirthKinematics(initial.px,initial.py,initial.pz,initial.energy,cuts),"birth pass retained despite exit PT below threshold");
  std::cout<<"PASS: "<<checks<<" birth-kinematic software assertions\n";
}
