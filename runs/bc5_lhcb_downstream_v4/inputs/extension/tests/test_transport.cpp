// Standalone software tests. No ROOT, generators, model events or sensitivity.
#include "../modules/DownstreamTransport.h"
#include <cstdlib>
#include <iostream>
#include <limits>
#include <random>
#include <string>
using namespace downstream;
int checks=0;
void check(bool ok,const std::string &name) {
  ++checks; if(!ok) { std::cerr<<"FAIL: "<<name<<'\n'; std::exit(1); }
}
void close(double a,double b,double eps,const std::string &name) {
  check(std::abs(a-b)<=eps*std::max(1.0,std::max(std::abs(a),std::abs(b))),name);
}
State initial(double px=0,double py=0,double pz=10,double z=2,double x=0.5) {
  State s; s.x=x;s.z=z;s.px=px;s.py=py;s.pz=pz;
  s.energy=std::hypot(s.momentum(),0.105);s.charge=1;s.ct=0.17;return s;
}
// Independent RK4 using z rather than path as the integration variable.
State derivative(const State &s,double by) {
  State d; const double kb=magneticConstant*s.charge*by;
  d.x=s.px/s.pz;d.y=s.py/s.pz;d.px=-kb;d.pz=kb*s.px/s.pz;
  d.path=s.momentum()/s.pz;d.ct=s.energy/s.pz;return d;
}
State add(State s,const State &d,double a) {
  s.x+=a*d.x;s.y+=a*d.y;s.px+=a*d.px;s.pz+=a*d.pz;
  s.path+=a*d.path;s.ct+=a*d.ct;return s;
}
State rk4(State s,double by,double zend) {
  const int n=4000;const double h=(zend-s.z)/n;
  for(int i=0;i<n;++i) {
    const State a=derivative(s,by),b=derivative(add(s,a,h/2),by);
    const State c=derivative(add(s,b,h/2),by),d=derivative(add(s,c,h),by);
    s=add(s,a,h/6);s=add(s,b,h/3);s=add(s,c,h/3);s=add(s,d,h/6);s.z+=h;
  }return s;
}
int main() {
  Config cfg;cfg.validate();
  std::mt19937 rng(271828);std::uniform_real_distribution<double> u(-1,1);
  for(int i=0;i<80;++i) {
    State input=initial(4*u(rng),3*u(rng),8+20*(u(rng)+1)/2,3);
    input.charge=i%2 ? -1:1;const double by=i%3 ? 1:-1;
    State exact=input;check(fieldToZ(exact,by,7),"forward crossing");
    State numeric=rk4(input,by,7);
    close(exact.x,numeric.x,2e-10,"RK x");close(exact.y,numeric.y,2e-10,"RK y");
    close(exact.px,numeric.px,2e-10,"RK px");close(exact.pz,numeric.pz,2e-10,"RK pz");
    close(exact.path,numeric.path,2e-10,"RK path");close(exact.ct,numeric.ct,2e-10,"RK ct");
    close(exact.momentum(),input.momentum(),2e-14,"momentum invariant");
    close(exact.py,input.py,0,"py invariant");close(exact.energy,input.energy,0,"energy invariant");
  }
  State positive=initial(0,0,20,3,0),negative=positive;
  negative.charge=-1;check(fieldToZ(positive,1,7)&&fieldToZ(negative,1,7),"polarity propagation");
  check(positive.px<0 && positive.x<0,"positive charge bends negative x");
  close(positive.x,-negative.x,1e-13,"charge flip x");close(positive.px,-negative.px,1e-13,"charge flip px");
  State flipped=initial(0,0,20,3,0);check(fieldToZ(flipped,-1,7),"B flip");close(flipped.x,negative.x,1e-13,"B/charge symmetry");
  State zero=initial(1,2,10,3),neutral=zero,tiny=zero;neutral.charge=0;
  check(fieldToZ(zero,0,7)&&fieldToZ(neutral,1,7)&&fieldToZ(tiny,1e-14,7),"zero/neutral/tiny field");
  close(zero.x,0.9,1e-13,"straight x");close(zero.y,0.8,1e-13,"straight y");
  close(neutral.x,zero.x,1e-13,"neutral straight");close(tiny.x,zero.x,1e-12,"tiny B stable");
  for(double by:{1e-100,1e-300,1e-320}) {
    State extremelyTiny=initial(1,2,10,3);
    check(fieldToZ(extremelyTiny,by,7),"subnormal field crossing");
    close(extremelyTiny.x,zero.x,1e-13,"subnormal field straight limit");
  }
  State turn=initial(0,0,0.1,3);check(!fieldToZ(turn,1,7),"turn before field exit");
  State tangent=initial(0,0,1,3);check(!fieldToZ(tangent,1,3+1/magneticConstant),"tangent exit rejected");
  State backward=initial(0,0,-10);check(transport(backward,cfg).status==Status::NotForward,"backward rejected");
  check(transport(initial(0,0,10,7.8),cfg).status==Status::BornAfterFront,"born behind front rejected");
  State nan=initial();nan.x=std::numeric_limits<double>::quiet_NaN();check(!transport(nan,cfg).accepted(),"NaN rejected");
  State unphysical=initial();unphysical.energy=1;check(!transport(unphysical,cfg).accepted(),"E below p rejected");
  cfg.by=0;
  Result straight=transport(initial(),cfg);check(straight.accepted(),"straight accepted");
  close(straight.front.z,7.7,0,"front plane");close(straight.end.z,9.4,0,"back plane");
  close(straight.end.path,7.4,1e-14,"total path");
  close(straight.end.ct,0.17+7.4*initial().energy/10,1e-14,"massive ct");
  check(transport(initial(0,0,10,2,0),cfg).status==Status::BeamHole,"front hole");
  check(transport(initial(0,0,10,2,3.3),cfg).status==Status::OutsideSciFi,"rectangle reject");
  check(transport(initial(0,0,10,2,3.24),cfg).accepted(),"rectangle boundary inclusive");
  State crossing=initial(-2,0,10,7.7,0.17);
  check(transport(crossing,cfg).status==Status::BeamHole,"intermediate hole crossing despite good endpoints");
  cfg.mode=AcceptanceMode::FrontFace;check(transport(crossing,cfg).accepted(),"front-only explicit mode");
  close(transport(crossing,cfg).end.z,7.7,0,"front-only endpoint");
  cfg.mode=AcceptanceMode::FullEnvelope;
  check(transport(initial(1,0,10,7.7,3.2),cfg).status==Status::OutsideSciFi,"back rectangle fails");
  cfg.by=1;
  Result bornIn=transport(initial(0,0,20,5),cfg);check(bornIn.accepted(),"born in field");
  close(bornIn.end.px,-2*magneticConstant,1e-12,"field exposure from birth");
  Result bornAfter=transport(initial(0,0,20,7),cfg);check(bornAfter.accepted(),"born at exit");
  close(bornAfter.end.px,0,0,"no field after exit");
  Config invalid=cfg;invalid.fieldZMax=8;bool threw=false;try{invalid.validate();}catch(const std::invalid_argument&){threw=true;}
  check(threw,"invalid geometry fails fast");
  std::cout<<"PASS: "<<checks<<" deterministic software assertions; 80 independent RK4 comparisons\n";
}
