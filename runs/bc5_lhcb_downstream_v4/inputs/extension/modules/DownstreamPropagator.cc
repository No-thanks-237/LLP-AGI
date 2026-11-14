// SPDX-License-Identifier: GPL-3.0-or-later
#include "modules/DownstreamPropagator.h"
#include "classes/DelphesClasses.h"
#include "TIterator.h"
#include "TObjArray.h"
#include "TString.h"
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

DownstreamPropagator::DownstreamPropagator():fIterator(0),fInputArray(0),fOutputArray(0) {}
DownstreamPropagator::~DownstreamPropagator() { delete fIterator; }
void DownstreamPropagator::Init() {
  fBirthCuts.minEnergy=GetDouble("MinBirthEnergy",0.0);
  fBirthCuts.minPT=GetDouble("MinBirthPT",0.0);
  fBirthCuts.validate();
  fConfig.by=GetDouble("By",1.0);
  fConfig.fieldZMin=GetDouble("FieldZMin",3.0);
  fConfig.fieldZMax=GetDouble("FieldZMax",7.0);
  fConfig.frontZ=GetDouble("SciFiZFront",7.7);
  fConfig.backZ=GetDouble("SciFiZBack",9.4);
  fConfig.halfX=GetDouble("SciFiHalfWidthX",3.24);
  fConfig.halfY=GetDouble("SciFiHalfWidthY",2.415);
  fConfig.holeRadius=GetDouble("BeamHoleRadius",0.09);
  const std::string mode=GetString("AcceptanceMode","full_envelope");
  if(mode=="full_envelope") fConfig.mode=downstream::AcceptanceMode::FullEnvelope;
  else if(mode=="front_face") fConfig.mode=downstream::AcceptanceMode::FrontFace;
  else throw std::invalid_argument("AcceptanceMode must be full_envelope or front_face");
  fConfig.validate();
  fInputArray=ImportArray(GetString("InputArray","Delphes/stableParticles"));
  fIterator=fInputArray->MakeIterator();
  fOutputArray=ExportArray(GetString("OutputArray","tracks"));
}
void DownstreamPropagator::Process() {
  fIterator->Reset();
  Candidate *particle;
  while((particle=static_cast<Candidate *>(fIterator->Next()))) {
    if(particle->Charge==0) continue;
    downstream::State s;
    const TLorentzVector &v=particle->Position, &p=particle->Momentum;
    // Read only the incoming truth four-vector; Track.PT below is an exit value.
    if(!downstream::passesBirthKinematics(p.Px(),p.Py(),p.Pz(),p.E(),fBirthCuts)) continue;
    s.x=v.X()*1e-3; s.y=v.Y()*1e-3; s.z=v.Z()*1e-3; s.ct=v.T()*1e-3;
    s.px=p.Px(); s.py=p.Py(); s.pz=p.Pz(); s.energy=p.E(); s.charge=particle->Charge;
    const downstream::Result result=downstream::transport(s,fConfig);
    if(!result.accepted()) continue;
    const downstream::State &end=result.end, &front=result.front;
    Candidate *track=static_cast<Candidate *>(particle->Clone());
    // Clone may copy an existing association list. TreeWriter requires index 0
    // to reference this original stable truth particle, never another daughter.
    track->GetCandidates()->Clear();
    track->AddCandidate(particle);
    track->InitialPosition=v;
    track->Position.SetXYZT(end.x*1e3,end.y*1e3,end.z*1e3,end.ct*1e3);
    track->Momentum.SetPxPyPzE(end.px,end.py,end.pz,end.energy);
    track->XFirstHit=front.x*1e3; track->YFirstHit=front.y*1e3; track->ZFirstHit=front.z*1e3;
    track->L=end.path*1e3;
    track->P=track->Momentum.P(); track->PT=track->Momentum.Pt();
    track->Phi=track->Momentum.Phi();
    track->CtgTheta=track->PT>0 ? end.pz/track->PT : 1e10;
    // Solenoidal perigee/curvature fields are undefined for this transport.
    // NaN is deliberate; do not interpret them as a reconstructed vertex.
    const float undefined=std::numeric_limits<float>::quiet_NaN();
    track->D0=track->DZ=track->C=undefined;
    track->Xd=track->Yd=track->Zd=undefined;
    fOutputArray->Add(track);
  }
}
void DownstreamPropagator::Finish() { delete fIterator; fIterator=0; }
