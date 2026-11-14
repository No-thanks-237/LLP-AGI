// Delphes interface unit test with manually constructed Candidate objects.
// No event generator, model sample, HepMC reader or sensitivity calculation.
#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "TApplication.h"
#include "TBranch.h"
#include "TClass.h"
#include "TClonesArray.h"
#include "TMemFile.h"
#include "TObjArray.h"
#include "TROOT.h"
#include "TTree.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
void require(bool ok,const char *name) { if(!ok) throw std::runtime_error(name); }
void near(double a,double b,double tol,const char *name) { require(std::abs(a-b)<tol,name); }
int main(int argc,char **argv) {
  if(argc!=2) { std::cerr<<"usage: test_module card.tcl\n";return 2; }
  try {
    require(TClass::GetClass("DownstreamPropagator")!=nullptr,"ROOT module dictionary missing");
    gROOT->SetBatch();
    int appargc=1;char *appargv[]={argv[0]};TApplication app("test_module",&appargc,appargv);
    TMemFile memory("test_module.root","RECREATE");
    ExRootTreeWriter writer(&memory,"Delphes");
    ExRootConfReader reader;reader.ReadFile(argv[1]);
    Delphes delphes("Delphes");delphes.SetConfReader(&reader);delphes.SetTreeWriter(&writer);
    TObjArray *all=delphes.ExportArray("allParticles"), *stable=delphes.ExportArray("stableParticles");
    delphes.ExportArray("partons");delphes.InitTask();
    Candidate *p=delphes.GetFactory()->NewCandidate();
    p->PID=-13;p->Status=1;p->Charge=1;p->Mass=0.105;
    p->Position.SetXYZT(500,200,2000,170);
    p->Momentum.SetPxPyPzE(0,0,20,std::hypot(20.0,0.105));
    all->Add(p);stable->Add(p);
    // Pre-existing association deliberately tests Clone's copied-list handling.
    Candidate *dummy=delphes.GetFactory()->NewCandidate();p->AddCandidate(dummy);
    delphes.ProcessTask();
    TObjArray *accepted=delphes.ImportArray("Downstream/tracks");
    require(accepted->GetEntriesFast()==1,"accepted track count");
    Candidate *track=static_cast<Candidate *>(accepted->At(0));
    require(track!=p,"must clone input");
    require(track->GetCandidates()->GetEntriesFast()==1 && track->GetCandidates()->At(0)==p,"original truth association");
    near(p->Position.Z(),2000,1e-12,"input position unchanged");near(p->Momentum.Px(),0,1e-12,"input momentum unchanged");
    require(p->GetCandidates()->At(0)==dummy,"input association unchanged");
    near(track->InitialPosition.Z(),2000,1e-12,"birth mm");
    near(track->ZFirstHit,7700,1e-3,"front hit mm");near(track->Position.Z(),9400,1e-9,"endpoint mm");
    near(track->Momentum.Px(),-4*0.299792458,1e-10,"field exit momentum");
    require(track->Position.T()>170 && track->L>7400,"ct/path accumulated in mm");
    require(std::isnan(track->D0)&&std::isnan(track->C),"undefined perigee fields");
    TBranch *branch=writer.GetTree()->GetBranch("Track");require(branch!=nullptr,"Track branch");
    TClonesArray *tracks=*reinterpret_cast<TClonesArray **>(branch->GetAddress());
    require(tracks && tracks->GetEntriesFast()==1,"TreeWriter track count");
    Track *entry=static_cast<Track *>(tracks->At(0));
    near(entry->Z,2000,1e-3,"TreeWriter birth");near(entry->ZOuter,9400,1e-3,"TreeWriter endpoint");
    near(entry->P,20,1e-5,"TreeWriter momentum");near(entry->ZFirstHit,7700,1e-3,"TreeWriter first hit");
    require(entry->TOuter>entry->T,"TreeWriter time conversion");
    // In-memory TRef may resolve to Candidate; on disk the same UID identifies
    // the Particle GenParticle branch. Do not cast an in-memory Candidate to it.
    TObject *truth=entry->Particle.GetObject();
    require(truth && truth->GetUniqueID()==p->GetUniqueID(),"TreeWriter truth TRef UID");
    writer.Fill();writer.Clear();delphes.Clear();
    // A zero-accepted-track entry must still be written, preserving the denominator.
    delphes.ProcessTask();
    require(accepted->GetEntriesFast()==0 && tracks->GetEntriesFast()==0,"no stale tracks after clear");
    writer.Fill();require(writer.GetTree()->GetEntries()==2,"zero-track entry retained");
    delphes.FinishTask();
    std::cout<<"PASS: synthetic Candidate module/TreeWriter interface checks and zero-track entry\n";
  }catch(const std::exception &e) { std::cerr<<"FAIL: "<<e.what()<<'\n';return 1; }
}
