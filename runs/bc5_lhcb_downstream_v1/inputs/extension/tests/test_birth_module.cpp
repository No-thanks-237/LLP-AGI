// In-memory Candidate interface unit test; no generated/HepMC events.
#include "modules/Delphes.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "ExRootAnalysis/ExRootConfReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "TApplication.h"
#include "TMemFile.h"
#include "TObjArray.h"
#include "TROOT.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
void require(bool ok,const char *name) { if(!ok) throw std::runtime_error(name); }
int main(int argc,char **argv) {
  if(argc!=2) return 2;
  try {
    gROOT->SetBatch();int n=1;char *a[]={argv[0]};TApplication app("test_birth_module",&n,a);
    TMemFile memory("test_birth_module.root","RECREATE");ExRootTreeWriter writer(&memory,"Delphes");
    ExRootConfReader reader;reader.ReadFile(argv[1]);
    require(reader.GetDouble("::Downstream::MinBirthEnergy",-1)==5,"formal E threshold");
    require(reader.GetDouble("::Downstream::MinBirthPT",-1)==.5,"formal birth PT threshold");
    Delphes delphes("Delphes");delphes.SetConfReader(&reader);delphes.SetTreeWriter(&writer);
    TObjArray *all=delphes.ExportArray("allParticles"),*stable=delphes.ExportArray("stableParticles");
    delphes.ExportArray("partons");delphes.InitTask();
    Candidate *wanted=nullptr;
    for(int i=0;i<4;++i) {
      Candidate *p=delphes.GetFactory()->NewCandidate();p->PID=-13;p->Status=1;p->Charge=1;
      p->Position.SetXYZT(500,200,2000,170);
      const double px=i==0 ? 0 : i==2 ? .5 : 1;
      const double pz=i==3 ? 4 : 20;
      const double energy=i==3 ? 5 : std::sqrt(px*px+pz*pz+.105*.105);
      p->Momentum.SetPxPyPzE(px,0,pz,energy);all->Add(p);stable->Add(p);
      if(i==1) wanted=p;
    }
    delphes.ProcessTask();TObjArray *tracks=delphes.ImportArray("Downstream/tracks");
    require(tracks->GetEntriesFast()==1,"zero PT / exact PT / exact E rejected");
    Candidate *track=static_cast<Candidate *>(tracks->At(0));
    require(track->GetCandidates()->At(0)==wanted,"birth-passing original truth selected");
    require(track->Momentum.Pt()<.5 && wanted->Momentum.Pt()==1,"birth cut survives magnetic cancellation");
    require(all->GetEntriesFast()==4,"cut does not discard truth particles");
    writer.Fill();require(writer.GetTree()!=nullptr,"output tree exists");
    writer.Clear();delphes.Clear();delphes.ProcessTask();
    require(tracks->GetEntriesFast()==0,"zero track entry after clear");writer.Fill();
    delphes.FinishTask();std::cout<<"PASS: formal-card birth cuts applied before transport; truth preserved\n";
  }catch(const std::exception &e) { std::cerr<<"FAIL: "<<e.what()<<'\n';return 1; }
}
