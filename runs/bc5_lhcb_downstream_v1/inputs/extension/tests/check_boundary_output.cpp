// Read serialized ROOT output. No tested transport code is linked as an oracle.
#include "classes/DelphesClasses.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
constexpr double cLight=299792458.0;
struct ParticleSpec {
  int index,pid,status,charge,m1,m2,d1,d2,accepted;
  double birth[4],initialMomentum[4],front[4],end[4],momentum[4],path;
};
struct EventSpec { int number,nParticles,nTracks;std::string name;double weight;std::vector<ParticleSpec> particles; };
void require(bool ok,const std::string &what) { if(!ok) throw std::runtime_error(what); }
void near(double a,double b,double absTol,double relTol,const std::string &what) {
  require(std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=absTol+relTol*std::abs(b),what+": actual="+std::to_string(a)+", expected="+std::to_string(b));
}
std::vector<EventSpec> readSpecs(const char *path) {
  std::ifstream in(path);require(bool(in),"cannot open expected TSV");
  std::vector<EventSpec> specs;std::string line;
  while(std::getline(in,line)) {
    std::istringstream row(line);char type;row>>type;
    if(type=='E') { EventSpec e;row>>e.number>>e.name>>e.weight>>e.nParticles>>e.nTracks;require(bool(row),"invalid E spec");specs.push_back(e); }
    else if(type=='P') {
      require(!specs.empty(),"P without E");ParticleSpec p;
      row>>p.index>>p.pid>>p.status>>p.charge>>p.m1>>p.m2>>p.d1>>p.d2;
      for(double &v:p.birth) row>>v;
      for(double &v:p.initialMomentum) row>>v;
      row>>p.accepted;
      for(double &v:p.front) row>>v;
      for(double &v:p.end) row>>v;
      for(double &v:p.momentum) row>>v;
      row>>p.path;require(bool(row),"invalid P spec");specs.back().particles.push_back(p);
    }else require(false,"unknown spec record");
  }
  for(const EventSpec &s:specs) require(int(s.particles.size())==s.nParticles,"incomplete event spec");
  return specs;
}
int main(int argc,char **argv) {
  if(argc!=3 && argc!=4) { std::cerr<<"usage: check_boundary_output output.root expected.tsv [report.json]\n";return 2; }
  try {
    const auto specs=readSpecs(argv[2]);TFile file(argv[1],"READ");require(!file.IsZombie(),"invalid ROOT file");
    TTree *tree=dynamic_cast<TTree *>(file.Get("Delphes"));require(tree,"missing Delphes tree");
    require(tree->GetEntries()==static_cast<Long64_t>(specs.size()),"event denominator/entry count");
    TClonesArray *particles=nullptr,*tracks=nullptr,*events=nullptr,*weights=nullptr;
    for(const char *name:{"Particle","Track","Event","Weight"}) require(tree->GetBranch(name),std::string("missing branch ")+name);
    tree->SetBranchAddress("Particle",&particles);tree->SetBranchAddress("Track",&tracks);
    tree->SetBranchAddress("Event",&events);tree->SetBranchAddress("Weight",&weights);
    int totalTracks=0,zeroTrackEvents=0,etaConditionedTracks=0,etaSentinelTracks=0;
    double sumWeights=0,expectedSumWeights=0,zeroTrackWeightSum=0,expectedZeroTrackWeightSum=0;
    for(size_t entryIndex=0;entryIndex<specs.size();++entryIndex) {
      require(tree->GetEntry(entryIndex)>0,"entry read failed");const EventSpec &s=specs[entryIndex];
      const std::string context="event "+std::to_string(s.number)+" "+s.name+" ";
      require(events->GetEntriesFast()==1,context+"Event count");
      HepMCEvent *event=static_cast<HepMCEvent *>(events->At(0));require(event->Number==s.number,context+"event number");
      near(event->Weight,s.weight,1e-7,2e-7,context+"Event.Weight");
      require(weights->GetEntriesFast()==2,context+"all event weights retained");
      near(static_cast<Weight *>(weights->At(0))->Weight,s.weight,1e-7,2e-7,context+"Weight[0]");
      near(static_cast<Weight *>(weights->At(1))->Weight,2*s.weight,1e-7,2e-7,context+"Weight[1]");
      sumWeights+=event->Weight;expectedSumWeights+=s.weight;
      require(particles->GetEntriesFast()==s.nParticles,context+"Particle count");
      require(tracks->GetEntriesFast()==s.nTracks,context+"Track count");
      if(!s.nTracks) {
        ++zeroTrackEvents;
        zeroTrackWeightSum+=event->Weight;
        expectedZeroTrackWeightSum+=s.weight;
      }
      totalTracks+=s.nTracks;
      std::map<const TObject *,int> indices;
      for(int i=0;i<s.nParticles;++i) {
        const ParticleSpec &p=s.particles[i];GenParticle *truth=static_cast<GenParticle *>(particles->At(i));indices[truth]=i;
        require(p.index==i,context+"expected index");
        require(truth->PID==p.pid && truth->Status==p.status && truth->Charge==p.charge,context+"truth PID/status/charge");
        require(truth->M1==p.m1 && truth->M2==p.m2 && truth->D1==p.d1 && truth->D2==p.d2,context+"truth parent/daughter chain");
        const double position[]={truth->X,truth->Y,truth->Z};
        for(int k=0;k<3;++k) near(position[k],p.birth[k]*1000,5e-4,0,context+"unmodified truth position");
        near(truth->T,p.birth[3]/cLight,1e-14,0,context+"truth ct-to-time conversion");
        const double momentum[]={truth->Px,truth->Py,truth->Pz,truth->E};
        for(int k=0;k<4;++k) near(momentum[k],p.initialMomentum[k],3e-6,2e-7,context+"unmodified truth momentum/units");
      }
      std::set<int> found;
      for(int i=0;i<s.nTracks;++i) {
        Track *track=static_cast<Track *>(tracks->At(i));TObject *ref=track->Particle.GetObject();
        require(ref && dynamic_cast<GenParticle *>(ref),context+"serialized TRef resolves GenParticle");
        const auto index=indices.find(ref);require(index!=indices.end(),context+"TRef belongs to current Particle branch");
        const ParticleSpec &p=s.particles[index->second];require(p.accepted && found.insert(p.index).second,context+"correct unique accepted truth");
        require(track->PID==p.pid && track->Charge==p.charge,context+"track identity");
        const double birth[]={track->X,track->Y,track->Z};
        const double first[]={track->XFirstHit,track->YFirstHit,track->ZFirstHit};
        const double end[]={track->XOuter,track->YOuter,track->ZOuter};
        for(int k=0;k<3;++k) {
          near(birth[k],p.birth[k]*1000,5e-4,0,context+"track birth");
          near(first[k],p.front[k]*1000,5e-4,0,context+"first hit");
          near(end[k],p.end[k]*1000,5e-4,0,context+"endpoint RK4 oracle");
        }
        near(track->T,p.birth[3]/cLight,1e-14,0,context+"birth time");
        near(track->TOuter,p.end[3]/cLight,1e-14,0,context+"endpoint time");
        near(track->L,p.path*1000,8e-4,0,context+"path RK4 oracle");
        // The standard Track::P4 uses PT/eta and cannot represent the exact
        // PT=0 axis (TreeWriter writes eta=999.9). Reconstruct from stored P/PT
        // and the forward sign instead, so exact-axis boundary probes are valid.
        const double pz=std::sqrt(std::max(0.0,double(track->P)*track->P-double(track->PT)*track->PT));
        const double actual[]={track->PT*std::cos(track->Phi),track->PT*std::sin(track->Phi),pz,
                               std::hypot(double(track->P),double(track->Mass))};
        for(int k=0;k<4;++k) near(actual[k],p.momentum[k],1e-5,2e-6,context+"endpoint momentum RK4 oracle");
        const double expectedPT=std::hypot(p.momentum[0],p.momentum[1]);
        if(expectedPT>0) {
          const double expectedP=std::hypot(expectedPT,p.momentum[2]);
          const double sinSquared=(expectedPT/expectedP)*(expectedPT/expectedP);
          const double expectedEta=std::asinh(p.momentum[2]/expectedPT);
          const double baseTolerance=2e-6+2e-7*std::abs(expectedEta);
          // ROOT obtains eta from cos(theta). d eta / d cos(theta) =
          // 1/sin^2(theta), so double rounding is amplified close to the axis.
          // Four ulps conservatively cover the direction-cosine arithmetic.
          const double cosineRoundoff=4*std::numeric_limits<double>::epsilon();
          const double conditionedTolerance=cosineRoundoff/sinSquared;
          if(conditionedTolerance>baseTolerance) ++etaConditionedTracks;
          if(sinSquared<=cosineRoundoff && track->Eta==999.9f) {
            ++etaSentinelTracks;
          }else near(track->Eta,expectedEta,2e-6+conditionedTolerance,2e-7,context+"endpoint eta (conditioned)");
        }else {
          require(track->Eta==999.9f,context+"axis eta sentinel");
          ++etaSentinelTracks;
        }
        require(std::isnan(track->D0)&&std::isnan(track->DZ)&&std::isnan(track->C),context+"undefined solenoidal fields");
      }
      for(const ParticleSpec &p:s.particles) require(bool(found.count(p.index))==bool(p.accepted),context+"every selected/rejected particle accounted for");
    }
    near(sumWeights,expectedSumWeights,1e-5,0,"weighted denominator retained including zero-track events");
    near(zeroTrackWeightSum,expectedZeroTrackWeightSum,1e-5,0,"zero-track weighted denominator retained");
    std::ostringstream report;report<<"{\"passed\":true,\"events\":"<<specs.size()<<",\"accepted_tracks\":"<<totalTracks
      <<",\"zero_track_events\":"<<zeroTrackEvents<<",\"sum_event_weights\":"<<sumWeights
      <<",\"zero_track_weight_sum\":"<<zeroTrackWeightSum
      <<",\"eta_conditioned_tracks\":"<<etaConditionedTracks<<",\"eta_sentinel_tracks\":"<<etaSentinelTracks
      <<",\"eta_roundoff_bound\":\"4_double_epsilon_over_sin_theta_squared\""
      <<",\"serialized_truth_refs\":true,\"oracle\":\"independent_z_RK4\"}";
    std::cout<<report.str()<<'\n';if(argc==4) { std::ofstream out(argv[3]);require(bool(out),"cannot write report");out<<report.str()<<'\n'; }
  }catch(const std::exception &e) {
    std::cerr<<"FAIL: "<<e.what()<<'\n';
    if(argc==4) { std::ofstream out(argv[3]);out<<"{\"passed\":false,\"see_stderr\":true}\n"; }
    return 1;
  }
}
