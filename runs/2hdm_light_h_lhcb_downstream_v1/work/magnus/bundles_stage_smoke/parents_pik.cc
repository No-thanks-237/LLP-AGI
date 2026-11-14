// Parent sample for the 2HDM light-H LHCb downstream analysis: pi/K from inelastic pp.
//
// Adapted from inputs/bc5_reuse/parents.cc (mode minbias) with these deliberate changes:
//   * ONE output file with all of pi+-, K+-, K0/K0bar (PDG 211,-211,321,-321,311,-311);
//   * no eta pre-cut of any kind (the downstream selection is applied to the H, not the parent);
//   * 10-column format: event pdg px py pz E x y z t   (x,y,z in mm, t in mm/c as Pythia gives);
//   * per-generation summary with sigmaGen and the pi/K multiplicities (validation).
// Everything else (beams, seed handling, no decays: ParticleDecays:tau0Max = 0) is identical
// to the BC5 v3.5 parents.cc, so the fragmentation-level particle content is comparable.
//
// Usage: parents_pik <nEvents> <seed> <outdir>
#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <map>

using namespace Pythia8;

static const int PID[6] = {211, -211, 321, -321, 311, -311};
static inline bool keep(int id) {
  for (int k = 0; k < 6; ++k) if (PID[k] == id) return true;
  return false;
}

int main(int argc, char** argv) {
  if (argc < 4) { fprintf(stderr, "usage: %s <nEvents> <seed> <outdir>\n", argv[0]); return 2; }
  long nEvents = atol(argv[1]);
  int seed = atoi(argv[2]);
  std::string outdir = argv[3];

  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("SoftQCD:inelastic = on");
  pythia.readString("SoftQCD:elastic = off");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 0.0");   // keep fragmentation-level hadrons
  if (!pythia.init()) { fprintf(stderr, "Pythia init failed\n"); return 3; }

  FILE* f = fopen((outdir + "/minbias_pik.txt").c_str(), "w");
  if (!f) { fprintf(stderr, "cannot open output\n"); return 4; }
  fprintf(f, "# event pdg px py pz E x y z t\n");
  fprintf(f, "# Pythia 8.316 pp 13 TeV SoftQCD:inelastic=on status 91/92, no eta cut, seed %d\n", seed);

  long kept = 0, nEvWithMeson = 0, nEvWithMesonEta6 = 0, nRowsEta6 = 0;
  std::map<int, long> byPid;
  for (long i = 0; i < nEvents; ++i) {
    if (!pythia.next()) continue;
    long rowsThisEvent = 0; bool inEta6 = false;
    for (int j = 0; j < pythia.event.size(); ++j) {
      const Particle& p = pythia.event[j];
      if (p.status() != 91 && p.status() != 92) continue;
      int id = p.id();
      if (!keep(id)) continue;
      fprintf(f, "%ld %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", i, id,
              p.px(), p.py(), p.pz(), p.e(), p.xProd(), p.yProd(), p.zProd(), p.tProd());
      ++rowsThisEvent; ++kept; byPid[id]++;
      if (std::fabs(p.eta()) < 6.0) { ++nRowsEta6; inEta6 = true; }
    }
    if (rowsThisEvent > 0) ++nEvWithMeson;
    if (inEta6) ++nEvWithMesonEta6;
  }
  fclose(f);

  printf("{\"mode\": \"minbias_pik\", \"events\": %ld, \"rows\": %ld, "
         "\"events_with_pi_or_K\": %ld, \"events_with_pi_or_K_eta6\": %ld, "
         "\"rows_eta6\": %ld, \"sigmaGen_mb\": %.8g, \"sigmaErr_mb\": %.8g, "
         "\"pid_counts\": {", nEvents, kept, nEvWithMeson, nEvWithMesonEta6, nRowsEta6,
         pythia.info.sigmaGen(), pythia.info.sigmaErr());
  bool first = true;
  for (std::map<int, long>::const_iterator it = byPid.begin(); it != byPid.end(); ++it) {
    printf("%s\"%d\": %ld", first ? "" : ", ", it->first, it->second);
    first = false;
  }
  printf("}}\n");
  return 0;
}
