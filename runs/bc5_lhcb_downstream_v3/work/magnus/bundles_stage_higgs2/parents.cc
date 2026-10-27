// Parent samples for the BC5 LHCb downstream analysis, Pythia8 8.316.
//   mode minbias : SoftQCD:inelastic = on  -> B, D and K mesons (fragmentation level)
//   mode bbar    : gg/qqbar -> b bbar only -> B mesons (high statistics)
//   mode cbar    : gg/qqbar -> c cbar only -> D mesons (high statistics)
//   mode higgs   : all Pythia8 HiggsSM production processes -> on-shell h four-momenta
// Usage: parents <mode> <nEvents> <seed> <outdir>
#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <set>

using namespace Pythia8;

static const std::set<int> HEAVY = {521,-521,511,-511,531,-531,541,-541,411,-411,421,-421,431,-431};
static const std::set<int> BONLY = {521,-521,511,-511,531,-531,541,-541};
static const std::set<int> DONLY = {411,-411,421,-421,431,-431};
static const std::set<int> LIGHT = {321,-321,311,-311};

int main(int argc, char** argv) {
  if (argc < 5) { fprintf(stderr, "usage: %s <mode> <nEvents> <seed> <outdir>\n", argv[0]); return 2; }
  std::string mode = argv[1];
  long nEvents = atol(argv[2]);
  int seed = atoi(argv[3]);
  std::string outdir = argv[4];

  Pythia pythia;
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Print:quiet = on");
  pythia.readString("Next:numberShowEvent = 0");
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 0.0");   // keep fragmentation-level hadrons

  bool isHiggs = (mode == "higgs");
  if (mode == "minbias") {
    pythia.readString("SoftQCD:inelastic = on");
    pythia.readString("SoftQCD:elastic = off");
  } else if (mode == "bbar" || mode == "cbar") {
    pythia.readString("SoftQCD:inelastic = off");      // clean pp -> Q Qbar sample
    pythia.readString("HardQCD:all = off");
    char q = (mode == "bbar") ? 'b' : 'c';
    pythia.readString(std::string("HardQCD:gg2") + q + q + "bar = on");
    pythia.readString(std::string("HardQCD:qqbar2") + q + q + "bar = on");
    pythia.readString("PhaseSpace:pTHatMin = 0.");
    pythia.readString("PhaseSpace:pTHatMinDiverge = 0.5");
    pythia.readString("MultipartonInteractions:ecmPow = 0.20");   // keep UE/MPI active
  } else if (isHiggs) {
    pythia.readString("HiggsSM:all = on");
    pythia.readString("25:mayDecay = off");                      // keep the Higgs on shell
    pythia.readString("25:m0 = 125.25");
    pythia.readString("25:mWidth = 0.0041");
    pythia.readString("25:doForceWidth = on");
  } else {
    fprintf(stderr, "unknown mode\n"); return 2;
  }

  if (!pythia.init()) { fprintf(stderr, "Pythia init failed\n"); return 3; }

  FILE* fh = nullptr; FILE* fl = nullptr; FILE* fg = nullptr;
  if (mode == "minbias") {
    fh = fopen((outdir + "/" + mode + "_heavy.txt").c_str(), "w");
    fl = fopen((outdir + "/" + mode + "_light.txt").c_str(), "w");
    fprintf(fh, "# event pdg px py pz E x y z t\n");
    fprintf(fl, "# event pdg px py pz E x y z t\n");
  } else if (mode == "bbar" || mode == "cbar") {
    fh = fopen((outdir + "/" + mode + "_heavy.txt").c_str(), "w");
    fprintf(fh, "# event pdg px py pz E x y z t\n");
  } else {
    fg = fopen((outdir + "/higgs.txt").c_str(), "w");
    // v3.4: the 11th column is the Pythia STATUS of the record, so that the consumer can pick
    // the pre-decay copy (|status| = 62) instead of guessing from the record order.  The 2nd
    // column stays pythia.info.code() (the hard-process code).
    fprintf(fg, "# event code px py pz E x y z t status\n");
  }

  const double ETA_MIN = 1.5, ETA_MAX = 5.5;
  long kept = 0;
  for (long i = 0; i < nEvents; ++i) {
    if (!pythia.next()) continue;
    if (mode == "higgs") {
      for (int j = 0; j < pythia.event.size(); ++j) {
        const Particle& p = pythia.event[j];
        if (p.id() != 25) continue;
        if (std::abs(p.eta()) > 6.0) continue;
        fprintf(fg, "%ld %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %d\n",
                i, pythia.info.code(), p.px(), p.py(), p.pz(), p.e(),
                p.xProd(), p.yProd(), p.zProd(), p.tProd(), p.status());
        ++kept;
      }
      continue;
    }
    for (int j = 0; j < pythia.event.size(); ++j) {
      const Particle& p = pythia.event[j];
      if (p.status() != 91 && p.status() != 92) continue;
      int id = p.id();
      double eta = p.eta();
      if (!(eta > ETA_MIN && eta < ETA_MAX)) continue;
      if (mode == "bbar") {
        if (!BONLY.count(id)) continue;
        fprintf(fh, "%ld %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", i, id,
                p.px(), p.py(), p.pz(), p.e(), p.xProd(), p.yProd(), p.zProd(), p.tProd());
        ++kept;
      } else if (mode == "cbar") {
        if (!DONLY.count(id)) continue;
        fprintf(fh, "%ld %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", i, id,
                p.px(), p.py(), p.pz(), p.e(), p.xProd(), p.yProd(), p.zProd(), p.tProd());
        ++kept;
      } else {
        FILE* tgt = LIGHT.count(id) ? fl : (HEAVY.count(id) ? fh : nullptr);
        if (!tgt) continue;
        fprintf(tgt, "%ld %d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", i, id,
                p.px(), p.py(), p.pz(), p.e(), p.xProd(), p.yProd(), p.zProd(), p.tProd());
        ++kept;
      }
    }
  }
  if (fh) fclose(fh); if (fl) fclose(fl); if (fg) fclose(fg);
  printf("{\"mode\": \"%s\", \"events\": %ld, \"kept\": %ld, \"sigmaGen_mb\": %.8g, \"sigmaErr_mb\": %.8g}\n",
         mode.c_str(), nEvents, kept, pythia.info.sigmaGen(), pythia.info.sigmaErr());
  return 0;
}