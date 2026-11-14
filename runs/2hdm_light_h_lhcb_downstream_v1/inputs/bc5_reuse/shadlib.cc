// Pythia 8.316 hadronisation / decay library for the BC5 scalar decay channels.
//
// For each channel and each m_S on a log grid the S decay products are injected at rest with
// the correct colour flow and hadronised + decayed by Pythia; the stable charged particles are
// written in the S rest frame.  This replaces the analytic jet / leading-track models.
//   usage: shadlib <channel> <mMin> <mMax> <nMass> <nSample> <seed> <out> [mass_file]
//
// v3.4 (2026-10-24): the mass grid must be EXACTLY the scan grid, because the consumer
// (run_scan.py) keys every point by round(m_S, 6) and no longer borrows a neighbouring node.
// If the optional mass_file is given, its lines (one mass per line, GeV) are used instead of
// the log grid.
// The events are written in the S rest frame with the parton / tau pair along +-z; the consumer
// rotates every event by a random SO(3) matrix (the decay axis of a spin-0 particle must be
// isotropic with respect to the laboratory axes).
#include "Pythia8/Pythia.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
using namespace Pythia8;

int main(int argc, char** argv) {
  std::string ch = argv[1];
  double mMin = atof(argv[2]), mMax = atof(argv[3]);
  int nMass = atoi(argv[4]);
  long ns = atol(argv[5]);
  int seed = atoi(argv[6]);
  std::string out = argv[7];

  std::vector<double> masses;
  if (argc > 8 && std::string(argv[8]) != "-") {
    FILE* fm = fopen(argv[8], "r");
    if (!fm) { fprintf(stderr, "cannot open mass file %s\n", argv[8]); return 4; }
    char buf[256];
    while (fgets(buf, sizeof(buf), fm)) {
      if (buf[0] == '#' || buf[0] == '\n') continue;
      double v = atof(buf);
      if (v > 0.0) masses.push_back(v);
    }
    fclose(fm);
    nMass = (int)masses.size();
  } else {
    for (int im = 0; im < nMass; ++im) masses.push_back(mMin * std::pow(mMax / mMin, (double)im / (nMass - 1)));
  }

  Pythia pythia;
  pythia.readString("ProcessLevel:all = off");
  pythia.readString("PartonLevel:all = off");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("ParticleDecays:limitTau0 = off");
  pythia.readString("Beams:eCM = 13000.");
  pythia.readString("Check:epTolErr = 1e12");
  pythia.readString("Print:quiet = on");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  if (!pythia.init()) { fprintf(stderr, "init failed\n"); return 3; }

  FILE* f = fopen(out.c_str(), "w");
  fprintf(f, "# channel=%s mMin=%.6g mMax=%.6g nMass=%d nSample=%ld seed=%d\n", ch.c_str(), mMin, mMax, nMass, ns, seed);
  fprintf(f, "# mass sample pid charge E px py pz\n");

  long rows = 0;
  for (int im = 0; im < nMass; ++im) {
    double m0 = masses[im];
    if (ch == "tautau" && m0 < 2.0 * pythia.particleData.m0(15)) continue;
    for (long is = 0; is < ns; ++is) {
      Event& ev = pythia.event;
      ev.reset();
      double E = m0 / 2.0, p = E;
      if (ch == "gg") {
        int a = ev.append(21, 23, 0,0,0,0,0,0, 0.,0., p, E, 0.);
        int b = ev.append(21, 23, 0,0,0,0,0,0, 0.,0.,-p, E, 0.);
        ev[a].col(1); ev[a].acol(2); ev[b].col(2); ev[b].acol(1);
      } else if (ch == "ss" || ch == "cc" || ch == "bb") {
        int q = (ch == "ss") ? 3 : (ch == "cc" ? 4 : 5);
        if (m0 < 2.0 * pythia.particleData.m0(q)) continue;
        int a = ev.append(q, 23, 0,0,0,0,0,0, 0.,0., p, E, 0.);
        int b = ev.append(-q, 23, 0,0,0,0,0,0, 0.,0.,-p, E, 0.);
        ev[a].col(1); ev[b].acol(1);
      } else if (ch == "tautau") {
        double mt = pythia.particleData.m0(15);
        double pt = std::sqrt(std::max(m0*m0/4.0 - mt*mt, 0.0));
        ev.append(15,  23, 0,0,0,0,0,0, 0.,0., pt, m0/2.0, mt);
        ev.append(-15, 23, 0,0,0,0,0,0, 0.,0.,-pt, m0/2.0, mt);
      } else { fprintf(stderr, "bad channel\n"); return 2; }

      pythia.forceHadronLevel();
      for (int rep = 0; rep < 40; ++rep) if (!pythia.moreDecays()) break;

      for (int i = 0; i < ev.size(); ++i) {
        const Particle& q = ev[i];
        if (!q.isFinal() || q.charge() == 0) continue;
        fprintf(f, "%.8g %ld %d %d %.8g %.8g %.8g %.8g\n",
                m0, is, q.id(), (int)q.charge(), q.e(), q.px(), q.py(), q.pz());
        ++rows;
      }
    }
  }
  fclose(f);
  printf("{\"channel\": \"%s\", \"rows\": %ld}\n", ch.c_str(), rows);
  return 0;
}