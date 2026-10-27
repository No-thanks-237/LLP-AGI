#!/bin/bash
set -euo pipefail
source ./config.env
PFX=$(pythia8-config --prefix)
g++ -O2 -std=c++17 -o parents parents.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
export LD_LIBRARY_PATH="$PFX/lib"
mkdir -p out
./parents minbias "$NEVENTS_MINBIAS" "$SEED" out | tee -a out/summary.jsonl
for f in out/*.txt; do [ -e "$f" ] && gzip -f "$f"; done
ls -l out/
