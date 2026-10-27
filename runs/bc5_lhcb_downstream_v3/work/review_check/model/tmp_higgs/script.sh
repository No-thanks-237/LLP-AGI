#!/bin/bash
set -euo pipefail
source ./config.env
echo '=== environment ==='
hostname; date -u
g++ --version | head -1
pythia8-config --version
PFX=$(pythia8-config --prefix)
echo "pythia8 prefix: $PFX"
echo '=== compile ==='
g++ -O2 -std=c++17 -o parents parents.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
echo compiled
export LD_LIBRARY_PATH="$PFX/lib"
mkdir -p out
sum=out/summary.jsonl
: > "$sum"
run_one() {
  mode=$1; nev=$2; sd=$3
  echo "=== $mode : $nev events, seed $sd ==="
  ./parents "$mode" "$nev" "$sd" out | tee -a "$sum"
}
run_one minbias "$NEVENTS_MINBIAS" "$SEED"
run_one bbar "$NEVENTS_BBAR" "$((SEED+1))"
run_one cbar "$NEVENTS_CBAR" "$((SEED+2))"
run_one higgs "$NEVENTS_HIGGS" "$((SEED+3))"
echo '=== gzip ==='
for f in out/*.txt; do [ -e "$f" ] && gzip -f "$f"; done
ls -l out/
echo '=== done ==='
