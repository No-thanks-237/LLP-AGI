#!/bin/bash
set -euo pipefail
PFX=$(pythia8-config --prefix)
echo "prefix: $PFX"
g++ -O2 -std=c++17 -o shadlib shadlib.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
export LD_LIBRARY_PATH="$PFX/lib"
mkdir -p out
echo "=== mass list ($(wc -l < mass_list.txt) points) ==="
head -3 mass_list.txt; tail -2 mass_list.txt
i=0
for ch in gg ss cc bb tautau; do
  i=$((i+1))
  echo "=== $ch ==="
  ./shadlib $ch 0.021 62.49 52 1500 $((900+i)) out/hadron_lib_$ch.txt mass_list.txt
done
ls -l out/
echo "=== gzip ==="
for f in out/hadron_lib_*.txt; do gzip -f "$f"; done
ls -l out/
echo "=== done ==="
