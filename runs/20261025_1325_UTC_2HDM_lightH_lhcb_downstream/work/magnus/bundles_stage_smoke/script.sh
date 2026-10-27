#!/bin/bash
set -euo pipefail
echo "=== SMOKE TEST: compile + tiny runs of both generators; date=$(date -u +%FT%TZ) ==="
PFX=$(pythia8-config --prefix)
echo "pythia prefix: $PFX"
pythia8-config --version 2>/dev/null || true
t0=$(date +%s)
g++ -O2 -std=c++17 -o parents_pik parents_pik.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
g++ -O2 -std=c++17 -o shadlib shadlib.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
t1=$(date +%s)
echo "compile_seconds=$((t1-t0))"
export LD_LIBRARY_PATH="$PFX/lib"
mkdir -p out
t2=$(date +%s)
./parents_pik 20000 20261025 out | tee -a out/summary.jsonl
t3=$(date +%s)
echo "minbias_20000_seconds=$((t3-t2))"
echo "minbias_events_per_second=$(python3 -c "print(20000.0/$((t3-t2)) if $((t3-t2))>0 else 0)")"
head -2 mass_list_gg.txt > mass_smoke.txt
t4=$(date +%s)
./shadlib gg 0.021 62.49 2 50 999 out/hadron_smoke.txt mass_smoke.txt | tee -a out/summary.jsonl
t5=$(date +%s)
echo "shadlib_2nodes_x50_seconds=$((t5-t4))"
echo "=== smoke sample head ==="
head -3 out/minbias_pik.txt
wc -l out/minbias_pik.txt out/hadron_smoke.txt
echo "=== column check (expect 10 columns, first non-comment line) ==="
grep -v '^#' out/minbias_pik.txt | head -1 | awk '{print NF}'
head -3 out/hadron_smoke.txt
gzip -f out/*.txt
ls -l out/
echo "=== done ==="
