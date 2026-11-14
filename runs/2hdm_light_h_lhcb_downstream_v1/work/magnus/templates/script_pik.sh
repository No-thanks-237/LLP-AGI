#!/bin/bash
set -euo pipefail
source ./config.env
echo "=== host=$(hostname) date=$(date -u +%FT%TZ) ==="
PFX=$(pythia8-config --prefix)
echo "pythia prefix: $PFX"
t0=$(date +%s)
g++ -O2 -std=c++17 -o parents_pik parents_pik.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
echo "compile_seconds=$(( $(date +%s) - t0 ))"
export LD_LIBRARY_PATH="$PFX/lib"
mkdir -p out
echo "=== generating $NEVENTS events, seed $SEED ==="
t1=$(date +%s)
./parents_pik "$NEVENTS" "$SEED" out | tee -a out/summary.jsonl
echo "generation_seconds=$(( $(date +%s) - t1 ))"
ls -l out/
echo "=== gzip ==="
for f in out/*.txt; do [ -e "$f" ] && gzip -f "$f"; done
ls -l out/
echo "=== checks ==="
zcat out/minbias_pik.txt.gz | sed -n '1,3p'
echo -n "rows (incl. 2 header lines): "; zcat out/minbias_pik.txt.gz | wc -l
echo -n "columns of first data line: "; zcat out/minbias_pik.txt.gz | awk '!/^#/{print NF; exit}'
sha256sum out/minbias_pik.txt.gz
echo "=== done ==="
