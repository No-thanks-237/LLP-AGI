#!/bin/bash
set -euo pipefail
source ./config.env
echo "=== host=$(hostname) date=$(date -u +%FT%TZ) ==="
PFX=$(pythia8-config --prefix)
echo "pythia prefix: $PFX"
t0=$(date +%s)
g++ -O2 -std=c++17 -o shadlib shadlib.cc $(pythia8-config --cxxflags) $(pythia8-config --libs)
echo "compile_seconds=$(( $(date +%s) - t0 ))"
export LD_LIBRARY_PATH="$PFX/lib"
mkdir -p out
echo "=== per-channel mass lists (exact nodes of work/common/mass_list_*.txt) ==="
for ch in gg ss cc tautau; do
  case "$ch" in
    tautau) mf=mass_list_tau.txt ;;
    *)      mf=mass_list_$ch.txt ;;
  esac
  echo "$ch ($mf, $(wc -l < $mf) nodes): $(tr '\n' ' ' < $mf)"
done
i=0
for ch in gg ss cc tautau; do
  i=$((i+1))
  case "$ch" in
    tautau) mf=mass_list_tau.txt ;;
    *)      mf=mass_list_$ch.txt ;;
  esac
  echo "=== $ch seed $((910+i)) nSample $NSAMPLE ==="
  ./shadlib "$ch" 0.021 62.49 "$(wc -l < $mf)" "$NSAMPLE" $((910+i)) "out/hadron_lib_$ch.txt" "$mf" | tee -a out/summary.jsonl
done
ls -l out/
echo "=== gzip ==="
for f in out/hadron_lib_*.txt; do gzip -f "$f"; done
ls -l out/
echo "=== headers (mass nodes present per channel) ==="
for ch in gg ss cc tautau; do echo -n "$ch: "; zcat out/hadron_lib_$ch.txt.gz | sed -n '1p'; done
for ch in gg ss cc tautau; do echo -n "$ch distinct mass nodes: "; zcat out/hadron_lib_$ch.txt.gz | awk '!/^#/{print $1}' | sort -u | wc -l; done
sha256sum out/hadron_lib_*.txt.gz
echo "=== done ==="
