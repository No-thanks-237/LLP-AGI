#!/usr/bin/env bash
# build.sh -- build the batch executable from the PATCHED author source.
# Author source: <run>/inputs/decay_source (Light_scalar_decay, Combined_Code.zip)
# Patches (see PATCHES.md): (1) M_PI definition in include/Parameters.h;
#                           (2) main() replaced by src_patched/drv.cpp batch driver.
# No physics edits.  Toolchain: MSYS2 g++ 15.2.0 (ucrt64), -std=c++11.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$HERE/src_patched"
OUT="$HERE/build"
mkdir -p "$OUT"
CXX=${CXX:-g++}
# -DM_PI is kept on the command line too (belt and braces; Parameters.h guard makes it idempotent)
"$CXX" -std=c++11 -O2 -DM_PI=3.14159265358979323846 \
      -I"$SRC/include" \
      -o "$OUT/light_scalar_decay.exe" \
      "$SRC/drv.cpp" "$SRC"/src/*.cpp "$SRC"/src/data_base/*.cpp
echo "built: $OUT/light_scalar_decay.exe"
