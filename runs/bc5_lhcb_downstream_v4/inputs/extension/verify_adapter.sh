#!/usr/bin/env bash
set -euo pipefail
src=/build/Delphes
cd "$src"
# Also build the shared library when resuming a completed reader-only build.
make -j "${BUILD_JOBS:-4}" libDelphes.so
g++ $(root-config --cflags) -I"$src" -I"$src/external" \
  /extension/tests/test_module.cpp -L"$src" -lDelphes \
  $(root-config --libs) -lEG -o /build/test_module
LD_LIBRARY_PATH="$src:${LD_LIBRARY_PATH:-}" /build/test_module \
  "$src/cards/delphes_card_downstream_geometry_only.tcl"
g++ $(root-config --cflags) -I"$src" -I"$src/external" \
  /extension/tests/test_birth_module.cpp -L"$src" -lDelphes \
  $(root-config --libs) -lEG -o /build/test_birth_module
LD_LIBRARY_PATH="$src:${LD_LIBRARY_PATH:-}" /build/test_birth_module \
  "$src/cards/delphes_card_downstream.tcl"
mkdir -p /build/package/bin /build/package/lib /build/package/cards /build/package/source/modules
cp DelphesHepMC2 DelphesHepMC3 /build/package/bin/
cp libDelphes.so /build/package/lib/
find . -maxdepth 1 -type f \( -name '*_rdict.pcm' -o -name '*.rootmap' \) \
  -exec cp '{}' /build/package/lib/ \;
cp cards/delphes_card_downstream*.tcl /build/package/cards/
cp /extension/modules/* /build/package/source/modules/
cp modules/ModulesLinkDef.h /build/package/source/modules/
cp LICENSE /build/package/
cp /extension/README.md /build/package/
tar -czf /build/downstream-delphes.tar.gz -C /build/package .
