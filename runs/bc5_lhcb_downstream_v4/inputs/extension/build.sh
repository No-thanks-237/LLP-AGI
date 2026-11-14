#!/usr/bin/env bash
# Run ONLY in an ephemeral container based on the pinned collider image.
set -euo pipefail
test -d /build
test -d /extension/modules
src=/build/Delphes
if [[ ! -d "$src" ]]; then
  cp -a /opt/MG5_aMC_v3_7_0/Delphes "$src"
fi
cp /extension/modules/DownstreamPropagator.h "$src/modules/"
cp /extension/modules/DownstreamPropagator.cc "$src/modules/"
cp /extension/modules/DownstreamTransport.h "$src/modules/"
cp /extension/modules/BirthKinematics.h "$src/modules/"
cp /extension/cards/delphes_card_downstream*.tcl "$src/cards/"
python3 - "$src/modules/ModulesLinkDef.h" <<'PY'
import pathlib, sys
p = pathlib.Path(sys.argv[1])
s = p.read_text()
if '#include "modules/DownstreamPropagator.h"' not in s:
    s = s.replace('#include "modules/ParticlePropagator.h"',
                  '#include "modules/ParticlePropagator.h"\n#include "modules/DownstreamPropagator.h"')
if '#pragma link C++ class DownstreamPropagator+;' not in s:
    anchor = '#pragma link C++ class ParticlePropagator+;'
    if anchor not in s:
        raise RuntimeError('Unexpected ModulesLinkDef format')
    s = s.replace(anchor, anchor + '\n#pragma link C++ class DownstreamPropagator+;')
p.write_text(s)
PY
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -I/extension/modules \
  /extension/tests/test_transport.cpp -o /build/test_transport
/build/test_transport
g++ -std=c++17 -O2 -Wall -Wextra -pedantic -I/extension/modules \
  /extension/tests/test_birth_kinematics.cpp -o /build/test_birth_kinematics
/build/test_birth_kinematics
cd "$src"
sh ./configure
# Discard copied base build products only inside the isolated source copy.
make clean
make -j "${BUILD_JOBS:-4}" libDelphes.so DelphesHepMC2 DelphesHepMC3
bash /extension/verify_adapter.sh
