#!/bin/bash
set -euo pipefail
cd /opt/magnus/back_end
python /opt/llp/service.py &
backend=$!
cd /opt/magnus/front_end
node server.mjs -p 3011 -H 0.0.0.0 &
frontend=$!
trap 'kill -TERM "$backend" "$frontend" 2>/dev/null || true; wait || true' TERM INT EXIT
wait -n "$backend" "$frontend"
