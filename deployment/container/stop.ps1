$ErrorActionPreference = 'Stop'
# Stops the control plane only; first confirm there are no Running jobs.
& docker compose -f (Join-Path $PSScriptRoot 'compose.json') stop
exit $LASTEXITCODE
