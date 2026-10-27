$ErrorActionPreference = 'Stop'
& docker compose -f (Join-Path $PSScriptRoot 'compose.json') up -d --wait --wait-timeout 90
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
Write-Host 'Magnus: http://localhost:3011/jobs'
