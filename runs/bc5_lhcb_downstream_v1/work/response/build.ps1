# Build downstream_kernel.dll and run the Python validation suite.
#
# Idempotent: every run deletes the previous artifact, rebuilds from source and
# rewrites the log.  The script fails loudly (non-zero exit, thrown error) when
# compilation fails or when the shared library is not produced.
#
# Usage:
#   pwsh -File build.ps1
#   pwsh -File build.ps1 -Compiler <g++.exe> -Python <python.exe>

[CmdletBinding()]
param(
    [string]$Compiler = 'D:/software/msys2/ucrt64/bin/g++.exe',
    [string]$Python   = 'D:/py/python.exe'
)

$ErrorActionPreference = 'Stop'

# Never let a non-zero native exit code become a terminating error before we
# have captured and logged it.
if (Get-Variable -Name PSNativeCommandUseErrorActionPreference -ErrorAction SilentlyContinue) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$ResponseDir = $PSScriptRoot
$WorkDir     = Split-Path -Parent $ResponseDir
$RunDir      = Split-Path -Parent $WorkDir
$LogDir      = Join-Path $RunDir 'logs'
$LogPath     = Join-Path $LogDir 'response_build.log'

New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
Set-Content -LiteralPath $LogPath -Value '' -Encoding utf8

function Emit {
    param([string]$Text)
    Add-Content -LiteralPath $LogPath -Value $Text -Encoding utf8
    Write-Host $Text
}

$IncludeDir = 'E:/AAAwork/LLP-AGI/detector/downstream/modules'
$Source     = Join-Path $ResponseDir 'downstream_kernel.cpp'
$Target     = Join-Path $ResponseDir 'downstream_kernel.dll'

$compileArgs = @(
    '-std=c++17'
    '-O2'
    '-Wall'
    '-Wextra'
    '-shared'
    '-static'
    '-static-libgcc'
    '-static-libstdc++'
    '-I'
    $IncludeDir
    $Source
    '-o'
    $Target
)

Emit '== downstream_kernel build =='
Emit "timestamp : $(Get-Date -Format o)"
Emit "compiler  : $Compiler"
Emit "source    : $Source"
Emit "output    : $Target"
Emit "include   : $IncludeDir"
Emit "command   : & '$Compiler' $($compileArgs -join ' ')"
Emit ''

if (-not (Test-Path -LiteralPath $Compiler)) {
    Emit "ERROR: compiler not found: $Compiler"
    throw "compiler not found: $Compiler"
}
if (-not (Test-Path -LiteralPath $Source)) {
    Emit "ERROR: source not found: $Source"
    throw "source not found: $Source"
}

# Idempotency: start from a clean artifact so a stale DLL can never be mistaken
# for a successful build.
if (Test-Path -LiteralPath $Target) {
    Remove-Item -LiteralPath $Target -Force
}

$compileOutput = & $Compiler @compileArgs 2>&1 | Out-String
$compileExit = $LASTEXITCODE

Emit '-- compiler stdout+stderr --'
Emit ($compileOutput.TrimEnd())
Emit "-- compiler exit code: $compileExit --"
Emit ''

if ($compileExit -ne 0) {
    Emit 'BUILD FAILED: compiler returned a non-zero exit code.'
    throw "downstream_kernel.dll build FAILED (exit $compileExit); see $LogPath"
}
if (-not (Test-Path -LiteralPath $Target)) {
    Emit 'BUILD FAILED: compiler returned 0 but the DLL is missing.'
    throw "compiler returned 0 but $Target is missing; see $LogPath"
}

$dllItem = Get-Item -LiteralPath $Target
Emit "BUILD OK: $($dllItem.FullName) ($($dllItem.Length) bytes)"
Emit ''

# ---- run the Python validation suite -------------------------------------
$TestScript = Join-Path $ResponseDir 'test_kernel.py'

if (-not (Test-Path -LiteralPath $Python)) {
    $fallback = Get-Command python -ErrorAction SilentlyContinue
    if ($null -eq $fallback) {
        Emit "ERROR: no Python interpreter found (tried $Python and 'python')."
        throw "no Python interpreter found"
    }
    Emit "WARNING: $Python not found; falling back to $($fallback.Source)"
    $Python = $fallback.Source
}

$env:PYTHONDONTWRITEBYTECODE = '1'
Emit "python    : $Python"
Emit "test      : $TestScript"
Emit ''
Emit '== test_kernel.py =='

$testOutput = & $Python -B $TestScript 2>&1 | Out-String
$testExit = $LASTEXITCODE

Emit ($testOutput.TrimEnd())
Emit ""
Emit "== test_kernel.py exit code: $testExit =="

exit $testExit
