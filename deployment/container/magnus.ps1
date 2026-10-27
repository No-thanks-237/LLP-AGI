# Project-local CLI; use the SDK inside the new service container.
$ErrorActionPreference = 'Stop'
$taskRoot = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..\..')).TrimEnd('\')
$mapped = @()
foreach ($value in $args) {
    $text = [string]$value
    if ($text -match '^[A-Za-z]:[\\/]') {
        $full = [IO.Path]::GetFullPath($text)
        if (-not $full.StartsWith($taskRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
            throw 'External files must first be copied into this run inputs directory.'
        }
        $text = '/workspace/' + $full.Substring($taskRoot.Length + 1).Replace('\','/')
    }
    elseif ($text -match '^(runs|temp|configs|deployment)[\\/]') {
        $text = $text.Replace('\','/')
    }
    $mapped += $text
}
$dockerArgs = @('exec','-e','PYTHONDONTWRITEBYTECODE=1','-w','/workspace','llp-magnus','magnus') + $mapped
& docker @dockerArgs
exit $LASTEXITCODE
