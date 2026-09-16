[CmdletBinding()]
param(
    [string]$QuartoPath = '',
    [switch]$KeepTyp
)

$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
$slidesPath = Join-Path $projectRoot 'presentation/slides.qmd'

if (-not (Test-Path -LiteralPath $slidesPath -PathType Leaf)) {
    throw "Slide source was not found: $slidesPath"
}

if (-not $QuartoPath) {
    $quartoCommand = Get-Command quarto -ErrorAction SilentlyContinue
    if ($quartoCommand) {
        $QuartoPath = $quartoCommand.Source
    } else {
        $QuartoPath = Join-Path $env:LOCALAPPDATA 'Programs/Quarto/bin/quarto.exe'
    }
}
if (-not (Test-Path -LiteralPath $QuartoPath -PathType Leaf)) {
    throw 'Quarto 1.8 or later is required. Install Quarto or pass -QuartoPath.'
}

# Keep the pinned Typst package downloads inside the project. A first build
# needs network access; subsequent builds reuse this ignored cache.
$cachePath = Join-Path $projectRoot '.cache/typst-packages'
New-Item -ItemType Directory -Force -Path $cachePath | Out-Null
$previousPackageCache = $env:TYPST_PACKAGE_CACHE_PATH
$env:TYPST_PACKAGE_CACHE_PATH = $cachePath
# Honor the format and theme selected in slides.qmd.
$renderArguments = @('render', 'presentation/slides.qmd')
if ($KeepTyp) { $renderArguments += @('--metadata', 'keep-typ:true') }

Push-Location $projectRoot
try {
    & $QuartoPath @renderArguments
    if ($LASTEXITCODE -ne 0) {
        throw "Quarto rendering failed (exit code $LASTEXITCODE)."
    }
    $pdfPath = Join-Path $projectRoot 'presentation/slides.pdf'
    if (-not (Test-Path -LiteralPath $pdfPath -PathType Leaf)) {
        throw "Quarto completed but the expected PDF was not found: $pdfPath"
    }
    Write-Output "PDF: $pdfPath"
} finally {
    Pop-Location
    $env:TYPST_PACKAGE_CACHE_PATH = $previousPackageCache
}
