param(
    [ValidatePattern('^[A-Za-z0-9_-]+$')]
    [string]$AttemptName = 'attempt_03_final'
)

$ErrorActionPreference = 'Stop'

$root = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
$baseOutDir = Join-Path $root 'relicweave-synlab\tmp\review_20261003\implementation\pg\contact_spike'
$outDir = Join-Path $baseOutDir $AttemptName
$julia = Join-Path $root 'third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$depot = Join-Path $root 'third_party\pgfrag-runtime\depot'
$project = Join-Path $root 'pgfrag\src\pgfrag\resources\julia'
$runner = Join-Path $PSScriptRoot 'run_contact_spike.jl'
$plotter = Join-Path $PSScriptRoot 'render_contact_spike.py'

if (Test-Path -LiteralPath $outDir) {
    throw "Attempt output already exists; choose a new attempt name to preserve it: $outDir"
}
New-Item -ItemType Directory -Path $outDir | Out-Null
$env:JULIA_DEPOT_PATH = $depot
$env:JULIA_NUM_THREADS = '1'

Write-Output "CONTACT_SPIKE_COMMAND julia=$julia project=$project threads=1 output=$outDir"
& $julia '--compiled-modules=no' '--startup-file=no' '--history-file=no' "--project=$project" '-t' '1' $runner '--outdir' $outDir 2>&1 | Tee-Object -FilePath (Join-Path $outDir 'run.log')
$juliaExit = $LASTEXITCODE

$python = 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$reportPath = Join-Path $outDir 'report.json'
if ((Test-Path -LiteralPath $reportPath) -and
    ((Get-Content -LiteralPath $reportPath -Raw | ConvertFrom-Json).status -eq 'PREFLIGHT_FAILED')) {
    Write-Output 'PLOT_SKIPPED no fixed case was started because the native one-step preflight failed.'
    'PLOT_SKIPPED no case rows are expected after PREFLIGHT_FAILED.' | Set-Content -LiteralPath (Join-Path $outDir 'plot.log') -Encoding utf8
    $plotExit = 0
}
else {
    & $python $plotter '--outdir' $outDir 2>&1 | Tee-Object -FilePath (Join-Path $outDir 'plot.log')
    $plotExit = $LASTEXITCODE
}

if ($juliaExit -ne 0) {
    Write-Output "JULIA_RUNNER_EXIT=$juliaExit"
    exit $juliaExit
}
if ($plotExit -ne 0) {
    Write-Output "PLOTTER_EXIT=$plotExit"
    exit $plotExit
}
Write-Output "CONTACT_SPIKE_ARTIFACTS_READY output=$outDir"
