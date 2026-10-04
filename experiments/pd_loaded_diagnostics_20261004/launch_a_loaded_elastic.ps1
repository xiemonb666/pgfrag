$ErrorActionPreference = 'Stop'
if (Test-Path variable:PSNativeCommandUseErrorActionPreference) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$Workspace = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
$Experiment = Join-Path $Workspace 'pgfrag\experiments\pd_loaded_diagnostics_20261004'
$ManifestPath = Join-Path $Experiment 'manifest.json'
$Manifest = Get-Content -LiteralPath $ManifestPath -Raw | ConvertFrom-Json
foreach ($Record in $Manifest.source_files) {
    $Path = Join-Path $Workspace ([string]$Record.path -replace '/', '\')
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { throw "Missing pinned file: $($Record.path)" }
    $Actual = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($Actual -ne ([string]$Record.sha256).ToLowerInvariant()) { throw "Source hash mismatch: $($Record.path)" }
}
$Julia = Join-Path $Workspace 'third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$JuliaProject = Join-Path $Workspace 'pgfrag\src\pgfrag\resources\julia'
$JuliaDepot = Join-Path $Workspace 'third_party\pgfrag-runtime\depot'
$RuntimeHash = (Get-FileHash -LiteralPath $Julia -Algorithm SHA256).Hash.ToLowerInvariant()
if ($RuntimeHash -ne ([string]$Manifest.runtime.sha256).ToLowerInvariant()) { throw 'Julia runtime hash mismatch' }

$Runner = Join-Path $Workspace 'pgfrag\experiments\pd_point_event_probe_20261004\run_pd_point_event_probe.jl'
$Points = Join-Path $Workspace 'pgfrag\experiments\pd_point_event_probe_20261004\input\coupon_uniform_voxel_5cube_125.points.csv'
$Metadata = Join-Path $Workspace 'pgfrag\experiments\pd_point_event_probe_20261004\input\coupon_uniform_voxel_5cube_125.metadata.toml'
$Common = Join-Path $Experiment 'common_loaded_case.toml'
$Case = Join-Path $Experiment 'case_a_loaded_elastic.toml'
$RunId = 'loaded_elastic_a8_run01'
$ReviewRoot = Join-Path $Workspace 'pgfrag\tmp\review_20261004\pd_loaded_diagnostics'
$OutDirectory = Join-Path $ReviewRoot 'OUT_A'
$SnapshotDirectory = Join-Path (Join-Path $ReviewRoot 'snapshots') $RunId
$ReceiptDirectory = Join-Path $ReviewRoot 'receipts'
$ReceiptPath = Join-Path $ReceiptDirectory ($RunId + '.launch_receipt.json')
$StdoutPath = Join-Path $ReceiptDirectory ($RunId + '.stdout.log')
$StderrPath = Join-Path $ReceiptDirectory ($RunId + '.stderr.log')
foreach ($Path in @($Julia, $JuliaProject, $JuliaDepot, $Runner, $Points, $Metadata, $Common, $Case)) {
    if (-not (Test-Path -LiteralPath $Path)) { throw "Required launch path missing: $Path" }
}
if (Test-Path -LiteralPath $SnapshotDirectory) { throw "Refusing to overwrite snapshot: $SnapshotDirectory" }
foreach ($Path in @($ReceiptPath, $StdoutPath, $StderrPath)) {
    if (Test-Path -LiteralPath $Path) { throw "Refusing to overwrite run artifact: $Path" }
}
New-Item -ItemType Directory -Force -Path $OutDirectory, (Split-Path -Parent $SnapshotDirectory), $ReceiptDirectory | Out-Null
if (@(Get-ChildItem -LiteralPath $OutDirectory -Force).Count -ne 0) { throw 'OUT_A must be empty' }

$PreviousThreads = $env:JULIA_NUM_THREADS
$PreviousDepot = $env:JULIA_DEPOT_PATH
$PreviousPrecompileTasks = $env:JULIA_NUM_PRECOMPILE_TASKS
$Process = $null
$StartedUtc = [DateTime]::UtcNow.ToString('o')
$TimedOut = $false
$WaitCompleted = $false
$LaunchError = $null
$ChildExit = $null
try {
    $env:JULIA_NUM_THREADS = '1'
    $env:JULIA_DEPOT_PATH = $JuliaDepot
    $env:JULIA_NUM_PRECOMPILE_TASKS = '1'
    $Arguments = @("--project=$JuliaProject", '--startup-file=no', '--history-file=no',
                   $Runner, $SnapshotDirectory, $Points, $Metadata, $Common, $Case)
    $Process = Start-Process -FilePath $Julia -ArgumentList ($Arguments -join ' ') `
        -WorkingDirectory $OutDirectory -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $StdoutPath -RedirectStandardError $StderrPath
    $WaitCompleted = $Process.WaitForExit([int]$Manifest.runtime.timeout_seconds * 1000)
    if (-not $WaitCompleted) {
        $TimedOut = $true
        try { $Process.Refresh(); if (-not $Process.HasExited) { $Process.Kill() } }
        catch { $LaunchError = "Timeout termination error: $($_.Exception.Message)" }
        $null = $Process.WaitForExit(5000)
    }
    $Process.Refresh()
    if ($Process.HasExited) { $ChildExit = [int]$Process.ExitCode }
} catch {
    $LaunchError = $_.Exception.Message
} finally {
    if ($null -eq $PreviousThreads) { Remove-Item Env:JULIA_NUM_THREADS -ErrorAction SilentlyContinue } else { $env:JULIA_NUM_THREADS = $PreviousThreads }
    if ($null -eq $PreviousDepot) { Remove-Item Env:JULIA_DEPOT_PATH -ErrorAction SilentlyContinue } else { $env:JULIA_DEPOT_PATH = $PreviousDepot }
    if ($null -eq $PreviousPrecompileTasks) { Remove-Item Env:JULIA_NUM_PRECOMPILE_TASKS -ErrorAction SilentlyContinue } else { $env:JULIA_NUM_PRECOMPILE_TASKS = $PreviousPrecompileTasks }
}
$FinishedUtc = [DateTime]::UtcNow.ToString('o')
$Status = if ($TimedOut) { 'TIMED_OUT' } elseif ($null -ne $LaunchError) { 'LAUNCH_ERROR' } elseif ($ChildExit -eq 0) { 'COMPLETED' } else { 'FAILED_OR_EXIT_UNKNOWN' }
$Receipt = [ordered]@{
    schema = 'pgfrag.pd.loaded_diagnostic.launch_receipt.v1'
    run_id = $RunId
    status = $Status
    wrapper_exit_code = if ($TimedOut) { 124 } elseif ($null -ne $LaunchError) { 125 } elseif ($null -ne $ChildExit) { $ChildExit } else { 125 }
    child_exit_code = $ChildExit
    child_exit_observed = ($null -ne $ChildExit)
    timed_out = $TimedOut
    wait_completed = $WaitCompleted
    launch_error = $LaunchError
    timeout_seconds = [int]$Manifest.runtime.timeout_seconds
    julia_threads = 1
    runtime_sha256 = $RuntimeHash
    started_utc = $StartedUtc
    finished_utc = $FinishedUtc
    argv = @('--project=' + $JuliaProject, '--startup-file=no', '--history-file=no', $Runner, $SnapshotDirectory, $Points, $Metadata, $Common, $Case)
    snapshot_directory = $SnapshotDirectory
    stdout_log = $StdoutPath
    stderr_log = $StderrPath
    native_energy_ledger = 'UNKNOWN_NOT_A_GATE'
    physical_qualification = 'NONE_UNVALIDATED_PROTOTYPE'
}
$Receipt | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $ReceiptPath -Encoding utf8
if ($TimedOut) { exit 124 }
if ($null -ne $LaunchError) { exit 125 }
if ($null -eq $ChildExit) { exit 125 }
exit $ChildExit
