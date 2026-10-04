$ErrorActionPreference = 'Stop'
if (Test-Path variable:PSNativeCommandUseErrorActionPreference) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$Workspace = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
$Experiment = Join-Path $Workspace 'pgfrag\experiments\pd_loaded_diagnostics_20261004'
$ManifestPath = Join-Path $Experiment 'manifest.json'
$Manifest = Get-Content -LiteralPath $ManifestPath -Raw | ConvertFrom-Json
if ([string]$Manifest.B.status -ne 'FROZEN_PENDING_SINGLE_RUN') {
    throw "B is not frozen for its single authorized run: status=$($Manifest.B.status)"
}

$LauncherPath = $PSCommandPath
$LauncherHash = (Get-FileHash -LiteralPath $LauncherPath -Algorithm SHA256).Hash.ToLowerInvariant()
if ($LauncherHash -ne ([string]$Manifest.B.launcher_sha256).ToLowerInvariant()) {
    throw 'B launcher hash does not match the frozen manifest pin'
}

$PrePins = @()
foreach ($Record in $Manifest.source_files) {
    $Path = Join-Path $Workspace ([string]$Record.path -replace '/', '\')
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) { throw "Missing pinned file: $($Record.path)" }
    $Item = Get-Item -LiteralPath $Path
    $Actual = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($Actual -ne ([string]$Record.sha256).ToLowerInvariant() -or $Item.Length -ne [long]$Record.bytes) {
        throw "Pre-run source pin mismatch: $($Record.path)"
    }
    $PrePins += [ordered]@{ path=[string]$Record.path; bytes=[long]$Item.Length; sha256=$Actual }
}
if ($PrePins.Count -ne 26) { throw "Expected 26 frozen source pins, found $($PrePins.Count)" }

$Julia = Join-Path $Workspace 'third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$JuliaProject = Join-Path $Workspace 'pgfrag\src\pgfrag\resources\julia'
$JuliaDepot = Join-Path $Workspace 'third_party\pgfrag-runtime\depot'
$RuntimeHash = (Get-FileHash -LiteralPath $Julia -Algorithm SHA256).Hash.ToLowerInvariant()
if ($RuntimeHash -ne ([string]$Manifest.runtime.sha256).ToLowerInvariant()) { throw 'Julia runtime hash mismatch' }

$Driver = Join-Path $Experiment 'run_affine_critical_stretch_probe.jl'
$Points = Join-Path $Workspace 'pgfrag\experiments\pd_point_event_probe_20261004\input\coupon_uniform_voxel_5cube_125.points.csv'
$Metadata = Join-Path $Workspace 'pgfrag\experiments\pd_point_event_probe_20261004\input\coupon_uniform_voxel_5cube_125.metadata.toml'
$Common = Join-Path $Experiment 'common_loaded_case.toml'
$Case = Join-Path $Experiment 'case_b_affine_threshold.toml'
$DriverHash = (Get-FileHash -LiteralPath $Driver -Algorithm SHA256).Hash.ToLowerInvariant()
$CaseHash = (Get-FileHash -LiteralPath $Case -Algorithm SHA256).Hash.ToLowerInvariant()
if ($DriverHash -ne ([string]$Manifest.B.driver_sha256).ToLowerInvariant()) { throw 'B driver hash does not match frozen manifest pin' }
$PinnedCase = @($Manifest.source_files | Where-Object { $_.path -eq 'pgfrag/experiments/pd_loaded_diagnostics_20261004/case_b_affine_threshold.toml' })
if ($PinnedCase.Count -ne 1 -or $CaseHash -ne ([string]$PinnedCase[0].sha256).ToLowerInvariant()) {
    throw 'B case hash does not match frozen manifest source pin'
}

$RunId = 'affine_threshold_b_run01'
$ReviewRoot = Join-Path $Workspace 'pgfrag\tmp\review_20261004\pd_loaded_diagnostics'
$OutDirectory = Join-Path $ReviewRoot 'OUT_B'
$SnapshotDirectory = Join-Path (Join-Path $ReviewRoot 'snapshots') $RunId
$ReceiptDirectory = Join-Path $ReviewRoot 'receipts'
$ReceiptPath = Join-Path $ReceiptDirectory ($RunId + '.launch_receipt.json')
$StdoutPath = Join-Path $ReceiptDirectory ($RunId + '.stdout.log')
$StderrPath = Join-Path $ReceiptDirectory ($RunId + '.stderr.log')
$PreflightPath = Join-Path $ReceiptDirectory ($RunId + '.preflight.json')
foreach ($Path in @($Julia, $JuliaProject, $JuliaDepot, $Driver, $Points, $Metadata, $Common, $Case)) {
    if (-not (Test-Path -LiteralPath $Path)) { throw "Required launch path missing: $Path" }
}
if (Test-Path -LiteralPath $SnapshotDirectory) { throw "Refusing to overwrite B snapshot: $SnapshotDirectory" }
foreach ($Path in @($ReceiptPath, $StdoutPath, $StderrPath, $PreflightPath)) {
    if (Test-Path -LiteralPath $Path) { throw "Refusing to overwrite B run artifact: $Path" }
}
New-Item -ItemType Directory -Force -Path $OutDirectory, (Split-Path -Parent $SnapshotDirectory), $ReceiptDirectory | Out-Null
if (@(Get-ChildItem -LiteralPath $OutDirectory -Force).Count -ne 0) { throw 'OUT_B must be empty' }

$Arguments = @('--project=' + $JuliaProject, '--startup-file=no', '--history-file=no',
    $Driver, $Points, $Metadata, $Common, $Case, $SnapshotDirectory)
$ManifestHash = (Get-FileHash -LiteralPath $ManifestPath -Algorithm SHA256).Hash.ToLowerInvariant()
$Preflight = [ordered]@{
    schema = 'pgfrag.pd.loaded_diagnostic.B_preflight.v1'
    status = 'PREFLIGHT_PASS_SINGLE_RUN_ARMED'
    run_id = $RunId
    launcher_path = $LauncherPath
    launcher_sha256 = $LauncherHash
    manifest_path = $ManifestPath
    manifest_sha256 = $ManifestHash
    source_pin_count = $PrePins.Count
    source_pins = $PrePins
    runtime_sha256 = $RuntimeHash
    driver_sha256 = $DriverHash
    case_sha256 = $CaseHash
    julia_threads = 1
    timeout_seconds = [int]$Manifest.runtime.timeout_seconds
    argv = $Arguments
    output_directory = $OutDirectory
    snapshot_directory = $SnapshotDirectory
    states = @('rest: eta/epsilon_c=0', 'below: eta/epsilon_c=0.5', 'above: eta/epsilon_c=1.5')
    execution_contract = 'Three pristine handlers; native force/failure evaluation only; no time integration; no manual bond activity edits; one attempt only.'
    materials = 'Frozen UNVALIDATED_PROTOTYPE values; no calibration or physical qualification.'
}
$Utf8NoBom = [System.Text.UTF8Encoding]::new($false)
[System.IO.File]::WriteAllText($PreflightPath, ($Preflight | ConvertTo-Json -Depth 10) + [Environment]::NewLine, $Utf8NoBom)
$PreflightHash = (Get-FileHash -LiteralPath $PreflightPath -Algorithm SHA256).Hash.ToLowerInvariant()

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
$PostPins = @()
$PostPinError = $null
try {
    foreach ($Record in $Manifest.source_files) {
        $Path = Join-Path $Workspace ([string]$Record.path -replace '/', '\')
        $Item = Get-Item -LiteralPath $Path
        $Actual = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
        $PostPins += [ordered]@{ path=[string]$Record.path; bytes=[long]$Item.Length; sha256=$Actual; matches_expected=($Actual -eq ([string]$Record.sha256).ToLowerInvariant() -and $Item.Length -eq [long]$Record.bytes) }
    }
} catch { $PostPinError = $_.Exception.Message }
$PostPinsMatch = ($null -eq $PostPinError) -and (@($PostPins | Where-Object { -not $_.matches_expected }).Count -eq 0) -and ($PostPins.Count -eq 26)
$OutputArtifacts = @()
foreach ($Name in @('affine_native_bond_states.csv', 'affine_canonical_pair_states.csv', 'summary.toml')) {
    $Path = Join-Path $SnapshotDirectory $Name
    if (Test-Path -LiteralPath $Path -PathType Leaf) {
        $Item = Get-Item -LiteralPath $Path
        $OutputArtifacts += [ordered]@{ path=$Path; bytes=[long]$Item.Length; sha256=(Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant() }
    }
}
$StdoutInfo = if (Test-Path -LiteralPath $StdoutPath) { $Item = Get-Item -LiteralPath $StdoutPath; [ordered]@{bytes=[long]$Item.Length;sha256=(Get-FileHash -LiteralPath $StdoutPath -Algorithm SHA256).Hash.ToLowerInvariant()} } else { $null }
$StderrInfo = if (Test-Path -LiteralPath $StderrPath) { $Item = Get-Item -LiteralPath $StderrPath; [ordered]@{bytes=[long]$Item.Length;sha256=(Get-FileHash -LiteralPath $StderrPath -Algorithm SHA256).Hash.ToLowerInvariant()} } else { $null }
$Status = if ($TimedOut) { 'TIMED_OUT' } elseif ($null -ne $LaunchError) { 'LAUNCH_ERROR' } elseif (-not $PostPinsMatch) { 'POST_RUN_PIN_DRIFT' } elseif ($ChildExit -eq 0) { 'COMPLETED' } else { 'FAILED_OR_EXIT_UNKNOWN' }
$Receipt = [ordered]@{
    schema = 'pgfrag.pd.loaded_diagnostic.B_launch_receipt.v1'
    run_id = $RunId
    status = $Status
    wrapper_exit_code = if ($TimedOut) { 124 } elseif ($null -ne $LaunchError) { 125 } elseif (-not $PostPinsMatch) { 3 } elseif ($null -ne $ChildExit) { $ChildExit } else { 125 }
    child_exit_code = $ChildExit
    child_exit_observed = ($null -ne $ChildExit)
    timed_out = $TimedOut
    wait_completed = $WaitCompleted
    launch_error = $LaunchError
    timeout_seconds = [int]$Manifest.runtime.timeout_seconds
    julia_threads = 1
    runtime_sha256 = $RuntimeHash
    launcher_path = $LauncherPath
    launcher_sha256 = $LauncherHash
    manifest_sha256_before = $ManifestHash
    preflight_path = $PreflightPath
    preflight_sha256 = $PreflightHash
    started_utc = $StartedUtc
    finished_utc = $FinishedUtc
    elapsed_seconds = ([DateTime]::Parse($FinishedUtc).ToUniversalTime() - [DateTime]::Parse($StartedUtc).ToUniversalTime()).TotalSeconds
    argv = $Arguments
    snapshot_directory = $SnapshotDirectory
    stdout_log = $StdoutPath
    stderr_log = $StderrPath
    stdout = $StdoutInfo
    stderr = $StderrInfo
    source_pin_count_before = $PrePins.Count
    source_pins_before = $PrePins
    source_pin_count_after = $PostPins.Count
    source_pins_after = $PostPins
    post_pin_error = $PostPinError
    post_source_pins_match = $PostPinsMatch
    output_artifacts = $OutputArtifacts
    no_time_integration = $true
    no_manual_bond_activity_edits = $true
    physical_qualification = 'NONE_UNVALIDATED_PROTOTYPE'
    energy = 'UNKNOWN_NOT_A_GATE'
    contact = 'NOT_RUN'
}
[System.IO.File]::WriteAllText($ReceiptPath, ($Receipt | ConvertTo-Json -Depth 12) + [Environment]::NewLine, $Utf8NoBom)
Write-Output "B_STATUS=$Status wrapper=$($Receipt.wrapper_exit_code) child=$ChildExit timeout=$TimedOut pins=$($PostPins.Count)/$($PostPinsMatch)"
if ($TimedOut) { exit 124 }
if ($null -ne $LaunchError) { exit 125 }
if ($null -eq $ChildExit) { exit 125 }
if (-not $PostPinsMatch) { exit 3 }
exit $ChildExit
