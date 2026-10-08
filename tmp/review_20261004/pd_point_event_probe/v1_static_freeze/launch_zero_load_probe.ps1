$ErrorActionPreference = 'Stop'
if (Test-Path variable:PSNativeCommandUseErrorActionPreference) {
    $PSNativeCommandUseErrorActionPreference = $false
}

$Workspace = 'C:\Users\xiemo\Desktop\ShardForge'
$Julia = 'C:\Users\xiemo\Desktop\ShardForge\third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$JuliaProject = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\src\pgfrag\resources\julia'
$JuliaDepot = 'C:\Users\xiemo\Desktop\ShardForge\third_party\pgfrag-runtime\depot'
$Runner = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\pd_point_event_probe_20261004\run_pd_point_event_probe.jl'
$ManifestPath = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\pd_point_event_probe_20261004\manifest.json'
$PointsCsv = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\pd_point_event_probe_20261004\input\coupon_uniform_voxel_5cube_125.points.csv'
$PointMetadata = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\pd_point_event_probe_20261004\input\coupon_uniform_voxel_5cube_125.metadata.toml'
$CommonCase = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\pd_point_event_probe_20261004\common_case.toml'
$Case = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\pd_point_event_probe_20261004\case_zero_load_control.toml'
$RunId = 'zero_load_control_run01'
$OutDirectory = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\tmp\review_20261004\pd_point_event_probe\OUT'
$SnapshotDirectory = Join-Path 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\tmp\review_20261004\pd_point_event_probe\snapshots' $RunId
$ReceiptDirectory = 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\tmp\review_20261004\pd_point_event_probe\receipts'
$ReceiptPath = Join-Path $ReceiptDirectory ($RunId + '.launch_receipt.json')
$StdoutPath = Join-Path $ReceiptDirectory ($RunId + '.stdout.log')
$StderrPath = Join-Path $ReceiptDirectory ($RunId + '.stderr.log')

foreach ($RequiredPath in @($Julia, $JuliaProject, $JuliaDepot, $Runner, $ManifestPath,
                             $PointsCsv, $PointMetadata, $CommonCase, $Case)) {
    if (-not (Test-Path -LiteralPath $RequiredPath)) {
        throw "Required frozen path is missing: $RequiredPath"
    }
}

$Manifest = Get-Content -LiteralPath $ManifestPath -Raw | ConvertFrom-Json
foreach ($Record in $Manifest.file_hashes) {
    $RelativePath = [string]$Record.path
    $ExpectedHash = ([string]$Record.sha256).ToLowerInvariant()
    $AbsolutePath = Join-Path $Workspace ($RelativePath -replace '/', '\')
    if (-not (Test-Path -LiteralPath $AbsolutePath -PathType Leaf)) {
        throw "Pinned source file is missing: $RelativePath"
    }
    $ActualHash = (Get-FileHash -LiteralPath $AbsolutePath -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($ActualHash -ne $ExpectedHash) {
        throw "Frozen source hash mismatch for $RelativePath"
    }
}
$RuntimeHash = (Get-FileHash -LiteralPath $Julia -Algorithm SHA256).Hash.ToLowerInvariant()
if ($RuntimeHash -ne ([string]$Manifest.runtime.sha256).ToLowerInvariant()) {
    throw 'Frozen Julia executable hash does not match manifest.json'
}

if (-not (Test-Path -LiteralPath $OutDirectory -PathType Container)) {
    New-Item -ItemType Directory -Path $OutDirectory | Out-Null
}
if (@(Get-ChildItem -LiteralPath $OutDirectory -Force).Count -ne 0) {
    throw "OUT must be empty before launch: $OutDirectory"
}
if (Test-Path -LiteralPath $SnapshotDirectory) {
    throw "Snapshot directory already exists; refusing overwrite: $SnapshotDirectory"
}
if ((Test-Path -LiteralPath $ReceiptPath) -or (Test-Path -LiteralPath $StdoutPath) -or
    (Test-Path -LiteralPath $StderrPath)) {
    throw "Receipt or process log already exists for $RunId; refusing overwrite"
}
New-Item -ItemType Directory -Force -Path (Split-Path -Parent $SnapshotDirectory) | Out-Null
New-Item -ItemType Directory -Force -Path $ReceiptDirectory | Out-Null

$PreviousThreads = $env:JULIA_NUM_THREADS
$PreviousDepot = $env:JULIA_DEPOT_PATH
$PreviousPrecompileTasks = $env:JULIA_NUM_PRECOMPILE_TASKS
$env:JULIA_NUM_THREADS = '1'
$env:JULIA_DEPOT_PATH = $JuliaDepot
$env:JULIA_NUM_PRECOMPILE_TASKS = '1'
$StartedUtc = [DateTime]::UtcNow.ToString('o')
$CommandArguments = @("--project=$JuliaProject", '--startup-file=no', '--history-file=no',
                      $Runner, $SnapshotDirectory, $PointsCsv, $PointMetadata,
                      $CommonCase, $Case)
if (@($CommandArguments | Where-Object { $_ -match '[\s"]' }).Count -ne 0) {
    throw 'Frozen Windows launch arguments must not contain spaces or quotes'
}
$TimeoutSeconds = [int]$Manifest.runtime_gate.timeout_seconds
if ($TimeoutSeconds -lt 1) { throw 'Manifest runtime timeout must be positive' }
$TimedOut = $false
try {
    $Process = Start-Process -FilePath $Julia -ArgumentList ($CommandArguments -join ' ') `
        -WorkingDirectory $OutDirectory -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput $StdoutPath -RedirectStandardError $StderrPath
    $ExitedInTime = $Process.WaitForExit($TimeoutSeconds * 1000)
    if (-not $ExitedInTime) {
        $TimedOut = $true
        if (-not $Process.HasExited) { $Process.Kill() }
        $Process.WaitForExit()
    }
    $JuliaExitCode = [int]$Process.ExitCode
} finally {
    if ($null -eq $PreviousThreads) { Remove-Item Env:JULIA_NUM_THREADS -ErrorAction SilentlyContinue }
    else { $env:JULIA_NUM_THREADS = $PreviousThreads }
    if ($null -eq $PreviousDepot) { Remove-Item Env:JULIA_DEPOT_PATH -ErrorAction SilentlyContinue }
    else { $env:JULIA_DEPOT_PATH = $PreviousDepot }
    if ($null -eq $PreviousPrecompileTasks) { Remove-Item Env:JULIA_NUM_PRECOMPILE_TASKS -ErrorAction SilentlyContinue }
    else { $env:JULIA_NUM_PRECOMPILE_TASKS = $PreviousPrecompileTasks }
}
$FinishedUtc = [DateTime]::UtcNow.ToString('o')
$OutEmptyAfter = @(Get-ChildItem -LiteralPath $OutDirectory -Force).Count -eq 0
$Receipt = [ordered]@{
    schema = 'pgfrag.pd.point_event.launch_receipt.v1'
    run_id = $RunId
    case_id = 'zero_load_control'
    status = $(if ($TimedOut) { 'TIMED_OUT' } elseif ($JuliaExitCode -eq 0) { 'COMPLETED' } else { 'FAILED' })
    julia_process_exit_code = $JuliaExitCode
    timed_out = $TimedOut
    timeout_seconds = $TimeoutSeconds
    started_utc = $StartedUtc
    finished_utc = $FinishedUtc
    command_executable = $Julia
    command_arguments = @("--project=$JuliaProject", '--startup-file=no', '--history-file=no',
                          $Runner, $SnapshotDirectory, $PointsCsv, $PointMetadata,
                          $CommonCase, $Case)
    project_path = $JuliaProject
    depot_path = $JuliaDepot
    runtime_sha256 = $RuntimeHash
    julia_threads = 1
    out_directory = $OutDirectory
    out_empty_before_launch = $true
    out_empty_after_launch = $OutEmptyAfter
    snapshot_directory = $SnapshotDirectory
    stdout_log = $StdoutPath
    stderr_log = $StderrPath
}
$ReceiptText = ($Receipt | ConvertTo-Json -Depth 8) + "`n"
$Utf8NoBom = New-Object System.Text.UTF8Encoding($false)
[System.IO.File]::WriteAllText($ReceiptPath, $ReceiptText, $Utf8NoBom)

if (-not $OutEmptyAfter) {
    throw "OUT gained files during launch; inspect the sibling receipt and logs: $OutDirectory"
}
exit $JuliaExitCode
