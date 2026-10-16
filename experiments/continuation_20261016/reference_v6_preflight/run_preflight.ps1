param()
$ErrorActionPreference = 'Stop'

$experiment = $PSScriptRoot
$pgfrag = Split-Path (Split-Path (Split-Path $experiment -Parent) -Parent) -Parent
$data = Join-Path $pgfrag 'data\continuation_20261016\reference_v6_preflight'
$controlPath = Join-Path $data 'run_control.json'
$budgetPath = Join-Path $data 'preflight_budget.json'
$nativeRoot = Join-Path $experiment 'native_preflight'
$run = Join-Path $nativeRoot 'run01'
$stdout = Join-Path $nativeRoot 'stdout.log'
$stderr = Join-Path $nativeRoot 'stderr.log'
$watchdogPath = Join-Path $data 'process_watchdog.json'
$executionPath = Join-Path $data 'execution_record.json'

function Get-FileDigest($path) {
    if (-not (Test-Path -LiteralPath $path)) { return $null }
    $file = Get-Item -LiteralPath $path
    $hash = (Get-FileHash -Algorithm SHA256 -LiteralPath $path).Hash.ToLowerInvariant()
    return @{ bytes = $file.Length; sha256 = $hash }
}

if (-not (Test-Path -LiteralPath $controlPath) -or -not (Test-Path -LiteralPath $budgetPath)) {
    throw 'Prepare the frozen inputs and pass the full-output budget forecast before running.'
}
$control = Get-Content -Raw -LiteralPath $controlPath | ConvertFrom-Json
$budget = Get-Content -Raw -LiteralPath $budgetPath | ConvertFrom-Json
if ($control.status -ne 'READY_FOR_SINGLE_NATIVE_INITIALIZATION_PREFLIGHT' -or $budget.status -ne 'WHOLE_OUTPUT_FORECAST_PASS_RUN_NOT_STARTED' -or -not $budget.native_preflight_output_forecast.pass) {
    throw 'The frozen preflight is not in the ready state or its 256 MiB forecast did not pass.'
}
$moduleNames = @('cohesive_vessel.py', 'cohesive_mesh.py', 'cohesive_stability.py')
$repoPg = $pgfrag
$currentSourceHashesBefore = [ordered]@{}
foreach ($name in $moduleNames) {
    $currentPath = Join-Path $repoPg "src\pgfrag\$name"
    $currentSourceHashesBefore[$name] = Get-FileDigest $currentPath
    if ($currentSourceHashesBefore[$name].sha256 -ne $control.input_hashes.$name) {
        throw "Current native source differs from frozen case before execution: $name"
    }
}
if (Test-Path -LiteralPath $run -or Test-Path -LiteralPath $watchdogPath -or Test-Path -LiteralPath $executionPath) {
    throw 'A run attempt or execution record already exists. This package permits one attempt and no retry.'
}
$runnerHashesBefore = [ordered]@{}
foreach ($name in $control.runner_script_hashes.PSObject.Properties.Name) {
    $scriptPath = Join-Path $experiment $name
    $actualHash = if (Test-Path -LiteralPath $scriptPath) { (Get-FileHash -Algorithm SHA256 -LiteralPath $scriptPath).Hash.ToLowerInvariant() } else { $null }
    $runnerHashesBefore[$name] = $actualHash
    if ($actualHash -ne $control.runner_script_hashes.$name) {
        throw "Runner script hash mismatch before native execution: $name"
    }
}

New-Item -ItemType Directory -Force -Path $nativeRoot | Out-Null
$root = (Resolve-Path -LiteralPath (Join-Path $experiment '..\..\..\..')).Path
$bashRunner = '/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/continuation_20261016/reference_v6_preflight/run_native_preflight.sh'
$runtime = '/mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python'
$runtimeProbe = @'
import sys, numpy, akantu
print(f"python={sys.version.split()[0]} numpy={numpy.__version__} akantu={akantu.__version__}")
'@
$probeResult = & wsl.exe -d Ubuntu -- $runtime -c $runtimeProbe
if ($LASTEXITCODE -ne 0 -or ($probeResult -join '') -notmatch 'python=3\.12\.15 numpy=2\.5\.3 akantu=5\.0\.7-post1') {
    throw "Expected WSL Akantu runtime did not verify: $($probeResult -join ' ')"
}

$started = [DateTime]::UtcNow
$clock = [System.Diagnostics.Stopwatch]::StartNew()
$process = Start-Process -FilePath 'wsl.exe' -ArgumentList @('-d', 'Ubuntu', '--', 'bash', $bashRunner) `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr -PassThru -WindowStyle Hidden
Write-Output "Started one v6 initialization-only attempt (WSL host PID $($process.Id)); watchdog 330 s."

$timedOut = $false
$watchdogSeconds = [int]$control.external_watchdog_s
while (-not $process.WaitForExit(2000)) {
    if ($clock.Elapsed.TotalSeconds -ge $watchdogSeconds) {
        $timedOut = $true
        try { $process.Kill($true) } catch { Stop-Process -Id $process.Id -Force -ErrorAction SilentlyContinue }
        break
    }
}
if ($timedOut) {
    $process.WaitForExit(10000) | Out-Null
}
$clock.Stop()
$process.Refresh()
$exitCode = if ($process.HasExited -and -not $timedOut) { $process.ExitCode } else { $null }

function Write-JsonNoBom($path, $value) {
    $json = $value | ConvertTo-Json -Depth 14
    [System.IO.File]::WriteAllText($path, $json + [Environment]::NewLine, [System.Text.UTF8Encoding]::new($false))
}

$inputDir = Join-Path $experiment 'inputs'
$frozenHashes = [ordered]@{}
Get-ChildItem -LiteralPath $inputDir -File | Sort-Object Name | ForEach-Object {
    $frozenHashes[$_.Name] = Get-FileDigest $_.FullName
}
$currentSourceHashesAfter = [ordered]@{}
foreach ($name in $moduleNames) {
    $currentPath = Join-Path $repoPg "src\pgfrag\$name"
    $currentSourceHashesAfter[$name] = Get-FileDigest $currentPath
}
$runnerHashesAfter = [ordered]@{}
$runnerScriptsUnchanged = $true
foreach ($name in $control.runner_script_hashes.PSObject.Properties.Name) {
    $scriptPath = Join-Path $experiment $name
    $runnerHashesAfter[$name] = if (Test-Path -LiteralPath $scriptPath) { (Get-FileHash -Algorithm SHA256 -LiteralPath $scriptPath).Hash.ToLowerInvariant() } else { $null }
    if ($runnerHashesAfter[$name] -ne $runnerHashesBefore[$name] -or $runnerHashesAfter[$name] -ne $control.runner_script_hashes.$name) { $runnerScriptsUnchanged = $false }
}
$moduleBytesUnchanged = $true
foreach ($name in $moduleNames) {
    if ($currentSourceHashesBefore[$name].sha256 -ne $currentSourceHashesAfter[$name].sha256 -or
        $currentSourceHashesAfter[$name].sha256 -ne $frozenHashes[$name].sha256) { $moduleBytesUnchanged = $false }
}
$casePath = Join-Path $inputDir 'case.json'
$stdoutDigest = Get-FileDigest $stdout
$stderrDigest = Get-FileDigest $stderr
$runControlSnapshot = @{
    schema = 'pgfrag.reference_v6_native_initialization_process_watchdog.v1'
    status = if ($timedOut) { 'EXTERNAL_WATCHDOG_TIMEOUT_PRESERVED' } elseif ($exitCode -eq 0) { 'PROCESS_COMPLETED' } else { 'PROCESS_FAILED_PRESERVED' }
    started_utc = $started.ToString('o')
    ended_utc = [DateTime]::UtcNow.ToString('o')
    wall_elapsed_s = [Math]::Round($clock.Elapsed.TotalSeconds, 6)
    external_watchdog_s = $watchdogSeconds
    timed_out = $timedOut
    exit_code = $exitCode
    wsl_host_pid = $process.Id
    runtime_probe = ($probeResult -join '').Trim()
    runtime_python = $runtime
    command_argv = $control.command_argv
    command_shell = 'bash run_native_preflight.sh; timeout --foreground --signal=TERM --kill-after=5s 300s; taskset -c 0; external PowerShell watchdog 330s'
    stdout = $stdoutDigest
    stderr = $stderrDigest
    frozen_input_hashes_after_process = $frozenHashes
    frozen_case_sha256_after_process = (Get-FileDigest $casePath).sha256
    current_source_hashes_before_process = $currentSourceHashesBefore
    current_source_hashes_after_process = $currentSourceHashesAfter
    runner_script_hashes_before = $runnerHashesBefore
    runner_script_hashes_after = $runnerHashesAfter
    runner_scripts_unchanged = $runnerScriptsUnchanged
    current_module_bytes_match_frozen_copies = $moduleBytesUnchanged
    run_directory_created = Test-Path -LiteralPath $run
    native_result_present = Test-Path -LiteralPath (Join-Path $run 'result.json')
    retry_forbidden = $true
}
Write-JsonNoBom $watchdogPath $runControlSnapshot

$execution = [ordered]@{
    schema = 'pgfrag.reference_v6_native_initialization_preflight_execution.v1'
    status = $runControlSnapshot.status
    run_mode = 'NATIVE_INITIALIZATION_ONLY_ZERO_SOLVER_STEPS'
    command_argv = $control.command_argv
    command_shell = $runControlSnapshot.command_shell
    runtime = $control.runtime
    runtime_probe = $runControlSnapshot.runtime_probe
    native_cpu_affinity = 'CPU 0 via taskset -c 0'
    thread_environment = $control.thread_environment
    external_watchdog_s = $watchdogSeconds
    wrapper_wall_time_s = $control.wrapper_wall_time_s
    elapsed_wall_s = $runControlSnapshot.wall_elapsed_s
    exit_code = $exitCode
    watchdog_timed_out = $timedOut
    case_sha256_before = $control.case_sha256
    case_sha256_after = $runControlSnapshot.frozen_case_sha256_after_process
    input_hashes_before = $control.input_hashes
    input_hashes_after = $frozenHashes
    source_geometry_hashes = $control.source_geometry_hashes
    current_source_hashes_before_process = $currentSourceHashesBefore
    current_source_hashes_after_process = $currentSourceHashesAfter
    runner_script_hashes_before = $runnerHashesBefore
    runner_script_hashes_after = $runnerHashesAfter
    runner_scripts_unchanged = $runnerScriptsUnchanged
    source_modules_unchanged = $moduleBytesUnchanged
    stdout = $stdoutDigest
    stderr = $stderrDigest
    process_watchdog_path = 'process_watchdog.json'
    process_watchdog_sha256 = (Get-FileDigest $watchdogPath).sha256
    native_run_directory = 'native_preflight/run01'
    native_result_present = $runControlSnapshot.native_result_present
    retry_count = 0
    retry_or_retune_performed = $false
}
Write-JsonNoBom $executionPath $execution
Write-Output "Attempt finished: status=$($execution.status), exit=$exitCode, elapsed=$($execution.elapsed_wall_s)s, modules_unchanged=$moduleBytesUnchanged."
if ($timedOut) { exit 124 }
if ($exitCode -ne 0) { exit 2 }
