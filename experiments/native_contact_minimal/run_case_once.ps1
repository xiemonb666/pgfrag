param(
    [Parameter(Mandatory = $true)]
    [ValidateSet('preflight', 'contact_baseline', 'contact_halfdt', 'nocontact_control')]
    [string]$CaseName
)

$ErrorActionPreference = 'Stop'
$workspace = (Resolve-Path (Join-Path $PSScriptRoot '..\..\..')).Path
$julia = Join-Path $workspace 'third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$project = 'pgfrag/src/pgfrag/resources/julia'
$nativeSource = Join-Path $PSScriptRoot 'run_native_contact_minimal.jl'
$resultRoot = Join-Path $workspace 'relicweave-synlab\tmp\review_20261004\implementation\native_contact_minimal'
$outputDir = Join-Path $resultRoot $CaseName
New-Item -ItemType Directory -Force -Path $outputDir | Out-Null

$argv = @(
    "--project=$project",
    '--compiled-modules=no',
    '--startup-file=no',
    '--history-file=no',
    '-t1',
    $nativeSource,
    $CaseName,
    $outputDir
)

$psi = [System.Diagnostics.ProcessStartInfo]::new()
$psi.FileName = $julia
$psi.WorkingDirectory = $workspace
$psi.UseShellExecute = $false
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.EnvironmentVariables['JULIA_DEPOT_PATH'] = 'third_party/pgfrag-runtime/depot'
$psi.Arguments = ($argv | ForEach-Object {
    if ($_ -match '[\s"]') { '"' + $_.Replace('"', '\"') + '"' } else { $_ }
}) -join ' '

$process = [System.Diagnostics.Process]::new()
$process.StartInfo = $psi
$clock = [System.Diagnostics.Stopwatch]::StartNew()
if (-not $process.Start()) { throw 'Julia process failed to start.' }
$stdoutTask = $process.StandardOutput.ReadToEndAsync()
$stderrTask = $process.StandardError.ReadToEndAsync()
$exited = $process.WaitForExit(120000)
if (-not $exited) {
    $process.Kill($true)
    $process.WaitForExit()
}
$clock.Stop()
$stdout = $stdoutTask.GetAwaiter().GetResult()
$stderr = $stderrTask.GetAwaiter().GetResult()
$exitCode = if ($exited) { $process.ExitCode } else { -1 }
[System.IO.File]::WriteAllText((Join-Path $outputDir 'stdout.log'), $stdout, [System.Text.UTF8Encoding]::new($false))
[System.IO.File]::WriteAllText((Join-Path $outputDir 'stderr.log'), $stderr, [System.Text.UTF8Encoding]::new($false))

function ConvertTo-TomlString([string]$value) {
    return '"' + $value.Replace('\', '\\').Replace('"', '\"').Replace("`r", '\r').Replace("`n", '\n') + '"'
}

$argvToml = ($argv | ForEach-Object { ConvertTo-TomlString $_ }) -join ', '
$timedOutToml = if ($exited) { 'false' } else { 'true' }
$manifestLines = @(
    "case = $(ConvertTo-TomlString $CaseName)",
    "native_source = $(ConvertTo-TomlString $nativeSource)",
    "julia_executable = $(ConvertTo-TomlString $julia)",
    "working_directory = $(ConvertTo-TomlString $workspace)",
    "depot_path = $(ConvertTo-TomlString 'third_party/pgfrag-runtime/depot')",
    "argv = [$argvToml]",
    "exit_code = $exitCode",
    "timed_out_at_120_seconds = $timedOutToml",
    "runtime_seconds = $([Math]::Round($clock.Elapsed.TotalSeconds, 6))",
    "stdout_file = $(ConvertTo-TomlString (Join-Path $outputDir 'stdout.log'))",
    "stderr_file = $(ConvertTo-TomlString (Join-Path $outputDir 'stderr.log'))"
)
[System.IO.File]::WriteAllText((Join-Path $outputDir 'execution.toml'), ($manifestLines -join "`n") + "`n", [System.Text.UTF8Encoding]::new($false))
$process.Dispose()

Write-Output "case=$CaseName exit_code=$exitCode runtime_seconds=$([Math]::Round($clock.Elapsed.TotalSeconds, 3)) timed_out=$(-not $exited) output=$outputDir"
if ($stdout) { Write-Output $stdout.TrimEnd() }
if ($stderr) { Write-Output $stderr.TrimEnd() }
if (-not $exited) { exit 124 }
exit $exitCode
