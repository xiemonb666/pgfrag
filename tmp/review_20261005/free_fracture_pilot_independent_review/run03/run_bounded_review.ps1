param(
    [Parameter(Mandatory = $true)]
    [switch]$RootAuthorized,
    [Parameter(Mandatory = $true)]
    [ValidatePattern('^[0-9a-fA-F]{64}$')]
    [string]$ExpectedFinalSummarySha256
)
$ErrorActionPreference = 'Stop'
if (-not $RootAuthorized) {
    Write-Error 'Explicit parent authorization is required.'
    exit 3
}
$reviewDir = (Resolve-Path -LiteralPath $PSScriptRoot).Path
$repoRoot = [IO.Path]::GetFullPath((Join-Path $reviewDir '..\..\..\..\..'))
$python = Join-Path $repoRoot 'third_party\vesselforge-runtime\venv\Scripts\python.exe'
$adapter = Join-Path $reviewDir 'pilot_independent_review.py'
$outputDir = Join-Path $reviewDir 'audit_output'
[void][IO.Directory]::CreateDirectory($outputDir)
$stdoutPath = Join-Path $outputDir 'bounded_runner_stdout.txt'
$stderrPath = Join-Path $outputDir 'bounded_runner_stderr.txt'
$runnerReceiptPath = Join-Path $reviewDir 'bounded_runner_receipt.json'
if (-not (Test-Path -LiteralPath $python -PathType Leaf) -or
    -not (Test-Path -LiteralPath $adapter -PathType Leaf)) {
    throw 'Required bundled Python runtime or adapter is missing.'
}
$start = [Diagnostics.Stopwatch]::StartNew()
$psi = [Diagnostics.ProcessStartInfo]::new()
$psi.FileName = $python
$psi.WorkingDirectory = $repoRoot
$psi.UseShellExecute = $false
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.CreateNoWindow = $true
[void]$psi.ArgumentList.Add($adapter)
[void]$psi.ArgumentList.Add('--root-authorized')
[void]$psi.ArgumentList.Add('--expected-final-summary-sha256')
[void]$psi.ArgumentList.Add($ExpectedFinalSummarySha256.ToLowerInvariant())
$process = [Diagnostics.Process]::new()
$process.StartInfo = $psi
if (-not $process.Start()) { throw 'Could not start the read-only review adapter.' }
$stdoutTask = $process.StandardOutput.ReadToEndAsync()
$stderrTask = $process.StandardError.ReadToEndAsync()
$completed = $process.WaitForExit(110000)
if (-not $completed) {
    try { $process.Kill($true) } catch { try { $process.Kill() } catch {} }
    $process.WaitForExit()
    $stdoutText = $stdoutTask.GetAwaiter().GetResult()
    $stderrText = $stderrTask.GetAwaiter().GetResult()
    [IO.File]::WriteAllText($stdoutPath, $stdoutText, [Text.UTF8Encoding]::new($false))
    [IO.File]::WriteAllText($stderrPath, $stderrText, [Text.UTF8Encoding]::new($false))
    $start.Stop()
    $receipt = [ordered]@{
        schema = 'pgfrag.free_fracture_pilot.bounded_review_runner.v1'
        status = 'TIMEOUT_KILLED'
        runtime_budget_seconds = 110
        elapsed_ms = $start.ElapsedMilliseconds
        adapter_exit_code = $null
        expected_summary_sha256 = $ExpectedFinalSummarySha256.ToLowerInvariant()
        stdout_path = $stdoutPath
        stderr_path = $stderrPath
    }
    $receipt | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $runnerReceiptPath -Encoding utf8
    exit 124
}
$process.WaitForExit()
$stdoutText = $stdoutTask.GetAwaiter().GetResult()
$stderrText = $stderrTask.GetAwaiter().GetResult()
[IO.File]::WriteAllText($stdoutPath, $stdoutText, [Text.UTF8Encoding]::new($false))
[IO.File]::WriteAllText($stderrPath, $stderrText, [Text.UTF8Encoding]::new($false))
$start.Stop()
$receipt = [ordered]@{
    schema = 'pgfrag.free_fracture_pilot.bounded_review_runner.v1'
    status = if ($process.ExitCode -eq 0) { 'COMPLETED' } else { 'COMPLETED_NONZERO' }
    runtime_budget_seconds = 110
    elapsed_ms = $start.ElapsedMilliseconds
    adapter_exit_code = $process.ExitCode
    expected_summary_sha256 = $ExpectedFinalSummarySha256.ToLowerInvariant()
    stdout_path = $stdoutPath
    stderr_path = $stderrPath
}
$receipt | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $runnerReceiptPath -Encoding utf8
exit $process.ExitCode
