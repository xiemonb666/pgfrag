$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'

$outDir = $PSScriptRoot
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $outDir '..\..\..\..'))
$python = Join-Path $repoRoot 'relicweave-synlab\.venv\Scripts\python.exe'
$sourcePath = 'pgfrag/experiments/material_visibility_native_bridge_20261005/bridge_contract.py'
$testPath = 'pgfrag/experiments/material_visibility_native_bridge_20261005/test_bridge_contract.py'
$expectedSource = '67b8f3090abebb093b7a1e099059c35c2ae09890b4f6b78453dd44d8a6ad4085'
$expectedTest = 'cb9c9b62214c6be2f5d9368dff9498f90a16c16537e135d9dd4a7c1e9cd09d07'
$stdoutPath = Join-Path $outDir 'bridge_contract.stdout.txt'
$stderrPath = Join-Path $outDir 'bridge_contract.stderr.txt'
$junitPath = Join-Path $outDir 'bridge_contract.junit.xml'
$exitPath = Join-Path $outDir 'exit_code.txt'
$nodeIdsPath = Join-Path $outDir 'distinct_nodeids.txt'
$manifestPath = Join-Path $outDir 'run_manifest.json'

function Get-Sha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

if (-not (Test-Path -LiteralPath $python -PathType Leaf)) {
    throw "SynLab Python not found: $python"
}
$sourceAbsolute = Join-Path $repoRoot $sourcePath
$testAbsolute = Join-Path $repoRoot $testPath
if (-not (Test-Path -LiteralPath $sourceAbsolute -PathType Leaf) -or
    -not (Test-Path -LiteralPath $testAbsolute -PathType Leaf)) {
    throw 'Frozen source or test file is missing.'
}

$preHashes = @{
    bridge_contract_py = Get-Sha256 $sourceAbsolute
    test_bridge_contract_py = Get-Sha256 $testAbsolute
}
if ($preHashes.bridge_contract_py -ne $expectedSource -or
    $preHashes.test_bridge_contract_py -ne $expectedTest) {
    throw "Pre-run SHA mismatch; refusing test execution. Observed: $($preHashes | ConvertTo-Json -Compress)"
}

$argv = @(
    '-B', '-m', 'pytest', '-q', '-rA',
    "--junitxml=$junitPath",
    $testPath
)
$startedUtc = [DateTimeOffset]::UtcNow
$clock = [System.Diagnostics.Stopwatch]::StartNew()
Push-Location -LiteralPath $repoRoot
try {
    & $python @argv 1> $stdoutPath 2> $stderrPath
    $testExitCode = $LASTEXITCODE
}
finally {
    Pop-Location
    $clock.Stop()
}
$finishedUtc = [DateTimeOffset]::UtcNow

[System.IO.File]::WriteAllText($exitPath, "$testExitCode`n", [System.Text.UTF8Encoding]::new($false))
$postHashes = @{
    bridge_contract_py = Get-Sha256 $sourceAbsolute
    test_bridge_contract_py = Get-Sha256 $testAbsolute
}

$caseIds = @()
$summary = @{ tests = 0; failures = 0; errors = 0; skipped = 0 }
if (Test-Path -LiteralPath $junitPath -PathType Leaf) {
    [xml]$junit = Get-Content -LiteralPath $junitPath -Raw
    foreach ($suite in @($junit.SelectNodes('//testsuite'))) {
        $summary.tests += [int]$suite.tests
        $summary.failures += [int]$suite.failures
        $summary.errors += [int]$suite.errors
        $summary.skipped += [int]$suite.skipped
    }
    foreach ($case in @($junit.SelectNodes('//testcase'))) {
        $caseIds += ('{0}::{1}' -f $case.classname, $case.name)
    }
}
$distinctIds = @($caseIds | Sort-Object -Unique)
[System.IO.File]::WriteAllLines($nodeIdsPath, $distinctIds, [System.Text.UTF8Encoding]::new($false))

$receipt = [ordered]@{
    test_file = $testPath
    source_file = $sourcePath
    python_executable = $python
    argv = $argv
    invocation_count = 1
    started_utc = $startedUtc.ToString('O')
    finished_utc = $finishedUtc.ToString('O')
    elapsed_seconds = [Math]::Round($clock.Elapsed.TotalSeconds, 6)
    exit_code = $testExitCode
    pre_sha256 = $preHashes
    post_sha256 = $postHashes
    source_hashes_unchanged = (($preHashes.bridge_contract_py -eq $postHashes.bridge_contract_py) -and
                               ($preHashes.test_bridge_contract_py -eq $postHashes.test_bridge_contract_py))
    junit = if (Test-Path -LiteralPath $junitPath -PathType Leaf) {
        @{ path = [System.IO.Path]::GetFileName($junitPath); sha256 = Get-Sha256 $junitPath }
    } else { $null }
    stdout = @{ path = [System.IO.Path]::GetFileName($stdoutPath); sha256 = Get-Sha256 $stdoutPath }
    stderr = @{ path = [System.IO.Path]::GetFileName($stderrPath); sha256 = Get-Sha256 $stderrPath }
    summary = $summary
    distinct_node_ids = $distinctIds
    distinct_node_id_count = $distinctIds.Count
}
$json = ConvertTo-Json -InputObject $receipt -Depth 8
[System.IO.File]::WriteAllText($manifestPath, "$json`n", [System.Text.UTF8Encoding]::new($false))

Write-Output (ConvertTo-Json -InputObject $receipt -Depth 8)
exit $testExitCode
