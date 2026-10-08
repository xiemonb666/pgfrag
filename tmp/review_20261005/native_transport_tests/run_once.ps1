$ErrorActionPreference = 'Stop'
$ProgressPreference = 'SilentlyContinue'

$outDir = $PSScriptRoot
$repoRoot = [System.IO.Path]::GetFullPath((Join-Path $outDir '..\..\..\..'))
$python = Join-Path $repoRoot 'relicweave-synlab\.venv\Scripts\python.exe'
$sourcePath = 'pgfrag/experiments/material_visibility_native_bridge_20261005/native_mask_transport.py'
$testPath = 'pgfrag/experiments/material_visibility_native_bridge_20261005/test_native_mask_transport.py'
$docPath = 'pgfrag/experiments/material_visibility_native_bridge_20261005/TRANSPORT.md'
$expectedSource = '6779b3242ef87c13d5626a1e6d89414e1cdd7d6ed8f7e264fb01ed83edef7ea9'
$expectedTest = '605cc23b38616fdcb848898c907d3e3677df7e11bdbd137c927a0c0342fc6643'
$expectedDoc = '1037dae0b326f76e62591a3bf244f527a724380631db2225e371b6006c77d685'
$stdoutPath = Join-Path $outDir 'native_transport.stdout.txt'
$stderrPath = Join-Path $outDir 'native_transport.stderr.txt'
$junitPath = Join-Path $outDir 'native_transport.junit.xml'
$exitPath = Join-Path $outDir 'exit_code.txt'
$nodeIdsPath = Join-Path $outDir 'distinct_testcase_ids.txt'
$manifestPath = Join-Path $outDir 'run_manifest.json'

function Get-Sha256([string]$Path) {
    return (Get-FileHash -Algorithm SHA256 -LiteralPath $Path).Hash.ToLowerInvariant()
}

if (-not (Test-Path -LiteralPath $python -PathType Leaf)) {
    throw "SynLab Python not found: $python"
}
$sourceAbsolute = Join-Path $repoRoot $sourcePath
$testAbsolute = Join-Path $repoRoot $testPath
$docAbsolute = Join-Path $repoRoot $docPath
if (-not (Test-Path -LiteralPath $sourceAbsolute -PathType Leaf) -or
    -not (Test-Path -LiteralPath $testAbsolute -PathType Leaf) -or
    -not (Test-Path -LiteralPath $docAbsolute -PathType Leaf)) {
    throw 'Frozen source, test, or documentation file is missing.'
}

$preHashes = @{
    native_mask_transport_py = Get-Sha256 $sourceAbsolute
    test_native_mask_transport_py = Get-Sha256 $testAbsolute
    transport_md = Get-Sha256 $docAbsolute
}
if ($preHashes.native_mask_transport_py -ne $expectedSource -or
    $preHashes.test_native_mask_transport_py -ne $expectedTest -or
    $preHashes.transport_md -ne $expectedDoc) {
    throw "Pre-run SHA mismatch; refusing test execution. Observed: $($preHashes | ConvertTo-Json -Compress)"
}

$argv = @(
    '-B', '-m', 'pytest', '-rA', '-p', 'no:cacheprovider',
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
    native_mask_transport_py = Get-Sha256 $sourceAbsolute
    test_native_mask_transport_py = Get-Sha256 $testAbsolute
    transport_md = Get-Sha256 $docAbsolute
}

$caseIds = @()
$caseCount = 0
$summary = @{ tests_attribute = 0; failures = 0; errors = 0; skipped = 0 }
if (Test-Path -LiteralPath $junitPath -PathType Leaf) {
    [xml]$junit = Get-Content -LiteralPath $junitPath -Raw
    foreach ($suite in @($junit.SelectNodes('//testsuite'))) {
        $summary.tests_attribute += [int]$suite.tests
        $summary.failures += [int]$suite.failures
        $summary.errors += [int]$suite.errors
        $summary.skipped += [int]$suite.skipped
    }
    $testcaseNodes = @($junit.SelectNodes('//testcase'))
    $caseCount = $testcaseNodes.Count
    foreach ($case in $testcaseNodes) {
        $caseIds += ('{0}::{1}' -f $case.classname, $case.name)
    }
}
$distinctIds = @($caseIds | Sort-Object -Unique)
[System.IO.File]::WriteAllLines($nodeIdsPath, $distinctIds, [System.Text.UTF8Encoding]::new($false))

$receipt = [ordered]@{
    test_file = $testPath
    source_file = $sourcePath
    documentation_file = $docPath
    python_executable = $python
    argv = $argv
    invocation_count = 1
    started_utc = $startedUtc.ToString('O')
    finished_utc = $finishedUtc.ToString('O')
    elapsed_seconds = [Math]::Round($clock.Elapsed.TotalSeconds, 6)
    exit_code = $testExitCode
    pre_sha256 = $preHashes
    post_sha256 = $postHashes
    all_three_hashes_unchanged = (($preHashes.native_mask_transport_py -eq $postHashes.native_mask_transport_py) -and
                                  ($preHashes.test_native_mask_transport_py -eq $postHashes.test_native_mask_transport_py) -and
                                  ($preHashes.transport_md -eq $postHashes.transport_md))
    junit = if (Test-Path -LiteralPath $junitPath -PathType Leaf) {
        @{ path = [System.IO.Path]::GetFileName($junitPath); sha256 = Get-Sha256 $junitPath }
    } else { $null }
    stdout = @{ path = [System.IO.Path]::GetFileName($stdoutPath); sha256 = Get-Sha256 $stdoutPath }
    stderr = @{ path = [System.IO.Path]::GetFileName($stderrPath); sha256 = Get-Sha256 $stderrPath }
    summary = $summary
    actual_testcase_xml_element_count = $caseCount
    distinct_testcase_ids = $distinctIds
    distinct_testcase_id_count = $distinctIds.Count
    scope_note = 'Serialization and declared-ledger validation only; raw-horizon/native application remains unverified.'
}
$json = ConvertTo-Json -InputObject $receipt -Depth 8
[System.IO.File]::WriteAllText($manifestPath, "$json`n", [System.Text.UTF8Encoding]::new($false))

Write-Output (ConvertTo-Json -InputObject $receipt -Depth 8)
exit $testExitCode
