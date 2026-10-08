$ErrorActionPreference = 'Stop'

$workspace = 'C:\Users\xiemo\Desktop\ShardForge'
$scratchRoot = Join-Path $workspace 'pgfrag\tmp\initial_bond_policy_body_pilot_20261005'
$attemptRoot = Join-Path $scratchRoot 'run_native_attempt_01'
$auditDir = Join-Path $attemptRoot 'audit'
$bundleRoot = Join-Path $attemptRoot 'bundles'
$freshDepot = Join-Path $attemptRoot 'depot'
$pilotExp = Join-Path $workspace 'pgfrag\experiments\initial_bond_policy_body_pilot_20261005'
$nativeExp = Join-Path $workspace 'pgfrag\experiments\initial_bond_policy_native_20261005'
$run03 = Join-Path $workspace 'pgfrag\tmp\initial_bond_policy_native_20261005\fork_run_03'
$sharedDepot = Join-Path $workspace 'third_party\pgfrag-runtime\depot'
$vendorRoot = Join-Path $workspace 'pgfrag\src\pgfrag\resources\julia\vendor\Peridynamics.jl'
$appProject = Join-Path $workspace 'pgfrag\Project.toml'
$appManifest = Join-Path $workspace 'pgfrag\Manifest.toml'
$julia = Join-Path $workspace 'third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$python = 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe'
$generator = Join-Path $pilotExp 'generate_fixture_policy.py'
$driver = Join-Path $pilotExp 'native_body_pilot.jl'
$pilotPins = Join-Path $pilotExp 'source_reference_pins.md'
$nativePins = Join-Path $nativeExp 'source_pin_manifest.md'
$utf8NoBom = [System.Text.UTF8Encoding]::new($false)

if (Test-Path -LiteralPath $attemptRoot) {
    throw "Refusing to reuse an existing pilot attempt path: $attemptRoot"
}
New-Item -ItemType Directory -Path $auditDir -Force | Out-Null
New-Item -ItemType Directory -Path $freshDepot -Force | Out-Null

function Get-FileSha256([string] $Path) {
    $algorithm = [System.Security.Cryptography.SHA256]::Create()
    $stream = [System.IO.File]::OpenRead($Path)
    try {
        return [Convert]::ToHexString($algorithm.ComputeHash($stream)).ToLowerInvariant()
    }
    finally {
        $stream.Dispose()
        $algorithm.Dispose()
    }
}

function Get-TreeSnapshot([string] $Label, [string] $Path, [string] $Phase) {
    if (-not (Test-Path -LiteralPath $Path -PathType Container)) {
        throw "Required tree is missing: $Path"
    }
    $rootFull = [System.IO.Path]::GetFullPath($Path).TrimEnd('\')
    $files = @(Get-ChildItem -LiteralPath $rootFull -File -Recurse -Force | Sort-Object FullName -CaseSensitive)
    $records = [System.Collections.Generic.List[string]]::new()
    [long] $bytes = 0
    foreach ($file in $files) {
        $relative = $file.FullName.Substring($rootFull.Length + 1).Replace('\', '/')
        $sha = Get-FileSha256 $file.FullName
        $records.Add("$relative`t$($file.Length)`t$sha")
        $bytes += $file.Length
    }
    $manifestPath = Join-Path $auditDir "${Phase}_tree_${Label}.tsv"
    [System.IO.File]::WriteAllLines($manifestPath, $records, $utf8NoBom)
    return [pscustomobject]@{
        label = $Label
        path = $rootFull
        file_count = $files.Count
        total_bytes = $bytes
        manifest_path = $manifestPath
        manifest_sha256 = Get-FileSha256 $manifestPath
    }
}

function Resolve-NativePinPath([string] $RelativePath) {
    if ($RelativePath.StartsWith('pgfrag/', [System.StringComparison]::Ordinal) -or
        $RelativePath.StartsWith('third_party/', [System.StringComparison]::Ordinal)) {
        return [System.IO.Path]::GetFullPath((Join-Path $workspace ($RelativePath.Replace('/', '\'))))
    }
    return [System.IO.Path]::GetFullPath((Join-Path $nativeExp ($RelativePath.Replace('/', '\'))))
}

function Read-TablePins([string] $ManifestPath, [string] $StopHeading = '') {
    $pins = [System.Collections.Generic.List[object]]::new()
    foreach ($line in [System.IO.File]::ReadAllLines($ManifestPath)) {
        if (-not [string]::IsNullOrEmpty($StopHeading) -and
            $line.StartsWith($StopHeading, [System.StringComparison]::Ordinal)) { break }
        if ($line -match '^\|\s*[^|]+\|\s*`([^`]+)`\s*\|\s*`([0-9a-f]{64})`\s*\|' -or
            $line -match '^\|\s*`([^`]+)`\s*\|\s*`([0-9a-f]{64})`\s*\|') {
            $pins.Add([pscustomobject]@{ relative_path = $Matches[1]; expected_sha256 = $Matches[2] })
        }
    }
    return $pins.ToArray()
}

function Resolve-PilotPinPath([string] $RelativePath) {
    $candidate = [System.IO.Path]::GetFullPath((Join-Path $pilotExp ($RelativePath.Replace('/', '\'))))
    return $candidate
}

function Start-CapturedProcess([string] $Name, [string] $FileName, [string[]] $Arguments,
                               [string] $WorkingDirectory, [hashtable] $Environment,
                               [int] $TimeoutMilliseconds) {
    $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
    $startInfo.FileName = $FileName
    $startInfo.WorkingDirectory = $WorkingDirectory
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    foreach ($argument in $Arguments) { $startInfo.ArgumentList.Add($argument) }
    foreach ($entry in $Environment.GetEnumerator()) { $startInfo.Environment[$entry.Key] = [string]$entry.Value }

    $process = [System.Diagnostics.Process]::new()
    $process.StartInfo = $startInfo
    $startedUtc = [DateTimeOffset]::UtcNow
    $timer = [System.Diagnostics.Stopwatch]::StartNew()
    if (-not $process.Start()) { throw "Could not start $Name process: $FileName" }
    $processId = $process.Id
    $stdoutTask = $process.StandardOutput.ReadToEndAsync()
    $stderrTask = $process.StandardError.ReadToEndAsync()
    $timedOut = -not $process.WaitForExit($TimeoutMilliseconds)
    if ($timedOut) {
        try { $process.Kill($true) } catch { }
        [void]$process.WaitForExit(5000)
    }
    else {
        $process.WaitForExit()
    }
    $timer.Stop()
    $stdout = $stdoutTask.GetAwaiter().GetResult()
    $stderr = $stderrTask.GetAwaiter().GetResult()
    $endedUtc = [DateTimeOffset]::UtcNow
    $exitCode = if ($process.HasExited) { $process.ExitCode } else { $null }
    $process.Dispose()
    return [pscustomobject]@{
        name = $Name
        executable = $FileName
        argv = @($FileName) + @($Arguments)
        working_directory = $WorkingDirectory
        environment_overrides = $Environment
        started_utc = $startedUtc.ToString('o')
        ended_utc = $endedUtc.ToString('o')
        duration_ms = $timer.ElapsedMilliseconds
        timeout_ms = $TimeoutMilliseconds
        timed_out = $timedOut
        process_id = $processId
        exit_code = $exitCode
        stdout = $stdout
        stderr = $stderr
    }
}

function Write-Json([string] $Path, $Value) {
    $json = ConvertTo-Json -InputObject $Value -Depth 32
    [System.IO.File]::WriteAllText($Path, $json + "`n", $utf8NoBom)
}

$result = [ordered]@{
    attempt_label = 'GEOMETRY_FILTER_FIXTURE_ONLY'
    status = 'PRECHECK'
    generator = $null
    julia = $null
    protected_pin_status = $null
    tree_integrity_status = $null
    artifact_bytes = $null
    artifact_limit_bytes = 20971520
}

$nativePinRows = Read-TablePins $nativePins
if ($nativePinRows.Count -ne 23) { throw "Expected 23 frozen overlay/source pins; found $($nativePinRows.Count)" }
$pilotPinRows = Read-TablePins $pilotPins '## Existing isolated'
if ($pilotPinRows.Count -ne 8) { throw "Expected 8 frozen pilot source/review pins; found $($pilotPinRows.Count)" }

$pinnedFiles = [System.Collections.Generic.List[object]]::new()
foreach ($pin in $nativePinRows) {
    $path = Resolve-NativePinPath $pin.relative_path
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Pinned source is missing: $path" }
    $actual = Get-FileSha256 $path
    if ($actual -cne $pin.expected_sha256) { throw "Frozen native source pin mismatch: $($pin.relative_path)" }
    $pinnedFiles.Add([pscustomobject]@{ path = $path; expected_sha256 = $pin.expected_sha256; actual_sha256 = $actual })
}
foreach ($pin in $pilotPinRows) {
    $path = Resolve-PilotPinPath $pin.relative_path
    if (-not (Test-Path -LiteralPath $path -PathType Leaf)) { throw "Pilot input is missing: $path" }
    $actual = Get-FileSha256 $path
    if ($actual -cne $pin.expected_sha256) { throw "Frozen pilot source pin mismatch: $($pin.relative_path)" }
    $pinnedFiles.Add([pscustomobject]@{ path = $path; expected_sha256 = $pin.expected_sha256; actual_sha256 = $actual })
}
$expectedNativeManifest = 'fd9f1278577ee8266d1a9767c6f4cedcc17dff716c2debbb69ed7921c844466f'
if ((Get-FileSha256 $nativePins) -cne $expectedNativeManifest) { throw '23-entry native pin manifest changed' }

$protectedSingleFiles = @(
    [pscustomobject]@{ path = $appProject; expected = 'eebff26de1f99a56b29a43ad767d2877bba5b6a0e4a03de1f909eab00f59f4ae' },
    [pscustomobject]@{ path = $appManifest; expected = 'eb7d8902b84cfc64faaa0d45b0e6a3ae8534a252d06399f3db074da5a738e594' },
    [pscustomobject]@{ path = $julia; expected = '29a5cb5bde6a3f4a151f6e714591ca9344638303d9ed2a893ab31ada3de5d103' },
    [pscustomobject]@{ path = $python; expected = (Get-FileSha256 $python) }
)
foreach ($item in $protectedSingleFiles) {
    if (-not (Test-Path -LiteralPath $item.path -PathType Leaf)) { throw "Protected file is missing: $($item.path)" }
    $actual = Get-FileSha256 $item.path
    if ($actual -cne $item.expected) { throw "Protected file hash mismatch: $($item.path)" }
}

$treeSpecs = @(
    [pscustomobject]@{ label = 'original_vendor'; path = $vendorRoot },
    [pscustomobject]@{ label = 'shared_depot'; path = $sharedDepot },
    [pscustomobject]@{ label = 'run03_project'; path = (Join-Path $run03 'project') },
    [pscustomobject]@{ label = 'run03_depot'; path = (Join-Path $run03 'depot') }
)
$beforeTrees = foreach ($spec in $treeSpecs) { Get-TreeSnapshot $spec.label $spec.path 'before' }
$beforePinnedPath = Join-Path $auditDir 'before_pinned_files.tsv'
$pinRecords = foreach ($item in $pinnedFiles) { "$( $item.path )`t$( $item.actual_sha256 )" }
[System.IO.File]::WriteAllLines($beforePinnedPath, $pinRecords, $utf8NoBom)
$beforeProtectedFiles = foreach ($item in $protectedSingleFiles) {
    [pscustomobject]@{ path = $item.path; sha256 = Get-FileSha256 $item.path }
}
Write-Json (Join-Path $auditDir 'before_snapshot.json') ([ordered]@{
    captured_utc = [DateTimeOffset]::UtcNow.ToString('o')
    overlay_source_pin_count = $nativePinRows.Count
    overlay_source_pins = $pinnedFiles | Select-Object -First $nativePinRows.Count
    pilot_source_pins = $pinnedFiles | Select-Object -Skip $nativePinRows.Count
    protected_single_files = $beforeProtectedFiles
    protected_trees = $beforeTrees
    protected_pin_manifest_sha256 = Get-FileSha256 $nativePins
})

$result.status = 'BASELINE_CAPTURED'
$generatorArgs = @('-B', $generator, $bundleRoot)
$generatorProcess = Start-CapturedProcess 'fixture_generator' $python $generatorArgs $workspace @{} 10000
$result.generator = [ordered]@{
    argv = $generatorProcess.argv
    executable_sha256 = Get-FileSha256 $python
    working_directory = $generatorProcess.working_directory
    started_utc = $generatorProcess.started_utc
    ended_utc = $generatorProcess.ended_utc
    duration_ms = $generatorProcess.duration_ms
    timeout_ms = $generatorProcess.timeout_ms
    timed_out = $generatorProcess.timed_out
    process_id = $generatorProcess.process_id
    child_exit_code = $generatorProcess.exit_code
    stdout_path = Join-Path $auditDir 'generator_stdout.txt'
    stderr_path = Join-Path $auditDir 'generator_stderr.txt'
}
[System.IO.File]::WriteAllText($result.generator.stdout_path, $generatorProcess.stdout, $utf8NoBom)
[System.IO.File]::WriteAllText($result.generator.stderr_path, $generatorProcess.stderr, $utf8NoBom)
$result.generator.stdout_sha256 = Get-FileSha256 $result.generator.stdout_path
$result.generator.stderr_sha256 = Get-FileSha256 $result.generator.stderr_path

$positiveBundle = Join-Path $bundleRoot 'positive_mask'
$blockedBundle = Join-Path $bundleRoot 'zero_degree_block'
if ($generatorProcess.timed_out -or $generatorProcess.exit_code -ne 0) {
    $result.status = if ($generatorProcess.timed_out) { 'GENERATOR_TIMEOUT_NO_RETRY' } else { 'GENERATOR_FAILED_NO_RETRY' }
}
else {
    $generatorReceipt = $generatorProcess.stdout | ConvertFrom-Json
    $fixtureReceiptsPath = Join-Path $auditDir 'generator_external_pins_receipts.json'
    Write-Json $fixtureReceiptsPath $generatorReceipt
    $positiveJson = Join-Path $positiveBundle 'audit_contract.json'
    $positiveToml = Join-Path $positiveBundle 'native_mask.toml'
    $blockedJson = Join-Path $blockedBundle 'audit_contract.json'
    $blockedToml = Join-Path $blockedBundle 'native_mask.toml'
    $goodJsonSha = [string]$generatorReceipt.positive_mask.audit_contract_json_sha256
    $goodTomlSha = [string]$generatorReceipt.positive_mask.native_toml_sha256
    $blockedJsonSha = [string]$generatorReceipt.zero_degree_block.audit_contract_json_sha256
    $blockedTomlSha = [string]$generatorReceipt.zero_degree_block.native_toml_sha256
    $fixtureBindings = @(
        [pscustomobject]@{ label='positive JSON'; path=$positiveJson; expected=$goodJsonSha },
        [pscustomobject]@{ label='positive TOML'; path=$positiveToml; expected=$goodTomlSha },
        [pscustomobject]@{ label='blocked JSON'; path=$blockedJson; expected=$blockedJsonSha },
        [pscustomobject]@{ label='blocked TOML'; path=$blockedToml; expected=$blockedTomlSha }
    )
    foreach ($binding in $fixtureBindings) {
        $actual = Get-FileSha256 $binding.path
        if ($actual -cne $binding.expected) { throw "Fixture external receipt did not bind actual bytes: $($binding.label)" }
    }
    $result.fixture_external_pins = $fixtureBindings | ForEach-Object {
        [pscustomobject]@{ label=$_.label; path=$_.path; expected_sha256=$_.expected; actual_sha256=(Get-FileSha256 $_.path) }
    }
    $nativeGraphReceipt = Join-Path $positiveBundle 'native_graph_readback.txt'
    $juliaArgs = @(
        '--startup-file=no', '--history-file=no', '--threads=1', '--compiled-modules=no',
        "--project=$(Join-Path $run03 'project')", $driver,
        $positiveJson, $positiveToml, $goodJsonSha, $goodTomlSha,
        $blockedJson, $blockedToml, $blockedJsonSha, $blockedTomlSha,
        $nativeGraphReceipt
    )
    $juliaEnv = @{
        JULIA_DEPOT_PATH = "$freshDepot;$sharedDepot"
        JULIA_PKG_OFFLINE = 'true'
        JULIA_PKG_PRECOMPILE_AUTO = '0'
    }
    $juliaProcess = Start-CapturedProcess 'native_body_pilot' $julia $juliaArgs $workspace $juliaEnv 60000
    $result.julia = [ordered]@{
        argv = $juliaProcess.argv
        executable_sha256 = Get-FileSha256 $julia
        working_directory = $juliaProcess.working_directory
        environment_overrides = $juliaProcess.environment_overrides
        started_utc = $juliaProcess.started_utc
        ended_utc = $juliaProcess.ended_utc
        duration_ms = $juliaProcess.duration_ms
        timeout_ms = $juliaProcess.timeout_ms
        timed_out = $juliaProcess.timed_out
        process_id = $juliaProcess.process_id
        child_exit_code = $juliaProcess.exit_code
        stdout_path = Join-Path $auditDir 'julia_stdout.txt'
        stderr_path = Join-Path $auditDir 'julia_stderr.txt'
        graph_receipt_path = $nativeGraphReceipt
    }
    [System.IO.File]::WriteAllText($result.julia.stdout_path, $juliaProcess.stdout, $utf8NoBom)
    [System.IO.File]::WriteAllText($result.julia.stderr_path, $juliaProcess.stderr, $utf8NoBom)
    $result.julia.stdout_sha256 = Get-FileSha256 $result.julia.stdout_path
    $result.julia.stderr_sha256 = Get-FileSha256 $result.julia.stderr_path
    if ($juliaProcess.timed_out) {
        $result.status = 'JULIA_TIMEOUT_NO_RETRY'
    }
    elseif ($juliaProcess.exit_code -ne 0) {
        $result.status = 'JULIA_FAILED_NO_RETRY'
    }
    elseif (-not (Test-Path -LiteralPath $nativeGraphReceipt -PathType Leaf)) {
        $result.status = 'JULIA_EXITED_ZERO_WITHOUT_GRAPH_RECEIPT_NO_RETRY'
    }
    else {
        $result.status = 'NATIVE_BODY_GRAPH_READBACK_PASS_FIXTURE_ONLY'
    }
}

$afterTrees = foreach ($spec in $treeSpecs) { Get-TreeSnapshot $spec.label $spec.path 'after' }
$treeComparisons = foreach ($before in $beforeTrees) {
    $after = $afterTrees | Where-Object label -EQ $before.label | Select-Object -First 1
    [pscustomobject]@{
        label = $before.label
        before_sha256 = $before.manifest_sha256
        after_sha256 = $after.manifest_sha256
        equal = ($before.manifest_sha256 -ceq $after.manifest_sha256)
        before_file_count = $before.file_count
        after_file_count = $after.file_count
        before_total_bytes = $before.total_bytes
        after_total_bytes = $after.total_bytes
    }
}
$afterPinnedPath = Join-Path $auditDir 'after_pinned_files.tsv'
$afterPinRecords = foreach ($item in $pinnedFiles) { "$( $item.path )`t$( Get-FileSha256 $item.path )" }
[System.IO.File]::WriteAllLines($afterPinnedPath, $afterPinRecords, $utf8NoBom)
$pinnedFilesEqual = (Get-FileSha256 $beforePinnedPath) -ceq (Get-FileSha256 $afterPinnedPath)
$afterProtectedFiles = foreach ($item in $protectedSingleFiles) {
    [pscustomobject]@{ path = $item.path; sha256 = Get-FileSha256 $item.path; equal_to_expected = ((Get-FileSha256 $item.path) -ceq $item.expected) }
}
$result.protected_pin_status = if ($pinnedFilesEqual -and ($afterProtectedFiles | Where-Object { -not $_.equal_to_expected }).Count -eq 0) { 'ALL_PINNED_FILES_UNCHANGED' } else { 'MISMATCH' }
$result.tree_integrity_status = if (($treeComparisons | Where-Object { -not $_.equal }).Count -eq 0) { 'ALL_PROTECTED_TREES_UNCHANGED' } else { 'MISMATCH' }
$result.post_run_integrity = [ordered]@{
    captured_utc = [DateTimeOffset]::UtcNow.ToString('o')
    pinned_files_equal = $pinnedFilesEqual
    protected_single_files = $afterProtectedFiles
    protected_tree_comparisons = $treeComparisons
    fresh_scratch_depot_path = $freshDepot
    fresh_scratch_depot_file_count = @(Get-ChildItem -LiteralPath $freshDepot -File -Recurse -Force).Count
}
$result.execution_status = $result.status
$result.completed_utc = [DateTimeOffset]::UtcNow.ToString('o')
$result.audit_notes = @(
    'Fixture-only two-box policy; external JSON and TOML receipt hashes were checked against actual bytes before Julia.',
    'Native body and BondSystem only; no BodyChunk, ThreadsBodyDataHandler, solver, force evaluation, CAD, or production job.',
    'The 20 MiB limit is total scratch artifact bytes, not process RSS.',
    'Any generator or Julia failure is preserved; this runner performs no retry.'
)
Write-Json (Join-Path $auditDir 'after_snapshot.json') ([ordered]@{
    captured_utc = $result.completed_utc
    protected_pin_status = $result.protected_pin_status
    tree_integrity_status = $result.tree_integrity_status
    protected_single_files = $afterProtectedFiles
    protected_tree_comparisons = $treeComparisons
})
$result.runner_script_path = Join-Path $scratchRoot 'run_once.ps1'
$result.runner_script_sha256 = Get-FileSha256 $result.runner_script_path
$reportPath = Join-Path $auditDir 'pilot_execution_result.json'
$result.artifact_bytes = 0
$result.artifact_file_count = 0
for ($iteration = 0; $iteration -lt 6; $iteration++) {
    if ($result.protected_pin_status -ne 'ALL_PINNED_FILES_UNCHANGED' -or
        $result.tree_integrity_status -ne 'ALL_PROTECTED_TREES_UNCHANGED') {
        $result.status = 'POSTRUN_PROTECTED_INTEGRITY_FAILURE'
    }
    else {
        $result.status = $result.execution_status
    }
    Write-Json $reportPath $result
    $allScratchFiles = @(Get-ChildItem -LiteralPath $scratchRoot -File -Recurse -Force)
    $artifactBytes = [long](($allScratchFiles | Measure-Object -Property Length -Sum).Sum)
    $artifactFileCount = $allScratchFiles.Count
    $result.artifact_bytes = $artifactBytes
    $result.artifact_file_count = $artifactFileCount
    $result.artifact_cap_status = if ($artifactBytes -le 20MB) { 'WITHIN_20_MIB' } else { 'OVER_20_MIB' }
    if ($result.artifact_cap_status -ne 'WITHIN_20_MIB' -and
        $result.protected_pin_status -eq 'ALL_PINNED_FILES_UNCHANGED' -and
        $result.tree_integrity_status -eq 'ALL_PROTECTED_TREES_UNCHANGED') {
        $result.status = 'ARTIFACT_CAP_EXCEEDED'
    }
}
Write-Json $reportPath $result
$allScratchFiles = @(Get-ChildItem -LiteralPath $scratchRoot -File -Recurse -Force)
$artifactBytes = [long](($allScratchFiles | Measure-Object -Property Length -Sum).Sum)
$result.artifact_bytes = $artifactBytes
$result.artifact_file_count = $allScratchFiles.Count
$result.artifact_cap_status = if ($artifactBytes -le 20MB) { 'WITHIN_20_MIB' } else { 'OVER_20_MIB' }
Write-Json $reportPath $result

Write-Output "PILOT_STATUS=$($result.status)"
Write-Output "PROTECTED_PINS=$($result.protected_pin_status)"
Write-Output "PROTECTED_TREES=$($result.tree_integrity_status)"
Write-Output "ARTIFACT_BYTES=$artifactBytes"
if ($result.generator) { Write-Output "GENERATOR_EXIT=$($result.generator.child_exit_code);TIMED_OUT=$($result.generator.timed_out);MS=$($result.generator.duration_ms)" }
if ($result.julia) { Write-Output "JULIA_EXIT=$($result.julia.child_exit_code);TIMED_OUT=$($result.julia.timed_out);MS=$($result.julia.duration_ms)" }
