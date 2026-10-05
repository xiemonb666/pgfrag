# Source-only one-shot launcher. Run only after the root agent authorizes execution.
$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$ExperimentDir = [System.IO.Path]::GetFullPath($PSScriptRoot).TrimEnd('\')
$RepoRoot = [System.IO.Path]::GetFullPath((Join-Path $ExperimentDir '..\..\..')).TrimEnd('\')
$JuliaPath = Join-Path $RepoRoot 'third_party\pgfrag-runtime\extract\julia-1.10.12\bin\julia.exe'
$JuliaProject = Join-Path $RepoRoot 'pgfrag\src\pgfrag\resources\julia'
$DriverPath = Join-Path $ExperimentDir 'run_free_fracture_pilot.jl'
$CasePath = Join-Path $ExperimentDir 'fixed_case.toml'
$FreezePath = Join-Path $ExperimentDir 'source_freeze.toml'
$AttemptDir = Join-Path $ExperimentDir 'run_attempt_01'
$OutputDir = Join-Path $AttemptDir 'results'
$summaryPath = Join-Path $OutputDir 'summary.toml'
$ReceiptPath = Join-Path $AttemptDir 'launcher_receipt.json'
$StdoutPath = Join-Path $AttemptDir 'native_stdout.log'
$StderrPath = Join-Path $AttemptDir 'native_stderr.log'
$ClaimPath = Join-Path $ExperimentDir 'run_attempt_01.claim'
$MaximumAttemptBytes = [int64]0
$MaximumRunMilliseconds = [int]0

function Get-Sha256([string]$Path) {
    return (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
}

function Get-RelativeRepoPath([string]$Path) {
    $fullPath = [System.IO.Path]::GetFullPath($Path)
    $prefix = $RepoRoot + '\'
    if (-not $fullPath.StartsWith($prefix, [System.StringComparison]::OrdinalIgnoreCase)) {
        throw "Path is outside the pinned repository: $fullPath"
    }
    return $fullPath.Substring($prefix.Length).Replace('\', '/')
}

function Get-PinnedSha256([string]$RelativePath, [string]$FreezeText) {
    $pinBlocks = [regex]::Split($FreezeText, '(?m)^\s*\[\[pins\]\]\s*$')
    foreach ($pinBlock in $pinBlocks) {
        $pathMatch = [regex]::Match($pinBlock, '(?m)^\s*path\s*=\s*"([^"]+)"\s*$')
        if ($pathMatch.Success -and $pathMatch.Groups[1].Value -ceq $RelativePath) {
            $shaMatch = [regex]::Match($pinBlock, '(?m)^\s*sha256\s*=\s*"([0-9a-fA-F]{64})"\s*$')
            if (-not $shaMatch.Success) { throw "Frozen pin has no valid SHA-256: $RelativePath" }
            return $shaMatch.Groups[1].Value.ToLowerInvariant()
        }
    }
    throw "Source freeze does not include required pin: $RelativePath"
}

function Get-TreeBytes([string]$Directory) {
    if (-not (Test-Path -LiteralPath $Directory -PathType Container)) { return [int64]0 }
    $totalBytes = [int64]0
    foreach ($filePath in [System.IO.Directory]::EnumerateFiles(
            $Directory, '*', [System.IO.SearchOption]::AllDirectories)) {
        $totalBytes += [System.IO.FileInfo]::new($filePath).Length
    }
    return $totalBytes
}

function Get-CaseInteger([string]$CaseText, [string]$Key) {
    $match = [regex]::Match($CaseText, "(?m)^\s*$([regex]::Escape($Key))\s*=\s*(\d+)\s*$")
    if (-not $match.Success) { throw "Fixed case lacks integer field $Key" }
    return [int64]$match.Groups[1].Value
}

foreach ($requiredFile in @($JuliaPath, $DriverPath, $CasePath, $FreezePath)) {
    if (-not (Test-Path -LiteralPath $requiredFile -PathType Leaf)) {
        throw "Required pinned file is missing: $requiredFile"
    }
}
if (-not (Test-Path -LiteralPath $JuliaProject -PathType Container)) {
    throw "Pinned Julia project is missing: $JuliaProject"
}
if (Test-Path -LiteralPath $AttemptDir) { throw "One-shot output already exists: $AttemptDir" }
if (Test-Path -LiteralPath $ClaimPath) { throw "One-shot claim already exists: $ClaimPath" }

$CaseText = [System.IO.File]::ReadAllText($CasePath)
if ($CaseText -notmatch '(?m)^run_status\s*=\s*"SOURCE_ONLY_NOT_RUN"\s*$') {
    throw 'Fixed case is no longer marked SOURCE_ONLY_NOT_RUN.'
}
$ExpectedPointCount = Get-CaseInteger $CaseText 'expected_point_count'
$MaximumPointCount = Get-CaseInteger $CaseText 'max_points'
$MaximumStepCount = Get-CaseInteger $CaseText 'max_steps'
$MaximumOutputBytes = Get-CaseInteger $CaseText 'max_output_bytes'
$MaximumAttemptBytes = Get-CaseInteger $CaseText 'max_attempt_bytes'
$MaximumRunMilliseconds = Get-CaseInteger $CaseText 'max_runtime_ms'
if ($ExpectedPointCount -gt $MaximumPointCount -or $MaximumPointCount -gt 5000) {
    throw 'Fixed case exceeds its point cap.'
}
if ($MaximumStepCount -gt 180 -or $MaximumStepCount -le 0) {
    throw 'Fixed case exceeds its 180-step source cap.'
}
if ($MaximumOutputBytes -gt (190 * 1024 * 1024) -or $MaximumOutputBytes -le 0) {
    throw 'Fixed case exceeds its 190 MiB experiment-output budget.'
}
$MaximumBinaryBytes = [int64](21 * $ExpectedPointCount * ($MaximumStepCount + 1) * 8)
$MaximumGraphLabelBytes = [int64](4 * $ExpectedPointCount * ($MaximumStepCount + 1))
if ($MaximumBinaryBytes + $MaximumGraphLabelBytes + (48 * 1024 * 1024) -gt $MaximumOutputBytes) {
    throw 'Maximum point-state and graph-label tensors plus bond/event reserve exceed the fixed-case budget.'
}
if ($MaximumAttemptBytes -gt (200 * 1024 * 1024) -or
    $MaximumAttemptBytes -lt $MaximumOutputBytes -or $MaximumAttemptBytes -le 0) {
    throw 'Fixed case exceeds its 200 MiB total-attempt budget.'
}
if ($MaximumRunMilliseconds -gt 180000 -or $MaximumRunMilliseconds -le 0) {
    throw 'Fixed case exceeds its 180-second runtime cap.'
}

$FreezeText = [System.IO.File]::ReadAllText($FreezePath)
$FreezeSha256 = Get-Sha256 $FreezePath
$PinnedInputPaths = @($DriverPath, $CasePath, $JuliaPath, $PSCommandPath)
$PinnedInputSha256 = @{}
foreach ($inputPath in $PinnedInputPaths) {
    $relativePath = Get-RelativeRepoPath $inputPath
    $actualSha256 = Get-Sha256 $inputPath
    $expectedSha256 = Get-PinnedSha256 $relativePath $FreezeText
    if ($actualSha256 -cne $expectedSha256) {
        throw "Frozen input changed: $relativePath"
    }
    $PinnedInputSha256[$relativePath] = $actualSha256
}

$driveRoot = [System.IO.Path]::GetPathRoot($ExperimentDir)
$driveInfo = [System.IO.DriveInfo]::new($driveRoot)
if ($driveInfo.AvailableFreeSpace -lt ($MaximumAttemptBytes + 32 * 1024 * 1024)) {
    throw 'Insufficient free disk space for the fixed 200 MiB attempt cap and safety margin.'
}

$claimStream = [System.IO.File]::Open(
    $ClaimPath, [System.IO.FileMode]::CreateNew, [System.IO.FileAccess]::Write,
    [System.IO.FileShare]::None)
try {
    $claimBytes = [System.Text.Encoding]::UTF8.GetBytes(
        "one-shot claim created $([DateTimeOffset]::UtcNow.ToString('o'))`n")
    $claimStream.Write($claimBytes, 0, $claimBytes.Length)
} finally {
    $claimStream.Dispose()
}
[System.IO.Directory]::CreateDirectory($AttemptDir) | Out-Null

$nativeProcess = $null
$nativeProcessStarted = $false
$nativeExitCode = $null
$nativeProcessId = $null
$nativeProcessKilled = $false
$timedOut = $false
$diskCapExceeded = $false
$outputCapExceeded = $false
$launcherError = $null
$summaryDeclaredStatus = $null
$status = 'START_FAILED'
$startedUtc = [DateTimeOffset]::UtcNow
$endedUtc = $null
$stopwatch = $null
$durationMilliseconds = -1
$attemptBytesAtStop = [int64]0
$experimentOutputBytes = [int64]0
$arguments = @(
    '--startup-file=no',
    '--compiled-modules=no',
    '--threads=1',
    "--project=$JuliaProject",
    $DriverPath,
    $CasePath,
    $OutputDir)
foreach ($argument in $arguments) {
    if ($argument -match '\s') { throw "This pinned launcher requires whitespace-free arguments: $argument" }
}

$environmentNames = @('JULIA_NUM_THREADS', 'JULIA_DEPOT_PATH', 'JULIA_PKG_OFFLINE')
$previousEnvironment = @{}
foreach ($name in $environmentNames) {
    $previousEnvironment[$name] = [System.Environment]::GetEnvironmentVariable($name, 'Process')
}
try {
    [System.Environment]::SetEnvironmentVariable('JULIA_NUM_THREADS', '1', 'Process')
    [System.Environment]::SetEnvironmentVariable(
        'JULIA_DEPOT_PATH', (Join-Path $RepoRoot 'third_party\pgfrag-runtime\depot'), 'Process')
    [System.Environment]::SetEnvironmentVariable('JULIA_PKG_OFFLINE', '1', 'Process')
    $nativeProcess = Start-Process -FilePath $JuliaPath -ArgumentList $arguments -PassThru `
        -WindowStyle Hidden -WorkingDirectory $RepoRoot `
        -RedirectStandardOutput $StdoutPath -RedirectStandardError $StderrPath
    $nativeProcessStarted = $true
    $nativeProcessId = [int]$nativeProcess.Id
    $stopwatch = [System.Diagnostics.Stopwatch]::StartNew()
    while (-not $nativeProcess.HasExited) {
        $remainingMs = $MaximumRunMilliseconds - [int]$stopwatch.ElapsedMilliseconds
        if ($remainingMs -le 0) {
            $timedOut = $true
            $nativeProcess.Kill()
            $nativeProcessKilled = $true
            $nativeProcess.WaitForExit()
            break
        }
        $waitMs = [int][Math]::Min(1000, $remainingMs)
        [void]$nativeProcess.WaitForExit($waitMs)
        $attemptBytesAtStop = Get-TreeBytes $AttemptDir
        if ($attemptBytesAtStop + 4096 -gt $MaximumAttemptBytes) {
            $diskCapExceeded = $true
            if (-not $nativeProcess.HasExited) {
                $nativeProcess.Kill()
                $nativeProcessKilled = $true
                $nativeProcess.WaitForExit()
            }
            break
        }
    }
    # A null ExitCode must remain unavailable; PowerShell's [int]$null becomes zero.
    $nativeProcess.WaitForExit()
    $nativeProcess.Refresh()
    if ($nativeProcess.HasExited) {
        $rawNativeExitCode = $nativeProcess.ExitCode
        if ($null -ne $rawNativeExitCode) {
            $nativeExitCode = [int]$rawNativeExitCode
        }
    }
    $stopwatch.Stop()
    $endedUtc = [DateTimeOffset]::UtcNow
    $attemptBytesAtStop = Get-TreeBytes $AttemptDir
    $experimentOutputBytes = Get-TreeBytes $OutputDir
    $outputCapExceeded = $experimentOutputBytes -gt $MaximumOutputBytes
    $summaryPath = Join-Path $OutputDir 'summary.toml'
    if (Test-Path -LiteralPath $summaryPath -PathType Leaf) {
        $summaryText = [System.IO.File]::ReadAllText($summaryPath)
        $summaryStatusMatch = [regex]::Match($summaryText, '(?m)^status\s*=\s*"([^"]+)"\s*$')
        if ($summaryStatusMatch.Success) { $summaryDeclaredStatus = $summaryStatusMatch.Groups[1].Value }
    }
    if ($null -eq $nativeExitCode) {
        $status = 'PROCESS_EXIT_CODE_UNAVAILABLE'
    } elseif ($timedOut) {
        $status = 'TIMED_OUT_KILLED'
    } elseif ($diskCapExceeded -or $outputCapExceeded) {
        $status = 'DISK_CAP_EXCEEDED'
    } elseif ($nativeExitCode -eq 0) {
        if ($summaryDeclaredStatus -ceq 'ENGINEERING_DIAGNOSTIC_COMPLETED') {
            $status = 'COMPLETED_EXIT_0'
        } else {
            $status = 'EXIT_0_SUMMARY_NOT_COMPLETED'
        }
    } else {
        $status = 'COMPLETED_NONZERO_EXIT'
    }
    $durationMilliseconds = [int64]$stopwatch.ElapsedMilliseconds
} catch {
    $launcherError = $_.Exception.Message
    $endedUtc = [DateTimeOffset]::UtcNow
    if ($nativeProcess) {
        try {
            if (-not $nativeProcess.HasExited) {
                try {
                    $nativeProcess.Kill()
                    $nativeProcessKilled = $true
                    $nativeProcess.WaitForExit()
                } catch { }
            }
            if ($nativeProcess.HasExited) {
                $nativeProcess.WaitForExit()
                $nativeProcess.Refresh()
                if ($nativeProcess.HasExited) {
                    $rawNativeExitCode = $nativeProcess.ExitCode
                    if ($null -ne $rawNativeExitCode) {
                        $nativeExitCode = [int]$rawNativeExitCode
                    }
                }
            }
        } catch { }
    }
    if ($null -eq $nativeExitCode) {
        $status = 'PROCESS_EXIT_CODE_UNAVAILABLE'
    } else {
        $status = 'LAUNCHER_ERROR'
    }
    $attemptBytesAtStop = Get-TreeBytes $AttemptDir
    $experimentOutputBytes = Get-TreeBytes $OutputDir
    $outputCapExceeded = $experimentOutputBytes -gt $MaximumOutputBytes
    if ($stopwatch) { $durationMilliseconds = [int64]$stopwatch.ElapsedMilliseconds }
} finally {
    foreach ($name in $environmentNames) {
        [System.Environment]::SetEnvironmentVariable($name, $previousEnvironment[$name], 'Process')
    }
}

$stdoutBytes = [int64]0
$stderrBytes = [int64]0
$stdoutSha256 = $null
$stderrSha256 = $null
if (Test-Path -LiteralPath $StdoutPath -PathType Leaf) {
    $stdoutBytes = [System.IO.FileInfo]::new($StdoutPath).Length
    $stdoutSha256 = Get-Sha256 $StdoutPath
}
if (Test-Path -LiteralPath $StderrPath -PathType Leaf) {
    $stderrBytes = [System.IO.FileInfo]::new($StderrPath).Length
    $stderrSha256 = Get-Sha256 $StderrPath
}
$summarySha256 = $null
if (Test-Path -LiteralPath $summaryPath -PathType Leaf) { $summarySha256 = Get-Sha256 $summaryPath }
$sourcePinsJson = ($PinnedInputSha256.GetEnumerator() | Sort-Object Name |
                   ForEach-Object { "$($_.Name)=$($_.Value)" }) -join ';'
$endedUtcText = $null
if ($endedUtc) { $endedUtcText = $endedUtc.ToString('o') }
$durationReceipt = $null
if ($durationMilliseconds -ge 0) { $durationReceipt = [int64]$durationMilliseconds }

# This receipt contains only primitive values. Raw native stdout/stderr remain in sibling files.
$receipt = [ordered]@{
    schema = 'pgfrag.free_fracture_pilot.launch_receipt.v1'
    status = $status
    process_started = [bool]$nativeProcessStarted
    process_id = $nativeProcessId
    native_process_exit_code = $nativeExitCode
    native_process_killed = [bool]$nativeProcessKilled
    timed_out = [bool]$timedOut
    disk_cap_exceeded = [bool]$diskCapExceeded
    output_cap_exceeded = [bool]$outputCapExceeded
    cap_exceeded = [bool]($diskCapExceeded -or $outputCapExceeded)
    max_runtime_ms = [int]$MaximumRunMilliseconds
    max_attempt_bytes = [int64]$MaximumAttemptBytes
    expected_point_count = [int64]$ExpectedPointCount
    max_point_count = [int64]$MaximumPointCount
    max_step_count = [int64]$MaximumStepCount
    max_output_bytes = [int64]$MaximumOutputBytes
    max_point_state_binary_bytes = [int64]$MaximumBinaryBytes
    max_graph_label_binary_bytes = [int64]$MaximumGraphLabelBytes
    started_utc = $startedUtc.ToString('o')
    ended_utc = $endedUtcText
    duration_ms = $durationReceipt
    julia_path = $JuliaPath
    julia_sha256 = $PinnedInputSha256[(Get-RelativeRepoPath $JuliaPath)]
    driver_sha256 = $PinnedInputSha256[(Get-RelativeRepoPath $DriverPath)]
    fixed_case_sha256 = $PinnedInputSha256[(Get-RelativeRepoPath $CasePath)]
    source_freeze_sha256 = $FreezeSha256
    required_source_pins = $sourcePinsJson
    output_directory = $OutputDir
    experiment_output_bytes = [int64]$experimentOutputBytes
    attempt_bytes_before_receipt = [int64]$attemptBytesAtStop
    stdout_path = $StdoutPath
    stdout_bytes = [int64]$stdoutBytes
    stdout_sha256 = $stdoutSha256
    stderr_path = $StderrPath
    stderr_bytes = [int64]$stderrBytes
    stderr_sha256 = $stderrSha256
    summary_sha256 = $summarySha256
    summary_status = $summaryDeclaredStatus
    result_authoritative = $false
    authority_note = 'Final launcher receipt process and budget fields govern; summary.toml is provisional until this receipt confirms COMPLETED_EXIT_0 and all caps remain clear.'
    launcher_error = $launcherError
}
$receiptJson = $receipt | ConvertTo-Json -Compress -Depth 2
[System.IO.File]::WriteAllText(
    $ReceiptPath, $receiptJson + "`n", [System.Text.UTF8Encoding]::new($false))
$attemptBytesFinal = Get-TreeBytes $AttemptDir
if ($attemptBytesFinal -gt $MaximumAttemptBytes) {
    $receipt.status = 'FINAL_ATTEMPT_DISK_CAP_EXCEEDED'
    $receipt.disk_cap_exceeded = $true
    $receipt.cap_exceeded = $true
    $receipt.attempt_bytes_before_receipt = [int64]$attemptBytesAtStop
    $receipt.final_attempt_bytes = [int64]$attemptBytesFinal
    $receipt.result_authoritative = $false
    $receiptJson = $receipt | ConvertTo-Json -Compress -Depth 2
    [System.IO.File]::WriteAllText(
        $ReceiptPath, $receiptJson + "`n", [System.Text.UTF8Encoding]::new($false))
    exit 1
}
$receipt.final_attempt_bytes = [int64]$attemptBytesFinal
$receipt.final_budget_ok = [bool](-not $receipt.cap_exceeded -and
    $experimentOutputBytes -le $MaximumOutputBytes -and
    $attemptBytesFinal -le $MaximumAttemptBytes)
$receipt.result_authoritative = [bool]($status -eq 'COMPLETED_EXIT_0' -and
    -not $receipt.cap_exceeded -and -not $receipt.timed_out -and
    -not $receipt.native_process_killed -and $receipt.native_process_exit_code -eq 0 -and
    $receipt.final_budget_ok -and
    $receipt.summary_status -ceq 'ENGINEERING_DIAGNOSTIC_COMPLETED')
$receiptJson = $receipt | ConvertTo-Json -Compress -Depth 2
[System.IO.File]::WriteAllText(
    $ReceiptPath, $receiptJson + "`n", [System.Text.UTF8Encoding]::new($false))
$attemptBytesAfterFinalReceipt = Get-TreeBytes $AttemptDir
if ($attemptBytesAfterFinalReceipt -gt $MaximumAttemptBytes) {
    $receipt.status = 'FINAL_ATTEMPT_DISK_CAP_EXCEEDED'
    $receipt.disk_cap_exceeded = $true
    $receipt.cap_exceeded = $true
    $receipt.final_budget_ok = $false
    $receipt.result_authoritative = $false
    $receipt.final_attempt_bytes = [int64]$attemptBytesAfterFinalReceipt
    $receiptJson = $receipt | ConvertTo-Json -Compress -Depth 2
    [System.IO.File]::WriteAllText(
        $ReceiptPath, $receiptJson + "`n", [System.Text.UTF8Encoding]::new($false))
    exit 1
}

if ($status -eq 'COMPLETED_EXIT_0') { exit 0 }
exit 1
