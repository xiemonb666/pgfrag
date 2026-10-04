param(
    [Parameter(Mandatory = $true)]
    [switch]$ConfirmedNativeProcessExited
)
$ErrorActionPreference = 'Stop'
$source = '\\wsl.localhost\Ubuntu\home\xiemo\cohesive_bilinear_patch_20261004\controlled_intrinsic_bilinear_modeI_201_run03'
$destination = 'C:\Users\xiemo\Desktop\ShardForge\relicweave-synlab\tmp\review_20261004\implementation\cohesive_bilinear_patch_20261004\controlled_intrinsic_bilinear_modeI_201_run03'
if (-not $ConfirmedNativeProcessExited) { throw 'Confirm the native WSL process has exited before copying closed outputs.' }
if (-not (Test-Path -LiteralPath $source -PathType Container)) { throw "Native run output missing: $source" }
if (Test-Path -LiteralPath $destination) { throw "Review destination already exists; refusing overwrite: $destination" }
if (-not (Test-Path -LiteralPath (Join-Path $source 'run_started.json') -PathType Leaf)) { throw 'Runner did not create its own run_started.json; no output is copied.' }
Copy-Item -LiteralPath $source -Destination $destination -Recurse
if (-not (Test-Path -LiteralPath (Join-Path $destination 'run_started.json') -PathType Leaf)) { throw 'Copied run lacks runner-created run_started.json.' }
$files = Get-ChildItem -LiteralPath $destination -File -Recurse -Force
$inventory = foreach ($file in $files) {
    $relative = $file.FullName.Substring($destination.Length).TrimStart('\')
    [pscustomobject]@{ path = $relative; bytes = $file.Length; sha256 = (Get-FileHash -LiteralPath $file.FullName -Algorithm SHA256).Hash.ToLowerInvariant() }
}
$inventory | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath (Join-Path $destination 'copy_inventory.json') -Encoding UTF8