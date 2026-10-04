$ErrorActionPreference = 'Stop'
$InputDirectory = Join-Path $PSScriptRoot 'input'
New-Item -ItemType Directory -Force -Path $InputDirectory | Out-Null

$CellCount = 5
[int]$JitterPatternSeed = 17
[double]$CellSizeM = 0.001
[double]$CellVolumeM3 = 1.0e-9
$GeneratorHash = (Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash.ToLowerInvariant()
$Utf8NoBom = New-Object System.Text.UTF8Encoding($false)

function Format-Number([double]$Value) {
    return $Value.ToString('R', [System.Globalization.CultureInfo]::InvariantCulture)
}

foreach ($Layout in @('uniform_voxel', 'cell_jitter_candidate')) {
    $BaseName = if ($Layout -eq 'uniform_voxel') {
        'coupon_uniform_voxel_5cube_125'
    } else {
        'coupon_cell_jitter_5cube_seed17'
    }
    $CsvPath = Join-Path $InputDirectory ($BaseName + '.points.csv')
    $MetadataPath = Join-Path $InputDirectory ($BaseName + '.metadata.toml')
    $Rows = New-Object 'System.Collections.Generic.List[string]'
    $Rows.Add('point_id,x_m,y_m,z_m,volume_m3')
    [double]$WeightedX = 0.0
    [double]$WeightedY = 0.0
    [double]$WeightedZ = 0.0
    $PointId = 0

    for ($k = 0; $k -lt $CellCount; $k++) {
        for ($j = 0; $j -lt $CellCount; $j++) {
            for ($i = 0; $i -lt $CellCount; $i++) {
                $PointId++
                [double]$X = ($i + 0.5) * $CellSizeM
                [double]$Y = ($j + 0.5) * $CellSizeM
                [double]$Z = ($k + 0.5) * $CellSizeM
                if ($Layout -eq 'cell_jitter_candidate') {
                    # Deterministic offsets in [-0.2, 0.2] of each original cell width.
                    $X += ((($i * $JitterPatternSeed + $j * 3 + $k * 4) % 5) - 2) * 0.1 * $CellSizeM
                    $Y += ((($i * 4 + $j * $JitterPatternSeed + $k * 2) % 5) - 2) * 0.1 * $CellSizeM
                    $Z += ((($i * 3 + $j * 2 + $k * $JitterPatternSeed) % 5) - 2) * 0.1 * $CellSizeM
                }
                $Rows.Add(($PointId.ToString([System.Globalization.CultureInfo]::InvariantCulture) + ',' +
                           (Format-Number $X) + ',' + (Format-Number $Y) + ',' +
                           (Format-Number $Z) + ',' + (Format-Number $CellVolumeM3)))
                $WeightedX += $X * $CellVolumeM3
                $WeightedY += $Y * $CellVolumeM3
                $WeightedZ += $Z * $CellVolumeM3
            }
        }
    }

    $CsvText = ($Rows -join "`n") + "`n"
    [System.IO.File]::WriteAllText($CsvPath, $CsvText, $Utf8NoBom)
    $CsvHash = (Get-FileHash -LiteralPath $CsvPath -Algorithm SHA256).Hash.ToLowerInvariant()
    $VolumeSum = $CellCount * $CellCount * $CellCount * $CellVolumeM3
    $ComX = $WeightedX / $VolumeSum
    $ComY = $WeightedY / $VolumeSum
    $ComZ = $WeightedZ / $VolumeSum
    $BottomIds = (1..($CellCount * $CellCount)) -join ', '
    $TopStart = ($CellCount - 1) * $CellCount * $CellCount + 1
    $TopEnd = $CellCount * $CellCount * $CellCount
    $TopIds = ($TopStart..$TopEnd) -join ', '

    if ($Layout -eq 'uniform_voxel') {
        $InputId = 'coupon_uniform_voxel_5cube_125_v1'
        $SourceDescription = 'Procedurally generated cell-center material points for a 5 mm cube; one positive 1 mm^3 voxel quadrature weight per cell.'
        $QuadratureMethod = 'uniform_voxel_cell_center'
        $VolumeSemantics = 'Known 1 mm^3 cell volume assigned to each cell-center quadrature point.'
    } else {
        $InputId = 'coupon_cell_jitter_5cube_seed17_v1'
        $SourceDescription = 'Procedurally generated 5 mm cube voxel quadrature candidate; fixed integer pattern seed 17 jitters each sample inside its original cell by at most 0.2 cell widths while retaining that cell volume.'
        $QuadratureMethod = 'uniform_voxel_quadrature_with_deterministic_in_cell_jitter_candidate'
        $VolumeSemantics = 'The original 1 mm^3 voxel cell quadrature weight is retained after positional jitter; this is not an exact arbitrary Voronoi volume.'
    }

    $Metadata = @"
schema = "pgfrag.pd.pointset.v1"
input_id = "$InputId"
point_data_file = "$([System.IO.Path]::GetFileName($CsvPath))"
point_data_sha256 = "$CsvHash"
coordinate_units = "m"
volume_units = "m^3"
coordinate_frame = "Right-handed Cartesian xyz in metres; origin is the lower corner of the declared 5 mm coupon."
point_id_semantics = "Positive source point identity; CSV row order defines native Body index and is preserved separately in outputs."
point_count = 125
declared_volume_sum_m3 = $(Format-Number $VolumeSum)
declared_com_m = [$(Format-Number $ComX), $(Format-Number $ComY), $(Format-Number $ComZ)]

[source]
kind = "declared_internal_fixture_generation"
description = "$SourceDescription"
reference = "make_point_inputs.ps1 in this experiment directory; not a measured specimen or imported shell."
source_sha256 = "$GeneratorHash"
source_artifact_sha256 = "NOT_APPLICABLE_GENERATED_FIXTURE"

[quadrature]
method = "$QuadratureMethod"
volume_semantics = "$VolumeSemantics"
voronoi_claim = false
real_shell_claim = false

[load_sets]
axis = 3
bottom_point_ids = [$BottomIds]
top_point_ids = [$TopIds]
"@
    [System.IO.File]::WriteAllText($MetadataPath, $Metadata, $Utf8NoBom)
}
