$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Drawing
$root = 'C:\Users\xiemo\Desktop\ShardForge'
$evidence = Join-Path $root 'relicweave-synlab\tmp\review_20261004\implementation\cohesive_capability'
$png = Join-Path $evidence 'cohesive_capability_preview.png'
$traces = @{
  dt = @(Import-Csv (Join-Path $evidence 'modeI_dt\trace.csv'))
  dt2 = @(Import-Csv (Join-Path $evidence 'modeI_dt2\trace.csv'))
}
$vis = Get-Content (Join-Path $evidence 'native_visualization_data.json') -Raw | ConvertFrom-Json
$bmp = [System.Drawing.Bitmap]::new(1125, 765)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.ScaleTransform(0.75,0.75)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAliasGridFit
$white = [System.Drawing.Brushes]::White
$g.Clear([System.Drawing.Color]::White)
$font = [System.Drawing.Font]::new('Arial', 10)
$small = [System.Drawing.Font]::new('Arial', 8)
$title = [System.Drawing.Font]::new('Arial', 19, [System.Drawing.FontStyle]::Bold)
$panelTitle = [System.Drawing.Font]::new('Arial', 12, [System.Drawing.FontStyle]::Bold)
$black = [System.Drawing.Brushes]::Black
$grayPen = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(210,215,220), 1)
$axisPen = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(70,75,80), 1.3)
$blue = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(40,110,190), 2)
$orange = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(225,115,35), 2)
$green = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(45,145,85), 2)
$red = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(190,45,55), 2)
$purple = [System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(130,70,165), 2)
$g.DrawString('Akantu 5.0.7.post1 · tiny dynamic cohesive capability evidence', $title, $black, 30, 20)
$g.DrawString('UNVALIDATED_TOY · SI · central 4 mm² API contract surface only · not a natural fracture / ceramic qualification', $font, $black, 32, 52)

function Draw-Panel {
  param($x,$y,$w,$h,$heading,$xlabel,$ylabel,$xmin,$xmax,$ymin,$ymax)
  $plot = [System.Drawing.RectangleF]::new($x+58,$y+42,$w-75,$h-87)
  $g.DrawString($heading,$panelTitle,$black,$x+4,$y+3)
  for($i=0;$i -le 4;$i++) {
    $yy=$plot.Y+$plot.Height*$i/4
    $g.DrawLine($grayPen,$plot.X,$yy,$plot.Right,$yy)
    $v=$ymax-($ymax-$ymin)*$i/4
    $g.DrawString(('{0:0.###}' -f $v),$small,$black,$x+2,$yy-7)
    $xx=$plot.X+$plot.Width*$i/4
    $g.DrawLine($grayPen,$xx,$plot.Y,$xx,$plot.Bottom)
    $vx=$xmin+($xmax-$xmin)*$i/4
    $g.DrawString(('{0:0.###}' -f $vx),$small,$black,$xx-15,$plot.Bottom+4)
  }
  $g.DrawRectangle($axisPen,$plot.X,$plot.Y,$plot.Width,$plot.Height)
  $g.DrawString($xlabel,$small,$black,$plot.X+$plot.Width/2-45,$y+$h-22)
  $g.DrawString($ylabel,$small,$black,$x+3,$y+25)
  return $plot
}
function Add-Series {
  param($plot,$xs,$ys,$xmin,$xmax,$ymin,$ymax,$pen)
  $pts = [System.Collections.Generic.List[System.Drawing.PointF]]::new()
  $n=[Math]::Min($xs.Count,$ys.Count)
  for($i=0;$i -lt $n;$i+=5) {
    $xx=$plot.X+($xs[$i]-$xmin)/($xmax-$xmin)*$plot.Width
    $yy=$plot.Bottom-($ys[$i]-$ymin)/($ymax-$ymin)*$plot.Height
    if($xx -ge $plot.X -and $xx -le $plot.Right -and $yy -ge $plot.Y -and $yy -le $plot.Bottom) {$pts.Add([System.Drawing.PointF]::new([single]$xx,[single]$yy))}
  }
  if($n -gt 0) {
    $i=$n-1; $xx=$plot.X+($xs[$i]-$xmin)/($xmax-$xmin)*$plot.Width; $yy=$plot.Bottom-($ys[$i]-$ymin)/($ymax-$ymin)*$plot.Height
    if($xx -ge $plot.X -and $xx -le $plot.Right -and $yy -ge $plot.Y -and $yy -le $plot.Bottom) {$pts.Add([System.Drawing.PointF]::new([single]$xx,[single]$yy))}
  }
  if($pts.Count -gt 1) {$g.DrawLines($pen,$pts.ToArray())}
}
function Legend {
  param($x,$y,$items)
  $xx=$x
  foreach($item in $items) { $g.DrawLine($item.Pen,$xx,$y+7,$xx+20,$y+7); $g.DrawString($item.Label,$small,$black,$xx+24,$y); $xx += $item.Width }
}

$disps=@(); $forces=@()
foreach($case in @('dt','dt2')) { foreach($r in $traces[$case]) { $disps += ([double]$r.imposed_top_displacement_m*1e6); $forces += [double]$r.actuator_force_estimate_N } }
$xmax=[Math]::Max(1,[Math]::Ceiling((($disps | Measure-Object -Maximum).Maximum)/5)*5)
$ymax=[Math]::Max(1,[Math]::Ceiling((($forces | Measure-Object -Maximum).Maximum)/5)*5)
$p1=Draw-Panel 35 90 700 390 'Boundary response (actual actuator estimator)' 'imposed top displacement (µm)' 'force (N)' 0 $xmax 0 $ymax
foreach($case in @('dt','dt2')) { $rows=$traces[$case]; $xs=@($rows | ForEach-Object {[double]$_.imposed_top_displacement_m*1e6}); $ys=@($rows | ForEach-Object {[double]$_.actuator_force_estimate_N}); Add-Series $p1 $xs $ys 0 $xmax 0 $ymax $(if($case -eq 'dt'){$blue}else{$orange}) }
Legend 115 423 @(@{Pen=$blue;Label='dt';Width=75},@{Pen=$orange;Label='dt/2';Width=90})

$p2=Draw-Panel 765 90 700 390 'Energy history (mJ; same end time)' 'time (ms)' 'energy (mJ)' 0 0.1 -0.1 0.7
foreach($case in @('dt','dt2')) {
  $rows=$traces[$case]; $xs=@($rows | ForEach-Object {[double]$_.time_s*1e3})
  foreach($series in @(@{Col='external_work_J';Pen=$(if($case -eq 'dt'){$blue}else{$orange})},@{Col='kinetic_J';Pen=$(if($case -eq 'dt'){$green}else{$purple})},@{Col='bulk_potential_J';Pen=$(if($case -eq 'dt'){$red}else{$grayPen})})) {
    $ys=@($rows | ForEach-Object {[double]$_.$($series.Col)*1e3}); Add-Series $p2 $xs $ys 0 0.1 -0.1 0.7 $series.Pen
  }
}
$g.DrawString('work / kinetic / bulk potential; color pairs = dt, dt/2',$small,$black,825,423)

$p3=Draw-Panel 35 510 700 390 'Native cohesive dissipated output (signed; threshold is only an energy proxy)' 'time (ms)' 'dissipated (mJ)' 0 0.1 -0.08 0.04
foreach($case in @('dt','dt2')) { $rows=$traces[$case]; $xs=@($rows | ForEach-Object {[double]$_.time_s*1e3}); $ys=@($rows | ForEach-Object {[double]$_.cohesive_dissipated_J*1e3}); Add-Series $p3 $xs $ys 0 0.1 -0.08 0.04 $(if($case -eq 'dt'){$blue}else{$orange}) }
$targetY=$p3.Bottom-(0.0304-(-0.08))/(0.04-(-0.08))*$p3.Height
$dashPen=[System.Drawing.Pen]::new([System.Drawing.Color]::DarkRed,1.5); $dashPen.DashStyle=[System.Drawing.Drawing2D.DashStyle]::Dash
$g.DrawLine($dashPen,$p3.X,$targetY,$p3.Right,$targetY)
$g.DrawString('Gc·A = +0.0304 mJ; native signed final values are negative; target not reached',$small,[System.Drawing.Brushes]::DarkRed,90,866)

$p4=Draw-Panel 765 510 700 390 'Actual native cohesive interface projection (mode-I dt, deformed x-y; full scale)' 'x (mm)' 'y (mm)' -1.1 1.1 -1.1 1.1
$caseData=$vis.modeI_dt
foreach($facet in $caseData.cohesive_facets) {
  foreach($side in @('side_a','side_b')) {
    $nodes=@($facet.$side.PSObject.Properties | ForEach-Object {$_.Value})
    if($nodes.Count -eq 3) {
      $pts=[System.Drawing.PointF[]]@($nodes | ForEach-Object {
        $xx=$p4.X+(([double]$_.deformed_xyz_m[0]*1000+1.1)/2.2)*$p4.Width
        $yy=$p4.Bottom-(([double]$_.deformed_xyz_m[1]*1000+1.1)/2.2)*$p4.Height
        [System.Drawing.PointF]::new([single]$xx,[single]$yy)
      })
      $g.DrawPolygon($(if($side -eq 'side_a'){$blue}else{$orange}),$pts)
    }
  }
}
Legend 825 866 @(@{Pen=$blue;Label='native side A';Width=125},@{Pen=$orange;Label='native side B';Width=125})
$g.DrawString('8 paired facets · rest area/side 4e-6 m² · plan projection; z-opening not shown',$small,$black,790,887)
$g.DrawString('NOT full-debond proof · damage/fullfailure UNKNOWN · roughness NOT_VALIDATED · contact NOT_RUN · GPU UNKNOWN',$font,$black,38,951)
$g.DrawString('Source/API review uses official v5.0.7 source tag; exact wheel build correspondence UNKNOWN. See COHESIVE_CAPABILITY_20261004.md',$small,$black,38,976)

$bmp.Save($png,[System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose(); $bmp.Dispose()
Get-Item $png | Select-Object FullName,Length
