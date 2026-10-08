$ErrorActionPreference='Stop'
$runOut='C:\Users\xiemo\Desktop\ShardForge\relicweave-synlab\tmp\review_20261004\implementation\cohesive_bilinear_patch_20261004\controlled_intrinsic_bilinear_modeI_201_run01'
$data=Get-Content -LiteralPath (Join-Path $runOut 'plot_data.json') -Raw | ConvertFrom-Json
try { Add-Type -AssemblyName System.Drawing.Common -ErrorAction Stop } catch { Add-Type -AssemblyName System.Drawing -ErrorAction Stop }
$bmp=[System.Drawing.Bitmap]::new(1600,1080)
$gfx=[System.Drawing.Graphics]::FromImage($bmp)
$gfx.SmoothingMode=[System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$gfx.Clear([System.Drawing.Color]::White)
$font=[System.Drawing.Font]::new('Arial',12)
$small=[System.Drawing.Font]::new('Arial',10)
$title=[System.Drawing.Font]::new('Arial',17,[System.Drawing.FontStyle]::Bold)
$black=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(45,45,45),2)
$grid=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(225,225,225),1)
$blue=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(32,96,180),3)
$orange=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(220,112,38),3)
$green=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(42,145,80),3)
$red=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(190,54,54),2)
$purple=[System.Drawing.Pen]::new([System.Drawing.Color]::FromArgb(128,74,160),3)
$orange.DashStyle=[System.Drawing.Drawing2D.DashStyle]::Dash
function Draw-Panel($rect,$heading,$xlabel,$ylabel,$xmin,$xmax,$ymin,$ymax) {
  $gfx.DrawString($heading,$title,[System.Drawing.Brushes]::Black,$rect.X,$rect.Y-29)
  for($i=0;$i -le 4;$i++){
    $x=$rect.X+$i*$rect.Width/4.0
    $y=$rect.Y+$i*$rect.Height/4.0
    $gfx.DrawLine($grid,[single]$x,[single]$rect.Y,[single]$x,[single]($rect.Y+$rect.Height))
    $gfx.DrawLine($grid,[single]$rect.X,[single]$y,[single]($rect.X+$rect.Width),[single]$y)
  }
  $gfx.DrawRectangle($black,$rect)
  $gfx.DrawString($xlabel,$font,[System.Drawing.Brushes]::Black,$rect.X+$rect.Width/2-80,$rect.Y+$rect.Height+8)
  $gfx.DrawString($ylabel,$font,[System.Drawing.Brushes]::Black,8,$rect.Y+$rect.Height/2)
  $gfx.DrawString($xmin.ToString('0.###'),$small,[System.Drawing.Brushes]::Black,$rect.X-8,$rect.Y+$rect.Height+2)
  $gfx.DrawString($xmax.ToString('0.###'),$small,[System.Drawing.Brushes]::Black,$rect.X+$rect.Width-30,$rect.Y+$rect.Height+2)
  $gfx.DrawString($ymin.ToString('0.###'),$small,[System.Drawing.Brushes]::Black,$rect.X-42,$rect.Y+$rect.Height-10)
  $gfx.DrawString($ymax.ToString('0.###'),$small,[System.Drawing.Brushes]::Black,$rect.X-42,$rect.Y-8)
}
function Map-X($x,$rect,$xmin,$xmax){ return [single]($rect.X+($x-$xmin)/($xmax-$xmin)*$rect.Width) }
function Map-Y($y,$rect,$ymin,$ymax){ return [single]($rect.Y+$rect.Height-($y-$ymin)/($ymax-$ymin)*$rect.Height) }
function Draw-Series($x,$y,$rect,$xmin,$xmax,$ymin,$ymax,$pen) {
  for($i=1;$i -lt $x.Count;$i++){
    $x1=Map-X ([double]$x[$i-1]) $rect $xmin $xmax
    $x2=Map-X ([double]$x[$i]) $rect $xmin $xmax
    $y1=Map-Y ([double]$y[$i-1]) $rect $ymin $ymax
    $y2=Map-Y ([double]$y[$i]) $rect $ymin $ymax
    $gfx.DrawLine($pen,$x1,$y1,$x2,$y2)
  }
}
function Draw-MapTriangle($tri,$rect,$pen,$xmin,$xmax,$ymin,$ymax) {
  for($i=0;$i -lt 3;$i++){
    $j=($i+1)%3
    $x1=Map-X ([double]$tri[$i][0]) $rect $xmin $xmax
    $x2=Map-X ([double]$tri[$j][0]) $rect $xmin $xmax
    $y1=Map-Y ([double]$tri[$i][1]) $rect $ymin $ymax
    $y2=Map-Y ([double]$tri[$j][1]) $rect $ymin $ymax
    $gfx.DrawLine($pen,$x1,$y1,$x2,$y2)
  }
}
$gfx.DrawString('Controlled intrinsic bilinear Mode I patch — native cell-average fields',$title,[System.Drawing.Brushes]::Black,45,25)
$r1=[System.Drawing.Rectangle]::new(95,105,1430,245)
Draw-Panel $r1 'Traction-opening: native mean norm and independent bilinear reference' 'Opening (nm)' 'Traction (MPa)' 0 ([double](($data.commanded_nm | Measure-Object -Maximum).Maximum)) 0 80
Draw-Series $data.opening_nm $data.traction_mean_MPa $r1 0 ([double](($data.commanded_nm | Measure-Object -Maximum).Maximum)) 0 80 $blue
Draw-Series $data.commanded_nm $data.analytic_traction_MPa $r1 0 ([double](($data.commanded_nm | Measure-Object -Maximum).Maximum)) 0 80 $red
$gfx.DrawLine($blue,1120,75,1160,75);$gfx.DrawString('native |T| mean',$small,[System.Drawing.Brushes]::Black,1165,67)
$gfx.DrawLine($red,1320,75,1360,75);$gfx.DrawString('analytic T',$small,[System.Drawing.Brushes]::Black,1365,67)
$maxE=[double](($data.dissipated_uJ | Measure-Object -Maximum).Maximum);$maxE=[Math]::Max($maxE,[double]$data.GcA_uJ)*1.15
$r2=[System.Drawing.Rectangle]::new(95,450,1430,245)
Draw-Panel $r2 'Native cohesive energies (diagnostic; all DOFs prescribed)' 'Commanded opening (nm)' 'Energy (µJ)' 0 ([double](($data.commanded_nm | Measure-Object -Maximum).Maximum)) 0 $maxE
Draw-Series $data.commanded_nm $data.dissipated_uJ $r2 0 ([double](($data.commanded_nm | Measure-Object -Maximum).Maximum)) 0 $maxE $green
Draw-Series $data.commanded_nm $data.reversible_uJ $r2 0 ([double](($data.commanded_nm | Measure-Object -Maximum).Maximum)) 0 $maxE $purple
$gcY=Map-Y ([double]$data.GcA_uJ) $r2 0 $maxE;$gfx.DrawLine($orange,[single]$r2.X,$gcY,[single]($r2.X+$r2.Width),$gcY)
$gfx.DrawLine($green,1075,420,1115,420);$gfx.DrawString('D native',$small,[System.Drawing.Brushes]::Black,1120,412)
$gfx.DrawLine($purple,1200,420,1240,420);$gfx.DrawString('R native',$small,[System.Drawing.Brushes]::Black,1245,412)
$gfx.DrawLine($orange,1325,420,1365,420);$gfx.DrawString('Gc*A',$small,[System.Drawing.Brushes]::Black,1370,412)
$r3=[System.Drawing.Rectangle]::new(690,815,200,200)
$gfx.DrawString('Rest interface at true scale: two coincident native sides; no opening exaggeration',$title,[System.Drawing.Brushes]::Black,180,775)
Draw-Panel $r3 'Rest interface XY projection' 'x (mm)' 'y (mm)' 0 2 0 2
foreach($tri in $data.interface_side_A_xy_mm){Draw-MapTriangle $tri $r3 $blue 0 2 0 2}
foreach($tri in $data.interface_side_B_xy_mm){Draw-MapTriangle $tri $orange 2 0 2 0}
$gfx.DrawLine($blue,965,870,1005,870);$gfx.DrawString('native side A',$small,[System.Drawing.Brushes]::Black,1010,862)
$gfx.DrawLine($orange,965,900,1005,900);$gfx.DrawString('native side B (coincident)',$small,[System.Drawing.Brushes]::Black,1010,892)
$path=Join-Path $runOut 'bilinear_patch_review.png'
$bmp.Save($path,[System.Drawing.Imaging.ImageFormat]::Png)
$gfx.Dispose();$bmp.Dispose()
Write-Output "PNG=$path"