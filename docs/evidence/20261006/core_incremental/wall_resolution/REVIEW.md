# Thin-wall resolution guard review

Status: **local analytic fixture qualification only**. The utility is isolated under `pgfrag/experiments/wall_resolution_guard_20261006/`; it does not modify or import the old ingress worker, invoke CAD, or launch a solver.

The old source is the bounded 96-cell ROI ingress experiment at `pgfrag/experiments/forge_volume_ingress_20261005/`. Its recovered local ledger reconciles 8.01676271854445 mm³ across the 96 cells and 1.68352017089433e-05 kg at 2100 kg/m³. The independent audit reports exact R=8/10 mm annular-cylinder ROI integration and local point-weight checks. Those records do not establish global minimum wall thickness or thickness-direction layer count; the producer also remained blocked at report finalization and did not promote a candidate.

This checker's exact local fixture is a hollow cylinder with `t=2 mm`, `Rin=8 mm`, `Rout=10 mm`, and a 2 mm axial patch. `h=0.5 mm` uses 64 exact annular-wedge cells. Every declared local normal sample traverses 4 actual cell intervals and contains 4 distinct projected centroid layers. Its result is **LOCAL_RESOLUTION_PRECONDITION_MET**. The `h=2 mm` fixture has 16 cells but only 1 actual interval and centroid layer, so its result is **UNDERRESOLVED**. The weights are exact sector volumes `0.5 × (Rout²−Rin²) × Δtheta × Δz`; no `h³` substitute is accepted.

The historical `h=8 mm` against the source's nominal 2 mm wall parameter is explicitly **UNDERRESOLVED** (`h/t=4.0`). This is a coarse nominal-parameter comparison, not a recovered local mesh measurement or a global wall certificate.

The declared three-interval and three-centroid thresholds are fixture-specific engineering preconditions. They are not a universal crack-resolution rule. Process-zone adequacy is a separate required check. Even the fine fixture reports `physical_status=NOT_QUALIFIED`, `physical_pass=false`, and `global_pass=false`; missing full-domain coverage can only produce a global UNKNOWN. The line-segment adapter uses exact quadratic intersections with inner and outer cylinder radii plus axial bounds; it rejects a long segment whose two material endpoints are joined through the hollow cavity.

The qualification run took 0.0036 s (limit 30 s), with 64 / 16 cells and 48 / 12 segments. See `thin_wall_qualification_report.json`, `thin_wall_qualification_inputs.json`, and `actual_section_samples.png` for raw review artifacts.

Run tests with:

```powershell
third_party/vesselforge-runtime/venv/Scripts/python.exe -m pytest -q pgfrag/experiments/wall_resolution_guard_20261006/test_thin_wall_qualification.py
```
