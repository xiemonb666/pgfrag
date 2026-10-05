# Thin-wall resolution guard

This is an isolated analytic review component for the wall-resolution failure mode noted after `forge_volume_ingress_20261005`. It does not modify that experiment or the `pgfrag` production API, and it imports no CAD kernel or solver.

The input schema requires unit-tagged cylinder wall samples, exact annular-wedge cell boundary geometry and local max-edge values, centroid positions with positive volume pointweights and source references, and unit-tagged graph segments. Cell volume comes from the analytic annular-wedge formula, not from `h³`. Cell boundary intervals crossed by every declared wall-normal sample and distinct volume-centroid layers are reported as separate evidence. Segment validation checks the whole line against the hollow-cylinder radii and axial bounds using exact quadratic breakpoints.

`make_analytic_hollow_cylinder_fixture()` creates the bounded review fixtures: `t=2 mm`, `Rin=8 mm`, `Rout=10 mm`, `h=0.5 mm` (64 cells, 48 segments) and `h=2 mm` (16 cells, 12 segments). Both are local partial patches. The three-interval/three-centroid rule is an explicitly declared engineering precondition for this review, not a general crack-resolution guarantee; process-zone checks remain separate. No report can set `physical_pass` or `global_pass` true.

Run the tests and materialize the JSON report plus actual section-samples figure with the bundled runtime:

```powershell
third_party/vesselforge-runtime/venv/Scripts/python.exe -m pytest -q pgfrag/experiments/wall_resolution_guard_20261006/test_thin_wall_qualification.py
third_party/vesselforge-runtime/venv/Scripts/python.exe pgfrag/experiments/wall_resolution_guard_20261006/write_review_artifacts.py
```

The generated evidence lives in `pgfrag/tmp/review_20261006/wall_resolution_guard_20261006/`. The fixture tests exercise mm-to-m conversion, bad units, positive but incorrect `h³` pointweights, a long segment through the cavity with material endpoints, empty segment evidence, missing or duplicated local coverage rays, volume/mass reconciliation, and explicit dispositions for `h=2 mm` and historical `h=8 mm`.
