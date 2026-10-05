# Discretization bias descriptor evaluator

This experiment-local component pairs two JSON run descriptors, checks what changed, transports geometry into one declared rest-material frame, and emits metrics. It is **not** a solver and does not prove mesh bias is absent. Fixtures always yield `FIXTURE_ONLY`; reports always say `NOT_ADMITTED`.

## Run the fixture example

Use a Python environment with matplotlib installed:

```powershell
$py = 'third_party\vesselforge-runtime\venv\Scripts\python.exe'
& $py pgfrag\experiments\discretization_bias_metrics_20261006\make_fixtures.py
& $py pgfrag\experiments\discretization_bias_metrics_20261006\bias_metrics.py `
  --study pgfrag\experiments\discretization_bias_metrics_20261006\fixtures\fixture_study.json `
  --output-dir pgfrag\tmp\review_20261006\discretization_bias_metrics_20261006
& $py -m unittest discover -s pgfrag\experiments\discretization_bias_metrics_20261006\tests -v
```

The fixture set includes a 90-degree local-frame remesh, a translated placement, a resolution pair with a changed material-field seed, a no-fracture descriptor, and graph-only components. Expected report scope is FIXTURE_ONLY.

## Descriptor contract

See `fixtures/*.json` for complete examples. Each descriptor records geometry reference/hash, load hash/units, material parameter and field hashes plus random-field seed/units, solver/build/law/correction, mesh hash, h/delta/dt with units, orientation and translation IDs, common frame ID, 4x4 local-to-common rigid transform, and evidence scope. The transform maps local coordinates to common coordinates; its translation is in the descriptor's `physical_frame.coordinate_units`. Supported units are mm/m/um for coordinates, s/ms/us for time, kg/g for mass, and J/mJ/kJ for energy.

A study manifest pairs run IDs and names exactly one changed factor: orientation, translation, h, dt, or delta. Geometry/load/material field/seed/material parameters/solver build/failure law/surface correction and crack-surface extractor/path-sampling-policy hashes must match. An h pair also holds dt and delta fixed; a dt pair holds h and delta fixed. Mismatches are INCOMPARABLE. Missing provenance, incomplete histories, or no fracture are UNKNOWN.

Crack paths are unordered physical-coordinate samples generated under one fixed, hashed sampling policy; the report computes symmetric nearest-sample p95 and Hausdorff distances after frame transport. Surface area is recomputed from triangle geometry and vertex IDs, so node IDs/order do not define the result. Energy terms are compared by name in joules. Fragment mass distributions count only descriptors with closed, nonoverlapping geometric surfaces; bond-graph components are UNKNOWN and are never relabelled as physical fragments. The mass distribution comparison is 1-D Wasserstein distance over fragment mass shares, with no required fragment count.

Caller tolerances are optional for metric generation but required for a physical assessment. Every provided tolerance must contain a nonnegative `max_abs` and a nonempty `basis`; there are no defaults. Even a fully attested native pair inside tolerances is reported as `WITHIN_TOLERANCES_NOT_BIAS_ELIMINATED`, never as proof of zero bias or A12 eligibility. Physical qualification claims are descriptor attestations for later review, not inferred from fixtures or geometry alone.

