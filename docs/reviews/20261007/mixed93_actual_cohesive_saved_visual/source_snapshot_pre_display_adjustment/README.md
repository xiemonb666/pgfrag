# Mixed93 actual cohesive saved-result visual

This is a source-only visual package for one future accepted diagnostic result. The planned input is `pgfrag/tmp/review_20261006/mixed93_actual_cohesive_bridge_20261006/numeric_runs/run01/result.json`; its SHA-256 is intentionally `null` because the actual cohesive worker result does not exist yet. The case and source pins keep render authorization false. No result payload was opened and no figure was rendered while preparing this package.

The reader follows the worker's frozen result schema `mixed93_actual_cohesive_bridge_numeric_result.v1` from `run_actual_cohesive_bridge.py` SHA-256 `36ac1710286b0abf607bda5b1c6ffab893f8191830eae091517e2acd41d62236`: saved interface coordinates/weights are in `evidence.interface`, the two 16-point ledgers are in `evidence.baseline.points` and `evidence.coupled.points`, and the actual q-change/force summary is in `solution`. Point IDs are joined to the interface's saved `quadrature_id`; the renderer does not build new geometry, query a field, or reevaluate an interface point. It refuses any result that is not `ACCEPTED_DIAGNOSTIC_ONLY` with all gates passed, has a different point count/order, or leaves the declared virgin elastic opening interval.

The future PNG has four panels. Two x/z views show baseline and coupled-final opening at the same 16 saved QPs, with equal spatial axes and one fixed shared color scale from zero to the toy `delta0` in metres. The traction-opening panel overlays the actual saved baseline/coupled values on the virgin line `t=K0*delta` only from zero through `delta0`, marks `delta0`, and annotates the actual saved `||q_c-q_0||` and internal-force norm. It does not invent a post-`delta0` softening curve. The sign panel uses saved point vectors and normal to show the internal bookkeeping share `+t n` beside physical upper-face `-t n` and lower-face `+t n`; the internal share is not mislabeled as an applied face traction.

Visible title labels are `FIXED_CRACK_TOY_ELASTIC_DIAGNOSTIC`, `NO_GROWTH`, and `NO_MATERIAL_CALIBRATION`. All plotted quantities carry SI units. The footer states the saved-data and no-growth boundary. The output is one PNG, bounded to 15 seconds and 2 MiB combined output; the runner refuses overwrite and never retries.

The runner is standard-library-only. Running it without flags performs source/hash/schema readiness checks and reports the pending result pin. Matplotlib is imported only inside `render_saved_bridge.py` after root authorization, matching result SHA, accepted status, saved-data checks, and exact GO-token checks. After the actual worker result exists, root must review it and then bind its path/SHA in both `case.json` and `SOURCE_PINS.json`, update those files' source pins, and explicitly authorize one render with a root token. The source check does not grant that authorization.

Use the configured VesselForge Python runtime for the future source check; do not render from this SOURCE_READY package:

```powershell
C:/Users/xiemo/Desktop/ShardForge/third_party/vesselforge-runtime/venv/Scripts/python.exe -B pgfrag/tmp/review_20261007/mixed93_actual_cohesive_saved_visual/run_saved_bridge_visual.py
```

This visual communicates only a fixed-crack toy elastic cohesive diagnostic. It does not establish real-material behavior, process-zone resolution, mixed-mode qualification, or physical crack growth.
