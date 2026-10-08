# Initial bond policy native body pilot — source only

This directory prepares a bounded native-graph pilot for the isolated initial-bond-policy fork. It has not been executed. The case is explicitly `GEOMETRY_FILTER_FIXTURE_ONLY`: it is a synthetic 16-point cloud made from two disconnected 2×2×2 point boxes with an empty x-gap. It does not use CAD, validate CAD accuracy, represent a physical specimen, establish a scientific material classifier, or demonstrate fracture behavior.

The fixed point payload is `fixedcase_points.csv`, in exact native Body row order. `fixedcase_components.csv` labels the two toy boxes. The point IDs are intentionally interleaved with respect to native row order so tests must map canonical external `(i,j)` pairs to native indices rather than assuming those orderings agree. Coordinates are millimetres and the native body frame uses the explicit identity transform plus one `0.001` mm-to-m conversion.

The 1.8 mm horizon includes all 28 within-box pairs in each box (the largest within-box distance is `sqrt(3)` mm), plus four across-gap pairs at 1.5 mm. The next across-gap distance is `sqrt(3.25)` mm, greater than the horizon. The expected raw set is therefore 60 canonical pairs / 120 directed candidates. The positive fixture marks the 56 within-box pairs included and four across-gap pairs excluded; each included point has degree 7. No pair lies exactly on the horizon. The native PointNeighbors implementation uses an inclusive squared-distance comparison, `distance2 <= search_radius^2`.

`generate_fixture_policy.py` writes a positive transport bundle and a second, deliberately invalid zero-degree probe bundle into a new disposable directory under `pgfrag/tmp/initial_bond_policy_body_pilot_20261005/`. The positive fixture labels are defined only by the two toy component labels. The blocked variant additionally forces the seven within-box candidate pairs incident to external ID 103 to `OUT_OR_VOID`, yielding the deliberate static gate state `BLOCKED_ZERO_DEGREE`; that state is a negative-test construction, not geometry truth. The generated audit uses `VALID` and `READY` only because the pinned transport schema requires those tokens. Its CAD-named path is the explicit non-CAD `fixedcase_source_marker.txt`; its zero “accuracy” is explicitly a schema placeholder. Neither receipt may be used for production or described as a CAD/compiler scientific receipt.

The generator invokes the pinned transport-only compiler from the sibling bridge experiment. It computes actual file SHA-256 values for the fixture marker, ordered CSV, generator source, classifier-result fixture, and interval-evidence fixture, then emits fixture-only receipt hashes. The native driver requires the JSON and TOML SHA values as separate CLI arguments; it does not trust the digests embedded only in TOML. These hashes establish fixture artifact identity, not classifier provenance.

The Julia driver constructs a real `OSBMaterial` Body with exact uniform 1.8 mm horizon and one `PointDecomposition` chunk. It first calls the native unmasked `find_bonds` search and compares every directed candidate against an independently reconstructed Euclidean set; it then constructs a real `BondSystem` with `InitialBondPolicy` attached and reads back the filtered native graph, degrees, neighbors, bond-ID ledger, and admission receipt. It checks both directed orientations of every excluded pair are absent and both orientations of every included pair remain. The negative policy must be rejected before native `BondSystem` admission.

OSB has no persistent `weighted_volume` field. The pilot calls the pinned native `calc_weighted_volume` method with a fresh all-active scratch `OSBStorage`, then compares it to an independent sum over the expected post-mask pair set using the same OSB formula and the default `NoCorrection` (`β=1`). This is a read-only formula check on the admitted graph; it is not a `BodyChunk`/ThreadsBodyDataHandler initialization, force evaluation, time integration, solver run, or physical event observation. The driver deliberately constructs no `Job`, `BodyChunk`, `ThreadsBodyDataHandler`, time solver, or CAD workflow. Consequently this source pilot proves the one-body/one-chunk native `BondSystem` seam if later run, but does not claim a complete Threads handler initialization.

## Planned invocation after source review

Use the already verified disposable fork at `pgfrag/tmp/initial_bond_policy_native_20261005/fork_run_03/project` only as a read-only Julia environment. Keep generated fixture files in a new scratch directory and do not run `Pkg.resolve`, package installation, tests, CAD, a `Job`, a solver, force evaluation, or fracture. The bridge compiler runs only to write the two fixture bundles. Then pass the four SHA values printed in its JSON receipt as the driver's explicit JSON/TOML pin arguments.

Example shape (replace `RUN_DIR` and the four hashes with the generated fixture-only receipt values):

```text
python pgfrag/experiments/initial_bond_policy_body_pilot_20261005/generate_fixture_policy.py pgfrag/tmp/initial_bond_policy_body_pilot_20261005/RUN_DIR
julia --project=pgfrag/tmp/initial_bond_policy_native_20261005/fork_run_03/project `
  pgfrag/experiments/initial_bond_policy_body_pilot_20261005/native_body_pilot.jl `
  pgfrag/tmp/initial_bond_policy_body_pilot_20261005/RUN_DIR/positive_mask/audit_contract.json `
  pgfrag/tmp/initial_bond_policy_body_pilot_20261005/RUN_DIR/positive_mask/native_mask.toml `
  GOOD_JSON_SHA GOOD_TOML_SHA `
  pgfrag/tmp/initial_bond_policy_body_pilot_20261005/RUN_DIR/zero_degree_block/audit_contract.json `
  pgfrag/tmp/initial_bond_policy_body_pilot_20261005/RUN_DIR/zero_degree_block/native_mask.toml `
  BLOCKED_JSON_SHA BLOCKED_TOML_SHA `
  pgfrag/tmp/initial_bond_policy_body_pilot_20261005/RUN_DIR/positive_mask/native_graph_readback.txt
```

The required external runner budget is at most 125 points, 20,000 directed candidates, 60 seconds, and 20 MiB for fixture artifacts/working output. The source asserts the point, directed-edge, elapsed-time, and result-file caps. A 20 MiB whole-Julia-process RSS cap is not established by this source; if that was intended, resolve the runner’s process-memory measurement before authorizing execution because Julia/package import overhead is outside the fixture’s bounded allocations. Stop after one attempt on any import, assertion, timeout, or resource failure.

## Status

- Source preparation: complete; not frozen until the source hashes in `source_reference_pins.md` are filled.
- Julia import / native Body / raw search / filtered graph / weighted-volume formula: `NOT_RUN`.
- Physical failures, fractures, force, solver, time integration, CAD, and actual ROI adapter: `NOT_RUN` / out of scope.
