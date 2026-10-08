# Fracture artifact contract review

## Result

New isolated contract files are under `pgfrag/experiments/fracture_artifact_contract_20261006/`. The actual tiny payload-backed sample is `fixtures/sample_fracture_artifact.json`; machine-readable checks and case outcomes are in `fracture_artifact_review.json` beside this note.

The sample is contract-valid and reports `FIXTURE_GEOMETRY_ONLY`, `physical_receipt_verified=false`, and `pretraining_qualified=false`. It carries two explicit closed fragment mesh references, mass/volume, one reciprocal two-sided patch in mother coordinates, reverse triangle winding and opposite unit normals, one zero-time fixture event with lineage, one geometry-only clean label, and two observation fragment references. The reference area is recomputed from the patch triangle coordinates and must match the declared positive area.

## Validation

The Draft 2020-12 schema passed `jsonschema.Draft202012Validator.check_schema` and validates the sample. The bounded focused suite ran 9 tests in 0.070 s and passed. It covers the complete fixture, missing labels, cross-fragment surface mismatch, units mismatch, strict non-finite rejection, GT leakage, known zero versus explicit unknown, receipt checkbox rejection, graph-components-only rejection, reciprocal triangle winding/area, observation/GT channel separation, artifact digest checks and referenced payload hash verification. The fixture table contains 7 countable rows; no solver was started and no vessel batch was generated.

The strict observation reader returns only observation rows and verified observed mesh/mask bytes (base64), with no package path to browse. It does not return mother identity, event/surface lineage or clean labels. Clean GT requires the separate `purpose="offline_audit"` reader and is tagged audit-only.

## Qualification and integration boundary

No checkbox can establish physical qualification. `PHYSICALLY_VALIDATED` requires calibrated material plus a structured receipt accepted by a caller-supplied trusted verifier; the default reader/exporter has no trust root and rejects promotion. The current sample's material is explicitly `UNCALIBRATED` and toughness is explicitly `UNKNOWN`; the zero fixture load/time and zero observation count/area remain `KNOWN` values.

The existing `pgfrag.result.v1` reader and SynLab chain dispatch were read only. Their documented and tested scope stays `PG_DIAGNOSTIC_ONLY`: endpoint point/component/bond rows do not produce closed fragments, crack surfaces, event histories or R2 labels. The new files provide an input contract for a later A12 adapter; no A12 or existing reader/chain integration is claimed here.
