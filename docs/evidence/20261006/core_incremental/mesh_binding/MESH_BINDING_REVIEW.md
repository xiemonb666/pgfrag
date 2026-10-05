# Planar fixture mesh-binding review

This review records one bounded geometry readback of the existing planar BOX split fixture. The new independent entrypoint is `mesh_binding_validator.py`; it is not wired into the existing PG reader, chain, production source, or solver.

## Inputs and binding

The fixture directory contains `explicit_fragment_fixture_output.json`, `interface_pair_ledger.json`, `piece_negative.obj`, and `piece_positive.obj`. It does **not** contain a separate `fragment_meshes.json`. The validator therefore binds the saved OBJ meshes to the piece vertices, triangles, face lineage, and declared metrics in the existing output JSON, then binds those actual face indices to the existing reciprocal pair ledger. No fixture asset was created or edited.

The small input manifest is [planar_box_mesh_binding_manifest.json](../../../experiments/fracture_artifact_contract_20261006/fixtures/planar_box_mesh_binding_manifest.json). It pins SHA-256 digests for both OBJs, the output JSON, and the pair ledger; specifies metres; caps the input at two fragments, 64 combined OBJ vertices, 64 combined triangles, and 256 KiB per input file; and names the actual OBJ object and fracture-face group. Bytes are size-checked before parsing. The fixed parser accepts only three-coordinate vertices and triangular absolute-index faces.

The validator derived closure from actual directed edge incidences, outward orientation from positive signed volume, and paired area/normals from the referenced OBJ triangles. It checked reciprocal vertex IDs and local indices against the same mother-frame coordinates and required reverse winding across each pair. Producer `checks` and `topology` booleans were ignored.

## Readback result

The saved [execution receipt](mesh_binding_execution_receipt.json) records `PASS`, `FIXTURE_GEOMETRY_ONLY`, `NOT_COUPLED`, and `physical_pack_count: 0`. The two actual OBJ meshes have signed volumes 4.44 m³ and 3.56 m³; their sum is 8.00 m³ against the fixture's declared 8.00 m³. The two paired triangles total 4.354859354789773 m² and their actual normals are opposite. The readback took 98.78 ms.

The geometry result is limited to closed, consistently oriented OBJ surfaces, actual paired-face binding, numeric area, and volume-sum conservation. It does not establish non-overlap or complete mother-body coverage: `intersection_test` is `NOT_RUN`, and both proof fields are `NOT_ESTABLISHED`. It does not couple mechanics or qualify a physical or pretraining pack. The contract's observation reader remains the only A12 inference-facing primitive and does not expose clean GT or mother geometry.

## Focused checks

Command:

```powershell
& 'relicweave-synlab/.venv/Scripts/python.exe' -m unittest -v pgfrag.experiments.fracture_artifact_contract_20261006.test_mesh_binding_validator
```

Result: 8 tests passed in 0.070 s. Negatives cover missing pair-face index, modified bytes with stale digest, coordinate changes after re-pinning, missing mesh, globally reversed mesh winding, unit mismatch, and a mesh file one byte over the cap. The positive case reads the saved fixture and checks actual OBJ binding and numerical results. The previously completed 9-test artifact-contract suite was not rerun in this readback pass.

## Preserved source pins

The four pre-read hashes and post-read hashes match exactly; all recorded fixture file lengths are unchanged.

| Existing fixture file | Bytes | Pre-read SHA-256 | Post-read SHA-256 |
| --- | ---: | --- | --- |
| `explicit_fragment_fixture_output.json` | 17,872 | `0b18a24b8f58d1df516743c55cc0cb4c4a3629e5288fbcd920be3c4afecebe16` | `0b18a24b8f58d1df516743c55cc0cb4c4a3629e5288fbcd920be3c4afecebe16` |
| `interface_pair_ledger.json` | 3,973 | `b5e2d578d7f56f28cdf0edf988f0755ff4002b2f3a3b9b6fe96a337672bfd9fe` | `b5e2d578d7f56f28cdf0edf988f0755ff4002b2f3a3b9b6fe96a337672bfd9fe` |
| `piece_negative.obj` | 812 | `5f5e8e3556173f7f214b98d0e1cf87ec3f0ff84898d127e60ecd4b78a3e4a083` | `5f5e8e3556173f7f214b98d0e1cf87ec3f0ff84898d127e60ecd4b78a3e4a083` |
| `piece_positive.obj` | 808 | `342e18222e59ce3e1863b606cab77841bd8035adb6da326573eba73512eec42f` | `342e18222e59ce3e1863b606cab77841bd8035adb6da326573eba73512eec42f` |

## New-file SHA-256 pins

| New file | SHA-256 |
| --- | --- |
| `pgfrag/experiments/fracture_artifact_contract_20261006/mesh_binding_validator.py` | `3cd8018456c8be1e8a17cdbac9cfb4d7b227f09bb3170f8a04083ede2db8d51f` |
| `pgfrag/experiments/fracture_artifact_contract_20261006/test_mesh_binding_validator.py` | `d21fd672bf799b24d6bd4a8fabc78346a4936be2036067552e47a23a085c4a0b` |
| `pgfrag/experiments/fracture_artifact_contract_20261006/fixtures/planar_box_mesh_binding_manifest.json` | `6c490acb9ac1f54b045373e633d1db441361b46669d9b1b559a65e6a19b058b3` |
| `pgfrag/tmp/review_20261006/fracture_artifact_contract_20261006/mesh_binding_execution_receipt.json` | `e14f8eae0ccdeab5b8e909d7748dd225187b122d5cdead2b967129e4a3db81f1` |
