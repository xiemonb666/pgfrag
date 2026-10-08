# PG-Frag native bond geometry-mask bridge: independent static review

Review date: 2026-10-05  
Scope: pgfrag/experiments/material_visibility_native_bridge_20261005/bridge_contract.py (SHA-256 67b8f3090abebb093b7a1e099059c35c2ae09890b4f6b78453dd44d8a6ad4085), the originally supplied test_bridge_contract.py (SHA-256 a0c2c2998f36502a6981d24e16e576348126304cb5140d0f9f9d7b88c0b26fc4), README, mask schema, and unapplied fork pseudopatch. Static inspection only; I did not run pytest, Julia, CAD, or generated/synthetic cases.

## Assessment

The builder binds the actual bytes it captures for the classifier result, ordered point CSV, CAD, classifier source, and sidecar. It checks per-file and aggregate size caps before opening inputs, reads each input twice with metadata checks, rejects duplicate JSON object keys and nonstandard JSON constants, checks the source CAD digest against both digest claims in the classifier result, and records the actual SHA-256 values. It verifies row order via the CSV digest and point-ID-order digest; every submitted pair has canonical increasing external IDs, known endpoint IDs, matching endpoint coordinates, and exactly two reciprocal direction IDs/states. The selected sidecar record must be unique and canonically equal to the complete result object, including raw_hits and open_intervals. The transform is finite homogeneous affine with a proper orthonormal rotation, the horizon is finite and positive, and the source/reference units and reference frame are checked.

For the supplied candidate ledger, any listed UNKNOWN blocks the whole asset; empty candidate set, empty included mask, and zero included degree are also blocked. initialization_authorized is hard-coded false; native gate status stays UNVERIFIED, raw candidate reproduction stays false, and graph readback stays NOT_RUN. Therefore this source cannot itself authorize native initialization.

## Findings

### P2 — The static builder does not prove that the submitted pair ledger is the Euclidean horizon candidate set

bridge_contract.py:222-296 validates each listed pair and its endpoints, but never computes pair distances or candidate-set completeness. bridge_contract.py:441-443 only validates that the caller-supplied horizon_m is finite and positive. The output nevertheless labels the rows raw_euclidean_candidate_pairs with candidate_rule=exact_pinned_Peridynamics_find_bonds_convention at lines 570-576, and lines 523-532 can assign CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED from the supplied states and degree counts alone.

Concrete counterexample accepted by these checks: two CSV points, external IDs 0 and 4, at 0 mm and 1 mm; horizon_m=0.0005; a self-consistent classifier result contains only pair [0,4] marked IN_MATERIAL, has both reciprocal directions and matching endpoint coordinates/counts, and the selected sidecar contains that exact result. With source/CAD claims matching the actual bytes, the builder returns static_contract_status=VALID and the positive-degree native-unverified candidate state even though the pair distance is 0.001 m, greater than the 0.0005 m horizon. Likewise, omitted Euclidean candidates—including omitted UNKNOWN candidates—cannot trigger the unknown block.

This is a false positive at the static-candidate stage, not evidence of a native pass: the output remains explicitly unverified and the pseudopatch says the native constructor must compare the complete solver candidate set before filtering. README lines 85 and 87 disclose this remaining native check. The static status/candidate fields should not be consumed as proof that the pair ledger was actually generated from the stated horizon.

### P2 — The JSON Schema does not encode the cross-field invariants enforced by the current builder

mask_contract.schema.json:15-20 permits the candidate eligibility string independently of counts and masks. Lines 35-41 independently type native gate flags. Lines 115-140 permit empty arrays and do not correlate pair arrays with counts or booleans; $defs.canonical_pair at lines 183-186 only requires two nonnegative integers and does not require i < j or unique pairs. Thus a hand-edited document can be schema-valid while claiming candidate eligibility with zero candidates, inconsistent nonzero flags, reversed/duplicate pairs, or counts that do not match its arrays. The current build_contract constructs consistent values, so this is not a demonstrated false output from that function. It is a downstream-consumer hazard: schema validation alone cannot validate a received/tampered contract. The pseudopatch’s validate_typed_transport and the JSON-to-TOML compiler/cross-artifact checks remain pseudocode/future work (patch lines 7, 51-62; README lines 61-63 and 83).

### P2 — Sidecar success and classifier provenance are self-asserted, not authenticated

bridge_contract.py:323-368 checks a literal pytest_exit_code == 0, matching algorithm text, a resolvable module_file path, unique selector, and exact selected-result equality. It does not authenticate the process, invocation, output capture, or that the reported result was produced by executing the bytes currently hashed at that path. A fabricated result and matching sidecar can satisfy these structural checks, then receive the candidate label if its counts/degrees are consistent. Actual hashes bind the submitted bytes and detect mismatches; they do not establish who produced the result or whether the classifier ran successfully. README line 84 states this limitation, so treat it as an explicit evidence boundary, not a native-gate failure. The pseudopatch’s future native check must retain process/runtime/argv/output evidence if trusted-run provenance is required.

### P3 — Raw-hit/interval evidence is equality-bound but not semantically validated

bridge_contract.py:245-269 only requires raw_hits and open_intervals to be lists and validates reciprocal IDs/states; it does not validate hit/interval item shapes or consistency between intervals and the classification. _validate_geometry_sidecar at lines 366-368 proves that the selected sidecar copy matches the classifier result, not that the hit/interval evidence is meaningful. An IN_MATERIAL row with empty lists and matching reciprocal evidence is structurally accepted. Whether empty intervals are invalid depends on the classifier’s evidence schema/semantics, which is not declared in this bridge contract; the bridge currently cannot make that determination.

## Checks and scope notes

- Strict JSON, actual-byte SHA binding, input size preflight, coordinate/frame validation, canonical pair ID/order checks, reciprocal directions, sidecar unique selection and result equality, and the three block states are implemented in the inspected builder.
- There is no binding from the supplied horizon_m to the classifier result, nor an independent candidate-set reconstruction in Python. The eventual native equality check described in the README/pseudopatch is essential.
- The UNKNOWN block applies to unknown rows present in the submitted ledger only; it does not establish ledger completeness.
- README, schema, and pseudopatch consistently describe the bridge as source-only/design-only; no native integration, native graph readback, or mechanical acceptance has been completed.
- The supplied test source was not executed. The parent reports correcting the expected zero-degree IDs in its working test file and will run the final focused suite separately; that update is outside this static review’s original test SHA baseline.
- Input paths are not explicitly required to be regular files before Path.open; special files such as FIFOs can defeat prompt completion even when their stat size appears within budget. This is an availability hardening point, not a path to a VALID contract.
