# Material visibility to native bond graph: source-only bridge candidate

Status: design and static review only. No Julia, CAD, native solver, new asset, or tests were run for this continuation. The frozen classifier/STEP evidence belongs to earlier runs and was not regenerated here. The production runner still does not consume a material-visibility mask; the native geometry gate remains **PENDING**.

## What this directory contains

- `bridge_contract.py` is a candidate packager. It computes the classifier-source SHA from the bytes at an explicit source path, binds the CAD, point CSV, classifier result, and geometry sidecar by actual file hashes, and emits a static mask contract. It does not call OCC or native code.
- `mask_contract.schema.json` is the transport schema. `static_contract_status=VALID` means only that this static package is structurally consistent. The native gate is always `UNVERIFIED`, the raw candidate set is not reproduced, graph readback is `NOT_RUN`, and initialization is never authorized by this packager.
- `test_bridge_contract.py` contains narrow regression tests. They cover source-byte hashing, sidecar raw-hit mismatch, unknown/empty/zero-degree blocking, source-accuracy basis rejection, pre-read size budgets, point/pair caps, and exact sidecar selection. These tests are written but **not run**.
- `pinned_fork_initial_bond_policy.patch` is an unapplied design patch for the next versioned Peridynamics fork. It is not applied to either vendor tree and is not yet an executable patch; it records the Body input field, a stdlib TOML transport, the graph-construction seam, and required recalculation.

The contract keeps classifier/global point IDs (which may be zero-based or sparse) separate from the 1-based native body row indices. `points.rows[*].native_index` maps each external point ID to the exact row in `Body.position`; canonical pairs and pair hashes always use the external IDs.

## Static audit of the pinned solver

The inspected source is `src/pgfrag/resources/julia/vendor/Peridynamics.jl` at the current pinned `0.5.4+pgfrag.1` source. The relevant API and construction path are:

1. `BodyChunk` calls `get_system(body, pd, chunk_id)` before storage and data-handler initialization (`src/discretization/body_chunk.jl:39-48`).
2. The bond-system `get_system` calls `BondSystem(body, pd, chunk_id)`; the constructor calls `get_bond_data`, then sizes correction storage and computes kernels (`src/discretization/bond_system.jl:55-80`).
3. `find_bonds` currently runs the native horizon search and calculates per-point counts, then applies only `filter_bonds!` for predefined cracks before building bond ranges (`src/discretization/bond_system.jl:96-108, 134-175`). No geometry policy is accepted.
4. `get_bond_data` derives halo membership, `ChunkHandler`, localizer, and bond IDs after `find_bonds` (`src/discretization/bond_system.jl:67-74`). A later graph edit would leave these derived values stale.
5. `Body` has no visibility field. Its current fields are internal to the pinned package (`src/discretization/body.jl:83-118`). The explicit constructor keyword proposed here therefore belongs in a versioned fork.
6. `precrack!` represents a physical crack and changes bond activity or initial damage. It cannot represent an `OUT_OR_VOID` geometry exclusion (`src/discretization/predefined_cracks.jl:12-19, 44-60, 63-81`).
7. Initial bond activity and neighbor counts derive from the graph, while energy correction allocates against it and computes its factors during handler initialization. OSB weighted volume and dilatation also iterate the graph (`src/discretization/bond_system.jl:349-359`, `src/discretization/bond_system_corrections.jl:43-101, 140-182`, `src/physics/ordinary_state_based.jl:207-235`). Filtering before halo/storage/handler construction makes these values derive from the admitted graph.

The available `run_fixture.jl` constructs a plain `uniform_box` `Body` and `Job`; it has no CAD or visibility-contract input. The production source has no call from `bond_visibility.py` to that Julia path. This proves an integration gap, not a native pass.

## Contract semantics and bounded inputs

The packager first stats and resolves all five distinct input files (classifier result, point payload, source CAD, raw-hit/interval sidecar, classifier source), checks every per-file cap and a 64 MiB aggregate cap, then performs two bounded reads per file and rejects metadata or byte changes. Current LocalTiny limits are: 8 MiB result, 2 MiB point CSV, 16 MiB CAD, 32 MiB sidecar, 2 MiB classifier source, at most 5,000 points, at most 250,000 candidate pairs, and at most 16 captured sidecar results. Point CSV parsing stops as soon as its point cap is exceeded.

The sidecar validator requires a supported wrapper schema, successful captured run, exact unique `(test_nodeid, call_index)` selection, a `module_file` resolving to the source file whose bytes were hashed, and exact canonical JSON equality between the selected report and classifier result. That equality includes every pair's `raw_hits` and `open_intervals`; a sidecar that only hashes successfully but contains unrelated or mismatched hit evidence is rejected.

The static contract has separate states:

- `static_contract_status=VALID`: paths, byte hashes, units, frames, coordinates, sidecar evidence, counts, and mask structure agree.
- `native_mask_eligibility=BLOCKED_UNKNOWN`: any unknown candidate exists; it blocks initialization of the complete asset.
- `BLOCKED_EMPTY_CANDIDATE_SET`, `BLOCKED_EMPTY_MASK`, and `BLOCKED_ZERO_DEGREE`: no empty candidate or zero-degree graph is treated as qualified.
- `CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED`: the static mask has included candidates and positive degree at every listed point. It remains an unverified candidate until native candidate-set and graph readback succeed.
- `native_gate.status=UNVERIFIED`, `raw_candidate_set_reproduced=false`, `native_graph_readback=NOT_RUN`, and `initialization_authorized=false` in every packager result.
- `initial_physics_events` is explicitly planned metadata with `evidence_status=PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN`. Its planned zero count is not a measurement or runtime receipt. Geometry exclusions never count as physical fracture events.

`source_accuracy_mm` and `source_accuracy_basis` preserve the CAD classifier's declared model precision. They do not assert measured physical accuracy. Input coordinates remain in mm, native coordinates are m, and the explicit rigid transform is recorded. The caller applies `m_per_input_unit=0.001` once, then the named reference-to-native transform.

The stable hash encodings are language-independent:

- Point order: SHA-256 of `pgfrag.point-ids.v1\n` followed by each canonical nonnegative external ID plus `\n` in payload row order.
- Full pair map: SHA-256 of `pgfrag.native-pair-map.v1\n` followed by sorted `i,j,STATE\n` rows, including included, excluded, and unknown decisions.
- Included mask: SHA-256 of `pgfrag.native-included-mask.v1\n` followed by sorted `i,j\n` rows. This never replaces the full pair-map hash.

## Proposed native API and graph-construction order

The patch draft chooses an explicit immutable `InitialBondPolicy` constructor input attached to `Body`, in a separately versioned pinned fork. The proposed call shape is:

```julia
policy = InitialBondPolicy(audit_json_path, native_toml_path) # rehashes both actual byte streams
body = Body(OSBMaterial{EnergySurfaceCorrection}(), position_m, volume;
            initial_bond_policy=policy)
```

The policy stores an immutable typed-TOML string and SHA pins, not a mutable global variable or a precrack. Python continues to emit the strict JSON audit contract. A future Python bridge step compiles its typed projection into a versioned TOML transport that contains the actual audit-JSON byte SHA. Julia reads both current files, verifies the JSON SHA embedded by the transport, computes the actual TOML byte SHA, and records both actual byte hashes in the native receipt. TOML contains typed point rows, external-ID/native-index mapping, candidate decisions, frame/unit transform, horizon, and all source-artifact pins needed by graph construction. The JSON remains the human/review audit artifact; the TOML is the native transport.

This avoids adding an external parser: the fork declares the Julia standard-library `TOML` and `SHA` dependencies directly. `TOML.parse` rejects malformed/duplicate keys; a recursive validator rejects non-finite floats and unsupported TOML values such as date/time objects, while fixed per-file/aggregate byte caps are checked before parsing. The policy retains only immutable strings and digests; a parsed TOML table is a local validation value and is not stored in the policy. The Python JSON→TOML compiler and cross-artifact equality checks are future bridge work and are not implemented in this directory.

At `get_bond_data`/`find_bonds`, the versioned fork must reject anything outside one body, one chunk, uniform horizon, and ordinary state-based (OSB) material for v1. It must then:

1. Check point IDs, native row positions, frame transform, horizon, and all source-file pins against the `InitialBondPolicy` before graph construction.
2. Run the solver's own Euclidean horizon search and collect the complete directed candidate ledger. Convert native row indices to external IDs, canonicalize pairs, and require exactly two reciprocal directed candidates per pair. The canonical set must equal the contract's complete candidate set before any pair is removed.
3. Reject the entire asset if any candidate is `UNKNOWN`, if the candidate set or included mask is empty, or if a point would have zero included degree. Do not drop unknown rows and proceed with the remainder.
4. Keep only `IN_MATERIAL` candidates. Count `OUT_OR_VOID` separately as geometric exclusions. Recompute directed `n_neighbors` from the kept edges before `find_bond_ids`; do not add exclusions to `bond_active`, damage, precrack state, or physical fracture-event output.
5. Only after this filter, build bond ranges, halo points, `ChunkHandler`/localizer, kernels, correction storage, and the native data handler. `EnergySurfaceCorrection` `mfactor`/`scfactor`, OSB weighted volume, and OSB dilatation must be computed from the filtered graph. Ordinary physical precrack behavior remains a later separate solver feature.

The current seam (`find_bonds` after native search and before `filter_bonds!`/`find_bond_ids`) is technically usable inside the pinned source, but it is not a supported extension hook today. For broader solver support, thread the explicit policy through `threads_data_handler → BodyChunk → get_system → BondSystem`, rather than relying on a global. Multi-chunk, MPI, other material systems, contact neighborhoods, and other graph-derived systems are out of v1 scope.

## Native readback and later mechanical acceptance

An eventual native receipt must prove, from the actual initialized system, that the imported fork/version/path and contract SHA are expected; all source and point bytes still match; point IDs/order/positions match; the independently enumerated Euclidean candidates match the entire contract set; and the canonical native graph equals only the included mask in both directions. It must independently read back `n_neighbors`, `bond_ids`, halo/localization, all-true initial active state, zero initial damage, and finite correction/OSB graph-derived values. Receipt fields must distinguish Euclidean candidates, geometric exclusions, and simulation-created physical broken events.

Mechanical acceptance is separate and remains future work: homogeneous undeformed coupon has zero strain energy/force within a preregistered bound; small affine strain has finite positive energy and force with consistent work/energy derivative; axis permutation and rigid rotation preserve energy and rotate force; reciprocal graph entries and global internal-force resultant pass; hollow-source cross-cavity links are absent while same-wall links remain; and free-surface horizon truncation causes neither geometry-excluded pairs nor phase fracture. No statement here claims any such native result.

## Known risks / decisions for review

- The JSON→typed-TOML compiler and checks proving its projection matches the strict audit JSON still need implementation and focused tests. The transport must carry the exact audit JSON byte SHA; the native receipt must record actual JSON and TOML byte SHAs.
- SHA binding and byte-for-byte sidecar consistency prove identity of the supplied artifacts, not that `pytest_exit_code=0` came from a fresh, trusted process. The wrapper fields are unsigned claims; a future native receipt must preserve process/runtime/argv and captured output, and must reread the current source and payload paths before graph construction.
- The current classifier result records point count and pair-endpoint coordinates, but does not include the original ordered point-CSV SHA. The packager checks each candidate endpoint against the supplied rows; a candidate mask that could proceed must also give every listed point positive included degree. The native fork still has to read back every row and independently reproduce the candidate set before initialization.
- Absolute paths and SHA pins do not guarantee artifact freshness by themselves. They identify bytes; native initialization must independently rehash each current file and stop if any path now resolves to different bytes.
- Exact candidate-set equality can reject a mesh if Python and Julia differ on horizon-boundary floating-point conventions. That is a fail-closed result; do not widen tolerance to make the sets pass.
- The one-chunk restriction is a scope assertion, not evidence that the current runner actually chooses one chunk. Native readback must record and validate it.
- A valid static contract does not establish the native CAD/body correspondence until the body points, CAD hash, source, and initialized graph are read back in the fork.
- Changes to `pgfrag/vendor/Peridynamics.jl` or `src/pgfrag/resources/julia/vendor/Peridynamics.jl` are outside this handoff and were not made.
