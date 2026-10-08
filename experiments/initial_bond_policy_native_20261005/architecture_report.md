# Initial-bond policy native overlay

Status: source-only staged fork overlay for root review. No vendor, production, current `Project.toml`/`Manifest.toml`, or Git state was modified. No Julia, tests, CAD, native solver, or synthetic geometry was executed.

## Frozen inputs and staged files

The overlay is based on the pinned Peridynamics source at `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl`. Exact base file SHA-256 values, including the copy-only license, are in `source_pin_manifest.md`. The staged source copies are under `overlay/`; the unapplied unified patch is `initial_bond_policy_overlay.patch`.

The patch was reconstructed in memory against each pinned base file and compared line-for-line with the staged overlay. This checks patch shape only; it does not apply the patch to a fork or compile Julia. It contains changes to six existing source files plus the new policy helper. Its isolated `Project.toml` hunk adds only the Julia stdlib UUIDs for `SHA` and `TOML`; no manifest is included.

The transport compiler read during this implementation was `pgfrag/experiments/material_visibility_native_bridge_20261005/native_mask_transport.py`, SHA-256 `6779b3242ef87c13d5626a1e6d89414e1cdd7d6ed8f7e264fb01ed83edef7ea9`. Its TOML root and nested table keys match the helper's exact-key schema. The compiler still needs its own frozen-source review and trusted job/compiler receipt before a native caller supplies expected JSON/TOML hashes.

## Policy inputs and validation

`InitialBondPolicy(json_path, toml_path; expected_json_sha256=..., expected_toml_sha256=...)` requires both expected digests from a root-authorized job/compiler receipt. It bounded-reads actual JSON and TOML bytes, requires both expected digests to match, and also checks the transport's JSON path and embedded JSON digest. The TOML text and actual digests are retained in immutable struct fields. Before bond construction, it rereads and rehashes both files and the CAD, point CSV, classifier-result JSON, classifier source, and raw-hit sidecar.

There is no JSON parser dependency in the Julia fork. Strict JSON syntax/schema validation remains the compiler job's responsibility. The native helper pins the actual JSON bytes and validates the complete typed TOML projection (exact tables and keys, typed values, caps, finite values, sorted full pair ledger, external/native row maps, all state partitions, recomputed counts and hashes, zero-degree report, gate state, and planned-only physical-event metadata). Hashes bind bytes; they do not prove who ran the compiler or that a declared test process succeeded. Consequently the native API requires both expected digests as explicit keywords instead of trusting the TOML's self-declared JSON SHA.

The point CSV is independently parsed from its actual pinned bytes with the expected four-column header, canonical decimal IDs, unique IDs, finite coordinates, and row cap. Parsed IDs/order and all three binary64 coordinates must equal the TOML rows exactly. These values are then compared to `Body.position` after one `0.001` mm-to-m conversion and the declared proper homogeneous transform, using `atol=1e-12 m`, `rtol=0`. The transform must be a finite rigid orientation-preserving transform; the horizon must match exactly at every point.

The helper validates the native transport as `TRANSPORT_ONLY_NATIVE_UNVERIFIED`, with native readback `NOT_RUN` and initialization authorization false. The raw-hit sidecar is reread and SHA-bound, but its contents and any declared pytest exit code are not independently reinterpreted as a test-run proof by Julia.

## Native construction path

The pinned lifecycle is `BodyChunk → get_system → BondSystem → get_bond_data → find_bonds → find_halo_points/ChunkHandler/localize/find_bond_ids` (`overlay/src/discretization/bond_system.jl:56-137`; the unchanged construction-order source is pinned in the manifest). The geometry filter is inserted directly after the complete solver-native `find_bonds!` directed search and before physical `filter_bonds!` (`bond_system.jl:123-137`); therefore halo discovery, localization, bond IDs, kernels, correction storage, and later OSB weighted-volume setup consume the admitted graph. Physical precrack filtering remains later and separate.

The policy type is included before `body.jl`, while helper signatures avoid the later-defined `PointDecomposition` and `Bond` names so the include order remains valid (`overlay/src/Peridynamics.jl:111`; helper signatures at `initial_bond_policy.jl:776,868`). `Body` has one policy field and a backward-compatible keyword default (`overlay/src/discretization/body.jl:100-128`). For policy-bearing bodies, v1 accepts only `OSBMaterial`, exactly one native chunk covering rows `1:N`, and a positive uniform horizon equal to the contract value (`initial_bond_policy.jl:776-822`). Policy-bearing bodies are rejected at the Threads handler and all MPI body/multibody entry points (`threads_multibody_data_handler.jl:37-40`, `mpi_body_data_handler.jl:81-83`, `mpi_multibody_data_handler.jl:2-4`). The MPI single-body guard runs before point decomposition, so even a single-rank/one-chunk MPI launch cannot enter this policy path. MPI multibody handling remains unimplemented in this pinned source and reports that explicitly after its policy guard. The supported native policy path is Threads CPU, one body, and one chunk.

The independent critic identified that guarding only MPI multibody setup left the `mpi_data_handler(body, solver)` overload able to reach `PointDecomposition` first. The updated overlay rejects any policy-bearing body at that overload's first statement. The focused source test checks the guard accepts a legacy body and rejects a policy-bearing body; it is staged and not executed. See the appended review addendum for the exact static evidence and remaining execution boundary.

Before filtering, every native directed edge is mapped from native row IDs to external point IDs. The implementation requires exact equality with both directions of every complete Euclidean candidate pair: no missing, extra, duplicate, or nonreciprocal direction is accepted. It then filters only `OUT_OR_VOID`, rebuilds per-row `n_neighbors`, checks positive degree, and independently checks the complete directed included-mask set (`initial_bond_policy.jl:868-965`). `UNKNOWN`, empty candidate set, empty included mask, and any zero-degree point reject construction. Geometry exclusions do not update damage, `bond_active`, precrack state, or physical-event counts.

`BondSystem` retains an immutable `InitialBondAdmissionReceipt` with actual artifact hashes, candidate/mask hashes and counts, the measured native directed candidate/included counts, and explicit unverified graph-readback / unauthorized initialization fields. The receipt describes admission at the construction seam; it is not a native post-initialization readback or mechanical acceptance result.

## Review status and remaining boundary

The Julia regression source is in `tests/test_initial_bond_policy.jl`; it was not run. The root batch is the only authorized test execution. `isolated_fork_compile_test_plan.md` describes a later disposable-copy source test, capped at 180 seconds, with no solver or production execution. This branch is source-only and does not claim the native application has been compiled, applied, exercised, or independently read back.

The classifier result and sidecar are bound to the transport through actual SHA-256 pins. The Python compiler/job receipt is a required trust input because Julia does not parse JSON. Static byte identity alone does not independently establish classifier execution, freshness of a test run, or scientific truth of the mask. Those statements remain outside this overlay's proof boundary.
