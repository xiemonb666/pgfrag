# Post-run acceptance scope addendum — 2026-10-05

This note supplements the source-frozen `static_review.md`. That review records the pre-run state; its pending-run statements are superseded by the authorized execution receipt below. The original review and source-pin manifests are left unchanged.

## Runtime result and interpreter note

Root reviewed the actual generator/Julia receipt and the native graph text, and accepts the narrow fixture-only `Body`/`BondSystem` assertions plus the direct OSB weighted-volume formula parity. The generator exited 0 in 172 ms and Julia 1.10.12 exited 0 in 27,004 ms, within the 60-second external limit. Total pilot scratch artifacts were 940,943 bytes, within the 20 MiB artifact limit. Full argv, actual child return codes, stdout/stderr, input-byte digests, and before/after source/depot integrity records are preserved at `pgfrag/tmp/initial_bond_policy_body_pilot_20261005/run_native_attempt_01/audit/pilot_execution_result.json`.

The executed generator used the bundled Codex app Python executable at `C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe` with `-B`, rather than the planned VF Python executable. It ran only the checked-in standard-library generator; the executable SHA-256 is `372c2eae555b344520bf147be0096e009069aeca4e7f78d6aecea6d53158056a`. Root accepted this recorded interpreter deviation; no rerun is needed. The native driver SHA-256 was `038f289af08b91f5d4cbee599aaab0e1e37061d5f1615316c62d15343faae4bf`.

The successful readback establishes this fixture's raw candidate equality, filtered `BondSystem` graph assertions, and native `calc_weighted_volume` formula agreement. It remains tagged `GEOMETRY_FILTER_FIXTURE_ONLY`: no CAD input, CAD classifier validation, `BodyChunk`/Threads handler, solver, force, physical event, or production job was run. The static transport receipt intentionally still says `native_graph_readback=NOT_RUN`; the separate fixture-run receipt documents the observed run.

## Interpretation limits

The fixture supplies a uniform nodal volume of `1e-9 m³` to each of the eight points in each box. Each box's coordinate bounding cube is `1 mm × 1 mm × 1 mm = 1e-9 m³`, while those eight nodal weights sum to `8e-9 m³`. These are formula-test fixture weights, not a mass-qualified or volume-conservative representation of either box. Weighted-volume parity validates the native formula on the chosen graph and weights only; it does not establish physical mass, specimen density, or CAD-qualified geometry. Do not change the fixture weights as part of this annotation.

The saved graph text at `pgfrag/tmp/initial_bond_policy_body_pilot_20261005/run_native_attempt_01/bundles/positive_mask/native_graph_readback.txt` records aggregate counts and assertion status, but does not export the full native directed-edge, point-ID, bond-ID, and degree ledger. The driver compared the complete sets in memory during the successful run; the saved receipt alone cannot support a separate edge-by-edge rebuild or visualization. A later full initializer probe should export those rows, together with ordered points and volumes, for independent read-only reconstruction.

The runtime's 20 MiB limit was enforced against total scratch artifact bytes, not Julia process RSS. No Julia retry, further fixture run, or production promotion is authorized by this addendum.
