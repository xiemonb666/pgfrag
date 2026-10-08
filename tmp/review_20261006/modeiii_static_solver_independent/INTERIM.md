# Independent static-solver review — interim snapshot

**Status:** Interim only. The authoring source was still changing when inspected; this is not a frozen-source review and must be rechecked after the root agent sends an explicit source freeze.

## Immediate blocker

At the captured snapshot, `static_solve.py` pins `static_solve_case.json` to SHA-256 `7a0311fa8fa8718b85cf1a70f8cfb9b843c923a734b49a1e85cd8f1540469c3d` (`static_solve.py:46-48`), but the on-disk case file hashes to `93c51463d717ad2aa8e78bc163fb8915aaf26164f411c2b46a20d80492bda42d`. `_verify_source_pins()` compares these values and raises before numerical imports (`static_solve.py:101-115`, called at `1317`). Therefore this exact interim snapshot cannot proceed to assembly or solve.

## Snapshot identity

- `static_solve.py` SHA-256: `d5cbfa799e6a52135eef5059281e1d9729eedda379fb6e2b0a5dbbee97261286`
- `static_solve_case.json` SHA-256: `93c51463d717ad2aa8e78bc163fb8915aaf26164f411c2b46a20d80492bda42d`
- Observed local write times during capture: README 00:51:53, case JSON 00:50:59, solver 00:52:36 on 2026-10-06.
- Dependency pins inspected in the case/provider/contract/reference/MLS/operator inputs match their current expected values; the mismatch observed is the static solve case pin.

## Review boundary

This was source inspection only. No solver, tests, numerical assembly, integration, or Git operation was run. I checked the current implementation statically for semantic DOF mapping, F1 coefficient scaling, stiffness strain ordering, saddle scaling, load construction, cap handling, quadrature shapes, off-grid metric contraction, and receipt branches. The inspected code appears internally consistent on those points, but the final review is deferred until the stable source freeze because the captured case hash already disagrees with the solver's pin.
