# Energy ledger v2 pre-review packet

- Diff: `pgfrag/tmp/review_20261004/energy_ledger/v2.diff` (v1 frozen snapshot to v2 working source/test/input).
- Snapshot: `pgfrag/tmp/review_20261004/energy_ledger/freeze_v2/`; integrity metadata records live and snapshot SHA-256 values.
- Static check: Python AST parse for the module and tests plus JSON parse for the audit input succeeded.
- v2 pytest: not run, pending root review as instructed. The v1 7-pass receipt and frozen snapshot remain untouched.
- Scope remains only `energy_ledger.py`, `test_energy_ledger.py`, and this review directory.

The v2 changes reject bool numerics and non-bool cumulative flags, reject duplicate/empty/illegal required keys, reject SI conversion overflow, return `UNKNOWN_NUMERIC_RANGE` for non-finite audit arithmetic or overflow exceptions while preserving raw channel values, and prevent single-sample free-dynamics balance from passing. A known signed negative-dissipation failure remains `FAIL` even when a cumulative increment also has a numeric-range error; that uncertainty is reported separately. Existing caller-attestation, contact-channel, all-DOF constitutive-only, and material-calibration limitations remain.
