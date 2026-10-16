# Reference vessel zero-load native benchmark audit

Status: **ZERO_LOAD_NATIVE_EVOLUTION_NUMERICAL_PASS**.

This is a bounded full-size v5 diagnostic-mesh cost/numerical-control benchmark with zero applied force. The bulk/interface values are uncalibrated placeholders; exterior ancestry remains unbound. It is not a physical qualification, fracture result, training result, or full-run cost guarantee.

The native wrapper reported 46 steps at 2.173913043478261e-09 s, with 18.758411692000003 s inside the wrapper and 23.753231219 s wall time for the outer process. These are observations from this short zero-load run, not a projected duration for 50 microseconds.

The final byte count across the two benchmark folders is 148,756,872 / 268,435,456 bytes. Every check is recorded in [`native_evolution_audit.json`](native_evolution_audit.json).

The attempt is preserved as run_attempt_01. No retry, resume, loading change, or tuning was performed.
