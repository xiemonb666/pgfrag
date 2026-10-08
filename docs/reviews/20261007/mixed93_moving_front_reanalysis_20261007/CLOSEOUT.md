# Moving-front run01 closeout

The single authorized attempt `ROOT_20261007_MIXED93_MOVING_FRONT_RUN01` is closed. Do not rerun `capture_once.py` or reuse `run01`; any future attempt needs new authorization and a new output ID.

The capture wrapper returned **2** after one worker invocation. The worker returned **2**, did not time out, and wrote a failure-only checkpoint. The external tool call took **1.269 s**; the wrapper recorded **0.968 s** through its post-hash/output tally and a 59.938 s child timeout budget. The wrapper does not record actual child execution duration, so that duration is unknown.

The worker completed 1,344 new cap MLS rows in six batches `[256, 192, 256, 256, 256, 128]`; all statuses were `OK`. It assembled both raw 93-DOF stiffness matrices and passed the recorded 6-zero/87-positive/no-significant-negative spectrum and six-rigid checks for both fronts. The side-traction cross-check was 1.842e-16. The failure checkpoint ends at `candidate_front_rigid_checks_complete`. No boundary wrench, KKT solve, energy-release value, or primary result was produced.

The stop was a source contract mismatch, not a physics failure: the driver checks `audit.producer_run.rc`, while the frozen audit records `audit.producer_run.return_code = 1` and `result_json_present = false`. The failed check occurred before wrench and KKT work. No source fix or retry was made.

The capture receipt reports 1,802,359 bytes including the receipt index, below the 12,582,912-byte limit; independent directory tally agrees at 1,802,359 bytes. Pre/post pins matched. Raw stdout/stderr, the failure JSON, compressed checkpoint arrays, receipt and index remain in `capture_once01/` and `run01/`; the NPZ payload was not opened.

Failure JSON SHA-256: `9da8793fda4925014154931462a111fae0f71a35fbae6db48f9d6ac3447f2935`. Failure arrays NPZ SHA-256: `8f6d3ca994d44ea5003379736ae4cd91f9c4daa98af1366c8c59c81f757c0e2f`. Capture receipt SHA-256: `779b503781d711a6b06257cfbc72fd00938af34c315307bb2de948e6908240f3`. Outer addendum SHA-256: `4bb9e1307bcfcba3c7b51a8147a56c7e9293d55197b34768c2f721f21e5c6dbf`. The pre-GO disabled case and pins are preserved in `pre_GO_disabled_snapshot/`.

This diagnostic does not establish crack growth, fracture energy, or material qualification.
