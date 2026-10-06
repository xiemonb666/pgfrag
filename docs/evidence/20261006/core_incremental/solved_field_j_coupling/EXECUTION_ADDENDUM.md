# Execution addendum: quadrature metadata mapping

The single authorized readback produced a complete worker stage, but the parent
kept the promoted result `UNKNOWN_WORKER_RESULT_CONTRACT`. The receipt records
worker return code 0 and an outer wall time of 0.422 seconds. The 649,133-byte
stage is preserved at
`runs/run01/.solved_field_j_coupling.json.34676.worker.json` with SHA-256
`1377e428a12f53bbc9dc9537c4dbf5d0c3f67a937bcdeceee6ddc6c98dfdfac9`.

The rejection came from the metadata mapper: it read `quadrature_rule` from
each frozen solve level, where that field does not exist, so it emitted `null`.
The frozen source stores `theta_order` and `volume_point_count`. Its two
`OK` levels contain `(theta_order, volume_point_count)` values `(8, 128)` and
`(16, 256)`, which map structurally to `GL4x8x4` and `GL4x16x4`. The validator
reported `semantic_coefficient_binding_missing` because this metadata check
shares a compound rejection branch; the staged payload does contain the
semantic coefficient bindings.

`metadata_mapping_fix.py` is an 86-line, source-only candidate helper. It checks
the frozen input hash and identity, validates those source fields, and derives
the descriptors. It imports no provider or J-integrator module and does not
call the static solver. A read-only structure check returned both expected
descriptors. The helper has not been applied to the frozen adapter.

The parent result remains `UNKNOWN`; the original worker stage, run receipt,
frozen solve result, adapter source, and source-freeze record remain in place.
No second provider or J run was performed. Machine-readable hashes and run
details are in [EXECUTION_ADDENDUM.json](EXECUTION_ADDENDUM.json).
