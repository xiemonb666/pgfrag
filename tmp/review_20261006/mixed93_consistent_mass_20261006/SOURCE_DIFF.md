# Source-only delta for root review

The consistent-mass primitive is unchanged from the first source freeze. The saved-only audit now includes pointwise fitted-versus-physical checks for all six rigid velocity fields, with translation errors dimensionless and rotational errors normalized by the fixed geometry length scales. The maximum normalized error uses the declared `1e-10` gate and participates in the aggregate diagnostic result.

The audit now calls a standard-library-only recovery metadata gate before its lazy NumPy import. It verifies the independent recovery JSON pin, failed NPZ path/hash/byte count, 93-DOF and 99-order metadata, diagnostic-only recovery status, and the original producer's retained `rc=1` serialization failure and absent sibling `result.json`. The case's `result_path` and hash identify the independent recovery JSON as a readback alias; they do not claim that the producer emitted its missing result file.

The fixed case pins recovery audit `pgfrag/tmp/review_20261006/mixed93_saved_independent/audit.json` at SHA-256 `b622b771f092768fc1dd73ef5fc4a7360b0f17329ff592639db48e5600ecd495` and failed arrays `pgfrag/tmp/review_20261006/mixed93_cached_lefm_solve_20261006/run01/arrays.npz` at SHA-256 `3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41` (1,569,232 bytes). Execution remains unauthorized. No NPZ payload, NumPy import, mass assembly, solver, or pytest was run during this source-only update.

The prior frozen audit SHA-256 was `7f6698c05c7e7702b2d60498e93946d7609e70a96520a5659e8caa5358508309`. The core module SHA-256 remains `8887a9efbdf19a20b0d7890cff1acadeacc5923aecfb3b4f8695d24df33dfda2`.
