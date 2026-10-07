# Independent saved-result audit — source only

This folder contains a second-party checker for the saved output of the fixed-crack mixed-93 cohesive bridge worker. The producer has one accepted saved result and receipt; this folder only freezes and audits those files after a separate root review. The independent numerical checker remains disabled in `pins.json` until that review issues its own GO token.

The pinned producer result is SHA-256 `e0819a792789fe8fdada86a75983e7cb8465c5d0aa78140ca93a84c4440159f2`; its outer receipt is pinned separately in `pins.json`. After reviewing this source package, root may authorize one independent read-only check by setting this folder's own GO fields and providing that token. The checker verifies the saved JSON and receipt, hashes the already-produced NPZ before loading it read-only, and independently recomputes the stated algebraic gates. It never imports the author worker, MLS, material law, or assembly code, and it has no solve/eigenvalue/provider path.

The independent checks are limited to the fixed 16-point interface trace and the saved linear algebra: one-side area `R*L=1e-4 m²`; the `F1` y-opening column at semantic coefficient 82; `KG[82,82]=20000 N/m`; six saved rigid candidates; baseline/final elastic cohesive force and energy; the face-traction sign convention; and the saved coupled KKT residual and energy identity. Raw spectrum gates are reported as producer evidence and are not independently recomputed.

The fixture uses `K0=1e8 Pa/m`, `tpeak=1000 Pa`, and `Gc=0.01 J/m²` as a numerical toy. This is DIAGNOSTIC_ONLY: fixed crack, no growth, no incremental external-work claim, and no ceramic or other material qualification. The zero-tangent control is synthetic and is not a physical history.

`check_source_contract.py` is the only check permitted before root authorization. It performs stdlib-only path/hash/AST checks, does not read the NPZ, result, or receipt, and prints `SOURCE_READY` only when the source package is internally consistent.
